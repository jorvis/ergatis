#!/usr/bin/perl

=head1 NAME

enzyme2obo.pl - Reads an ExPASy release of Enzyme Commission data and converts to an OBO file.

=head1 SYNOPSIS

USAGE: enzyme2obo.pl 
            --dat_file=/path/to/enzyme.dat 
            --class_file=/path/to/enzclass.txt
            --output_file=/path/to/somefile.obo
          [ --old_obo=/path/to/old_ec.obo ]

=head1 OPTIONS

B<--dat_file>
    Input enzyme.dat file, provided by ExPASy.

B<--class_file>
    Input enzclass.txt file, provided by ExPASy.

B<--output_file>
    Output OBO file to be created.  If exists, it will be overwritten.

B<--old_obo>
    Optional.  If an existing OBO file for this exists the EC: IDs for terms will be maintained.

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

This script was written to provide a way of representing EC as an ontology, with the intent
of loading into a Chado database instance to support genome annotation.  Because of some
limitations in the cvterm table (or at least our usage of it) some modifications to the
source data had to be made.  For details see the OUTPUT section.

=head1  INPUT

The currently  used line  types, along with their respective line codes, are listed below:

   ID  Identification                         (Begins each entry; 1 per entry)
   DE  Description (official name)            (>=1 per entry)
   AN  Alternate name(s)                      (>=0 per entry)
   CA  Catalytic activity                     (>=1 per entry)
   CF  Cofactor(s)                            (>=0 per entry)
   CC  Comments                               (>=0 per entry)
   PR  Cross-references to PROSITE            (>=0 per entry)
   DR  Cross-references to Swiss-Prot         (>=0 per entry)
   //  Termination line                       (Ends each entry; 1 per entry)

Note that if you're pointing this at a very old ec.obo file you'll get errors like:

    A format problem has been detected (and ignored) in ...
    
This may have to do with incorrect encoding of the 'def:' lines in the file.  They should
each end with [] symbols.  If not, fix it like this:

    perl -i -pe 's|^(def\:.+)|$1 \[\]|' ec.obo.old

=head1 MAPPING

    DAT file code           OBO file attribute
    -------------           ------------------
        ID                      xref_analog
        DE                      name
        AN                      synonym
        CA                      def
        CC                      comment

All other code values are not currently mapped.

The full specification for the enzyme.dat file is provided here:

    ftp://ftp.expasy.org/databases/enzyme/enzuser.txt

The OBO Flat File format followed is defined here:

    http://www.geneontology.org/GO.format.obo-1_2.shtml#S.2.2.3

=head1  OUTPUT

The output conforms to version 1.2 of the OBO specification, with example entries like:

    [Term]
    id: EC:00000002
    name: UTP--xylose-1-phosphate uridylyltransferase.
    def: "UTP + alpha-D-xylose 1-phosphate = diphosphate + UDP-D-xylose." []
    synonym: "Xylose-1-phosphate uridylyltransferase." RELATED []
    xref: EC:2.7.7.11
    is_a: EC:00001637 ! Transferases. Transferring phosphorous-containing groups. Nucleotidyltransferases.


Note that the values for all synonym scopes are 'RELATED'.  Because the enzyme.dat file doesn't
have an equivalent value for this I choose what I thought to be the safest assertion.  In 
reality, the scope 'EXACT' is probably the correct value for most.

Note that parental relationships (via is_a) are not defined for deleted or transferred entries.
This is because the parent entries may be completely removed from the collection, as was done
with 3.4.4.-

=head2 Name value modifications

Because the Chado schema restricts cvterm.name to be unique within a given controlled
vocabulary and is_obsolete status, the values of the DE lines in the enzyme.dat need to be
manipulated somewhat.

Because their definitions are curated to be fully descriptive and unique, full precision 
terms have names exported that are identical to their DE lines.  Each level above this, on
the other hand, has a name formed by the concatenation of their parents.  This conforms to 
the encoding strategy (manually curated) of these in the Gene Ontology.

Entries annotated as 'Deleted entry.' within the input set have the original EC number
appended to the end of the name in parentheses.

Entries annotated as Transferred entry: Y are renamed like Transferred entry: (X -> Y), 
where X is the original EC number and Y is the new one.

=head1  CONTACT

    Joshua Orvis
    jorvis@som.umaryland.edu

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

use lib '/usr/local/packages/perllib/';
use OBO::Core::Ontology;
use OBO::Core::Term;
use OBO::Core::Relationship;
use OBO::Core::RelationshipType;
use OBO::Parser::OBOParser;

my %options = ();
my $results = GetOptions (\%options, 
                          'dat_file=s',
                          'class_file=s',
                          'output_file=s',
                          'old_obo=s',
                          'log|l=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

## how many digits to use (for zero-padding) in EC IDs? (EC:00002464, for example)
my $EC_ID_WIDTH = 8;

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">$options{log}") || die "can't create log file: $!";
}

my $old_obo;
my $old_ec_id_lookup;
my $new_ec_id_lookup;

## store the max EC ID so we know what to increment from
#   for IDs like: EC:00000008
my $max_ec_id = 0;

## store any pre-defined OBO information
if ( defined $options{old_obo} ) {
    my $obo_parser = OBO::Parser::OBOParser->new();
    
    $old_obo = $obo_parser->work( $options{old_obo} );
    
    ## store an EC number -> old EC ID lookup
    for my $term ( @{$old_obo->get_terms()} ) {
        for my $xref ( $term->xref_set_as_string ) {
            if ( $xref->name =~ /^EC\:/ ) {
                $$old_ec_id_lookup{ $xref->name } = $term->id;
                
                if ( $term->id =~ /^EC:(\d+)/ ) {
                    if ( $1 > $max_ec_id ) {
                        $max_ec_id = $1;
                    }
                } else {
                    die "expected previous EC ID to be in a format like: EC:00000006.  Instead, got" . $term->name;
                }
            }
        }
    }
}

## create a new OBO
my $ontology = OBO::Core::Ontology->new();

$ontology->remark('Created from ExPASy enzyme.dat by enzyme2obo.pl');
$ontology->default_namespace('EC');

## add any relationship types
$ontology->add_relationship_type_as_string( 'is_a', 'is_a' );

## load the enzyme class defs
open( my $classfh, "<$options{class_file}") || die "failed to open class file: $!";

my $entry;
my ($class_id, $class_value);

while ( my $line = <$classfh> ) {
    chomp $line;

    ## skip whitespace lines
    next if ( $line =~ /^\s*$/ );

    ## matches patterns like: 1. 1.99.-    With other acceptors.
    if ( $line =~ m|^([ \d\.\-]{1,3}\.[ \d\.\-]{1,3}\.[ \d\.\-]{1,3}\.[ \d\.\-]{1,3})\s+(.+)| ) {
        
        if ( defined $class_id ) {
            add_enzyme_to_obo( { id => $class_id, de => [ $class_value ] } );
        }
        
        ($class_id, $class_value) = ($1, $2);
        
        ## make sure the ID doesn't have any spaces
        $class_id =~ s|\s||g;
        
        ## unless this is top-level, prepend that name with the parent's
        my $parent_num = parent_ec_num( $class_id);
        
        if ( $parent_num ) {
            my $parent_term = $ontology->get_term_by_xref( 'EC', $parent_num );
            
            $class_value = $parent_term->name . " $class_value";
            #LEFT OFF HERE
        }
    
    } elsif ( $line =~ m|^ {9,}\s*(.+)| ) {
        $class_value .= " $1";
    }
}

## make sure to add the last one
add_enzyme_to_obo( { id => $class_id, de => [ $class_value ] } );

## open the input file
open(my $ifh, "<$options{dat_file}") || die "failed to read input DAT file: $!";

$entry = undef;

## $h->{old_ec_num} -> new_ec_num
my $replaced_id_lookup = {};

while ( my $line = <$ifh> ) {
    chomp $line;
    
    ## first check for the end of a record
    if ( $line =~ m|^//| && defined $$entry{id} ) {
        add_enzyme_to_obo( $entry );
        undef $entry;
    }
    
    ## all lines are structured like this:
    #    Characters    Content
    #    ---------     ----------------------------------------------------------
    #    1 to 2        Two-character line code. Indicates the type of information
    #                  contained in the line.
    #    3 to 5        Blank
    #    6 up to 78    Data
    my ($code, $data);

    if ( $line =~ /(..)...(.+)/ ) {
        ( $code, $data ) = ( $1, $2 );
    } else {
        next;
    }
    
    if ( $code eq 'ID' ) {

        $$entry{id} = $data;
        
    } elsif ( $code eq 'DE' ) {
    
        push @{$$entry{de}}, $data;
    
    } elsif ( $code eq 'AN' ) {
        
        push @{$$entry{an}}, $data;
    
    } elsif ( $code eq 'CA' ) {
    
        push @{$$entry{ca}}, $data;
    
    } elsif ( $code eq 'CC' ) {
        
        push @{$$entry{cc}}, $data;
    }
    
}

## some terms were replaced by others.  make those assertions now
#   (doing it now prevents a double parse/scan of the document)
for my $deprecated_num ( keys %$replaced_id_lookup ) {
    my ($replacement_id, $deprecated_id);
    my $replacement_num = $$replaced_id_lookup{$deprecated_num};

    if ( defined $$old_ec_id_lookup{$deprecated_num} ) {
        $deprecated_id = $$old_ec_id_lookup{$deprecated_num};
        
    } elsif ( defined $$new_ec_id_lookup{$deprecated_num} ) { 
        $deprecated_id = $$new_ec_id_lookup{$deprecated_num};
        
    } else {
        die "couldn't find $deprecated_num in my lookups.";
    }
    
    if ( defined $$old_ec_id_lookup{ $replacement_num } ) {
        $replacement_id = $$old_ec_id_lookup{$replacement_num};
        
    } elsif ( defined $$new_ec_id_lookup{$replacement_num} ) { 
        $replacement_id = $$new_ec_id_lookup{$replacement_num};
        
    } else {
        die "couldn't find $replacement_num in my lookups.";
    }
    
    my $ec_term = $ontology->get_term_by_id( $deprecated_id );
    
    if ( defined $ec_term ) {
        $ec_term->replaced_by( $replacement_id );
    } else {
        die "$deprecated_num was tagged for replacement but I couldn't find it in my ontology data structure.";
    }
}

## now add parental relationship
for my $term ( @{$ontology->get_terms()} ) {
    ## relationships are not defined for transferred and deleted entries (since their parents
    #   might have been completely removed (ex, 3.4.4.-)
    next if ( $term->name =~ /^Deleted entry\./ || $term->name =~ /^Transferred entry\: / );
    
    ## get the EC number
    my $ec_num;
    
    for my $xref ( $term->xref_set_as_string ) {
        if ( $xref->name =~ /^EC\:(.+)/ ) {
            $ec_num = $1;
            last;
        }
    }
    
    if (! defined $ec_num ) {
        die "there was a problem getting the EC xref from term with ID: " . $term->id;
    }
    
    my $parent_ec_num = parent_ec_num( $ec_num );
    my $parent_ec_id;
    
    ## next if we're at the top level
    next unless defined $parent_ec_num;
    
    if ( defined $$old_ec_id_lookup{ "EC:$parent_ec_num" } ) {
        $parent_ec_id = $$old_ec_id_lookup{ "EC:$parent_ec_num" };
    
    } elsif ( defined $$new_ec_id_lookup{ "EC:$parent_ec_num" } ) {
        $parent_ec_id = $$new_ec_id_lookup{ "EC:$parent_ec_num" };
    
    } else {
        die "failed to create a relationship because I couldn't find an EC ID for $parent_ec_num\n";
    }
    
   my $parent_term = $ontology->get_term_by_id($parent_ec_id);
   $ontology->create_rel( $term, 'is_a', $parent_term );
}

sub parent_ec_num {
    my $child = shift;
    my $parent;
    
    my @ec_num_parts = split('\.', $child);
    
    ## this is easier backwards
    my @new_parts = ();
    @ec_num_parts = reverse ( @ec_num_parts );
    
    while ( scalar @ec_num_parts ) {
        my $part = shift @ec_num_parts;
        
        ## if ec_num_parts is now empty, we have a top-level term
        if ( ! scalar @ec_num_parts ) {
            return undef;
        
        } elsif ( $part ne '-' ) {
            push @new_parts, '-';
            push @new_parts, @ec_num_parts;
            undef @ec_num_parts;
        } else {
            push @new_parts, '-';
        }
    }
    
    if ( scalar @new_parts ) {
        $parent = join('.', reverse @new_parts);
    }
    
    return $parent;
}

open(my $ofh, ">$options{output_file}") || die "failed to create output file: $!";

$ontology->export( $ofh );

exit(0);


sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
}


sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( dat_file class_file output_file );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
    
    ##
    ## you can do other things here, such as checking that files exist, etc.
    ##
    
    ## handle some defaults
    #$options{optional_argument2}   = 'foo'  unless ($options{optional_argument2});
}

sub add_enzyme_to_obo {
    my $atts = shift;
    
    my $term = OBO::Core::Term->new();
    my $id;
    my $is_obsolete = 0;
    
    if ( defined $$old_ec_id_lookup{ "EC:$$atts{id}" } ) {
        $id = $$old_ec_id_lookup{ "EC:$$atts{id}" }
        
    } else {
        $id = sprintf("EC:%0${EC_ID_WIDTH}d", $max_ec_id++);
        
        ## remember this one
        $$new_ec_id_lookup{"EC:$$atts{id}"} = $id;
    }
    
    $term->id( $id );
    
    ## put the EC number as an xref
    $term->xref_set_as_string("[EC:$$atts{id}]");
    
    ## handle names
    if ( ! defined $$atts{de} ) {
        die "failed to find a name (required) for term with id $$atts{id}\n";
    } 
    
    ## the spec only allows one name per term, so if there are multiple lines they
    #   should be concatenated.
    my $name = join(' ', @{$$atts{de}});
    $term->name( $name );
    
    ## the DE lines also contain info on obsoleted entries
    if ( $term->name eq 'Deleted entry.' ) {
        $is_obsolete = 1;

        ## rename it for uniqueness
        $term->name( "Deleted entry. ($$atts{id})" );

    } elsif ( $term->name =~ /Transferred entry\: (\d+\.\d+\.\d+\.\d+)/ ) {
        
        $$replaced_id_lookup{"EC:$$atts{id}"} = "EC:$1";
        $is_obsolete = 1;
        
        ## rename it for uniqueness
        $term->name( "Transferred entry: ($$atts{id} -> $1)" );
    }
    
    if ( $is_obsolete ) {
        $term->is_obsolete(1);
    }
    
    ## handle synonyms
    if ( defined $$atts{an} ) {
        for my $syn_str ( @{$$atts{an}} ) {
            $term->synonym_as_string( $syn_str, "[]", 'RELATED' );
        }
    }
    
    ## handle the definition (interactions)
    if ( defined $$atts{ca} ) {
        $term->def_as_string( join(' ', @{$$atts{ca}}), "[]" );
    }
    
    ## now any comments
    if ( defined $$atts{cc} ) {
        my $comment = join(' ', @{$$atts{cc}});
        
        ## collapse multiple spaces
        $comment =~ s/ +/ /g;
        $term->comment( $comment );
    }
    
    $ontology->add_term( $term );
}






























