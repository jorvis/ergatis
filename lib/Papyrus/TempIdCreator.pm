package Papyrus::TempIdCreator;

=head1 NAME

TempIdCreator.pm - A module for creating temporary IDs within BSML documents.

=head1 SYNOPSIS

    use Papyrus::TempIdCreator;

    my $creator = new Papyrus::TempIdCreator;

    print "the so type used for primary_transcript is " .
        $creator->so_used('primary_transcript') . "\n";

    print "the next ID for so_type primary_transcript is " .
        $creator->new_id(db => 'aa1', so_type => 'primary_transcript', prefix => '30') . "\n";

    my $id = $creator->next_id_num('transcript');
    print "next id for transcript is $id\n";

    $id = $creator->next_id_num('transcript');
    print "next id for transcript is $id\n";


    prints:

    the so type used for primary_transcript is transcript
    the next ID for so_type primary_transcript is ir.aa1.transcript.30.1
    next id for transcript is 2
    using id transcript 3

=head1 DESCRIPTION

Though not tied to anything BSML, this script is intended to both provide a
method for creating incrementing IDs based on SO type and give a method
to look up a SO type to see its current usage status, which will probably be done
most often while creating BSML files.

For example, we may not currently use the SO term "primary_transcript", instead
preferring the more general "transcript".  We can write each of our BSML creation
scripts to lookup our usage for "primary_transcript" here, and then use the
term returned (so_used method).  Then, when we decide to support "primary_transcript",
we only have to change its mapping here and all our scripts that build BSML will
be up to date.

The primary method for ID generation is the "new_id" method.  It requires 'db' and
'so_type' as arguments and will return a full temporary ID for use in your BSML
document.  Other methods in this module are provided but may not often be used
directly.

You can pass a SO_term to the "next_id_num" method, which will return an incremented
counter specific to that term.  It will croak if that SO term is not in its index,
which can be checked (in a non-croaking way) using the indexes method.

=head1 METHODS

=over 3

=item I<PACKAGE>->new()

Returns a newly created "Papyrus::TempIdCreator" object.  No arguments are necessary.
More than one object can be created if, for example, you are creating two different
BSML documents within the same script and want separate counters for each.

=item I<$OBJ>->new_id( db => I<'dbid'>, so_type => I<'so type'>, [ prefix => I<'prefix'> ])

Perhaps the only method you'll use in this package, this generates fully-qualified
temporary IDs according to our current naming scheme.  The required arguments are
'db' and 'so_type', along with the optional 'prefix'.  Names generated are returned
according to the following convention:

    ir.$DB.$SO_TYPE.$PREFIX.N
    
Where N is an incrementing integer for that so type.  For example:

    ir.aa1.gene.38556579.3
    
The prefix is optional, and is used to group the temporary IDs and help make them unique
across documents.  In practice the prefix is usually the Workflow command id for the
command generating the new features, such as a foo2bsml.pl script.  Thus, all of the temp
IDs within the same document will usually have the same prefix.

The so_type passed to this method will first be checked and changed to the currently
implemented so_type, if necessary.

=item I<$OBJ>->indexes(I<'so type'>)

Returns 1 (true) if the SO term passed is among those indexed in this module.  Returns
false otherwise and, unlike other methods, does not croak.

=item I<$OBJ>->next_id_num(I<'so type'>)

Returns an integer of the next available id for the SO term passed.  Croaks if the SO
term passed is not indexed by this module.

=item I<$OBJ>->so_used(I<'so term'>)

Returns the text string of the SO term currently used in place of the passed term.  In
many cases, these may be the same.  Croaks if the SO term passed is not indexed by
this module.

=back

=head1 AUTHOR

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;
use warnings;
use Carp;

## class data and methods
{
    ## here we put the mapping table of SO types with the type we
    ##  currently use for each.
    my %_so_used = &_load_terms();
    my %_ids = ();
    
    sub new {
        my ($class, %args) = @_;
        
        return bless { %_ids }, $class;
    }

    ## generate a full temporary id with format like ir.aa1.gene.50.234 or ir.aa1.gene.234
    ##  depending on whether a prefix is passed
    sub new_id {
        my ($self, %args) = @_;
        
        ## user should have passed a 'db' and 'so_type'.  a 'prefix' is optional
        if (! exists $args{db}) { croak "db is a required argument to the new_id method!" }
        if (! exists $args{so_type}) { croak "so_type is a required argument to the new_id method!" }
        
        ## look up the SO term we are actually using for this term
        my $so = $self->so_used($args{so_type});
         
        my $newid = "ir.$args{db}.$so.";
        
        if (exists $args{prefix}) {
            $newid .= "$args{prefix}." . $self->next_id_num( $so );
        } else {
            $newid .= $self->next_id_num( $so );
        }
        
        return $newid;
    }

    ## check if so term is indexed, return true or false.
    sub indexes {
        my ($self, $qry) = @_;
        
        return 1 if (exists $_so_used{$qry});
        return 0;
    }

    ## return the next id num for a given so term
    sub next_id_num {
        my ($self, $qry) = @_;
        
        ## make sure this is a term we have indexed
        if ( $self->indexes($qry) ) {
            return ++$self->{$qry};
        } else {
            croak ("so id $qry is not in my lookup table!");
        }
    }

    ## return the actual so term used for any given so term.
    sub so_used {
        my ($self, $qry) = @_;
        
        if (exists $_so_used{$qry}) {
            return $_so_used{$qry};
        } else {
            croak ("so id $qry is not in my lookup table!");
        }
    }
    
    ## load the so terms from the data table at the end of this module
    sub _load_terms {
        my %types;
        
        while (<DATA>) {
            if (/^\s*(.+?)\s+(.+)\s*$/) {
                $types{$1} = $2;
            }
        }
        
        return %types;
    }
}

1;

__DATA__

canonical_five_prime_splice_site        splice_site
canonical_splice_site                   splice_site
canonical_three_prime_splice_site       splice_site
CDS                                     CDS
cds                                     CDS
coding_exon                             exon
exon                                    exon
five_prime_UTR                          five_prime_UTR
gene                                    gene
intergenic_region                       intergenic_region
intron                                  intron
noncoding_exon                          exon
polyA_signal_sequence                   polyA_signal_sequence
primary_transcript                      transcript
promoter                                promoter
repeat_region                           repeat_region
signal_peptide                          signal_peptide
splice_site                             splice_site
start_codon                             start_codon
stop_codon                              stop_codon
tandem_repeat                           repeat_region
three_prime_UTR                         three_prime_UTR
transcript                              transcript
transcription_start_site                transcription_start_site
transcription_end_site                  transcription_end_site
transit_peptide                         transit_peptide
tRNA                                    tRNA
SNP										SNP
snp										SNP
nucleotide_insertion					nucleotide_insertion
nucleotide_deletion						nucleotide_deletion
indel									indel
gap										gap
