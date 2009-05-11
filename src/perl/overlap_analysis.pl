#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use XML::Twig;
use Pod::Usage;
use File::OpenFile qw (open_file);
use IntervalTree;
use Data::Dumper;
use BSML::FeatureRelationshipLookup;

###################### Globals #####################
my $logfh = \*STDOUT;
my $input_file;
my $output;
my @rna_bsml;
my @evidence_bsml;
my $flagfh = \*STDOUT;
my $cutoff = 60;
my $sequence_class = "assembly";
my %delete;
##################################################

&check_options();

#parse all the genes from the input file
print "Parsing input file\n";
my %genes = &get_genes( $input_file );

#get the relationships between feature types
my $feature_rel_lookup = new BSML::FeatureRelationshipLookup( 'bsml' => $input_file );

#parse all the rnas 
my @rnas;
print "Parsing RNA files\n";
foreach my $rna_file ( @rna_bsml ) {
    push(@rnas, &get_rnas( $rna_file ) );
}
#parse evidence files
print "Parsing evidence files\n";
my %evidence = &parse_evidence_files( \@evidence_bsml );

#create the IntervalTrees
my %iTrees;

foreach my $gene_id ( keys %genes ) {
    my $seq_id = $genes{$gene_id}->{'parent_seq'};
    $iTrees{$seq_id} = new IntervalTree if( !defined( $iTrees{$seq_id} ) );
    $iTrees{$seq_id}->addInterval( $gene_id, $genes{$gene_id}->{'left'}, $genes{$gene_id}->{'right'} );
}

print "Building Trees\n";
map { $iTrees{$_}->buildTree } keys %iTrees;

#searching for gene overlaps
foreach my $gene_id ( keys %genes ) {
    my $seq = $genes{$gene_id}->{'parent_seq'};
    my @overlaps = $iTrees{$seq}->searchInterval( $genes{$gene_id}->{'start'},
                                                  $genes{$gene_id}->{'stop'} );

    foreach my $ol ( @overlaps ) {

        #determine the amount of overlap
        my @ordered = sort ( $genes{$gene_id}->{'left'}, $genes{$gene_id}->{'right'}, $ol->[0], $ol->[1] );
        my $ol_amount = (( $genes{$gene_id}->{'right'} - $genes{$gene_id}->{'left'} ) + ( $ol->[1] - $ol->[0] )) -
            ( $ordered[3] - $ordered[0] );

        #does it pass the cutoff?
        if( $ol_amount > $cutoff ) {
            &handle_gene_overlap( $ol, [ $genes{$gene_id}->{'left'}, $genes{$gene_id}->{'right'}, $gene_id ] );
        } else {
            &_log("Overlap (size: $ol_amount) found between $gene_id and $ol->[2]. Skipping because it does not".
                  "meet the maximum cutoff (cutoff: $cutoff)");
        }
    }
}

#searching for RNA overlaps
foreach my $rna ( @rnas ) {
    my ($id, $start, $stop, $seq) = @{$rna};
    next if( !exists( $iTrees{$seq} ) );
    my @overlaps = $iTrees{$seq}->searchInterval( $start, $stop );

    foreach my $ol ( @overlaps ) {
        &handle_rna_overlap( $rna, $ol );
    }

}
 
my $outfh = open_file( $output, 'out' );
my $twig = new XML::Twig( 'twig_roots' => {
    'Sequence' => sub {
        my $id = $_[1]->att('id');
        if( $id =~ /^(.*)_seq$/ ) {
            my $real_id = $1;
            my $gene = $feature_rel_lookup->lookup( $real_id, 'gene' );
            if( exists( $delete{$gene} ) ) {
                $_[1]->delete;
            } else {
                $_[1]->print;
            }
        } else {
            $_[1]->print;
        }
        $_[0]->purge;
    },
    'Feature' => sub {
        my $id = $_[1]->att('id');
        my $gene = $feature_rel_lookup->lookup( $id, 'gene' );
        if( exists( $delete{$gene} ) ) {
            $_[1]->delete;
        }
    },
    'Feature-group' => sub {
        my $id = $_[1]->att('group-set');
        my $fgm = $_[1]->first_child('Feature-group-member');
        my $gene_id = $feature_rel_lookup->lookup( $fgm->att('featref'), 'gene' );
        if( exists( $delete{$gene_id} ) ) {
            $_[1]->delete;
        }
    },
    'Analysis' => sub {
        $_[1]->set_att('id', "overlap_analysis_analysis");
        map { $_->delete } $_[1]->children('Attribute');
        my $att = new XML::Twig::Elt('Attribute');
        $att->set_atts( { 'name' => 'sourcename',
                         'content' => $output });
        $att->paste('first_child', $_[1] );
        $att = new XML::Twig::Elt('Attribute');
        $att->set_atts( {'name' => 'program',
                        'content' => 'overlap_analysis' });
        $att->paste('first_child', $_[1] );
        $att = new XML::Twig::Elt('Attribute');
        $att->set_atts( {'name' => 'algorithm',
                         'content' => 'overlap_analysis' } );
        $att->paste('first_child', $_[1] );
                        
        $_[0]->flush;
    }
},
                          'twig_print_outside_roots' => $outfh,
                          'pretty_print' => 'indented' );
my $in = open_file( $input_file, 'in' );
$twig->parse( $in );
close($in);
close($outfh);



################################### SUBS ###################################
sub handle_gene_overlap {
    my ($gene1, $gene2) = @_;
    
    #does gene1 have evidence?
    my $ev1 = 0;
    my $ev2 = 0;
    $ev1 = 1 if( $evidence{$gene1->[2]} );
    $ev2 = 1 if( $evidence{$gene2->[2]} );

    if( $ev1 && $ev2 ) {
        &flag_overlap( $gene1->[2], 'gene', $gene1->[0], $gene1->[1],
                       $gene2->[2], 'gene', $gene2->[0], $gene2->[1] );
    } elsif( $ev1 || $ev2 ) {
        if( $ev1 ) {
            &mark_for_delete( $gene1->[2] );
            &_log( "Marking $gene1->[2] for delete because it has no evidence. Overlaps $gene2->[2] ".
                   "($gene1->[0] - $gene1->[1] and $gene2->[0] - $gene2->[1] )" );
        } elsif( $ev2 ) {
            &mark_for_delete( $gene2->[2] );
            &_log( "Marking $gene2->[2] for delete because it has no evidence. Overlaps $gene1->[2] ".
                   "($gene2->[0] - $gene2->[1] and $gene1->[0] - $gene1->[1] )" );
        }
    } else {
        &flag_overlap( $gene1->[2], 'gene', $gene1->[0], $gene1->[1],
                       $gene2->[2], 'gene', $gene2->[0], $gene2->[1], "no_evidence" );
    }
}

sub handle_rna_overlap {
    my ($rna, $gene) = @_;
    
    #does the gene have evidence?
    my $gene_id = $gene->[2];
    
    if( exists( $evidence{$gene_id} ) ) {
        &flag_overlap( $gene_id, 'gene', $gene->[0], $gene->[1], $rna->[0],
                       'rna', $rna->[1], $rna->[2] );
    } else {
        &mark_for_delete( $gene_id );
        &_log("Marking $gene_id for deletion because it overlaps an RNA prediction and has no evidence ".
              "(RNA gene: $rna->[0]) ( $gene->[0] - $gene->[1] and $rna->[1] - $rna->[2] )");
    }
    
}

sub mark_for_delete {
    my ($gene_id) = @_;
    $delete{$gene_id} = 1;
}

sub flag_overlap {
    my (@fields) = @_;
    my $string = join("\t", @fields);
    print $flagfh $string."\n";
}

sub parse_evidence_files {
    my ($files) = @_;
    my %retval;

    my $twig = new XML::Twig( 'twig_roots' => {
        'Seq-pair-alignment' => sub {
            my $refseq = $_[1]->att('refseq');
            my $gene_id;
            eval {
                $gene_id = $feature_rel_lookup->lookup( $refseq, 'gene' );
            };
            return unless( defined( $gene_id ) );
            $retval{$gene_id} = 1;
        }
    } );
    
    foreach my $file (@{$files}) {
        my $in = open_file( $file, "in" );
        $twig->parse( $in );
        close($in);
    }

    return %retval;
}

sub get_rnas {
    my ($file) = @_;
    my @retval;
    
    my $twig = new XML::Twig( 'twig_roots' => {
        'Sequence[@class="'.$sequence_class.'"]' => sub {
            my $seq_id = $_[1]->att('id');
            my @features = $_[1]->find_nodes('//Feature[@class="gene"]');
            foreach my $feature ( @features ) {
                my $gene_id = $feature->att('id');
                my $int_loc = $feature->first_child('Interval-loc');
                my ($start, $stop) = ($int_loc->att('startpos'), $int_loc->att('endpos') );
                push( @retval, [ $gene_id, $start, $stop, $seq_id ] );
            }
        } } );

    my $in = open_file( $file, "in" );
    $twig->parse( $in );
    close($in);

    return @retval;
}

sub get_genes {
    my ($file) = @_;
    my %retval;

    my $twig = new XML::Twig( 'twig_roots' => {
        'Sequence[@class="'.$sequence_class.'"]' => sub {
            my $seq_id = $_[1]->att('id');
            my @features = $_[1]->find_nodes('//Feature[@class="gene"]');
            foreach my $feature ( @features ) {
                my $feat_id = $feature->att('id');
                my $int_loc = $feature->first_child('Interval-loc');
                my ($left,$right,$comp) = ( $int_loc->att('startpos'),
                                            $int_loc->att('endpos'),
                                            $int_loc->att('complement') );
                $retval{$feat_id} = {
                    'left' => $left,
                    'right' => $right,
                    'comp' => $comp,
                    'parent_seq' => $seq_id,
                    };
            }
            
        } } );

    my $in = open_file( $file, 'in' );
    $twig->parse( $in );
    close($in);
    return %retval;
}

sub check_options {
    my %options;
    my $results = GetOptions (\%options,
                              'input_file|i=s',
                              'output|o=s',
                              'rna_bsml|r=s',
                              'evidence_bsml|e=s',
                              'overlap_cutoff|c=s',
                              'flagged_overlaps_file|f=s',
                              'input_sequence_class|s=s',
                              'log|l=s',
                              'help|h'
                              );

    &_pod if( $options{'help'} );
    $logfh = open_file( $options{'log'}, "out" ) if( $options{'log'} );

    my @reqs = qw(input_file output rna_bsml evidence_bsml);
    foreach my $req ( @reqs ) {
        die("Option $req is required") unless( $options{$req} );
    }

    $input_file = $options{'input_file'};
    $output = $options{'output'};
    @rna_bsml = &get_input_files( $options{'rna_bsml'} );
    @evidence_bsml = &get_input_files( $options{'evidence_bsml'} );
    $flagfh = open_file( $options{'flagged_overlaps_file'}, 'out' ) if( $options{'flagged_overlaps_file'} );
    $sequence_class = $options{'input_sequence_class'} if( $options{'input_sequence_class'} );
}

sub get_input_files {
    my ($input_string) = @_;
    my @retval;
    my @tokens = split(/[,\s]+/, $input_string);

    foreach my $token ( @tokens ) {
        if( $token =~ /.bsml(\.gz)?$/ ) {
            push( @retval, $token );
        } elsif( $token =~ /.list(\.gz)?$/ ) {
            my $in = open_file( $token, "in" );
            chomp( my @tmp = <$in> );
            close($in);
            push( @retval, @tmp );
        } else {
            die("Don't recognize extension of file $token.  Looking for .bsml or .list files");
        }
    }

    return @retval;
}

sub _log {
    my ($msg) = @_;
    print $logfh "$msg\n";
}

sub _pod {
    pod2usage( { 
        -exitval => 0, 
        -verbose => 2, 
        -output => \*STDERR} );
}
=head1 NAME

overlap_analysis.pl - Analyzes overlaps and attempts to resolve them

=head1 SYNOPSIS

 USAGE: overlap_analysis.pl
       --input_file=/path/to/some/gene_describing.bsml
       --output=/path/to/output.bsml
       --rna_bsml=/path/to/file1.bsml,/path/to/bsml.list
       --evidence_bsml=/path/to/evidence1.bsml,/path/to/evidence_bsml.list
     [ --flagged_overlaps_file=/path/to/flagged.file
       --overlap_cutoff=60
       --input_sequence_class=assembly
       --log=/path/to/file.log
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>
    A gene describing bsml document.  Will parse out Feature elements.

B<--output,-o>
    The output bsml file

B<--rna_bsml,-r>
    Bsml describing RNA predictions/models.

    Can take in a comma separated list of bsml files or bsml lists.  Will determine the
    format of the file based on the extension (.bsml or .list). Can be gzipped
    and therefore contain .gz extension.

B<--evidence_bsml,-e>
    Bsml that holds alignment evidence (Blast/Ber and HMM)

    Takes a comma separated list of bsml files or bsml lists. Will determine the 
    format of the file based on the extension (.bsml or .list). Can be gzipped
    and therefore contain .gz extension.

B<--flagged_overlaps_file,-f>
    Output file containing information about the overlaps that could not be resolved.

B<--overlap_cutoff,-c>
    The maximum allowed overlap (in nucleotides). Default = 60 nucleotides.

B<--input_sequence_class,-s>
    Optional. Default = assembly.  The class of the parent sequence type

B<--log,-l>
    Logfile.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

 Finds overlaps between genes and between genes and RNAs. Will try to resolve the overlap:

 1. If two genes overlap more than the cutoff
   a. And one gene does not have evidence (as determined by input evidence bsml),
      the gene with no evidence is removed.
   b. If both genes have evidence, a line is printed to the flagged_overlaps_file
   c. If both genes have no evidence, a line is printed to the flagged_overlaps file
 2. If a gene overlaps an RNA (with any overlap)
   a. Remove gene if it has no evidence
   b. Print to flagged_overlaps_file if gene has evidence
 
=head1  INPUT

  --input_file:
     Should be a gene describing bsml document. Will parse out Feature[@class="gene"] elements
     and compare the Interval-loc coordinates with other genes/RNAs with the same parent
     Sequence (based on id).

  --rna_bsml:
     An RNA describing bsml document. Parses Feature[@class="gene"] and looks for overlaps
     with genes

  --evidence_bsml:
     Bsml files describing sequence alignments from blast, ber, and hmm components
       

=head1 OUTPUT

  Will produce a bsml file much like that provided in the --input_file, only with some
  genes removed due to overlaps.

  The flagged_overlaps_file will contain information about unresolved overlaps in the
  following format:

    feature_id1  type1  left1  right1  feature_id2  type2  left2  right2  

  Where type is either gene or RNA

=head1  CONTACT

    Kevin Galens
    kgalens@som.umaryland.edu

=cut
