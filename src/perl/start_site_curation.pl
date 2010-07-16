#!/usr/bin/perl

=head1 NAME

start_site_curation.pl - Looks at the start sites of genes in a bsml file and will
    move the start site based on alignment evidence

=head1 SYNOPSIS

 USAGE: start_site_curation.pl
       --input_file=/path/to/some/glimmer3.bsml
       --output_bsml=/path/to/output.all.bsml
       --changed_features_bsml=/path/to/output.new.bsml
       --evidence=/path/to/ber.bsml.list
       --char_db=/path/to/tchar.db
       --username=user
       --password=pass
       --ber_extension=300
       --percent_identity_cutoff=60
       --p_value_cutoff=1e-30
       --characterized_vote_bonus
       --min_vote_cutoff=3
       --rbs_sliding_window=4
       --rbs_ag_percent_cutoff=75
     [ --log=/path/to/file.log
       --debug=4
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>
    REQUIRED. Path to input gene describing bsml file.

B<--output_file,-o>
    REQUIRED. Path to output bsml file.  All features (those that have not changed
    and those features which have) are included in this file.

B<--changed_features_bsml,-O>
    OPTIONAL. Path to output bsml file which will contain only those features that
    have changed.  If not specified will not print to this file.

B<--evidence,-e>
    REQUIRED. BER evidence list for genes (bsml).

B<--char_db,-c>
    OPTIONAL. Stored hash of characterized proteins. If not specified, will not score
    characterized proteins more

B<--ber_extension,-b>
    OPTIONAL. Default = 300.  The number of nucleotides used in the ber extension

B<--percent_identity_cutoff,-p>
    OPTIONAL. Default = 60.  The percent identity of BER hits.  Filter out
    all those hits with a lower percent_identity.

B<--p_value_cutoff,-P>
    OPTIONAL. Default = 1e-30.  The p_value cutoff for BER alignments.

B<--characterized_vote_bonus,-C>
    OPTIONAL. Default = 4.  The bonus number of votes a start site will receive when 
    confirmed by a characterized protein. (1 non-characterized match = 1 vote)

B<--min_vote_cutoff,-m>
    OPTIONAL. Default = 2.  The number of votes a start_site should receive before being
    considered. Basically, gives a bonus to current start sites.

B<--rbs_sliding_window_size,-r>
    OPTIONAL. Default = 6.  The window size over which percent AG is calculated when
    searching for ribosomal binding site.

B<--rbs_ag_percent_cutoff,-a>
    OPTIONAL. Default = 75.  The percent presence of As and Gs in the sliding window. If
    percent AG is higher than this cutoff, the sequence is considered to have a RBS.

B<--username>
    Database username for querying for characterized status of proteins

B<--password>
    Database password for querying for characterized status of proteins

B<--log,-l>
    Logfile.

B<--debug,-D>
    Debug level

B<--help,-h>
    Print this message

=head1  DESCRIPTION

    A voting based start site curation algorithm.  The algorithm uses evidence generated
    by BER alignments, as well as taking into consideration ribosomal binding sites
    (shine delgarno consensus sequence) as well as start site probabilities.  Reads gene
    information from gene describing BSML input (such as that output from promote_gene_prediction
    or glimmer3) and also evidence (a list of BSML output files from BER alignment run).  
username
    A start site will be changed according to the following algorithm:

    1. BER Evidence
       - All alignments not passing the percent_identity and p_value cutoffs are eliminated
       - Any alignment which begins at a start site within the query sequence gets a vote
       - A bonus is added if the subject sequence is a characterized protein (characterized_
         vote_bonus).
       - The start site with the highest number of votes is chosen as new start if the
         number of votes is greater than or equal to the min_vote_cutoff.  If there is a 
         tie (i.e. two start sites with the same vote count):
    2. Ribosomal Binding Site
       - The sequence upstream of the start site is searched (20 nucleotides upstream
         to 5 nucleotides upstream [coords -20 -> -5 relative to start site] ).
       - Lengths of sequences of size sliding_window_size are measured for AG content. If
         the AG percentage is greater than or equal to the percent_ag_cutoff, the start
         site is considered to have an RBS and is given one vote.
    3. Start Site probabilities
       - The occurence of all the start sites is counted and the relative frequencies are
         calculated.  The relative frequency is added to the vote total.  This will never be
         greater than one vote, and therefore will only be useful in deciding ties at this 
         point.
    

=head1  INPUT

    Input is a gene describing bsml document.  It should contain a Sequence element 
    representing the genomic sequence with Feature children representing the genes.
    Each gene should have a Feature element of class polypeptide, CDS, gene, transcript,
    and exon and their relationships should be described with Feature-groups.  The Sequence
    element must also have a Seq-data-import linking to fasta.  Attributes of the Sequence
    element that are exptected are:
      Sequence[@molecule]
      Sequence[@class]
      Sequence[@id]

    Ids for the Feature elements are assumed to be generated by Ergatis::IdGenerator and
    in the format proj.feature_type.12345.1 . 

    The BER evidence input is a list of BER bsml files which are generated using the 
    ergatis component.  The BSML is formatted as our standard usage for alignment
    analyses.  

=head1 OUTPUT
    
    The output is a gene describing bsml document.  The version number on the genes and
    related features where start site changes have occurred are incremented to indicate 
    a new version of the gene.  

    The script can print two sets of bsml files.  The first set includes the all genes
    and is representative of the working models for the sequence.  The second set is a 
    collection of bsml files which only list the genes which have changed.  The second
    set bsml files (only changed features) is optional and will not be printed if the 
    --changed_features_bsml option is not used.

=head1  CONTACT

    Kevin Galens
    kgalens@som.umaryland.edu

=cut

use strict;
use warnings;
use Ergatis::Logger;
use Data::Dumper;
use BSML::BsmlAlignmentFilter;
use BSML::BsmlBuilder;
use File::Find;
use File::Basename;
use File::OpenFile qw(open_file);
use UnirefClusters::Database;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use MLDBM 'DB_File';
use DB_File;

## required for proper operation on NFS
##  see SF.net bug 2142533 - https://sourceforge.net/tracker2/?func=detail&aid=2142533&group_id=148765&atid=772583
$File::Find::dont_use_nlink = 1;

############# PARAMETERS ########################
## See perldoc for explanation of parameters.
## Defaults are set here.
my $ber_extension = 300;
my $percent_identity_cutoff = 60;
my $p_value_cutoff = 1e-30;
my $characterized_vote_bonus = 4;
my $min_vote_cutoff = 2;

##rbs finding params
my $sliding_window_size = 6;
my $ag_percent_cutoff = 75;

## name of the analysis (for bsml printing)
my $analysis_name = 'start_site_curation';
################################################# 

########### GLOBALS ##################
my %options;
my $input_file;
my $output_file;
my $changed_features_bsml;
my @evidence_files;
my %evidence;
my %feature_relationships;
my %characterized;
my %changed_start_sites;
my $use_db = 0;
my $ur_db;
#######################################

my $results = GetOptions (\%options,
                          'input_file|i=s',
                          'output_bsml|o=s',
                          'changed_features_bsml|O=s',
                          'evidence|e=s',
                          'ber_extension|b=s',
                          'char_db|c=s',
                          'percent_identity_cutoff|p=s',
                          'p_value_cutoff|P=s',
                          'characterized_vote_bonus|C=s',
                          'min_vote_cutoff|m=s',
                          'rbs_sliding_window|r=s',
                          'rbs_ag_percent_cutoff|C=s',
                          'log|l=s',
                          'help|h=s',
                          'debug|D=s',
                          );

## Setup the logger
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE' =>$logfile,
                                 'LOG_LEVEL'=>$options{'debug'});
$logger = Ergatis::Logger::get_logger();

## Check the options
&check_options(\%options);

%evidence = &parse_evidence_files( \@evidence_files );

## Get all the sequences (and features)
my $sequences = &parse_input_bsml( $input_file );

SEQUENCE:
foreach my $sequence_id ( keys %{$sequences} ) {

    my $genes = $sequences->{$sequence_id}->{'genes'};
    $logger->info( "Analyzing ".scalar( @{$genes} )." genes on sequence $sequence_id");

    ## Get the start site frequencies.
    my $start_site_relative_frequencies = &get_start_site_relative_frequencies( $genes );

  GENE:
    foreach my $gene ( sort { $a->{'endpos'} <=> $b->{'endpos'} } @{$genes} ) {
        
        my $gene_id = $gene->{'id'};

        my $evidence = &get_evidence($gene_id);

        #Skip if we have no evidence
        if( @{$evidence} == 0 ) {
            next GENE;
            $logger->info("Skipped $gene_id because it had no evidence");
        } 

        # Grab the sequence
        my $gene_sequence = &get_sequence_interval( $gene->{'parent_seq'},
                                                    $gene->{'startpos'}, 
                                                    $gene->{'endpos'}, );

        # Reverse complement those on the reverse strand
        $gene_sequence = &reverse_complement( $gene_sequence ) if( $gene->{'complement'} );

        # Array ref of in frame start sites and initialization of all in frame start site
        # vote counts to zero.
        my $starts = 
            &find_in_frame_start_sites( $gene_sequence, $ber_extension%3 ); ## The ber extension represents
                                                                            ## the current start site and 
                                                                            ## therefore which frame the
                                                                            ## the gene starts in

        my $votes = {};
        map { $votes->{$_} = 0 } @{$starts};

        # This does the voting for the evidence hits and also a bonus for characterized
        # matches
        foreach my $evidences ( @{$evidence} ) {
            foreach my $ev ( @{$evidences} ) {

                my $is_characterized = &is_characterized( $ev->{'compseq'} );

                foreach my $start_site ( @{$starts} ) {
                    
                    if( $ev->{'left'} == $start_site + $ber_extension ) {
                        $votes->{$ev->{'left'}}++;
                        $votes->{$ev->{'left'}} += $characterized_vote_bonus
                            if( $is_characterized );
                    }

                }
            }
        }

        # Array ref of top scoring start sites
        my $top_start_sites = &get_top_start_sites( $votes );

        # Deal with ties (i.e. two start sites got the same number of evidence votes).
        my $new_start_site;
        my $ret_string;
        if( @{$top_start_sites} > 1 ) {
            ($new_start_site, $ret_string) = &tie_breaker( $top_start_sites, $gene,
                                                           $start_site_relative_frequencies);
        } else {
            $new_start_site = $top_start_sites->[0];
            $ret_string = "Chose start site based on evidence";
        }


        # If new_start_site is equal to the extension, this means that it's the original
        # start site and we don't have to do anything.
        if( $new_start_site != $ber_extension ) { 
            $changed_start_sites{$sequence_id}->{$gene_id} = $new_start_site;
            $logger->info("Changed start site for $gene_id to $new_start_site");
            $logger->info("$ret_string");
        } 

    } ## END GENES

    if( !exists( $changed_start_sites{$sequence_id} ) ) {
        $changed_start_sites{$sequence_id} = {};
    }

} ## END SEQUENCES


&write_output( $input_file, $output_file, $sequences, $changed_features_bsml );

############################## SUBS ############################
sub is_characterized {
    my ($compseq) = @_;
    my $retval = 0;

    #Are we using the lookup file or the database
    if( $use_db ) {
        my $cluster_id = $ur_db->get_cluster_id_by_acc( $compseq );
        $retval = $ur_db->cluster_is_trusted( $cluster_id ) if( defined( $cluster_id ) );
    } else {
        my $portion = "$1|$2" if( $compseq =~ /^([^\_]+)\_([^\_]+\_?[^\_]+)\_/ );
        $retval = defined( $portion ) ? exists( $characterized{$portion} ) : exists( $characterized{$compseq} );
    }

    return $retval;
}

sub write_output {
    my ($in_file, $out_file, $sequences, $changed_features_out ) = @_;

    &write_bsml( $in_file, $out_file, $sequences );
    print STDERR "$out_file\n";

    if( $changed_features_out ) {
        &write_bsml( $in_file, $changed_features_out, $sequences, 1 );
        print STDERR "$changed_features_out\n";
    }
}

sub write_bsml {
    my ($in_file, $out_file, $sequences, $only_new_features) = @_;

    my $doc = new BSML::BsmlBuilder;

  SEQUENCE:
    foreach my $seq_id ( keys %{$sequences} ) {
        #Add sequence
        my $seq = $doc->createAndAddSequence( $seq_id, $seq_id, $sequences->{$seq_id}->{'length'},
                                              $sequences->{$seq_id}->{'molecule'}, 
                                              $sequences->{$seq_id}->{'class'} );

        #Add seq-data-import
        my $sdi = $doc->createAndAddSeqDataImport( $seq, 'fasta', $sequences->{$seq_id}->{'source'},
                                                   '', $sequences->{$seq_id}->{'identifier'} );

        #Add link to the analysis
        my $link = $doc->createAndAddLink( $seq, 'analysis',  '#'.$analysis_name.'_analysis',
                                           'input_of' );

        #Add the feature table
        my $feature_table = $doc->createAndAddFeatureTable( $seq );

        #tmp
        my $count = 0;

        #Cycle through the genes\
      GENE:
        foreach my $gene ( @{$sequences->{$seq_id}->{'genes'}} ) {

            if( !exists( $changed_start_sites{ $seq_id } ) ) {
                $logger->logdie("$seq_id doesn't exist in changed_start_sites hash");
            }
            
            #Don't add the gene if we are only printing genes with new start sites
            #and the start site wasn't changed.
            next GENE if( $only_new_features && !exists( $changed_start_sites{$seq_id}->{$gene->{'id'}} ) );

            #Otherwise, go ahead and add the gene (and related features) to the document
            &add_gene_to_doc( $doc, $seq, $feature_table, $gene );
            
        }
    }

    #Add analysis
    my $source_dir = dirname( $in_file );
    my $analysis = $doc->createAndAddAnalysis( 'id'         => $analysis_name."_analysis",
                                               'sourcename' => $source_dir,
                                               'algorithm'  => $analysis_name,
                                               'program'    => $analysis_name
                                               );

    my $oh = open_file( $out_file, 'out' );
    eval {
        $doc->write( $oh );
    };
    if( $@ ) {
        $logger->logdie("$@");
    }
    close( $oh );
}

sub add_gene_to_doc {
    my ($doc, $seq, $ft, $gene) = @_;
    my $new_ss_flag = 0;
    my $gene_id = $gene->{'id'};
    my $interval_loc = {
        'startpos' => $gene->{'startpos'},
        'endpos' => $gene->{'endpos'},
        'complement' => $gene->{'complement'},
    };

    my $new_ss =  $changed_start_sites{ $gene->{'parent_seq'} }->{ $gene->{'id'} }
    if( exists( $changed_start_sites{ $gene->{'parent_seq'} }->{ $gene->{'id'} } ) );

    #was the start site changed?
    if( $new_ss ) { 
        $new_ss_flag = 1;
        $gene_id = &get_next_id_version( $gene->{'id'} );
        
        #Stored are start sites from the ber hits. These are not genomic coords
        my $relative_start_site = $changed_start_sites{$gene->{'parent_seq'}}->{$gene->{'id'}};
        my $delta_start_site = $relative_start_site - $ber_extension;

        #If it's on the reverse strand, must substract
        if( $interval_loc->{'complement'} ) {
            $interval_loc->{'endpos'} = $interval_loc->{'endpos'} - $delta_start_site;
        } else {
            $interval_loc->{'startpos'} = $interval_loc->{'startpos'} + $delta_start_site;
        }

    }

    #Add the gene to the feature table
    my $bsml_gene = $doc->createAndAddFeature( $ft, $gene_id, $gene_id, 'gene' );
    $doc->createAndAddIntervalLoc( $bsml_gene, $interval_loc->{'startpos'}, $interval_loc->{'endpos'},
                                   $interval_loc->{'complement'} );
    $doc->createAndAddLink( $bsml_gene, 'analysis', '#'.$analysis_name."_analysis", 'computed_by' );

    #Create the feature group
    my $fg = $doc->createAndAddFeatureGroup( $seq, '', $gene_id );

    my $fgm = $doc->createAndAddFeatureGroupMember( $fg, $gene_id, 'gene' );

    #Add the other features
    foreach my $class ( qw(transcript exon polypeptide CDS) ) {
        my $feat_id = &get_feature_relationship( $gene->{'id'}, $class );
        $logger->logdie("Could not find relationship to $gene->{'id'} for $class")
            unless( $feat_id );

        $feat_id = &get_next_id_version( $feat_id ) if( $new_ss_flag );

        my $feat = $doc->createAndAddFeature( $ft, $feat_id, $feat_id, $class );
        $doc->createAndAddIntervalLoc( $feat, $interval_loc->{'startpos'}, $interval_loc->{'endpos'},
                                       $interval_loc->{'complement'} );
        $doc->createAndAddLink( $feat, 'analysis', '#'.$analysis_name."_analysis", 'computed_by' );
        $doc->createAndAddFeatureGroupMember( $fg, $feat_id, $class );
    }
    
}

#expects ids in the format db.feature_type.number.version
#and will return the same id with the version number incremented
sub get_next_id_version {
    my ($id) = @_;
    my $retval;

    if( $id =~ /(.*\.)(\d+)$/ ) {
        my $new_version = $2 + 1;
        $retval = $1.$new_version;
    } else {
        $logger->logdie("Could not update the version of feature id $id");
    }

    return $retval;
}

sub tie_breaker {
    my ($top_start_sites, $gene, $start_site_relative_frequencies ) = @_;
    my $votes = {};
    my $ret_string = "";

  START_SITE:
    foreach my $start_site (@{$top_start_sites}) {

        # transform $start_site (relative to alignment) into genomic coords
        my $region = &get_upstream_region( $gene, $start_site );

        #Pass in upstream region from -20 -> -5 (relative to start codon position 0 )
        if( &has_rbs( substr( $region, 0 , 15 ) ) ) {
            $votes->{$start_site} = 1;
            $ret_string = "Chose start site based on evidence, presence of rbs site";
        }

        #Grab the start codon
        my $start_codon = substr( $region, length($region) - 3, length($region) );
        my $freq = 0;
        $freq = $start_site_relative_frequencies->{$start_codon} 
        if( exists( $start_site_relative_frequencies->{$start_codon} ) );
        $votes->{$start_site} += $freq;
        $ret_string .= " and start site [$start_site] has a frequency of $freq";

    }

    my $new_top_start_sites = &get_top_start_sites( $votes );
    my $new_start_site;
    
    #if we still have a tie
    if( @{$new_top_start_sites} > 1 ) {
        
        #If the original start site is still in contention use that
        $new_start_site = $ber_extension if( grep( /$ber_extension/, @{$new_top_start_sites} ));
        $ret_string = "tied after rbs and start site frequency, using original start site";
        
        #If not, use the longest prediction
        unless( $new_start_site ) {
            ($new_start_site) = sort { $a <=> $b } @{$new_top_start_sites};
            $ret_string = "tied after rbs and start site frequency consideration, original start site ".
                "didn't get enough votes, so using longest prediction";
        }
                                  

    } else {
        $new_start_site = $new_top_start_sites->[0];
    }

    return ($new_start_site, $ret_string);    
}

sub get_upstream_region {
    my ($gene, $start_site) = @_;
    
    my ($genomic_left, $genomic_right) = ($gene->{'startpos'}, $gene->{'endpos'});
    my $genomic_start = $genomic_left - $ber_extension + $start_site;

    my ($rbs_start, $rbs_stop) = ($genomic_start - 20, 
                                  $genomic_start + 3 );
    if( $gene->{'complement'} ) {
        $genomic_start = $genomic_right - ( $start_site - $ber_extension );
        ($rbs_start, $rbs_stop) = ($genomic_start + 20,
                                   $genomic_start - 3 );
    }
    
    #Get the upstream region of the gene (include the start codon).
    my $region = &get_sequence_interval( $gene->{'parent_seq'},
                                         $rbs_start, 
                                         $rbs_stop );

    return $region;
}

#Takes sequence 15 na long 
sub has_rbs {
    my ($sequence) = @_;
    $logger->logdie("sequence passed to has_rbs must be 23 nas long") 
        unless( length($sequence) == 15 );

    my $retval = 0;
    for( my $i = 0; $i < length( $sequence ) - $sliding_window_size; $i++ ) {
        my $window = substr( $sequence, $i, $i + $sliding_window_size );
        my $count = 0;
        $count++ while( $window =~ /[AG]/ig );
        my $percent = ($count/length($window)) * 100;
        if( $percent >= $ag_percent_cutoff ) {
            $retval = 1;
            last;
        }
    }

    return $retval;
    
    
}
sub get_start_site_relative_frequencies {
    my ($genes) = @_;
    my $ss_frequencies = {}; #Start site frequencies
    my $total = 0;

    foreach my $gene ( @{$genes} ) {
        my ($ss_start, $ss_end) = ($gene->{'startpos'}, $gene->{'startpos'}+3);
        if( $gene->{'complement'} ) {
            ($ss_start, $ss_end) = ($gene->{'endpos'}, $gene->{'endpos'}-3);
        }
        
        my $region = &get_sequence_interval( $gene->{'parent_seq'},
                                             $ss_start,
                                             $ss_end );

        $ss_frequencies->{$region}++;
        $total++;
    }
    
    my %ss_rel_freqs = map{ $_ => ($ss_frequencies->{$_})/$total } keys %{$ss_frequencies};
    return \%ss_rel_freqs;

}

sub get_top_start_sites {
    my $votes = shift;
    my @values;
    
    my ($max) = sort{ $b <=> $a } values( %{$votes} );
    
    #If the start site with the most votes doesn't beat the minimum
    #vote cutoff, then none of them will.  So just push on the extension
    #which indicates the original start site (i.e. the alignments weren't good enough
    #to determine new start site, so default to gene caller).
    if( $max < $min_vote_cutoff ) {
        push(@values, $ber_extension );
    } else {

        map { push(@values, $_) if( $votes->{$_} == $max ); } keys %{$votes};

    }

    return \@values;
}

sub find_in_frame_start_sites {
    my ($sequence, $frame) = @_;
    my $retval = [];

    if( $frame < 0 || $frame > 2 ) {
        $logger->logdie("Frame passed into find_in_frame_start_sites should be a value ".
                        "between 0 and 2.  [current value: $frame]");
    }

    #Data structure of start Codons
    my $codons ={
        'ATG' => 'M',
        'GTG' => 'V',
        'TTG' => 'L' 
        };
    
    foreach my $codon ( keys %{$codons} ) {

        my $pos = -1; #Make sure we search from index 0
        while(($pos = index($sequence, $codon, $pos+1)) > -1) {

            #only those in frame
            next unless( $pos%3 == $frame );
            push(@{$retval}, $pos);
        }
    }
    
    return $retval;
    
}
sub reverse_complement {
    my ($seq) = @_;
    $seq =~ tr/ACGT/TGCA/;
    return reverse($seq);
}

sub parse_input_bsml {
    my ($bsml_file) = @_;
    my @genes = ();
    my %sequences;

    $logger->logdie("File does not exist [$bsml_file]") unless( -e $bsml_file );

    my $twig = new XML::Twig( 'twig_roots' => {
        'Sequence' => sub {
            my ($xtwig, $seq_elem) = @_;
            my $seq_id = $seq_elem->att('id');
            my $seq_class = $seq_elem->att('class');
            my $seq_data_import = $seq_elem->first_child('Seq-data-import');

            foreach my $feat_elem ( $seq_elem->find_nodes('//Feature[@class="gene"]') ) {

                my $int_loc = $feat_elem->first_child('Interval-loc');
                my ($startpos,$endpos,$complement) = ($int_loc->att('startpos'),
                                                      $int_loc->att('endpos'),
                                                      $int_loc->att('complement'),);

                my $id = $feat_elem->att('id');
                my $feat = {
                    'id' => $id,
                    'startpos' => $startpos,
                    'endpos' => $endpos,
                    'complement' => $complement,
                    'parent_seq' => $seq_id,
                };

                push(@genes, $feat);
            }
            

            #Create the sequence data structure if the sequence had features.
            #Otherwise, it's just a CDS or polypeptide Sequence.
            if( $seq_class ne 'CDS' && $seq_class ne 'polypeptide' ) {
                my $sequence = { 
                    'source' => $seq_data_import->att('source'),
                    'identifier' => $seq_data_import->att('identifier'),
                    'genes' => \@genes,
                    'class' => $seq_elem->att('class'),
                    'molecule' => $seq_elem->att('molecule') 
                    };
                $sequences{$seq_id} = $sequence;
                
            }

        },
        'Feature-group' => sub {
            my ($xtwig, $fg_elem) = @_;
            my $added = &store_feature_relationship( $fg_elem->children('Feature-group-member') );
            $logger->logdie("Didn't store any feature relationships")
                unless( $added > 0 );
        },
    });

    my $fh;
    eval {
        $fh = open_file( $bsml_file, 'in' );
    };
    if( $@ ) {
        $logger->logdie("Could not open $bsml_file [$!]");
    }

    $twig->parse($fh);
    close($fh);

    return \%sequences;
            
}

#Name: get_evidence
#Desc: will retrieve parsed evidence stored in the global var %evidence.
#      This evidence was parsed in subroutine 'parse_evidence' and the
#      hash is keyed using whatever id was stored in the 
#      Seq-pair-alignment[@refseq] attribute. (ex. CDS for BER and polypeptide for HMM).
#Args: A feature_id
#Rets: Reference to an array of parsed evidence objects (as stored by sub parse_evidence)
#      and undef if there is no evidence.
sub get_evidence {
    my ($feature_id) = @_;
    my $evidence = [];

  CLASSES:
    foreach my $class ( qw/polypeptide CDS gene transcript exon/ ) {
        my $class_id = &get_feature_relationship( $feature_id, $class );
        $logger->logdie("Could not find id for $class when using id $feature_id")
            unless( defined( $class_id ) );
        if( exists( $evidence{$class_id} ) ) {
            push(@{$evidence}, $evidence{$class_id});
        }
    }
    return $evidence;
}

#Name: parse_evidence
#Desc: will parse a list of evidence files into a data structure
#Args: reference to an array of evidence file names
#Rets: a hash in the following format:
sub parse_evidence_files {
    my $evidence_files = shift;

    # return value
    my %evidence = ();

    my $total = scalar(@{$evidence_files});
    my $count = 0;

    #Parse each of the evidence files
    foreach my $evidence_file ( @{$evidence_files} ) {
        print "\r$count/$total";
        $count++;

        my $af = new BSML::BsmlAlignmentFilter( { 
            'file' => $evidence_file
            });

        $af->add_filter( 'p_value', "-".$p_value_cutoff );
        $af->add_filter( 'percent_identity', "+".$percent_identity_cutoff );
        $af->add_filter( 'trusted_cutoff', 1 );
        
        $af->parse_alignment_file();

        foreach my $refseq ( $af->refseqs() ) {
            $evidence{$refseq} = $af->get_alignment_interval( $refseq );
        }
    }

    return %evidence;
    
}

#interbase numbering
sub get_sequence_interval {
    my ($seq_id, $left, $right) = @_;

    my $reverse_complement = 0;
    if( $left > $right ) {
        ($right, $left) = ($left,$right);
        $reverse_complement = 1;
    } elsif( $left == $right ) {
        return "";
    }

    my $seq = "";

    #Have we already parsed the sequence?
    if( exists( $sequences->{$seq_id}->{'sequence'} ) ) {
        $seq = $sequences->{$seq_id}->{'sequence'};
    } else {

        #Get the fasta file
        my $fasta_file;
        if( exists( $sequences->{$seq_id} ) ) {
            $fasta_file = $sequences->{$seq_id}->{'source'};
        } else {
            $logger->logdie("Could not find sequence information for sequence $seq_id".
                            " in subroutine get_sequence_interval");
        }
        my $fh = open_file("$fasta_file", 'in' );

        my $parse = 0;
        while( <$fh> ) {

            if( $parse ) {
                next if( /^\s+$/ );
                chomp;
                $seq .= $_;            
            }

            #Find the correct sequence (in case there are multiple)
            if( /^>(\S+)/ ) {
                if( $1 eq $sequences->{$seq_id}->{'identifier'} ) {
                    $parse = 1;
                }
            }
        }

        $sequences->{$seq_id}->{'sequence'} = $seq;
        
    }

    #Store the length of the sequence
    $sequences->{$seq_id}->{'length'} = length($seq);
    
    if( $right > length( $seq ) ) {
        $right = length( $seq );
    }

    my $retval = substr($seq, $left, $right - $left);
    $retval = &reverse_complement( $retval ) if( $reverse_complement );
    return $retval;

}

sub store_feature_relationship {
    my (@feature_group_members) = @_;

    my $lookup = {};
    foreach my $fgm ( @feature_group_members ) {
        my $featref = $fgm->att('featref');
        my $feature_type = $fgm->att('feature-type');

        $lookup->{$feature_type} = $featref;
        $feature_relationships{$featref} = $lookup;
    }

    return scalar(keys %{$lookup});
}

sub get_feature_relationship {
    my ($id, $class) = @_;
    my $return_id = undef;

    if( exists( $feature_relationships{$id} ) ) {
        $return_id = $feature_relationships{$id}->{$class};
    }

    return $return_id;

}


## Checks the input options.  See pod for specifics.
sub check_options {
   my $opts = shift;

   if( $opts->{'help'} ) {
       &_pod;
   }

   if( $opts->{'input_file'} ) {
       $input_file = $opts->{'input_file'};
   } else {
       $logger->logdie("Option input_file is required");
   }

   if( $opts->{'output_bsml'} ) {
       $output_file = $opts->{'output_bsml'};
   } else {
       $logger->logdie("Option output_bsml is required");
   }

   if( $opts->{'changed_features_bsml'} ) {
       $changed_features_bsml = $opts->{'changed_features_bsml'};
   } else {
       undef $changed_features_bsml;
   }

   if( $opts->{'evidence'} ) {
       my @lists = split(/[,\s]+/, $opts->{'evidence'} );

       foreach my $list ( @lists ) {
           my $fh = open_file( $list, 'in' );
           chomp( my @files = <$fh> );
           push(@evidence_files, @files);
       }
   } else {
       $logger->logdie("Option evidence is required");
   }

   if( $opts->{'char_db'} ) {
       tie(%characterized, 'MLDBM', $opts->{'char_db'}, O_RDONLY )
           or $logger->logdie("Could not tie $opts->{'char_db'} to hash");
   }

   if( $opts->{'username'} && $opts->{'password'} ) {
       $use_db = 1;
       $ur_db = new UnirefClusters::Database( "username" => $opts->{'username'},
                                              "password" => $opts->{'password'} );
   }

   if( $opts->{'ber_extension'} ) {
       $ber_extension = $opts->{'ber_extension'};
   }

   if( $opts->{'percent_identity_cutoff'} ) {
       $percent_identity_cutoff = $opts->{'percent_identity_cutoff'};
   }

   if( $opts->{'p_value_cutoff'} ) {
       $p_value_cutoff = $opts->{'p_value_cutoff'};
   }

   if( $opts->{'characterized_vote_bonus'} ) {
       $characterized_vote_bonus = $opts->{'characterized_vote_bonus'};
   }

   if( $opts->{'min_vote_cutoff'} ) {
       $min_vote_cutoff = $opts->{'min_vote_cutoff'};
   }

   if( $opts->{'rbs_sliding_window'} ) {
       $sliding_window_size = $opts->{'rbs_sliding_window'};
   }

   if( $opts->{'rbs_ag_percent_cutoff'} ) {
       $ag_percent_cutoff = $opts->{'rbs_ag_percent_cutoff'};
   }

}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => *STDERR} );
}
