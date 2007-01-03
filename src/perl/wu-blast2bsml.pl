#!/usr/local/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

wu-blast2bsml.pl  - convert raw blast output to BSML documents

=head1 SYNOPSIS

USAGE:  wu-blast2bsml.pl -i blastp.raw -o blastp.bsml -q /path/to/query_file.fsa -l /path/to/logfile.log

=head1 OPTIONS

B<--input,-i>
    The input raw blast output from the wu-blast suite.

B<--output,-o>
    The output BSML file name.

B<--query_file_path,-q>
    The full path to the query FASTA file (used to populate refxref)
    [When used in workflow should be passed ITER_FILE_PATH]

B<--gzip_output,-g>
    Compress the bsml output.  If there is not .gz extension on the the output file
    one will be added (optional).
    
B<--max_hsp_count,-m> Maximum number of HSPs stored per alignment (optional)

B<--pvalue,-p> Maximum P-value at and above which to exclude HSPs (optional)

B<--class,-c> The ref/comp sequence type.  Default is 'assembly'

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1   DESCRIPTION

blastbtab2bsml.pl takes raw output from WU-BLAST and translates it to BSML.
In the case of an empty results set, a sequence stub and analysis link will
still be created for the query sequence. 

NOTE:  

Calling the script name with NO flags/options or --help will display the syntax requirement.

=head1  CONTACT
    Brett Whitty
    bwhitty@tigr.org

=cut

use strict;
use Bio::SearchIO;
use English;
use File::Basename;
use File::Path;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
BEGIN {
use Ergatis::Logger;
use BSML::BsmlRepository;
use BSML::BsmlBuilder;
use BSML::BsmlParserTwig;
}
my $defline;
my %options = ();
my $results = GetOptions (\%options, 
                          'input|i=s',
                          'output|o=s',
                          'query_file_path|q=s',
                          'gzip_output|g=s',
                          'log|l=s',
                          'debug|d=s',
                          'analysis_id|a=s',
                          'bsml_dir|d=s', ## deprecated.  keeping for backward compat (for now)
                          'max_hsp_count|m=s',
                          'pvalue|p=s', 
                          'debug=s',
                          'class|c=s',
                          'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

## display documentation
if( $options{'help'} || scalar keys(%options) == 0 ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

if($options{'pvalue'} eq ""){
    $options{'pvalue'} = 10;
}

my $ref_molecule = {
                        'wu-blastn' => 'na',
                        'wu-blastp' => 'aa',
                        'wu-blastx' => 'na',
                        'wu-tblastn' => 'aa',
                        'wu-tblastx' => 'na',
                    };
my $comp_molecule = {
                        'wu-blastn' => 'na',
                        'wu-blastp' => 'aa',
                        'wu-blastx' => 'aa',
                        'wu-tblastn' => 'na',
                        'wu-tblastx' => 'na',
                    };


## make sure everything passed was peachy
&check_parameters(\%options);

my $class;
if (!defined($options{'class'})){
    $logger->logdie("class was not defined");
} else {
    $class = $options{'class'};
}

##determine wu-blast program
my $blast_program;
if ($options{'analysis_id'} =~ /^([^_]+)_analysis/) {
    $blast_program = $1;
} else {
    $logger->logdie("analysis_id '$options{analysis_id}' has unexpected structure");
}

##check that we've got defined molecule types for the analysis program
unless($ref_molecule->{$blast_program} && $comp_molecule->{$blast_program}) {
    $logger->logdie("Molecule types for '$options{analysis_id}' are undefined");
}

## get a filehandle on the input
open(my $ifh, "<$options{input}") || $logger->logdie("can't read the input sequence: $!");

my $in = new Bio::SearchIO(-format => 'blast', 
                           -fh     => $ifh);

## open the output file:
##open (my $ofh, ">$options{output}") || $logger->logdie("can't create output file for BLAST report: $!");
my $doc = new BSML::BsmlBuilder();
$doc->makeCurrentDocument();

my $hsplookup = {};
my @hsp_ref_array;

# parse each blast record:
while( my $result = $in->next_result ) {
    my $hsp_counter = 0;
    # parse each hit per record.
    while( my $hit = $result->next_hit ) {
        # a hit consists of one or more HSPs
        while( my $hsp = $hit->next_hsp ) {
            my @x;
            $x[0] = $result->query_name();
            # date
            $x[2] = $result->query_length();
            $x[3] = $hsp->algorithm();
            
            ## database name will get parsed with whitespace if its path is long
            $x[4] = $result->database_name();
            $x[4] =~ s/\s//g;
            
            $x[5] = $hit->name();
            $x[6] = $hsp->start('query');
            $x[7] = $hsp->end('query');
            my $queryStrand = $hsp->strand('query');
            if ($queryStrand == -1) {
                ($x[6], $x[7]) = ($x[7], $x[6]);
            }

            $x[8] = $hsp->start('hit');
            $x[9] = $hsp->end('hit');
            my $hitStrand = $hsp->strand('hit');
            if ($hitStrand == -1) {
                ($x[8], $x[9]) = ($x[9], $x[8]);
            }

            $x[10] = sprintf ("%.1f", $hsp->percent_identity());   

            my $similarity = $hsp->frac_conserved('total') * 100; 
            $x[11] = sprintf("%.1f", $similarity);
            $x[12] = $hsp->score();
            $x[13] = $hsp->bits();
            $x[15] = $hit->description();
            $x[16] = ( ($hsp->query->frame + 1) * $hsp->query->strand); #blast frame (1, 2, 3, -1, -2, -3).

            my $strandDescript = "null";
            if ($queryStrand == 1) {
                $strandDescript = "Plus";
            } elsif ($queryStrand == -1) {
                $strandDescript = "Minus";
            }

            $x[17] = $strandDescript;
            $x[18] = $hit->length();
            $x[19] = $hsp->evalue();
            $x[20] = $hsp->pvalue();

            if(($x[20] < $options{'pvalue'}) && ($x[0] ne "") && ($x[5] ne "")){
                ## pvalue is less than cutoff parameter
                ##      so process btab line
                $hsp_counter++;
                
                if(($x[20] <= $options{'pvalue'}) && ($x[0] ne "") && ($x[5] ne "")){
                    if(!(exists $hsplookup->{$x[0]}->{$x[5]})){
                        $hsplookup->{$x[0]}->{$x[5]} = [];
                    }

                    push @{$hsplookup->{$x[0]}->{$x[5]}},{'pvalue'=>$x[20], 'line'=>\@x};

                } else {
                    $logger->debug("Skipping $x[0] $x[5] pvalue $x[20] above pvalue cutoff of $options{'pvalue'}") if($logger->is_debug());
                }
            } else {
                ## else we're gonna skip the line
                $logger->debug("Skipping $x[0] $x[5] pvalue $x[20] above pvalue cutoff of $options{'pvalue'}") if($logger->is_debug());
            }

        }
    }
    if ($hsp_counter == 0) {
        my $align = &createAndAddNullResult(
                                            doc             => $doc,
                                            query_name      => $result->query_name(),
                                            query_length    => $result->query_length(),
                                            class           => $class,
                                           );
    }
}

    my $qual_stats = {};
    foreach my $query (keys %{$hsplookup}){
        #   print "Query: $query\n";
        foreach my $subject (keys %{$hsplookup->{$query}}) {
        #   print "Subject: $subject\n";
            my @hsps = sort {$a->{'pvalue'} <=> $b->{'pvalue'}} @{$hsplookup->{$query}->{$subject}};
            my $maxhsp;
            if ($options{'max_hsp_count'} ne "") {
                $maxhsp = ($options{'max_hsp_count'}<scalar(@hsps)) ? $options{'max_hsp_count'} : scalar(@hsps);
            } else {
                $maxhsp = scalar(@hsps);
            }
            my $queryid;
            
                ######## Calculate Coverage info
                #
                #

            ## calculate coverage stuff
            my $new_subject;
            @hsp_ref_array = ();
            for (my $i=0; $i<$maxhsp; $i++) {
                my @btab = @{$hsps[$i]->{'line'}};
                $queryid = $btab[0] if ($btab[0] && (!$queryid));
                for (my $i=0; $i<scalar(@btab); $i++) {
                    if ($btab[$i] eq 'N/A'){
                        $btab[$i] = undef;
                    }
                }
                my $orig_dbmatch_accession = $btab[5];
                #$btab[5] =~ s/[^a-z0-9\_\.\-]/_/gi;

                $new_subject = $btab[5];
                
                my $qfmin = $btab[6];
                my $qfmax = $btab[7];
                my $qstrand = 0;
                my $tfmin = $btab[8];
                my $tfmax = $btab[9];
                my $tstrand = 0;
                    
                ## if query positions are on the reverse strand
                if ($btab[6] > $btab[7]) {
                        $qfmin = $btab[7];
                        $qfmax = $btab[6];
                        $qstrand = 1;
                }
                
                ## if target positions are on the reverse strand
                if ($btab[8] > $btab[9]) {
                        $tfmin = $btab[9];
                        $tfmax = $btab[8];
                        $tstrand = 1;
                }
                
                ## transform the start positions to interbase 
                $qfmin = $qfmin - 1;
                $tfmin = $tfmin - 1;    
                    
                my $hsp_ref = { 
                                'query_protein_id'  => $btab[0],
                                'target_protein_id' => $btab[5],
                                'significance'      => $btab[20],
                                'percent_identity'  => $btab[10],
                                'percent_similarity'=> $btab[11],
                                'query_seqlen'      => $btab[2],
                                'target_seqlen'     => $btab[18],
                                'query_fmin'        => $qfmin,
                                'query_fmax'        => $qfmax,
                                'query_strand'      => $btab[17],
                                'target_fmin'       => $tfmin,
                                'target_fmax'       => $tfmax,
                                'target_strand'     => $tstrand, ## target strand is not captured in btab 
                              };
                push (@hsp_ref_array, $hsp_ref);
            }
            if (!defined($qual_stats->{$queryid})) {
               $qual_stats->{$queryid} = {};
            }
            
            my $coverage_arr_ref = &getAvgBlastPPctCoverage(\@hsp_ref_array);
            my $id_sim_arr_ref = &getAvgBlastPIdSim(\@hsp_ref_array); 
            
            $qual_stats->{$queryid}->{$new_subject} = {
                'percent_coverage_refseq'   =>  sprintf("%.1f",$coverage_arr_ref->[0]),
                'percent_coverage_compseq'  =>  sprintf("%.1f",$coverage_arr_ref->[1]),
                'percent_identity'          =>  sprintf("%.1f",$id_sim_arr_ref->[0]),
                'percent_similarity'        =>  sprintf("%.1f",$id_sim_arr_ref->[1]),
                                                        };
                #
                #
                #################################

                
            for (my $i=0; $i<$maxhsp; $i++) {
                my @btab = @{$hsps[$i]->{'line'}};
                $logger->debug("Storing HSP $btab[0] $btab[5] $btab[20]") if ($logger->is_debug());
                
                $queryid = $btab[0] if ($btab[0] && (!$queryid));
                for (my $i=0;$i<scalar(@btab);$i++){
                    if ($btab[$i] eq 'N/A'){
                        $btab[$i] = undef;
                    }
                }
            
                ## dbmatch_accession needs to be alphanumeric or _-.
                ##  but the original needs to be passed to createAndAddBtabLine so it can
                ##  be recognized and parsed
                my $orig_dbmatch_accession = $btab[5];
                #$btab[5] =~ s/[^a-z0-9\_\.\-]/_/gi;

                my $align = &createAndAddBlastResultLine(
                              doc                => $doc,
                              query_name         => $btab[0],
                              date               => $btab[1],
                              query_length       => $btab[2],
                              blast_program      => $btab[3],
                              search_database    => $btab[4],
                              dbmatch_accession  => $btab[5],
                              start_query        => $btab[6],
                              stop_query         => $btab[7],
                              start_hit          => $btab[8],
                              stop_hit           => $btab[9],
                              percent_identity   => $btab[10],
                              percent_similarity => $btab[11],
                              raw_score          => $btab[12],
                              bit_score          => $btab[13],
                              dbmatch_header     => $btab[15],
                              blast_frame        => $btab[16],
                              qry_strand         => $btab[17],
                              hit_length         => $btab[18],
                              e_value            => $btab[19],
                              p_value            => $btab[20],
                              percent_coverage_refseq    => $qual_stats->{$btab[0]}->{$btab[5]}->{'percent_coverage_refseq'},
                              percent_coverage_compseq   => $qual_stats->{$btab[0]}->{$btab[5]}->{'percent_coverage_compseq'},
                              percent_identity_total     => $qual_stats->{$btab[0]}->{$btab[5]}->{'percent_identity'},
                              percent_similarity_total   => $qual_stats->{$btab[0]}->{$btab[5]}->{'percent_similarity'},
                              class                      => $class,
                              orig_dbmatch_accession     => $orig_dbmatch_accession
                                                    );

                my $seq = $doc->returnBsmlSequenceByIDR($btab[5]);

            }
            if ($queryid) {
                my $seq = $doc->returnBsmlSequenceByIDR($queryid);
            }
        }
    }
    
$doc->createAndAddAnalysis(
                            id => $options{analysis_id},
                            sourcename => $options{'output'},
                          );

$doc->write($options{'output'}, '', $options{'gzip_output'});

exit(0);

##Adds BSML tags for the case where 
##the query sequence returned no hits
sub createAndAddNullResult {
    my %args = @_;
    my $doc = $args{'doc'};
    
    if( !( $doc->returnBsmlSequenceByIDR( "$args{'query_name'}")) ){
        my $seq = $doc->createAndAddSequence(
                            "$args{'query_name'}", 
                            "$args{'query_name'}",
                            $args{'query_length'}, 
                            $ref_molecule->{$blast_program}, 
                            $args{'class'}
                                            );
        $doc->createAndAddSeqDataImport(
                            $seq, 
                            'fasta', 
                            $options{'query_file_path'}, 
                            '', 
                            $args{'query_name'}
                                       );
        $seq->addBsmlLink(
                            'analysis', 
                            '#' . $options{'analysis_id'}, 
                            'input_of'
                         );
        $doc->createAndAddBsmlAttribute( 
                            $seq, 
                            'defline', 
                            $defline
                         );
    }
}

sub createAndAddBlastResultLine {
    
    my %args = @_;
    my $doc = $args{'doc'};

#    my $orig_dbmatch_accession = $args{dbmatch_accession};
#    $args{dbmatch_accession} =~ s/[^a-z0-9\_]/_/gi;

    #determine if the query name and the dbmatch name are a unique pair in the document
    my $alignment_pair_list = BSML::BsmlDoc::BsmlReturnAlignmentLookup( $args{'query_name'}, $args{'dbmatch_accession'} );

    my $alignment_pair = '';
    if( $alignment_pair_list ) {
        $alignment_pair = $alignment_pair_list->[0];
    }

    if( $alignment_pair ) {
        #add a new BsmlSeqPairRun to the alignment pair and return
        my $seq_run = $alignment_pair->returnBsmlSeqPairRunR( $alignment_pair->addBsmlSeqPairRun() );

        $doc->createAndAddBsmlAttribute($seq_run, 'class', 'match_part');

        if( $args{'start_query'} > $args{'stop_query'} ) {
            $seq_run->setattr( 
                                'refpos', 
                                $args{'stop_query'}-1 
                             );
            $seq_run->setattr( 
                                'runlength', 
                                $args{'start_query'} - $args{'stop_query'} + 1 
                             );
            $seq_run->setattr(
                                'refcomplement', 
                                1 
                             );
        } else {
            $seq_run->setattr( 
                                'refpos', 
                                $args{'start_query'}-1
                             );
            $seq_run->setattr(
                                'runlength', 
                                $args{'stop_query'} - $args{'start_query'} + 1
                             );
            $seq_run->setattr( 
                                'refcomplement', 
                                0 
                             );
        }

        #the database sequence is always 5' to 3'
        $seq_run->setattr(
            'comppos', 
            $args{'start_hit'} -1
                         ) if (defined ($args{'start_hit'}));
        $seq_run->setattr(
            'comprunlength', 
            $args{'stop_hit'} - $args{'start_hit'} + 1
                         ) if ((defined ($args{'start_hit'})) and (defined ($args{'stop_hit'})));
        $seq_run->setattr( 
            'compcomplement', 
            0 
                         );
        $seq_run->setattr(
            'runscore', 
            $args{'bit_score'} 
                         ) if (defined ($args{'bit_score'}));
        $seq_run->setattr(
            'runprob', 
            $args{'e_value'} 
                         ) if (defined ($args{'e_value'}));
        $seq_run->addBsmlAttr( 
            'percent_identity', 
            $args{'percent_identity'}
                             ) if (defined ($args{'percent_identity'}));   
        $seq_run->addBsmlAttr( 
            'percent_similarity', 
            $args{'percent_similarity'}
                             ) if (defined ($args{'percent_similarity'}));

        ## add hsp percent coverage stats
        $seq_run->addBsmlAttr( 
            'percent_coverage_refseq', 
            sprintf("%.1f", $args{'runlength'} / $args{'query_length'})
                             ) if (defined ($args{'runlength'}) && defined($args{'query_length'}));
                             
        $seq_run->addBsmlAttr(
            'percent_coverage_compseq', 
            sprintf("%.1f", $args{'comprunlength'} / $args{'hit_length'})
                             ) if (defined ($args{'comprunlength'}) && defined($args{'hit_length'}));
        ## ^^
        
        $seq_run->addBsmlAttr(
            'p_value', 
            $args{'p_value'}
                             ) if (defined ($args{'p_value'}));

        return $alignment_pair;
    }

    #no alignment pair matches, add a new alignment pair and sequence run
    #check to see if sequences exist in the BsmlDoc, if not add them with basic attributes
    my $seq;
    
    if( !( $doc->returnBsmlSequenceByIDR( "$args{'query_name'}")) ){
        $seq = $doc->createAndAddSequence( 
                        "$args{'query_name'}", 
                        "$args{'query_name'}", 
                        $args{'query_length'}, 
                        $ref_molecule->{$blast_program}, 
                        $args{'class'}
                                         );
        $doc->createAndAddSeqDataImport(
                        $seq, 
                        'fasta', 
                        $options{'query_file_path'}, 
                        '', 
                        $args{'query_name'}
                                       );
        $seq->addBsmlLink(
                        'analysis', 
                        '#' . $options{analysis_id}, 
                        'input_of'
                         );
    }
    
    if( !( $doc->returnBsmlSequenceByIDR( "$args{'dbmatch_accession'}")) ){
        $seq = $doc->createAndAddSequence(
                        "$args{'dbmatch_accession'}", 
                        "$args{'dbmatch_header'}", 
                        ($args{'hit_length'} || 0), 
                        $comp_molecule->{$blast_program}, 
                        $args{'class'} 
                                         );
        $doc->createAndAddSeqDataImport(
                        $seq, 
                        'fasta', 
                        $args{'search_database'}, 
                        '', 
                        $args{'dbmatch_accession'}
                                       );
    }

    ## see if the dbmatch_header format is recognized.  if so, add some cross-references
    if (defined $args{'dbmatch_header'}) {
        $doc->createAndAddCrossReferencesByParse( 
                    sequence => $seq, 
                    string => $args{orig_dbmatch_accession}
                                                );
    }
    
    $alignment_pair = $doc->returnBsmlSeqPairAlignmentR( $doc->addBsmlSeqPairAlignment() );
    
    ## to the alignment pair, add a Link to the analysis
    $alignment_pair->addBsmlLink(
                            'analysis', 
                            '#' . $options{analysis_id},
                            'computed_by'
                                );

    $alignment_pair->setattr( 
                            'refseq', 
                            "$args{'query_name'}"
                            ) if (defined ($args{'query_name'}));
    $alignment_pair->setattr(
                            'compseq', 
                            "$args{'dbmatch_accession'}"
                            ) if (defined ($args{'dbmatch_accession'}));
    
    BSML::BsmlDoc::BsmlSetAlignmentLookup( 
                            "$args{'query_name'}", 
                            "$args{'dbmatch_accession'}", 
                            $alignment_pair
                                         );

    $alignment_pair->setattr(
                            'refxref', 
                            $options{'query_file_path'}.':'.$args{'query_name'}
                            ) if (defined ($args{'query_name'}));                     
    $alignment_pair->setattr(
                            'refstart',
                            0
                            );
    $alignment_pair->setattr(
                            'refend',
                            $args{'query_length'}
                            ) if (defined ($args{'query_length'}));
    $alignment_pair->setattr(
                            'reflength',
                            $args{'query_length'}
                            ) if (defined ($args{'query_length'}));
    $alignment_pair->setattr(
                            'method',
                            $args{'blast_program'}
                            ) if (defined ($args{'blast_program'}));
    $alignment_pair->setattr(
                            'compxref',
                            $args{'search_database'}.':'.$args{'dbmatch_accession'} 
                            ) if ((defined ($args{'search_database'})) and (defined ($args{'dbmatch_accession'})));
    $alignment_pair->setattr(
                            'class',
                            'match'
                            );

    ## add average percent coverage values to alignment pair
    $alignment_pair->addBsmlAttr(
                            'percent_coverage_refseq', 
                            $args{'percent_coverage_refseq'}
                                ) if (defined ($args{'percent_coverage_refseq'}));   
    $alignment_pair->addBsmlAttr(
                            'percent_coverage_compseq', 
                            $args{'percent_coverage_compseq'}
                                ) if (defined ($args{'percent_coverage_compseq'}));   
    
    ## add percent identity and similarity values to alignment pair
    $alignment_pair->addBsmlAttr(
                            'percent_identity', 
                            $args{'percent_identity_total'}
                                ) if (defined ($args{'percent_identity_total'}));   
    $alignment_pair->addBsmlAttr(
                            'percent_similarity', 
                            $args{'percent_similarity_total'}
                                ) if (defined ($args{'percent_similarity_total'}));   

    my $seq_run = $alignment_pair->returnBsmlSeqPairRunR( $alignment_pair->addBsmlSeqPairRun() );

    $doc->createAndAddBsmlAttribute($seq_run, 'class', 'match_part');
    
    if( $args{'start_query'} > $args{'stop_query'} ) {
        $seq_run->setattr(
                        'refpos', 
                        $args{'stop_query'}-1
                         );
        $seq_run->setattr(
                        'runlength', 
                        $args{'start_query'} - $args{'stop_query'} + 1
                         );
        $seq_run->setattr(
                        'refcomplement', 
                        1 
                         );
    } else {
        $seq_run->setattr(
                        'refpos', 
                        $args{'start_query'} -1
                         );
        $seq_run->setattr(
                        'runlength', 
                        $args{'stop_query'} - $args{'start_query'} + 1
                         );
        $seq_run->setattr(
                        'refcomplement', 
                        0
                         );
    }

    #the database sequence is always 5' to 3'
    $seq_run->setattr(
                    'comppos', 
                    $args{'start_hit'} -1
                     ) if (defined  ($args{'start_hit'}));
    $seq_run->setattr(
                    'comprunlength',
                    ($args{'stop_hit'} - $args{'start_hit'} + 1)
                     ) if ((defined ($args{'start_hit'})) and (defined ($args{'stop_hit'})));
    $seq_run->setattr(
                    'compcomplement', 
                    0
                     );
    $seq_run->setattr(
                    'runscore',
                    $args{'bit_score'}
                     ) if (defined  ($args{'bit_score'}));
    $seq_run->setattr( 
                    'runprob',
                    $args{'e_value'}
                     ) if (defined  ($args{'e_value'}));
    $seq_run->addBsmlAttr(
                    'percent_identity', 
                    $args{'percent_identity'}
                         ) if (defined  ($args{'percent_identity'}));
    $seq_run->addBsmlAttr(
                    'percent_similarity', 
                    $args{'percent_similarity'}
                         ) if (defined  ($args{'percent_similarity'}));
    
    ## add hsp percent coverage stats
    $seq_run->addBsmlAttr( 
        'percent_coverage_refseq', 
        sprintf("%.1f", $args{'runlength'} / $args{'query_length'})
                         ) if (defined ($args{'runlength'}) && defined($args{'query_length'}));
                             
    $seq_run->addBsmlAttr(
        'percent_coverage_compseq', 
        sprintf("%.1f", $args{'comprunlength'} / $args{'hit_length'})
                         ) if (defined ($args{'comprunlength'}) && defined($args{'hit_length'}));
    ## ^^
        
    $seq_run->addBsmlAttr(
                    'chain_number',
                    $args{'chain_number'}
                         ) if (defined  ($args{'chain_number'}));
    $seq_run->addBsmlAttr(
                    'segment_number',
                    $args{'segment_number'}
                         ) if (defined  ($args{'segment_number'}));
    $seq_run->addBsmlAttr(
                    'p_value',
                    $args{'p_value'}
                         ) if (defined  ($args{'p_value'}));

    return $alignment_pair;
}

sub getAvgBlastPIdSim {
    my($hsps) = @_;
    my $sim_sum=0;
    my $id_sum=0;
    my $q_len_sum=0;
    my $numHsps = 0;
    
    # Group by query and target id
    my $hspsByQuery = &groupByMulti($hsps, ['query_protein_id', 'target_protein_id']);

    foreach my $queryId (keys %$hspsByQuery) {
        my $hspsByTarget = $hspsByQuery->{$queryId};

        foreach my $subjId (keys %$hspsByTarget) {
            my $shsps = $hspsByTarget->{$subjId};
#            my $querySeqLen = $shsps->[0]->{'query_seqlen'};
#            my $targetSeqLen = $shsps->[0]->{'target_seqlen'};
            
            foreach my $hsp(@{$shsps}) {
                ++$numHsps;
                my $q_seg_len = $hsp->{'query_fmax'} - $hsp->{'query_fmin'};
                $q_len_sum += $q_seg_len;
                $sim_sum += $hsp->{'percent_similarity'} * $q_seg_len;
                $id_sum += $hsp->{'percent_identity'} * $q_seg_len;
            }
        }
    }
    
    if ($numHsps == 0) {
        return undef;
    } else {
        return [($id_sum/$q_len_sum), ($sim_sum/$q_len_sum)];
    }
            
    
}

## Returns an array reference where 
##     [0] = query percent coverage, 
## and [1] = target percent coverage
sub getAvgBlastPPctCoverage {
    my($hsps) = @_;
    my $qsum=0;
    my $tsum=0;
    my $numHsps=0;

    # Group by query and target id
    my $hspsByQuery = &groupByMulti($hsps, ['query_protein_id', 'target_protein_id']);

    foreach my $queryId (keys %$hspsByQuery) {
        my $hspsByTarget = $hspsByQuery->{$queryId};

        foreach my $subjId (keys %$hspsByTarget) {
            ++$numHsps;
            my $shsps = $hspsByTarget->{$subjId};
            my $querySeqLen = $shsps->[0]->{'query_seqlen'};
            my $targetSeqLen = $shsps->[0]->{'target_seqlen'};

            my @queryIntervals = map { {'fmin' => $_->{'query_fmin'}, 'fmax' => $_->{'query_fmax'}, 'strand' => $_->{'query_strand'}} } @$shsps;
            my @targetIntervals = map { {'fmin' => $_->{'target_fmin'}, 'fmax' => $_->{'target_fmax'}, 'strand' => $_->{'target_strand'}} } @$shsps;

            my $mergedQueryIntervals = &mergeOverlappingIntervals(\@queryIntervals);
            my $mergedTargetIntervals = &mergeOverlappingIntervals(\@targetIntervals);

            my $queryHitLen = 0;
            my $targetHitLen = 0;

            map { $queryHitLen += ($_->{'fmax'} - $_->{'fmin'}); } @$mergedQueryIntervals;
            map { $targetHitLen += ($_->{'fmax'} - $_->{'fmin'}); } @$mergedTargetIntervals;

            $qsum += $queryHitLen / $querySeqLen;
            $tsum += $targetHitLen / $targetSeqLen;
        }
    }

    if ($numHsps == 0) {
        return undef;
    } else {
        return [($qsum/$numHsps*100.0), ($tsum/$numHsps*100.0)];
    }
}

# Generalized version of groupBy 
sub groupByMulti {
    my($arrayref, $keyFields) = @_;
    my $nKeys = scalar(@$keyFields);
    my $groups = {};

    foreach my $a (@$arrayref) {
        my @keyValues = map { $a->{$_} } @$keyFields;
        my $hash = $groups;

        for (my $i = 0;$i < $nKeys;++$i) {
            my $kv = $keyValues[$i];

            if ($i < ($nKeys-1)) {
                $hash->{$kv} = {} if (!defined($hash->{$kv}));
                $hash = $hash->{$kv};
            } else {
                $hash->{$kv} = [] if (!defined($hash->{$kv}));
                push(@{$hash->{$kv}}, $a);
            }
        }
    }
    return $groups;
}

# Generate a new set of intervals by merging any that overlap in the original set.
#
sub mergeOverlappingIntervals {
    my($intervals) = @_;

    # result set of intervals
    my $merged = [];

    # sort all intervals by fmin
    my @sorted = sort { $a->{'fmin'} <=> $b->{'fmin'} } @$intervals;
    
    # current interval
    my $current = undef;

    foreach my $i (@sorted) {
        if (!defined($current)) {
            # case 1: no current interval
            $current = $i;
        } else {
            # case 2: compare current interval to interval $i
            if ($i->{'fmin'} > $current->{'fmax'}) {   
                # case 2a: no overlap
                push(@$merged, $current);
                $current = $i;
            } elsif ($i->{'fmax'} > $current->{'fmax'}) {
                # case 2b: overlap, with $i ending to the right of $current
                $current->{'fmax'} = $i->{'fmax'};
            }
        }
    }
    push(@$merged, $current) if (defined($current));

    return $merged;
}


sub check_parameters {
    my $options = shift;
    
    unless (-e $options{'input'}) {
        $logger->logdie("input option not passed or does not exist!");
        exit(1);
    }

    unless (defined $options{'output'}) {
        $logger->logdie("output option not passed");
        exit(1);
    }

    ## handle some defaults
    if (! $options{'analysis_id'}) {
        $options{'analysis_id'} = 'unknown_analysis';
    }
    if (! $options{'query_file_path'}) {
        $options{'query_file_path'} = '';
    } else {
        open(IN, "< $options{query_file_path}") or
            $logger->logdie("Unable to open query_file $options{query_file_path} ($!)");
        while(<IN>) {
            chomp;
            if(/^>(.*)/) {
                $defline = $1;
                last;
            }
        }
        close(IN);
    }
}


