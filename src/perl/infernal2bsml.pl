#!/usr/bin/perl

=head1 NAME

infernal.pl - Turns infernal raw output into bsml. 

=head1 SYNOPSIS

USAGE: template.pl
            --input_file=/path/to/some/transterm.raw
            --output=/path/to/transterm.bsml
            --project=aa1
            --id_repository=/some/id_repository/dir
            --query_file_path=/input/file.fsa
            --gzip_output=1
          [ --infernal_v1
            --annot_bsml_list=/path/to/annot_bsml.list
            --stockholm_path=/path/to/stockholm/files
            --default_type=ncRNA
            --log=/path/to/file.log
            --debug=4
            --help
          ]

=head1 OPTIONS

B<--input_file,-i>
    The input file (should be prosite scan output)

B<--output,-o>
    Where the output bsml file should be

B<--project,-p>
    [DEPRECATED] The script will parse the project from the sequence id.
    Used in id generation.  It's the first token in the id.  (Ex. project.class.number.version)

B<--id_repository,-r>
    Used to make the ids (See Ergatis::IdGenerator for details)

B<--query_file_path,-g>
    Path to the query file (input fasta file) for infernal.

B<--gzip_output,-g>
    A non-zero value will result in compressed bsml output.  If no .gz is on the end of the bsml output name, one will
    be added.

B<--infernal_v1,-v1>
    Whether or not to expect infernal 1.0 or later output. If left off then the assumption is that the input is 
    pre-infernal 1.0

B<--annot_bsml_list,-v1>
    A list file containing the existing annotation bsml. This will activate the appending of the genome tag to the output bsml.

B<--stockholm_path,-sp>
    The path to the stockholm files that contain the alignments for the Rfam CMs. This is required if you want a good type to 
    be assigned to your hits.

B<--default_type,-dt>
    The default SO type to be assigned to either all hits or only those hits that don't have a valid SO mapping in their 
    stockholm file (if --stockholm_path is specified).

B<--log,-l>
    In case you wanted a log file.

B<--debug,-d>
    There are no debug statements in this program.  Sorry.

B<--help,-h>
    Displays this message.

=head1  DESCRIPTION

    Reads in infernal output and produces a bsml representation of the matches.  

=head1  INPUT

    The raw output of infernal.  

=head1  OUTPUT

    Generates a BSML document representing the ps_scan matches.

=head1  CONTACT

    Kevin Galens
    kgalens@tigr.org

=cut


use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

use strict;
use warnings;
use BSML::GenePredictionBsml;
use BSML::BsmlBuilder;
use Ergatis::IdGenerator;
use Ergatis::Logger;
use Chado::Gene;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use XML::Twig;
#use Pod::Usage;

###############GLOBALS######################
my %options = ();                     #Options hash
my $input;
my $output;
my $gene_pred;

my $idMaker;                          #Id generator object
my $project;
my $gzip;
my $defline;

# The following is a mapping of all the supported RFAM types 
# and their associated SO type. These types were derived by 
# manually matching the TP line data from RFAM stockholm files
# with appropriate SO terms. If the TP line of the associated
# stockholm file does not have a value in SO then the behaviour 
# will default to whetever type is specified in --default_type.
# If no value is specified then the hit will be ignored.

my $RFAM_TP_TO_SO_TYPE = {
    'snoRNA' => 'snoRNA',
    'rRNA' => 'rRNA',
    'tRNA' => 'tRNA',
    'riboswitch' => 'riboswitch',
    'miRNA' => 'miRNA',
    'IRES' => 'internal_ribosome_entry_site',
    'ribozyme' => 'ribozyme',
    'snRNA' => 'snRNA'
    };


############################################


#Get the options.
my $results = GetOptions (\%options, 
                          'input|i=s',
                          'output|o=s',
                          'project|p=s',
                          'id_repository|r=s',
                          'query_file_path|q=s',
                          'gzip_output|g=s',
                          'infernal_v1|v1',
                          'stockholm_path|sp:s',
                          'annot_bsml_list|al:s',
                          'default_type|dt:s',
                          'log|l=s',
                          'command_id=s',       ## passed by workflow
                          'logconf=s',          ## passed by workflow (not used)
                          'debug=s',
                          'help|h') || pod2usage();


#Make the logger
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

#Create the gene prediction bsml document
$gene_pred = new BSML::GenePredictionBsml( 'infernal', $options{'output'});
my $doc = $gene_pred->{'doc'}; #new BSML::BsmlBuilder();    #Bsml Builder object
#If the help option is used, spit out documentation
&_pod if( $options{'help'} );

#Check in input parameters and set up some variables.
&check_parameters(\%options);

my $rawStructure;

#Parse the input file
#if($options{'infernal_v1'}) {
    print STDERR "Parsing a V1 infernal output\n";
    $rawStructure = &parseInfernalRawV1($input);
#}
#else {
#    $rawStructure = &parseInfernalRawPreV1($input);    
#}

#Create a bsml document out of the data
&createBsml($rawStructure);

#Write the bsml file.
$gene_pred->writeBsml($output, '', $gzip);

# If we have been passed a annotation file list then we'll need to add the genome tag to 
# the output bsml files.

# Note that this is a total HACK
if($options{'annot_bsml_list'}) {
    &appendGenomeLink(&getInputFiles($options{'annot_bsml_list'}),[$options{'output'}]);
}

exit(0);

###################################### SUBROUTINES ################################

#This is the parsing function.  Kinda gross.
sub parseInfernalRawPreV1 {
    my $inputFile = shift;
    my @prevLines;
    my $infHit;
    my $retval;

    open(IN, "< $inputFile") or &_die("Unable to open input $inputFile");

    while(<IN>) {
        #A crappy (har, har) regex.
        if(/^\shit\s.*?(\d+).*?:\s+(\d+)\s+(\d+)\s+([.\d]+)/) {
            my ($hitNo, $start, $stop) = ($1,$2,$3);
            $infHit->{'bit_score'} = $4;

            $start--; #For interbase numbering

            #Set all the stuff and things
            if($prevLines[1] =~ /sequence\:\s+(.*)\:\:(\d+)-(\d+)/) {
                $infHit->{'seqId'} = $1;
                $project = "";
                $project = $1 if($infHit->{'seqId'} =~ /^([^\.]+)/);
                $logger->logdie("Unable to parse project from $infHit->{'seqId'}") unless($project);
                $start+=($2-1);
                $stop+=($2-1);
            } 

            $infHit->{'cmId'} = $1 if($prevLines[0] =~ /cmsearch\s.*\/(.*?\.cm)\s+/);
            $infHit->{'strand'} = 0;  #0 means forward strand.

            $infHit->{'domain_num'} = $hitNo;
            
            if($start > $stop) {
                my $tmp = $start;
                $start = $stop;
                $stop = $tmp;
                $infHit->{'strand'} = 1;
            }            
            $infHit->{'start'} = $start;
            $infHit->{'stop'} = $stop;

            $infHit->{'runlength'} = $stop-$start;

            my $curLine=$_;
            my ($cmStart, $cmStop);
            my $cmLine = 1;
            
            HIT: while(<IN>) {
                if(/^\s+(\d+)/ && !defined($cmStart) && $cmLine) {
                    $cmStart = $1-1;
                } 
                if(/^\s+\d+.*?(\d+)\s*$/ && $cmLine) {
                    $cmStop = $1;
                } 
                if(/\sCPU\stime/) {
                    last HIT;
                }
                $cmLine = (($cmLine-1)**2) if(/^\s+\d/);
                
            }

            $_ = $curLine;
            $infHit->{'comprunlength'} = $cmStop-$cmStart;

            
            #print Dumper($infHit);
            my %tmpHash = %{$infHit};
            push(@{$retval}, \%tmpHash);

        } else {
            if(@prevLines == 2) {
                shift @prevLines;
                push(@prevLines, $_);
            } elsif(@prevLines < 2) {
                push(@prevLines, $_);
            } else {
                &_die("Collecting too many lines");
            }
        }

     }
    close(IN);

#    $doc->createAndAddAnalysis(
#                               id => "infernal_analysis",
#                               sourcename => $options{'input'},
#                               );

    return $retval;

}

sub parseInfernalRawV1 {
    my $inputFile = shift;
    my @prevLines;
    my $infHit;
    my $retval;

    open(IN, "< $inputFile") or &_die("Unable to open input $inputFile");

    my $seqstart;
    my $seqstop;
    my $qstart;
    my $qstop;
    my $tstart;
    my $tstop;
    my $strand = 0;
    while(<IN>) {

        # HACK-O-RAMA: Assumes we name our cm files the Rfam_ID.cm.
        if(/\# command:.*\/.*\/(.*)\.cm/) {
            $infHit->{'cmId'} = $1;
        }
        elsif(/^ CM\: (.*)/) {
            $infHit->{'name'} = $1;
        }
        elsif(/ \>(.*)\:\:(\d+)-(\d+)/) {
            ($infHit->{'seqId'}, $seqstart, $seqstop) = ($1,$2,$3);
        }
        elsif(/ Query = (\d+) - (\d+), Target = (\d+) - (\d+)/) {
            ($qstart, $qstop, $tstart, $tstop) = ($1,$2,$3,$4);
        }
        elsif(/Score = (\S+),/) {
            my($score) = ($1);
            $infHit->{'bit_score'} = $score;
            $project = $1 if($infHit->{'seqId'} =~ /^([^\.]+)/);
            
            $tstart += $seqstart;
            $tstop += $seqstop;
            
            $infHit->{'strand'} = 0;
            if($tstart > $tstop) {
                my $tmp = $tstart;
                $tstart = $tstop;
                $tstop = $tmp;
                $infHit->{'strand'} = 1;
            }           
            $infHit->{'start'} = $tstart;
            $infHit->{'stop'} = $tstop;
            $infHit->{'runlength'} = $tstop - $tstart;
            $infHit->{'comprunlength'} = abs($qstop - $qstart);
            my %tmpHash = %{$infHit};

            push(@{$retval}, \%tmpHash);
        }
    }

    close(IN);

#    $doc->createAndAddAnalysis(
#                               id => "infernal_analysis",
#                               sourcename => $options{'input'},
#                               version => '1.0',
#                               );
    return $retval;
}

#######GENERATE BSML (helper functions follow)#############
sub createBsml {    
    my $infernalRaw = shift;
    my %querySeqs;
    my %cmSeqs;
    my %spas;
    my %fts;
    my %countHits;
    
    foreach my $match (@{$infernalRaw}) {

        #See if we can parse out a useful type. If not we're moving on.
        my($type,$desc) = &getTypeAndDesc($match->{'cmId'});
        if(!$type) {
            next;
        }

        #Has the query seq been added?
#        &addQuerySeq(\%querySeqs, \%fts, $match);
        
            
        #Make sure the covariance model has been added as a sequence
        &addCoModel(\%cmSeqs, $match);
        

        #Add the Sequence Pairs Currently not doing this because it seems redundant.
#        &addSeqPair(\%spas, $match);
        


        # First, create the gene feature
        my $currGene = new Chado::Gene ( $idMaker->next_id( 'type' => 'gene',
                                                              'project' => $project ),
                                         $match->{'start'}, $match->{'stop'},
                                         $match->{'strand'},
                                         $match->{'seqId'}
                                         );
        
        ## Next, with the same coords, create the RNA feature
        my $rna_feat_id = $currGene->addFeature( $idMaker->next_id ( 'type' => $type,
                                                                        'project' => $project ),
                                                 $match->{'start'}, $match->{'stop'},
                                                 $match->{'strand'},                                                 
                                                  $type
                                                  );
        ## add name attribute to feature
        $currGene->addFeatureAttribute(
                                       $rna_feat_id,
                                       'gene_product_name',
                                       $desc,
                                       );
        
        ## add name attribute to feature
        $currGene->addFeatureAttribute(
                                       $rna_feat_id,
                                       'gene_product_name_source',
                                       $match->{'cmId'},
                                       );

        ## add score attribute to feature
        $currGene->addFeatureAttribute(
                                       $rna_feat_id,
                                       'score',
                                       $match->{'bit_score'},
                                       );
        foreach my $t (qw(exon CDS)) {
            
            # Add the exon or CDS to the gene model object
            $currGene->addFeature( $idMaker->next_id( 'type' => $t,
                                                  'project' => $project ),
                                   $match->{'start'}, 
                                   $match->{'stop'},
                                   $match->{'strand'},         
                                   $t
                                   );
            
        }
        # Handle Group now:
        my $count = $currGene->addToGroup( $currGene->getId, { 'all' => 1} );
        $gene_pred->addGene($currGene);
        #print "About to addBsmlLink\n";
#        $feat->addBsmlLink('analysis', '#infernal_analysis', 'computed_by');
        #print "Adding intervalLoc\n";
#        $feat->addBsmlIntervalLoc($match->{'start'}, $match->{'stop'}, $match->{'strand'});
        #print "Add the $type to the feature group\n";
        
        #print "Finished and writing.\n";
        
    }
}

sub addQuerySeq {
    my ($querySeqs, $fts, $match) = @_;
    
    unless($querySeqs->{$match->{'seqId'}}) {
        $querySeqs->{$match->{'seqId'}}  = 
            $doc->createAndAddSequence($match->{'seqId'},undef, 
                                       undef, 'dna', 'assembly');
        $querySeqs->{$match->{'seqId'}}->addBsmlLink('analysis', '#infernal_analysis', 'input_of');
        $querySeqs->{$match->{'seqId'}}->addBsmlAttr('defline', $defline);
        $fts->{$match->{'seqId'}} = 
            $doc->createAndAddFeatureTable($querySeqs->{$match->{'seqId'}});
    }
    
}

sub addCoModel {
    my ($cmSeqs, $match) = @_;

    unless($cmSeqs->{$match->{'cmId'}}) {
        $cmSeqs->{$match->{'cmId'}} = 
            $doc->createAndAddSequence($match->{'cmId'}, undef,
                                       undef, 'rna', 'covariance_model');
        $cmSeqs->{$match->{'cmId'}}->addBsmlLink('analysis', '#infernal_analysis', 'input_of');
    }

}

sub addSeqPair {
    my($spas, $match) = @_;

    unless($spas->{$match->{'seqId'}.$match->{'cmId'}}) {
        $spas->{$match->{'seqId'}.$match->{'cmId'}} = 
            my $aln = $doc->createAndAddSequencePairAlignment( refseq   => $match->{'seqId'},
                                                     compseq  => $match->{'cmId'},
                                                     complength => "",
                                                     class => "match",
                                                     refstart => 0,);
            $aln->addBsmlLink('analysis', '#infernal_aalysis', 'computed_by');
    }
    
    my $spr = $doc->createAndAddSequencePairRun(alignment_pair => 
                                                $spas->{$match->{'seqId'}.$match->{'cmId'}},
                                                runscore => $match->{'score'},
                                                runlength => $match->{'runlength'},
                                                comprunlength => $match->{'comprunlength'},
                                                refpos => $match->{'start'},
                                                refcomplement => $match->{'strand'},
                                                comppos => 0,
                                                compcomplement => 0,
                                                );

    $doc->createAndAddBsmlAttribute($spr, 'class', 'match_part');
    $doc->createAndAddBsmlAttribute($spr, bit_score => $match->{'bit_score'});


}

#Given our current default location for RFAM HMMs and CMs
# were searching for a more informative title.
sub getTypeAndDesc {
    my $cmId = shift;
    my $type = "";
    my $desc = "";

    my $cm = $cmId;
    $cm = $1 if($cmId =~ /^([^.]+)/);
    
    opendir(DIR, $options{'stockholm_path'}) || warn "can't opendir $options{'stockholm_path'}: $!";
    my @files = grep { /$cm/ } readdir(DIR);
    closedir(DIR);
    if(@files > 0) {
        open(IN, "<$options{'stockholm_path'}/$files[0]" ) or &_warn("Unable to open file $files[0]");
        while(<IN>) {
            if(/^\#=GF\sTP\s+(.*)/) {
                my @types = split(/;/, $1);
                map {
                    s/\s//g;
                    if($RFAM_TP_TO_SO_TYPE->{$_}) {
                        $type = $RFAM_TP_TO_SO_TYPE->{$_};
                    }
                }@types;
                $type = $options{'default_type'} if(!$type && $options{'default_type'});
                if(!$type) {
                    print STDERR "WARN: ignoring $cmId hit because it's type could not be determined\n";
                }
            }
            elsif(/^\#=GF\sDE\s+(.*)/) {
                $desc = $1;
            }
            
        }
        close(IN);
    }
    return ($type,$desc);
}
#############END OF BSML GENERATION#####################

sub check_parameters {
    
    #input must be passed and must exist
    if($options{'input'}) {
        &_die("input option ($options{input}) does not exist") unless(-e $options{'input'});
    } else {
        &_die("input option was not passed and is required");
    }
    $input = $options{'input'};

    #output must be passed
    &_die("output option was not passed and is required") unless($options{'output'});
    $output = $options{'output'};

    #the id_repository should be passed (perhaps it should also exist, but whatever).
    unless($options{'id_repository'}) {
        &_die("option id_repository is required (for id generation).".
              " See Ergatis::IdGenerator for more information");
    }
    $idMaker = new Ergatis::IdGenerator( id_repository => "$options{'id_repository'}" );
    
    #If the query file path was given, parse out the definition line.
    if($options{'query_file_path'}) {
        open(IN, "<$options{query_file_path}") or
            &_die("Unable to open $options{query_file_path}");
        
        while(<IN>) {
            chomp;
            if(/^>(.*)/) {
                $defline = $1;
                last;
            }
        }
        close(IN);
    }
    if(!$options{'stockholm_path'} && !$options{'default_type'}) {
        &_die("You must specify either --stockholm_path or --default_type or both");
    }

}

#------------------------------------------------------------------------------------------#
# From here down is a super duper hack. This was at one point another script that has been #
# copied into this script because we cannot do conditionals in ergatis.                    #

sub appendGenomeLink {

    my $annotFiles = shift;
    my $analFiles = shift;

    my @annotationFiles = @$annotFiles;
    my @analysisFiles = @$analFiles;

    # First make a mapping of sequence_id to genome tag.
    my $sequence_id_to_genome = {};
    foreach my $file (@annotationFiles) {
        
        my $genome_tag;
        my $genome_id;
        my $found_sequence = 0;
        my $twig = new XML::Twig( twig_handlers => {
            'Genomes' => sub {
                my($t, $genome) = @_;
                $t->set_pretty_print('indented');
                $genome_tag = $genome->sprint(0,0);
                $genome_id = $genome->child('Genome')->{'att'}->{'id'};
            },
            'Sequence[@class="assembly"]' => sub {
                my($t, $seq) = @_;
                if($seq->{'att'}->{'class'} eq 'assembly') {
                    my $ident = $seq->first_child('Seq-data-import')->{'att'}->{'identifier'};
                    if($ident =~ /\|/) {
                        print STDERR "Found some pipes \| in identifier $ident. Replacing them with '_'. This is a no smoking script.\n";
                        $ident =~ s/\|/\_/g;
                    }
                    $sequence_id_to_genome->{$ident} = {
                        'tag' => $genome_tag,
                        'id' => $genome_id
                        };
                    $found_sequence = 1;
                }
            }
        },
        ignore_elts => {
            'Feature-tables' => 1,
            'Feature-group' => 1,
            'Feature' => 1,
            'Sequence[@class!="assembly"]' => 1});
        
        # Parse the file
        my $openFile;
        if( -e $file ) {
            open($openFile, "< $file");
        } elsif( -e $file.".gz" ) {
            open($openFile, "<:gzip", $file.".gz");
        } else {
        }
        $twig->parse($openFile);
        
        print STDERR "Couldn't find genome tag or Sequence tag for $file\n" if !$found_sequence;
    }
    
    print STDERR "Done making sequence_id->genome_tag mapping\n";
    
    # Next make a mapping of filename to sequence_id
    my $analysis_file_name_to_sequence_id = {};
    foreach my $file (@analysisFiles) {
        
        my $genome_tag;
        my $genome_id;
        my $twig = new XML::Twig( twig_handlers => {
            'Sequence' => sub {
                my($t, $seq) = @_;
                if($seq->{'att'}->{'class'} eq 'assembly') {
                    print STDERR $seq->{'att'}->{'id'}."\n";
                    $analysis_file_name_to_sequence_id->{$file} = $seq->{'att'}->{'id'};
                }
            }  
        });
        
        # Parse the file
        my $openFile;
        if( -e $file ) {
            open($openFile, "< $file");
        } elsif( -e $file.".gz" ) {
            open($openFile, "<:gzip", $file.".gz");
        } else {
        }
        $twig->parse($openFile);
        
    }
    print STDERR "Done making filename->sequence_id\n";
    
    # Lastly we'll add the genome text and the link.
    foreach my $file (@analysisFiles) {
        
        # Parse the file
        my $openFile;
        if( -e $file ) {
            open($openFile, " $file");
        } elsif( -e $file.".gz" ) {
            open($openFile, ":gzip", $file.".gz");
        } else {
        }
        
        my $in_sequence = 0;
        my $fname = $file;
        $fname =~ s/.*\/([^\/]+)/$1/;
        open(my $outputFile, ">".$file.".new");
        while(<$openFile>) {
            chomp;
            my $line = $_;
            if($line =~ /<Definitions>/) {
                if($analysis_file_name_to_sequence_id->{$file}) {
                    print $outputFile $line.($sequence_id_to_genome->{$analysis_file_name_to_sequence_id->{$file}}->{'tag'}."\n");
                }
                else {
                    print STDERR "Wasn't able to add a genome tag to $file\n";
                }
            }
            elsif($line =~ /<Sequence.*assembly.*/) {
                $in_sequence = 1;
            }
            elsif($line =~ /(.*)<\/Sequence>/ && $in_sequence) {
                $in_sequence = 0;
                if($analysis_file_name_to_sequence_id->{$file}) {
                    print $outputFile "$1<Link rel=\"genome\" href=\"\#".
                        ($sequence_id_to_genome->{$analysis_file_name_to_sequence_id->{$file}}->{'id'})."\"></Link>\n";
                    print $outputFile $line."\n";
                }
            }
            else {
                print $outputFile "$line\n";
            }
        }
        my $er =  `cp $file\.new $file`;
        if(!$er) {
            print STDERR `rm $file\.new`;
        }
        else {
            print STDERR "Issues copying $file\.new to $file\n";
        }
    }
}

sub getInputFiles {
    my $inputStr = shift;
    my @files;

    #comma unseparate
    my @tokens = split( /,/, $inputStr);

    foreach my $token ( @tokens ) {
        my $tfh;
        $token =~ s/\.gz$//;
        if( -e $token ) {
            open($tfh, "< $token") or &_die("Can't open $token ($!)");
        } elsif( -e $token.".gz" ) {
            my $compFile = $token.".gz";
            open($tfh, "<:gzip", "$compFile") or &_die("Can't open $compFile ($!)");
        } else {
            die("Can't find file $token");
        }
        

        my $isList = 0;

        while( my $line = <$tfh> ) {
            next if($line =~ /^\s+$/);
            
            #It's a bsml file.
            if($line =~ /</ && !$isList) {
                push(@files, $token);
                last;
            #It's a list.
            } elsif($line =~ m|\.bsml|) {
                chomp $line;
                &_die("$line doesn't exist (from list $token) ($!)")
                    unless(-e $line || -e $line.".gz");
                push(@files, $line);
                $isList = 1;
            }
            
        }

    }

    return \@files;
    
}

sub _pod {
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

#DIE!
sub _die {
    my $msg = shift;
    $logger->logdie($msg);
}
######EOF###############################
