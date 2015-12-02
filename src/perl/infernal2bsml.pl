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
    The table output from Infernal's 'cmsearch' run.

B<--output,-o>
    Where the output bsml file should be

B<--project,-p>
    [DEPRECATED] The script will parse the project from the sequence id.
    Used in id generation.  It's the first token in the id.  (Ex. project.class.number.version)

B<--id_repository,-r>
    Used to make the ids (See Ergatis::IdGenerator for details)

B<--query_file_path,-q>
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
    be assigned to your hits.  Can be a directory with the CM IDs in the file names, or a file with all CMs. 

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

    The tbl output of infernal. Can be obtained by appending --tblout <filename> to the 'cmsearch' command 

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
my $DEBUG = 4;						  # Default DEBUG number
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


$options{'debug'} = $DEBUG if ! $options{'debug'};

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
print STDERR "Parsing a V1 infernal output\n";
$rawStructure = &parseInfernalRawV1($input);

print STDERR "Building a hash of Rfam types and descriptions\n";
my $cm_hash = build_cm_type_desc_hash();

#Create a bsml document out of the data
&createBsml($rawStructure, $cm_hash);
print STDERR "Finished creating BSML\n";
#Write the bsml file.
$gene_pred->writeBsml($output, '', $gzip);
print STDERR "Finished writing BSML\n";
# If we have been passed a annotation file list then we'll need to add the genome tag to 
# the output bsml files.

# Note that this is a total HACK
if($options{'annot_bsml_list'}) {
    &appendGenomeLink(&getInputFiles($options{'annot_bsml_list'}),[$options{'output'}]);
	print STDERR "Finished creating annotated BSML list\n";
}

exit(0);

###################################### SUBROUTINES ################################

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
    my $strand;
    while(<IN>) {
		#next if (/^#target/) || (/^#---/);
		next if (/^#/);
		my @fields = split(/\s+/);

		$infHit->{'seqId'} 		= $fields[0];
		$infHit->{'name'} 		= $fields[2];
		$infHit->{'cmId'} 		= $fields[3];
		$qstart 				= $fields[5];
		$qstop 					= $fields[6];
		$tstart 				= $fields[7];
		$tstop 					= $fields[8];
		$infHit->{'bit_score'}	= $fields[13];

        # Currently assuming the Rfam.cm general file is being used.  Store the path
		#if(/\# query CM file:\s+(\/.*\.cm)/) {
		#    $infHit->{'cmPath'} = $1;
		#}
        
		# Hit start should always be less than hit stop
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
		$project = $1 if($infHit->{'seqId'} =~ /^([^\.]+)/);

        my %tmpHash = %{$infHit};

        push(@{$retval}, \%tmpHash);
		%$infHit = ();	# Reset hash


    }

    close(IN);

    return $retval;
}

#######GENERATE BSML (helper functions follow)#############
sub createBsml {    
    my $infernalRaw = shift;
	my $cm_hash = shift;
    my %querySeqs;
    my %cmSeqs;
    my %spas;
    my %fts;
    my %countHits;
    
	my ($type, $desc);

    foreach my $match (@{$infernalRaw}) {

        #See if we have a type for the given Rfam CM. If not we're moving on.
        $type = $cm_hash->{$match->{'cmId'}}->{type};
        if(!$type) {
            next;
        } else {
			$desc = $cm_hash->{$match->{'cmId'}}->{desc};
		}

        #Make sure the covariance model has been added as a sequence
        &addCoModel(\%cmSeqs, $match);

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

# Build a hash consisting of the type and decription for each CM
sub build_cm_type_desc_hash {
    my $cmId = '';
    my $type = "";
    my $desc = "";

	my $cm_hash = {};
    my $file;

	# Determine if stockholm input is a file or a directory
	if (-d $options{'stockholm_path'}){
    	opendir(DIR, $options{'stockholm_path'}) || warn "can't opendir $options{'stockholm_path'}: $!";
		$cm_hash = parse_file($file, $cm_hash) while ($file = readdir(DIR));
		closedir(DIR);
	} elsif (-f $options{'stockholm_path'}) {
		$file = $options{'stockholm_path'};
		$cm_hash = parse_file($file, $cm_hash);
	} else {
		&_die("Stockholm path specified was neither a directory nor a file");
	}

	return $cm_hash;
}

# Parse the Stockholm Rfam file and push relevant parts into a hash-ref
sub parse_file {
	my ($file, $cm_hash) = @_;
	open(IN, $file ) or &_warn("Unable to open file $file for reading");
    
	my $cm = '';
	my $type = '';
	my $desc = '';
	
	while(<IN>) {
		if(/^\#=GF\sAC\s+(.+)/) {
			$cm = $1;
		}
		
        if(/^\#=GF\sTP\s+(.+)/) {
            my @types = split(/;/, $1);
            map {
                s/\s//g;
                if($RFAM_TP_TO_SO_TYPE->{$_}) {
                    $type = $RFAM_TP_TO_SO_TYPE->{$_};
                }
            }@types;
            $type = $options{'default_type'} if(!$type && $options{'default_type'});
            if(!$type) {
                print STDERR "WARN: ignoring $cm hit because it's type could not be determined\n";
            }
        }
        elsif(/^\#=GF\sDE\s+(.*)/) {
            $desc = $1;
        }

		if (m{^//}) {
			$cm_hash->{$cm}->{type} = $type;
			$cm_hash->{$cm}->{desc} = $desc;
		}
    }
        close(IN);
    return ($cm_hash);
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
