#!/usr/bin/perl

#./new_auto_gene_curation --input_file mdg.scaffold.1.glimmer3.bsml --output_dir ./ --hmm_input /usr/local/annotation/MOORE/output_repository/hmmpfam/5925_protein/ --ber_input /usr/local/annotation/MOORE/output_repository/ber/5925_default/ --names_dump /usr/local/annotation/MOORE/data_dbs/tax_dump/names.dmp --nodes_dump /usr/local/annotation/MOORE/data_dbs/tax_dump/nodes.dmp --base_name tmp.base.name --debug 5

####################################

### Use these things. ###########
use strict;
use warnings;
use IntervalTree;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use XML::Twig;
use File::Find;
use Ergatis::Logger;
use Time::localtime;
use BSML::BsmlBuilder;
use MLDBM 'DB_File';
$|++;
#################################

### Globals ##############
my @inputFiles;
my $sequenceSource;
my $tmpDirectory;           
my $outDirectory;
my $berFiles;
my $hmmFiles;
my %options = ();
my $names;
my $nodes;
my $baseName = "";
my $logger;
my $features;
my $genes;
my $startSites;
my $familyLookup;
my %hmmInfo;
my $changedFeatures;
my $lengths;
my $featureSeqElems;
my $deletedFeatures;
###########################

######## Parameters #################
my $overlapCutoff = 30;
my $btabEvalCutoff1 = 1e-10;
my $btabEvalCutoff2 = 1e-5;
my $moveCutoff = 15;
my $delPerOverlapCutoff = 0.6;
my $minGeneLength = 90;
my $strongEvCutoff = 5;
my $weakEvCutoff = 10;
my $berExtension = 300;
#####################################

####### Option Gathering and Logger Setup ########
my $results = GetOptions (\%options, 
                          'input_list|i=s',
                          'input_file|f=s',
                          'output_dir|o=s',
                          'hmm_input|h=s',
                          'hmm_info_db|m=s',
                          'ber_input|b=s',
                          'names_dump|n=s',
                          'nodes_dump|d=s',
                          'tmp_directory|t=s',
                          'ber_extension|e=s',
                          'base_name|b=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || &_pod;

&checkParametersAndLog(\%options);
##################################################

############## MAIN ##############################
foreach my $inputFile (@inputFiles) {
    
    #First we need to prepare.  Parse taxonomy information
    $familyLookup = new FamilyLookup( ( 'names' => $names, 'nodes' => $nodes ) );
    print "preprocessing taxonomy information\n" if($logger->is_debug);
    my $status = $familyLookup->preprocessTaxonomyInformation();
    $logger->logdie($status) unless($status == 1);

    #Make new BSML file
    my $outBsml = "$outDirectory/$baseName.all.bsml";

    my $out;
    open($out, "> $outBsml") or 
        $logger->logdie("Unable to open $outBsml for writing ($!)");

    #This is used to make sure we gather all the sequences before we find any features
    #which is dealt with in the next twig.
    my $seqTwig = new XML::Twig( twig_handlers => {
        'Sequence' => \&collectSeqs });
    $seqTwig->parsefile($inputFile);
                                 
    
    #Create a twig object
    my $twig = new XML::Twig( twig_roots => 
                              { 
                                  'Sequence'  => \&handleSequence
                              } ,
                              twig_print_outside_roots => $out,
                              pretty_print => 'indented',
                              );

    #Parse input file
    print STDOUT "parsing $inputFile\n" if($logger->is_debug);
    $twig->parsefile($inputFile);
    print STDOUT "wrote new file to $outBsml\n" if($logger->is_debug);

    &printChangedFeatures;
}
##################################################

######## SUB ROUTINES (alphabetically sorted) ############################
sub addFeature {
    my ($feature, $featureTable, $bsmlObj) = @_;

    $logger->logdie("Feature does not have new id (old id: ".$feature->{'old_id'})
        unless($feature->{'new_id'});

    my $newId = $feature->{'new_id'};

    my ( $title, $class) = ( $feature->{'title'}, $feature->{'class'});  

    my $bsmlFeature = $bsmlObj->createAndAddFeature( $featureTable, $newId, $title, $class );

    my $intLoc = $feature->{'interval-loc'};
    $logger->logdie("Unable to parse Interval-loc element from $newId") unless($intLoc);
    my ($startpos, $endpos, $complement) = ($intLoc->att('startpos'), $intLoc->att('endpos'), 
                                            $intLoc->att('complement') );

    $logger->logdie("cannot parse startpos, endpos or complement from $newId\n") 
        unless(defined($startpos) &&
               defined($endpos)   &&
               defined($complement) );

    $bsmlFeature->addBsmlIntervalLoc( $startpos, $endpos, $complement );

    $bsmlFeature->addBsmlLink( 'analysis', '#auto_gene_curation_analysis', 'computed_by' );

    return $bsmlFeature;

}

sub analyzeBerEvidence {
    my ($berFile, $feature) = @_;
    
    #The e-value cutoff score is dependent on the length of the gene
    my $evalCutoff = $feature->{'endpos'} - $feature->{'startpos'} > 900  ?
        $btabEvalCutoff1 : $btabEvalCutoff2;

    my $twig = new XML::Twig( twig_handlers => {
        'Definitions' => sub {&definitionHandler(@_,$feature)} } );
    $twig->parsefile($berFile);

  
    
}

sub analyzeHmmEvidence {
    my ($hmm, $feature) = @_;

    $logger->logdie("$hmm does not exist") unless(-e $hmm);

    #Get the poly id so we can find the HMM file.
    my $cds_id = $feature->{'old_id'};
    my $poly_id = $genes->{$cds_id}->{'polypeptide'};
    $logger->logdie("Could not find polypeptide id for $cds_id") unless($poly_id);

    #Create a twig and a return array
    my $twig = new XML::Twig( twig_handlers => 
                              { 
                                  'Seq-pair-alignment'  => 
                                      sub { &handleHmmSPA(@_, $feature) }
                              });
    $twig->parsefile($hmm);
}

sub changeFeature {
    my ($changeFeature, $newStart) = @_;

    #The genes object holds the ids for the features, not the features themselves;
    foreach my $type ( keys %{$genes->{$changeFeature->{'old_id'}}} ) {
        next if($type eq 'feature-group' || $type eq 'fgm');

        my $featureId = $genes->{$changeFeature->{'old_id'}}->{$type};

        #First change the bsml id
        my $feature = $features->{$featureId};
        $feature->{'new_id'} = $1.($2+1) if($features->{$featureId}->{'old_id'} =~ /(.*\.)(\d+)/);
        $logger->logdie("Could not change the old id. ".$feature->{'old_id'}) unless($feature->{'new_id'});

        $logger->logdie("Could not find bsml_feature_object for $featureId") unless($feature->{'bsml_feature_object'});
        $feature->{'bsml_feature_object'}->set_att('id', $feature->{'new_id'});

        if(exists($featureSeqElems->{$featureId})) {
            #Update the sequence id
            $featureSeqElems->{$featureId}->{'seqId'} = $feature->{'new_id'};
        }

        #Now change the feature group
        $genes->{$changeFeature->{'old_id'}}->{'fgm'}->{$type}->set_att('featref', $feature->{'new_id'});
        
        #For debugging, keep the old start and stop.
        my ($oldStart,$oldStop) = $feature->{'complement'} ? ($feature->{'endpos'}, $feature->{'startpos'}) :
            ($feature->{'startpos'}, $feature->{'endpos'});

        #And now the coordinate.
        if($feature->{'complement'}) {
            $feature->{'endpos'} = $newStart;
            $feature->{'interval-loc'}->set_att('endpos', $newStart);
        } else {
            $feature->{'startpos'} = $newStart;
            $feature->{'interval-loc'}->set_att('startpos', $newStart);
        }

        #Just a check.
        if($feature->{'startpos'} > $feature->{'endpos'}) {
            print STDOUT $feature->{'old_id'}."\n";
            print STDOUT $feature->{'new_id'}."\n";
            print STDOUT "Old Start: $oldStart\n";
            print STDOUT "Old Stop:  $oldStop\n";
            print STDOUT "Startpos: ".$feature->{'startpos'}."\n";
            print STDOUT "Endpos:   ".$feature->{'endpos'}."\n";
            print STDOUT $feature->{'complement'}."\n";
            
            $logger->logdie("I have no clue what happened but the endpos is now lower than the startpos.  Good luck");
        }

        #Keeps a records of features that have changed.
        $changedFeatures->{$featureId} = 1;

        #Change the sequence element that is related to the feature that has changed.
        if(exists($featureSeqElems->{$feature->{'old_id'}})) {
            
            #Update the sequence id.
            $featureSeqElems->{$feature->{'old_id'}}->{'seqId'} = $feature->{'new_id'}."_seq";

            #We don't have fasta sequence for it.  So we should get that.  And since we are only ever
            #making genes shorter, we can use the existing fasta and just take off some nucleotides/amino acids.
            my $source = $featureSeqElems->{$feature->{'old_id'}}->{'source'};
            $logger->logdie("Attribute source does not exist in Seq-data-import element of Sequence ".
                            $featureSeqElems->{$feature->{'old_id'}}->{'seqId'}) unless($source);
            
            my $ident = $featureSeqElems->{$feature->{'old_id'}}->{'identifier'};
            $logger->logdie("Unable to parse identifier from Seq-data-import element of Sequence ".
                            $featureSeqElems->{$feature->{'old_id'}}->{'seqId'}) unless($ident);

            $logger->logdie("Sorry, I only deal with fasta input for Seq-data-import elements ".
                            $featureSeqElems->{$feature->{'old_id'}}->{'seqId'}) 
                unless($featureSeqElems->{$feature->{'old_id'}}->{'format'} eq 'fasta');

            my $difference = $feature->{'complement'} ? $oldStart - $newStart : $newStart - $oldStart;
          
            $logger->logdie("The gene has gotten bigger and therefore cannot make a new fasta file from $source")
                if($difference <= 0);

            my $fsaFile = &makeFsaFile($source, $difference/3, $ident, $feature->{'new_id'});

            $featureSeqElems->{$feature->{'old_id'}}->{'source'} = $fsaFile;
            $featureSeqElems->{$feature->{'old_id'}}->{'identifier'} = $feature->{'new_id'};

            print STDOUT "Printing the seq elem for ".$feature->{'old_id'}."\n" if($logger->is_debug);
            print STDOUT $featureSeqElems->{$feature->{'old_id'}}->{'seqId'} if($logger->is_debug);
            
        } 
    }

}

sub checkParametersAndLog {
    my $opts = shift;

    #Setup Logger
    my $logfile = $opts->{'log'} || Ergatis::Logger::get_default_logfilename();
    $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$opts->{'debug'});
    $logger = $logger->get_logger();
    print "setup logger\n" if($logger->is_debug);

    #If -h, then pod
    &_pod if($opts->{'help'});

    #Input is required.
    unless($opts->{'input_list'} || $opts->{'input_file'}) {
        $logger->logdie("Either input_list or input_file option is required");
    } else {
        if($opts->{'input_list'}) {
            print "retrieving input from list\n" if($logger->is_debug);
            open(IN, "< $opts->{'input_list'}") or
                $logger->logdie("Unable to open input_list $opts->{'input_list'} ($!)");
            @inputFiles = <IN>;
            chomp(@inputFiles);
            close(IN);
            
        }
        if($opts->{'input_file'}) {
            print "using file as input\n" if($logger->is_debug);
            push(@inputFiles, $opts->{'input_file'});
        }
    }

    #Output directory is required
    unless($opts->{'output_dir'}) {
        $logger->logdie("Option output_dir is required");
    } else {
        $logger->logdie("$opts->{'output_dir'} is not a directory") unless(-d $opts->{'output_dir'});
    }
    $outDirectory = $opts->{'output_dir'};
    $outDirectory =~ s/\/$//;

    #Option temp directory
    unless($opts->{'tmp_dir'}) {
        $tmpDirectory = $outDirectory;
    } else {
        $logger->logdie("$opts->{'tmp_dir'} is not a directory") unless(-d $opts->{'tmp_dir'});
        $tmpDirectory = $opts->{'tmp_dir'};
    }
    $tmpDirectory =~ s/\/$//;

    #Evidence Input (HMM and BER).
    if($opts->{'hmm_input'}) {
        if(-d $opts->{'hmm_input'}) {
            print STDOUT "gathering hmm input files\n" if($logger->is_debug);
            #Input was a directory, get all the bsml files
            find( sub {  $hmmFiles->{$1} = $File::Find::name 
                             if($File::Find::name =~ /.*\.bsml/ && /polypeptide\.(\d+)\./) },
                  $opts->{'hmm_input'});
                
        } elsif(-e $opts->{'hmm_input'}) {
            open(IN, "< $opts->{'hmm_input'}") or
                $logger->logdie("Unable to open hmm_input file $opts->{'hmm_input'} ($!)");
            my $line;
            while($line = <IN>) {
                chomp $line;
                $hmmFiles->{$1} = $line if($line =~ /polypeptide.(\d+)\./);
            }
            close(IN);
            
        } else {
            $logger->logdie("Option hmm_input ($opts->{'hmm_input'}) is not a file or a directory");
        }
    } else {
        $logger->logdie("Option hmm_input is required");
    }

    if($opts->{'ber_input'}) {
        if(-d $opts->{'ber_input'}) {
            #Input was a directory, get all the bsml files
            print STDOUT "gathering ber input files\n" if($logger->is_debug);
            find( sub {  $berFiles->{$1} = $File::Find::name 
                             if($File::Find::name =~ /.*\.bsml/ && /polypeptide\.(\d+)\./) },
                  $opts->{'ber_input'});

        } elsif(-e $opts->{'ber_input'}) {
            open(IN, "< $opts->{'ber_input'}") or
                $logger->logdie("Unable to open ber_input file $opts->{'ber_input'} ($!)");
            my $line;
            while($line = <IN>) {
                chomp $line;
                $berFiles->{$1} = $line if($line =~ /polypeptide.(\d+)\./);
            }
            close(IN);
            
        } else {
            $logger->logdie("Option ber_input ($opts->{'ber_input'}) is not a file or a directory");
        }
    } else {
        $logger->logdie("Option ber_input is required");
    }

    #Names and nodes (NCBI file dump) names.bcp and nodes.bcp
    #(Can also use bcp dump of SYBPANDA, database prometheus, gb_names and gb_nodes).
    if($opts->{'names_dump'}) {
        $logger->logdie("$opts->{'names_dump'} does not exist") unless( -e $opts->{'names_dump'} );
        $names = $opts->{'names_dump'};
    } else {
        $logger->logdie("Option names_dump is required");
    }

    if($opts->{'nodes_dump'}) {
        $logger->logdie("$opts->{'nodes_dump'} does not exist") unless( -e $opts->{'nodes_dump'} );
        $nodes = $opts->{'nodes_dump'};
    } else {
        $logger->logdie("Option nodes_dump is required");
    }

    if($opts->{'ber_extension'}) {
        $berExtension = $opts->{'ber_extension'};
    }

    if($opts->{'base_name'}) {
        $baseName = $opts->{'base_name'};
    } else {
        $logger->logdie("Option --base_name is required");
    }

    
    ## MLDBM file stored with this structure:
    ##  $h->{hmm_id} = {
    ##                      trusted_cutoff => N,
    ##                      
    ##                 };
     
    if($opts->{'hmm_info'}) {
        tie(%hmmInfo, 'MLDBM', $opts->{'hmm_info'}, mode => 0_RDONLY )
            or $logger->logdie("Unable to tie hash to $opts->{'hmm_info'}");
    } elsif( $opts->{'hmm_info_db'} ) {
        my $file = $opts->{'hmm_info_db'};
        tie(%hmmInfo, 'MLDBM', $file, mode => O_RDONLY )
            or $logger->logdie("Unable to tie hash to $file");
    } else {
        $logger->logdie("option hmm_info_db is required");
    }

} #End checkParameters

sub collectSeqs {
    my ($twig, $seq) = @_;

    
    my $seqId = $seq->att('id');
    $logger->logdie("Unable to parse id from sequence") unless($seqId);

    my $class = $seq->att('class');
    unless($class) {
        $class = $1 if($seqId =~ m|^[^\.]+\.([^\.]+)\.|);
        $logger->logdie("Can't parse class from $seqId\n") unless($class);
    }

    #If it's the assembly sequence I don't care about it (yet).
    return if($class eq 'assembly');
    
    my $featId = $1 if($seqId =~ /^(.*)_seq/);
    $logger->logdie("Unable to parse the feature id from $seqId") unless($featId);

    my $title = $seq->att('title');
    unless($title) {
        $title = $seqId;
    }
    $featureSeqElems->{$featId}->{'title'} = $title;

    $featureSeqElems->{$featId}->{'class'} = $class;
    $featureSeqElems->{$featId}->{'seqId'} = $seqId;

    my $molecule = $seq->att('molecule');
    my $length = $seq->att('length');

    $featureSeqElems->{$featId}->{'moluecule'} = $molecule if($molecule);
    $featureSeqElems->{$featId}->{'length'} = $length if($length);

    #Gather the Seq-data-import information
    my $sdi = $seq->first_child('Seq-data-import');
    $logger->logdie("$seqId does not have Seq-data-import") unless($sdi);
    
    my $source = $sdi->att('source');
    my $format = $sdi->att('format');
    my $ident  = $sdi->att('identifier');

    $logger->logdie("$seqId does not have a complete Seq-data-import (missing either source, format or identifier)")
        unless($source && $format && $ident);

    $featureSeqElems->{$featId}->{'source'} = $source;
    $featureSeqElems->{$featId}->{'format'} = $format;
    $featureSeqElems->{$featId}->{'identifier'} = $ident;
    
}

sub definitionHandler {
    my ($twig, $definitionsElem, $feature) = @_;
    my $bsmlId2FamilyLookup = {};

    #Get the Sequence elements.
    my $seqsElem = $definitionsElem->first_child('Sequences');
    my @sequences;
    @sequences = $seqsElem->children('Sequence') if(defined($seqsElem));
    
    #We go through all the sequences first so that we can collect family information.
    foreach my $seq ( @sequences ) {

        #Parse the sequence id
        my $id = $seq->att('id');
        $logger->logdie("Unable to parse id from sequence bsml element") unless($id);
        
        #Make sure this sequence was not the input sequence
        my $link = $seq->first_child('Link');
        if($link) {
            next if($link->att('role') eq 'input_of');
        }

        #Get the title to parse out organism
        # ex { organism };
        my $title = $seq->att('title');
        $logger->logdie("Unable to parse title from sequence bsml element ($id)") unless($title);
        my $organism = $1 if($title =~ /\{([^\;]+)/g);

        #if we couldn't parse the organism, move on
        next unless($organism);
        
        #Lookup the family
        my $family = $familyLookup->lookupFamilyByName($organism);
        $bsmlId2FamilyLookup->{$id} = $family;
        
    }

    #Get the Seq-pair-alignment elements.
    my $tablesElem = $definitionsElem->first_child('Tables');
    my @spas;
    @spas = $tablesElem->children('Seq-pair-alignment') if($tablesElem);

	my $familiesPresent;

    #Go through each Seq-pair-alignment (spa).
    foreach my $spa ( @spas ) {
        my $refseqId = $spa->att('compseq');
        $logger->logdie("Could not parse refseq id out of Seq-pair-alignment") unless($refseqId);

        my $compseqId = $spa->att('refseq');
        $logger->logdie("Could not parse compseq id out of Seq-pair-alignment") unless($compseqId);
        
        my $fam;
        $fam = $bsmlId2FamilyLookup->{$refseqId};
        next if( !defined($fam) || $fam eq "");
        
        #Insert the start and stop codons into the hashref->arrayref
        my $spr = $spa->first_child('Seq-pair-run');
        $logger->logdie("Unable to retrieve Seq-pair-run element from Seq-pair-alignment element ( $refseqId )")
            unless($spr);
        my ($start, $len) = ($spr->att('refpos'), $spr->att('runlength'));

        #Adjust gene to genomic coordinates
        #Find the start of the evidence
        #Then find the extension that was done for the ber step.  We can't just subtract the $berExtension
        # because the extension could have been truncated if the gene was predicted near the begning or
        # end of the sequence.


        my ($bsmlStart,$bsmlStop) = ($feature->{'startpos'}, $feature->{'endpos'});

        my $ext = ($bsmlStart < $berExtension) ? $bsmlStart : $berExtension;

        my $evidenceStart = ($bsmlStart - $ext) + $start;
        my $evidenceStop = $evidenceStart + $len;
	
        push(@{$familiesPresent->{$fam}}, [$evidenceStart, $evidenceStop]);
    }

	#How many families are represented?
	my $phyloSpread = scalar (keys %{$familiesPresent});

	#The spread of phylogenetic families represents the strength of the evidence.
    print STDOUT "did not pass weak cutoff score\n" if($phyloSpread < $weakEvCutoff && $logger->is_debug);
	return if($phyloSpread < $weakEvCutoff); 

	my ($starts, $stops);

	#The most conserved is the median... ?
	#GO through each family and find the median start and stop site
	foreach my $family (keys %{$familiesPresent}) {
		my @sStart = sort { $b->[0] <=> $a->[0] } @{$familiesPresent->{$family}};
		$starts->{$family} = $sStart[@sStart/2]->[0];
		my @sStop = sort { $a->[1] <=> $b->[1] } @{$familiesPresent->{$family}};
		$stops->{$family} = $sStop[@sStop/2]->[1];	
	}

    print STDOUT "Found no start sites\n" if(scalar(keys %{$starts}) == 0 && $logger->is_debug);
	return if(scalar(keys %{$starts}) == 0);

	#If we have less than five conserved starts, just set the weak evidence. 
	#Otherwise, set both strong and weak.
	if(scalar(keys %{$starts}) < 5) {
		my $ev_start = [sort {$a <=> $b} values %{$starts}]->[scalar(keys %{$starts}) - 1];
		my $ev_stop = [sort {$b <=> $a} values %{$stops}]->[scalar(keys %{$stops}) - 1];
		push(@{$feature->{'evidence'}}, {
			'type' => 'BER',
			'accession' => 'BER',
			'start'     => $ev_start,
			'stop'      => $ev_stop,
			'quality'   => $phyloSpread });
		
	} else {
		my $ev_start = [sort {$a <=> $b} values %{$starts}]->[4];
		my $ev_stop = [sort {$b <=> $a} values %{$stops}]->[4];
		push(@{$feature->{'evidence'}}, {
			'type' => 'BER',
			'accession' => 'BER',
			'start'     => $ev_start,
			'stop'      => $ev_stop,
			'quality'   => 5});

		$ev_start = [sort {$a <=> $b} values %{$starts}]->[4];
		$ev_stop = [sort {$b <=> $a} values %{$stops}]->[4];
		push(@{$feature->{'evidence'}}, {
			'type' => 'BER',
			'accession' => 'BER',
			'start'     => $ev_start,
			'stop'      => $ev_stop,
			'quality'   => $phyloSpread});
	} 


} #End definitionHandler

#Must be going in the same direction. Maybe.
sub findNewStart {
    my ($upstream, $downstream) = @_;

    $logger->logdie("downstream feat (".$downstream->{'old_id'}.") does not have any clear starts")
        unless($downstream->{'clearStarts'});

    my ($ufEnd5, $ufEnd3) = ($upstream->{'complement'}) ? ($upstream->{'endpos'}, $upstream->{'startpos'}) :
        ($upstream->{'startpos'}, $upstream->{'endpos'});

    my ($dfEnd5, $dfEnd3) = ($downstream->{'complement'}) ? ($downstream->{'endpos'}, $downstream->{'startpos'}) :
        ($downstream->{'startpos'}, $downstream->{'endpos'});

    #Find the reading frame
    my $readingFrame = $dfEnd3 % 3; 

    #If we find a start that doesn't overlap with any of the upstream feature's evidence, that's a winner.
    #Return it immediately.  If somehow we make it all the way through, return negative 1.

  STARTS:
    foreach my $possibleStart ( @{$downstream->{'clearStarts'}} ) {
        next unless($possibleStart % 3 == $readingFrame);
        foreach my $upEvidence ( @{$upstream->{'evidence'}} ) {

            my $evEnd5 = ($upstream->{'complement'}) ? $upEvidence->{'stop'} : $upEvidence->{'start'};
            
            #Determine if these coordinates overlap.
            my @sorted = sort {$a<=>$b} ($possibleStart, $dfEnd3, 
                               $evEnd5, $ufEnd3);
            my $overlap = (abs($possibleStart - $dfEnd3) + 
                        abs($evEnd5 - $ufEnd3)) -
                        (abs($sorted[0] - $sorted[3]));

            next STARTS if($overlap > 0);

        }

        return $possibleStart;
    }

    return undef;
    
    
}#end findNewStartSite

#Will find all the start sites in a given fasta file with identifier.
sub findStartSites {
    my ($identifier, $fastaFile) = @_;
    my $retval = {};

    #Data structure of start Codons
    my $codons ={
        'forward' => {
            'ATG' => 'M',
            'GTG' => 'V',
            'TTG' => 'L' },
        'reverse' => {
            'CAT' => 'M',
            'CAC' => 'V',
            'CAA' => 'L' }
    };
    
    my $found = 0; #Flag set to 1 after encountering the identifier
    my $lengthOffset = 0; #Parse each line as we find it, this holds how much sequence we've seen

    open(IN, "< $fastaFile") or $logger->logdie("Unable to open $fastaFile ($!)");

    my $prevTwoNucs = "";
    $lengths->{$identifier} = 0;

    while(<IN>) {
        chomp;
        if(/^>$identifier/) {
            $found = 1;
        } elsif(/^>/ && $found) {
            $found = 0;
        } elsif($found) {
            $lengths->{$identifier}+=length($_);
            foreach my $strand ( keys %{$codons} ) {
                foreach my $codon ( keys %{$codons->{$strand}} ) {
                    my $pos = -1; #Make sure we search from index 0
                    while(($pos = index($_, $codon, $pos+1)) > -1) {
                        my $frame = ( ($pos+$lengthOffset) % 3 ) + 1;
                        $pos+=3 if($strand eq 'reverse');
                        push(@{$retval->{$strand}->{$frame}}, $pos+$lengthOffset);
                    }
                    
                    #Do any start/stops span the lines?
                    next if($prevTwoNucs eq "");
                    my $check = $prevTwoNucs. substr($_, 0, 2);
                    $pos = -1; 
                    while(($pos = index($check, $codon, $pos+1)) > -1) {
                        my $frame = ( ($pos+$lengthOffset-2) % 3) + 1;
                        $pos += 3 if($strand eq 'reverse');
                        push(@{$retval->{$strand}->{$frame}}, $pos-2+$lengthOffset);
                    }
                }
            }
            $prevTwoNucs = substr($_, -2, 2);
            $lengthOffset+=length($_);
            
        }
    }
    return $retval;
    
} #End findStartSites

#I don't know what this does.
sub getClearStartSites {
    my $feature = shift;

    my @clearStarts;

    my ($fEnd5, $fEnd3) = ($feature->{'complement'}) ? ($feature->{'endpos'}, $feature->{'startpos'}) :
        ($feature->{'startpos'}, $feature->{'endpos'});

    #Now find the no further downstream point.
    my $nfd = $fEnd3;
    foreach my $evidence ( @{$feature->{'evidence'}} ) {
        next unless($evidence->{'quality'} >= $strongEvCutoff);
        my $evEnd5 = ($feature->{'complement'}) ? $evidence->{'stop'} : $evidence->{'start'};
        $nfd = $evEnd5 if(abs($evEnd5 - $fEnd3) > abs($nfd - $fEnd3));
    }

    #Now find all the start sites between the feature start and the nfd.
    foreach my $strand ( keys %{$startSites->{$feature->{'parentseq'}}}) {
        #Make sure we only look for starts on the correct strand.
        next unless( ($feature->{'complement'} == 1 && $strand eq 'reverse') ||
                     ($feature->{'complement'} == 0 && $strand eq 'forward'));
        foreach my $frame ( keys %{$startSites->{$feature->{'parentseq'}}->{$strand}} ) {

            foreach my $start( @{$startSites->{$feature->{'parentseq'}}->{$strand}->{$frame}} ) {
                if($start >= $fEnd5 && $start <= $nfd && !($feature->{'complement'}) && $start <= $fEnd3 ) {
                    push (@clearStarts, $start);
                } elsif($start <= $fEnd5 && $start >= $nfd &&  $feature->{'complement'} && $start >= $fEnd3) {
                    push @clearStarts, $start;
                }
            }

        }
    }
    
    return \@clearStarts;
    
}#end getClearStartSites

#Will retrieve evidence for a feature object.
sub getEvidence {
    my $feature = shift;

    #First look for ber and hmm file.
    my $polypeptideId = $genes->{$feature->{'old_id'}}->{'polypeptide'};
    $logger->logdie("Unable to find polypeptide id from ".$feature->{'old_id'})
        unless($polypeptideId);
    my $berFile = $berFiles->{$1} if($polypeptideId =~ /polypeptide\.(\d+?)\./);
    $logger->logdie("Unable to find ber analysis for polypeptide $polypeptideId (".
                    $feature->{'old_id'}.")") unless($berFile);
    my $hmmFile = $hmmFiles->{$1} if($polypeptideId =~ /polypeptide\.(\d+?)\./);
    $logger->logdie("Unable to find hmm analysis for polypeptide $polypeptideId (".
                    $feature->{'old_id'}.")") unless($hmmFile);

    #Parse BER file
    &analyzeBerEvidence($berFile, $feature);

    #Parse HMM File
    &analyzeHmmEvidence($hmmFile, $feature);
}#end getEvidence

#Twig function for analyzeHmmEvidence.  Will collect information about hmm hits and store
#them in the array ref found in $feature->{'evidence'}.  Consults the tied has hmmInfo
#for trusted cutoff scores.
sub handleHmmSPA {
    my ($twig, $spaElem, $feature) = @_;

    #Get the seq pair runs
    my @sprs = $spaElem->children('Seq-pair-run');

    $logger->logdie("SPA (".$spaElem->att('compseq').") did not have any spr's") 
    
    ## you would only do this if your HMM lib had duplicate IDs and the BSML wasn't formed correctly
    ##next if(@sprs == 0);

    foreach my $spr ( @sprs ) {
        
        my ($start, $length, $compseq) = ($spr->att('refpos'),$spr->att('runlength'), $spaElem->att('compseq'));
        
        ## if compseq ID doesn't exist and starts with an underscore, this may be because it really started with
        ##  an integer and BSML couldn't allow that within an id attribute.  Check for this.
        if ( ! exists $hmmInfo{$compseq} ) {

            if ( $compseq =~ /^_(\d.+)/ && exists $hmmInfo{$1} ) {
                $logger->warn("replaced ID of $compseq with $1 to correct for the BSML id attribute limitations");
                $compseq = $1;
            } else {
                $logger->logdie("Could not find hmm $compseq in hmm_info");
            }
        }
        
        unless($hmmInfo{$compseq}->{'trusted_cutoff'}) {
            print STDOUT Dumper($hmmInfo{$compseq});
        }
            
        $logger->logdie("Hmm does not have trusted cutoff ($compseq)") 
            unless( $hmmInfo{$compseq}->{'trusted_cutoff'} );
        $logger->logdie("current spr does not have a run length ($compseq)") unless($spr->att('runlength'));
        next if($hmmInfo{$compseq}->{'trusted_cutoff'} >= $spr->att('runlength') );

        my $featureId = $feature->{'old_id'};
        #These are gene coordinates, lets change them to genomic coordinates
        $start = $feature->{'startpos'};

        $logger->logdie("Couldn't parse start position for $featureId")
            unless($start);
        
        #Store the necesary information from the hit.
        push(@{$feature->{'evidence'}}, {
            'type'    => 'HMM',
            'start'    => $start,
            'stop'    => $start+$length,
            'quality' => 50, #Arbitrarily set quality to 50.
        });  
    }
}# end handleHmmSPA

#We have already determined that ther exists an overlap between these two feature elements.
# Conditions:
#   - If the coordinates are the same, remove one (the second one passed in).
#   - If they features overlap in the same frame, remove the shorter one.
#   - Check evidence
#      - If neither have evidence, don't remove either (FUTURE: check with metagene results)
#      - If genes are in opposite directions, try to move the end5s clear of the other's evidence.
#           ex. (Represent feature evidence [BER and hmmpfam])
#            3'             5'                      3'            5'
#             |-------------|                        |------------|
#                                          ---->
#                       |--------------|                            |-------------|
#                      5'              3'                           5'            3'      
#
#      - If genes are going in same direction, move the end5 of downstream gene away from upstream genes
#        evidence.  Make sure that the gene still makes the length cutoff.
sub handleOverlap {
    my ($feat1, $feat2) = @_;

    #Check to see if coordinates are the same
    if($feat1->{'startpos'} == $feat2->{'startpos'} && $feat1->{'endpos'} == $feat2->{'endpos'}) {
        &removeFeature($feat2);
        print STDOUT "Removed because the coordinates are the same\n" if($logger->is_debug);
        return;
    }
    
    #If overlapping features are in the same frame, get rid of the shorter
    if($feat1->{'startpos'} % 3 == $feat2->{'startpos'} % 3 &&
       $feat1->{'endpos'} % 3   == $feat2->{'endpos'}   % 3 ) {
        my $removeMe = (abs($feat1->{'startpos'} - $feat1->{'endpos'})  <
                        abs($feat2->{'startpos'} - $feat2->{'endpos'})) ? 
                        $feat1 : $feat2;
        &removeFeature($removeMe);
        print STDOUT "Removed because in same frame\n" if($logger->is_debug);
        return;
    }

    #If the genes are going in the same direction 
    if($feat1->{'complement'} == $feat2->{'complement'}) {

        #Put this in terms of upstream and downstream
        my ($upstream, $downstream) = ($feat1->{'startpos'} < $feat2->{'startpos'}) ? 
            ($feat1, $feat2) : ($feat2, $feat1);   

        #Get evidence for these two features
        &getEvidence($upstream);
        &getEvidence($downstream);

        #Start sites clear of the most downstream evidence.
        $downstream->{'clearStarts'} = &getClearStartSites($downstream);

        #Now try and move the current start site. 
        #move end5 of downstream gene within cutoff tolerance of upstream end3.
        my $newStart = &findNewStart($upstream, $downstream);

        my ($dfEnd5, $dfEnd3) = $downstream->{'complement'} ? ($downstream->{'endpos'},$downstream->{'startpos'}) :
            ($downstream->{'startpos'},$downstream->{'endpos'});
       
        if( defined($newStart) && $newStart == $dfEnd3 && scalar @{$upstream->{'evidence'}} == 0 ) {
            print STDOUT "removing ".$upstream->{'old_id'}." because it has no evidence\n" if($logger->is_debug);
            &removeFeature($upstream);
        } elsif( defined($newStart) && abs($newStart-$dfEnd3) < $minGeneLength ) {
            print STDOUT "removing ".$downstream->{'old_id'}." because the new start makes it too short\n" if($logger->is_debug);
            &removeFeature($downstream);
        } elsif( defined($newStart) ) {
            print STDOUT "found a new start for ".$downstream->{'old_id'}.". Changing from $dfEnd5 to $newStart\n" if($logger->is_debug);
            my $oldStart = $downstream->{'complement'} ? $downstream->{'endpos'} : $downstream->{'startpos'};
            &changeFeature($downstream, $newStart) unless($newStart == $oldStart);
        } else {
            print STDOUT "newStart was undefined\n" if($logger->is_debug);
        }
        

    #If the features are going in opposite directions (different strands).
    } else {
        print STDOUT "They are not on the same strand.  Calling evidence\n" if($logger->is_debug);
        &getEvidence( $feat1 );
        &getEvidence( $feat2 );

        $feat1->{'clearStarts'} = &getClearStartSites($feat1);
        $feat2->{'clearStarts'} = &getClearStartSites($feat2);

        my $f1_start = &findNewStart($feat2, $feat1);
        my $f2_start = &findNewStart($feat1, $feat2);

        #We know that they are opposite.
        my ($f1End5, $f1End3, $f2End5, $f2End3);
        if($feat1->{'complement'}) {
            ($f1End5, $f1End3) = ('endpos', 'startpos');
            ($f2End5, $f2End3) = ('startpos', 'endpos');
        } else {
            ($f2End5, $f2End3) = ('endpos', 'startpos');
            ($f1End5, $f1End3) = ('startpos', 'endpos');
        }

        if(!$f1_start) {
            if(scalar @{$feat1->{'evidence'}} == 0) {
                &removeFeature($feat1);
            } else {
                print STDOUT "No new start site but there is evidence\n" if($logger->is_debug);
            }
        } elsif( abs($f1_start - $feat1->{$f1End3}) < $minGeneLength ) {
            &removeFeature($feat1);
        } elsif($feat1->{$f1End5} != $f1_start) {
            &changeFeature($feat1, $f1_start);
        } else {
            print STDOUT "New start was the same.  Do nothing\n" if($logger->is_debug);
        }

        if(!$f2_start) {
            if(scalar @{$feat2->{'evidence'}} == 0) {
                print STDOUT "feat 2 removed because there was no evidence and no new start site.\n" if($logger->is_debug);
                &removeFeature($feat2);
            } else {
                print STDOUT "No new start site but there is evidence\n"if($logger->is_debug) ;
            }
        } elsif( abs($f2_start - $feat2->{$f2End3}) < $minGeneLength ) {
            print STDOUT "feat 2 removed because the new start site did not meet the minimum gene length\n" if($logger->is_debug);
            &removeFeature($feat2);
        } elsif($feat2->{$f2End5} != $f2_start) {
            print STDOUT "feat 2 is being changed\n" if($logger->is_debug);
            &changeFeature($feat2, $f2_start);
        } else {
            print STDOUT "New start was the same.  Do nothing\n" if($logger->is_debug);
        }

            
    }

    
}

sub handleSequence {
    my ($twig, $seqElem) = @_;

    my $intervalTree = new IntervalTree;

    my $id = $seqElem->att('id');
    $logger->logdie("Can't parse id from input bsml file") unless($id);

    my $class = $seqElem->att('class');
    unless($class) {
        $class = $1 if($id =~ /^[^\.]+\.([^\.]+)\./);
        $logger->logdie("Can't parse class out of $id") unless($class);
    }
    return unless($class eq 'assembly');
       

    print STDOUT "handling $id\n" if($logger->is_debug);

    #We need to find all the start sites for this sequence to use later on (when finding no further downstream
    #sites, etc.  
    my $seqDataImport = $seqElem->first_child('Seq-data-import');
    $logger->logdie("Unable to parse the Seq-data-import element for sequence $id")
        unless($seqDataImport);
    my $identifier = $seqDataImport->att('identifier');
    $logger->logdie("Unable to parse identifier attribute from Seq-data-import element inside Sequence $id")
        unless($identifier);
    my $fsaFile = $seqDataImport->att('source');
    $logger->logdie("Unable to parse source from Seq-data-import element inside Sequence $id")
        unless($fsaFile);
    $startSites->{$id} = &findStartSites( $identifier, $fsaFile );
    $sequenceSource->{'fsaFile'} = $fsaFile;
    $sequenceSource->{'id'} = $identifier;
    

    #Find overlaps.  First go through all the bsml cds feature elements, adding
    #the start and stop coordinates to the interval tree.  Then go through the
    #features again

    #Get feature elements
    my $fTables = $seqElem->first_child('Feature-tables');
    my @bsmlFeatures = ();
    if($fTables) {
        @bsmlFeatures = $fTables->first_child('Feature-table')->children('Feature');
    }

    my $total = scalar @bsmlFeatures;
    my $count = 0;

    #Loop through feature elements
    foreach my $feature ( @bsmlFeatures ) {
        
        print STDOUT "\r$count/$total " unless($logger->is_debug);
        my $percentage = int(($count/$total)*10000)/100;
        print STDOUT "$percentage %" unless($logger->is_debug);
        $count++;
        
        
        #Get the feature id, die if we don't parse one
        my $featId = $feature->att('id');
        $logger->logdie("Unable to parse id from feature element (in $id)")
            unless($featId);


        #Make a new feature object and store for later use.
        $features->{$featId} = new Feature();
        $features->{$featId}->init($feature);
        $features->{$featId}->{'parentseq'} = $id;

        #Add the coordinates to the interval tree.
        $intervalTree->addInterval( $featId, $features->{$featId}->{'startpos'}, $features->{$featId}->{'endpos'} )
            if($features->{$featId}->{'class'} eq 'CDS');

    }

    print STDOUT "\nBuilding the tree\n";
    $intervalTree->buildTree;

    #Go through all the feature groups to find out which features are part of which genes.
    my @featureGroups;
    @featureGroups = $fTables->children('Feature-group') if($fTables);

    #Make a gene out of each featureGroup and index them by CDS id for easy access later.
    foreach my $featureGroup ( @featureGroups ) {
        my $tmpGene = new Gene;
        $tmpGene->init($featureGroup);
        $logger->logdie("Unable to find CDS feature id in gene") unless($tmpGene->{'CDS'});
        $genes->{$tmpGene->{'CDS'}} = $tmpGene;
    }


    $count = 0;
    #Now look for overlaps
    
    print STDOUT "Starting to look for overlaps\n";
    foreach my $feature ( values %{$features} ) {
        print STDOUT "\r$count/$total ";
        my $percentage = int(($count/$total)*10000)/100;
        print STDOUT "$percentage %";
        $count++;
        
        next unless($feature->{'class'} eq 'CDS');

        my @overlaps = $intervalTree->searchInterval( $feature->{'startpos'},
                                                      $feature->{'endpos'} );


        #Make sure that we haven't seen this overlap before.
        my @filteredOverlaps = grep ( &isNewOverlap($feature->{'old_id'}, $_->[2]),
                                     @overlaps);

        
        
        foreach my $overlap ( @filteredOverlaps ) {
            $feature->{'overlaps'}->{ $overlap->[2] } = 1;
            &handleOverlap( $feature, $features->{$overlap->[2]} );
        }   

    }

    #Print the other sequences
    print STDOUT "Found ".(scalar (keys %{$featureSeqElems}))." non-assembly sequences\n";
    
    my $seqParent = $seqElem->parent();
    while( my ($featId, $sequenceInfo) = each(%{$featureSeqElems}) ) {

        my $title = $sequenceInfo->{'title'} ? $sequenceInfo->{'title'} : $featId;

         my $intLoc = $seqParent->new( 'Seq-data-import', 
                                       { 'source' => $sequenceInfo->{'source'}, 
                                         'format' => $sequenceInfo->{'format'},
                                         'identifier' => $sequenceInfo->{'identifier'},
                                     }, ('') );

         my $link = $seqParent->new( 'Link', { 'rel' => 'analysis', 
                                               'href' => '#auto_gene_curation_analysis', 
                                               'role' => 'input_of' }, ('') );
        
        unless($sequenceInfo->{'molecule'}) {
            $sequenceInfo->{'molecule'} = 'dna';
            $sequenceInfo->{'molecule'} = 'aa' if($sequenceInfo->{'class'} eq 'polypeptide');
        }

        

        my $newSeq = $seqParent->new( 'Sequence', 
                                      { 'id' => $sequenceInfo->{'seqId'},
                                        'class' => $sequenceInfo->{'class'},
                                        'title' => $title,
                                        'molecule' => $sequenceInfo->{'molecule'} },
                                      ($intLoc, $link));
        $newSeq->print;
    
    }

    $seqElem->print; 
}

#Pass in two feature ids.  
#Returns 1 if we've seen this overlap or 0 otherwise.
#Checks the 'overlaps' hashref of feature elements.
#Also filters out identity hits.
sub isNewOverlap {
    my ($featId1, $featId2) = @_;
    my $retval = 1; #assume it's new

    #Check to see if it's already in the overlap hash.
    $retval = 0 if($features->{$featId1}->{'overlaps'}->{$featId2});
    $retval = 0 if($features->{$featId2}->{'overlaps'}->{$featId1});

    #Filter out identity hits.
    $retval = 0 if($featId1 eq $featId2);

    return $retval;
    
}
sub makeFsaFile {
    my ($source, $difference, $ident, $newId) = @_;
    my $flag = 0;

    #Open the input source file.
    open(IN, "< $source") or $logger->logdie("Unable to open $source ($!)");

    my $seq = "";

    #Search for the correct identifier in the fasta file.
    while(my $line = <IN>) {
        chomp $line;
        if($line =~ /^>$ident/) {
            $flag = 1;
        } elsif($flag && $line !~ /^\s+$/ && $line !~ /^>/) {
            $seq .= $line;
        } elsif($flag && $line =~ /^>/) {
            last;
        }
    }
    close(IN);

    if(length($seq) < ($difference - length($seq))) {
        print STDOUT "$newId\n";
        print STDOUT "$source\n";
        print STDOUT "diff: $difference\n";
        print STDOUT length($seq)."\n";
        die("Bad\n");
    }

    #$difference contains the number of nucleotides we are chopping off from the beginning.  This is under the
    #assumption that the start site will create a smaller protein (because we are only handling overlaps, 
    #changing a start site that results in the increase in size of a gene will only increase the overlap). 
    #Therefore we just need to chop a few nucleotides off the beginning of the sequence.
    my $newSeq = substr($seq, ($difference - length($seq)), (length($seq) - $difference));

    #Open the new output file.
    my $newOutFsa = "$outDirectory/$newId.fsa";
    open(OUT, "> $newOutFsa") or
        $logger->logdie("Unable to open $newOutFsa for writing ($!)");

    print OUT ">$newId\n";
    
    my $lineLength = 60;
    my $buffer;

    for (my $i = 0; $i < length($seq); $i+=$lineLength) {
        print OUT substr($seq, $i, $lineLength);
        print OUT "\n";
    }

    close(OUT);

    return $newOutFsa;

}

sub printChangedFeatures {
    print STDOUT "building the new bsml output file\n";

    my $out = "$outDirectory/$baseName.new.bsml";

    my $bsmlObj = new BSML::BsmlBuilder;

    my $asmblId = $sequenceSource->{'id'};
    my $length = $lengths->{ $asmblId };
    my $seq = $bsmlObj->createAndAddSequence( $asmblId, $asmblId, $length, 'dna', 'assembly' );
    $seq->addBsmlLink( 'analysis', '#auto_gene_curation_analysis', 'input_of' );
    my $sdi = $bsmlObj->createAndAddSeqDataImport( $seq, 'fasta', $sequenceSource->{'fsaFile'}, '', $asmblId);
    my $fTable = $bsmlObj->createAndAddFeatureTable( $seq );

    #Cycle through the features
    foreach my $featureId ( keys %{$changedFeatures} ) {
        my $feature = $features->{$featureId};
        next unless($feature->{'class'} eq 'CDS');

        #Create a feature group
        my $featureGroup = $bsmlObj->createAndAddFeatureGroup( $seq, '', $feature->{'new_id'} );

        foreach my $type ( keys %{$genes->{$feature->{'old_id'}}} ) {
            next if($type eq 'fgm' || $type eq 'feature-group');
	    my $featureGroupMemberId = $genes->{$feature->{'old_id'}}->{$type};
            my $other = $features->{$featureGroupMemberId};
            my $fgmId = $other->{'new_id'};
            $logger->logdie("Unable to find the feature object for $featureGroupMemberId") unless($other);

            my $otherBsmlFeature = &addFeature( $other, $fTable, $bsmlObj );

            my $otherFeatType = $1 if($featureGroupMemberId =~ /^[^\.]+\.([^\.]*)\./);
            my $other_fgm = $bsmlObj->createAndAddFeatureGroupMember( $featureGroup, 
                                                                      $fgmId, 
                                                                      $otherFeatType );

        }
    }

    $bsmlObj->createAndAddAnalysis();
    print STDOUT "Wrote to $out\n" if($logger->is_debug);
    $bsmlObj->write($out);
    

}

#Deletes the bsml_feature_object (so it will not remain in the bsml file)
#Also removes the feature from the $features hash (which should contain all current features).
#Also removes it from the genes object.
#Also removes the sequence element if its there.
sub removeFeature {
    my $feature = shift;

    foreach my $type ( keys %{$genes->{$feature->{'old_id'}}} ) {
        next if($type eq 'fgm' || $type eq 'feature-group');
        my $featId = $genes->{$feature->{'old_id'}}->{$type};
        $features->{$featId}->{'bsml_feature_object'}->delete 
            if($features->{$featId}->{'bsml_feature_object'});
        $deletedFeatures->{$featId} = 1;
        $logger->logdie("No fgm. $type $featId") unless($genes->{$feature->{'old_id'}}->{'fgm'}->{$type});

    }

    #remove the sequence elements (for polypeptides and cds')
    delete($featureSeqElems->{$feature->{'old_id'}}) if(exists($featureSeqElems->{$feature->{'old_id'}}));
    delete($featureSeqElems->{ $genes->{$feature->{'old_id'}}->{'polypeptide'} } ) 
        if( exists( $featureSeqElems->{ $genes->{ $feature->{'old_id'} }->{'polypeptide'} } ) );

    $genes->{$feature->{'old_id'}}->{'feature-group'}->delete;

    
}

#Just in case I want to print out the documentation.
sub _pod {
    &pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDOUT} );
}
######################## END SUBROUTINES #################################

######################## PACKAGES ########################################

package Feature;

sub new {
    my $class = shift;
    my %params = @_;
    my $self = \%params;

    return( bless( $self, $class ) );
}

sub init {
    my ($self, $bsmlFeature) = @_;
    
    $self->{'old_id'} = $bsmlFeature->att('id');
    $self->{'class'} = $bsmlFeature->att('class');
    $self->{'interval-loc'} = $bsmlFeature->first_child('Interval-loc');

    ($self->{'startpos'}, $self->{'endpos'}, $self->{'complement'}) = 
        ($self->{'interval-loc'}->att('startpos'),
         $self->{'interval-loc'}->att('endpos'), 
         $self->{'interval-loc'}->att('complement'));

    $self->{'bsml_feature_object'} = $bsmlFeature;   
}

package Gene;

sub new {
    my $class = shift;
    my %params = @_;
    my $self = {};

    return( bless( $self, $class ) );
}

sub init {
    my ($self, $featureGroup) = @_;
    my @fgms = $featureGroup->children('Feature-group-member');

    foreach  my $fgm ( @fgms ) {
        my $type = $fgm->att('feature-type');
        $logger->logdie("Unable to open parse type from feature group member") unless($type);
        $self->{$type} = $fgm->att('featref');
        $self->{'fgm'}->{$type} = $fgm;
        
    }
    
    $self->{'feature-group'} = $featureGroup;
}

package FamilyLookup;

sub new {
    my $class = shift;
    my %params = @_;

    return undef unless(exists($params{'names'}) && exists($params{'nodes'}));

    my $self = {
        'nodes' => $params{'nodes'},
        'names' => $params{'names'}
    };

    return( bless( $self, $class ) );
}

#Reads in BCP dumps of tables 'names' and 'nodes'
#from NCBI's taxonomy database (or TIGR's prometheus on SYBPANDA).
#
# Column headers, delimited by white space and '|'.
#
# nodes.dmp
# ---------
# This file represents taxonomy nodes. The description for each node includes 
# the following fields:
#   tax_id					           -- node id in GenBank taxonomy database
#  	parent tax_id				       -- parent node id in GenBank taxonomy database
#  	rank					           -- rank of this node (superkingdom, kingdom, ...) 
#  	embl code				           -- locus-name prefix; not unique
#  	division id				           -- see division.dmp file
#  	inherited div flag  (1 or 0)	   -- 1 if node inherits division from parent
#  	genetic code id				       -- see gencode.dmp file
#  	inherited GC  flag  (1 or 0)	   -- 1 if node inherits genetic code from parent
#  	mitochondrial genetic code id	   -- see gencode.dmp file
#  	inherited MGC flag  (1 or 0)	   -- 1 if node inherits mitochondrial gencode from parent
#  	GenBank hidden flag (1 or 0)       -- 1 if name is suppressed in GenBank entry lineage
#  	hidden subtree root flag (1 or 0)  -- 1 if this subtree has no sequence data yet
#  	comments				           -- free-text comments and citations
#
# names.dmp
# ---------
# Taxonomy names file has these fields:
# 	tax_id					-- the id of node associated with this name
# 	name_txt				-- name itself
# 	unique name				-- the unique variant of this name if name not unique
# 	name class				-- (synonym, common name, ...)
#
# (ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_readme.txt)
#
#Function returns 1 on success and an error string on failure
sub preprocessTaxonomyInformation {
    my $self = shift;

    #Here we will create a hash linking species to their related families using above 
    #described database dump files.
    my ($names, $nodes) = ($self->{'names'}, $self->{'nodes'});

    open(NAMES, "< $names") or return "Could not open $names ($!)";
    
    #First go through the names file and record all the tax_ids that contain a scientific name.
    while(<NAMES>) {
        
        #Column separator is '|'
        my @cols = split(/\s*\|\s*/,$_);
        next unless($cols[3] eq 'scientific name');

        my $tmpHash = { 'id' => $cols[0],
                        'name' => $cols[1] };
        
        $self->{'byId'}->{$cols[0]} = $tmpHash;
        $self->{'byName'}->{$cols[1]} = $tmpHash;

    }

    close(NAMES);


    #Now go through nodes.  We only need to look at the rows that are scientific names
    #(the information retrieved from the names table).
    open(NODES, "< $nodes") or return "Could not open $nodes ($!)";

    while(<NODES>) {
        my @cols = split(/\s*\|\s*/,$_);
        next unless(exists($self->{'byId'}->{$cols[0]}));
        $self->{'byId'}->{$cols[0]}->{'parent'} = $cols[1];
        $self->{'byId'}->{$cols[0]}->{'rank'} = $cols[2];
    }

    close(NODES);

    return 1;

}

sub lookupFamilyByName {
    my ($self,$name) = @_;
    my $retval;

    unless(exists($self->{'byName'}->{$name})) {
        return "";
    }

    if($self->{'byName'}->{$name}->{'rank'} eq 'family') {
        $retval = $name;
    } elsif( $name eq 'root') {
        $retval = "";
    } else {
        $retval = $self->lookupFamilyById($self->{'byName'}->{$name}->{'parent'});
    }

    return $retval;
    
}

sub lookupFamilyById {
    my ($self, $id) = @_;
    my $retval;

    unless(exists($self->{'byId'}->{$id})) {
        return "";
    }

    if($self->{'byId'}->{$id}->{'rank'} eq 'family') {
        $retval = $self->{'byId'}->{$id}->{'name'};
    } elsif( $self->{'byId'}->{$id}->{'name'} eq 'root') {
        $retval = "";
    } else {
        $retval = $self->lookupFamilyById($self->{'byId'}->{$id}->{'parent'});
    }

    return $retval;
}
