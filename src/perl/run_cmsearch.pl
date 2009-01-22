#!/usr/bin/perl

=head1 NAME

run_cmsearch.pl - Takes hmmpfam files as input and runs the hits against
    cmsearch.  Optimization for cmsearch.

=head1 SYNOPSIS

USAGE: run_cmsearch.pl 
              --input_list=/path/hmmpfam/results.raw.list
              --input_file=/path/hmmpfam/result.raw
              --tmp_dir=/path/to/tmpdir/
              --output_file=/path/to/infernal.raw
              --flanking_seq=50
              --hmm_cm_table=/path/to/some/file.table
              --cmsearch_bin=/path/to/cmsearch
              --cm_dir=/dir/with/covariance/models/
          [   --other_opts=cmsearch options
              --log=/path/to/some/file.log
              --debug= > 2 for verbose
          ]

=head1 OPTIONS

B<--input_list,-l>
    A list of hmmpfam raw output files to be formatted and run with cmsearch

B<--input_file,-i>
    An input file (hmmpfam) raw output

B<--tmp_dir,-t>
    Directory to write temporary files.  Will be cleaned up (not yet).

B<--output_file,-o>
    The output file for cmsearch results to go into.  Will concate all results to this file.

B<--flanking_seq,-f>
    The number of nucleotides on either side of the hmmpfam hit to parse out of database for 
    cmsearch run. Default: 50.

B<--hmm_cm_table,-c>
    File containing a lookup for covariance models given an hmm model.  See input section of perldoc
    more specific details on this file.

B<--cmsearch_bin,-b>
    Path to the cmsearch binary. If not it will be assumed that the binary is in the PATH.

B<--cm_dir,-m>
    If you don't feel like making a hmm_cm_table, you can always just put the directory where all
    the covariance models are.  The program will try and guess which covariance model is related to
    each hmm.  Needless to say, hmm_cm_table is a better option.  Turns out, computers are bad
    guessers.

B<--other_opts,-e>
    Other options to be passed into the cmsearch program.  
    *Note: The -W option is used by this script and therefore will be parsed out of the other_opts
    string prior to running cmsearch.  This program uses the length of the sequence found by the
    hmmpfam run (plus extraSeq) as the window size.

B<--debug> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

The program infernal is used to search covariance models against biological sequences
(see infernal userguide: oryx.ulb.ac.be/infernal.pdf).  More specifically infernal uses
profiles of RNA secondary structure for sequence analysis.  The program is computationally
intensive when used with large genomic sequences.  

Run_cmsearch.pl will take in a set of hmmpfam results as an initial screening to infernal.  
RNA HMMs can quickly identify regions where RNA is located, and these sections can then be further
refined by running this section through cmsearch (infernal).  With a set of HMMs and CMs created from
the same alignments, infernal can be used in a high-throughput manner.

=head1  INPUT

The main input for run_cmsearch is either a list of hmmpfam raw results or one hmmpfam raw result.
The sequences are then parsed from these hmmpfam results and the database queryed.  The sequence
identifier is taken from the name of the file (which is reliable in the ergatis naming scheme, but
probably could be made more general in the future).  The name of the HMMs are also parsed from these
files.  Important to note that the HMM names are not the same as the file names (since hmmpfam can
use a multi-hmm file for searching).  This becomes important in the creation of the hmm_cm_table.

If a hit was found with an HMM file, the related CM file will be used to search that portion of the
sequence.  To do this, the program must know which HMMs relate to which CMs.  There are two ways
to provide this information.  This by providing a tab-delimited file with HMM name (not file name,
but the name found in the HMM header and the one used in hmmpfam output) and CM (actual path to file)
on each line.  For example:

    RF00001.HMM     /usr/local/db/RFAM/CMs/RF00001.cm
    RF00002.HMM     /usr/local/db/RFAM/CMs/RF00002.cm
    ...

The other option is to provide the program with the directory where all the covariance models are
stored.  The program will then look in that directory for a covariance model that contains the first
word of the HMM name (not file name, but actual id name found in hmmpfam results).  

=head1  OUTPUT

The program will create individual fasta files for each hmmpfam hit provided in the input and is 
cleaned up later.  All the cmsearch results (all the stdout from the program run) is written to
output_file.  See infernal userguide: oryx.ulb.ac.be/infernal.pdf for more information on output format.

=head1  CONTACT

    Kevin Galens
    kgalens@tigr.org

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Ergatis::Logger;
use Data::Dumper;
use IPC::Open3;

my %options = ();
my $results = GetOptions (\%options, 
                          'input_list|l=s',
                          'input_file|i=s',
                          'tmp_dir|t=s',
                          'output_file|o=s',
                          'sequence_list|s=s',
                          'flanking_seq|f=i',
                          'hmm_cm_table|c=s',
                          'cmsearch_bin|b=s',
                          'other_opts|e=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || &_pod();

## display documentation
&_pod if( $options{'help'} );

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();


############ GLOBALS AND CONSTANTS #############
my $PROG_NAME = 'cmsearch';

my @files;
my $input;
my $hmmCmTable;
my %hmmCmLookup;
my $cmDir;
my $outputDir;
my $outputFile;
my $debug;
my @seqList;
my $extraSeqLen = 50;
my $other_opts = "";
my @dirsMade;
################################################

## make sure everything passed was peachy
&check_parameters(\%options);


################## MAIN #########################
&createHmmCmLookup($hmmCmTable) if(defined($hmmCmTable));

#All the input files if it's a list.
foreach my $file(@files) {
    my $stats = &processHmmpfamFile($file);
    $stats = &getSequencesAndPrint($stats);

    #Make sure that we make an output file even if there aren't any results to read in.
    if( (scalar(keys %{$stats})) == 0 ) {
        system("touch $outputFile");
    }
    
    #Okay this may seem to be overkill, but I'm almost sure it has to be this way (almost).
    #It allows for the case where one HMM can have a hit against the same sequence twice.  
    #Therefore one cmsearch run should be identified not only by the hmm/cm and query sequence, 
    #but also where the hit has occured on the query sequence.
    foreach my $querySeq(keys %{$stats}) {
        foreach my $hmm(keys %{$stats->{$querySeq}}) {
            foreach my $boundaries(keys %{$stats->{$querySeq}->{$hmm}}) {
                my $fileName = $stats->{$querySeq}->{$hmm}->{$boundaries};
                unless($fileName) {
                    print Dumper($stats->{$querySeq}->{$hmm});
                    &_die("Could not find file name");
                }
                my $exitVal = &runProg($fileName, $hmm, $other_opts, $outputFile);
                &handleExitVal($exitVal);
            }
        }
    }
    
}

&cleanUp;

exit(0);
################## END MAIN #####################


##################################### SUB ROUTINES ############################################

#Name: check_parameters
#Desc: checks the input parameters of the program
#Args: Hash Ref - holding the options (from GetOptions routine)
#Rets: Nothing.
sub check_parameters {
    my $options = shift;

    #I wish this function could be uglier.  So I'm not commenting it.
    if($options{input_list} && $options{input_list} ne "") {
        &_die("input_list [$options{input_list}] does not exist") unless( -e $options{input_list});
        open(IN, "< $options{input_list}") or &_die("Unable to open $options{input_list}");
        while(<IN>) {
            push(@files, $_);
        }
        close(IN);
    } elsif($options{input_file} && $options{input_file} ne "") {
        &_die("input_file [$options{input_file}] does not exist") unless( -e $options{input_file});
        push(@files,$options{input_file});
    } else {
        &_die("Either input_list or input_file must be provided");
    }

    if($options{tmp_dir}) {
        &_die("$options{tmp_dir} (tmp_dir) does not exist") unless( -d $options{tmp_dir});
    } else {
        &_die("Option tmp_dir must be provided");
    }
    $outputDir = $options{tmp_dir};
    $outputDir =~ s|/$||;

    &_die("output_file option is required") unless($options{output_file});
    $outputFile = $options{output_file};
    system("rm -f $outputFile") if(-e $options{output_file});

    $debug = $options{debug} if($options{debug});

    if($options{sequence_list}) {
        &_die("sequence_list $options{sequence_list} does not exist") unless(-e $options{sequence_list});
    } else {
        &_die("Option sequence_list must be provided");
    }
    open(IN, "<$options{sequence_list}") or &_die("Unable to open $options{sequence_list}");
    @seqList = <IN>;
    close(IN);

    if($options{hmm_cm_table} && $options{hmm_cm_table} ne "") {
        if(-f $options{hmm_cm_table}) {
            $hmmCmTable = $options{hmm_cm_table};
        } elsif(-d $options{hmm_cm_table}) {
            $cmDir = $options{hmm_cm_table};
        } else {
            &_die("$options{hmm_cm_dir} does not exist");
        }
    } else {
        &_die("Option --hmm_cm_table is required");
    }

    $other_opts = $options{other_opts} if($options{other_opts});
    $extraSeqLen = $options{flanking_seq} if($options{flanking_seq});
    
    if($options{cmsearch_bin} && $options{cmsearch_bin} ne "") {
        &_die("Could not locate cmsearch binary at $options{cmsearch_bin}") unless(-e $options{cmsearch_bin});
        $PROG_NAME = $options{cmsearch_bin};
    }
}

#Name: cleanUp
#Desc: Will remove tmp files and directories made.
#Args: None (Uses @dirsMade array)
#Rets: Nothing
sub cleanUp {
    foreach my $dir(@dirsMade) {
        print STDERR "DEBUG: removing $dir\n";
        system("rm -rf $dir");
    }
}

#Name: createHmmCmLookup
#Desc: Will create a hash lookup from a file that coordinates hmm -> cmm files from the 
#      supplied table file. (Tab delimited matching an hmm with a cm).
#Args: String: File name of the table
#Rets: Nothing
sub createHmmCmLookup {
    my $lookupInput = shift;
    print "In createHmmCmLookup\n" if($debug > 2);
    open(IN, "< $lookupInput") or 
        &_die("Unable to open $lookupInput for reading");
    while(<IN>) {
        my @tmpList = split(/\s/);
        $hmmCmLookup{$tmpList[0]} = $tmpList[1];
        &_die("Check the format of the HMM->CM lookup file.  The HMM should be the HMM name, ".
              "and the CM should be a full path to the CM file.  See perldoc for more details.")
            unless(-e $tmpList[1]);
    }
    close(IN);
  
}

#Name: findFile
#Desc: Given a sequence, will search for the file name in the @seqList array
#Args: Part of the name of the file
#Rets: The sequence in that file.
sub findFile {
    my $seqID = shift;
    my ($retval,$fileFound);
    my $fileFlag = 0; 
    $seqID =~ s/\|/\_/g;
    print "Searching with $seqID\n";
    foreach my $file (@seqList) {
        if($file =~ /$seqID\./) {
            $fileFound = $file;
            print $fileFound."\n";
            $fileFlag++;
        }
    }

    if($fileFlag > 1) {
        &_die("Found $fileFlag possible files");
    } elsif($fileFlag == 0) {
        &_die("Didn't find a match for $seqID in the sequence_list");
    }

    open(IN, "< $fileFound") or &_die("Unable to open $fileFound to get sequence information");
    
    my $header = "";
    my %sequence;
    while(<IN>) {
        my $line = $_;
        chomp($line);

        if($line =~ />(.*)\n*$/) {
            $header = $1;
            $sequence{$header} = "";
        } else {
            $sequence{$header}.=$line;
        }
    }

    if(scalar(keys %sequence) > 1) {
        &_die("More than one sequence found in file $fileFound.  ".
              "Sorry, but only single sequence fasta format files are accepted (for now at least).");
    } elsif($header eq "") {
        &_die("Couldn't find sequence in file $fileFound.  Perhaps it's not fasta format?");
    }
    
    return $sequence{$header};
}

#Name: getInputFiles
#Desc: Opens the passed in file name and makes an array of all the listed files
#Args: String: the file name fo the input list
#Rets: List: An array of the contained files.
sub getInputFiles {
    my $fileName = shift;
    my @retval;
    my $inH;

    open($inH, "< $fileName") or &_die("Could not open $fileName for reading");

    while(<$inH>) {
        push(@retval,$_);
    }
    close($inH);
    
    chomp(@retval);

    return \@retval;

}

#Name: getSequenceAndPrint
#Desc: Retrieves sequence from the database (by looking up id's stored in a hash ref created by
#      processHmmpfamFile.
#Args: Hash Ref: Stats returned from processHmmpfamFile
#Rets: Hash Ref: The same hash ref, just with sequences included.
sub getSequencesAndPrint {
    my $seqHash = shift;
    

    #Loop through sequences
    foreach my $querySeq (keys %{$seqHash}) {

        my $tmpSeq = "";
        $tmpSeq = &findFile($querySeq);
        &_die("Couldn't find file matching $querySeq") unless($tmpSeq ne "");

        #Loop through all the hits for a query sequence
        # (one hit is defined by a unique combination of query sequence,
        # hmm file, start and stop locations).
        foreach my $hmm (keys %{$seqHash->{$querySeq}}) {

            foreach my $startEnd (keys %{$seqHash->{$querySeq}->{$hmm}}) {

                #Since the start and end coordinates are going to change, remove the old entry.
                delete($seqHash->{$querySeq}->{$hmm}->{$startEnd});
                
                #Get the start and end of the hit
                my ($start, $end) = split(/::/, $startEnd);
                
                #Parse sequence information
                my ($tmpSection, $start, $end) = &parseSeqAndExtra($tmpSeq, $start, $end);
                
                #Store file name information and print sequence
                $seqHash->{$querySeq}->{$hmm}->{"${start}::$end"} = &printSeqToFile($querySeq, $hmm, $start, 
                                                                                  $end, $tmpSection);

            
              

            }

        }

    }

    #Return the result
    return $seqHash;
}

#Name: handleExitVal
#Desc: Will die with a correct message depending on the exit value of runProg 
#Args: Hash Ref :: Run information (returned from run prog)
#Rets: Nothing
sub handleExitVal {
    my $exitVal = shift;
    &_die("command:\n$exitVal->{cmd}\nError: @{$exitVal->{err}}") 
        unless(@{$exitVal->{err}} == 0);
}

#Name: lookupCM
#Desc: Looks up a certain cm file from an $hmm name.
#Args: Scalar: HMM name.  Must be in the HMM->CM table.
#Rets: Scalar: The CM related to the HMM.
sub lookupCM {
    my $hmm = shift;
    my $retval;
    
    if($cmDir) {
        $hmm = $1 if($hmm =~ /([^_]+)/);
        opendir(IN,$cmDir) or &_die("Could not open directory cm_dir: [$cmDir] ($!)");
        print "Opening $cmDir\n" if($debug > 2);
        my @posCM = grep { /$hmm/ && -f "$cmDir/$_" } readdir(IN);
        closedir(IN);

        if(@posCM == 0) {
            &_die("Could not find a match for hmm $hmm in $cmDir");
        } elsif(@posCM > 1) {
            $" = ", ";
            &_die("Found more than one possible match for $hmm in $cmDir:\n@posCM");
            $" = " ";
        } else {
            $retval = "$cmDir/$posCM[0]";
        }
    } else {
        
        $retval = $hmmCmLookup{$hmm};
    }

    return $retval;
}

#Name: parseSeqAndExtra
#Desc: Parses the sequence information out of the database.  Will take $extraSeqLen nucleotides
#      on either side of the boudaries provided
#Args: Scalar: Full sequence
#      Scalar: The start boundary to be parsed
#      Scalar: The end boundary
#Rets: Returns the new start and end coordinates (after $extraSeqLen) and the newly parsed
#      sequence
sub parseSeqAndExtra {
    my ($seq, $start, $end) = @_;

    if($start < $extraSeqLen) {
        $start = 0;
    } else {
        $start -= $extraSeqLen;
    }

    #A little sloppy here.  If it's over the length, it will just
    #take it all, so doesn't matter if $end > length($seq)
    $end += $extraSeqLen;

    my $retval = substr($seq, $start, $end-$start);

    return ($retval, $start, $end);
}

#Name: printSeqToFile
#Desc: Will print a sequence to file in the output directory in fasta format
#Args: Scalarx5: QuerySequence ID, the hmm file used in the analysis, the start and end boudnaries
#                and the sequence.
#Rets: Nothing
sub printSeqToFile {
    my ($querySeq, $hmm, $start, $end, $tmpSeq) = @_;
    my $outFileName;
    $querySeq =~ s/\|/\_/g;
    unless(-d "$outputDir/$querySeq") {
        system("mkdir $outputDir/$querySeq");
        push(@dirsMade, "$outputDir/$querySeq");
    }

    my $tmpHmm = $1 if($hmm =~ /.*::(.*)/);

    $outFileName = "$outputDir/$querySeq/$tmpHmm.$start.$end.fsa";
    open(OUT, "> $outFileName") 
        or &_die("Unable to open output file $outFileName ($!)");
    
    print OUT ">${querySeq}::$start-$end\n$tmpSeq";
    close(OUT);

    return $outFileName;
    

}

#Name: processHmmpfamFile
#Desc: Takes in an hmmpfam output file and returns a hash of full of the hits (summarized)
#Args: Scalar: Hmmpfam name
#Rets: Hash Ref: $hashRef->{seqId}->{hmmFile::hmmName}->{start::end}
sub processHmmpfamFile {
    my $fileName = shift;
    my $hmmpfamStats;  #Hash Ref
    my ($inH,$hmmFile,$querySeq);
    my $alignmentFlag = 0;

    open($inH, "<$fileName") or &_die("Unable to open $fileName ($!)");
    print "Parsing $fileName\n" if($debug > 2);
    
    while(<$inH>) {
        if(/^HMM file:\s+(.*)/) {
            $hmmFile = $1;
            chomp $hmmFile;
            print "$hmmFile was used\n" if($debug > 2);
        } elsif(/^Query sequence:\s(.*)/) {
            $querySeq = $1;
            chomp $querySeq;
            print "Query Seq :: $querySeq\n" if($debug > 2);
        } elsif(/^Alignments/) {
            $alignmentFlag = 1;
            print "Found the alignments section\n" if($debug > 2);
        } elsif($alignmentFlag && /^(\S+):.*from\s(\d+)\sto\s(\d+)/) { 
            print "Found start: $2\n" if($debug > 3);
            print "Found end: $3\n" if($debug > 3);
            $hmmpfamStats->{$querySeq}->{"${hmmFile}::$1"}->{"$2::$3"} = "";
        } 
    }
    close($inH);
    return $hmmpfamStats;

}

#Name: runProg
#Desc: Runs the cmsearch program
#Args: Scalar :: Fasta file path
#      Scalar :: Covariance model path
#      Scalar :: Other options
#Rets: Hash Ref :: std - stdout of program
#                  err - stderr of program
#                  exitVal - the exitVal
#                  cmd - the command that was run
sub runProg {
    my ($fsaFile, $hmm, $oOpts, $outputFile) = @_;
    my $retval;

    my $cm = &lookupCM($1) if($hmm =~ /.*::(.*)/);
    &_die("Could not parse HMM name from line $hmm") unless($cm);
    #Make sure the -W option isn't used in the other opts.  
    $oOpts =~ s/--window\s\S+//;

    #Get the length of the fasta file
    my $length;
    open(IN, "< $fsaFile") or &_die("Unable to open $fsaFile to determine length ($!)");
    while(<IN>) {
        if(/^[^>]/) {
            $length = length($_);
        }
    }
    &_die("Could not determine length of sequence in file $fsaFile") unless($length);

    #set up the cmsearch command
    my $cmd = $PROG_NAME." $oOpts $cm $fsaFile";
    print " running [$cmd]\n" if($debug > 2);

    #Run the command and store it's std out and err and exitval.
    open3(undef, \*STD, \*ERR, $cmd);
    my $exitVal = $? << 8;
    $retval->{exitVal} = $exitVal;
    $retval->{std} = [<STD>];
    $retval->{err} = [<ERR>];
    $retval->{cmd} = $cmd;

    #Print the std out to a file.  The output file.  If it exists, concatenate, otherwise
    #just make it.
    my $openOut = "> $outputFile";
    $openOut = ">".$openOut if(-e $outputFile);
    open(OUT, $openOut) or &_die("Unable to open output file $outputFile ($!) [$openOut]");
    print OUT "$cmd\n";
    print OUT "@{$retval->{std}}\n";
    close(OUT);

    #Return the info about the cmsearch run
    return $retval;
}


#Name: _pod
#Desc: Because I'm too lazy to put this line twice.
#Args: none
#Rets: none
sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

#Name: _die
#Desc: Because I'm too lazy to write $logger->logdie(blah) everywhere
#Args: Scalar: The program's last words
#Rets: Nothing, this function will never exit successfully. Shame...
sub _die {
    my $msg = shift;
    &cleanUp;
    $logger->logdie($msg);
}
##EOF
