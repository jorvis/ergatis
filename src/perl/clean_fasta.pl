#!/usr/bin/perl -w
# Copyright @ 2002 - 2010 The Institute for Genomic Research (TIGR).
# All rights reserved.
# 
# This software is provided "AS IS".  TIGR makes no warranties, express or
# implied, including no representation or warranty with respect to the
# performance of the software and derivatives or their safety,
# effectiveness, or commercial viability.  TIGR does not warrant the
# merchantability or fitness of the software and derivatives for any
# particular purpose, or that they may be exploited without infringing the
# copyrights, patent rights or property rights of others.
# 
# This software program may not be sold, leased, transferred, exported or
# otherwise disclaimed to anyone, in whole or in part, without the prior
# written consent of TIGR.

use strict;
use lib '/usr/local/packages/perllib';
use TIGR::Foundation;
use TIGR::FASTAreader;
use File::Copy;
use File::Basename;


our @DEPEND = qw(TIGR::Foundation TIGR::FASTAreader);
our $REVISION = (qw$Revision: 5284 $)[-1];
our $VERSION = '1.26';
our $VERSION_STRING = "$VERSION (Build $REVISION)"; 

my $APP = basename($0,"");

my $wtime = 15; # input waiting time in seconds for usage #2

my $__HELP_INFO = qq~
NAME
$APP - a script that rearranges a FASTA file to make it correct according
         to the FASTA specification.
USAGE

   1. $APP [OPTIONS] <input_file> [-o <output_file>]

   2. cat|zcat <fasta|compressed file|> | $APP [OPTIONS] -o <output_file>
	Note: $APP waits $wtime seconds and quit if no input during that time.

   <input_file> - the FASTA file to be modified.
   -o|output <output_file> - the output file for the clean FASTA information.  

OPTIONS  
   -c|check check if the FASTA file is parse-able and write a message to stdout.
   -u|uniq  make the identifiers of all the records unique.
   -p|progress <number> display clean progress as a progress bar.
   -w|wait <seconds> assign input waiting time for usage #2, default is $wtime.

SUMMARY
The $APP program is used to modify a FASTA file according to the FASTA
specification so that it parses. If the -check option is specified, the result
of the parsing is written to stdout. if the -check option is not specified, 
the file is modified if required and parsed. If the user does not specify a 
filename with the -o option, the original file is replaced with the "clean" 
file and the original file is copied into a file with a "_orig" suffix. 

The script makes the following changes to a FASTA file:
1. All the data lines in the FASTA file are made equal to 60 bases.
2. Non-graphical characters are removed from the data.
3. Spaces, forward slashes and <CTRL-A>s are removed from the identifier.
4. Non-graphical characters (other than <CTRL-A> from the non identifier 
   portion) are removed from the headers.
5. Empty lines are removed from the file.
6. If the -uniq option is specified, the record identifiers are made unique.

EXIT CODES
   0   The FASTA file parsed after cleaning.
   1   There was an error in program execution
   2   The FASTA file did not parse(only with the check option)
   3   The FASTA file parsed(only with the check option)
   4   The FASTA file did not parse after cleaning.

SEE ALSO
  createFasta, alignFasta, mergeFasta, clrFasta, cmpFasta, splitfasta, cutFasta
~;

# Adding help and depend information using TIGR Foundation
my $tf_obj = new TIGR::Foundation;
my $help = 0; # this variable is 1 when the help information is required.
my $depend = 0;  # this variable is 1 when program dependency information is 
                 # required.
my $check = 0; # this variable is 1 when the user wants to check if the file 
               # parsed.
my $uniq = 0; # this variable is 1 when the record identifiers are to be made 
              #unique.
my $output_file = undef; #The output file for writing the clean FASTA 
                         #information.
my $pbar = 0; # don't display progress bar by default

$tf_obj->addDependInfo(@DEPEND);
$tf_obj->setVersionInfo($VERSION_STRING);
$tf_obj->setHelpInfo($__HELP_INFO);

# a flag which specifies that an error should be printed to stderr as well as 
# the error file
my $ERR_FLAG = 1;
# program default parameters
my $LINE_LENGTH = 60;

# Declaring a hash which contains (record identifier, frequency) pairs
my %idenTable = ();

# This function returns a string representation of the FASTA record.
# It takes a header and a data component. It makes changes to the data 
# and the header to a make the record more "fasta_like".
# This string representation conforms to the TIGR definition of a
# FASTA record.  On failure, this function returns undefined.
sub cleanRecord($$) {
   my $header = shift;
   my $data = shift;
   my $identifier = undef;
   
   if(defined $header) {
      #removing all spaces after > in the header.
      $header =~ s/^> +/>/g;
      #retrieving the identifer for the header.
      ($identifier) = $header =~ /^>(\S+).*$/;
   }  else {
	$tf_obj->logError("the header is not defined");}

   if(defined ($identifier)) {
      print STDERR "\n" if $pbar;
      print("Cleaning the record $identifier\n");
      #removing forward slashes and CTRL characters from the identifier.
      $identifier =~ s/[[:cntrl:]\/]//g;
  
      if($uniq == 1) { #the identifiers have to be made unique.
         if(!defined($idenTable{$identifier})) {
            $idenTable{$identifier} = 1;}
         else {
	    my $freq = $idenTable{$identifier};
            $idenTable{$identifier} = $freq + 1;
            $identifier.="_$freq";}
      }
      $header =~ s/^>\S+\s*/>$identifier /g;
      $tf_obj->logLocal("Cleaned the identifier for $identifier", 3);
      
      #removing all non-graphical characters from the header except CTRL-A
      my @ctrl_list = $header =~ /([[:cntrl:]])/g;
      my $list_len = @ctrl_list;
      my $ctrl_char = "";
      my $ctrl_str = "";
      for(my $i = 0; $i<$list_len; $i++) {
         if($ctrl_list[$i] eq "\cA") {
            $ctrl_list[$i] = "";
         }
         $ctrl_str.=$ctrl_list[$i];
      }
      $header =~ s/[$ctrl_str\t\n\f\a\e]//g;
      $header =~ s/[^[:graph:][:cntrl:]\s]//g;
      $tf_obj->logLocal("Cleaned the header for $identifier", 3);
   
      # removing all non-graphical characters from the data.
      if(defined $data) {
	$$data =~ s/[[:cntrl:]\s\a\e]//g;
	#removing all no-nucleotide and non-peptide bases from the data.
	$$data =~ s/[^ATUGCMRWSYKVHDBN\.\-QEZILFPX\*]//gi;
	$tf_obj->logLocal("Cleaned the data for $identifier", 3);
      }
      else {
         $tf_obj->logError("the data is not defined for identifier");}
  
      if(!defined($identifier)) {
         $identifier = "<undef>";
      }
   }
   $header .= "\n";
   $tf_obj->logLocal("Printing the header information in the output file", 2);
   if((defined $header) && ($header !~ /^\s*$/) && ($header ne "") 
	&& (defined $data) && ($$data ne "")) {
        # && ($$data !~ /^\s*$/) # -> that causes segmentation fault
        # printing the header information in the output file
        print FASTAOUT $header;
    }
   return $data;
}

sub makeEqual($) {
    # making all data lines equal.
    my $data = shift;
    my $formatted_data = undef;
    my $seg = undef;
    my $err_condition = 1;

    $err_condition = undef if(!defined $data);

    while ( ( defined ( $seg = substr $$data, 0, $LINE_LENGTH, '' ) ) 
	    && ( $seg ne '' ) ) { 
		#$formatted_data .= $seg . "\n"; 
		print FASTAOUT $seg ."\n";}

    $tf_obj->logLocal("Made all data line lengths equal", 2);
    # printing formatted data to output_file
    #if(defined $formatted_data) {
    #	print FASTAOUT $formatted_data;
    #}
    return ( defined ( $err_condition ) ) ? 1 : undef;
}

sub alarm_h () {
    alarm(0);
    $tf_obj->bail("No input during ".$wtime."s");
}

MAIN:
{
   my @errs_list_bef_clean = ();
   my @errs_list_aft_clean = ();
   my $EXEC_ERROR = 1; #exit code if there is an error in program execution
   my $FILE_INCORRECT = 2; #exit code if the file did not parse.
   my $FILE_CORRECT = 3; #exit code if the file parsed.
   my $FILE_INCORRECT_AFTER_CLEANING = 4; #exit code if the file did not parse
                                          #after cleaning.
   my $FILE_CORRECT_AFTER_CLEANING = 0; #exit code if the file parsed after
                                        #cleaning.
   my $input_file = undef; # the FASTA file to be cleaned
   my $err_condition = 1; # this variable keeps track of execution errors

   $tf_obj->logLocal("taking in the options from the command line",1);
   my $result = $tf_obj->TIGR_GetOptions (
    'c|check', \$check, 
    'u|uniq', \$uniq, 
    'p|progress:i', \$pbar, 
    'w|wait:i', \$wtime,
    'o|output=s' => \$output_file);

   $tf_obj->logLocal("performing the option dependent checks", 2);
   
   # incorrect execution command                               
   if ($result == 0) {
      $tf_obj->logError("Command line parsing failed\n", $ERR_FLAG);
      $err_condition = undef; 
   }

  $pbar = abs($pbar);
  # the warnings file where data content errors are stored
  my $warningsfile = "$APP.warnings";

  # grab input file from command line
  if ( scalar (@ARGV) == 1 ) {
      $input_file = shift @ARGV;
      $tf_obj->logLocal("Got input file = \'$input_file\'", 2);
      if((!defined $input_file) || (!-r $input_file) ||
         (-z $input_file)) {
	    $tf_obj->logError("the inputfile cannot be used", $ERR_FLAG);
	    $err_condition = undef;
	    $warningsfile = undef;}

    # Exit after notifying the user whether the file is parseable
    if(($check == 1) && (defined $err_condition)) {
      my $fasta_reader_bef_clean = undef; 
      if((defined $err_condition) && (defined $tf_obj)) {
         # create a new FASTAreader object to check if file is parseable.
         $fasta_reader_bef_clean = new TIGR::FASTAreader $tf_obj,
    	    \@errs_list_bef_clean,"$input_file";
      }

      if((!defined ($fasta_reader_bef_clean)) && (defined $err_condition)) { 
         #the file did not parse.
         print STDOUT "$input_file did not parse before cleaning\n";
      
         if ( scalar(@errs_list_bef_clean) > 0 ) { # are there parse errors?
            if ( ! open (WARNINGS, ">$warningsfile") ) {
               $tf_obj->logError
	       ("Cannot open the warnings file: '$warningsfile'");
               $err_condition = undef;
            }
            if(defined $err_condition) {
               print STDOUT "look at $warningsfile for details\n";
	       print WARNINGS 
	        "These are the errors in the file before cleaning\n";
               while ( @errs_list_bef_clean ) { # get the messages from the list
                  my $message = shift @errs_list_bef_clean; 
                  print WARNINGS $message, "\n";
               }
               print WARNINGS "\n";
            }
	 }
         exit($FILE_INCORRECT);
      }
      elsif(defined $err_condition) {
         print("$input_file parsed before cleaning\n");
         exit($FILE_CORRECT); 
      }
    }
  
    if (defined $err_condition) {
      # rename the input file in the same directory 
      # with the "_orig" extension if the output file is not specified
      my $input = $input_file;
      if(!defined $output_file) { 
         $input = "$input_file"."_orig";
         if((rename ($input_file, $input)) == 1) {
            $tf_obj->logLocal("renamed the $input_file to $input", 4);
	 }
         else {
            $tf_obj->logError
	    ("Cannot make a copy of the input_file: \'$input_file\'");
            $err_condition = undef;
         }
      }
      # opening the input file 
      if ( ! open (FASTAIN, $input) ) {
         $tf_obj->logError("Cannot open input_file: '$input'");
         $err_condition = undef;
     }
    }
   
  } elsif ($output_file) {

    $SIG{'ALRM'} = \&alarm_h;
    alarm($wtime);
    open (FASTAIN, "-"); # open STDIN

  } else {
 
    $tf_obj->printHelpInfoAndExit();
  } # if ( scalar (@ARGV) = 1 ) 

   #setting the default output file to be the same as the input file
   if((!defined($output_file)) && (defined $err_condition)) {
      $output_file = "$input_file";}
   
   if((defined $output_file) && ( ! open (FASTAOUT, ">$output_file") )) {
      $tf_obj->logError("Cannot open $output_file for writing");
      $err_condition = undef;}
  
    my $header = "";
    my $data = "";
    my $cleanData = "";
    my $formatted_data = "";
    my $line;

    while ((defined $err_condition ) && (defined($line = <FASTAIN>))) {

	chomp $line;

	#if the line is blank ignore it.
	next if($line eq "");

	if( $line =~ /^>/ ) { #if the line is a header.
         if( $header ne '' ) {
            #clean the record and print header
	    if(defined $err_condition ) {
               
               if(!defined($cleanData = cleanRecord($header, \$data))) {
                  $tf_obj->logError("the record was not cleaned");
                  $err_condition = undef;
		}
		$data = $$cleanData;}

            if((defined $err_condition) && (!defined makeEqual(\$data))) { 
               $tf_obj->logError("the record data could not be split ".
        	"into sixty character lines");
               $err_condition = undef;} 
	}
        $header = $line;
        $data = "";
	print STDERR ">" if $pbar;
      }
      else {
    	    $data .= $line; #concatenate all data lines.
	    print STDERR ">" if ($pbar && !($.%$pbar));
      }
    }

   # print last case
   if( $header ne '' ) {
      #clean the record and print header
      if(defined $err_condition ) {
         $cleanData = cleanRecord($header, \$data);

         if(!defined $cleanData) {
            $tf_obj->logError("the record was not cleaned");
            $err_condition = undef;
	 }
         $data = $$cleanData;
      }

      if((defined $err_condition) && (!defined makeEqual(\$data))) { 
         $tf_obj->logError("the record data could not be split ".
    	    "into sixty character lines");
         $err_condition = undef;
      } 
   }
   close(FASTAOUT);
   close(FASTAIN);
   
   if(defined $err_condition) {
      print("\nParsing the file after cleaning\n"); 
      # create a new FASTAreader object to see if file 
      # is parseable after cleaning.
      my $fasta_reader_aft_clean = new TIGR::FASTAreader $tf_obj, 
              \@errs_list_aft_clean, $output_file;
       
      if(!defined ($fasta_reader_aft_clean)) { 
         #the file did not parse after cleaning
         print STDOUT "$input_file did not parse after cleaning\n";
      
         if ( scalar(@errs_list_aft_clean) > 0 ) { # are there parse errors?
            if ( ! open (WARNINGS, ">$warningsfile") ) {
               $tf_obj->logError("Cannot open the warnings file: ".
                           "\'$warningsfile\'");
               $err_condition = undef;
            }
            if(defined $err_condition) {
               print STDOUT "look at $warningsfile for details\n";
	       print WARNINGS 
	        "These are the errors in the file after cleaning\n";
               while (@errs_list_aft_clean) { # get the messages from the list
                  my $message = shift @errs_list_aft_clean; 
                  print WARNINGS $message, "\n";
               }
               print WARNINGS "\n";
               close(WARNINGS);
            }
	 } 
         exit($FILE_INCORRECT_AFTER_CLEANING);
      }
      else {
         print("$input_file parsed after cleaning\n") if $input_file;
         close(WARNINGS);
         exit($FILE_CORRECT_AFTER_CLEANING);
      }
   }
   else { #there was an error in program execution
      exit($EXEC_ERROR);
   }
}

