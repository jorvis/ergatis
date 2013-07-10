#!/usr/local/bin/perl -w
##############################################################################
### This program breaks down a multi fastx file into multiple fastx files
### FastQ --> FastQ
##############################################################################

use strict;
use Cwd;
use Carp;
use Getopt::Long;
use Pod::Usage;
use File::Spec;
use File::Temp qw/ tempfile /;

##############################################################################
### Constants
##############################################################################

use constant FALSE => 0;
use constant TRUE  => 1;

use constant VERSION => '2.0.0';
use constant PROGRAM => eval { ($0 =~ m/(\w+\.pl)$/) ? $1 : $0 };

##############################################################################
### Globals
##############################################################################

my %hCmdLineOption = ();
my $sHelpHeader = "\nThis is ".PROGRAM." version ".VERSION."\n";

my (@aInputFiles);
my ($sOutDir, $sOutFile, $sFile);
my ($bDebug, $bVerbose);

##############################################################################
### Main
##############################################################################

GetOptions( \%hCmdLineOption,
            'input|i=s', 'outdir|o=s', 'numseq|n=i',
            'exclude|e=s', 'verbose|v',
            'debug',
            'help',
            'man') or pod2usage(2);
            
if ($hCmdLineOption{'debug'}) {
	$hCmdLineOption{'input'} = "";
	$hCmdLineOption{'outdir'} = "";
	$hCmdLineOption{'numseq'} = "";
}

if ($hCmdLineOption{'help'} || 
	(! defined $hCmdLineOption{'input'})) {
    pod2usage( -msg => $sHelpHeader, -exitval => 1);
}

pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};

$bDebug   = (defined $hCmdLineOption{'debug'}) ? TRUE : FALSE;
$bVerbose = (defined $hCmdLineOption{'verbose'}) ? TRUE : FALSE;

$sOutDir = File::Spec->curdir();
if (defined $hCmdLineOption{'outdir'}) {
    $sOutDir = $hCmdLineOption{'outdir'};

    if (! -e $sOutDir) {
        mkdir($hCmdLineOption{'outdir'}) ||
            croak "ERROR! Cannot create output directory\n";
    }
    elsif (! -d $hCmdLineOption{'outdir'}) {
            croak "ERROR! $hCmdLineOption{'outdir'} is not a directory\n";
    }
}
$sOutDir = File::Spec->canonpath($sOutDir);

if ((!defined $hCmdLineOption{'numseq'}) || !( $hCmdLineOption{'numseq'} )) {
	$hCmdLineOption{'numseq'} = 20000000;
}
	print 	$hCmdLineOption{'numseq'};
  

# process the input file

($bDebug || $bVerbose) ? print STDERR "Split fastq files .....\n" : ();
Split_Fastq_File(\%hCmdLineOption, $hCmdLineOption{'input'}, $sOutDir);

exit;


##############################################################################
### Subroutines
##############################################################################


# Split_Fastq_File()
#
# Purpose
#   split single multi-fastq file into multiple multi-fastx files.
#
# Required Parameters
#   phCmdLineOption = pointer to hash of command line options
#   sInFile         = /path/to/fastq_file
#	sOutDir			= /path/to/output_directory
#   
# Optional Parameters
#   none
#
# Returns
#   nothing
#
# Side Effects
#   none
#
# Assumptions
#   none
#
# Notes
#
sub Split_Fastq_File {
	my $phCmdLineOption = shift;
    my $sInFile         = shift;
    my $sOutDir			= shift;

    my $sSubName = (caller(0))[3];

	# check for required parameters
    if (! ((defined $sInFile) && (defined $sOutDir))) {
        croak "$sSubName - required parameters missing\n";
    }
    
    # Local variables 
       
    my $bDebug   = (defined $phCmdLineOption->{'debug'}) ? TRUE : FALSE;
    my $bVerbose = (defined $phCmdLineOption->{'verbose'}) ? TRUE : FALSE;

    my ($fpIN, $fnOUT, $fpOUT);
    my ($sSubFileBase, $sSubFile, $sNumFile, $sOutFile);
    my ($sHeader, $sRangeSpec, $sSequence, $sRegex);
    my ($nFileNum, $nN, $nX, $nT);
    my ($bFlag);

	# Start
	
    ($bDebug) ? print STDERR "In $sSubName\n" : ();

    ($bDebug || $bVerbose) ? 
        print STDERR "\nSplitting $sInFile\n" : ();
  
    $nFileNum = 1;
    if ($sInFile =~ m/\.txt$/|| $sInFile =~ m/\.txt\.gz$/ ) {
	    $sSubFileBase = Init_OutFileName(\%$phCmdLineOption, $sOutDir, $sInFile, '.txt');

	}
	elsif ($sInFile =~ m/\.fastq$/ || $sInFile =~ m/\.fastq\.gz$/   ) {
	    $sSubFileBase = Init_OutFileName(\%$phCmdLineOption, $sOutDir, $sInFile, '.fastq');
		
	}
	else {
	    print STDERR "\tERROR : Incorrect file extension for $sInFile\n";
	    exit;
	}

	$sSubFileBase =~ s/sequence/dequence/;
	$sSubFile = $sSubFileBase.sprintf("_%08d.%s", $nFileNum, 'fastq');
	$sNumFile =  $sSubFileBase.".num.txt";
	$sOutFile = (File::Spec->splitpath($sSubFile))[2];
	   
    $nN = $nT = 0; 
    $nX = 1;
    
    ($bDebug || $bVerbose) ? print STDERR "\r\tSubFileName : $sOutFile\tSequence Number : $nX" : ();

    # initialize input file objects
	if($sInFile =~/\.gz$/){
		open($fpIN, "gunzip -c $sInFile |") or die "\t Cannot open $sInFile for reading\n";
	}
	else {
		open($fpIN, "<$sInFile") or die "\t Cannot open $sInFile for reading\n";
	}
    open($fpOUT, ">$sSubFile") or die "\t Cannot open $sSubFile for writing\n";
    open($fnOUT, ">$sNumFile") or die "\t Cannot open $sNumFile for writing\n"; 

	while (<$fpIN>) {
		$_ =~ s/\s+$//;
		if (/^\@([A-Za-z0-9#\/_.:-]+.*)$/) {
		
				  $sHeader = $1;
			
				  $bFlag = 0;
				  if (defined $hCmdLineOption{'exclude'}) {
					  $sRegex = $hCmdLineOption{'exclude'};
					  $bFlag = 1 if ($sHeader =~ m/$sRegex/);
				  }
			
				  if ($nX > $phCmdLineOption->{'numseq'}) {
					  close($fpOUT);
					  ($bDebug || $bVerbose) ? print STDERR "\r\tSubFileName : $sOutFile\tSequence Number : ".($nX-1)."\n" : ();
				
					  $nFileNum++;
					  $sSubFile = $sSubFileBase.sprintf("_%08d.%s", $nFileNum, 'fastq');
					  $sOutFile = (File::Spec->splitpath($sSubFile))[2];
					  open($fpOUT, ">$sSubFile") or die "\t Cannot open $sSubFile for writing\n";
					  $nX = 1;
					  ($bDebug || $bVerbose) ? print STDERR "\r\tSubFileName : $sOutFile\tSequence Number : $nX" : ();
				  }
				  
				  if (!($bFlag)) {
					  
					  print $fpOUT "\@$sHeader\n";
				  }
				  
				  while (<$fpIN>) {
					  last if (/^\+$/ || /^\+\Q$sHeader\E$/);
					  $_ =~ s/\s+$//;
					  print $fpOUT uc($_)."\n" if (!($bFlag));
					  $nN++;
				  }
			
				  if (!($bFlag)) {
					  print $fpOUT "+$sHeader\n";
				  }
				  
				  while ($nN > 0) {
					  $_ = <$fpIN>;
					  $_ =~ s/\s+$//;
					  if (!($bFlag)) {
						  print $fpOUT $_."\n";
					  }
					  $nN--;
				  }
				  
#			($bDebug || $bVerbose) ? print STDERR "\r\tSubFileName : $sOutFile\tSequence Number : $nX" : ();
			
				  $nX++;
				  $nT++;
		  }
	   }

			print $fnOUT "$sSubFileBase\t$nFileNum\n";	
			close ($fnOUT);
			close($fpOUT);
			close($fpIN);
			
			($bDebug || $bVerbose) ? print STDERR "\r\tSubFileName : $sOutFile\tSequence Number : ".($nX-1)."\n" : ();
	
			($bDebug || $bVerbose) ? print STDERR "\nTotal number of reads processed : $nT\n\n" : ();
        
			($bDebug) ? print STDERR "\nLeaving $sSubName\n" : ();



			return;
 }
	

# Init_OutFileName()
#
# Purpose
#   generates a new filename based on existing filename by removing the
#   extension, and prepending an output directory
#
# Required Parameters
#   phCmdLineOption = pointer to hash containing command line options
#   sOutDir         = output directory
#   sFileName       = existing filename
#   sExtension      = filename extension
#
# Optional Parameters
#   none
#
# Returns
#   sOutFile        = new filename
#
# Side Effects
#   none
#
# Assumptions
#
# Notes
#
		sub Init_OutFileName {
			my $phCmdLineOption = shift;
			my $sOutDir         = shift;
			my $sFileName       = shift;
			my $sExtension      = shift;
			
			my $sSubName = (caller(0))[3];
			
			if (! ((defined $sFileName) &&
           (defined $sExtension) &&
				   (defined $sOutDir) &&
				   (defined $phCmdLineOption))) {
				croak "$sSubName - required parameters missing\n";
			}

			# Local variables
			my $bDebug   = (defined $phCmdLineOption->{'debug'}) ? TRUE : FALSE;
			my $bVerbose = (defined $phCmdLineOption->{'verbose'}) ? TRUE : FALSE;
 
			my ($sFile, $sDir);
			my $sOutFile;
			
			# Start

#    ($bDebug) ? print STDERR "In $sSubName\n" : ();

#    ($bDebug || $bVerbose) ? print STDERR "\n\tInitializing filename...\n" : ();

			($_, $sDir, $sFile) = File::Spec->splitpath($sFileName);
			if ($sOutDir ne File::Spec->curdir()) {
				$sDir = '';
			}
#	print "\$1: $1\n";
			$sFile =~ /^(\S+)$sExtension/;
			if (defined $1) {
				$sOutFile = $sOutDir.'/'.$sDir.'/'.$1;
			}
			else {
				$sOutFile = $sOutDir.'/'.$sDir.'/'.$sFile;
			}
#	print "$sOutFile\n";
#	exit;
			$sOutFile = File::Spec->canonpath($sOutFile);
			
#    ($bDebug || $bVerbose) ? print STDERR "\n" : ();
			
#    ($bDebug) ? print STDERR "Leaving $sSubName\n" : ();

    
			return $sOutFile
			}


##############################################################################
### POD Documentation
##############################################################################

__END__

=head1 NAME

split_fastq.pl - program to split_fastq file.

=head1 SYNOPSIS

    split_fasq.pl --i <fastq_file> [--o <output_dir>] 
                 [--n < max number of sequences per split file>] [--e <exclude_regex>] [--v]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

   
    --i <fastq_file>                         = fastq file to split.

    --o <output_dir>                               = output directory. Optional.

    --n < max number of sequences per split file>  = maximum number of reads in per output file. Optional.

    --e <exclude regex>                            = regex to exclude reads that match the regex. Optional.

    --v                                            = generate runtime messages. Optional

=head1 DESCRIPTION

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT


location

=head1 DEPENDENCIES

Bio::Perl

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

There are no known bugs in this module. Please report problems to Amol Shetty
(ashetty@som.umaryland.edu). Patches are welcome.

=head1 AUTHOR

 Weizhong Chang, Amol Carl Shetty
 Bioinformatics Software Engineer
 Institute of Genome Sciences
 University of Maryland
 Baltimore, Maryland 21201

=head1 LICENSE AND COPYRIGHT

Copyright (c) 2010 Amol Carl Shetty (<ashetty@som.umaryland.edu>). All rights
reserved.

This program is free software; you can distribute it and/or modify it under
GNU-GPL licenscing terms.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or FITNESS
FOR A PARTICULAR PURPOSE. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=cut

##############################################################################
### Internal Documentation
##############################################################################

# Variable/Subroutine Name Conventions
#
# Type Prefixes
#
#    b = boolean
#    c = character
#    i = integer
#    n = number (real or integer)
#    p = pointer aka reference
#    s = string
#    a = array
#    h = hash
#    r = rec (aka hash)
#    f = file (typically combined with p, e.g. $fpIN)
#    o = object (e.g. for vars created/initialized via module->new() methods)
#
#    Type prefixes can be combined where appropriate, e.g. as = array of
#    strings or ph = pointer to a hash, but not cn (character of numbers???)
#
# Variable Names
#
#    Variable names start with a type prefix followed by capitalized words
#    descriptive of the purpose of the variable. Variable names do not have
#    underscores. E.g. sThisIsAString
#
# Subroutine Names
#
#    Subroutine names do not have a type prefix, have capitalized words
#    descriptive of the purpose of the subroutine, and the words are
#    separated by an underscore. E.g. This_Is_A_Subroutine
#
#    Internal utility subroutines begin with an underscrore and follow the
#    naming conventions of subroutines.
#

