#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell
##############################################################################
### This program generates BedGraph and a BigWig file from a BAM file
##############################################################################

use strict;
use Cwd;
use Carp;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Spec;

##############################################################################
### Constants
##############################################################################

use constant FALSE => 0;
use constant TRUE  => 1;

use constant VERSION => '0.1.0';
use constant PROGRAM => eval { ($0 =~ m/(\w+\.pl)$/) ? $1 : $0 };

use constant UCSC_UTIL_DIR => '/usr/local/packages/ucsc_utils';
use constant BEDTOOLS_BIN_DIR => '/usr/local/packages/bedtools';
use constant SAMTOOLS_BIN_DIR => '/usr/local/bin';

##############################################################################
### Globals
##############################################################################

my %hCmdLineOption = ();
my $sHelpHeader = "\nThis is ".PROGRAM." version ".VERSION."\n";

my (@aRegions);
my ($sRefFile, $sBamFile, $sAnnotationFile);
my ($sGenicBedFile, $sExonicBedFile, $sIntronicBedFile, $sIntergenicBedFile);
my ($sOutDir, $sSizeFile, $sOutFile, $sSortedFile, $sSummaryFile);
my ($sSenseFile, $sAntiSenseFile, $sBedGraph, $sBigWig);
my ($sSenseBedGraph, $sAntisenseBedGraph, $sSenseBigWig, $sAntisenseBigWig);
my ($sPrefix, $sSampleId, $sRegion, $nTotalMappedReads, $nUniqueMappedReads);
my ($sCmd);
my ($bDebug, $bVerbose);

##############################################################################
### Main
##############################################################################

GetOptions( \%hCmdLineOption,
            'reffile|r=s', 'infile|i=s', 'stranded|s=s', 'outdir|o=s',
            'ucsc_util_dir|uu=s', 'bedtools_bin_dir|bb=s', 'samtools_bin_dir|sb=s', 
            'single_end|e', 'include_wig|w', 'verbose|v',
            'debug',
            'help',
            'man') or pod2usage(2);

if ($hCmdLineOption{'help'} || 
	(! defined $hCmdLineOption{'reffile'}) || 
	(! defined $hCmdLineOption{'infile'})) {
    pod2usage( -msg => $sHelpHeader, -exitval => 1);
}

if ((defined $hCmdLineOption{'stranded'}) &&
	($hCmdLineOption{'stranded'} !~ m/^firststrand$/i) &&
	($hCmdLineOption{'stranded'} !~ m/^secondstrand$/i)) {
    die "\tERROR: Incorrect library type specified. (firststrand or secondstrand)";
}

pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};

$bDebug   = (defined $hCmdLineOption{'debug'}) ? TRUE : FALSE;
$bVerbose = (defined $hCmdLineOption{'verbose'}) ? TRUE : FALSE;

$sOutDir = File::Spec->curdir();
if (defined $hCmdLineOption{'outdir'}) {
    $sOutDir = $hCmdLineOption{'outdir'};

    if (! -e $sOutDir) {
        mkdir($hCmdLineOption{'outdir'}) ||
            die "ERROR! Cannot create output directory\n";
    }
    elsif (! -d $hCmdLineOption{'outdir'}) {
            die "ERROR! $hCmdLineOption{'outdir'} is not a directory\n";
    }
}
$sOutDir = File::Spec->canonpath($sOutDir);

if (! (defined $hCmdLineOption{'ucsc_util_dir'}) ) {
	$hCmdLineOption{'ucsc_util_dir'} = UCSC_UTIL_DIR;
}

if (! (defined $hCmdLineOption{'bedtools_bin_dir'}) ) {
	$hCmdLineOption{'bedtools_bin_dir'} = BEDTOOLS_BIN_DIR;
}

if (! (defined $hCmdLineOption{'samtools_bin_dir'}) ) {
	$hCmdLineOption{'samtools_bin_dir'} = SAMTOOLS_BIN_DIR;
}

# process reference file

$hCmdLineOption{'reffile'} = File::Spec->rel2abs($hCmdLineOption{'reffile'});
($_, $_, $sRefFile) = File::Spec->splitpath($hCmdLineOption{'reffile'});

symlink $hCmdLineOption{'reffile'}, "$sOutDir/$sRefFile";

if ( -e "$hCmdLineOption{'reffile'}.fai" ) {
	symlink "$hCmdLineOption{'reffile'}.fai", "$sOutDir/$sRefFile.fai";
}
else {
	($bDebug || $bVerbose) ? 
		print STDERR "\nIndexing $sRefFile .....\n" : ();
	
	$sCmd = $hCmdLineOption{'samtools_bin_dir'}."/samtools faidx".
			" ".$sOutDir."/".$sRefFile;
	
	exec_command($sCmd);
	
	($bDebug || $bVerbose) ? 
		print STDERR "Indexing $sRefFile ..... done\n" : ();
}

# process chromosome size file

$sPrefix = $sRefFile;
$sPrefix =~ s/\.([A-za-z])+$//;

$sSizeFile = File::Spec->rel2abs("$sOutDir/$sPrefix.chromosome.sizes");

($bDebug || $bVerbose) ? 
	print STDERR "\nGenerating $sSizeFile .....\n" : ();

$sCmd = "cut -f1,2 $sOutDir/$sRefFile.fai > $sSizeFile";

exec_command($sCmd);

($bDebug || $bVerbose) ? 
	print STDERR "Generating $sSizeFile ..... done\n" : ();

$hCmdLineOption{'infile'} = File::Spec->rel2abs($hCmdLineOption{'infile'});
($_, $_, $sBamFile) = File::Spec->splitpath($hCmdLineOption{'infile'});

$sSampleId = $sBamFile;
$sSampleId =~ s/\.bam$//;

# generate sense and antisense BAM file

if (defined $hCmdLineOption{'stranded'}) {
	($sSenseFile, $sAntiSenseFile) = Generate_Stranded_BAM(\%hCmdLineOption, $hCmdLineOption{'infile'}, $sSampleId, $sOutDir);
}

# generate bedGraph file

if (defined $hCmdLineOption{'stranded'}) {
	$sSenseBedGraph = File::Spec->rel2abs("$sOutDir/$sSampleId.sense.bedGraph");
	Bam2BedGraph(\%hCmdLineOption, $sSenseFile, $sSizeFile, $sSenseBedGraph);
	
	$sOutFile = File::Spec->rel2abs("$sOutDir/$sSampleId.antisense.tmp");
	Bam2BedGraph(\%hCmdLineOption, $sAntiSenseFile, $sSizeFile, $sOutFile);
	
	$sAntisenseBedGraph = File::Spec->rel2abs("$sOutDir/$sSampleId.antisense.bedGraph");
	Reverse_BedGraph(\%hCmdLineOption, $sOutFile, $sAntisenseBedGraph);
}
else {
	$sBedGraph = File::Spec->rel2abs("$sOutDir/$sSampleId.bedGraph");
	Bam2BedGraph(\%hCmdLineOption, $hCmdLineOption{'infile'}, $sSizeFile, $sBedGraph);
}

# generate BigWig file

if (defined $hCmdLineOption{'stranded'}) {
	$sSenseBigWig = File::Spec->rel2abs("$sOutDir/$sSampleId.sense.bw");
	BedGraph2BigWig(\%hCmdLineOption, $sSenseBedGraph, $sSizeFile, $sSenseBigWig);
	
	$sAntisenseBigWig = File::Spec->rel2abs("$sOutDir/$sSampleId.antisense.bw");
	BedGraph2BigWig(\%hCmdLineOption, $sAntisenseBedGraph, $sSizeFile, $sAntisenseBigWig);
}
else {
	$sBigWig = File::Spec->rel2abs("$sOutDir/$sSampleId.bw");
	BedGraph2BigWig(\%hCmdLineOption, $sBedGraph, $sSizeFile, $sBigWig);
}

exit;


##############################################################################
### Subroutines
##############################################################################


# exec_command()
#
# Purpose
#   executes shell command
#
# Required Parameters
#   sCmd			= /path/to/output_directory
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
sub exec_command {
	my $sCmd = shift;
	
	if ((!(defined $sCmd)) || ($sCmd eq "")) {
		die "\nSubroutine::exec_command : ERROR! Incorrect command!\n";
	}
	
	my $nExitCode;
	
	print STDERR "$sCmd\n";
	$nExitCode = system("$sCmd");
	if ($nExitCode != 0) {
		die "\tERROR! Command Failed!\n\t$!\n";
	}
	print STDERR "\n";
	
	return;
}

# samtools_view()
#
# Purpose
#   executes samtools view command
#
# Required Parameters
#   phCmdLineOption		= pointer to hash containing command line options
#   sInFile				= path to input alignment SAM/BAM file
#   sOutFile			= path to output file
#	sParams				= parameters
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
sub samtools_view {
    my $phCmdLineOption		= shift;
    my $sInFile				= shift;
    my $sOutFile			= shift;
    my $sParams				= shift;
	
	my $sSubName = (caller(0))[3];

    if (! ((defined $sInFile) &&
           (defined $sParams) &&
           (defined $sOutFile) &&
           (defined $phCmdLineOption))) {
        croak "$sSubName - required parameters missing\n";
    }
	
	my $sCmd;
	$sCmd = $phCmdLineOption->{'samtools_bin_dir'}."/samtools view".
			" ".$sParams.
			" -o $sOutFile $sInFile";
	
	exec_command($sCmd);
	
	return;
}

# samtools_cat()
#
# Purpose
#   executes samtools cat command
#
# Required Parameters
#   phCmdLineOption		= pointer to hash containing command line options
#   sInFile				= comma-separated list of alignment SAM/BAM files 
#   sOutFile			= path to output file
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
sub samtools_cat {
    my $phCmdLineOption		= shift;
    my $sInFile				= shift;
    my $sOutFile			= shift;
    
	my $sSubName = (caller(0))[3];

    if (! ((defined $sInFile) &&
           (defined $sOutFile) &&
           (defined $phCmdLineOption))) {
        croak "$sSubName - required parameters missing\n";
    }
	
	my $sCmd;
	
	$sInFile =~ s/,/ /g;
	
	$sCmd = $phCmdLineOption->{'samtools_bin_dir'}."/samtools cat".
			" -o $sOutFile $sInFile";
	
	exec_command($sCmd);
	
	return;
}

# samtools_sort()
#
# Purpose
#   executes samtools sort command
#
# Required Parameters
#   phCmdLineOption		= pointer to hash containing command line options
#   sInFile				= comma-separated list of alignment SAM/BAM files 
#   sOutPrefix			= path to output prefix
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
sub samtools_sort {
    my $phCmdLineOption		= shift;
    my $sInFile				= shift;
    my $sOutPrefix			= shift;
    
	my $sSubName = (caller(0))[3];

    if (! ((defined $sInFile) &&
           (defined $sOutPrefix) &&
           (defined $phCmdLineOption))) {
        croak "$sSubName - required parameters missing\n";
    }
	
	my $sCmd;
	
	$sCmd = $phCmdLineOption->{'samtools_bin_dir'}."/samtools sort".
			" $sInFile $sOutPrefix";
	
	exec_command($sCmd);
	
	return;
}

# Generate_Stranded_BAM()
#
# Purpose
#   generate the sense and antisense BAM files for an alignment BAM file
#
# Required Parameters
#   phCmdLineOption = pointer to hash containing command line options
#   sBamFile        = input BAM alignment file
#   sSampleId       = Prefix for output files
#   sOutDir			= output directory
#
# Optional Parameters
#   none
#
# Returns
#   sSenseFile			= path to sense strand BAM file
#	sAntiSenseFile		= path to antisense strand BAM file
#
# Side Effects
#   none
#
# Assumptions
#
# Notes
#
sub Generate_Stranded_BAM {
    my $phCmdLineOption = shift;
    my $sBamFile		= shift;
    my $sSampleId       = shift;
    my $sOutDir			= shift;

    my $sSubName = (caller(0))[3];

    if (! ((defined $sBamFile) &&
           (defined $sSampleId) &&
           (defined $sOutDir) &&
           (defined $phCmdLineOption))) {
        croak "$sSubName - required parameters missing\n";
    }

    # Local variables
    my $bDebug   = (defined $phCmdLineOption->{'debug'}) ? TRUE : FALSE;
    my $bVerbose = (defined $phCmdLineOption->{'verbose'}) ? TRUE : FALSE;
 
    my ($sSenseFile, $sAntiSenseFile);
    my $sOutFile;

    # Start

#    ($bDebug) ? print STDERR "In $sSubName\n" : ();
	
	($bDebug || $bVerbose) ? 
		print STDERR "\nGenerating alignment BAM file for sense reads .....\n" : ();
	
	if (defined $phCmdLineOption->{'single_end'}) {
		$sOutFile = File::Spec->rel2abs("$sOutDir/$sSampleId.sense.bam");
		if ($phCmdLineOption->{'librarytype'} =~ m/^firststrand$/i) {
			samtools_view($phCmdLineOption, $sBamFile, $sOutFile, "-bh -f 16");
		}
		else {
			samtools_view($phCmdLineOption, $sBamFile, $sOutFile, "-bh -F 16");
		}
		
		samtools_sort($phCmdLineOption, $sOutFile, "$sOutDir/$sSampleId.sense.sorted");
		
		$sSenseFile = File::Spec->rel2abs("$sOutDir/$sSampleId.sense.sorted.bam");
	}
	else {
		$sOutFile = File::Spec->rel2abs("$sOutDir/$sSampleId.sense.1.bam");
		if ($phCmdLineOption->{'stranded'} =~ m/^firststrand$/i) {
			samtools_view($phCmdLineOption, $sBamFile, $sOutFile, "-bh -f 128 -F 16");
		}
		else {
			samtools_view($phCmdLineOption, $sBamFile, $sOutFile, "-bh -f 64 -F 16");
		}
		
		$sOutFile = File::Spec->rel2abs("$sOutDir/$sSampleId.sense.2.bam");
		if ($phCmdLineOption->{'stranded'} =~ m/^firststrand$/i) {
			samtools_view($phCmdLineOption, $sBamFile, $sOutFile, "-bh -f 80");
		}
		else {
			samtools_view($phCmdLineOption, $sBamFile, $sOutFile, "-bh -f 144");
		}
		
		$sOutFile = File::Spec->rel2abs("$sOutDir/$sSampleId.sense.bam");
		samtools_cat($phCmdLineOption,
					 "$sOutDir/$sSampleId.sense.1.bam,$sOutDir/$sSampleId.sense.2.bam",
					 $sOutFile);
		samtools_sort($phCmdLineOption, $sOutFile, "$sOutDir/$sSampleId.sense.sorted");
		
		$sSenseFile = File::Spec->rel2abs("$sOutDir/$sSampleId.sense.sorted.bam");
		
		if ( -e $sSenseFile ) {
			unlink("$sOutDir/$sSampleId.sense.1.bam");
			unlink("$sOutDir/$sSampleId.sense.2.bam");
		}
	}
	
	($bDebug || $bVerbose) ? 
		print STDERR "\nGenerating alignment BAM file for sense reads ..... done\n" : ();
	
	# generate BAM file for antisense reads
	
	($bDebug || $bVerbose) ? 
		print STDERR "\nGenerating alignment BAM file for antisense reads .....\n" : ();
	
	if (defined $phCmdLineOption->{'single_end'}) {
		$sOutFile = File::Spec->rel2abs("$sOutDir/$sSampleId.antisense.bam");
		if ($phCmdLineOption->{'librarytype'} =~ m/^firststrand$/i) {
			samtools_view($phCmdLineOption, $sBamFile, $sOutFile, "-bh -F 16");
		}
		else {
			samtools_view($phCmdLineOption, $sBamFile, $sOutFile, "-bh -f 16");
		}
		
		samtools_sort(\%hCmdLineOption, $sOutFile, "$sOutDir/$sSampleId.antisense.sorted");
		
		$sAntiSenseFile = File::Spec->rel2abs("$sOutDir/$sSampleId.antisense.sorted.bam");
	}
	else {
		$sOutFile = File::Spec->rel2abs("$sOutDir/$sSampleId.antisense.1.bam");
		if ($phCmdLineOption->{'stranded'} =~ m/^firststrand$/i) {
			samtools_view($phCmdLineOption, $sBamFile, $sOutFile, "-bh -f 144");
		}
		else {
			samtools_view($phCmdLineOption, $sBamFile, $sOutFile, "-bh -f 80");
		}
		
		$sOutFile = File::Spec->rel2abs("$sOutDir/$sSampleId.antisense.2.bam");
		if ($phCmdLineOption->{'stranded'} =~ m/^firststrand$/i) {
			samtools_view($phCmdLineOption, $sBamFile, $sOutFile, "-bh -f 64 -F 16");
		}
		else {
			samtools_view($phCmdLineOption, $sBamFile, $sOutFile, "-bh -f 128 -F 16");
		}
		
		$sOutFile = File::Spec->rel2abs("$sOutDir/$sSampleId.antisense.bam");
		samtools_cat($phCmdLineOption,
					 "$sOutDir/$sSampleId.antisense.1.bam,$sOutDir/$sSampleId.antisense.2.bam",
					 $sOutFile);
		samtools_sort($phCmdLineOption, $sOutFile, "$sOutDir/$sSampleId.antisense.sorted");
		
		$sAntiSenseFile = File::Spec->rel2abs("$sOutDir/$sSampleId.antisense.sorted.bam");
		
		if ( -e $sAntiSenseFile ) {
			unlink("$sOutDir/$sSampleId.antisense.1.bam");
			unlink("$sOutDir/$sSampleId.antisense.2.bam");
		}
	}
	
	($bDebug || $bVerbose) ? 
		print STDERR "\nGenerating alignment BAM file for antisense reads ..... done\n" : ();
	
#    ($bDebug) ? print STDERR "Leaving $sSubName\n" : ();
    
    return $sSenseFile, $sAntiSenseFile;
}

# Bam2BedGraph()
#
# Purpose
#   generate bedgraph file from alignment BAM file
#
# Required Parameters
#   phCmdLineOption		= pointer to hash containing command line options
#   sInFile				= path to alignment SAM/BAM file
#   sSizeFile			= path to genome size file
#   sOutFile			= path to output bedgraph file
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
#
# Notes
#
sub Bam2BedGraph {
    my $phCmdLineOption		= shift;
    my $sInFile				= shift;
    my $sSizeFile			= shift;
    my $sOutFile			= shift;
    
    my $sSubName = (caller(0))[3];

    if (! ((defined $sInFile) &&
           (defined $sSizeFile) &&
           (defined $sOutFile) &&
           (defined $phCmdLineOption))) {
        croak "$sSubName - required parameters missing\n";
    }

    # Local variables
    my $bDebug   = (defined $phCmdLineOption->{'debug'}) ? TRUE : FALSE;
    my $bVerbose = (defined $phCmdLineOption->{'verbose'}) ? TRUE : FALSE;
    
    my $sCmd;
    
    # Start

#    ($bDebug) ? print STDERR "In $sSubName\n" : ();
	
	($bDebug || $bVerbose) ? 
		print STDERR "\nGenerating bedGraph $sOutFile from alignment BAM $sInFile .....\n" : ();
	
	$sCmd = $phCmdLineOption->{'bedtools_bin_dir'}."/bedtools genomecov".
			" -split -bg".
			" -ibam ".$sInFile.
			" -g ".$sSizeFile.
			" > ".$sOutFile;
	
	exec_command($sCmd);
	
	($bDebug || $bVerbose) ? 
		print STDERR "Generating bedGraph $sOutFile from alignment BAM $sInFile ..... done\n" : ();
	
#    ($bDebug) ? print STDERR "Leaving $sSubName\n" : ();
    
    return;
}

# Reverse_BedGraph()
#
# Purpose
#   reverse depth value of input bedgraph file
#
# Required Parameters
#   phCmdLineOption		= pointer to hash containing command line options
#   sInFile				= path to alignment SAM/BAM file
#   sOutFile			= path to output bedgraph file
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
#
# Notes
#
sub Reverse_BedGraph {
    my $phCmdLineOption		= shift;
    my $sInFile				= shift;
    my $sOutFile			= shift;
    
    my $sSubName = (caller(0))[3];

    if (! ((defined $sInFile) &&
           (defined $sOutFile) &&
           (defined $phCmdLineOption))) {
        croak "$sSubName - required parameters missing\n";
    }

    # Local variables
    my $bDebug   = (defined $phCmdLineOption->{'debug'}) ? TRUE : FALSE;
    my $bVerbose = (defined $phCmdLineOption->{'verbose'}) ? TRUE : FALSE;
    
    my ($sRefId, $nStart, $nEnd, $nDepth);
    my $sCmd;
    my ($fpIN, $fpOUT);
    
    # Start

#    ($bDebug) ? print STDERR "In $sSubName\n" : ();
	
	($bDebug || $bVerbose) ? 
		print STDERR "\nReversing depth values for bedGraph $sInFile .....\n" : ();
	
	open($fpIN, "<$sInFile") or die "\tError! Cannot open $sInFile for reading .....\n";
	open($fpOUT, ">$sOutFile") or die "\tError! Cannot open $sInFile for writing .....\n";
	
	while (<$fpIN>) {
		$_ =~ s/\s+$//;
		
		if ($_ =~ m/^#/) {
			print $fpOUT "$_\n";
			next;
		}
		
		($sRefId, $nStart, $nEnd, $nDepth) = split(/\t/, $_);
		
		$nDepth *= -1;
		
		print $fpOUT "$sRefId\t$nStart\t$nEnd\t$nDepth\n";
	}
	
	close($fpOUT);
	close($fpIN);
	
	($bDebug || $bVerbose) ? 
		print STDERR "Reversing depth values for bedgraph $sInFile ..... done\n" : ();
	
#    ($bDebug) ? print STDERR "Leaving $sSubName\n" : ();
    
    return;
}

# BedGraph2BigWig()
#
# Purpose
#   generate BigWig file from bedGraph file
#
# Required Parameters
#   phCmdLineOption		= pointer to hash containing command line options
#   sInFile				= path to bedGraph file
#   sSizeFile			= path to genome size file
#   sOutFile			= path to output BigWig file
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
#
# Notes
#
sub BedGraph2BigWig {
    my $phCmdLineOption		= shift;
    my $sInFile				= shift;
    my $sSizeFile			= shift;
    my $sOutFile			= shift;
    
    my $sSubName = (caller(0))[3];

    if (! ((defined $sInFile) &&
           (defined $sSizeFile) &&
           (defined $sOutFile) &&
           (defined $phCmdLineOption))) {
        croak "$sSubName - required parameters missing\n";
    }

    # Local variables
    my $bDebug   = (defined $phCmdLineOption->{'debug'}) ? TRUE : FALSE;
    my $bVerbose = (defined $phCmdLineOption->{'verbose'}) ? TRUE : FALSE;
    
    my ($sCmd, $sWigFile);
    
    # Start

#    ($bDebug) ? print STDERR "In $sSubName\n" : ();
	
	($bDebug || $bVerbose) ? 
		print STDERR "\nGenerating BigWig $sOutFile from bedGraph $sInFile .....\n" : ();
	
	$sCmd = $phCmdLineOption->{'ucsc_util_dir'}."/bedGraphToBigWig".
			" ".$sInFile.
			" ".$sSizeFile.
			" ".$sOutFile;
	
	exec_command($sCmd);
	
	($bDebug || $bVerbose) ? 
		print STDERR "Generating BigWig $sOutFile from bedGraph $sInFile ..... done\n" : ();
	
	if (defined $phCmdLineOption->{'include_wig'}) {
		
		$sWigFile = $sOutFile;
		$sWigFile =~ s/.bw$/.wig/;
		
		BedGraph2Wig($phCmdLineOption, $sInFile, $sWigFile);
	}
	
#    ($bDebug) ? print STDERR "Leaving $sSubName\n" : ();
    
    return;
}

# BedGraph2Wig()
#
# Purpose
#   generate Wig file from BedGraph file
#
# Required Parameters
#   phCmdLineOption		= pointer to hash containing command line options
#   sInFile				= path to BedGraph file
#   sOutFile			= path to output Wig file
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
#
# Notes
#
sub BedGraph2Wig {
    my $phCmdLineOption		= shift;
    my $sInFile				= shift;
	my $sWigFile			= shift;
    
    my $sSubName = (caller(0))[3];

    if (! ((defined $sInFile) &&
           (defined $sWigFile) &&
           (defined $phCmdLineOption))) {
        croak "$sSubName - required parameters missing\n";
    }

    # Local variables
    my $bDebug   = (defined $phCmdLineOption->{'debug'}) ? TRUE : FALSE;
    my $bVerbose = (defined $phCmdLineOption->{'verbose'}) ? TRUE : FALSE;
    
    my ($sPrevRefId, $nPos);
    my ($sRefId, $nStart, $nEnd, $nDepth);
    my ($fpIN, $fpWIG);
    
    # Start

#    ($bDebug) ? print STDERR "In $sSubName\n" : ();
	
	($bDebug || $bVerbose) ? 
		print STDERR "\nGenerating Wig $sWigFile from bedGraph $sInFile .....\n" : ();
	
	open($fpIN, "<$sInFile") or die "\tError! Cannot open $sInFile for reading .....\n";
	open($fpWIG, ">$sWigFile") or die "\tError! Cannot open $sWigFile for writing .....\n";
	
	$sPrevRefId = "undef";
	while (<$fpIN>) {
		$_ =~ s/\s+$//;
		
		($sRefId, $nStart, $nEnd, $nDepth) = split(/\t/, $_);
		
		if ($sRefId ne $sPrevRefId) {
			print $fpWIG "variableStep chrom=".$sRefId."\n";
			$sPrevRefId = $sRefId;
		}
		
		for ($nPos = ($nStart + 1); $nPos <= $nEnd; $nPos++) {
			print $fpWIG "$nPos\t$nDepth\n";
		}			
	}
	
	close($fpWIG);
	close($fpIN);
	
	($bDebug || $bVerbose) ? 
		print STDERR "Generating Wig $sWigFile from bedGraph $sInFile ..... done\n" : ();
	
#    ($bDebug) ? print STDERR "Leaving $sSubName\n" : ();
    
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
	$sFile =~ m/^(\S+)$sExtension$/;
	if (defined $1) {
	    $sOutFile = $sOutDir.'/'.$sDir.'/'.$1;
	}
	else {
	    $sOutFile = $sOutDir.'/'.$sDir.'/'.$sFile;
	}
    $sOutFile = File::Spec->canonpath($sOutFile);

#    ($bDebug || $bVerbose) ? print STDERR "\n" : ();

#    ($bDebug) ? print STDERR "Leaving $sSubName\n" : ();
    
    return $sOutFile;
}


##############################################################################
### POD Documentation
##############################################################################

__END__

=head1 NAME

bam2bigwig.pl - program to generate BedGraph and BigWig files from alignment BAM file.

=head1 SYNOPSIS

    bam2bigwig.pl --r <reffile> --i <alignment_bam_file> [--s <strandedness>] [--o <output_dir>] [--v]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

   
    --r <reference FastA file>                     = /path/to/reference FastA file.

    --i <alignment_BAM_file>                       = /path/to/alignment BAM file sorted by position.

    --s <strandedness>                             = strandedness (firststrand or secondstrand). Optional.

    --o <output_dir>                               = output directory. Optional.

    --uu <ucsc_util_dir>                           = /path/to/ucsc_utility_directory. Optional.

    --bb <bedtools_bin_dir>                        = /path/to/bedtools_bin_directory. Optional.

    --sb <samtools_bin_dir>                        = /path/to/samtools_bin_directory. Optional.

    --v                                            = generate runtime messages. Optional

=head1 DESCRIPTION

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

PERL5LIB environment variable should be set if Bio::Perl is in non-standard
location

=head1 DEPENDENCIES

Bio::Perl

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

There are no known bugs in this module. Please report problems to Amol Shetty
(ashetty@som.umaryland.edu). Patches are welcome.

=head1 AUTHOR

 Amol Carl Shetty
 Bioinformatics Software Engineer II
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
