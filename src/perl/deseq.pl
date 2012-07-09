#!/usr/bin/perl

################################################################################
### POD Documentation
################################################################################

=head1 NAME

deseq.pl - script to analyze differential gene expression from HTSeq read counts.

=head1 SYNOPSIS

    deseq.pl --i input_samples_file [--c htseq_counts_listfile] [--a annotation_file] 
    		 [--o outdir] [--rb R_bin_dir] [--rp R_parameters] [--ra R_script] [--v] 

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --i <input_samples_file>       = /path/to/input tab-delimited file with the following sample information.
                                     #Sample_Name    Group_Name    [HTSeq_Counts_File]

    --c <htseq_counts_listfile>    = /path/to/listfile of htseq counts files for each sample.
                                     The counts filename should be in the format /path/to/<sample_name>.*.counts.
                                     Required if neither alignment nor counts file specified in sample.info file.

    --a <annotation_file>          = /path/to/annotation file in tab-delimited. Optional
                                     First column of annotation file must be same as ID used for read counting.

    --o <output dir>               = /path/to/output directory. Optional. [present working directory]

    --rb <R_bin_dir>               = /path/to/R binary. Optional. [/usr/local/bin]

    --rp <R_parameters>            = additonal R parameters. Optional. [--slave --vanilla]

    --rs <R_script>                = /path/to/R script for DESeq analysis. Optional. [deseq.R]

    --v                            = generate runtime messages. Optional

=head1 DESCRIPTION

The script executes the deseq.R script using the R Statistical Package.

=head1 AUTHOR

 Amol Carl Shetty
 Bioinformatics Software Engineer II
 Institute of Genome Sciences
 University of Maryland
 Baltimore, Maryland 21201

=cut

################################################################################

use strict;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Spec;
use FindBin qw($RealBin);

##############################################################################
### Constants
##############################################################################

use constant FALSE => 0;
use constant TRUE  => 1;

use constant R_BIN_DIR => '/usr/local/packages/R-2.12.0/bin';

use constant VERSION => '1.0.0';
use constant PROGRAM => eval { ($0 =~ m/(\w+\.pl)$/) ? $1 : $0 };

##############################################################################
### Globals
##############################################################################

my %hCmdLineOption = ();
my $sHelpHeader = "\nThis is ".PROGRAM." version ".VERSION."\n";

GetOptions( \%hCmdLineOption,
            'infile|i=s', 'listfile|c=s', 'annotation|a=s', 'outdir|o=s', 
            'r_bin_dir|rb=s', 'r_params|rp=s', 'r_script|rs=s', 
            'verbose|v',
            'debug',
            'help',
            'man') or pod2usage(2);

## display documentation
pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};
pod2usage( -msg => $sHelpHeader, -exitval => 1) if $hCmdLineOption{'help'};

## make sure everything passed was peachy
check_parameters(\%hCmdLineOption);

my (%hSamples, @aSampleNames);
my ($sOutDir, $sSampleName, $sGroupName, $sCountsFile, $sFile, $sRScript);
my ($fpSMPL, $fpLST, $fpOUT);
my ($sCmd, $nI);
my $bDebug   = (defined $hCmdLineOption{'debug'}) ? TRUE : FALSE;
my $bVerbose = (defined $hCmdLineOption{'verbose'}) ? TRUE : FALSE;

################################################################################
### Main
################################################################################

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'infile'} ...\n" : ();

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

$sRScript = File::Spec->rel2abs($hCmdLineOption{'r_script'});

($bDebug || $bVerbose) ? 
	print STDERR "\nGenerating sample info for DESeq Analysis ...\n" : ();

open ($fpSMPL, "$hCmdLineOption{'infile'}") or die "Error! Cannot open $hCmdLineOption{'infile'} for reading !!!\n";

while (<$fpSMPL>) {
	$_ =~ s/\s+$//;
	
	next if ($_ =~ /^#/);
		
	next if ($_ =~ /^$/);
	
	($sSampleName, $sGroupName, $sCountsFile) = (split(/\t/, $_))[0, 1, 2];
	
	$hSamples{$sSampleName} = [$sGroupName, undef];
	
	$hSamples{$sSampleName}[1] = $sCountsFile if (defined $sCountsFile);
	
	push @aSampleNames, $sSampleName;
}

close($fpSMPL);

if ((defined $hCmdLineOption{'listfile'}) && ($hCmdLineOption{'listfile'} !~ m/^$/)) {
	open ($fpLST, "$hCmdLineOption{'listfile'}") or die "Error! Cannot open $hCmdLineOption{'listfile'} for reading !!!\n";
	
	while (<$fpLST>) {
		$_ =~ s/\s+$//;
		
		next if ($_ =~ /^#/);
		
		next if ($_ =~ /^$/);
		
		$sCountsFile = $_;
		
		($_, $_, $sFile) = File::Spec->splitpath($sCountsFile);
		
		foreach $sSampleName (sort keys %hSamples) {
			if ($sFile =~ m/^$sSampleName/) {
				$hSamples{$sSampleName}[1] = $sCountsFile;
				last;
			}
		}
	}
	
	close($fpLST);
}

open ($fpOUT, ">$sOutDir/deseq.sample.info.txt") or die "Error! Cannot open $sOutDir/deseq.sample.info.txt for writing !!!\n";

for ($nI = 0; $nI < @aSampleNames; $nI++) {
	$sSampleName = $aSampleNames[$nI];
	die "Error! Count file not specified for $sSampleName!!!\n" if (! (defined $hSamples{$sSampleName}[1]));
	print $fpOUT "$sSampleName\t".
				 "$hSamples{$sSampleName}[0]\t".
				 "$hSamples{$sSampleName}[1]\n";
}

close($fpOUT);

($bDebug || $bVerbose) ? 
	print STDERR "\nInitiating DESeq Analysis ...\n" : ();

$sCmd  = $hCmdLineOption{'r_bin_dir'}."/R".
		 " ".$hCmdLineOption{'r_params'}.
		 " --args ".$sOutDir."/deseq.sample.info.txt".
		 " ".$sOutDir;

$sCmd .= " ".$hCmdLineOption{'annotation'} if ((defined $hCmdLineOption{'annotation'}) && ($hCmdLineOption{'annotation'} !~ m/^$/));

$sCmd .= " < ".$sRScript;

exec_command($sCmd);

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'infile'} ... done\n" : ();


################################################################################
### Subroutines
################################################################################

sub check_parameters {
    my $phOptions = shift;
    
    ## make sure input parameters are provided
    if (! (defined $phOptions->{'infile'}) ) {
		pod2usage( -msg => $sHelpHeader, -exitval => 1);
	}
	
    ## handle some defaults
    $phOptions->{'r_bin_dir'} = R_BIN_DIR if (! (defined $phOptions->{'r_bin_dir'}) );
    $phOptions->{'r_params'} = '--slave --vanilla' if (! (defined $phOptions->{'r_params'}) );
    $phOptions->{'r_script'} = $RealBin."/deseq_1.2.1.R" if (! (defined $phOptions->{'r_script'}) );
    
    # set environment variables
	set_environment($phOptions);
}

sub set_environment {
	my $phOptions = shift;
	
	umask 0000;
	
	# add R path to global path
	$ENV{PATH} = $phOptions->{'r_bin_dir'}.":".$ENV{PATH};
}

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
}

################################################################################
