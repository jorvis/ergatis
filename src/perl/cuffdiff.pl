#!/usr/bin/env perl

eval 'exec /usr/bin/env perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

################################################################################
### POD Documentation
################################################################################

=head1 NAME

cuffdiff.pl - script to generate Cuffdiff transcript differential expression for Cufflinks/Cuffcompare results.

=head1 SYNOPSIS

    cuffdiff.pl --i transcripts.gtf --l list_of_SAM_alignment_files [--o outdir] 
                [--labels sample#1,sample#2] [--num-threads threads] [--fdr FDR cut-off]
                [--cb cufflinks_bin_dir] [--args other_parameters] [--v] 

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --i <transcripts.gtf>          = /path/to/transcripts in GTF format from Cufflinks or Cuffcompare.

    --l <list_SAM_alignment_files> = /path/to/listfile of SAM alignment file(s) for each sample.
                                     Replicates' alignment files can be comma-separated for a single sample.
                                     eg: Sample#1_Replicate#1,Sample#1_Replicate#2,Sample#1_Replicate#3
                                         Sample#2_Replicate#1,Sample#2_Replicate#2,Sample#2_Replicate#3

    --o <output dir>               = /path/to/output directory. Optional. [present working directory]

    --labels <sample#1,sample#2>   = Label Ids for each sample to be analyzed. Optional.

    --num-threads <# threads>      = Use # threads to align reads. Optional. [1]

    --fdr <FDR cut-off>            = The allowed false discovery rate. Optional. [0.05]

    --cb <cufflinks_bin_dir>       = /path/to/cufflinks bin directory. Optional. [/usr/local/bin]

    --args <other_params>          = additional Cuffdiff parameters. Optional. Refer Cufflinks manual.

    --v                            = generate runtime messages. Optional

=head1 DESCRIPTION

The script executes the cuffdiff script from the Cufflinks RNA-Seq Analysis package

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

##############################################################################
### Constants
##############################################################################

use constant FALSE => 0;
use constant TRUE  => 1;

use constant BIN_DIR => '/usr/local/bin';

use constant VERSION => '1.0.0';
use constant PROGRAM => eval { ($0 =~ m/(\w+\.pl)$/) ? $1 : $0 };

##############################################################################
### Globals
##############################################################################

my %hCmdLineOption = ();
my $sHelpHeader = "\nThis is ".PROGRAM." version ".VERSION."\n";

GetOptions( \%hCmdLineOption,
            'gtffile|i=s', 'listfile|c=s', 'outdir|o=s', 
            'labels|L=s', 'num-threads|p=i', 'min-alignment-count=i', 'fdr=f', 
            'library-type=s', 'max-mle-iterations=i', 'min-isoform-fraction=f', 'max-bundle-frags=i', 
            'cufflinks_bin_dir|cb=s', 'args=s', 
            'verbose|v',
            'debug',
            'help',
            'man') or pod2usage(2);

## display documentation
pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};
pod2usage( -msg => $sHelpHeader, -exitval => 1) if $hCmdLineOption{'help'};

## make sure everything passed was peachy
check_parameters(\%hCmdLineOption);

my (@aSamFiles, @aFiles);
my ($sOutDir, $sFile, $sSamFiles, $sSamFileList, $sPrefix);
my ($fpLST);
my ($sCmd, $sKey, $nI);
my $bDebug   = (defined $hCmdLineOption{'debug'}) ? TRUE : FALSE;
my $bVerbose = (defined $hCmdLineOption{'verbose'}) ? TRUE : FALSE;

################################################################################
### Main
################################################################################

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'listfile'} ...\n" : ();

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

($bDebug || $bVerbose) ? 
	print STDERR "\nExecuting Cuffdiff Differential Isoform Expression Analysis for input $hCmdLineOption{'listfile'} ...\n" : ();

@aSamFiles = ();

open ($fpLST, "$hCmdLineOption{'listfile'}") or die "Error! Cannot open $hCmdLineOption{'listfile'} for reading !!!\n";

while (<$fpLST>) {
	$_ =~ s/\s+$//;
	
	next if ($_ =~ /^#/);
		
	next if ($_ =~ /^$/);
	
	@aFiles = split(/,/, $_);
	
	$sSamFiles = undef;
	
	for ($nI = 0; $nI < @aFiles; $nI++) {
		($_, $_, $sFile) = File::Spec->splitpath($aFiles[$nI]);
		$sCmd = "ln -sf $aFiles[$nI] $sOutDir/$sFile";
		
		exec_command($sCmd);
		
		$sSamFiles .= ((defined $sSamFiles) ? ",$sOutDir/$sFile" : "$sOutDir/$sFile");
	}
	
	push @aSamFiles, $sSamFiles;
}

close($fpLST);

$sSamFileList = join(" ", @aSamFiles);

($bDebug || $bVerbose) ? 
	print STDERR "\nInitiating CuffDiff Analysis ...\n" : ();

$sCmd  = $hCmdLineOption{'cufflinks_bin_dir'}."/cuffdiff".
		 " -o ".$sOutDir;

foreach $sKey ( keys %hCmdLineOption) {
	next if (($sKey eq "gtffile") || ($sKey eq "listfile") || ($sKey eq "outdir") || 
			 ($sKey eq "cufflinks_bin_dir") || ($sKey eq "library-type") || ($sKey eq "fdr") || ($sKey eq "args") || 
			 ($sKey eq "verbose") || ($sKey eq "debug") || ($sKey eq "help") || ($sKey eq "man") );
	
	$sCmd .= " --".$sKey." ".$hCmdLineOption{$sKey} if ((defined $hCmdLineOption{$sKey}) && ($hCmdLineOption{$sKey} !~ m/^$/i));
}

$sCmd .= " --library-type ".$hCmdLineOption{'library-type'} if ((defined $hCmdLineOption{'library-type'}) && ($hCmdLineOption{'library-type'} !~ m/^$/i));
$sCmd .= " --FDR ".$hCmdLineOption{'fdr'} if (defined $hCmdLineOption{'fdr'});
if ((! defined $hCmdLineOption{'labels'}) || ($hCmdLineOption{'labels'} =~ m/^$/i)) {
	($_, $_, $sFile) = File::Spec->splitpath($hCmdLineOption{'listfile'});
	if ($sFile =~ m/(.+)vs(.+)\.sample.list/) {
		$sCmd .= " --labels ".$2.",".$1;
	}
}
$sCmd .= " ".$hCmdLineOption{'args'} if (defined $hCmdLineOption{'args'});
$sCmd .= " ".$hCmdLineOption{'gtffile'};
$sCmd .= " ".$sSamFileList;

exec_command($sCmd);

($bDebug || $bVerbose) ? 
	print STDERR "\nRenaming CuffDiff output files ...\n" : ();

($_, $_, $sPrefix) = File::Spec->splitpath($hCmdLineOption{'listfile'});

$sPrefix =~ s/.list$//;

rename_files($sOutDir, ".diff", $sPrefix);

rename_files($sOutDir, ".fpkm_tracking", $sPrefix);

rename_files($sOutDir, ".count_tracking", $sPrefix);

rename_files($sOutDir, ".read_group_tracking", $sPrefix);

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'listfile'} ... done\n" : ();


################################################################################
### Subroutines
################################################################################

sub check_parameters {
    my $phOptions = shift;
    
    ## make sure input parameters are provided
    if ((! (defined $phOptions->{'gtffile'}) ) ||
    	(! (defined $phOptions->{'listfile'}) )) {
		pod2usage( -msg => $sHelpHeader, -exitval => 1);
	}
	
    ## handle some defaults
    $phOptions->{'cufflinks_bin_dir'} = BIN_DIR if (! (defined $phOptions->{'cufflinks_bin_dir'}) );
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

sub rename_files {
    my $sOutDir = shift;
    my $sSuffix	= shift;
    my $sPrefix = shift;
    
    ## make sure input parameters are provided
    if ((! (defined $sOutDir) ) ||
    	(! (defined $sSuffix) ) ||
    	(! (defined $sPrefix) )) {
		print STDERR "\nSubroutine::rename_files\n\tIncomplete parameter list!!!!!\n";
	}
	
	my (@aInputFiles);
	my ($sBasename);
	
	@aInputFiles = glob("$sOutDir/*$sSuffix");
	
	foreach $sFile (@aInputFiles) {
		($_, $_, $sBasename) = File::Spec->splitpath($sFile);
		
		rename ($sFile, "$sOutDir/$sPrefix.$sBasename") or printf STDERR "\tError! Cannot rename $sFile !!!\n";
	}
}

################################################################################
