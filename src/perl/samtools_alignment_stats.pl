#!/usr/bin/perl

################################################################################
### POD Documentation
################################################################################

=head1 NAME

samtools_alignment_stats.pl - script to generate alignment stats from BAM file

=head1 SYNOPSIS

    samtools_alignment_stats.pl --i input_BAM_file [--o outdir] [--sb samtools_bin_dir] [--v]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --i <input_BAM_file>           = /path/to/alignment BAM file.

    --o <output dir>               = /path/to/output directory. Optional. [present working directory]

    --sb <samtools_bin_dir>        = /path/to/samtools_build binary. Optional. [/usr/local/bin]

    --v                            = generate runtime messages. Optional

=head1 DESCRIPTION

The script executes the samtools flagstat script from the samtools software package.

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
            'infile|i=s', 'outdir|o=s', 'samtools_bin_dir|sb=s', 'verbose|v',
            'debug',
            'help',
            'man') or pod2usage(2);

## display documentation
pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};
pod2usage( -msg => $sHelpHeader, -exitval => 1) if $hCmdLineOption{'help'};

## make sure everything passed was peachy
check_parameters(\%hCmdLineOption);

my ($sOutDir, $sInFile);
my ($sCmd, $sPrefix);
my $bDebug   = (defined $hCmdLineOption{'debug'}) ? TRUE : FALSE;
my $bVerbose = (defined $hCmdLineOption{'verbose'}) ? TRUE : FALSE;

################################################################################
### Main
################################################################################

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'infile'} ...\n" : ();

($_, $_, $sInFile) = File::Spec->splitpath($hCmdLineOption{'infile'});

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

$sCmd = "ln -sf ".$hCmdLineOption{'infile'}." ".$sOutDir."/".$sInFile;

exec_command($sCmd);

$sInFile = $sOutDir."/".$sInFile;

($bDebug || $bVerbose) ? 
	print STDERR "\nGenerating alignment stats for $sInFile ...\n" : ();

($_, $_, $sPrefix) = File::Spec->splitpath($sInFile);
$sPrefix =~ s/.bam$//;

$sCmd  = $hCmdLineOption{'samtools_bin_dir'}."/samtools flagstat ".$sInFile." > $sOutDir/$sPrefix.mapstats.txt";

exec_command($sCmd);

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'infile'} ... done\n" : ();


################################################################################
### Subroutines
################################################################################

sub check_parameters {
    my $phOptions = shift;
    
    ## make sure input fastx is provided
    if (! (defined $phOptions->{'infile'}) ) {
		pod2usage( -msg => $sHelpHeader, -exitval => 1);
	}
	
    ## handle some defaults
    $phOptions->{'samtools_bin_dir'} = BIN_DIR if (! (defined $phOptions->{'samtools_bin_dir'}) );
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
