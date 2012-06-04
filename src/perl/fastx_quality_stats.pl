#!/usr/bin/perl

################################################################################
### POD Documentation
################################################################################

=head1 NAME

fastx_quality_stats.pl - script to execute fastx_quality_stats on FastA/Q file

=head1 SYNOPSIS

    fastx_quality_stats.pl --i fastx file [--o outdir] [--b fastx_bin_dir] 
                           [--Q quality(33 or 64)] [--z] [--v]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --i <fastx file>     = /path/to/fastx read file to analyze.

    --o <output dir>     = /path/to/output directory. Optional. [present working directory]

    --b <fastx_bin_dir>  = /path/to/fastx toolkit binary. Optional. [/usr/local/bin]

    --Q <quality format> = quality string format. Optional. [64]

    --z                  = indicates fastx_file is gzipped. Optional

    --v                  = generate runtime messages. Optional

=head1 DESCRIPTION

The script executes the fastx_quality_stats script from the fastx toolkit.

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
            'infile|i=s', 'outdir|o=s', 'quality|Q=i',
            'fastx_bin_dir|b=s', 'gzip|z', 'verbose|v',
            'debug',
            'help',
            'man') or pod2usage(2);

## display documentation
pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};
pod2usage( -msg => $sHelpHeader, -exitval => 1) if $hCmdLineOption{'help'};

## make sure everything passed was peachy
check_parameters(\%hCmdLineOption);

my ($sOutDir);
my ($sCmd, $sPrefix);
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

($_, $_, $sPrefix) = File::Spec->splitpath($hCmdLineOption{'infile'});

if (defined $hCmdLineOption{'gzip'}) {
	$sCmd = "zcat ".$hCmdLineOption{'infile'}." | ".
			$hCmdLineOption{'fastx_bin_dir'}."/fastx_quality_stats";
	
	$sPrefix =~ s/.(\w+).gz$/.base_quality_stats/;
}
else {
	$sCmd = $hCmdLineOption{'fastx_bin_dir'}."/fastx_quality_stats".
			" -i ".$hCmdLineOption{'infile'};
	
	$sPrefix =~ s/.(\w+)/.base_quality_stats/;
}


$sCmd .= " -Q ".$hCmdLineOption{'quality'};
$sCmd .= " -o ".$sOutDir."/$sPrefix.txt";
$sCmd .= " -v" if ($bDebug || $bVerbose);

exec_command($sCmd);

($bDebug || $bVerbose) ? 
	print STDERR "\nGenerating Quality Score Box Plot for $hCmdLineOption{'infile'} ...\n" : ();

$sCmd = $hCmdLineOption{'fastx_bin_dir'}."/fastq_quality_boxplot_graph.sh".
		" -i ".$sOutDir."/$sPrefix.txt".
		" -o ".$sOutDir."/$sPrefix.png";

exec_command($sCmd);

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'infile'} ... done\n" : ();


################################################################################
### Subroutines
################################################################################

sub check_parameters {
    my $phOptions = shift;
    
    ## make sure input fastx file provided
    if (! (defined $phOptions->{'infile'}) ) {
		pod2usage( -msg => $sHelpHeader, -exitval => 1);
	}
	
    ## handle some defaults
    $phOptions->{'quality'} = 64 if (! (defined $phOptions->{'quality'}) );
    $phOptions->{'fastx_bin_dir'} = BIN_DIR if (! (defined $phOptions->{'fastx_bin_dir'}) );
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
