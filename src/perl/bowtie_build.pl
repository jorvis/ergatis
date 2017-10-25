#!/usr/bin/env perl

eval 'exec /usr/bin/env perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

################################################################################
### POD Documentation
################################################################################

=head1 NAME

bowtie_build.pl - script to create Bowtie index for a reference FastA file

=head1 SYNOPSIS

    bowtie_build.pl --r reference_fasta_file [--o outdir] [--b bowtie_bin_dir] 
                    [--p index_prefix] [--a other_parameters] [--v]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --r <reference_fasta_file>     = /path/to/reference fasta file to index.

    --o <output dir>               = /path/to/output directory. Optional. [present working directory]

    --b <bowtie_bin_dir>           = /path/to/bowtie_build binary. Optional. [/usr/local/bin]

    --p <index_prefix>             = prefix for index files. Optional. [ref]

    --a <other_parameters>         = additonal bowtie_build parameters. Optional

    --v                            = generate runtime messages. Optional

=head1 DESCRIPTION

The script executes the bowtie_build script from the bowtie software package.

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
            'reffile|r=s', 'outdir|o=s', 'prefix|p=s',
            'bowtie_bin_dir|b=s', 'args|a=s', 'verbose|v',
            'debug',
            'help',
            'man') or pod2usage(2);

## display documentation
pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};
pod2usage( -msg => $sHelpHeader, -exitval => 1) if $hCmdLineOption{'help'};

## make sure everything passed was peachy
check_parameters(\%hCmdLineOption);

my ($sOutDir, $sOutFile);
my ($sCmd, $sPrefix);
my $bDebug   = (defined $hCmdLineOption{'debug'}) ? TRUE : FALSE;
my $bVerbose = (defined $hCmdLineOption{'verbose'}) ? TRUE : FALSE;

################################################################################
### Main
################################################################################

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'reffile'} ...\n" : ();

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
	print STDERR "\nGenerating Bowtie index for $hCmdLineOption{'reffile'} ...\n" : ();

$sCmd = "ln -sf ".$hCmdLineOption{'reffile'}." ".$sOutDir."/".$hCmdLineOption{'prefix'}.".fa";

exec_command($sCmd);

$sCmd  = $hCmdLineOption{'bowtie_bin_dir'}."/bowtie-build";
$sCmd .= " ".$hCmdLineOption{'args'} if (defined $hCmdLineOption{'args'});
$sCmd .= " -q" if (!($bDebug || $bVerbose));
$sCmd .= " -f ".$sOutDir."/".$hCmdLineOption{'prefix'}.".fa".
		 " ".$sOutDir."/".$hCmdLineOption{'prefix'};

exec_command($sCmd);

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'reffile'} ... done\n" : ();


################################################################################
### Subroutines
################################################################################

sub check_parameters {
    my $phOptions = shift;
    
    ## make sure input fastx is provided
    if (! (defined $phOptions->{'reffile'}) ) {
		pod2usage( -msg => $sHelpHeader, -exitval => 1);
	}
	
    ## handle some defaults
    $phOptions->{'prefix'} = "ref" if (! (defined $phOptions->{'prefix'}) );
    $phOptions->{'bowtie_bin_dir'} = BIN_DIR if (! (defined $phOptions->{'bowtie_bin_dir'}) );
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
