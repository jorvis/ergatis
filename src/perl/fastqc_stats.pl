#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

################################################################################
### POD Documentation
################################################################################

=head1 NAME

fastqc_stats.pl             - script to execute fastqc_quality_stats on Fastq files.

=head1 SYNOPSIS

    fastqc_stats.pl         --s1 fastq seq file1 [--s2 fastq seq file 2][--o outdir] [--b fastqc_bin_dir] 
                           [--v]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --s1 <fastq file>    = /path/to/fastq read file to analyze.
    
    --s2 <fastq file>    = /path/to/fastq paired read file to analyze. Optional

    --o <output dir>     = /path/to/output directory. Optional. [present working directory]

    --b <fastqc_bin_dir> = /path/to/fastx toolkit binary. Optional. [/usr/local/bin]

    --v                  = generate runtime messages. Optional

=head1 DESCRIPTION

The script executes the fastqc statistics script from the fastqc.

=head1 AUTHOR

 Priti Kumari
 Bioinformatics Software Engineer I
 Institute for Genome Sciences
 University of Maryland
 Baltimore, Maryland 21201

=cut

################################################################################

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Spec;

##############################################################################
### Constants
##############################################################################

use constant FALSE => 0;
use constant TRUE  => 1;

use constant BIN_DIR => '/usr/local/bin';

use constant VERSION => 'v0.10.0';
use constant PROGRAM => eval { ($0 =~ m/(\w+\.pl)$/) ? $1 : $0 };

##############################################################################
### Globals
##############################################################################

my %hCmdLineOption = ();
my $sHelpHeader = "\nThis is ".PROGRAM." version ".VERSION."\n";

GetOptions( \%hCmdLineOption,
            'seq1file|s1=s', 
			'seq2file|s2=s',
			'outdir|o=s', 
            'fastqc_bin_dir|b=s',
			'verbose|v',
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
	print STDERR "\nProcessing $hCmdLineOption{'seq1file'} ...\n" : ();

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

#($_, $_, $sPrefix) = File::Spec->splitpath($hCmdLineOption{'seq1file'});

if ((defined $hCmdLineOption{'seq2file'}) && ($hCmdLineOption{'seq2file'} !~ m/^$/)) {
    $sCmd = $hCmdLineOption{'fastqc_bin_dir'}."/fastqc"." ".$hCmdLineOption{'seq1file'}." ".$hCmdLineOption{'seq2file'}." -o ".$sOutDir;
}
else{
    $sCmd = $hCmdLineOption{'fastqc_bin_dir'}."/fastqc"." ".$hCmdLineOption{'seq1file'}." -o ".$sOutDir;
} 	


exec_command($sCmd);
out_format($sOutDir,$hCmdLineOption{'seq1file'});
if (defined $hCmdLineOption{'seq2file'}){
    out_format($sOutDir,$hCmdLineOption{'seq2file'});
}

################################################################################
### Subroutines
################################################################################

sub check_parameters {
    my $phOptions = shift;
    
    ## make sure input fastx file provided
    if (! (defined $phOptions->{'seq1file'}) ) {
		pod2usage( -msg => $sHelpHeader, -exitval => 1);
	}
	
    ## handle some defaults
    #$phOptions->{'quality'} = 64 if (! (defined $phOptions->{'quality'}) );
    $phOptions->{'fastqc_bin_dir'} = BIN_DIR if (! (defined $phOptions->{'fastqc_bin_dir'}) );
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

sub out_format {
	my $outdir = shift;
	my $s_file = shift;

	my($pref,@arr,$cmd,$new_dir,@plots,$temp);
	
	if ($s_file =~m/(.+)\.(\w+)(\.gz)$/){
		if ($2 eq 'fastq' or $2 eq 'txt'){
			$pref = $1;
		}
		else {
			$pref =$1.".".$2;
		}
	}
	
	@arr=split (/\//,$pref);
	$pref=$arr[scalar @arr -1];

	$new_dir=$pref."_"."fastqc";
	
	@plots=('per_base_quality.png','sequence_length_distribution.png');
	
	foreach $temp (@plots) {
		$cmd= "mv ".$outdir."/".$new_dir."/Images/".$temp." ".$outdir."/".$new_dir."/Images/".$pref.".".$temp ;   
		exec_command($cmd);
	}
	
}



################################################################################
