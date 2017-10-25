#!/usr/bin/env perl

eval 'exec /usr/bin/env perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

################################################################################
### POD Documentation
################################################################################

=head1 NAME

fastx_trimming.pl - script to execute fastx_trimmer on FastA/Q file

=head1 SYNOPSIS

    fastx_trimming.pl --i fastx file [--o outdir] [--b fastx_bin_dir] 
                      [--Q quality(33 or 64)] [--f first_base] 
                      [--l last_base] [--qs] [--z] [--v]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --i <fastx file>     = /path/to/fastx read file to analyze.

    --o <output dir>     = /path/to/output directory. Optional. [present working directory]

    --b <fastx_bin_dir>  = /path/to/fastx toolkit binary. Optional. [/usr/local/bin]

    --Q <quality format> = quality string format. Optional. [64]

    --f <first_base>     = first base to keep. Optional. [1]

    --l <last_base>      = last base to keep. Optional. [Read Length]

    --qs                 = derive 'last_base' trim length from quality statistics. Optional.

    --z                  = indicates fastx_file is gzipped. Optional

    --v                  = generate runtime messages. Optional

=head1 DESCRIPTION

The script executes the fastx_trimmer script from the fastx toolkit.

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
            'fastx_bin_dir|b=s', 'first_base|f=i', 'last_base|l=i',
            'qual_stats|qs', 'gzip|z', 'verbose|v',
            'debug',
            'help',
            'man') or pod2usage(2);

## display documentation
pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};
pod2usage( -msg => $sHelpHeader, -exitval => 1) if $hCmdLineOption{'help'};

## make sure everything passed was peachy
check_parameters(\%hCmdLineOption);

my ($sOutDir, $sOutFile);
my ($nReadLength);
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

if (defined $hCmdLineOption{'qual_stats'}) {
	($bDebug || $bVerbose) ? 
		print STDERR "\nGenerating Quality Statistics for $hCmdLineOption{'infile'} ...\n" : ();
	
	if (defined $hCmdLineOption{'gzip'}) {
		$sCmd = "zcat ".$hCmdLineOption{'infile'}." | ".
				$hCmdLineOption{'fastx_bin_dir'}."/fastx_quality_stats";
		
		$sPrefix =~ s/.(\w+).gz$/.base_quality_stats/;
	}
	else {
		$sCmd = $hCmdLineOption{'fastx_bin_dir'}."/fastx_quality_stats".
				" -i ".$hCmdLineOption{'infile'};
		
		$sPrefix =~ s/.(\w+)$/.base_quality_stats/;
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
		print STDERR "\nDetermining trim length for $hCmdLineOption{'infile'} ...\n" : ();
		
	($nReadLength, $hCmdLineOption{'last_base'}) = trim_length(\%hCmdLineOption, "$sOutDir/$sPrefix.txt");
	
	($bDebug || $bVerbose) ? 
		print STDERR "\tDerived Trim Length\t: $hCmdLineOption{'last_base'}\n" : ();
}

($bDebug || $bVerbose) ? 
	print STDERR "\nTrimming reads for $hCmdLineOption{'infile'} ...\n" : ();

($_, $_, $sOutFile) = File::Spec->splitpath($hCmdLineOption{'infile'});

if (defined $hCmdLineOption{'gzip'}) {
	$sCmd = "zcat ".$hCmdLineOption{'infile'}." | ".
			$hCmdLineOption{'fastx_bin_dir'}."/fastx_trimmer -z";
	
	$sOutFile =~ s/.(\w+).gz$/.trimmed.\1.gz/;
}
else {
	$sCmd = $hCmdLineOption{'fastx_bin_dir'}."/fastx_trimmer".
			" -i ".$hCmdLineOption{'infile'};
	
	$sOutFile =~ s/.(\w+)$/.trimmed.\1/;
}

$sCmd .= " -Q ".$hCmdLineOption{'quality'};
$sCmd .= " -f ".$hCmdLineOption{'first_base'} if (defined $hCmdLineOption{'first_base'});
$sCmd .= " -l ".$hCmdLineOption{'last_base'} if (defined $hCmdLineOption{'last_base'});
$sCmd .= " -o ".$sOutDir."/".$sOutFile;
$sCmd .= " -v" if ($bDebug || $bVerbose);

if ((defined $hCmdLineOption{'last_base'}) && (defined $nReadLength)) {
	if (($hCmdLineOption{'last_base'} - $hCmdLineOption{'first_base'} + 1) == $nReadLength) {
		$sCmd = "ln -sf $hCmdLineOption{'infile'} $sOutDir/$sOutFile";
	}
}

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
    $phOptions->{'quality'} = 64 if (! (defined $phOptions->{'quality'}) );
    $phOptions->{'fastx_bin_dir'} = BIN_DIR if (! (defined $phOptions->{'fastx_bin_dir'}) );
    $phOptions->{'first_base'} = 1 if (! (defined $phOptions->{'first_base'}) );
}

sub trim_length {
	my $phOptions	= shift;
	my $sStatsFile	= shift;
	
	my (@aMedQuals, @aSortedMedQuals);
	my ($nMedQual, $nCount, $nMedianOfMedians, $nTrimLength);
	my ($nI);
	my ($fpIN);
	
	## read quality statistics file and extract median quality at each base position
	@aMedQuals = ();
	open($fpIN, "$sStatsFile") or die "\tCannot open $sStatsFile for reading .....\n";
	while(<$fpIN>) {
		$_ =~ s/\s+$//;
		next if ($_ !~ /^\d/);
		
		$nMedQual = (split(/\t/, $_))[7];
		push @aMedQuals, $nMedQual;
	}
	close($fpIN);
	
	$nCount = scalar @aMedQuals;
	@aSortedMedQuals = sort { $a <=> $b } @aMedQuals;
	
	if ($nCount % 2) {
		$nMedianOfMedians = $aSortedMedQuals[int($nCount/2)];
	}
	else {
		$nMedianOfMedians = ($aSortedMedQuals[int($nCount/2)] + $aSortedMedQuals[(int($nCount/2) - 1)]) / 2;
	}
	
	$nTrimLength = $nCount;
	for ($nI = $#aMedQuals; $nI >= 0; $nI--) {
		last if ($aMedQuals[$nI] >= ($nMedianOfMedians - 10));
		$nTrimLength = $nI;
	}
	
	return ($nCount, $nTrimLength);
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
