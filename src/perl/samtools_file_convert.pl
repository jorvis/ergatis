#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

################################################################################
### POD Documentation
################################################################################

=head1 NAME

samtools_file_convert.pl - script to convert alignment files between SAM and BAM file formats.

=head1 SYNOPSIS

    samtools_file_convert.pl --i input_file --f input_file_format --opt options (1|2|3|4|)
                             [--r reference_fasta_file] [--o output_dir] [--sb samtools_bin_dir]
                             [--sv samtools_view_options] [--ss samtools_sort_options] [--v]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --i <input_file>               = /path/to/input SAM or BAM alignment file.

    --f <input_file_format>        = format of input alignment file (SAM or BAM).

    --opt <options (1|2|3|4)>      = String of options indicating conversion.
                                     1 : BAM        --> sorted BAM
                                     2 : sorted BAM --> indexed BAM
                                     3 : BAM        --> SAM
                                     4 : SAM        --> BAM

    --r <reference_fasta_file>     = /path/to/reference FastA file. Optional.

    --o <output_dir>               = /path/to/output directory. Optional. [input file directory]

    --sb <samtools_bin_dir>        = /path/to/samtools binary. Optional. [/usr/local/bin]

    --sv <samtools_view_options>   = additonal samtools view parameters. Optional.

    --ss <samtools_sort_options>   = additonal samtools sort parameters. Optional

    --v                            = generate runtime messages. Optional

=head1 DESCRIPTION

The script executes the samtools script from the samtools software package.

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
            'infile|i=s', 'infile_format|f=s', 'options|opt=s',
            'reffile|r=s', 'outdir|o=s', 'samtools_bin_dir|sb=s',
            'samtools_view_options|sv=s', 'samtools_sort_options|ss=s',
            'verbose|v',
            'debug',
            'help',
            'man') or pod2usage(2);

## display documentation
pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};
pod2usage( -msg => $sHelpHeader, -exitval => 1) if $hCmdLineOption{'help'};

## make sure everything passed was peachy
check_parameters(\%hCmdLineOption);

my ($sOutDir, $sInFile, $sFormat);
my (@aOptions, $nOption);
my ($sCmd, $nI);
my $bDebug   = (defined $hCmdLineOption{'debug'}) ? TRUE : FALSE;
my $bVerbose = (defined $hCmdLineOption{'verbose'}) ? TRUE : FALSE;

################################################################################
### Main
################################################################################

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'infile'} ...\n" : ();

($_, $sOutDir, $sInFile) = File::Spec->splitpath($hCmdLineOption{'infile'});

if (defined $hCmdLineOption{'outdir'}) {
    $sOutDir = $hCmdLineOption{'outdir'};

    if (! -e $sOutDir) {
        mkdir($hCmdLineOption{'outdir'}) ||
            die "ERROR! Cannot create output directory\n";
    }
    elsif (! -d $hCmdLineOption{'outdir'}) {
            die "ERROR! $hCmdLineOption{'outdir'} is not a directory\n";
    }
    
    $sCmd = "ln -sf ".$hCmdLineOption{'infile'}." ".$sOutDir."/".$sInFile;
    exec_command($sCmd);
	
	$sInFile = $sOutDir."/".$sInFile;
}
$sOutDir = File::Spec->canonpath($sOutDir);

$sFormat = $hCmdLineOption{'infile_format'};

@aOptions = split(//, $hCmdLineOption{'options'});

for ($nI = 0; $nI < @aOptions; $nI++) {
	$nOption = $aOptions[$nI];
	
	if ($nOption == 1) {
		die "Error! Input file format is not BAM.\n" if ($sFormat !~ m/^BAM$/i);
		($sInFile, $sFormat) = bam2sorted_bam(\%hCmdLineOption, $sInFile, $sOutDir);
	}
	elsif ($nOption == 2) {
		die "Error! Input file format is not BAM.\n" if ($sFormat !~ m/^BAM$/i);
		($sInFile, $sFormat) = sorted_bam2indexed_bam(\%hCmdLineOption, $sInFile, $sOutDir);
	}
	elsif ($nOption == 3) {
		die "Error! Input file format is not BAM.\n" if ($sFormat !~ m/^BAM$/i);
		($sInFile, $sFormat) = bam2sam(\%hCmdLineOption, $sInFile, $sOutDir);
	}
	elsif ($nOption == 4) {
		die "Error! Input file format is not SAM.\n" if ($sFormat !~ m/^SAM$/i);
		($sInFile, $sFormat) = sam2bam(\%hCmdLineOption, $sInFile, $sOutDir);
	}
	else {
		die "Error! Invalid option $nOption. Value has to been 1, 2, 3 or 4.\n";
	}
}

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'infile'} ... done\n" : ();


################################################################################
### Subroutines
################################################################################

sub check_parameters {
    my $phOptions = shift;
    
    ## make sure input fastx is provided
    if ((! (defined $phOptions->{'infile'}) ) ||
    	(! (defined $phOptions->{'infile_format'}) ) ||
    	(! (defined $phOptions->{'options'}) )) {
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

sub bam2sorted_bam {
    my $phOptions	= shift;
    my $sFile		= shift;
    my $sOutDir		= shift;
    
    ## make sure input file is provided
    if ((! (defined $sFile) ) ||
    	(! (defined $sOutDir) )) {
		die "Error! Undefined parameters in BAM --> Sorted BAM.\n";
	}
	
	my ($sPrefix, $sOutFile, $sFormat);
	my $bDebug   = (defined $phOptions->{'debug'}) ? TRUE : FALSE;
	my $bVerbose = (defined $phOptions->{'verbose'}) ? TRUE : FALSE;
	
	($_, $_, $sPrefix) = File::Spec->splitpath($sFile);
	$sPrefix =~ s/.bam$//;
	
	$sOutFile = $sOutDir."/".$sPrefix.".sorted_by_position";
	if (defined $phOptions->{'samtools_sort_options'}) {
		if ($phOptions->{'samtools_sort_options'} =~ m/-n/) {
			$sOutFile = $sOutDir."/".$sPrefix.".sorted_by_name";
		}
	}
	
    ## samtools execution
    ($bDebug || $bVerbose) ? 
		print STDERR "\nConverting $sFile to $sOutFile.bam ...\n" : ();
	
    $sCmd  = $phOptions->{'samtools_bin_dir'}."/samtools sort";
	$sCmd .= " ".$phOptions->{'samtools_sort_options'} if (defined $phOptions->{'samtools_sort_options'});
	$sCmd .= " ".$sFile." ".$sOutFile;
	
	exec_command($sCmd);
	
	($bDebug || $bVerbose) ? 
		print STDERR "\nConverting $sFile to $sOutFile.bam ... done\n" : ();
	
	$sFormat = "BAM";
	
	return ("$sOutFile.bam", $sFormat);
}

sub sorted_bam2indexed_bam {
    my $phOptions	= shift;
    my $sFile		= shift;
    my $sOutDir		= shift;
    
    ## make sure input file is provided
    if ((! (defined $sFile) ) ||
    	(! (defined $sOutDir) )) {
		die "Error! Undefined parameters in Sorted BAM --> Indexed BAM.\n";
	}
	
	my ($sPrefix, $sOutFile, $sFormat);
	my $bDebug   = (defined $phOptions->{'debug'}) ? TRUE : FALSE;
	my $bVerbose = (defined $phOptions->{'verbose'}) ? TRUE : FALSE;
	
	($_, $_, $sPrefix) = File::Spec->splitpath($sFile);
	$sPrefix =~ s/.bam$//;
	
	$sOutFile = $sOutDir."/".$sPrefix.".bam.bai";
	
    ## samtools execution
    ($bDebug || $bVerbose) ? 
		print STDERR "\nConverting $sFile to $sOutFile ...\n" : ();
	
    $sCmd  = $phOptions->{'samtools_bin_dir'}."/samtools index ".$sFile;
	
	exec_command($sCmd);
	
	($bDebug || $bVerbose) ? 
		print STDERR "\nConverting $sFile to $sOutFile ... done\n" : ();
	
	$sFormat = "BAM";
	
	return ($sFile, $sFormat);
}

sub bam2sam {
    my $phOptions	= shift;
    my $sFile		= shift;
    my $sOutDir		= shift;
    
    ## make sure input file is provided
    if ((! (defined $sFile) ) ||
    	(! (defined $sOutDir) )) {
		die "Error! Undefined parameters in BAM --> SAM.\n";
	}
	
	my ($sPrefix, $sOutFile, $sFormat);
	my $bDebug   = (defined $phOptions->{'debug'}) ? TRUE : FALSE;
	my $bVerbose = (defined $phOptions->{'verbose'}) ? TRUE : FALSE;
	
	($_, $_, $sPrefix) = File::Spec->splitpath($sFile);
	$sPrefix =~ s/.bam$//;
	
	$sOutFile = $sOutDir."/".$sPrefix.".sam";
	
    ## samtools execution
    ($bDebug || $bVerbose) ? 
		print STDERR "\nConverting $sFile to $sOutFile ...\n" : ();
	
    $sCmd  = $phOptions->{'samtools_bin_dir'}."/samtools view -h";
    $sCmd .= " ".$phOptions->{'samtools_view_options'} if (defined $phOptions->{'samtools_view_options'});
    $sCmd .= " -t ".$phOptions->{'reffile'}.".fai"
    	if ((defined $phOptions->{'reffile'}) && ( -e "$phOptions->{'reffile'}.fai"));
	$sCmd .= " -o ".$sOutFile." ".$sFile;
	
	exec_command($sCmd);
	
	($bDebug || $bVerbose) ? 
		print STDERR "\nConverting $sFile to $sOutFile ... done\n" : ();
	
	$sFormat = "SAM";
	
	return ($sOutFile, $sFormat);
}

sub sam2bam {
    my $phOptions	= shift;
    my $sFile		= shift;
    my $sOutDir		= shift;
    
    ## make sure input file is provided
    if ((! (defined $sFile) ) ||
    	(! (defined $sOutDir) )) {
		die "Error! Undefined parameters in BAM --> SAM.\n";
	}
	
	my ($sPrefix, $sOutFile, $sFormat);
	my $bDebug   = (defined $phOptions->{'debug'}) ? TRUE : FALSE;
	my $bVerbose = (defined $phOptions->{'verbose'}) ? TRUE : FALSE;
	
	($_, $_, $sPrefix) = File::Spec->splitpath($sFile);
	$sPrefix =~ s/.sam$//;
	
	$sOutFile = $sOutDir."/".$sPrefix.".bam";
	
    ## samtools execution
    ($bDebug || $bVerbose) ? 
		print STDERR "\nConverting $sFile to $sOutFile ...\n" : ();
	
    $sCmd  = $phOptions->{'samtools_bin_dir'}."/samtools view -bhS";
    $sCmd .= " ".$phOptions->{'samtools_view_options'} if (defined $phOptions->{'samtools_view_options'});
    $sCmd .= " -t ".$phOptions->{'reffile'}.".fai"
    	if ((defined $phOptions->{'reffile'}) && ( -e "$phOptions->{'reffile'}.fai"));
	$sCmd .= " -o ".$sOutFile." ".$sFile;
	
	exec_command($sCmd);
	
	($bDebug || $bVerbose) ? 
		print STDERR "\nConverting $sFile to $sOutFile ... done\n" : ();
	
	$sFormat = "BAM";
	
	return ($sOutFile, $sFormat);
}

################################################################################
