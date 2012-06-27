#!/usr/bin/perl

################################################################################
### POD Documentation
################################################################################

=head1 NAME

cuffcompare.pl - script to compare Cufflinks transcript GTF files against a reference or within samples.

=head1 SYNOPSIS

    cuffcompare.pl --i cufflinks_gtf] [--a annotation_file] [--p outprefix]
    		       [--o outdir] [--cb cufflinks_bin_dir] [--v] 

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --i <cufflinks_gtf>            = /path/to/cufflinks gtf file for a single sample or 
                                     a listfile of cufflinks gtf files for multiple samples.
                                     The gtf filename should be in the format /path/to/<sample_name>.*.transcripts.gtf

    --a <annotation_file>          = /path/to/annotation file in GTF format. Optional

    --p <prefix>                   = Output files prefix. Optional. [cuffcmp]

    --o <output dir>               = /path/to/output directory. Optional. [present working directory]

    --cb <cufflinks_bin_dir>       = /path/to/cufflinks binary directory. Optional. [/usr/local/bin]

    --args <other_params>          = additional Cuffcompare parameters. Optional. Refer Cufflinks manual.

    --v                            = generate runtime messages. Optional

=head1 DESCRIPTION

The script executes the cuffcompare script from the Cufflinks RNA-Seq Analysis package.

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
            'gtffile|i=s', 'annotation|a=s', 'prefix|p=s', 'outdir|o=s', 
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

my (@aGtfFile, $sGtfFile, $sGtfFileList);
my ($sOutDir, $sSampleName, $sGroupName, $sFile, $sPrefix);
my ($fpLST);
my ($sCmd);
my $bDebug   = (defined $hCmdLineOption{'debug'}) ? TRUE : FALSE;
my $bVerbose = (defined $hCmdLineOption{'verbose'}) ? TRUE : FALSE;

################################################################################
### Main
################################################################################

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'gtffile'} ...\n" : ();

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
	print STDERR "\nComparing Cufflinks Transcript Analysis ...\n" : ();

($_, $_, $sFile) = File::Spec->splitpath($hCmdLineOption{'gtffile'});

if ($sFile =~ m/\.gtf$/) {
	$sCmd = "ln -sf $hCmdLineOption{'gtffile'} $sOutDir/$sFile";
	
	exec_command($sCmd);
	
	push @aGtfFile, "$sOutDir/$sFile";
	
	$sPrefix = $hCmdLineOption{'prefix'};
}
else {
	open ($fpLST, "$hCmdLineOption{'gtffile'}") or die "Error! Cannot open $hCmdLineOption{'gtffile'} for reading !!!\n";
	
	while (<$fpLST>) {
		$_ =~ s/\s+$//;
		
		next if ($_ =~ /^#/);
		
		next if ($_ =~ /^$/);
		
		$sGtfFile = $_;
		
		($_, $_, $sFile) = File::Spec->splitpath($sGtfFile);
		
		$sCmd = "ln -sf $sGtfFile $sOutDir/$sFile";
	
		exec_command($sCmd);
		
		push @aGtfFile, "$sOutDir/$sFile";
	}
	
	close($fpLST);
	
	($_, $_, $sPrefix) = File::Spec->splitpath($hCmdLineOption{'gtffile'});
	
	$sPrefix =~ s/.list$//; 
}

$sGtfFileList = join(" ", @aGtfFile);

($bDebug || $bVerbose) ? 
	print STDERR "\nInitiating Cuffcompare Analysis ...\n" : ();

$sCmd  = $hCmdLineOption{'cufflinks_bin_dir'}."/cuffcompare".
		 " -o ".$sOutDir."/".$sPrefix;
$sCmd .= " -r ".$hCmdLineOption{'annotation'} if ((defined $hCmdLineOption{'annotation'}) && ($hCmdLineOption{'annotation'} !~ m/^$/));
$sCmd .= " ".$hCmdLineOption{'args'} if (defined $hCmdLineOption{'args'});
$sCmd .= " ".$sGtfFileList;

exec_command($sCmd);

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'gtffile'} ... done\n" : ();


################################################################################
### Subroutines
################################################################################

sub check_parameters {
    my $phOptions = shift;
    
    ## make sure input parameters are provided
    if (! (defined $phOptions->{'gtffile'}) ) {
		pod2usage( -msg => $sHelpHeader, -exitval => 1);
	}
	
    ## handle some defaults
    $phOptions->{'cufflinks_bin_dir'} = BIN_DIR if (! (defined $phOptions->{'cufflinks_bin_dir'}) );
    $phOptions->{'prefix'} = "cuffcmp" if ((!(defined $hCmdLineOption{'prefix'})) || ($hCmdLineOption{'prefix'} !~ m/^$/));
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
