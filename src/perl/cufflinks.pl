#!/usr/bin/perl

################################################################################
### POD Documentation
################################################################################

=head1 NAME

cufflinks.pl - script to generate Cufflinks transcript detection from sorted SAM/BAM alignment file

=head1 SYNOPSIS

    cufflinks.pl --i alignment_file [--g annotation_file] [--o outdir] 
                 [--frag-len-mean mean_length] [--frag-len-std-dev stdev_length] 
                 [--label sample_id] [--min-isoform-fraction 0.0-1.0] [--pre-mrna-fraction 0.0-1.0] 
                 [--min-intron-length length] [--max-intron-length length] [--junc-alpha 0.0-1.0] 
                 [--small-anchor-fraction 0.0-1.0] [--min-frags-per-transfrag threshold] 
                 [--num-threads threads] [--cb cufflinks_bin_dir] [--args other_parameters] [--v] 

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --i <alignment_file>                   = /path/to/sorted-by-position alignment SAM or BAM file.

    --a <annotation_file>                  = /path/to/annotation file in GTF or GFF format.

    --o <output dir>                       = /path/to/output directory. Optional. [present working directory]

    --frag-len-mean <mean_length>          = the expected (mean) fragment length. Optional.

    --frag-len-std-dev <stdev_length>      = standard deviation for the distribution on fragment lengths. Optional.

    --label <sample_id>                    = prefix for transfrags in GTF format. Optional. [CUFF]

    --min-isoform-fraction <0.0-1.0>       = filters out transcripts with very low abundance. Optional. [0.1]

    --pre-mrna-fraction <0.0-1.0>          = parameter to filter out alignments that lie within the intronic intervals. Optional. [0.15]

    --min-intron-length <length>           = minimum intron size allowed. Optional. [50]

    --max-intron-length <length>           = maximum intron length. Optional. [300,000]

    --junc-alpha <0.0-1.0>                 = alpha value for the binomial test (spliced alignment filtration). Optional. [0.001]

    --small-anchor-fraction <0.0-1.0>      = Spliced reads with small anchors are filtered out. Optional. [0.09]

    --min-frags-per-transfrag <threshold>  = Assembled transfrags supported by minimum fragments. Optional. [10]

    --num-threads <# threads>              = Use # threads to align reads. Optional. [1]

    --cb <cufflinks_bin_dir>               = /path/to/cufflinks bin directory. Optional. [/usr/local/bin]

    --args <other_params>                  = additional Cufflinks parameters. Optional. Refer Cufflinks manual.

    --v                                    = generate runtime messages. Optional

=head1 DESCRIPTION

The script executes the cufflinks script from the Cufflinks RNA-Seq Analysis package.

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
            'infile|i=s', 'annotation|g=s', 'outdir|o=s', 'library-type=s', 
            'frag-len-mean|m=i', 'frag-len-std-dev|s=i', 'max-mle-iterations=i', 'max-bundle-frags=i', 
            'label|L=s', 'min-isoform-fraction|F=f', 'pre-mrna-fraction|j=f', 
            'min-intron-length|l=i', 'max-intron-length|I=i', 'junc-alpha|a=f', 
            'small-anchor-fraction|A=f', 'min-frags-per-transfrag|f=i', 'overhang-tolerance=i', 
            'max-bundle-length=i', 'trim-3-avgcov-thresh=i', 'trim-3-dropoff-frac=f', 'num-threads|p=i', 
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

my ($sOutDir, $sOutFile);
my ($sCmd, $sKey, $sPrefix);
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

($bDebug || $bVerbose) ? 
	print STDERR "\nExecuting Cufflinks Isoform Identification for input $hCmdLineOption{'infile'} ...\n" : ();

$sCmd  = $hCmdLineOption{'cufflinks_bin_dir'}."/cufflinks".
		 " --output-dir ".$sOutDir;

$sCmd .= " --GTF ".$hCmdLineOption{'annotation'} if (defined $hCmdLineOption{'annotation'});

foreach $sKey ( keys %hCmdLineOption) {
	next if (($sKey eq "infile") || ($sKey eq "annotation") || ($sKey eq "outdir") || 
			 ($sKey eq "cufflinks_bin_dir") || ($sKey eq "library-type") || ($sKey eq "args") || 
			 ($sKey eq "verbose") || ($sKey eq "debug") || ($sKey eq "help") || ($sKey eq "man") );
	
	$sCmd .= " --".$sKey." ".$hCmdLineOption{$sKey} if ((defined $hCmdLineOption{$sKey}) && ($hCmdLineOption{$sKey} !~ m/^$/i));
}

$sCmd .= " --library-type ".$hCmdLineOption{'library-type'} if ((defined $hCmdLineOption{'library-type'}) && ($hCmdLineOption{'library-type'} !~ m/^$/i));

$sCmd .= " ".$hCmdLineOption{'args'} if (defined $hCmdLineOption{'args'});

$sCmd .= " ".$hCmdLineOption{'infile'};

exec_command($sCmd);

($_, $_, $sPrefix) = File::Spec->splitpath($hCmdLineOption{'infile'});
$sPrefix =~ s/\.[bs]am$//;

if ( -e "$sOutDir/transcripts.gtf" ) {
	$sCmd = "mv $sOutDir/transcripts.gtf $sOutDir/$sPrefix.transcripts.gtf";
	exec_command($sCmd);
}

if ( -e "$sOutDir/isoforms.fpkm_tracking" ) {
	$sCmd = "mv $sOutDir/isoforms.fpkm_tracking $sOutDir/$sPrefix.isoforms.fpkm_tracking";
	exec_command($sCmd);
}

if ( -e "$sOutDir/genes.fpkm_tracking" ) {
	$sCmd = "mv $sOutDir/genes.fpkm_tracking $sOutDir/$sPrefix.genes.fpkm_tracking";
	exec_command($sCmd);
}

if ( -e "$sOutDir/skipped.gtf" ) {
	$sCmd = "mv $sOutDir/skipped.gtf $sOutDir/$sPrefix.skipped.gtf";
	exec_command($sCmd);
}

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'infile'} ... done\n" : ();


################################################################################
### Subroutines
################################################################################

sub check_parameters {
    my $phOptions = shift;
    
    ## make sure input alignment is provided
    if (! (defined $phOptions->{'infile'})) {
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

################################################################################
