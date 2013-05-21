#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

################################################################################
### POD Documentation
################################################################################

=head1 NAME

tophat.pl - script to execute reference based alignment for input sequence file(s).

=head1 SYNOPSIS

    tophat.pl --r1 mate_1_sequence_file(s) [--r2 mate_2_sequence_file(s)] 
              --bi bowtie_index_dir --p index_prefix [--o outdir] 
              [--tb tophat_bin_dir] [--bb bowtie_bin_dir] [--sb samtools_bin_dir] 
              [--mate-inner-dist] [--mate-std-dev] [--min-anchor-length] [--splice-mismatches] 
              [--min-intron-length] [--max-intron-length] [--max-insertion-length] [--max-deletion-length]
              [--num-threads] [--max-multihits] [--library-type]
              [--bowtie_mode] [--initial-read-mismatches] [--segment-mismatches] [--segment-length] 
              [--read-gap-length] [--read-edit-dist]
              [--min-coverage-intron] [--max-coverage-intron] [--min-segment-intron] [--max-segment-intron]
              [--raw-juncs] [--GTF] [--transcriptome-index] [--insertions] [--deletions]
              [--a other_parameters] [--v]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --r1 <mate_1_sequence_file(s)> = /path/to/input sequence file(s) for the first mates.
    
    --r2 <mate_2_sequence_file(s)> = /path/to/input sequence file(s) for the second mates. Optional.

    --bi <bowtie_index_dir>        = /path/to/bowtie_index directory.

    --p <index_prefix>             = prefix for index files.

    --o <output dir>               = /path/to/output directory. Optional. [present working directory]

    --tb <tophat_bin_dir>          = /path/to/tophat binary. Optional. [/usr/local/bin]

    --bb <bowtie_bin_dir>          = /path/to/bowtie binary. Optional. [/usr/local/bin]

    --sb <samtools_bin_dir>        = /path/to/samtools binary. Optional. [/usr/local/bin]

    Tophat parameters              = refer to Tophat User Manual accessible at 
                                     http://tophat.cbcb.umd.edu/manual.html for specific details and 
                                     defaults of the above mentioned tophat parameters. Optional.

    --a <other_parameters>         = additonal tophat parameters. Optional. Refer Tophat manual.

    --v                            = generate runtime messages. Optional

=head1 DESCRIPTION

The script executes the tophat script from the tophat software package.

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

use constant BOWTIE_BIN_DIR => '/usr/local/bin';
use constant SAMTOOLS_BIN_DIR => '/usr/local/bin';
use constant TOPHAT_BIN_DIR => '/usr/local/bin';

use constant VERSION => '1.0.0';
use constant PROGRAM => eval { ($0 =~ m/(\w+\.pl)$/) ? $1 : $0 };

##############################################################################
### Globals
##############################################################################

my %hCmdLineOption = ();
my $sHelpHeader = "\nThis is ".PROGRAM." version ".VERSION."\n";

GetOptions( \%hCmdLineOption,
            'seq1file|r1=s', 'seq2file|r2=s', 'bowtie_index_dir|bi=s', 'prefix|p=s',
            'outdir|o=s', 'tophat_bin_dir|tb=s', 'bowtie_bin_dir|bb=s', 'samtools_bin_dir|sb=s',
            'mate-inner-dist=i', 'mate-std-dev=i', 'min-anchor-length=i', 'splice-mismatches=i', 
            'min-intron-length=i', 'max-intron-length=i', 'max-insertion-length=i', 'max-deletion-length=i', 
            'num-threads=i', 'max-multihits=i', 'library-type=s', 
            'bowtie_mode=s', 'initial-read-mismatches=i', 'segment-mismatches=i', 'segment-length=i', 
            'read-gap-length=i', 'read-edit-dist=i',  
            'min-coverage-intron=i', 'max-coverage-intron=i', 'min-segment-intron=i', 'max-segment-intron=i', 
            'raw-juncs=s', 'GTF=s', 'transcriptome-index=s', 'insertions=s', 'deletions=s', 
            'args|a=s', 'verbose|v',
            'debug',
            'help',
            'man') or pod2usage(2);

## display documentation
pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};
pod2usage( -msg => $sHelpHeader, -exitval => 1) if $hCmdLineOption{'help'};

## make sure everything passed was peachy
check_parameters(\%hCmdLineOption);

my ($sOutDir, $sReadFile);
my ($sCmd, $sKey, $sFile, $sPrefix);
my $bDebug   = (defined $hCmdLineOption{'debug'}) ? TRUE : FALSE;
my $bVerbose = (defined $hCmdLineOption{'verbose'}) ? TRUE : FALSE;

################################################################################
### Main
################################################################################

if ($bDebug || $bVerbose) { 
	print STDERR "\nProcessing $hCmdLineOption{'seq1file'} ";
	print STDERR "& Processing $hCmdLineOption{'seq2file'} " if (defined $hCmdLineOption{'seq2file'});
	print STDERR "...\n";
}

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
	print STDERR "\nExecuting Tophat reference based alignment for input sequence file(s) ...\n" : ();

$sCmd  = $hCmdLineOption{'tophat_bin_dir'}."/tophat".
		 " --output-dir ".$sOutDir;

foreach $sKey ( keys %hCmdLineOption) {
	next if (($sKey eq "seq1file") || ($sKey eq "seq2file") || ($sKey eq "bowtie_index_dir") || ($sKey eq "prefix") ||
			 ($sKey eq "outdir") || ($sKey eq "tophat_bin_dir") || ($sKey eq "bowtie_bin_dir") ||
			 ($sKey eq "samtools_bin_dir") || ($sKey eq "library-type") || ($sKey eq "bowtie_mode") || ($sKey eq "args") || 
			 ($sKey eq "GTF") || ($sKey eq "transcriptome-index") || 
			 ($sKey eq "verbose") || ($sKey eq "debug") || ($sKey eq "help") || ($sKey eq "man") );
	
	$sCmd .= " --".$sKey." ".$hCmdLineOption{$sKey} if ((defined $hCmdLineOption{$sKey}) && ($hCmdLineOption{$sKey} !~ m/^$/i));
}

$sCmd .= " --library-type ".$hCmdLineOption{'library-type'} if ((defined $hCmdLineOption{'library-type'}) && ($hCmdLineOption{'library-type'} !~ m/^$/i));
$sCmd .= " --bowtie-n" if ((defined $hCmdLineOption{'bowtie_mode'}) && ($hCmdLineOption{'bowtie_mode'} =~ m/^n$/i));

if ((defined $hCmdLineOption{'GTF'}) && ($hCmdLineOption{'GTF'} !~ m/^$/i)) {
	$sCmd .= " --GTF ".$hCmdLineOption{'GTF'};
	$sCmd .= " --transcriptome-index ".$hCmdLineOption{'transcriptome-index'}
		if ((defined $hCmdLineOption{'transcriptome-index'}) && ($hCmdLineOption{'transcriptome-index'} !~ m/^$/i));
}

$sCmd .= " ".$hCmdLineOption{'args'} if (defined $hCmdLineOption{'args'});
$sCmd .= " ".$hCmdLineOption{'bowtie_index_dir'}."/".$hCmdLineOption{'prefix'};

$sReadFile = $hCmdLineOption{'seq1file'};
$sReadFile .= " ".$hCmdLineOption{'seq2file'} if ((defined $hCmdLineOption{'seq2file'}) && ($hCmdLineOption{'seq2file'} !~ m/^$/));
$sCmd .= " ".$sReadFile;

exec_command($sCmd);

$sFile = (split(/,/, $hCmdLineOption{'seq1file'}))[0];
($_, $_, $sPrefix) = File::Spec->splitpath($sFile);
$sPrefix =~ s/.1_1_sequence.*//;
$sPrefix =~ s/.sequence.*//;
$sPrefix =~ s/.fastq.*$//;
$sPrefix =~ s/.fq.*$//;

if ( -e "$sOutDir/accepted_hits.bam" ) {
	$sCmd = "mv $sOutDir/accepted_hits.bam $sOutDir/$sPrefix.accepted_hits.bam";
	exec_command($sCmd);
}

if ( -e "$sOutDir/accepted_hits.sam" ) {
	$sCmd = "mv $sOutDir/accepted_hits.sam $sOutDir/$sPrefix.accepted_hits.sam";
	exec_command($sCmd);
}

if ( -e "$sOutDir/deletions.bed" ) {
	$sCmd = "mv $sOutDir/deletions.bed $sOutDir/$sPrefix.deletions.bed";
	exec_command($sCmd);
}

if ( -e "$sOutDir/insertions.bed" ) {
	$sCmd = "mv $sOutDir/insertions.bed $sOutDir/$sPrefix.insertions.bed";
	exec_command($sCmd);
}

if ( -e "$sOutDir/junctions.bed" ) {
	$sCmd = "mv $sOutDir/junctions.bed $sOutDir/$sPrefix.junctions.bed";
	exec_command($sCmd);
}

if ( -e "$sOutDir/unmapped_left.fq.z" ) {
	$sCmd = "mv $sOutDir/unmapped_left.fq.z $sOutDir/$sPrefix.unmapped_left.fq.z";
	exec_command($sCmd);
}

if ( -e "$sOutDir/unmapped_right.fq.z" ) {
	$sCmd = "mv $sOutDir/unmapped_right.fq.z $sOutDir/$sPrefix.unmapped_right.fq.z";
	exec_command($sCmd);
}

if ($bDebug || $bVerbose) { 
	print STDERR "\nProcessing $hCmdLineOption{'seq1file'} ";
	print STDERR "& Processing $hCmdLineOption{'seq2file'} " if (defined $hCmdLineOption{'seq2file'});
	print STDERR "... done\n";
}


################################################################################
### Subroutines
################################################################################

sub check_parameters {
    my $phOptions = shift;
    
    ## make sure input fastx is provided
    if ( (!(defined $phOptions->{'seq1file'})) ||
    	 (!(defined $phOptions->{'bowtie_index_dir'})) ||
    	 (!(defined $phOptions->{'prefix'})) ) {
		pod2usage( -msg => $sHelpHeader, -exitval => 1);
	}
	
    ## handle some defaults
    $phOptions->{'tophat_bin_dir'} = TOPHAT_BIN_DIR if (! (defined $phOptions->{'tophat_bin_dir'}) );
    
    if (! (defined $phOptions->{'bowtie_bin_dir'}) ) {
    	$phOptions->{'bowtie_bin_dir'} = BOWTIE_BIN_DIR;
	}
	
	if (! (defined $phOptions->{'samtools_bin_dir'}) ) {
    	$phOptions->{'samtools_bin_dir'} = SAMTOOLS_BIN_DIR;
	}
	
	# set environment variables
		set_environment($phOptions);
}

sub set_environment {
	my $phOptions = shift;
	
	umask 0000;
	
	# adding bowtie executible path to user environment
	$ENV{PATH} = $phOptions->{'bowtie_bin_dir'}.":".$ENV{PATH};
	
	# adding samtools executible path to user environment
	$ENV{PATH} = $phOptions->{'samtools_bin_dir'}.":".$ENV{PATH};
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
