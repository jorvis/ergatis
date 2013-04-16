#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

################################################################################
### POD Documentation
################################################################################

=head1 NAME

bowtie.pl - script to execute reference based alignment for input sequence file(s).

=head1 SYNOPSIS

    bowtie.pl --r1 mate_1_sequence_file(s) [--r2 mate_2_sequence_file(s)] 
              --bi bowtie_index_dir --p index_prefix [--o outdir] [--bb bowtie_bin_dir] 
              [--mode v or n] [--num_mismatches mismatches_permitted] [--seedlen length] 
              [--minins min_insert_size] [--maxins max_insert_size] [--library-type fr|rf|ff] 
              [--k num_valid_alignments] [--m allowed_multihits] [--M allowed_multihits] 
              [--un /path/to/unmapped_reads_file] [--threads num_threads] [--file-type fastq|fasta] 
              [--a other_parameters] [--gz] [--v]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --r1 <mate_1_sequence_file(s)> = /path/to/input sequence file(s) for the first mates.
    
    --r2 <mate_2_sequence_file(s)> = /path/to/input sequence file(s) for the second mates. Optional.

    --bi <bowtie_index_dir>        = /path/to/bowtie_index directory.

    --p <index_prefix>             = prefix for index files.

    --o <output dir>               = /path/to/output directory. Optional. [present working directory]

    --bb <bowtie_bin_dir>          = /path/to/bowtie binary. Optional. [/usr/local/bin]

    --sb <samtools_bin_dir>        = /path/to/samtools binary. Optional. [/usr/local/bin]

    Bowtie parameters              = refer to Bowtie User Manual accessible at 
                                     http://bowtie-bio.sourceforge.net/manual.shtml for specific details and 
                                     defaults of the above mentioned tophat parameters. Optional.

    --a <other_parameters>         = additonal tophat parameters. Optional. Refer Bowtie manual.

    --gz                           = Required if input sequence files are gzipped files.

    --v                            = generate runtime messages. Optional

=head1 DESCRIPTION

The script executes the bowtie script from the bowtie software package.

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

use constant VERSION => '1.0.0';
use constant PROGRAM => eval { ($0 =~ m/(\w+\.pl)$/) ? $1 : $0 };

##############################################################################
### Globals
##############################################################################

my %hCmdLineOption = ();
my $sHelpHeader = "\nThis is ".PROGRAM." version ".VERSION."\n";

GetOptions( \%hCmdLineOption,
            'seq1file|r1=s', 'seq2file|r2=s', 'bowtie_index_dir|bi=s', 'prefix|p=s',
            'outdir|o=s', 'bowtie_bin_dir|bb=s', 'samtools_bin_dir|sb=s', 'mode=s', 'num_mismatches=i',
            'file-type|ft=s', 'seedlen|l=i', 'minins|I=i', 'maxins|X=i', 'library-type|lt=s',
            'k=i', 'm=i', 'M=i', 'un=s', 'threads=i', 
            'args|a=s', 'gzip|gz', 'verbose|v',
            'debug',
            'help',
            'man') or pod2usage(2);

## display documentation
pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};
pod2usage( -msg => $sHelpHeader, -exitval => 1) if $hCmdLineOption{'help'};

## make sure everything passed was peachy
check_parameters(\%hCmdLineOption);

## hash of parameters which need single '-' tag
my %hParams = ("k" => 1, "m" => 1, "M" => 1);

my ($sOutDir, $sReadFile);
my (@aSeqFiles, $sSeq1File, $sSeq2File);
my ($sCmd, $sKey, $sFile, $sPrefix, $nI);
my $bDebug   = (defined $hCmdLineOption{'debug'}) ? TRUE : FALSE;
my $bVerbose = (defined $hCmdLineOption{'verbose'}) ? TRUE : FALSE;

################################################################################
### Main
################################################################################

if ($bDebug || $bVerbose) { 
	print STDERR "\nProcessing $hCmdLineOption{'seq1file'} ";
	print STDERR "& Processing $hCmdLineOption{'seq2file'} " if ((defined $hCmdLineOption{'seq2file'}) && ($hCmdLineOption{'seq2file'} !~ m/^$/));
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
	print STDERR "\nExecuting Bowtie reference based alignment for input sequence file(s) ...\n" : ();

if (defined $hCmdLineOption{'gzip'}) {
	@aSeqFiles = split(/,/, $hCmdLineOption{'seq1file'});
	
	$sSeq1File = "";
	for ($nI = 0; $nI < @aSeqFiles; $nI++) {
		($_, $_, $sFile) = File::Spec->splitpath($aSeqFiles[$nI]);
		
		$sCmd = "ln -sf $aSeqFiles[$nI] $sOutDir/$sFile";
		exec_command($sCmd);
		
		$sFile =~ s/\.gz$//;
		$sCmd = "zcat $sOutDir/$sFile.gz > $sOutDir/$sFile";
		exec_command($sCmd);
		
		$sSeq1File .= (($sSeq1File =~ m/^$/) ? "$sOutDir/$sFile" : ",$sOutDir/$sFile");
	}
	
	if ((defined $hCmdLineOption{'seq2file'}) && ($hCmdLineOption{'seq2file'} !~ m/^$/)) {
		@aSeqFiles = split(/,/, $hCmdLineOption{'seq2file'});
		
		$sSeq2File = "";
		for ($nI = 0; $nI < @aSeqFiles; $nI++) {
			($_, $_, $sFile) = File::Spec->splitpath($aSeqFiles[$nI]);
			
			$sCmd = "ln -sf $aSeqFiles[$nI] $sOutDir/$sFile";
			exec_command($sCmd);
			
			$sFile =~ s/\.gz$//;
			$sCmd = "zcat $sOutDir/$sFile.gz > $sOutDir/$sFile";
			exec_command($sCmd);
			
			$sSeq2File .= (($sSeq2File =~ m/^$/) ? "$sOutDir/$sFile" : ",$sOutDir/$sFile");
		}
	}
}
else {
	$sSeq1File = $hCmdLineOption{'seq1file'};
	$sSeq2File = $hCmdLineOption{'seq2file'} if ((defined $hCmdLineOption{'seq2file'}) && ($hCmdLineOption{'seq2file'} !~ m/^$/));
}

$sCmd  = $hCmdLineOption{'bowtie_bin_dir'}."/bowtie";
$sCmd .= " -".$hCmdLineOption{'mode'}." ".$hCmdLineOption{'num_mismatches'};

foreach $sKey ( keys %hCmdLineOption) {
	next if (($sKey eq "seq1file") || ($sKey eq "seq2file") || ($sKey eq "bowtie_index_dir") || ($sKey eq "prefix") ||
			 ($sKey eq "outdir") || ($sKey eq "bowtie_bin_dir") || ($sKey eq "samtools_bin_dir") || ($sKey eq "mode") || ($sKey eq "num_mismatches") ||
			 ($sKey eq "library-type") || ($sKey eq "file-type") || ($sKey eq "un") || ($sKey eq "args") || 
			 ($sKey eq "gzip") || ($sKey eq "verbose") || ($sKey eq "debug") || ($sKey eq "help") || ($sKey eq "man") );
	
	if (defined $hParams{$sKey}) {
		$sCmd .= " -".$sKey." ".$hCmdLineOption{$sKey} if ($hCmdLineOption{$sKey} !~ m/^$/);
	}
	else {
		$sCmd .= " --".$sKey." ".$hCmdLineOption{$sKey} if ($hCmdLineOption{$sKey} !~ m/^$/);
	}
}

$sCmd .= " --un $sOutDir/".$hCmdLineOption{'un'} if (defined $hCmdLineOption{'un'});
$sCmd .= " --".$hCmdLineOption{'library-type'} if (defined $hCmdLineOption{'library-type'});
$sCmd .= " -q" if ($hCmdLineOption{'file-type'} =~ m/^fastq$/);
$sCmd .= " -f" if ($hCmdLineOption{'file-type'} =~ m/^fasta$/);
$sCmd .= " ".$hCmdLineOption{'args'} if (defined $hCmdLineOption{'args'});
$sCmd .= " ".$hCmdLineOption{'bowtie_index_dir'}."/".$hCmdLineOption{'prefix'};

if ((defined $hCmdLineOption{'seq2file'}) && ($hCmdLineOption{'seq2file'} !~ m/^$/) && ($hCmdLineOption{'seq2file'} !~ m/^""$/)) {
	$sReadFile = "-1 ".$sSeq1File.
				 " -2 ".$sSeq2File;
}
else {
	$sReadFile = $sSeq1File;
}

$sFile = (split(/,/, $hCmdLineOption{'seq1file'}))[0];
($_, $_, $sPrefix) = File::Spec->splitpath($sFile);
$sPrefix =~ s/.1_1_sequence.*//;
$sPrefix =~ s/.sequence.*//;
$sPrefix =~ s/.fastq.*$//;
$sPrefix =~ s/.fq.*$//;

$sCmd .= " ".$sReadFile." $sOutDir/$sPrefix.bowtie.sam";

exec_command($sCmd);

if ($bDebug || $bVerbose) { 
	print STDERR "\nProcessing $hCmdLineOption{'seq1file'} ";
	print STDERR "& Processing $hCmdLineOption{'seq2file'} " if ((defined $hCmdLineOption{'seq2file'}) && ($hCmdLineOption{'seq2file'} !~ m/^$/));
	print STDERR "... done\n";
}

if (defined $hCmdLineOption{'gzip'}) {
	@aSeqFiles = split(/,/, $sSeq1File);
	foreach $sFile (@aSeqFiles) {
		unlink $sFile or warn "Could not unlink $sFile: $!";
	}
	if ((defined $hCmdLineOption{'seq2file'}) && ($hCmdLineOption{'seq2file'} !~ m/^$/) && ($hCmdLineOption{'seq2file'} !~ m/^""$/)) {
		@aSeqFiles = split(/,/, $sSeq2File);
		foreach $sFile (@aSeqFiles) {
			unlink $sFile or warn "Could not unlink $sFile: $!";
		}
	}
}

($bDebug || $bVerbose) ? 
	print STDERR "\nConverting $sPrefix.bowtie.sam to $sPrefix.bowtie.bam ...\n" : ();

die "Error! $sOutDir/$sPrefix.bowtie.sam doesn't exist.\n" if (! -e "$sOutDir/$sPrefix.bowtie.sam");
sam2bam(\%hCmdLineOption, "$sOutDir/$sPrefix.bowtie.sam", $sOutDir);

($bDebug || $bVerbose) ? 
	print STDERR "\nConverting $sPrefix.bowtie.sam to $sPrefix.bowtie.bam ... done\n" : ();

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
    $phOptions->{'bowtie_bin_dir'} = BOWTIE_BIN_DIR if (! (defined $phOptions->{'bowtie_bin_dir'}) );
    $phOptions->{'samtools_bin_dir'} = SAMTOOLS_BIN_DIR if (! (defined $phOptions->{'samtools_bin_dir'}) );
    $phOptions->{'mode'} = 'v' if (! (defined $phOptions->{'mode'}) );
    $phOptions->{'num_mismatches'} = 2 if (! (defined $phOptions->{'num_mismatches'}) );
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
    $sCmd .= " -o ".$sOutFile." ".$sFile;
	
	exec_command($sCmd);
	
	($bDebug || $bVerbose) ? 
		print STDERR "\nConverting $sFile to $sOutFile ... done\n" : ();
	
	$sFormat = "BAM";
	
	return;
}

################################################################################
