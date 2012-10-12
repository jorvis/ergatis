#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

################################################################################
### POD Documentation
################################################################################

=head1 NAME

create_paired_list_file.pl - script to generate paired list file from two list files.

=head1 SYNOPSIS

    create_paired_list_file.pl --l1 listfile1 --l2 listfile2 [--s sample_info_file [--o outdir] [--v]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --l1 <listfile1>               = /path/to/list file for input file 1s.

    --l2 <listfile2>               = /path/to/list file for input file 2s.

    --s <sample_info_file>         = /path/to/samples file with information on all samples to be analyzed.

    --o <output dir>               = /path/to/output directory. Optional. [present working directory]

    --v                            = generate runtime messages. Optional

=head1 DESCRIPTION

The script combines two separate list files of input files into a list file of paired files.

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

use constant VERSION => '1.0.0';
use constant PROGRAM => eval { ($0 =~ m/(\w+\.pl)$/) ? $1 : $0 };

##############################################################################
### Globals
##############################################################################

my %hCmdLineOption = ();
my $sHelpHeader = "\nThis is ".PROGRAM." version ".VERSION."\n";

GetOptions( \%hCmdLineOption,
            'listfile1|l1=s', 'listfile2|l2=s', 'samplefile|s=s', 
            'outdir|o=s', 'verbose|v',
            'debug',
            'help',
            'man') or pod2usage(2);

## display documentation
pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};
pod2usage( -msg => $sHelpHeader, -exitval => 1) if $hCmdLineOption{'help'};

## make sure everything passed was peachy
check_parameters(\%hCmdLineOption);

my %hSamples = ();
my ($sOutDir, $sOutFile);
my (@aSampleNames, $sSampleName, $sGroupName, $sRead1File, $sRead2File, $sReadsFile, $sFile, $nI);
my ($fpSMPL, $fpIN1, $fpIN2, $fpOUT);
my $bDebug   = (defined $hCmdLineOption{'debug'}) ? TRUE : FALSE;
my $bVerbose = (defined $hCmdLineOption{'verbose'}) ? TRUE : FALSE;

################################################################################
### Main
################################################################################

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'listfile1'} and $hCmdLineOption{'listfile2'} ...\n" : ();

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

$sOutFile = $sOutDir."/paired_input_file.list";

($bDebug || $bVerbose) ? 
	print STDERR "\nCombining input list files into paired list file $sOutFile ...\n" : ();

if ((defined $hCmdLineOption{'samplefile'}) && ($hCmdLineOption{'samplefile'} =~ m/^$/)) {
	open ($fpSMPL, "<$hCmdLineOption{'samplefile'}") or die "Error! Cannot open $hCmdLineOption{'samplefile'} for reading !!!\n";
	
	while (<$fpSMPL>) {
		$_ =~ s/\s+$//;
		
		next if ($_ =~ /^#/);
		
		next if ($_ =~ /^$/);
		
		($sSampleName, $sGroupName, $sRead1File, $sRead2File) = split(/\t/, $_);
		
		die "Error! Missing Reads 1 File !!!\n" if (! (defined $sRead1File));
		
		my @aRead1 = ();
		my @aRead2 = ();
		
		$hSamples{$sSampleName} = [\@aRead1, \@aRead2];
		
		push @aSampleNames, $sSampleName;
	}
	
	close($fpSMPL);
	
	open ($fpIN1, "<$hCmdLineOption{'listfile1'}") or die "Error! Cannot open $hCmdLineOption{'listfile1'} for reading !!!\n";
	while (<$fpIN1>) {
		$_ =~ s/\s+$//;
		
		next if ($_ =~ /^#/);
		
		next if ($_ =~ /^$/);
		
		$sReadsFile = $_;
		
		($_, $_, $sFile) = File::Spec->splitpath($sReadsFile);
		
		foreach $sSampleName (sort keys %hSamples) {
			if ($sFile =~ m/^$sSampleName\_0/) {
				push @{$hSamples{$sSampleName}[0]}, $sReadsFile;
				last;
			}
			elsif ($sFile =~ m/^$sSampleName\_1/) {
				push @{$hSamples{$sSampleName}[0]}, $sReadsFile;
				last;
			}
			elsif ($sFile =~ m/^$sSampleName\_2/) {
				push @{$hSamples{$sSampleName}[1]}, $sReadsFile;
				last;
			}
		}
	}
	close($fpIN1);
	
	open ($fpIN2, "<$hCmdLineOption{'listfile2'}") or die "Error! Cannot open $hCmdLineOption{'listfile2'} for reading !!!\n";
	while (<$fpIN2>) {
		$_ =~ s/\s+$//;
		
		next if ($_ =~ /^#/);
		
		next if ($_ =~ /^$/);
		
		$sReadsFile = $_;
		
		($_, $_, $sFile) = File::Spec->splitpath($sReadsFile);
		
		foreach $sSampleName (sort keys %hSamples) {
			if ($sFile =~ m/^$sSampleName\_0/) {
				push @{$hSamples{$sSampleName}[0]}, $sReadsFile;
				last;
			}
			elsif ($sFile =~ m/^$sSampleName\_1/) {
				push @{$hSamples{$sSampleName}[0]}, $sReadsFile;
				last;
			}
			elsif ($sFile =~ m/^$sSampleName\_2/) {
				push @{$hSamples{$sSampleName}[1]}, $sReadsFile;
				last;
			}
		}
	}
	close($fpIN2);
	
	open ($fpOUT, ">$sOutFile") or die "Error! Cannot open $sOutFile for writing !!!\n";
	
	for ($nI = 0; $nI < @aSampleNames; $nI++) {
		$sSampleName = $aSampleNames[$nI];
		die "Error! Reads file not specified for $sSampleName!!!\n" if (@{$hSamples{$sSampleName}[0]} == 0);
		
		my @aRead1 = sort {$a eq $b} @{$hSamples{$sSampleName}[0]};
		my @aRead2 = sort {$a eq $b} @{$hSamples{$sSampleName}[1]};
			
		print $fpOUT "".join(",", @aRead1)."\t".
					 "".join(",", @aRead2)."\n";
	}
	
	close($fpOUT);
}
else {
	open ($fpIN1, "<$hCmdLineOption{'listfile1'}") or die "Error! Cannot open $hCmdLineOption{'listfile1'} for reading !!!\n";
	open ($fpIN2, "<$hCmdLineOption{'listfile2'}") or die "Error! Cannot open $hCmdLineOption{'listfile2'} for reading !!!\n";
	open ($fpOUT, ">$sOutFile") or die "Error! Cannot open $sOutFile for writing !!!\n";
	
	while (<$fpIN1>) {
		$_ =~ s/\s+$//;
		print $fpOUT "$_\t";
		
		$_ = <$fpIN2>;
		$_ =~ s/\s+$//;
		print $fpOUT "$_\n";
	}
	
	close($fpIN1);
	close($fpIN2);
	close($fpOUT);
}

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'listfile1'} and $hCmdLineOption{'listfile2'} ... done\n" : ();


################################################################################
### Subroutines
################################################################################

sub check_parameters {
    my $phOptions = shift;
    
    ## make sure input fastx is provided
    if ((! (defined $phOptions->{'listfile1'}) ) || 
    	(! (defined $phOptions->{'listfile2'}) )) {
		pod2usage( -msg => $sHelpHeader, -exitval => 1);
	}
}

################################################################################
