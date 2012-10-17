#!/usr/bin/perl 

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

################################################################################
### POD Documentation
################################################################################

=head1 NAME

align_bowtie_stats.pl      -script to calculate alignment statistics for Bowtie output.

=head1 SYNOPSIS

    align_bowtie_stats.pl  --i <path to mapstat list file>  
                          [--o outdir] [--v]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --i <mapstat file list> = /path/to mapstats list file

    --o <output dir>        = /path/to/output directory. Optional. [present working directory]

    --v                     = generate runtime messages. Optional

=head1 DESCRIPTION

The script generates summary of alignment statistics from Bowtie output.

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

#use constant BIN_DIR => '/usr/local/bin';

#use constant VERSION => 'v0.10.0';
use constant PROGRAM => eval { ($0 =~ m/(\w+\.pl)$/) ? $1 : $0 };

##############################################################################
### Globals
##############################################################################

my %hCmdLineOption = ();
my $sHelpHeader = "\nThis is ".PROGRAM."\n";

GetOptions( \%hCmdLineOption,
            'infile|i=s','outdir|o=s', 
            'verbose|v',
            'debug',
            'help',
            'man') or pod2usage(2);

## display documentation
pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};
pod2usage( -msg => $sHelpHeader, -exitval => 1) if $hCmdLineOption{'help'};

## make sure everything passed was peachy
check_parameters(\%hCmdLineOption);

my ($sOutDir,$prefix);
my ($bowtie_list,$bowtie_file, $mapstat_list, $f1, $path, $key);
my %bam;
my ($cfile,$mfile,$lfile,$rfile,$readcount,@arr,$p_paired,$tot_reads,$p_mapped,$fout,$out_all);
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

open ($mapstat_list, "<$hCmdLineOption{infile}") or die "Error! Cannot open the mapstat list file";

while (<$mapstat_list>) {
    chomp($_);
    ($f1,$path,$prefix) = File::Spec->splitpath($_);
    @arr = split(/\./,$prefix);
    $prefix = $arr[0];
    if (!exists $bam{$prefix}) {
	$bam{$prefix}{"mapstats"} = $_;
	$_ =~ s/mapstats.txt/mapped_reads.count/;
	$bam{$prefix}{"count"} = $_;
    }
}


close $mapstat_list; 

open ($out_all, ">$sOutDir/All_Samples.txt") or die "Error Cannot open output file";


foreach $key (keys (%bam)) {
    open($cfile,"<$bam{$key}{'count'}") or die "Error! Cannot open the read count file.";
    open($mfile,"<$bam{$key}{'mapstats'}") or die "Error! Cannot open the mapstat file.";

    ###Reading read count file::
    $readcount=<$cfile>;
    chomp($readcount);

    ###Reading mapstat file to obtain properly paired an total reads:
    while (<$mfile>) {
	chomp ($_);
	if ($_ =~ m/properly paired/) {
	    $p_paired =(split(/[(:]/,$_))[1] ; 
	    $p_paired =~s/%//;
        }
        if ($_ =~ m/in total/) {
            $tot_reads =(split(/[(+]/,$_))[0] ; 
		     
	}
    }

    ###Percent mapped..
    $p_mapped = sprintf("%.2f",eval(($readcount/$tot_reads ) * 100)); 


    ###Write to Output Directory..
    open ($fout, ">$sOutDir/$key.txt");

    print $fout "\#Sample Id\t$key\n";
    print $fout "\#Total.reads\tMapped.reads\tPercent.Mapped\tPercent.Properly.Paired\n";
    print $fout "$tot_reads\t$readcount\t$p_mapped\t$p_paired\n";

    print $out_all "\#Sample Id\t$key\n";
    print $out_all "\#Total.reads\tMapped.reads\tPercent.Mapped\tPercent.Properly.Paired\n";
    print $out_all "$tot_reads\t$readcount\t$p_mapped\t$p_paired\n";
    
 
    close $fout;
    close $mfile;
    close $cfile;
     
}
close $out_all;
################################################################################
### Subroutines
################################################################################

sub check_parameters {
    my $phOptions = shift;
    
    ## make sure input files provided
    if (! (defined $phOptions->{'infile'}) ) {
		pod2usage( -msg => $sHelpHeader, -exitval => 1);
	}
  

}


################################################################################
