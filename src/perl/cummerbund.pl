#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

################################################################################
### POD Documentation
################################################################################

=head1 NAME

    cummerbund.pl  -       Uses R cummerbund package to generate plots from cuffdiff 
                           analysis.

=head1 SYNOPSIS

    cummerbund.pl          --i <input_list> [--rs R script]  
                           [--o outdir] [--v] [--rb R_bin_dir]
                           [--rp R parameters]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS
    
    --i  <input>           = List of cuffdiff output file.

    --o <output dir>       = /path/to/output directory. Optional.[PWD]

    --rs [R script]        = /path/to/Rscript. Optional.

    --rb <R_bin_dir>       = /path/to/R binary. Optional. [/usr/local/bin]

    --rp <R_parameters>    = additonal R parameters. Optional. [--vanilla]

    --v                    = generate runtime messages. Optional

=head1 DESCRIPTION

Uses R cummerbund package to generate plots from cuffdiff analysis.

=head1 AUTHOR

 Priti Kumari
 Bioinformatics Software Engineer 
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
use FindBin qw($RealBin);

##############################################################################
### Constants
##############################################################################

use constant FALSE => 0;
use constant TRUE  => 1;

use constant R_BIN_DIR => '/usr/local/packages/R-2.15.2/bin';


use constant PROGRAM => eval { ($0 =~ m/(\w+\.pl)$/) ? $1 : $0 };

##############################################################################
### Globals
##############################################################################

my %hCmdLineOption = ();
my $sHelpHeader = "\nThis is ".PROGRAM."\n";

GetOptions( \%hCmdLineOption,
            'outdir|o=s', 'infile|i=s',
            'verbose|v',
            'debug','Rscript|rs=s',
            'help','r_bin_dir|rb=s', 'r_params|rp=s',
            'man') or pod2usage(2);

## display documentation
pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};
pod2usage( -msg => $sHelpHeader, -exitval => 1) if $hCmdLineOption{'help'};

## make sure everything passed was peachy
check_parameters(\%hCmdLineOption);

## Define variables
my ($sOutDir, $raw_dir, $result_dir, $cuff_dir);
my $bDebug   = (defined $hCmdLineOption{'debug'}) ? TRUE : FALSE;
my $bVerbose = (defined $hCmdLineOption{'verbose'}) ? TRUE : FALSE;
my ($fh,$dh, $path) ;
my ($temp, $cuff_file, $files, $sample, $new_file, $old_file, $sCmd, $sRscript);
my @arr;

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

$sRscript = File::Spec->rel2abs($hCmdLineOption{'Rscript'});

$sOutDir = File::Spec->canonpath($sOutDir);
$raw_dir = $sOutDir."/raw_data";
$result_dir = $sOutDir."/figures";

if (! -e $raw_dir) {
    mkdir($raw_dir) || die "ERROR! Cannot create output directory\n";
}


if (! -e $result_dir) {
    mkdir($result_dir) || die "ERROR! Cannot create output directory\n";
}

open($fh, "<$hCmdLineOption{'infile'}") or die "Cannot open the list file.";

while (<$fh>) {
    chomp($_);
    next if($_=~ /^\#/);
    $_ =  File::Spec->rel2abs($_);
    ($temp, $temp, $cuff_file) =  File::Spec->splitpath($_);
    @arr = split(/\./, $cuff_file) ;
    $cuff_dir = $sOutDir."/raw_data/".$arr[0].".".$arr[1] ;
    if (! -e $cuff_dir) {
	mkdir ($cuff_dir) || die "Error cannot create sample directory.";
    }
    $temp =  File::Spec->rel2abs($temp);
    opendir($dh, $temp) || die "Error cannot open the ergatis output directory.";
    while ($files = readdir($dh)) {
	$new_file = $files;
	$sample = $arr[0].".".$arr[1].".";
	$new_file =~ s/$sample//;
	$old_file = $temp."/".$files;
	$new_file = $cuff_dir."/".$new_file;
	symlink $old_file, $new_file;
    }
    close $dh;
    $path = $cuff_dir."/cuffData.db";
    if (-e $path) {unlink $path ;}
   

}
close $fh;


$sCmd  = $hCmdLineOption{'r_bin_dir'}."/R".
		 " ".$hCmdLineOption{'r_params'}.
                 " --args ".$raw_dir." " .$result_dir ;


$sCmd .= " < ".$sRscript;

#$sCmd = $sRscript." ".$raw_dir." ".$result_dir ;
exec_command($sCmd);

print "You have pretty figures now!!!!!!!\n";


################################################################################
### Subroutines
################################################################################

sub check_parameters {
    my $phOptions = shift;
    
    ## make sure input files provided
    if (! (defined $phOptions->{'infile'}) ) {
	pod2usage( -msg => $sHelpHeader, -exitval => 1);
    }
    if (! (defined $phOptions->{'Rscript'})) {
	$phOptions->{'Rscript'} = $RealBin."/cummerbund_isoform_analysis.R";
    }

    $phOptions->{'r_bin_dir'} = R_BIN_DIR if (! (defined $phOptions->{'r_bin_dir'}) );
    $phOptions->{'r_params'} = '--vanilla' if (! (defined $phOptions->{'r_params'}) );

}

sub set_environment {
	my $phOptions = shift;
	
	umask 0000;
	
	# add R path to global path
	$ENV{PATH} = $phOptions->{'r_bin_dir'}.":".$ENV{PATH};
}    
sub exec_command {
        my $sCmd = shift;

        if ((!(defined $sCmd)) || ($sCmd eq "")) {
                die "\nSubroutine::exec_command : ERROR! Incorrect command!\n";
        }

        my $nExitCode;

        print STDERR "\n$sCmd\n";
        $nExitCode = system("$sCmd");
        if ($nExitCode != 0) {
                die "\tERROR! Command Failed!\n\t$!\n";
        }
        print STDERR "\n";

        return;
}  

################################################################################
