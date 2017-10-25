#!/usr/bin/env perl

eval 'exec /usr/bin/env perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

################################################################################
### POD Documentation
################################################################################

=head1 NAME

    expression_plots.pl  - Creates plots for deseq/cuffdiff or rpkm analysis.             

=head1 SYNOPSIS

    expression_plots.pl    --i <infile> --a <analysis> [--lfc] [--fdr]   
                           [--gene_col] [--rpkm_col][--o outdir] [--v]
                           [--out_prefix] [--rs Rscript] [--rb R_bin_dir] 
                           [--rp R parameters]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS
    
    --i <infile>           = Sample file.<For cuffdiff/deseq/rpkm : list file>
                             
    
    --a <analysis>         = cuffdiff or deseq or rpkm.

    --lfc [lfc cutoff]     = For cuffdiff/deseq analysis. Default 1.

    --fdr [fdr cutoff]     = For cuffdiff/deseq analysis. Default 0.05.

    --gene_col [geneId col]= For rpkm analysis. Default 4.

    --rpkm_col [rpkm col]  = For rpkm analysis. Default 11.

    --o [output dir]       = /path/to/output directory. Optional.[PWD]

    --rs [Rscript]         = /path/to/Rscript. Optional.

    --rb [R_bin_dir]       = /path/to/R binary. Optional. [/usr/local/bin]

    --rp [R_parameters]    = additonal R parameters. Optional. [--vanilla]

    --out_prefix           = Output prefix. Optional

    --v                    = generate runtime messages. Optional

=head1 DESCRIPTION

Creates plots for deseq/cuffdiff or rpkm analysis.

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
            'lfc|l=s','fdr|f=s','analysis|a=s',
	    'gene_col|g=s','rpkm_col|r=s','verbose|v','debug',
            'Rscript|rs=s','out_prefix=s','help',
            'man','r_bin_dir|rb=s', 'r_params|rp=s') or pod2usage(2);

## display documentation
pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};
pod2usage( -msg => $sHelpHeader, -exitval => 1) if $hCmdLineOption{'help'};

## make sure everything passed was peachy
check_parameters(\%hCmdLineOption);

## Define variables
my ($sOutDir, $sRscript);
my $cmd;
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
$sRscript = File::Spec->rel2abs($hCmdLineOption{'Rscript'});



if( $hCmdLineOption{'analysis'} eq 'rpkm' ) {
    
    $cmd =  $hCmdLineOption{'r_bin_dir'}."/R".
	    " ".$hCmdLineOption{'r_params'}.
            " --args ".$hCmdLineOption{'infile'}." ".$hCmdLineOption{'gene_col'}." ".$hCmdLineOption{'rpkm_col'}." ".$sOutDir." ".$hCmdLineOption{'out_prefix'};


    $cmd .= " < ".$sRscript;


   #$cmd =  $sRscript." ".$hCmdLineOption{'infile'}." ".$hCmdLineOption{'gene_col'}." ".$hCmdLineOption{'rpkm_col'}." ".$sOutDir." ".$hCmdLineOption{'out_prefix'}; 
}
else {

    $cmd = $hCmdLineOption{'r_bin_dir'}."/R".
	    " ".$hCmdLineOption{'r_params'}.
            " --args ".$hCmdLineOption{'infile'}." ".$sOutDir." ".$hCmdLineOption{'lfc'}." ".$hCmdLineOption{'fdr'}." ".$hCmdLineOption{'analysis'};


    $cmd .= " < ".$sRscript;

#   $cmd =  $sRscript." ".$hCmdLineOption{'infile'}." ".$sOutDir." ".$hCmdLineOption{'lfc'}." ".$hCmdLineOption{'fdr'}." ".$hCmdLineOption{'analysis'}; 
}

exec_command($cmd);



################################################################################
### Subroutines
################################################################################

sub check_parameters {
    my $phOptions = shift;
    
    ## make sure input files provided
    if (! (defined $phOptions->{'infile'}) ) {
	pod2usage( -msg => $sHelpHeader, -exitval => 1);
    }
    if (! (defined $phOptions->{'out_prefix'}) ) {
	$phOptions->{'out_prefix'} = 'plot' ;
    }
    if (! (defined $phOptions->{'analysis'})) {
	pod2usage( -msg => $sHelpHeader, -exitval => 1);
    }
    if (! (defined $phOptions->{'lfc'})) {
	$phOptions->{'lfc'} = 1;
    }
    if (! (defined $phOptions->{'fdr'})) {
	$phOptions->{'fdr'} = 0.05;
    }
    if (! (defined $phOptions->{'Rscript'})) {
	$phOptions->{'Rscript'} = $RealBin."/expression_plots.R";
    }
    if ($phOptions->{'analysis'} eq 'rpkm' ) {
	$phOptions->{'Rscript'} = $RealBin."/rpkm_density_plot.r";
	if (! (defined $phOptions->{'gene_col'})) {
	    $phOptions->{'gene_col'} = 4;
	}
	if (! (defined $phOptions->{'rpkm_col'})) {
	    $phOptions->{'rpkm_col'} = 11 ;
	}
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
