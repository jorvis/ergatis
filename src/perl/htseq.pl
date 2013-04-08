#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

################################################################################
### POD Documentation
################################################################################

=head1 NAME

htseq.pl - script to generate HTSeq read counts from sorted-by-name SAM alignment file

=head1 SYNOPSIS

    htseq.pl --i input_SAM_file --a annotation_file [--o outdir] 
    		 [--m union, intersection-strict, or intersection-nonempty] 
    		 [--t feature type from annotation file to count over] 
    		 [--rq minimum quality of alignment read] [--s strand-specific] 
    		 [--id attribute id to be used as feature id] [--p other_parameters] [--v] 

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --i <input_SAM_file>           = /path/to/sorted-by-name alignment SAM file.

    --a <annotation_file>          = /path/to/annotation file in GTF or GFF format.

    --o <output dir>               = /path/to/output directory. Optional. [present working directory]

    --m <count_mode>               = Mode to handle reads overlapping more than one feature.
                                     Possible values are union, intersection-strict or intersection-nonempty. Optional. [union]

    --t <feature_type>             = feature type (3rd column in GFF file) to be used. Optional. [exon]

    --rq <minimum read quality>    = skip all reads with alignment quality lower than the given minimum value. Optional. [0]

    --s <strand-specific>          = whether the data is from a strand-specific assay. Optional. [yes]

    --id <attribute id>            = GFF attribute ID to be used as feature ID. Optional. [gene_id]

    --p <other_params>             = additonal htseq parameters. Optional. Refer HTSeq manual.

    --v                            = generate runtime messages. Optional

=head1 DESCRIPTION

The script executes the htseq-count (count.py python) script from the htseq python package.

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

use constant PYTHON_BIN_DIR => '/usr/local/packages/Python-2.6.4/bin';
use constant PYTHON_LIB_DIR => '/usr/local/packages/pythonlib/2.6';

use constant VERSION => '1.0.0';
use constant PROGRAM => eval { ($0 =~ m/(\w+\.pl)$/) ? $1 : $0 };

##############################################################################
### Globals
##############################################################################

my %hCmdLineOption = ();
my $sHelpHeader = "\nThis is ".PROGRAM." version ".VERSION."\n";

GetOptions( \%hCmdLineOption,
            'infile|i=s', 'annotation|a=s', 'outdir|o=s', 'mode|m=s', 'type|t=s',
            'minqual|rq=i', 'stranded|s=s', 'idattr|id=s', 'params|p=s',
            'python_bin_dir|pb=s', 'python_lib_dir|pl=s', 'verbose|v',
            'debug',
            'help',
            'man') or pod2usage(2);

## display documentation
pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};
pod2usage( -msg => $sHelpHeader, -exitval => 1) if $hCmdLineOption{'help'};

## make sure everything passed was peachy
check_parameters(\%hCmdLineOption);

my ($sOutDir, $sOutFile);
my ($sCmd);
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

($_, $_, $sOutFile) = File::Spec->splitpath($hCmdLineOption{'infile'});
$sOutFile =~ s/.sam$//;
$sOutFile .= ".".$hCmdLineOption{'type'}.".counts";

($bDebug || $bVerbose) ? 
	print STDERR "\nGenerating $hCmdLineOption{'type'} counts for $hCmdLineOption{'infile'} ...\n" : ();

$sCmd  = $hCmdLineOption{'python_bin_dir'}."/python -m HTSeq.scripts.count".
		 " -m ".$hCmdLineOption{'mode'}.
		 " -t ".$hCmdLineOption{'type'}.
		 " -i ".$hCmdLineOption{'idattr'}.
		 " -s ".$hCmdLineOption{'stranded'};
$sCmd .= " ".$hCmdLineOption{'params'} if (defined $hCmdLineOption{'params'});
$sCmd .= " ".$hCmdLineOption{'infile'}.
		 " ".$hCmdLineOption{'annotation'}.
		 " > ".$sOutDir."/".$sOutFile;

exec_command($sCmd);

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'infile'} ... done\n" : ();


################################################################################
### Subroutines
################################################################################

sub check_parameters {
    my $phOptions = shift;
    
    ## make sure input fastx is provided
    if ((! (defined $phOptions->{'infile'}) ) ||
    	(! (defined $phOptions->{'annotation'}) )) {
		pod2usage( -msg => $sHelpHeader, -exitval => 1);
	}
	
    ## handle some defaults
    $phOptions->{'mode'} = 'union' if (! (defined $phOptions->{'mode'}) );
    $phOptions->{'type'} = 'exon' if (! (defined $phOptions->{'type'}) );
    $phOptions->{'minqual'} = 0 if (! (defined $phOptions->{'minqual'}) );
    $phOptions->{'stranded'} = 'yes' if (! (defined $phOptions->{'stranded'}) );
    $phOptions->{'idattr'} = 'gene_id' if (! (defined $phOptions->{'idattr'}) );
    
    $phOptions->{'python_bin_dir'} = PYTHON_BIN_DIR if (! (defined $phOptions->{'python_bin_dir'}) );
    $phOptions->{'python_lib_dir'} = PYTHON_LIB_DIR if (! (defined $phOptions->{'python_lib_dir'}) );
    
    # set environment variables
	set_environment($phOptions);
}

sub set_environment {
	my $phOptions = shift;
	
	umask 0000;
	
	# add python path to global path
	$ENV{PATH} = $phOptions->{'python_bin_dir'}.":".$phOptions->{'python_lib_dir'}.":".$ENV{PATH};
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
