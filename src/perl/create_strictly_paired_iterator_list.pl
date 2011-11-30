#!/usr/bin/perl

################################################################################
### POD Documentation
################################################################################

=head1 NAME

create_paired_iterator_list.pl - generates list file from various input sources of filenames.

=head1 SYNOPSIS

    create_paired_iterator_list.pl --input_file_list <path/to/listfile> --input_file1 <path/to/listfile> 
                                   --input_file2 <path/to/listfile> input_directory <path/to/directory of input files> 
                                   --input_directory_extension <input file extension> --output_iter_list <path/to/output_iter_file>
                                   [--log </path/to/logfile>] [--timestamp] [--v]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --input_file_list <path/to/listfile>               = /path/to/listfile for paired input files.
    
    --input_file1 <path/to/listfile>                   = /path/to/input file 1.

    --input_file2 <path/to/listfile>                   = /path/to/input file 2.

    --input_directory <path/to/file directory>         = /path/to/directory containing input files.

    --input_directory_extension <input file extension> = input file extension. Required with --input_directory.

    --output_iter_list <path/to/output_iter_file>      = /path/to/output iteration file.

    --log </path/to/logfile>                           = /path/to/output log file. Optional.

    --timestamp                                        = sybase timestamp. Optional.

    --v                                                = generate runtime messages. Optional.

=head1   DESCRIPTION
	
	Modified from ./create_file_iterator_list.pl.
	
    This script is used to accept a selection of inputs from either an input list or multiple input 
    fields. These inputs will then be distributed randomly across a certain number of groups.

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
            'input_file_list|l=s', 'input_file1|f1=s', 'input_file2|f2=s',
            'input_directory|d=s', 'input_directory_extension|e=s', 'output_iter_list|o=s', 
            'verbose|v', 
            'help',
            'man') or pod2usage(2);

## display documentation
pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};
pod2usage( -msg => $sHelpHeader, -exitval => 1) if $hCmdLineOption{'help'};

## play nicely
umask(0000);

## make sure everything passed was peachy
check_parameters(\%hCmdLineOption);

my (@aInputFiles1, @aInputFiles2);
my ($sInFile1, $sInFile2, $sOutFile);
my ($fpLST, $fpOUT);
my ($nI, $sPrefix);
my $bDebug   = (defined $hCmdLineOption{'debug'}) ? TRUE : FALSE;
my $bVerbose = (defined $hCmdLineOption{'verbose'}) ? TRUE : FALSE;

################################################################################
### Main
################################################################################

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing input parameters to generate list file ...\n" : ();

$sOutFile = $hCmdLineOption{'output_iter_list'};

open($fpOUT, ">$sOutFile") or die "Error: Cannot open $sOutFile for writing:\n$!\n";

print $fpOUT "\$;I_FILE_BASE\$;\t\$;INPUT_FILE_1\$;\t\$;INPUT_FILE_2\$;\n";

if (defined $hCmdLineOption{'input_file_list'}) {
	open($fpLST, "<$hCmdLineOption{'input_file_list'}")
		or die "Error: Cannot open $hCmdLineOption{'input_file_list'} for reading:\n$!\n";
	
	while (<$fpLST>) {
		$_ =~ s/\s+$//;
		
		next if ($_ =~ m/^#/);
		
		($sInFile1, $sInFile2) = split(/\t/, $_);
		
		if ((defined $sInFile1) && ($sInFile1 !~ m/^$/)) {
			($_, $_, $sPrefix) = File::Spec->splitpath($sInFile1);
			
			print $fpOUT "$sPrefix\t$sInFile1";
			
			print $fpOUT "\t$sInFile2" if ((defined $sInFile2) && ($sInFile2 !~ m/^$/));
		}
		
		print $fpOUT "\n";
	}
	
	close($fpLST);
}

if ((defined $hCmdLineOption{'input_file1'}) && (!(defined $hCmdLineOption{'input_file2'}))) {
	@aInputFiles1 = split(/,/, $hCmdLineOption{'input_file1'});
	
	for ($nI = 0; $nI<@aInputFiles1; $nI++) {
		($_, $_, $sPrefix) = File::Spec->splitpath($aInputFiles1[$nI]);
		print $fpOUT "$sPrefix\t$aInputFiles1[$nI]\n";
	}
}

if ((defined $hCmdLineOption{'input_file1'}) && (defined $hCmdLineOption{'input_file2'})) {
	@aInputFiles1 = split(/,/, $hCmdLineOption{'input_file1'});
	@aInputFiles2 = split(/,/, $hCmdLineOption{'input_file2'});
	
	die "Error: The number of Input files 2 (".@aInputFiles2.") needs to be equal to ".
		"the number of Input files 1 (".@aInputFiles1.") !!!!!\n" if (@aInputFiles2 != @aInputFiles1);
	
	for ($nI = 0; $nI<@aInputFiles1; $nI++) {
		($_, $_, $sPrefix) = File::Spec->splitpath($aInputFiles1[$nI]);
		print $fpOUT "$sPrefix\t$aInputFiles1[$nI]\t$aInputFiles2[$nI]\n";
	}
}

if ((defined $hCmdLineOption{'input_directory'}) && (defined $hCmdLineOption{'input_directory_extension'})) {
	die "Error: Cannot find directory $hCmdLineOption{'input_directory'} !!!!!\n" if (!(-e $hCmdLineOption{'input_directory'}));
	
	@aInputFiles1 = glob("$hCmdLineOption{'input_directory'}/*.$hCmdLineOption{'input_directory_extension'}");
	
	for ($nI = 0; $nI<@aInputFiles1; $nI++) {
		($_, $_, $sPrefix) = File::Spec->splitpath($aInputFiles1[$nI]);
		print $fpOUT "$sPrefix\t$aInputFiles1[$nI]\n";
	}
}

close($fpOUT);

exit;

################################################################################
### Subroutines
################################################################################

sub check_parameters {
    my $phOptions = shift;
    
    ## make sure input file(s) is provided
    my $bFlag = 0;
    
    undef $phOptions->{'input_file_list'} if ($phOptions->{'input_file_list'} =~ /^$/);
    undef $phOptions->{'input_file1'} if ($phOptions->{'input_file1'} =~ /^$/);
    undef $phOptions->{'input_file2'} if ($phOptions->{'input_file2'} =~ /^$/);
    undef $phOptions->{'input_directory'} if ($phOptions->{'input_directory'} =~ /^$/);
    undef $phOptions->{'input_directory_extension'} if ($phOptions->{'input_directory_extension'} =~ /^$/);
    
    $bFlag = 1 if (defined $phOptions->{'input_file_list'});
    
    $bFlag = 1 if (defined $phOptions->{'input_file1'});
    
    $bFlag = 1 if ((defined $phOptions->{'input_file2'}) &&
    			   (defined $phOptions->{'input_file1'}));
    
    $bFlag = 1 if ((defined $phOptions->{'input_directory'}) &&
    			   (defined $phOptions->{'input_directory_extension'}));
    
    pod2usage( -msg => $sHelpHeader, -exitval => 1) if (!($bFlag));
}

################################################################################
