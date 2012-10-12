#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

################################################################################
### POD Documentation
################################################################################

=head1 NAME

generate_list_file.pl - script to generate a list of files in a directory that
                        satisfy the given regular expression.

=head1 SYNOPSIS

    generate_list_file.pl --d search_directory --r regex 
                          [--o output_listfile] [--gz] [--v]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --d <search directory>         = /path/to/directory to begin the search at

    --r <regex>                    = regular expression criterion to filter files

    --o <output_listfile>          = /path/to/output list of files found. Optional. [stdout]

    --se                           = do not die on error. Optional

    --gz                           = include gzipped files in the search output. Optional

    --v                            = generate runtime messages. Optional

=head1 DESCRIPTION

The script uses the Unix find command to search for specific files within a directory tree.

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
use File::Find;

## required for proper operation on NFS (as suggested by Joshua Orvis)
## see SF.net bug 2142533 - https://sourceforge.net/tracker2/?func=detail&aid=2142533&group_id=148765&atid=772583
$File::Find::dont_use_nlink = 1;

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
            'directory|d=s', 'regex|r=s', 'output_list|o=s', 
            'skip_error|se', 'gzip|gz', 'verbose|v',
            'debug',
            'help',
            'man') or pod2usage(2);

## display documentation
pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};
pod2usage( -msg => $sHelpHeader, -exitval => 1) if $hCmdLineOption{'help'};

## make sure everything passed was peachy
check_parameters(\%hCmdLineOption);

my (%hFiles);
my ($sFile, $nNumFiles);
my ($fpOUT, $nI);
my $bDebug   = (defined $hCmdLineOption{'debug'}) ? TRUE : FALSE;
my $bVerbose = (defined $hCmdLineOption{'verbose'}) ? TRUE : FALSE;

################################################################################
### Main
################################################################################

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'directory'} ...\n" : ();

($bDebug || $bVerbose) ? 
	print STDERR "\nSearching $hCmdLineOption{'directory'} for files matching $hCmdLineOption{'regex'} ...\n" : ();

search_regex_files(\%hCmdLineOption, \%hFiles);

$nNumFiles = keys %hFiles;

if (!($nNumFiles)) {
	if (defined $hCmdLineOption{'skip_error'}) {
		print STDERR "Error! No files matching $hCmdLineOption{'regex'} found in $hCmdLineOption{'directory'} !!!\n";
		exit;
	}
	else {
		die "Error! No files matching $hCmdLineOption{'regex'} found in $hCmdLineOption{'directory'} !!!\n";
	}
}

if (defined $hCmdLineOption{'output_list'}) {
	open($fpOUT, ">$hCmdLineOption{'output_list'}") or die "Error! Cannot open $hCmdLineOption{'output_list'} for writing !!!\n";
}

foreach $sFile ( sort keys %hFiles ) {
	if (defined $hCmdLineOption{'output_list'}) {
		print $fpOUT $hFiles{$sFile}."\n";
	}
	else {
		print STDOUT $sFile."\t".$hFiles{$sFile}."\n";
	}
}

if (defined $hCmdLineOption{'output_list'}) {
	close($fpOUT);
}

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'directory'} ... done\n" : ();


################################################################################
### Subroutines
################################################################################

sub check_parameters {
    my $phOptions = shift;
    
    ## make sure input directory and regex is provided
    if ((! (defined $phOptions->{'directory'}) ) || 
    	(! (defined $phOptions->{'regex'}) )) {
		pod2usage( -msg => $sHelpHeader, -exitval => 1);
	}
}

sub search_regex_files {
	my $phOptions 	= shift;
	my $phList		= shift;
	
	## make sure input directory and regex is provided
    if ((! (defined $phOptions) ) || 
    	(! (defined $phList) )) {
		die "Error! Insufficient parameters provided to subroutine search_regex_files !!!\n";
	}
	
	my $sDir = $phOptions->{'directory'};
	my $sRegex = $phOptions-> {'regex'};
	$sRegex .= "(\.gz)?" if (defined $phOptions->{'gzip'});
	
	my $sProcess = sub {
		return if ! -f;
        return if ! /$sRegex$/;
        $phList->{$_} = $File::Find::name;
    };
    
    find($sProcess, $sDir);
}


################################################################################
