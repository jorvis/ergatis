#!/usr/bin/perl

###########################################################
# POD DOCUMENTATION                                       #
###########################################################
=head1 NAME

create_hmms.pl - Program to create Pfam and TIGRFAM HMM library and MLDB

=head1 SYNOPSIS

	create_hmm.pl --hmm_list <list_file> --tigrinfo_dir <TIGRFAMs_INFO> --tigrgo_link <TIGRFAMS_GO_LINK> --output_dir <output_dir> [--tigrrole_link <TIGRFAMS_ROLE_LINK> --log <log file> --debug <debug level> --help <usage>]
	
	parameters in [] are optional
    	do NOT type the carets when specifying options

=head1 OPTIONS
	
	--hmm_list	= /path/to/hmm_list_file.list. This is a file containing a list of different hmm
			  libraries to be merged. It contains the full paths of the untarred, ungzipped
			  recently downloaded HMM files. 

	--tigrinfo_dir	= /path/to/TIGRFAMs_INFO. This is the most current untarred, ungzipped TIGRFAMs_INFO
			  directory containing individual TIGRFAMs INFO files. [Latest : TIGRFAMs_10.1_INFO.tar.gz] 
 	
	--tigrgo_link	= /path/to/TIGRFAMS_GO_LINK. This is the most current TIGRFAMS_GO_LINK file
			  [Latest : 10.1_Release]
	
	--output_dir	= /path/to/output_dir. This is the output directory where created files will be stored.

      [ --tigrrole_link = /path/to/TIGRFAMS_ROLE_LINK. This is the most current TIGRFAMS_ROLE_LINK file
			  [Latest : 10.1_Release]. Optional

        --log		= /path/to/log_file.log. Log file. Optional

	--debug		= Debug level. Optional

	--help		= Help message, script usage. Optional
      ]

=head1 DESCRIPTION

This program creates Pfam and TIGRFAM lilbrary as well as MLDBM to be searched by polypeptide sequences. The current library is HMMER3 compatible
and uses HMMER3 suite for creating library files. Following steps are taken in library creation:
1.Merge the HMM library files mentioned in the list file supplied as arguement using merge_hmm_libraries.pl. 
  Duplicate HMMs are removed during the merge to create a single library of unique HMMs.
2.Convert the merged single library file to a binary format using hmmconvert (HMMER3 utility).
3.Create .h3f, .h3i, .h3m and .h3p files using hmmpress (HMMER3 utility) required by polypeptide sequences during search.
4.Create MLDBM from the merged library file using hmmlib_to_mldbm.pl.
5.Additional information from TIGRFAMs_INFO is incorporated in the above created MLDBM using tigrfam_info_mldbm.pl

=head1  CONTACT

Sonia Agrawal
sagrawal@som.umaryland.edu

=cut

use strict;
use File::Basename;
use FindBin qw($Bin);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use POSIX;

#################################################
# GLOBAL VARIABLES                              #
#################################################
my %options;
my ($results,$merged_lib,$merged_lib_bin,$mldbm_file, $logfh, $logfile, $log_dir);
my $hmmer_loc = "/usr/local/packages/hmmer-3.0/bin";
#################################################
# MAIN PROGRAM                                  #
#################################################
$results = GetOptions (\%options,
                'hmm_list|h=s',
                'tigrinfo_dir|i=s',
                'tigrgo_link|g=s',
		'tigrrole_link|r=s',
		'output_dir|o=s',
                'log|l=s',
                'debug|b=s',
                'help|h') || pod2usage();

## Display documentation
if( $options{'help'} ){
        pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## Make sure everything passed was correct
&check_parameters();

## Getting the log file
if (defined $logfile) {
	my ($log_base,$log_ext);
	($log_base,$log_dir,$log_ext) = fileparse($logfile,qr/\.[^.]*/);
	open($logfh, "> $logfile") or die "Could open $logfile log file for writing\n";
}

## Merged HMM library file
$merged_lib = $options{'output_dir'}."/coding_hmm.lib"; 

## Merged HMM library file in binary format
$merged_lib_bin = $merged_lib.".bin";

## MLDBM file from merged HMM library file
$mldbm_file = $merged_lib.".db";

## Merge HMM library files specified in the input list file
log_msg("INFO : Merging HMMs provided in $options{'hmm_list'}\n", "");
if (defined $logfile) {
	system("perl $Bin/merge_hmm_libraries.pl --input_list $options{'hmm_list'} --output_lib $merged_lib --log $log_dir/merge_hmm_libraries.log");
} else {
	system("perl $Bin/merge_hmm_libraries.pl --input_list $options{'hmm_list'} --output_lib $merged_lib");
}
if (-e $merged_lib) {
	log_msg("INFO : Created merged HMM library file $merged_lib\n", "");
} else {
	log_msg("ERROR : Could not create merged HMM library file $merged_lib. Check the log file $log_dir/merge_hmm_libraries.log for errors\n", "die");
}

## Convert the merged library file into binary format
log_msg("INFO : Converting merged library file $merged_lib into binary format $merged_lib_bin\n", "");
system("$hmmer_loc/hmmconvert -b $merged_lib > $merged_lib_bin");
if (-e $merged_lib_bin) {
	log_msg("INFO : Created merged library file in binary format $merged_lib_bin\n", "");
} else {
	log_msg("ERROR : Could not create merged HMM library file in binary format $merged_lib_bin.\n","die");
}

## Create .h3f, .h3i, .h3m and .h3p files
log_msg("INFO : Pressing merged binary library file $merged_lib_bin\n", "");
if (defined $logfile) {
	system("$hmmer_loc/hmmpress -f $merged_lib_bin > $logfile");
} else {
	system("$hmmer_loc/hmmpress -f $merged_lib_bin");
}
log_msg("\nINFO : Created .h3f, .h3i, .h3m and .h3p files\n", "");

## Create MLDBM from the merged library file
log_msg("INFO : Creating MLDBM from merged library file $merged_lib\n", "");
if (defined $logfile) {
	system("perl $Bin/hmmlib_to_mldbm.pl --hmm_file $merged_lib --output_file $mldbm_file --log $log_dir/hmmlib_to_mldbm.log");
} else {
	system("perl $Bin/hmmlib_to_mldbm.pl --hmm_file $merged_lib --output_file $mldbm_file");
}
if (-e $mldbm_file) {
	log_msg("INFO : Created MLDBM $mldbm_file\n", "");
} else {
	log_msg("ERROR : Could not create MLDBM file $mldbm_file. Check the log file $log_dir/hmmlib_to_mldbm.log for errors\n", "die");
}

## Add supplemental information to TIGRFAMs from TIGR_INFO
log_msg("INFO : Adding supplemental information to TIGRFAMs in MLDBM $mldbm_file\n", "");
if (defined $logfile) {
	system("perl $Bin/tigrfam_info_to_mldbm.pl --mldbm_file $mldbm_file --info_dir $options{'tigrinfo_dir'} --go_link $options{'tigrgo_link'} --role_link '$options{'tigrrole_link'}' --log $log_dir/tigrfam_info_to_mldbm.log");
} else {
	system("perl $Bin/tigrfam_info_to_mldbm.pl --mldbm_file $mldbm_file --info_dir $options{'tigrinfo_dir'} --go_link $options{'tigrgo_link'} --role_link '$options{'tigrrole_link'}'");
}
log_msg("INFO : Added supplemental information to TIGRFAMs in MLDBM $mldbm_file\nCompleted process\n", "");

close($logfh) if $logfh;

#########################################################
# SUBROUTINES                                           #
#########################################################

# Subroutine to check the supplied paramaters are correct
sub check_parameters {
	
	my @params = qw(hmm_list tigrinfo_dir tigrgo_link output_dir);
	foreach my $param_opt (@params) {
		## Check for all mandatory options
		unless ($options{$param_opt}) {
			pod2usage( {-exitval => 0 -verbose => 2, -output => \*STDERR} );
			log_msg("ERROR : Manadatory parameter $param_opt not passed\n","die");
		}
		## Check if the file exist
		unless (-e $options{$param_opt}) {
			log_msg("ERROR : Required file $param_opt $options{$param_opt} does not exist\n", "die");
		}
		## Check if the file is readable
		unless (-r $options{$param_opt}) {
			log_msg("ERROR : Required file $options{$param_opt} is not readable\n", "die");
		}
	}
	if(!($options{'tigrrole_link'})) {
		$options{'tigrrole_link'} = "";
		log_msg("INFO : TIGRFAMS_ROLE_LINK file not provided. Thus, TIGR_Role_Id information could not be incorporated into MLDBM file\n", "");
	}

	if($options{'log'}) {
		$logfile = $options{'log'};
	} 
}

sub log_msg {
	my ($msg,$action) = @_;
	if ($logfh) {
		print $logfh "$msg";
	} else {
		print "$msg\n";
	}
	if ($action eq "die") {
		die "";
	} 
}
__END__		
