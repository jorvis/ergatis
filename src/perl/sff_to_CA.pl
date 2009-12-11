#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

sff_to_CA.pl - run sffToCA from the wgs-assembler software package. This script acts as a wrapper to ensure input is handled correctly in the Ergatis framework.

=head1 SYNOPSIS

USAGE: sff_to_CA.pl

=head1 OPTIONS

B<--sff_file, -s> 
    A single sff file input.

B<--sff_file_list, -sl>
    A list of sff files to act as input to a single CA assembly run.
    
B<--sff_to_CA_exe, -exe>
    Path to the sffToCA executable.
    
B<--library, -lib>
    The UID of the library these reads are added to.
    
B<--clear, -c>
    Can be one of four options:
        * none 
        * soft
        * hard
        * chop

B<--linker, -lk>
    Search for linker, created mated reads.

B<--trim, -t>
    Trim back a specific amount of bases before the first N.

B<--insert_size, -is>
    Mates are on average i +- d bp apart.

B<--output, -o>
    Output path where fragment files will be placed.
    
B<--log, -l>
    Log file

B<--debug>
    Debug level
            
=head1 DESCRIPTION

Wrapper for the sffToCA executable. Ensures that input and parameters are correctly formatted for use in the Ergatis framework.

=head1 INPUT

Either a single or list of SFF files that will comprise a single assembly.

=head1 OUTPUT 

Fragment files that are ready to be supplied to the celera-assembler

=head1 CONTACT
    
    Cesar Arze
    carze@som.umaryland.edu

=cut                             

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Ergatis::Logger;

my %options = ();
my $results = GetOptions (\%options,
						  'sff_file|s=s',
						  'sff_file_list|sl=s',
						  'sff_to_CA_exe|exe=s',
						  'library|lib=s',
						  'linker|lk=s',
						  'trim|t=s',
						  'clear|c=s',
						  'insert_size|is=s',
						  'output|o=s',
						  'log|l=s',
						  'debug=s',
						  'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger( 'LOG_FILE' 	=>	$logfile,
								  'LOG_LEVEL'	=>	$options{'debug'} );
$logger = Ergatis::Logger::get_logger();
						  
my $SFF_TO_CA_EXE = $options{'sff_to_CA_exe'};						  
						  
my $sff_file = $options{'sff_file'};
my $sff_file_list = $options{'sff_file_list'};
my $library = $options{'library'};
my $linker = $options{'linker'};
my $trim = $options{'trim'};
my $clear = $options{'clear'};
my $insert_size = $options{'insert_size'};
my $output = $options{'output'};

my @exe_params = &setup_exe_params($output, $library, $linker, $trim, $clear, $insert_size);						  
my $params = join(' ', @exe_params);						  
my $input;		
						  
if ($options{'sff_file_list'}) {
	$input = &get_file_string($options{'sff_file_list'});
} elsif ($options{'sff_file'}) {
	$input = $sff_file;
}	

my $cmd = "$SFF_TO_CA_EXE $params $input";
run_system_cmd($cmd);

#####################################################
#													#
#				    SUBROUTINES						#
#													#
#####################################################

sub run_system_cmd {
	my $cmd = shift;
	my $res = system($cmd);
	$res = $res >> 8;
	
	unless ($res == 0) {
		$logger->logdie("Could not run $cmd");
	}
}

sub setup_exe_params {
	my ($out, $lib, $linker, $trim, $clear, $ins_size) = @_;
	my @params;
	
	if ($out) {
		push (@params, "-output $out");
	} else {
		$logger->logdie("Please provide a valid output");
	}
	
	if ($lib) {
		push (@params, "-libraryname $lib"); 
	} else {
		$logger->logdie("Please provide a library name");
	}
	
	if ($linker) {
		push (@params, "-linker $linker");
	}
	if ($trim) {
		push (@params, "-trim $trim");
	}
	if ($clear) {
		push (@params, "-clear $clear");
	}
	if ($ins_size) {
		push (@params, "-insertsize $ins_size");
	}
	
	return @params;
}

sub get_file_string {
	my $file = shift;
	my $sff_list_string;
	
	open (FILELIST, $file) or $logger->logdie("Could not open sff file list $file: $!");
	while (my $line = <FILELIST>) {
		chomp($line);
		verify_file($line);
		$sff_list_string .= " $line";
	}
	
	return $sff_list_string;
}				  

sub verify_file {
	my $file = shift;
	
	if ( ( defined($file) ) && (-e $file) && (-r $file) && (-s $file) ) {
		return 1;
	} else {
		
		if ( !defined($file) ) {
			$logger->logdie("FASTA file list $file was not defined");
		}
		
		if ( !-e $file ) {
			$logger->logdie("FASTA file list $file does not exist");
		}
		
		if ( !-r $file ) {
			$logger->logdie("FASTA file list $file does not have read permissions");
		}
		
		if ( !-s $file ) {
			$logger->logdie("FASTA file list $file has zero content");
		}
	}
}
