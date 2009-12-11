#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

mother.pl - Wrapper to facilitate the execution of mothur through ergatis.

=head1 SYNOPSIS

USAGE ./mothur.pl --mothur_exe=/path/to/mothur --input_file=/path/to/input --output_file=/path/to/output --args=<MOTHUR ARGS>

=head1 OPTIONS

B<--mothur_exe, -m>
	Path to the mothur executable.
	
B<--input_file, -i> 
	Input file that to be fed into mothur
	
B<--output_file, -o>
	Desired directory where mothur output should go

B<--args, -a> 
	Any command-line arguments that should be passed to mothur

B<--log, -l>
	OPTIONAL. Log file
	
=head1 DESCRIPTION

This wrapper is written to get around the funky way that mothur handles where it writes any output to. The output directory
is not user-controlled in the normal sense with output being written to the directory where the primary input file is house.
The primary input file is different depending on the mothur command invoked:

			* trim.seqs - fasta
			
=head1 INPUT

Input file for the mothur command to be invoked; this can be anything from a dist file to a FASTA file.

=head1 CONTACT

	Cesar Arze
	carze@som.umaryland.edu

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Ergatis::Logger;
use File::Copy;
use File::Basename;

my %options = ();
my $results = GetOptions (\%options,
						   'mothur_exe|m=s',
						   'input_file|i=s',
						   'output_file|o=s',
						   'args|a=s',
						   'log|l=s',
						   'help|h') || pod2usage();
						   
if ( $options{'help'} ) {
	pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}						   
						   					
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger( 'LOG_FILE' 	=>	$logfile,
							   	  'LOG_LEVEL'	=>	$options{'debug'} );							   	  
$logger = Ergatis::Logger::get_logger();

copy_input_files($options{input_file}, $options{output_file});
$options{'args'} =~ s/\,\s*\)$/\)/;

## Because mothur handles where output is generated to a bit wonky
## we not only need to copy our input file to the desired output
## directory but also change the working directory to the desired output
## directory
my ($name, $path) = basename($options{'output_file'});
chdir($path);

my $cmd = "$options{'mothur_exe'} \"$options{'args'}\""; 
run_system_cmd($cmd);

unlink($options{'output_file'});

#####################################################
#													#
#				    SUBROUTINES						#
#													#
#####################################################

sub copy_input_files {
	my ($file, $output_file) = @_;

	eval {
		copy($file, $output_file);
	};
	if ($@) {
		$logger->logdie("Could not copy input file to output directory: $@");
	}
}

sub run_system_cmd {
	my $cmd = shift;
	my $res = system($cmd);
	$res = $res >> 8;
	
	unless ($res == 0) {
		$logger->logdie("Could not run $cmd");
	}
	
	return $res;
}
