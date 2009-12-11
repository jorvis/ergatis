#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

velvetg_post_processing.pl - Renames default velvetg output to store each file with a unique filename

=head1 SYNOPSIS

USAGE: ./velvetg_post_processing.pl --input_dir /path/to/velvet/output --output_dir /path/to/output --output_prefix OUTPUT_PREFIX --log /path/to/log --debug DEBUG

=head1 OPTIONS

B<--input_dir, -i>
	Path to the velvetg output directory

B<--output_dir, -o>
	Path to output directory where renamed velvetg output files will be stored.
		
B<--output_prefix, -p>
	Output prefix to rename velvetg output files to. Files will be renamed to PREFIX.*

B<--log, -l>
	Log file
	
B<--debug, -d>
	Debug level

=head1 DESCRIPTION

The velvetg software does not uniquely name its output files a process is needed to give the output files unique names for easier retrieval
in later pipelines. This script uses the output prefix provided by the user to achieve this.

=head1 INPUT

Standard output from velvetg

=head1 OUPTUT

Renamed velvetg output files using the --output_prefix flag as the filenmame.

=head1 CONTACT

	Cesar Arze
	carze@som.umaryland.edu
	
=cut

use strict;
use File::Copy;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Ergatis::Logger;

my %options = ();
my $result = GetOptions (\%options,
						 'input_dir|i=s',
						 'output_dir|o=s',
						 'output_prefix|p=s',
						 'log|l=s',
						 'debug|d=s',
						 'help|h') || pod2usage();

if ( $options{'help'} ) {
	pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
						  
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger( 'LOG_FILE' 	=>	$logfile,
								 'LOG_LEVEL'	=>	$options{'debug'} );
$logger = Ergatis::Logger::get_logger();

my $input_dir = $options{'input_dir'};
my $prefix = $options{'output_prefix'};
my $output_dir = $options{'output_dir'};

my $velvet_output = &get_velvet_output_files($input_dir);

copy($velvet_output->{contig}, "$output_dir/$prefix.contigs") or $logger->logdie("Could not copy contigs file: $!");
copy($velvet_output->{stat}, "$output_dir/$prefix.stats") or $logger->logdie("Could not copy stats file: $!");
copy($velvet_output->{graph}, "$output_dir/$prefix.graph") or $logger->logdie("Could not copy LastGraph file: $!");
copy($velvet_output->{amos}, "$output_dir/$prefix.afg") or $logger->logdie("Could not copy amos file: $!");

#####################################################
#													#
#				    SUBROUTINES						#
#													#
#####################################################

sub get_velvet_output_files {
	my $dir = shift;
	my $output_files = ();
	
	opendir(VELVET, $dir) or $logger->logdie("Could not open velvet output director $dir: $!");
	my @files = readdir(VELVET);
	@files = map { $dir . "/" . $_ } @files;
	closedir VELVET;
	
	foreach my $file (@files) {
		if ($file =~ /contigs\.fa/) {
			$output_files->{contig} = $file;
		} elsif ($file =~ /stats\.txt/) {
			$output_files->{stats} = $file;
		} elsif ($file =~ /LastGraph/) {
			$output_files->{graph} = $file;
		} elsif ($file =~ /velvet_asm\.afg/) {
			$output_files->{amos} = $file;
		} elsif ($file ne "." || $file ne "..") {
			$logger->warn("File $file is not a velvet output file");
		}
	}	
	
	return $output_files;
}

						
						
