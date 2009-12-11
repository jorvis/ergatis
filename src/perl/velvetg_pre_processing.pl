#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

velvetg_pre_processing.pl - Renames velveth output to default naming schema.

=head1 SYNOPSIS

USAGE: ./velvetg_pre_processing.pl --input_dir /path/to/velvet/output --output_dir /path/to/output --prefix OUTPUT_PREFIX --log /path/to/log --debug DEBUG

=head1 OPTIONS

B<--input_dir, -i>
	Path to the velvet output directory containing the Logs, Roadmap, and Sequences files.

B<--output_dir, -o>
	Path to output directory where renamed velvet output files will be stored.
		
B<--prefix, -p>
	Output prefix to rename velvet output files to. Files will be named OUTPUT_PREFIX.sequences, OUTPUT_PREFIX.roadmaps, OUTPUT_PREFIX.log

B<--log, -l>
	Log file
	
B<--debug, -d>
	Debug level

=head1 DESCRIPTION

The velvet software does not uniquely name its output files a process is needed to give the output files unique names for easier retrieval
in later pipelines. This script takes the renamed unique output names and renames them back to the default schema so that the velvetg software
can run correctly.

=head1 INPUT

Standard output from velvet uniquely renamed; *.sequences, *.roadmaps, *.log

=head1 OUPTUT

Renamed velvet output files using the default naming schema -- Sequences, Roadmaps, Log

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
						 'prefix|p=s',
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
my $prefix = $options{'prefix'};
my $output_dir = $options{'output_dir'};

my $processed_velvet_output = &get_velvet_output_files($input_dir);

copy($processed_velvet_output->{sequence}, "$output_dir/Sequences") or $logger->logdie("Could not copy Sequence file: $!");
copy($processed_velvet_output->{roadmap}, "$output_dir/Roadmaps") or $logger->logdie("Could not copy Roadmaps file: $!");
copy($processed_velvet_output->{log}, "$output_dir/Log") or $logger->logdie("Could not copy Log file: $!");

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
		if ($file =~ /$prefix\.sequences/) {
			$output_files->{sequence} = $file;
		} elsif ($file =~ /$prefix\.roadmaps/) {
			$output_files->{roadmap} = $file;
		} elsif ($file =~ /$prefix\.log/) {
			$output_files->{log} = $file;
		} elsif ($file ne "." || $file ne "..") {
			$logger->warn("File $file is not a velvet output file");
		}
	}	
	
	return $output_files;
}
