#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

sffinfo.pl - Wrapper script for sffinfo to facilitate its execution via Ergatis.

=head1 SYNOPSIS

USAGE: sffinfo.pl --sffinfo_exec=/path/to/sffinfo --input=/path/to/input --output=/path/to/output --output_type=<TYPE OF OUTPUT>

=head1 OPTIONS

B<--sffinfo_exec, -e>
	Path to sffinfo executable

B<--config_opts, -c>
    Optional. Any command-line arguments to be passed to sffinfo.

B<--input, -i>
	Input file; An SFF file.

B<--output, -o>
	Output. Path to output not including a prefix i.e. /test/sff_test and not /test/sff_test.out
	
B<--output_type, -t>
	What types of output file to produce. Options include
		a -- Output just accessions
		s -- Output just the sequences
		q -- Ouput just the quality scores
		f -- Output just the flowgrams
		t -- Output the seq/qual/flow as tab-delimited lines
		
B<--log, -l>
	Optional. Log file.
	
B<--debug, -d>
	Optional. Debug level.
	
=head1 DESCRIPTION

A wrapper script that handles that prefixing of output depending on the output type selected.

=head1 INPUT

Sequence data in SFF format

=head1 OUTPUT

One of the selected output types (accessions, sequences, quality scores, etc.)

=head1 CONTACT

	Cesar Arze
	carze@som.umaryland.edu
	
=cut

use strict;
use warnings;
use Pod::Usage;
use Ergatis::Logger;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = 	GetOptions (\%options,
							'sffinfo_exec|e=s',
				 			'input|i=s',
                            'config_opts|c=s',
				 			'output|o=s',
				 			'output_type|t=s',
				 			'log|l=s',
				 			'debug|d=s') || pod2usage();

## Set logger
my $logfile = $options{'log'} || Ergatis::Logger::get_default_filename();
my $logger = new Ergatis::Logger( 'LOG_FILE' 	=>	$logfile,
        						  'LOG_LEVEL'	=>	$options{'debug'} );
$logger = Ergatis::Logger::get_logger();
			
my $output = $options{'output'};
my $input = $options{'input'};
my $type = $options{'output_type'};
my $sffinfo = $options{'sffinfo_exec'};
my $sffinfo_args = $options{'config_opts'};

## Check if we have multiple output types
my @types = ();
if ($type =~ /\,/) {
	@types = split(",", $type);
	
	foreach my $type (@types) {
		$type =~ s/^\s+|\s+$//g;
		my $prefix = &get_output_prefix($type);
		my $out = $output . $prefix;
		run_system_cmd("$sffinfo -$type $sffinfo_args $input > $out");
	}
} else {
	$type =~ s/^\s+|\s+$//g;
	my $prefix = &get_output_prefix($type);
	my $out = $output . $prefix;
	run_system_cmd("sffinfo -$type $sffinfo_args $input > $out");
}

#########################################################################
#                                                                       #
#                           SUBROUTINES                                 #
#                                                                       #
#########################################################################

## Because sffinfo can generate so many different types of file (fasta, qual, etc)
## we need to handle prefixing the output correctly.
sub get_output_prefix {
	my $output_type = shift;
	my $prefix;
	
	# TODO: Add support for other types?
	if ($output_type eq "a") {
		$prefix = ".acc";
	} elsif ($output_type eq "s") {
		$prefix = ".fsa";
	} elsif ($output_type eq "q") {
		$prefix = ".qual";
	} elsif ($output_type eq "f") {
		$prefix = ".flow";
	} elsif ($output_type eq "t") {
		$prefix = ".tab";
	} else {
		$logger->logdie("Invalid output type: $output_type");
	}
}

## Run system command and make sure return value indicates success.
sub run_system_cmd {
    my $cmd = shift;
    my $res = system($cmd);
        $res = $res >> 8;
        
    unless ($res == 0) {
        $logger->logdie("Could not run $cmd");
    }   
        
    return $res;
}				 
