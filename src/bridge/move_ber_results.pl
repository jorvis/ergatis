#!/usr/local/bin/perl

#########################################################################
## This script will move the BER result files to be used by Manatee,
## changing the Workflow style IDs into legacy style IDs along the way.
## It also fixes offsets in the BTAB files as needed.  Here's a sample
## usage:
##
## $progname -D /usr/local//usr/local/annotation/EHA2/output_repository/ber/1984_AllGroup.niaa/ -d EHA2
##
## The script expects both ber.bsml.list and ber.raw.list to reside in -D.
#########################################################################

use strict;
use warnings FATAL => 'all';

use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;
use IO::File;

my $ber_dir	= undef;
my $db		= undef;
my $log		= undef;

umask(0000);

parse_options();
move_results();

sub print_usage
{
	my $progname = basename($0);
	die << "END";
usage: $progname --ber_directory|-D <workflow_ber_results_directory>
	--database|-d <database> [--log|-l <log>] [--help|-h]
END
}

sub parse_options
{
	my %opts = ();
	GetOptions(\%opts, "ber_directory|D=s", "database|d=s",
			"log|l=s", "help|h");
	print_usage() if $opts{help};
	$ber_dir = $opts{ber_directory} if $opts{ber_directory};
	$db = $opts{database} if $opts{database};
	$log = new IO::File($opts{log}, "w") or
		die "Error writing log $opts{log}: $!"
		if $opts{log};
	print "No BER output directory provided\n" and print_usage()
		if not defined $ber_dir;
	print "No database provided\n" and print_usage()
		if not defined $db;
}

sub move_results
{
	my $btab_list = new IO::File("$ber_dir/ber.btab.list") or
		die "Error reading btab list $ber_dir/ber.btab.list: $!";
	while (my $file = <$btab_list>) {
		chomp $file;
		fix_and_move_btab($file);
	}
	my $raw_list = new IO::File("$ber_dir/ber.raw.list") or
		die "Error reading raw list $ber_dir/ber.raw.list: $!";
	while (my $file = <$raw_list>) {
		chomp $file;
		fix_and_move_raw($file);
	}
}

sub get_model_info
{
	my ($file) = @_;
	$file =~ /[\w\d]+\.model\.(\d+)_(\d+)/;
	return ($1, $2, $&);
}

sub fix_and_move_btab
{
	my ($file) = @_;
	my ($asmbl_id, $model_id, $orig_id) = get_model_info($file);
	my $new_id = "$asmbl_id.m$model_id";
	my $output_dir = "/usr/local/annotation/$db/asmbls/$asmbl_id/" .
		"BER_searches/CURRENT";
	my $offset_diff = length($orig_id) - length($new_id);
	if (! -d $output_dir) {
		system("mkdir -m 0777 -p $output_dir");
		die "Error creating directory $output_dir: $!"
			if $? >> 8;
	}
	my $output_file = "$output_dir/$new_id.nr.btab";
	print $log "Moving $file to $output_file\n" if $log;
	my $btab_in = new IO::File($file) or
		die "Error reading input BTAB $file: $!";
	my $btab_out = new IO::File($output_file, "w") or
		die "Error writing output BTAB $output_file: $!";
	my $line_num = 1;
	while (my $line = <$btab_in>) {
		chomp $line;
		my @tokens = split /\t/, $line;
		$tokens[0] = $new_id;
		$tokens[13] -= $offset_diff * $line_num;
		$tokens[14] -= $offset_diff * $line_num;
		print $btab_out join("\t", @tokens), "\n";
		++$line_num;
	}
}

sub fix_and_move_raw
{
	my ($file) = @_;
	my ($asmbl_id, $model_id, $orig_id) = get_model_info($file);
	my $new_id = "$asmbl_id.m$model_id";
	my $output_dir = "/usr/local/annotation/$db/asmbls/$asmbl_id/" .
		"BER_searches/CURRENT";
	my $offset_diff = length($orig_id) - length($new_id);
	if (! -d $output_dir) {
		system("mkdir -m 0777 -p $output_dir");
		die "Error creating directory $output_dir: $!"
			if $? >> 8;
	}
	my $output_file = "$output_dir/$new_id.nr.gz";
	print $log "Moving $file to $output_file\n" if $log;
	my $raw_in = new IO::File($file) or
		die "Error reading input raw file $file: $!";
	my $raw_out = new IO::File($output_file, ">:gzip") or
		die "Error writing output raw file $output_file: $!";
	while (my $line = <$raw_in>) {
		chomp $line;
		$line =~ s/$orig_id/$new_id/g;
		print $raw_out "$line\n";
	}
}
