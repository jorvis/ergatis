#!/usr/local/bin/perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use IO::File;
use File::Basename;

use Blast::BlastHitDataType;

my @input_files		= ();
my %ids			= ();
my $out			= *STDOUT;
my $require_both_ids	= 0;

&parse_options;
&extract_hits;

sub print_usage
{
	my $progname = basename($0);
	die << "END";
usage: $progname [--input|-i <hit_file> |
		  --input_list|-l <hit_file_list>]
		[--output|-o <output>] [--id|-I <id_to_extract>]
		[--require_both_ids|-r]
		[--help|-h]

	i:	hit_file (can be called multiple times)
	I:	id to look for (will extract from either query or subject -
		can be called multiple times)
	r:	require that both query and subject ids be in id list
		[default = false]
END
}

sub parse_options
{
	my %opts = ();
	my $out_fname = undef;
	my @id_list = ();
	my $help = 0;
	my $input_list = undef;
	GetOptions("input|i=s"	=> \@input_files,
		   "input_list|l=s" => \$input_list,
		   "output|o=s" => \$out_fname,
		   "id|I=s" => \@id_list,
		   "require_both_ids|r" => \$require_both_ids,
		   "help|h" => \$help) or &print_usage;
	&print_usage if $help;
	$out = new IO::File($out_fname, "w") or
		die "Error writing to $out_fname: $!"
		if $out_fname;
	foreach my $id (@id_list) {
		my @tokens = split /\s+/, $id;
		#++$ids{$id};
		++@ids{@tokens};
	}
	process_input_list($input_list) if $input_list;
}

sub extract_hits
{
	foreach my $file (@input_files) {
		my $fh = new IO::File($file) or
			die "Error reading hits file $file: $!";
		while (my $line = <$fh>) {
			chomp $line;
			my $hit = new Blast::BlastHitDataType($line);
			my $ok = 1;
			if (scalar(keys %ids)) {
				$ok = $require_both_ids ?
					exists($ids{$hit->GetQueryName()}) &&
					exists($ids{$hit->GetSubjectName()}) :
					exists($ids{$hit->GetQueryName()}) ||
					exists($ids{$hit->GetSubjectName()});
			}
			$out->printf("%s\n", $hit->ToString()) if $ok;
		}
	}
}

sub process_input_list
{
	my $input_list = shift;
	my $fh = new IO::File($input_list) or
		die "Error reading input_list $input_list: $!";
	while (my $file = <$fh>) {
		chomp $file;
		push @input_files, $file;
	}
}
