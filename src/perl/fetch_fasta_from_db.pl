#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

use warnings;
use strict;
use IO::File;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;

my $input_file	= 0;
my $input_list	= 0;
my @ids		= ();
my $out		= new IO::File("/dev/stdout", "w");
my $db		= 0;
my $format	= "xdformat";
my $protein	= 0;
my %BINS	= ( "cdbfasta"	=> "/usr/local/bin/cdbyank -w",
		    "xdformat"	=> "/usr/local/packages/blast2/xdget",
		    "formatdb"	=> "/usr/local/bin/fastacmd"
		  );

sub parse_options;
sub print_usage;
sub dump_fasta;

&parse_options;
&dump_fasta;

sub print_usage
{
	my $progname = basename($0);
	die << "END";
usage: $progname [-i <id> | -I <id_list>] -d <database>
	[-o <output>] -f <format> -p <[t,T]/[f,F]> [-h]

	f:	format [cdbfasta, xdformat, formatdb]
	p:	t/T if database contains protein data
END
}

sub parse_options
{
	my %opts	= ();
	GetOptions(\%opts,
		   "i|input_id=s",
		   "I|input_list=s",
		   "o|output=s",
		   "d|database=s",
		   "f|format=s",
		   "p|protein=s",
		   "h|help") || &print_usage;
	while (my ($key, $val) = each %opts) {
		if ($key eq "i") {
			my @tokens = split /\s+/, $val;
			push @ids, @tokens;
		}
		elsif ($key eq "I") {
			my $fh = new IO::File($val) or
				die "Error reading input $val: $!";
			while (<$fh>) {
				chomp;
				my @tokens = split /\s+/;
				push @ids, @tokens;
			}
		}
		elsif ($key eq "o") {
			$out->open($val, "w") or
				die "Error writing to $val: $!";
		}
		elsif ($key eq "d") {
			$db = $val;
		}
		elsif ($key eq "f") {
			$format = $val;
			print "Invalid format: $format" and print_usage
				unless exists $BINS{$format};
		}
		elsif ($key eq "p") {
			$protein = 1 if uc($val) eq "T";
		}
		elsif ($key eq "h") {
			&print_usage;
		}
	}
}

sub dump_fasta
{
	my $cmd = "$BINS{$format}";
	if ($format eq "cdbfasta") {
		$cmd .= " $db.cidx -a";
	}
	elsif ($format eq "xdformat") {
		$cmd .= ($protein ? " -p" : " -n") . " $db";
	}
	elsif ($format eq "formatdb") {
		$cmd .= " -d $db -p " . ($protein ? "T" : "F") . " -s";
	}
	foreach my $id (@ids) {
		print $out `$cmd "$id"` || die "Error fetching id $id\n";
	}
}
