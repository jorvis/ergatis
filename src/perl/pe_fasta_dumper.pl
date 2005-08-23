use strict;
use warnings;

use IO::File;
use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;

use Fasta::FastaIndexedReader;

my $indexer		= undef;
my $clusters		= new IO::File->fdopen(fileno(STDIN), "r");
my $output_dir		= undef;
my $chunk_size		= 50;

&parse_options;
&dump_fasta;

sub print_usage
{
	my $progname = basename($0);
	die << "END";
usage: $progname
	-clusters|c <tab_delimited_cluster_file>
	-output_dir|d <output_dir>
	{-fasta_file|f <fasta_file> ...}
	[-fasta_list|F <fasta_list>]
	[-help|h]
END
}

sub parse_options
{
	my $cluster_data	= undef;
	my $help	= undef;
	my @fasta_files	= ();
	my $fasta_list	= undef;
	GetOptions("clusters|c=s"	=> \$cluster_data,
		   "output_dir|d=s"	=> \$output_dir,
		   "help|h"		=> \$help,
		   "fasta_file|f=s"	=> \@fasta_files,
		   "fasta_list|F=s"	=> \$fasta_list);
	&print_usage if $help;
	$clusters->open($cluster_data) or
		die "Error reading cluster data $cluster_data: $!"
		if $cluster_data;
	&process_fasta_list(\@fasta_files, $fasta_list) if $fasta_list;
	print STDERR "No output directory provided\n" and &print_usage
		if !$output_dir;
	print STDERR "No FASTA file(s) provided\n" and &print_usage
		if !scalar(@fasta_files);
	$indexer = new Fasta::FastaIndexedReader(\@fasta_files, 1);
}

sub dump_fasta
{
	my $cid = 0;
	my $chunk = 0;
	my $dir = "$output_dir";
	while (my $line = <$clusters>) {
		if ($cid % $chunk_size == 0) {
			$dir = "$output_dir/$chunk";
			mkdir($dir);
			++$chunk;
		}
		chomp $line;
		my @ids = split /\t/, $line;
		my $fname = "$dir/cluster_$cid.fsa";
		my $out = new IO::File($fname, "w")
			or die "Error writing FASTA data to $fname: $!";
		++$cid;
		foreach my $id (@ids) {
			my $seq = $indexer->FetchSequence($id) or
				die "Error fetching sequence for $id\n";
			$seq->Print($out);
		}
		print $out "\n";
	}
}

sub process_fasta_list
{
	my ($fasta_files, $fasta_list) = @_;
	my $data = new IO::File($fasta_list) or
		die "Error reading FASTA list $fasta_list: $!";
	while (my $line = <$data>) {
		chomp $line;
		push @$fasta_files, $line if length $line;
	}
}
