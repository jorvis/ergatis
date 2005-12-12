use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use IO::File;

use Blast::BlastHitDataType;

my %p_values	=	();
my $praze_btab	=	new IO::File->fdopen(fileno(STDIN), "r");
my $out		=	new IO::File->fdopen(fileno(STDOUT), "w");

&parse_options;
&add_pvalues;

sub print_usage
{
	die << "END";
usage: add_pvalue_to_praze_output.pl -blast_btab|b <blast_btab_data>
	-praze_btab|p <praze_brab_data> [-output|o <output>]
	[-help|h]
END
}

sub parse_options
{
	my %opts = ();
	GetOptions(\%opts, "blast_btab|b=s", "praze_btab|p=s",
		   "output|o=s", "help|h");
	print_usage if $opts{help};
	initialize_pvalues($opts{blast_btab}) if $opts{blast_btab};
	$praze_btab->open($opts{praze_btab}) or
		die "Error reading praze btab data $opts{praze_btab}: $!"
		if $opts{praze_btab};
	$out->open($opts{output}, "w") or
		die "Error writing to $opts{output}: $!"
		if $opts{output};
}

sub initialize_pvalues
{
	my $fname = shift;
	my $fh = new IO::File($fname) or
		die "Error reading BLAST btab data $fname: $!";
	while (my $line = <$fh>) {
		chomp $line;
		my $hit = new Blast::BlastHitDataType($line);
		$p_values{$hit->GetSubjectName} = $hit->GetPValue;
	}
}

sub add_pvalues
{
	while (my $line = <$praze_btab>) {
		chomp $line;
		my @tokens = split /\t/, $line;
		my $subject_name = $tokens[5];
		my $p_value = $p_values{$subject_name} or
			die "No p-value found for $subject_name";
		print $out "$line\t$p_value\n";
	}
}
