#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;
use IO::File;
use Ergatis::IdGenerator;
use GenePredictionUtils::GeneMarkBsmlGenerator;

use constant
{
	EXON_LABELS	=>	{	"Initial"	=>	1,
					"Internal"	=>	1,
					"Terminal"	=>	1,
					"Single"	=>	1,
                    "CDS"       =>  1,
				}
};

my $in			= new IO::File->fdopen(fileno(STDIN), "r");
my $out			= "/dev/stdout";
my $fasta		= undef;
my $project		= "unknown";
my $gene_finder_name	= "GeneMark-ES";
my $sequence_id = undef;
my $id_repository = undef;

&parse_options;

my $id_generator = Ergatis::IdGenerator->new('id_repository' => $id_repository);

&create_bsml;

sub print_usage
{
	my $progname = basename($0);
	die << "END";
usage: $progname [-input|i <input_genscan_data>] [-output|o <output_bsml>]
	[-fasta_file|f <fasta_data>] [-project|p <project_name>] [-sequence_id|s <sequence_id>]
	[-help|h]
END
}

sub parse_options
{
	my %opts = ();
    GetOptions(\%opts, "input|i=s", "output|o=s", "fasta_file|f=s",
		   "project|p=s", "help|h", "sequence_id|s=s", "id_repository|d=s");
	&print_usage if $opts{help};
	$in->open($opts{input}, "r") or
		die "Error reading genscan input $opts{input}: $!"
		if $opts{input};
	$out = $opts{output} if $opts{output};
	$fasta = $opts{fasta_file} if $opts{fasta_file};
	$project = $opts{project} if $opts{project};
    $sequence_id = $opts{sequence_id} if $opts{sequence_id};
    $id_repository = $opts{id_repository} if $opts{id_repository};
}

sub create_bsml
{
    my $bsml_generator = new GenePredictionUtils::GeneMarkBsmlGenerator(
                            id_generator => $id_generator,
                         );
	my %gene_data = ();
	while (my $line = <$in>) {
		chomp $line;
		next if $line =~ /^Sequence name/;
		next if !length($line);
		$line =~ s/^\s+//g;
		my @tokens = split("\t", $line);
		my $gene_num;
        
        if ( $tokens[8] =~ /^gene_id \"(.+?)\"/ ) {
            $gene_num = $1;
        } else {
            die "failed to get a gene_id from this line: ($tokens[8])\n";
        }
        
        if (exists EXON_LABELS->{$tokens[2]}) {
			my $from = $tokens[3] - 1;
			my $to = $tokens[4];
            
            ## handle the direction and adjust for phase
			if ($tokens[6] eq '-') {
				$to -= $tokens[7];
                ($from, $to) = ($to, $from);
                
			} else {
                $to += $tokens[7];
            }
			
            $bsml_generator->AddExon($gene_num, $from, $to);
		}
	}
	$bsml_generator->WriteBsml($out, $sequence_id, $project,
				   $gene_finder_name, $fasta);
}
