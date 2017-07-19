#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
=head1  NAME

run_sighunt.pl - Run the SigHunt R package using the assembly input

=head1 SYNOPSIS

  USAGE: run_sighunt.pl
    --input_file=/path/to/assembly.fa
    --r_script=/path/to/run_sighunt.R
    --output_path=/path/to/output/
    [ --log=/path/to/some/log ]

=head1 OPTIONS

B<--input_file,-i>
	Path to an assembly fasta file

B<--cutoff, -c>
    DIAS value to exceed when filtering interval regions

B<--r_script,-r>
    The R script to run on the input pangenome_table

B<--output_path,-o>
    Path to which output files will be written.

B<--help,-h>
    This help message/documentation.

=head1   DESCRIPTION

    This script serves as a wrapper to run the R package SigHunt, which is used to detect
	lateral gene transfer (LGT) in eukaryotic organisms

=head1 INPUT

    The input should be an assembly fasta file.

=head1 OUTPUT

    1) A plot of the assemblies' Discrete Interval Accumulated Score (DIAS)
		- DIAS measures how many 4-mers deviate in their frequency
		from the local background of the genomic sequence and by how
		much
	2) List of positions for candidate regions (those with DIAS greater than chosen cutoff)

=cut

use Pod::Usage;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my $CUTOFF = 8;

my %options = ();
my $results = GetOptions( \%options,
                          'input_file|i=s',
						  'cutoff|c=i',
                          'r_script|r=s',
			              'r_exec_path|p=s',
                          'output_path|o=s',
                          'help|h') || pod2usage();

pod2usage if $options{'help'};


my $input_file = $options{'input_file'};
my $r_script = $options{'r_script'};
my $output_path = $options{'output_path'};
my $R_EXEC_PATH = $options{'r_exec_path'};
my $cutoff = defined $options{'cutoff'} ? $options{'cutoff'} : $CUTOFF;

# Just keep basename for R script template
my $r_filename = $r_script;
$r_filename =~  s/^.*\/([^\/]*)/$1/;

open (IN, "$r_script") || die "Couldn't open R script '$r_script': $!"; 

# Get the R script name for the output
my $input_r = "$output_path/$r_filename"."in";

# In the R-script template, customize for current job by subbing certain fields
open (OUT, ">$input_r");
while (<IN>) {
	s/###cutoff###/$cutoff/;
    s/###input_file###/$input_file/;
    s/###output_path###/$output_path/;
    print OUT;
}
close OUT;
close IN;

# Using 'CMD BATCH' in arguments allows for R to be executed in a wrapper script
# Adding --quiet just silences the R startup prompt, but not output.
if(system("$R_EXEC_PATH CMD BATCH --quiet $input_r $output_path/".$r_filename."out") !=0 ) {
    warn "Unable to run the R command $!\n";
}

print STDERR "done.\n";
exit();
