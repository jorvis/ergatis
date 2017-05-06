#!/usr/bin/env perl

=head1 NAME

perl_template.pl - Description

=head1 SYNOPSIS

 USAGE: perl_template.pl
       --input_file=/path/to/some/input.file
       --output=/path/to/transterm.file
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>

B<--output_file,-o>

B<--log,-l>
    Logfile.

B<--debug,-d>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

 DESCRIPTION
 
=head1  INPUT

    Describe the input

=head1 OUTPUT

    Describe the output

=head1  CONTACT

    Shaun Adkins
    sadkins@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

############# GLOBALS AND CONSTANTS ################

my $BLAST_BIN = "/usr/local/bin/";
my $MCL_BIN = "/usr/local/bin/";
my $MAX_WEIGHT = 316;
my $PROJECT = "lgt";
my $NUM_ALIGNMENTS = 200;

my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
####################################################

my @files;
my %options;

# Allow program to run as module for unit testing if necessary (change extension to .pm)
main() unless caller();

sub main {
    my $results = GetOptions (\%options,
                         "input_file|i=s",
						 "input_list|I=s",
                         "output_dir|o=s",
						 "max_eval|e:1-e5",
						 "min_pct_id|p:0.0",
						 "min_pct_match|m:0.0",
						 "inflation|f:1.5",
						 "num_cpus|n:1",
						 "mixed_genomes|m:1",
						 "blast_bin:$BLAST_BIN",
						 "mcl_bin:$MCL_BIN",
                         "log|l=s",
                         "debug|d=s",
                         "help|h"
                          );

    check_options(\%options);
	write_config_file(\%options, \@files);
}

sub write_config_file {
	my $opts = shift;
	my $files = shift;

	my $output_file = $opts->{'output_dir'} . "/fastortho.config";
	open OUT, ">$output_file" or &_log($ERROR, "Cannot open $output_file for writing: $!");
	
	print OUT "--project_name $PROJECT\n";
	print OUT "--working_directory " . $opts->{'output_dir'} . "\n";
	print OUT "--result_file " . $opts->{'output_dir'} . "/orthologs.txt\n";
	# Options suggest legacy BLAST but BLAST+ paths are default for FastOrtho
	print OUT "--formatdb_path " . $opts->{'blast_bin'} . "/makeblastdb\n";
	print OUT "--blastall_path " . $opts->{'blast_bin'} . "/blastp\n";
	print OUT "--mcl_path " . $opts->{'mcl_bin'} . "/mcl\n";
	print OUT "--pv_cutoff " . $opts->{'max_eval'} . "\n";
	print OUT "--pi_cutoff " . $opts->{'min_pct_id'} . "\n";
	print OUT "--pmatch_cutoff " . $opts->{'min_pct_match'} . "\n";
	print OUT "--maximum_weight $MAX_WEIGHT\n";
	print OUT "--inflation " . $opts->{'inflation'} . "\n";
	print OUT "--blast_cpus " . $opts->{'num_cpus'} . "\n";
	print OUT "--blast_b $NUM_ALIGNMENTS\n";
	print OUT "--blast_e $NUM_ALIGNMENTS\n";
	print OUT "--blast_e " . $opts->{'max_eval'} . "\n";

	# Add all the FASTA files to the config file
	my $genomes_option =  $opts->{'mixed_genomes'} ? "--mixed_genome_fasta" : "--single_genome_fasta";
	print OUT "$genomes_option $_\n" foreach (@$files);

	close OUT;
}

sub check_options {
    my $opts = shift;
    if( $opts->{'help'} ) {
        &_pod;
    }

    if( $opts->{'log'} ) {
        open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
    }

    $debug = $opts->{'debug'} if( $opts->{'debug'} );

    foreach my $req ( qw(output_dir) ) {
        &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
    }

	if($opts->{'input_file'}){
		push @files, $opts->{'input_file'};
	} elsif ($opts->{'input_list'}) {
		open LIST, $opts->{'input_list'} or &_log($ERROR, "Cannot open ". $opts->{'input_list'} . " for reading: $!");
		while (<LIST>) {
			chomp;
			push @files, $_;
		}
		close LIST;
	} else {
		&_log($ERROR, "Need to pass either the --input_file or --input_list option: $!");
	}
}

sub _log {
    my ($level, $msg) = @_;
    if( $level <= $debug ) {
        print STDOUT "$msg\n";
    }
    print $logfh "$msg\n" if( defined( $logfh ) );
    exit(1) if( $level == $ERROR );
}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
