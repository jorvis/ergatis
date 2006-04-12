#!/usr/local/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME

do_pangenome_analysis.pl - Performs the merging and analysis of BLAST results to create pangenome data.

=head1 SYNOPSIS

USAGE: pangenome_query_list.pl
        --input_list=/path/to/somefile.dat.list
		--output=/path/to/output.list
      [ --log=/path/to/some/log ]

=head1 OPTIONS

B<--input_list,-i>
    List of files containing serialized array of BLAST data.

B<--output_path,-o>
    Path to which output files will be written.

B<--output_file,-o>
    Full path and name of pangenome data output file.

B<--log,-d>
    optional. Will create a log file with summaries of all actions performed.

B<--debug>
    optional. Will display verbose status messages to STDERR if logging is disabled.

B<--help,-h>
    This help message/documentation.

=head1   DESCRIPTION

	The pangenome analysis script creates an array of BLAST results data which is then
	processed to create pangenome data.

=head1 INPUT

    The input should be a list of files containing serialized BLAST results array data.

=head1 OUTPUT

    There is no output unless you use the --log option.

=head1 CONTACT

    Brett Whitty
    bwhitty@tigr.org

=cut

use Pod::Usage;
use File::Basename;
use XML::Twig;
use Math::Combinatorics;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Storable;
use strict;

my @results = ();
my $dups = {};
my %options = ();
my $results = GetOptions (	\%options,
#							'filter_list|f=',
                          	'input_list|i=s',
              				'output_path|o=s',
							'output_file|f=s',
							'write_lists',
              				'debug|d=s',
              				'command_id=s',       ## passed by workflow
              				'logconf=s',          ## passed by workflow (not used)
              				'log|l=s',
              				'help|h',
		  				 ) || pod2usage();


my $output_path = $options{'output_path'};
$output_path =~ s/\/$//;
if ($output_path eq '') {
	$output_path = '.';
}

my $output_file;
if ($options{'output_file'}) {
	$output_file = $options{'output_file'};
} else {
	$output_file = 'pangenome.table.txt';
}

if (!-e $options{'input_list'}) {
	die "must specify an input list which exists";
}

open (IN, $options{'input_list'}) || die "couldn't open input list";

#if (

while (<IN>) {
	chomp;
	my $db;
	my $filename = basename($_);
	
	## parse the query database name from the filename
	if ($filename =~ /^([^_]+)_/) {
		$db = $1;
	} else {
		die "failed parsing db name from $_";
	}

	## unserialize the stored data
	my $temp_ref = retrieve($_) || die "failed unserializing $_";

	## pull the blast results data array out of the stored array
	my $results_ref = shift(@{$temp_ref});
	push(@results, @{$results_ref});

	## pull the dups hash out of the stored array
	my $dups_ref = shift(@{$temp_ref});
	foreach my $key(keys(%{$dups_ref})) {
		$dups->{$db}->{$key}=1;
	}

	$temp_ref = undef;
}

my %dbs;
my @genomes;
my $genes={};

foreach (@results) {
	if (!defined($genes->{$_->[0]})) {
		$genes->{$_->[0]} = {};
	}
	if (!defined($genes->{$_->[0]}->{$_->[1]})) {
		$genes->{$_->[0]}->{$_->[1]} = {};
	}
	$genes->{$_->[0]}->{$_->[1]}->{$_->[2]} = 1;
	$dbs{$_->[0]} = 1;
}	
@genomes = keys(%dbs);

my $max = scalar @genomes;

open(RESULT, ">".$output_path."/".$output_file) || die "couldn't open $output_file for writing";

foreach (@genomes) {
	my $gene_count = scalar keys(%{$genes->{$_}});
	print RESULT "1\t$gene_count\t$gene_count\n";
}

my $comp_counter = {};

for (my $i = 1; $i <= $max; $i++) {
	my $combinat = Math::Combinatorics->new(
											count => $i,
    	                                    data  => [@genomes],
        	                               );

	while(my @reference = $combinat->next_combination){
		my $ref_string = "(".join(",",@reference).")\n";
		my @comparison_set = @{array_diff(\@genomes, \@reference)};
		foreach my $comp_genome(@comparison_set) {
			print STDERR $ref_string."plus\n$comp_genome\n";
			my $new = 0;
			my $corecount = 0;
			$comp_counter->{$comp_genome}->{$i}++;
			my $genes_by_category = {};
			GENE: foreach my $gene(keys(%{$genes->{$comp_genome}})) {
				my $count = 0;

				foreach my $ref_genome(@reference) {
					if ($genes->{$comp_genome}->{$gene}->{$ref_genome}) {
						$count++;
						$genes_by_category->{$comp_genome}->{'shared'}->{$gene}=1;
					}
				}		
				if ($count == scalar @reference) {
					$corecount++;
					$genes_by_category->{$comp_genome}->{'core'}->{$gene}=1;
				}
				if ($count == 0) {
					$new++;
					$genes_by_category->{$comp_genome}->{'new'}->{$gene}=1;
				}
								
				#A## This block requires that a gene have one
				#A# or more hits to be considered core.
				#A# The changes below require that to be core it must
				#A# be in all
				#A#$count++;
				#A#foreach my $ref_genome(@reference) {
				#A#	if ($genes->{$comp_genome}->{$gene}->{$ref_genome}) {
				#A#		next GENE;
				#A#	}
				#A#}
				#A#$new++;
#				foreach my $ref_genome(@reference) {
#					if ($genes->{$comp_genome}->{$gene}->{$ref_genome}) {
#						$count++;
#					}
#				}
#				if ($count == scalar @reference) {
#					$corecount++;
#				}
#				if ($count == 0) {
#					$new++;
#				}
			}
#			#A#my $corecount = $count-$new;
			print STDERR "Core: $corecount New: $new\n";
			my $rgcount = (scalar @reference) + 1;
			print RESULT "$rgcount\t$corecount\t$new\n";

			if ($options{'write_lists'}) { 
			
				open (OUT, ">".$output_path."/".$comp_genome."_core_".$rgcount."_".$comp_counter->{$comp_genome}->{$i});
				print OUT "#".$comp_genome."\n";
				print OUT "#".join(" ", @reference)."\n";
				foreach my $g(keys(%{$genes_by_category->{$comp_genome}->{'core'}})) {
					print OUT $g."\n";
				}
				close OUT;
				open (OUT, ">".$output_path."/".$comp_genome."_shared_".$rgcount."_".$comp_counter->{$comp_genome}->{$i});
				print OUT "#".$comp_genome."\n";
				print OUT "#".join(" ", @reference)."\n";
				foreach my $g(keys(%{$genes_by_category->{$comp_genome}->{'shared'}})) {
					print OUT $g."\n";
				}
				close OUT;
				open (OUT, ">".$output_path."/".$comp_genome."_new_".$rgcount."_".$comp_counter->{$comp_genome}->{$i});
				print OUT "#".$comp_genome."\n";
				print OUT "#".join(" ", @reference)."\n";
				foreach my $g(keys(%{$genes_by_category->{$comp_genome}->{'new'}})) {
					print OUT $g."\n";
				}
				close OUT;
			}
				
		}
	}
}

exit(0);

sub array_diff {
	my ($superset_ref, $subset_ref) = @_;
	my %index;
	my @diff_array;
	foreach (@{$superset_ref}) {
		$index{$_} = 1;
	}
	foreach (@{$subset_ref}) {
		delete $index{$_};
	}
	@diff_array = keys(%index);
	return \@diff_array;
}
