#!/usr/local/packages/perl-5.8.5/bin/perl

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
							'write_lists:i',
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

print STDERR "Reading stored data...\n";
if (! -e $options{'output_path'}."/pangenome.blast.stored") {
	die "no stored blast data file found in output dir";
}

my $temp_ref = retrieve($options{'output_path'}."/pangenome.blast.stored");
## pull the blast results data array out of the stored array
my $results_ref = shift(@{$temp_ref});
@results = @{$results_ref};

## pull the dups hash out of the stored array
my $dups_ref = shift(@{$temp_ref});
$dups = $dups_ref;

my %dbs;
my @genomes;
my $genes={};

print STDERR "Processing results...\n";

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

# output the gene counts for single genomes
foreach (@genomes) {
	my $genome_dup_count = 0;
	my $gene_count = scalar keys(%{$genes->{$_}});
	foreach my $dup_set(keys(%{$dups->{$_}})) {
			my @dup_genes = split(" ", $dup_set);
			$genome_dup_count += (scalar(@dup_genes) - 1);
	}
	## our output table will take the form:	
	## genomecount\tcore\tshared\tnew\tcore dup count\tshared dup\tnew dup
	print RESULT "1\t0\t0\t$gene_count\t0\t0\t$genome_dup_count\n";
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
			my $dup_counts = {};
			$comp_counter->{$comp_genome}->{$i}++;
			my $genes_by_category = {};
			GENE: foreach my $gene(keys(%{$genes->{$comp_genome}})) {
				my $count = 0;

				foreach my $ref_genome(@reference) {
					## if the hash value == 1
					if ($genes->{$comp_genome}->{$gene}->{$ref_genome}) {
						$count++;
						## we have a hit
						$genes_by_category->{$comp_genome}->{'shared'}->{$gene}=1;
					}
				}		
				## if the number of genomes we have hits against is equal to the total
				if ($count == scalar @reference) {
					## then this is a core gene
					$genes_by_category->{$comp_genome}->{'core'}->{$gene}=1;
				}
				## if we didn't have any hits at all
				if ($count == 0) {
					## then this is a new (strain specific) gene
					$genes_by_category->{$comp_genome}->{'new'}->{$gene}=1;
				}
			}
			my $rgcount = (scalar @reference) + 1;
			
			## process lists to see how many duplicated genes are in each category
			foreach my $cat(('core', 'shared', 'new')) {
				$dup_counts->{$comp_genome}->{$cat} = 0;
				## we'll look at each set of duplicates one at a time
				foreach my $dup_set(keys(%{$dups->{$comp_genome}})) {
					my @dup_genes = split(" ", $dup_set);
					my $dup_count = scalar(@dup_genes);
					my $dup_counter=0;
					## for each gene of a duplicate set
					foreach my $dup(@dup_genes) {
						## check if it's in the category
						if ($genes_by_category->{$comp_genome}->{$cat}->{$dup}) {
							$dup_counter++;
						}
					}
					## check if all of the dups of a set weren't found in the same category
					if ($dup_counter > 0 && $dup_count != $dup_counter) {
						print STDERR "Only $dup_counter of the following dup set found:\n$dup_set\n***This could be a problem\n";
					## and if they are, then add the dup overcount to the total
					} elsif ($dup_counter == $dup_count) {
						$dup_counts->{$comp_genome}->{$cat} += ($dup_count - 1);
					}
				}
			}
			
			my $core_count     = scalar(keys(%{$genes_by_category->{$comp_genome}->{'core'}}));
			my $shared_count   = scalar(keys(%{$genes_by_category->{$comp_genome}->{'shared'}}));
			my $new_count      = scalar(keys(%{$genes_by_category->{$comp_genome}->{'new'}}));
			my $core_dup_count   = $dup_counts->{$comp_genome}->{'core'};
			my $shared_dup_count = $dup_counts->{$comp_genome}->{'shared'};
			my $new_dup_count    = $dup_counts->{$comp_genome}->{'new'};
			
			print RESULT "$rgcount\t$core_count\t$shared_count\t$new_count\t$core_dup_count\t$shared_dup_count\t$new_dup_count\n";

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
