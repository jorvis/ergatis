#!/usr/bin/perl

=head1  NAME

pangenome_make_table.pl - Creates counts of new and core genes per pangenome size

=head1 SYNOPSIS

USAGE: pangenome_make_table.pl
        --input_list=/path/to/somefile.dat.list
        --output=/path/to/output.list
      [ --log=/path/to/some/log ]

=head1 OPTIONS

B<--blast_stored_file,-b>
    List of files containing serialized array of BLAST data.

B<--output_path,-o>
    Path to which output files will be written.

B<--output_file,-o>
    Full path and name of pangenome data output file.

B<--multiplicity,-m>
    The multiplicity value to be used in data sampling.

B<--comparisons,-c>
    The number of comparisons to sample.

B<--num_threads,-t>
    Optional.  Enable multithreading, using the specified number of threads.

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

	Shaun Adkins
	sadkins@som.umaryland.edu

=cut

use Pod::Usage;
use File::Basename;
use XML::Twig;
use Math::Combinatorics;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Storable qw(nstore retrieve);
use Data::Random qw(:all);    ## temp
use Benchmark;
use bignum;
#use DBM::Deep;

use threads;    # SAdkins - 7/22/16 - Discouraged by Perl but let's try it!
use Thread::Queue;
use Thread::Semaphore;

use warnings;
use strict;

my $dups    = ();
my %options = ();
my @threads;

### Sampling analysis hashes
my $comp_counter = ();
my $dup_counts = ();
my $genes_by_category = ();

my $results = GetOptions(
    \%options,

    #                           'filter_list|f=',
    'db_list|dl:s',
    'blast_stored_file|b=s',
    'output_path|o=s',
    'output_file|f=s',
    'write_lists:i',
    'multiplicity|m:s',
    'comparisons|c:s',
    'num_threads|t:s',
    'debug|d=s',
    'log|l=s',
    'help|h',
) || pod2usage();

my $comparisons;
my $multiplicity;
my $db_filter = undef;

my $max_threads = $options{'num_threads'} ? $options{'num_threads'} : 1;
print STDERR "Max number of threads used: $max_threads\n";
#my @threads = init_threads($max_threads);

if ( $options{'db_list'} ) {
    &read_db_list();
}

unless ( $options{'comparisons'} || $options{'multiplicity'} ) {

    #    die "must provide flag --comparisons or --multiplicity";
    print STDERR "Will run analysis without sampling\n";
}

my $output_path = $options{'output_path'};
$output_path =~ s/\/$//;
if ( $output_path eq '' ) {
    $output_path = '.';
}

my $output_file;
if ( $options{'output_file'} ) {
    $output_file = $options{'output_file'};
} else {
    $output_file = 'pangenome.table.txt';
}

if ( !-e $options{'blast_stored_file'} ) {
    die "must specify an input list which exists";
}

##open (IN, $options{'input_list'}) || die "couldn't open input list";

print STDERR "Reading stored data...";
if ( !-e $options{'blast_stored_file'} ) {
    die "no stored blast data file found in output dir";
}

my $temp_ref = retrieve( $options{'blast_stored_file'} );
## pull the blast results data array out of the stored array
my $results_ref = shift( @{$temp_ref} );

## pull the dups hash out of the stored array
my $dups_ref = shift( @{$temp_ref} );
$dups = $dups_ref;

print STDERR "done.\n";

my %dbs;        #Keeps track of organism name entries
my @genomes;    #Also keeps track of org entries
my $genes = (); # keeps track of query_db -> query_prot -> subject_db

#my $genes = DBM::Deep->new( ".pangenome.temp.db" );

print STDERR "Processing results...";
# Create profile hash structure from the binary input structure
foreach my $qgenome (keys %$results_ref) {
    next if ($db_filter && !$db_filter->{$qgenome});    # skip genome if not in optional filter list
    $genes->{$qgenome} = ();
    $dbs{$qgenome} = 1;
    foreach my $qgene (keys %{$results_ref->{$qgenome}}) {
        $genes->{$qgenome}->{$qgene} = {};
        foreach my $sgenome (@{$results_ref->{$qgenome}->{$qgene}}) {
            next if ($db_filter && !$db_filter->{$sgenome});
            $genes->{$qgenome}->{$qgene}->{$sgenome} = 1;
        }
    }
}
@genomes = keys(%dbs);

open( RESULT, ">" . $output_path . "/" . $output_file )
  || die "couldn't open $output_file for writing";

my %genome_index = ();
my $i_counter    = 1;
foreach (@genomes) {
    $genome_index{$_} = $i_counter;

    #print RESULT "## $i_counter\t$_\n";
    $i_counter++;
}

# output the gene counts for single genomes
foreach (@genomes) {
    my $genome_dup_count = 0;
    my $gene_count       = scalar keys( %{ $genes->{$_} } );
    foreach my $dup_set ( keys( %{ $dups->{$_} } ) ) {
        my @dup_genes = split( " ", $dup_set );
        $genome_dup_count += scalar(@dup_genes);
    }
    ## our output table will take the form:
    ## genomecount\tcore\tshared\tnew\tcore dup count\tshared dup\tnew dup
    print RESULT "1\t0\t0\t$gene_count\t0\t0\t$genome_dup_count\t"
      . $genome_index{$_}
      . "\t$_\n";
}
print STDERR "done.\n";
if ( $options{'comparisons'} || $options{'multiplicity'} ) {
    &do_analysis_with_sampling();
} else {
    # Will take a lot longer to finish due to doing all possible combinations
    &do_analysis_no_sampling();
}

exit(0);

sub do_analysis_with_sampling {
    if ( $options{'comparisons'} ) {
        my ( $est_multiplicity, $t, $t2 ) =
          estimate_multiplicity( scalar(@genomes), $options{'comparisons'} );
        print STDERR "# of comparisons: $options{comparisons}\nest. multiplicity: $est_multiplicity\n";
        $options{'multiplicity'} = $est_multiplicity;
    } else {
        my ( $est_comparisons, $t ) =
          estimate_comparisons( scalar(@genomes), $options{'multiplicity'} );
        print STDERR "multiplicity: $options{multiplicity}\nest. # comparisons: $est_comparisons\n";
    }

    $comp_counter = ();

    my $max = scalar @genomes;

    print STDERR "Total of $max genomes in results set\n";
    print STDERR "Running analysis...\n";

    my %seen;

    my @i_genomes = ( 0 .. ( $max - 1 ) );

    # Iterate using combinations of reference genomes of size $i
    for ( my $i = 1; $i < $max; $i++ ) {

        # start timer
        my $start = new Benchmark;

        print STDERR "Running analysis with $i reference genomes...";

        ## new
        my $true_max =
          int( factorial( $max - 1 ) /
              ( factorial( ( $i + 1 ) - 1 ) * factorial( $max - ( $i + 1 ) ) )
          );
#        print STDERR "Maximum number of computations for $i on $max genomes: $true_max\n";

        # Iterate through each available genome
        for ( my $j = 0; $j < $max; $j++ ) {
            #Get comparison genome name and index number
            my $comp_genome = $genomes[$j];

            my @ref_genomes = @i_genomes;
            splice @ref_genomes, $j, 1;

            # Allocated concurrent thread limit
#           my $sem = Thread::Semaphore->new($max_threads);

            my $point_count = 0;

#           @threads = map { threads->create( sub {
#               $sem->down();   # Note resource is being used
                # Iterate until either the number of multiplicity rounds is reached, or all possible combinations of ref genomes is seen
                while ( $point_count < min($options{'multiplicity'}, $true_max) ) {
                    # Choose a random set of N reference genomes
                    my @reference_set =
                      rand_set( set => \@ref_genomes, size => $i, shuffle => 0 );

					# Create a vector of genome ids among the reference set
                    my @seen_vec = ();
                    for ( my $ii = 0; $ii < scalar(@reference_set); $ii++ ) {
                        push @seen_vec, $reference_set[$ii];
                    }
					# Sort the vector, numerically, and create a string
					my @sorted_vec = sort {$a <=> $b} @seen_vec;
					my $sort_vec_string = join(':', @sorted_vec);

                    unless ( defined $seen{$j}{$sort_vec_string} ) {
                        $dup_counts = ();
                        # Counter to track number of iterations of comparison genome in use per size-$i ref genomes
                        $comp_counter->{$comp_genome}->{$i}++;
                        $genes_by_category = ();

                        # Iterate through the genes in the query genome
                        GENE:
                        foreach my $gene ( keys( %{ $genes->{$comp_genome} } ) ) {
                            # Iterate through our N reference genome indexes
                            foreach my $i_ref_genome (@reference_set) {
                                # grabbing the subject db name
                                my $ref_genome = $genomes[$i_ref_genome];
                                my $hit = $genes->{$comp_genome}->{$gene}->{$ref_genome};
                                categorize_shared_gene($hit, $comp_genome, $gene);
                            }

                            my $count = (defined $genes_by_category->{$comp_genome}->{'shared'}->{$gene}) ? $genes_by_category->{$comp_genome}->{'shared'}->{$gene} : 0;
                            categorize_core_gene($comp_genome, $gene) if $count == scalar @reference_set;
                            categorize_new_gene($comp_genome, $gene) if $count == 0;
                          } #/GENE

                          count_dups($comp_genome);
                          write_results($comp_genome, \@reference_set, $i);

                          ## Ensure we don't use this set of genomes again
                          $seen{$j}{$sort_vec_string} = 1;
						  # Only increment point count if we had a unique reference genome set
						  $point_count++;
                    } ## end of unless
                } ## end of while 'point_count'
#                $sem->up(); # Free up resource
#            }); ## end of thread->create sub
#            } ## end of threads
#            $_->join() for @threads;
        }
        # end timer
        my $end = new Benchmark;

        # calculate difference
        my $diff = timediff( $end, $start );

        print STDERR "runtime: " . timestr( $diff, 'noc' ) . "\n";
    }
}

sub do_analysis_no_sampling {
    $comp_counter = ();
    my $max          = scalar @genomes;

    # Get combinations of reference genomes of size $i
    for ( my $i = 1; $i < $max; $i++ ) {
        # start timer
        my $start = new Benchmark;

        my $combinat = Math::Combinatorics->new(
            count => $i,
            data  => [@genomes],
        );

        # For each combination of reference genomes...
        while ( my @reference = $combinat->next_combination ) {
            my $ref_string = "(" . join( ",", @reference ) . ")\n";

            # Allocated concurrent thread limit
            my $sem = Thread::Semaphore->new($max_threads);            # Grab the genomes not in our reference set.
            my @comparison_set = @{ array_diff( \@genomes, \@reference ) };

            @threads = map { threads->create( sub {
                $sem->down();   # Note resource is being used
                # Each comparison genome will be compared to the reference genomes
                foreach my $comp_genome (@comparison_set) {
                    $dup_counts = ();
                    $comp_counter->{$comp_genome}->{$i}++;
                    $genes_by_category = ();

                    # Process each gene from the comparison genome
                    GENE:
                    foreach my $gene ( keys( %{ $genes->{$comp_genome} } ) ) {
                        # Keep track of how many ref genomes had a hit for this particular gene
                        foreach my $ref_genome (@reference) {
                            my $hit = $genes->{$comp_genome}->{$gene}->{$ref_genome};
                            categorize_shared_gene($hit, $comp_genome, $gene);
                        }
                        my $count = (defined $genes_by_category->{$comp_genome}->{'shared'}->{$gene}) ? $genes_by_category->{$comp_genome}->{'shared'}->{$gene} : 0;
                        categorize_core_gene($comp_genome, $gene) if $count == scalar @reference;
                        categorize_new_gene($comp_genome, $gene) if $count == 0;

                    } # /GENE

                   count_dups($comp_genome);
                   write_results($comp_genome, \@reference, $i);
                }
                $sem->up(); # Free up resource
            });
            }
            $_->join() for @threads;
        }
        # end timer
        my $end = new Benchmark;

        # calculate difference
        my $diff = timediff( $end, $start );

        print "runtime: " . timestr( $diff, 'noc' ) . "\n";
    }
}

# Return the minimum of the two values
sub min {
    my ($a, $b) = @_;
    return $a < $b ? $a : $b;
}

# Return an array of elements present in array A but not in array B
sub array_diff {
    my ( $superset_ref, $subset_ref ) = @_;
    my %index;
    my @diff_array;
    foreach ( @{$superset_ref} ) {
        $index{$_} = 1;
    }
    foreach ( @{$subset_ref} ) {
        delete $index{$_};
    }
    @diff_array = keys(%index);
    return \@diff_array;
}

no warnings;    # Guards against 'deep recursion' warnings
sub factorial {
    my $i = shift @_;
    if ( $i <= 1 ) {
        return 1;
    } else {
        return ( $i * factorial( $i - 1 ) );
    }
}
use warnings;

sub estimate_multiplicity {
    my ( $ngenomes, $req_comparisons ) = @_;
    my $i           = 0;
    my $lower_mult  = 0;
    my $lower_comp  = 0;
    my $lower_theor = 0;
    my $upper_mult  = 0;
    my $upper_comp  = 0;
    my $upper_theor = 0;
    my $ldiff       = 0;
    my $udiff       = 0;
  LOOP: for ( $i = 5; $i <= 5000; $i++ ) {
        my ( $a, $b ) = estimate_comparisons( $ngenomes, $i );
        if ( $a < $req_comparisons ) {
            $lower_mult  = $i;
            $lower_comp  = $a;
            $lower_theor = $b;
        } else {
            $upper_mult  = $i;
            $upper_comp  = $a;
            $upper_theor = $b;
            last LOOP;
        }
    }
    if ( $upper_mult == 5 ) {
        my ( $a, $b ) = estimate_comparisons( $ngenomes, 5 );
        return ( 5, $a, $b );
    } elsif ( $lower_mult == 5000 ) {
        return ( 5000, $lower_comp, $lower_theor );
    } else {
        $ldiff = $req_comparisons - $lower_comp;
        $udiff = $upper_comp - $req_comparisons;
        if ( $ldiff <= $udiff ) {
            return ( $lower_mult, $lower_comp, $lower_theor );
        } else {
            return ( $upper_mult, $upper_comp, $upper_theor );
        }
    }
}

sub estimate_comparisons {
    my ( $ngenomes, $multiplex ) = @_;

    my $i                 = 0;
    my $tot_comparisons   = 0;
    my $theor_comparisons = 0;
    my $theor             = 0;
    my $real              = 0;
    for ( $i = 2; $i <= $ngenomes; $i++ ) {
        my $theor = factorial($ngenomes) /
          ( factorial( $i - 1 ) * factorial( $ngenomes - $i ) );
        my $real = $multiplex * $ngenomes;
        $theor_comparisons += $theor;
        $tot_comparisons += min($theor, $real);
    }

    return ( $tot_comparisons, $theor_comparisons );
}

sub read_db_list {
    open IN, "<" . $options{'db_list'}
      || die "couldn't open '$options{db_list}' for reading: $!";
    while (<IN>) {
        chomp;
        $db_filter->{$_} = 1;
    }
    close IN;

}

# Subroutine to force a specific size for the max number of threads
sub init_threads {
    # An array to place our threads in
    my $num_of_threads = shift;
	my @initThreads;
	for(my $i = 1;$i<=$num_of_threads;$i++){
		push(@initThreads,$i);
	}
	return @initThreads;
}

# Determine if the given gene is also present in another genome up to this point.
sub categorize_shared_gene {
    my ($hit, $comp_genome, $gene) = @_;
    ## increment shared count for this gene if we have a hit
    $genes_by_category->{$comp_genome}->{'shared'}->{$gene}++ if $hit;
}

# Determine if the given gene is present in all other genomes up to this point.
sub categorize_core_gene {
    my ($comp_genome, $gene) = @_;
    $genes_by_category->{$comp_genome}->{'core'}->{$gene} = 1;
}

# Determine if the given gene is present in no other genomes up to this point.
sub categorize_new_gene {
    my ($comp_genome, $gene) = @_;
    $genes_by_category->{$comp_genome}->{'new'}->{$gene} = 1;
}

# Count numbers of core, shared, and new genes amongst duplicate hits
sub count_dups {
    my $comp_genome = shift;
    ## process lists to see how many duplicated genes are in each category
    foreach my $cat ( ( 'core', 'shared', 'new' ) ) {
        $dup_counts->{$comp_genome}->{$cat} = 0;
        ## we'll look at each set of duplicates one at a time
        foreach
          my $dup_set ( keys( %{ $dups->{$comp_genome} } ) )
        {
            my @dup_genes   = split( " ", $dup_set );
            my $dup_count   = scalar(@dup_genes);
            my $dup_counter = 0;

            ## for each gene of a duplicate set
            foreach my $dup (@dup_genes) {
                ## check if it's in the category
                if ( $genes_by_category->{$comp_genome}->{$cat}
                    ->{$dup} )
                {
                    $dup_counter++;
                }
            }
            ## check if all of the dups of a set weren't found in the same category
            if (   $dup_counter > 0
                && $dup_count != $dup_counter )
            {
                print STDERR
                  "Only $dup_counter of the following dup set found:\n$dup_set\n***This could be a problem: Comp_genome $comp_genome\n";
                ## and if they are, then add the dup overcount to the total
            } elsif ( $dup_counter == $dup_count ) {
                $dup_counts->{$comp_genome}->{$cat} +=
                  $dup_count;
            }
        }
    }
}

# Write the results to file
sub write_results {
    my $comp_genome = shift;
    my $ref = shift;
    my $size = shift;

    my $core_count =
      scalar(
        keys( %{ $genes_by_category->{$comp_genome}->{'core'} } ) );
    my $shared_count =
      scalar(
        keys( %{ $genes_by_category->{$comp_genome}->{'shared'} } ) );
    my $new_count =
      scalar(
        keys( %{ $genes_by_category->{$comp_genome}->{'new'} } ) );

    my $core_dup_count   = $dup_counts->{$comp_genome}->{'core'};
    my $shared_dup_count = $dup_counts->{$comp_genome}->{'shared'};
    my $new_dup_count    = $dup_counts->{$comp_genome}->{'new'};

    # This equals the number of references plus our one comparison genome
    my $rgcount = ( scalar @$ref ) + 1;

    print RESULT
      "$rgcount\t$core_count\t$shared_count\t$new_count\t$core_dup_count\t$shared_dup_count\t$new_dup_count\t"
      . $genome_index{$comp_genome}
      . "\t$comp_genome\n";

# Would love to make this its own subroutine eventually but don't feel like dealing with all the extra variables
    if ( $options{'write_lists'} ) {
        write_list("core", $size, $rgcount, $ref, $comp_genome);
        write_list("new", $size, $rgcount, $ref, $comp_genome);
        write_list("shared", $size, $rgcount, $ref, $comp_genome);
    }
}

sub write_list {
    my ($group, $size, $rgcount, $ref, $comp_genome) = @_;
    open( OUT,
            ">"
          . $output_path . "/"
          . $comp_genome . "_${group}_"
          . $rgcount . "_"
          . $comp_counter->{$comp_genome}->{$size} );
    print OUT "#" . $comp_genome . "\n";
    print OUT "#" . join( " ", @$ref ) . "\n";
    foreach my $g (keys(%{ $genes_by_category->{$comp_genome}->{$group} })){
        print OUT $g . "\n";
    }
    close OUT;
}
