#!/usr/bin/perl
=head1  NAME

pangenome_make_pangenome.pl - Reads a pangenome profile file and generates
a pangenome data.

=head1 SYNOPSIS

  USAGE: pangenome_do_R.pl
    --profile=/path/to/pangenome.table.txt
    --output_path=/path/to/output/
    --graph_title='Genus species'
    [ --log=/path/to/some/log ]

=head1 OPTIONS

B<--profile,-p>
    The pangenome profile from pangenome_make_profile

B<--comparisons,-c>
    The number of comparisons to make for any 1 value of n (sampling)

B<--multiplicity,-m>
    Another option for sampling based on a multiplicty factor( (sum(m*n) for n=[2...n])=number of comparisons).

B<--output_path,-o>
    Path to which output files will be written.

B<--graph_title,-g>
    Path to which output files will be written.

B<--help,-h>
    This help message/documentation.

=head1   DESCRIPTION

    The pangenome analysis script creates an array of BLAST results data which is then
    processed to create pangenome data.

=head1 INPUT

    The input should be a list of files containing serialized BLAST results array data.

=head1 OUTPUT

    There is no output unless you use the --log option.

=cut

use Pod::Usage;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Data::Random qw(:all);
use List::Util qw(min max);
my %options = ();
my $results = GetOptions( \%options,
                          'profile|p=s',
                          'comparisons|c:s',
                          'multiplicity|m:s',
                          'output_path|o=s',
                          'respect_order|r',
                          'help|h') || pod2usage();



pod2usage if $options{'help'};

# Reads a "profile" file from the pangenome pipeline, and determines distributions of values for
# the total number of genes pertaining to the pan-genome. The total number of permutations 
# evaluated for each number of genomes is maxed at 1000, permutations are randomly selected.

#
# READ first line of profile file, count genomes and assign names to genome columns
#
open IN, "<$options{'profile'}" or die "Unable to open profile file $options{'profile'}\n";
my $l=<IN>;
chop $l;
my @a = split /\t+/,$l;
my $num_genomes = scalar(@a)-2;

my $genome_name;
my %genome_number;
my @genome_names;
for (my $i=2;$i<scalar(@a);$i++){
	$genome_name->{$i-1}=$a[$i];
    push(@genome_names, $a[$i]);
	$genome_number{$a[$i]}=$i-1;
}
#
# READ all profile file, and store a matrix with gene and genome names and 0/1 values
#
my $matrix;
while($l=<IN>){
	chop $l;
	@a = split /\t+/,$l;
	for (my $i=2;$i<scalar(@a);$i++){
		$matrix->{$a[0]}->{$a[1]}->{$genome_name->{$i-1}}=$a[$i];
	}
}
close IN;

my $comparisons = $options{'comparisons'} ? $options{'comparisons'} :(($num_genomes)*1000);
my $multiplicity = $options{'multiplicity'};
if($multiplicity) {
    my ($est_comp, $tot_comp)  =  &estimate_comparisons($num_genomes,$multiplicity);
    $comparisons = $multiplicity ? $est_comp : $comparisons;
}

#
# GENERATE randomly permutations to be used as sequences of genomes for the pan-genome calculation
#
my %seen1;
my %seen2;

# HACK - basically generating a list of indices of the matrix.
# Could probably just use the actual names.
my $all_genomes = [2..($num_genomes+2)];

$options{'respect_order'} = defined($options{'respect_order'}) ? $options{'respect_order'} : 1;
print STDERR "Respecting order\n" if $options{'respect_order'};

open OUT, ">$options{'output_path'}/pangenome.output" or die "Unable to open output file $options{'output_path'}/pangenome.output";
for (my $n=1;$n<=$num_genomes;$n++){

    my $max1;

    # If we are respecting order (i.e. reordering a set counts as a new set) then
    # we'll use a different formula to calculate our max.
    if($options{'respect_order'}) {
        $max1 = int(factorial($num_genomes)/factorial($num_genomes-$n));

    }
    else {
        # Used to be calculating the max using the following but this assumes that order matters:
        $max1 = int(factorial($num_genomes) / (factorial($n) * factorial($num_genomes - ($n)))); 
    }

#	if ($n>2){$max1*=0.8}; # Not sure why we would do this

    # HACK - doing this so that we estimate the number of comparisons per n
	my $max2 = $comparisons/($num_genomes);

	my $iter=0;
#    print STDERR "$iter $max1 $max2\n";

    # Loop until we have the actual max (max1) or the sample max (max2)
	while($iter<$max1 && $iter<$max2){
 		my @genomes;
		my $string ;
		my %seen1;
		my $i=0;

        # Using the rand_set function from Math::Random
        # If we are respecting order than we want to allow for shuffle.
        # If not then we want to prevent shuffles.
        if($options{'respect_order'}) {
            @genomes = rand_set( set => \@genome_names, size => $n, shuffle => 1 );
        }
        else {
            @genomes = rand_set( set => \@genome_names, size => $n, shuffle => 0 );
        }

        $string = join("-", @genomes);
		unless ($seen2{$string}){
			$seen2{$string}=1;
#            my @genomes_names;
#            map {push(@genomes_names,($genome_name->{$_}))}@genomes;
			print OUT $n,"\t",calcpang(@genomes),"\t",join("-",@genomes),"\n";
			$iter++;
            print STDERR "\r".sprintf("%.0f",100*($iter/(min(($max1,$max2)))))."% done with $n";;
		}
        else {
#            print STDERR "Found a duplicate $string\n";
        }
	}
	print STDERR "\nperformed $iter permutations for N=$n genomes\n";
}
close OUT;
print STDERR "Done.\n";

sub calcpang{
#    my $matref = shift;
#    my $matrix = %$matref;
	my @ref_genomes = @_;
	my @done_genomes;
	my $pangenome_size = 0;
	my $pippo=0;
	my $genome1;
	my $genome2;
	my $gene;
	my $count;
    

	foreach $genome1 (@ref_genomes){
		foreach $gene (keys %{ $matrix->{$genome1} }){
			$count = 0;
			foreach $genome2 (@done_genomes){
                $count += $matrix->{$genome1}->{$gene}->{$genome2};
			}
			if ($count==0){$pangenome_size++}
		}
		push @done_genomes, $genome1;
	}
	return $pangenome_size;
}


sub factorial{
        my $n = shift;
        if ($n>1){
                my $facn = 1;
                my $i = 0;
                for ($i=$n;$i>1;$i--){
                        $facn*=$i;
                }
                return $facn;
        }
        else{
                return 1;
        }
}
sub estimate_multiplicity {
    my ($ngenomes, $req_comparisons) = @_;
    my $i = 0;
    my $lower_mult = 0;
    my $lower_comp = 0;
    my $lower_theor = 0;
    my $upper_mult = 0;
    my $upper_comp = 0;
    my $upper_theor = 0;
    my $ldiff = 0;
    my $udiff = 0;
LOOP:   for ($i=5; $i<=5000; $i++){
        my ($a, $b) = estimate_comparisons($ngenomes, $i); 
        if ($a < $req_comparisons){
            $lower_mult = $i;
            $lower_comp = $a;
            $lower_theor = $b;
        } else {
            $upper_mult = $i;
            $upper_comp = $a;
            $upper_theor = $b;
            last LOOP;
        }
    }
    if ($upper_mult == 5){
        my ($a,$b) = estimate_comparisons($ngenomes, 5);
        return (5,$a,$b);
    }
    elsif ($lower_mult == 5000){
        return (5000,$lower_comp,$lower_theor);
    }
    else {
        $ldiff = $req_comparisons - $lower_comp;
        $udiff = $upper_comp - $req_comparisons;
        if ($ldiff <= $udiff){
            return ($lower_mult, $lower_comp, $lower_theor);
        } else {
            return ($upper_mult, $upper_comp, $upper_theor);
        }
    }   
}

sub estimate_comparisons{
    my ($ngenomes, $multiplex) = @_;
    
    my $i = 0;
    my $tot_comparisons = 0;
    my $theor_comparisons = 0;
    my $theor = 0;
    my $real = 0;
    for ($i=2; $i<=$ngenomes; $i++){
        my $theor = int(factorial($ngenomes) / ( factorial($i - 1) * factorial($ngenomes - $i) ));
        my $real = $multiplex * $ngenomes;
        $theor_comparisons += $theor;
        if ($theor < $real){
            $tot_comparisons += $theor;
        } else {
            $tot_comparisons += $real;
        }
    }
    
    return ($tot_comparisons, $theor_comparisons);
}
