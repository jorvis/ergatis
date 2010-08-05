#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

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
use Storable qw(nstore retrieve);
use lib "/usr/local/scratch/bwhitty/lib"; ## temp
#use lib "/export/lib/5.8.8"; ## temp
use Data::Random qw(:all);   ## temp
use Benchmark;
#use DBM::Deep;

use strict;

$|++;

my @results = ();
my $dups = {};
my %options = ();
my $results = GetOptions (  \%options,
#                           'filter_list|f=',
                            'db_list|dl:s',
                            'blast_stored_file|b=s',
                            'output_path|o=s',
                            'output_file|f=s',
                            'write_lists:i',
                            'multiplicity|m:i',
                            'comparisons|c:i',
                            'debug|d=s',
                            'command_id=s',       ## passed by workflow
                            'logconf=s',          ## passed by workflow (not used)
                            'log|l=s',
                            'help|h',
                         ) || pod2usage();


my $comparisons;
my $multiplicity;
my $db_filter = undef;

if($options{'db_list'}) {
    &read_db_list();
}

unless ($options{'comparisons'} || $options{'multiplicity'}) {
}

my $output_path = $options{'output_path'};
$output_path =~ s/\/$//;
if ($output_path eq '') {
    $output_path = '.';
}

my $output_file;
if ($options{'output_file'}) {
    $output_file = $options{'output_file'};
} else {
    $output_file = 'pangenome.profile.txt';
}

print STDERR "Reading stored data...";
if (! -e $options{'blast_stored_file'}) {
    die "no stored blast data file found in output dir";
}

my $temp_ref = retrieve($options{'blast_stored_file'});
## pull the blast results data array out of the stored array
my $results_ref = shift(@{$temp_ref});
@results = @{$results_ref};

## pull the dups hash out of the stored array
my $dups_ref = shift(@{$temp_ref});
$dups = $dups_ref;

print STDERR "done.\n";

my %dbs;
my @genomes;
my $genes={};
#my $genes = DBM::Deep->new( ".pangenome.temp.db" );


print STDERR "Processing results...";

foreach (@results) {
    if(!$db_filter || ($db_filter->{$_->[0]} && $db_filter->{$_->[2]})) {
        if (!defined($genes->{$_->[0]})) {
            $genes->{$_->[0]} = {};
        }
        if (!defined($genes->{$_->[0]}->{$_->[1]})) {
            $genes->{$_->[0]}->{$_->[1]} = {};
        }
        $genes->{$_->[0]}->{$_->[1]}->{$_->[2]} = 1;
        $dbs{$_->[0]} = 1;
    }   
}
@genomes = keys(%dbs);

open(RESULT, ">".$output_path."/".$output_file) || die "couldn't open $output_file for writing";

my %genome_index = ();
my $i_counter = 0;
foreach (@genomes) {
    $genome_index{$_} = $i_counter;
    #print RESULT "## $i_counter\t$_\n";
    $i_counter++;
}

print RESULT "# GENOME\tGENE\t".join("\t",@genomes)."\n";

foreach my $genome (@genomes) {
    foreach my $gene (keys(%{$genes->{$genome}})) {
        print RESULT "$genome\t$gene";
        for (my $i=0; $i < scalar(@genomes); $i++) {
            my $bit;
            if ($genomes[$i] eq $genome) { 
                $bit = '1'; 
            } else {
                $bit = ($genes->{$genome}->{$gene}->{$genomes[$i]} == 1) ? '1' : '0';
            }
            print RESULT "\t$bit";
        }
        print RESULT "\n";
    }
}
exit(0);

sub read_db_list {
    open IN, "<".$options{'db_list'} || die "couldn't open '$options{db_list}' for reading: $!";
    while(<IN>) {
        chomp;
        $db_filter->{$_} = 1;
    }close IN;
    
}
