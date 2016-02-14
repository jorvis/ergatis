#!/usr/bin/perl

=head1  NAME

pangenome_merge_results_pl - Performs the merging of BLAST results to create pangenome data.

=head1 SYNOPSIS

USAGE: pangenome_merge_results.pl
        --input_list=/path/to/pangenome.stored.list
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
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Storable qw(nstore retrieve);
use List::MoreUtils qw(uniq);
use strict;

my $hit_results = {};
my $dups = {};
my %options = ();
my $results = GetOptions (  \%options,
#                           'filter_list|f=',
                            'input_list|i=s',
                            'output_path|o=s',
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

if (!-e $options{'input_list'}) {
    die "must specify an input list which exists";
}

open (IN, $options{'input_list'}) || die "couldn't open input list";

print STDERR "Reading stored data...\n";
my $results_count = 0;
my $file_count = 0;
while (<IN>) {
    chomp;
    my $db;
    my $filename = basename($_);
    
    ## parse the query database name from the filename
    if ($filename =~ /^([^\.]+)\./) {
        $db = $1;
    } else {
        die "failed parsing db name from $_";
    }

    ## unserialize the stored data
    my $temp_ref = retrieve($_) || die "failed unserializing $_";

    ## Add file's query gene-subject genome matches to big hash
	my $added_hits = 0;
	my $hits_before_add = 0;
	my $hits_after_add = 0;
    my $results_ref = shift(@{$temp_ref});
	foreach my $qgenome (keys %$results_ref) {
		foreach my $qgene (keys %{$results_ref->{$qgenome}}) {
			if (defined $hit_results->{$qgenome}->{$qgene}) {
				$hits_before_add += scalar @{$hit_results->{$qgenome}->{$qgene}};
			} else {
				$hits_before_add = 0;
			}
			# Since there are 2 files per query genome (blastp and tblastn) concatenate hits
			push @{$hit_results->{$qgenome}->{$qgene}}, @{$results_ref->{$qgenome}->{$qgene}};
			# Get uniq hits from concatenated gene hits and make that the new gene hits list
			my @uniq_hits = uniq @{$hit_results->{$qgenome}->{$qgene}};
			$hit_results->{$qgenome}->{$qgene} = @uniq_hits;
			$hits_after_add += scalar @uniq_hits;
		}
	}

	# Get updated number of added hits
	$added_hits = $hits_after_add - $hits_before_add;

	# print running total of results;
	$results_count += $added_hits;
	print STDERR "File - " . $file_count++ . "\tAddedHits - " . $added_hits . "\tRunningTotal - " . $results_count . "\n";

    ## pull the dups hash out of the stored array
    my $dups_ref = shift(@{$temp_ref});
    foreach my $org (keys(%{$dups_ref})) {
    	foreach my $dup_set (keys %{$dups_ref->{$org}}) {
        	$dups->{$db}->{$dup_set}=1;
        }
    }

    $temp_ref = undef;
}

nstore([$hit_results, $dups], $options{'output_path'}."/pangenome.blast.stored") || die "couldn't serialize results";

exit(0);

