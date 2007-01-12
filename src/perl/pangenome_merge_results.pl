#!/usr/local/bin/perl

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
use strict;

my @results = ();
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

nstore([\@results, $dups], $options{'output_path'}."/pangenome.blast.stored") || die "couldn't serialize results";

exit(0);

