#!/usr/local/bin/perl -w

=head1 NAME

merge_twinscan_raw.pl - pull together all raw twinscan files generated from the same source.

=head1 SYNOPSIS

    USAGE: merge_twinscan_raw.pl
            --input_list|-i=/path/to/twinwise.raw.list
            --chunk_size|-c [ size of chunks.]
            --overlap_size|-o [size of overlaps between chunks]

=head1 OPTIONS

B<--input_list, -i>
    A twinscan.raw.list file, containing paths to all of the files in the output
    of a single twinscan run.

B<--chunk_size, -c>
    Size of the largest chunks of sequence used as input.  The whole reason for having
    to make this silly program is that sometimes the sequences have been split.
    NOTE: Currently not implemented, since the initial coordinate is given in the name
    of the file, in this initial use case.

B<--overlap_size, -o>
    Size of the overlaps between adjacent chunks fo sequence.  These regions will be
    examined closely to remove the smaller of redundant features.
    NOTE:  Currently not implemented, since we're emulating the current script which
    resolves overlaps directly in the db.

=head1  DESCRIPTION

    merge_twinscan_raw.pl is intended to be part of a process that allows results of
    twinscan to be loaded into chado databases when run as previously done outside of
    ergatis. 

    The current method of running twinscan is to have the input be chunked into sections
    of some predetermined size, resulting in possibly several raw output files per
    original input sequence.  The initial coordinate is given in the raw output file name,
    and the coordinates of features are adjusted on the fly during the loading process.
    Once the data have been loaded, an additional script is run to remove the smaller of
    any 'redundant' features.

    This 3 step process is irritating in three ways:
    First, the initial process of chunking input, based on our experiences with aat, gives
    very little improvement with regards to actual time.  The cost of computing twice over
    the large areas that must overlap and the increased complexity of the process quickly
    eat into the time savings gained by the splitting of the sequences.

    Second, the manipulation of data on the fly as it is loaded created a disjointed data-
    set, with values appearing in the database which don't appear on the filesystem, which
    could create problems if the files were needed for use elsewhere.

    Thirdly, becuase the process was grown up around the legacy db schema, it is impossible
    to simply 'load' into chado, even after conversion to bsml.

    Because of all three of the above steps ,a second merge script, working on either the
    filesystem or the chado target db, is needed.  I've chosen to do the former, creating
    a set of files that can be loaded into the legacy db without further maniupulation, is
    representative of the data downstream, can easily be converted into bsml, and thus also
    loaded into chado with no additional manipulation.

=head1  INPUT

    Input is a listing of the output files from a twinscan run produced in the legacy
    fashion.  Note that this script will be rendered obsolete once an ergatis twinscan
    component is produced.

=head1  OUTPUT

    Output is a set of raw files where each file contains the complete and non-redundant
    set of predictions for a single input sequence, and a file listing members of the set.

=head1  CONTACT

    Jason Inman
    jinman@tigr.org

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
$|++;

my $input_list;
my $overlap;
my $output_dir;

my %options = ();
my $results = GetOptions(\%options,
                            'input_list|i=s',
                            'overlap|o=i',
                            'output_dir|d=s',
                        );

&check_paramaters;

# open the input list
open (LIST, "< $input_list") || die "Can't open $input_list: $!\n";

my $cur_base_file = '';
my $cur_output = '';

while (<LIST>) {

    chomp;
    my $raw_file = $_;

    # get the file name and true starting coordinate.
    my $basename;
    my $increment;
    my $output_file;

    if ($raw_file =~ /.*\/([^\/]+)$/) {
        $basename = $1;
        if ($basename =~ /(.*)\.(\d+)\.twinscan/) {
            $increment = $2;
            $output_file = "$output_dir/$1.twinscan";
        }
    }

    unless (defined $basename && defined $increment) {
        die "Unexpected filename found: $_\n";
    }

    # go ahead and put move the old results into the new ones.
    &merge_file($raw_file,$output_file,$increment);

    # If this is a new file, we'll need to possibly resolve overlapping predictions.
    # Even if there was only one file, we need to point back to the original source
    if ($cur_base_file ne $basename) {
        &resolve_overlaps($cur_output,$cur_base_file) unless ($cur_output eq '');
        $cur_base_file = $basename;
        $cur_output = $output_file;
    }

}

exit(0);


sub check_paramaters {

    if ($options{'input_list'}) {
        $input_list = $options{'input_list'};
        die "Can't find input list $input_list\n" unless (-f $input_list);
    } else {
        die "Please specify a list of input using --input_list\n";
    }

    if ($options{'overlap'}) {
        $overlap = $options{'overlap'};
        die "can't use negative overlap value.\n" if ($overlap < 0);
    } else {
#        die "Please specify an amount of overlap using --overlap\n";
    }

    if ($options{'output_dir'}) {
        $output_dir = $options{'output_dir'};
        die "$output_dir is not a directory/doesn't exist!\n" unless (-d $output_dir);
        die "$output_dir is not writeable\n" unless (-w $output_dir );
    } else {
        die "Please specify an output directory using --output_dir\n";
    }

}

sub merge_file {

    my ($raw,$out,$inc) = @_;

    # open the input and output handlers
    open (IN,"< $raw") || die "Can't open $raw: $!\n";
    if ($inc == 0) {
        open (OUT,"> $out") || die "Can't open $out for writing: $!\n";
    } else {
        open (OUT,">> $out") || die "Can't open $out for appending: $!\n";
    }

    # copy the file contents, with slight differences if $increment is zero or non-zero
    while (<IN>) {

        if (/^#/) {
            next if $inc;
        }

        if ($inc > 0) {
            my @fields = split(/\t/,$_);
            $fields[3] += $inc;
            $fields[4] += $inc;
            $_ = join ("\t",@fields);
        }

        print OUT;

    }

    close IN;
    close OUT;

}

sub resolve_overlaps {
# look into the raw twinscan file (now containing all of the results) and
# remove redundancies caused by the overlap of the split files.

    my $raw = shift;
    my $base = shift;

    # This will be a two pass process, the first pass creates a structure
    # similar to the one created by a query in the historical script.  This
    # structure is then examined and entries to be removed are flagged.  The
    # second pass simply copies every line unless the gene_id matches one
    # that has been flagged for deletion.

    my @data = ();
    my %prediction = ();
    my %deleted = ();
    my $cur_pred = '';

    open (RAW,"< $raw") || die "Problem opening $raw: $!\n";

    while (<RAW>) {

        next if /^#/;

        chomp;

        # get us some values to work with
        my @line = split(/\t/,$_);

        my ($lend, $rend, $ori, $group) = @line[3,4,6,8];

        # make sure the lend and rend are in order:
        ($lend, $rend) = sort {$a<=>$b} ($lend,$rend);
        my $length = $rend - $lend + 1;  # (not in interbase yet)

        # is it time to store this group? Do so and set up initial range for this
        # gene
        if ($cur_pred ne $group) {
            push (@data, {%prediction}) unless ($cur_pred eq '');
            %prediction = ();
            $cur_pred = $group;
            %prediction = ('group'  => $group,
                        'orient' => $ori,
                        'lend'   => $lend,
                        'rend'   => $rend,
                        'length' => $length,
                        'deletion_target' => 0,
                       );
        } else {
            # otherwise, alter lend, rend, and length as necessary
            if ($lend < $prediction{'lend'}) {
                $prediction{'length'} += ($prediction{'lend'} - $lend);
                $prediction{'lend'} = $lend;
            }
            if ($rend > $prediction{'rend'}) {
                $prediction{'length'} += ($rend - $prediction{'rend'});
                $prediction{'rend'} = $rend;
            }
        }
    }
    
    # handle that last group that won't get pushed on otherwise:
    push (@data, {%prediction});

    close RAW;

    # We're done building the structure, let's examine it now.

    # sort on length
    @data = sort {$a->{'lend'}<=>$b->{'lend'}} @data;


    # below is almost verbatim the code that does this in the historical script
    my $prev_data = $data[0];

    for (my $i = 1; $i <= $#data; $i++) {

        my $next_data = $data[$i];

        my ($prev_lend, $prev_rend, $prev_orient) = ($prev_data->{'lend'},
                                                     $prev_data->{'rend'},
                                                     $prev_data->{'orient'});

        my ($next_lend, $next_rend, $next_orient) = ($next_data->{'lend'},
                                                     $next_data->{'rend'},
                                                     $next_data->{'orient'});

        if ($prev_orient eq $next_orient && $next_lend < $prev_rend) {
        ## ah!  an overlap

            # mark the shorty for deletion
            if ($next_data->{'length'} <= $prev_data->{'length'}) {
                $deleted{$next_data->{'group'}} = 1;
            } else {
                $deleted{$prev_data->{'group'}} = 1;
            }

        }

        # Set up for the next comparison
        unless (exists $deleted{$next_data->{'group'}}) {
            $prev_data = $next_data;
        }

    }

    # Step three, the second pass, printing out only rows where the last field hasn't
    # been marked as deleted.

    open (RAW, "< $raw") || die "Couldn't open $raw: $!\n";
    open (TMP, "> $raw.tmp") || die "Couldn't open $raw.tmp: $!\n";

    while (<RAW>) {

        if (/^#/) {

            # all comments get printed
            print TMP;

        } else {

            # for other lines, screen against the 'deleted' gene predictions
            # also, make sure the source name is corrected. 
            chomp;
            my @line = split(/\t/,$_);

            # This is an evil hack that I'm too lazy to bother with doing as
            # it should be at the moment:
            ($line[0] = $base) =~ s/\.\d+\.twinscan//;
            $_ = join ("\t",@line);
            # On second thought, maybe it's just fine.

            print TMP "$_\n" unless (exists $deleted{$line[8]});

        }

    }

    close RAW;
    close TMP;

    rename "$raw.tmp", $raw;

}
