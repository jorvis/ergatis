#!/usr/bin/perl

=head1  NAME

pangenome_query_list.pl - Prepare a list of polypeptide and assembly sequence IDs from the set of pangenome BSML input files.

=head1 SYNOPSIS

USAGE: pangenome_query_list.pl
        --input_list=/path/to/somefile.bsml.list
        --output=/path/to/output.list
      [ --log=/path/to/some/log ]

=head1 OPTIONS

B<--input_list,-i>
    BSML list file containing sequences used in a pangenome analysis pipeline. 

B<--output,-o>
    Path to output file that will contain the list of query polypeptide and assembly sequence ids.

B<--log,-d>
    optional. Will create a log file with summaries of all actions performed.

B<--debug>
    optional. Will display verbose status messages to STDERR if logging is disabled.

B<--help,-h>
    This help message/documentation.

=head1   DESCRIPTION

    This script is used to extract the list of all query polypeptide and assembly sequence ids from a set of 
    input BSML files that were used in the pangenome analysis pipeline. You can provide a list of BSML
    documents representing a subset of the complete set of genomes used to perform the blast queries if desired.

=head1 INPUT

    The input should be a list of BSML files that were used as input to the pangenome analysis pipeline.

=head1 OUTPUT

    There is no output unless you use the --log option.

=head1 CONTACT

    Brett Whitty
    bwhitty@tigr.org

=cut

use Pod::Usage;
use Storable qw(nstore retrieve);
use XML::Twig;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use strict;

my %options = ();
my $results = GetOptions (  \%options,
                            'input_list|i:s',
                            'output_path|o=s',
                            'debug|d=s',
                            'command_id=s',       ## passed by workflow
                            'logconf=s',          ## passed by workflow (not used)
                            'log|l=s',
                            'help|h',
                         ) || pod2usage();

$options{'output_path'} =~ s/\/$//;

unless(-d $options{'output_path'}) {
    die "must provide output directory as argument";
}

if ($options{'input_list'} eq '') {
    ## when no input list is provided, we will not filter the results
    exit();
}

my %sequence_hash = ();

my $twig = XML::Twig->new(
                            twig_roots  => { 
                                            'Sequence' => \&processSequence
                                           }
                         );
                         
open (LIST, $options{'input_list'}) || die "couldn't open input list for reading";

while (<LIST>) {
    chomp;
    print STDERR "File: $_\n";
    if (-e $_) {
        $twig->parsefile($_);
    } else {
        die "specified input file does not exist";
    }
}

store(\%sequence_hash, $options{'output_path'}."/pangenome.filter.stored");

exit(0);

sub processSequence { 
    my ($twig, $feat) = @_;
    
    my $class = $feat->{'att'}->{'class'};
    my $id = $feat->{'att'}->{'id'};
    if ($class eq 'polypeptide') {
        my @seqdataimport = $feat->children('Seq-data-import');
        my @link = $feat->children('Link');
        if ($link[0]->{'att'}->{'rel'} ne 'analysis') {
            my $seq_id = $seqdataimport[0]->{'att'}->{'identifier'};
	    if ($seq_id =~ /^(([^\.]+)\.[a-z]+\.\d+\.\d+)$/) { 
		$sequence_hash{$2}{$1} = 1;
	    } else {
		$seq_id =~ /^([^_]+)_(.*)_[^_]+$/ || print STDERR "couldn't parse seq-data-import id for '$seq_id'\n";
		$sequence_hash{$1}{$2} = 1;
	    }
        }
    } elsif ($class eq 'assembly') {
	if ($id =~ /^(([^\.]+)\.assembly\.\d+\.\d+)$/) {
	    $sequence_hash{$2}{$1} = 1;
	} else {
	    $id =~ /^([^_]+)_(.*)_[^_]+$/ || print STDERR "couldn't parse db and seq id from '$id'\n";
	    $sequence_hash{$1}{$2} = 1;
	}
    }
}
