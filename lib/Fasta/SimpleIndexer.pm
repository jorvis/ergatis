#!/usr/local/bin/perl

package Fasta::SimpleIndexer;

=head1 NAME

  SimpleIndexer.pm - Provides lightweight object for indexed reading from a fasta file

=head1 Description

=head2 Overview

  Makes use of the BSML library BSML::Indexer::Fasta class to facilitate indexed reads
  from fasta files. This class is intentionally lightweight to be used wherever reading
  from fasta files is required outside of the scope of BSML document generation.

=head2 Constructor and initialization

  The object interface has been intentionally simplified to support and encourage wide use
  of this module in ergatis scripts.

  The constructor can be called in any of the following ways:

  my $fasta = new Fasta::SimpleIndexer( "/path/to/some/fasta_file.fsa" );
                       *or*
  my $fasta = new Fasta::SimpleIndexer( fasta_file => "/path/to/some/fasta_file.fsa" );
                       *or*                    
  my $fasta = new Fasta::SimpleIndexer( 
                                        fasta_file => "/path/to/some/fasta_file.fsa",
                                        index_dir  => "/path/to/another/dir"
                                      );

  If 'index_dir' is omitted, File::Basename::dirname(fasta_file) will be used.
                                      
=head2 Class and object methods

  get_record(accession)

    returns the formatted fasta record for the provided accession

  get_sequence(accession)
    
    returns the raw sequence string for the provided accession
  
=over 4

=cut

use warnings;
use strict;

use BSML::Indexer::Fasta;
use File::Basename qw( dirname );
use Carp;

my ($index_dir, $fasta_file);
my (%h, %e);

sub new {
    my ($class) = shift(@_);
    ## handle arguments, supporting use cases:
    ##
    ## new Fasta::SimpleIndexed("/path/to/some/fasta_file.fsa");
    ## new Fasta::SimpleIndexed(fasta_file => "/path/to/some/fasta_file.fsa");
    ## new Fasta::SimpleIndexed(fasta_file => "/path/to/some/fasta_file.fsa", index_dir => "/path/to/index/dir");
    ##
    my %args = ();
    if (scalar(@_) % 2 == 0) {
       (%args) = @_;
    } elsif (scalar(@_) == 1) {
       ($fasta_file) = @_;
       $args{'fasta_file'} = $fasta_file;
    } else {
        confess "Malformed arguments passed to constructor";
    }
    
    if (! $args{'fasta_file'}) {
        confess "Constructor requires 'fasta_file' as a parameter";
    } else {
        $fasta_file = $args{'fasta_file'};
    }
    if (! $args{'index_dir'}) {
        $index_dir = dirname($fasta_file);
    }
    unless (-d $index_dir) {
        confess "index_dir specified as '$index_dir' is not a directory!";
    }
    unless (-e $fasta_file) {
        confess "fasta_file specified as '$fasta_file' does not exist!";
    }

    my $indexer = BSML::Indexer::Fasta->new($fasta_file, $index_dir);

    # Check the health of the various indices for the data file.
    my @check = $indexer->check_indices;

    # Create the indices if necessary...
    if ($check[0] == 1) { $indexer->index_entries };
    if ($check[1] == 1) { $indexer->index_headers };

    # Get the name of the header index file.
    my $header_index = $indexer->header_index;
    my $entry_index = $indexer->entry_index;

    my $c = tie %h, 'CDB_File', $header_index or die "tie failed: $!\n";
       $c = tie %e, 'CDB_File', $entry_index or die "tie failed: $!\n";

    return bless {}, $class;  
}

## returns a formatted fasta record for a given accession 
sub get_record {
    my ($self, $accession) = @_;
    
    unless ($accession) {
        confess "No accession provided to get_record";
    }
    
    my $ref = _fetch($accession);
    
    $ref->[0] .= "\n";
    $ref->[1] =~ s/(.{1,60})/$1\n/g;
    
    return ">".join("", @{$ref});
}

## returns a string containing the raw sequence for a specified accession
sub get_sequence {
    my ($self, $accession) = @_;

    unless ($accession) {
        confess "No accession provided to get_sequence";
    }
    
    my $ref = _fetch($accession);

    return $ref->[1];
}

## returns an array containing header string and string containing the raw sequence for a specified accession
sub get_sequence_array {
    my ($self, $accession) = @_;

    unless ($accession) {
        confess "No accession provided to get_sequence_array";
    }
    
    my $ref = _fetch($accession);

    return ($ref->[0], $ref->[1]);
}

## private method that does the fetching
sub _fetch {
    
    my ($specified_header) = @_;
    
    my $entry = $h{$specified_header};

    my ($offset, $length) = split( ',', $e{$entry} );

    # if the sequence could not be loaded with an index, perform a grep operation on the
    # FASTA file.

    open (IN, $fasta_file) or die "Unable to open $fasta_file due to $!";

    if ($offset) {
        seek(IN, $offset, 0);
    }

    my $defline = '';
    my $line = <IN>;
    my $seq_ref = [];
    while (defined($line)) {
        unless($line =~ /^>([^\s]+)/) {
            $line = <IN>;
        } else {
            my $header = $1;
            
            $defline = $line;
            chomp($defline);
            $defline =~ s/^>//;
            
            if ($header eq $specified_header) {
                while(defined($line=<IN>) and $line !~ /^>/ ) {
                    next if($line =~/^\s+$/);                   #skip blank lines
                    chomp($line);
                    push(@$seq_ref, $line);
                }
                last;   #seq found, terminating fasta_file parasing
            } else { $line = <IN>; };  #wrong seq, keep looking
        }
    }
    close IN;

    my $final_seq = join("", @$seq_ref);

    return [$defline, $final_seq];
}

1;
