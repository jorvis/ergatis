#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

=head1 NAME

bowtie.pl - run the bowtie short-read aligner

=head1 SYNOPSIS

  USAGE: bowtie.pl [
            --bowtie_exec=full path to bowtie binary
            --reference=full path to bowtie reference index
            --reads=full path to reads fastq files, if paired-end then comma-separated paths to each mate
            --sam_output=full path to output sam file
            --max_insert=maximum length of insert for paired-end reads, default 300
            --max_mismatches=maxixmum number base-pairs that can mismatch, default 2
            --max_aligns=do not report reads with more than this many alignments, default 1
            --more_options=additional bowtie options appended to options area before reference
            --help
          ]
=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use POSIX;

my %options = ();
my $results = GetOptions (\%options, 
			  'bowtie_exec|bin=s',
                          'reference|r=s',
			  'reads|q=s',
                          'sam_output|o=s',
                          'max_insert|X=s',
                          'max_mismatches|v=s',
                          'max_aligns|m=s',
			  'more_options|more_options=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

my $bin = $options{'bowtie_exec'};
my $reference = $options{'reference'};
my $reads = $options{'reads'};
my $sam_output = $options{'sam_output'};
my $max_insert = $options{'max_insert'};
my $max_mismatches = $options{'max_mismatches'};
my $max_aligns = $options{'max_aligns'};
my $more_options = $options{'more_options'};

## Our list file does not directly reference the bowtie indices (these have extension of ebwt) but instead 
## contain a list of placeholder files (idx files) which can be parsed to produce the path to the reference
## that bowtie understands. 
##
## e.x. reference files = e_coli.ebwt, e_coli.1.ebwt, e_coli.2.ebwt, e_coli.rev.1.ebwt
##      idx file = e_coli.idx
##
##      (the reference passed to the bowtie executble would be /path/to/e_coli with no extension)
##
## e.x. bowtie -S /local/references/e_coli sequences_1.fastq /local/reads/sequences.fastq /local/output/sequences.sam
$reference =~ s/\.idx$//;

if($reads =~ /,/g){ #paired end reads
    my @p = split(",", $reads);
    system("$bin $more_options -X $max_insert -v $max_mismatches -m $max_aligns -S $reference -1 " . trim($p[0]) . " -2 " . trim($p[1]) . " $sam_output");

}
else{ #unpaired reads
    system("$bin $more_options -v $max_mismatches -m $max_aligns -S $reference $reads $sam_output");
}    

sub check_parameters {
    my $options = shift;
    
    ## make sure reference, reads, and output were passed
    unless ( $options{reference} && $options{reads} && $options{sam_output} ) {
        die "You must pass --reference --reads and --sam_output!";
    }
    
    ## handle some defaults
    $options{bowtie_exec} = "/usr/local/packages/bowtie/bowtie" unless ($options{bowtie});
    $options{max_insert}   = 300  unless ($options{output_subdir_size});
    $options{max_mismatches} = 2 unless ($options{output_subdir_prefix});
    $options{max_aligns}        = 1  unless ($options{seqs_per_file});
}

#-------------------------------------
# Trims whitespace from a string
#-------------------------------------
sub trim {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

