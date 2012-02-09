#!/usr/bin/perl

=head1 NAME

filter_sequences_by_id.pl - Filter a sequence file to include/exclude entries by their IDs.

=head1 SYNOPSIS

USAGE: split_multifasta.pl 
            --input_file=/path/to/some_file.fsa 
            --output_dir=/path/to/somedir
            --id_file=/path/to/somefile.tab
          [ --mode=include
            --format=fasta
            --id_file_list=/path/to/some.list
            --id_column_num=1
          ]

=head1 OPTIONS

B<--input_file,-i>
    The input file (parallelization available in Ergatis).

B<--output_file,-o>
    The output file to be created/overwritten.

B<--id_file,-f>
    This file contains the IDs on which the script will filter, and can either have
    one value per line or be tab-delimited, in which case you'll need to also pass
    the --id_column_num option.

B<--mode,-m>
    Optional.  Define whether to include/exclude only those sequences in the ID list in 
    the output file.  (default = include)

B<--format,-r>
    Optional.  Format of input files to filter.  Must be either 'fasta' or 'fastq' (default = fasta)

B<--id_file_list,-l>
    Optional.  Use this if you have more than one file containing filter IDs.

B<--id_column_num,-n>
    Optional.  In which column of the tab-delimited input files are the IDs?  (Default = 1)

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

This script is used to filter entries in a FASTA or FASTQ file by their IDs.  Depending on the
value of the --mode setting, the resulting output files can include or exclude those sequences
matching the ID list.

=head1  INPUT

This script operates on a single file.  In practice, it's called as part of an Ergatis component
by the same name if you have multiple input files to process.

The input is defined with --input_file and should be a single file in one of the supported formats.
File extensions are ignored.  The first *word* of each header line is accepted to be the ID

For example:

    >gi53791237 Tragulus javanicus p97bcnt gene for p97Bcnt
    ACAGGAGAAGAGACTGAAGAGACACGTTCAGGAGAAGAGCAAGAGAAGCCTAAAGAAATGCAAGAAGTTA
    AACTCACCAAATCACTTGTTGAAGAAGTCAGGTAACATGACATTCACAAACTTCAAAACTAGTTCTTTAA
    AAAGGAACATCTCTCTTTTAATATGTATGCATTATTAATTTATTTACTCATTGGCGTGGAGGAGGAAATG

    @HWI-ST180:243:C0589ACXX:1:1101:16202:57975 1:N:0:CGATGCA
    TTAGCCTTTTTTGCTTCCTTGGCAGCCCTGATAGCTTGTTCTCGTTGAGCCTTCCTAACTTCAGGTTTCTGATTCCTCTTGGCCATTATATCAGCAAGAGA
    +
    CCCFFFFAHGHHHJIIJJJJGEGGIJJJIGDHICFGDFDGIDGIFEHHIJJEGHIIGGIJIHCEE;AEHHHBECBC>DCEEEAEDCDCDEACAACCCCC?C

Whitespace is ignored within the input file.  

=head1  OUTPUT

A single output file will be created.  Depending on the value passed to the --mode option, 
the output file will either include/exclude only those with the IDs passed.

=head1  CONTACT

    Joshua Orvis
    jorvis@users.sf.net

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use POSIX;
use Ergatis::Logger;

my %options = ();
my $results = GetOptions (\%options, 
                          'input_file|i=s',
                          'output_file|o=s',
                          'id_file|f=s',
                          'mode|m=s',
                          'format|r=s',
                          'id_file_list|l=s',
                          'id_column_num|n=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

my $ids = load_ids( $options{id_file}, $options{id_file_list}, $options{id_column_num} );

my $sfh;

## load the sequence file
if ($options{'input_file'} =~ /\.(gz|gzip)$/) {
    open ($sfh, "<:gzip", $options{'input_file'})
      || $logger->logdie("can't open sequence file:\n$!");
} else {
    open ($sfh, "<$options{'input_file'}")
      || $logger->logdie("can't open sequence file:\n$!");
}

open(my $ofh, ">$options{output_file}") || $logger->logdie("failed to open output file: $!");

if ( $options{format} eq 'fasta' ) {
    process_fasta( $sfh, $ofh, $ids );
    
} elsif ( $options{format} eq 'fastq' ) {
    process_fastq( $sfh, $ofh, $ids );

} else {
    $logger->logdie("The only currently supported formats are fasta and fastq");
}




exit;


sub load_ids {
    my ( $id_file, $id_file_list, $id_column_num ) = @_;
    
    my $ids_loaded = {};
    
    if ( $id_file ) {
        load_id_file( $id_file, $ids_loaded, $id_column_num );
    }
    
    if ( $id_file_list ) {
        open(my $list_fh, $id_file_list) || $logger->logdie( "Failed to read input ID list file" );
        
        while ( my $line = <$list_fh> ) {
            chomp $line;
            next if $line =~ /^\s*$/;
            
            load_id_file( $line, $ids_loaded, $id_column_num );
        }
    }
    
    $logger->info("Loaded " . scalar(keys %$ids_loaded) . " IDs to filter\n");
    
    return $ids_loaded;
}

sub load_id_file {
    my ($file, $dmap, $col_num) = @_;
    
    open(my $ifh, $file) || $logger->logdie("Failed to read ID file ($file): $!");
    
    while (my $line = <$ifh>) {
        chomp $line;
        
        if ( $col_num ) {
            my @cols = split("\t", $line);
            $$dmap{ $cols[$col_num - 1] } = 1;
        } else {
            $$dmap{$line} = 1;
        }
    }
}

sub process_fasta {
    my ( $in, $out, $ids ) = @_;
    
    my $keep_this_seq = 0;
    my $seq_count = 0;
    
    while (my $line = <$in>) {

        ## if we find a header line ...
        if ($line =~ /^\>(\S+)/) {
            $seq_count++;
            
            if ( $options{mode} eq 'include' ) {
                if ( exists $$ids{$1} ) {
                    $keep_this_seq = 1;
                } else {
                    $keep_this_seq = 0;
                }
            } elsif ( $options{mode} eq 'exclude' ) {
                if ( exists $$ids{$1} ) {
                    $keep_this_seq = 0;
                } else {
                    $keep_this_seq = 1;
                }
            }
        }
        
        print $out $line if $keep_this_seq;
    }
}

sub process_fastq {
    my ( $in, $out, $ids ) = @_;
    
    my $keep_this_seq = 0;
    my $line_num = 0;
    
    while (my $line = <$in>) {
        $line_num++;
        
        ## sanity check.  there should be 4 lines per entry.
        #   the first of each should start with @, and the 3rd the + symbol
        if ( $line_num % 4 == 1 ) {
            if ( $line =~ /^\@(\S+)/ ) {
                
                if ( $options{mode} eq 'include' ) {
                    if ( exists $$ids{$1} ) {
                        $keep_this_seq = 1;
                    } else {
                        $keep_this_seq = 0;
                    }
                } elsif ( $options{mode} eq 'exclude' ) {
                    if ( exists $$ids{$1} ) {
                        $keep_this_seq = 0;
                    } else {
                        $keep_this_seq = 1;
                    }
                }
                
            } else {
                $logger->logdie("Error in source FASTQ file.  Expected @ symbol to begin line $line_num");
            }
            
        } elsif ( $line_num % 4 == 3 ) {
            if ( $line !~ /^\+/ ) {
                $logger->logdie("Error in source FASTQ file.  Expected + symbol to begin line $line_num");
            }
        }
        
        print $out $line if $keep_this_seq;
    }
}

sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( output_file );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }

    ## make sure input_file or input_file_list were passed
    unless ( $options{input_file} || $options{input_file_list} ) {
        $logger->logdie("You must pass either input_file or input_file_list");
    }

    ## make sure either id_file or id_file_list were passed
    unless ( $options{id_file} || $options{id_file_list} ) {
        $logger->logdie("You must pass either id_file or id_file_list");
    }
    
    ## make sure input_file exists
    if (! -e $options{input_file} ) {
        if ( -e "$options{input_file}.gz" ) {
            $options{input_file} .= '.gz';
        } else {
            $logger->logdie("the input file passed ($options{input_file}) cannot be read or does not exist");
        }
    }
    
    ## handle some defaults
    $options{mode} = 'include' unless ($options{mode});
    $options{id_column_num} = 1 unless ($options{id_column_num});
    $options{format} = 'fasta' unless ($options{format});

    ## mode should be either 'include' or 'exclude'
    if ( $options{mode} ne 'include' && $options{mode} ne 'exclude' ) {
        $logger->logdie("The value for --mode should be either 'include' or 'exclude'");
    }

    ## the only supported formats are fasta and fastq
    if ( $options{format} !~ /^fast[aq]$/ ) {
        $logger->logdie("The value for --format should be either 'fasta' or 'fastq'");
    }
}

