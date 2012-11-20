#!/usr/bin/env perl

=head1 NAME

create_fasta_bam_map.pl - Will create an output list mapping input query fasta to output bam files. Written to be
    used as part of the BWA component.

=head1 SYNOPSIS

 USAGE: create_fasta_bam_map.pl
    --query_list=/path/to/input.list
    --bam_list=/path/to/output.list
    --output=/path/to/output.map
    --strict=1
    [ 
    --log=/path/to/file.log
    --debug=3
    --help
    ]

=head1 OPTIONS

B<--query_list,-q>
    Input list of query and reference sequences used as input to BWA. If this option is blank, 
    the script will not run and just return 0. 

B<--bam_list,-b>
    The list of BAM files. The list which is created by bwa component.

B<--output,-o>
    The output map file. Can be used as input to other components. Such as mpileup, for example.

B<--strict,-s>
    Will require that all files in the query_list file map to files in the bam_list (and vice versa).

B<--log,-l>
    Logfile.

B<--debug,-d>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

    This was written to create an output map file for use as input to other components. Developed for use
    at the end of the BWA pipeline and the map file can be used as the input map file for mpileup. Might
    have other uses as well down the line.

=head1  CONTACT

    Kevin Galens
    kgalens@gmail.com

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Basename;
use File::OpenFile qw(open_file);

use Data::Dumper;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
my $strict = 1;
my @possible_query_extensions = qw(.fa .fsa .fasta .fna .fas .contigs);
####################################################

my %options;
my $results = GetOptions (\%options,
			  "query_list|q=s",
			  "bam_list|b=s",
			  "output_map|o=s",
			  "strict|s=s",
			  "log|l=s",
			  "debug|d=s",
			  "help|h"
                          );

&check_options(\%options);

my %queries = &parse_query_list( $options{'query_list'} );
my @pairs = &map_bam_files( $options{'bam_list'}, \%queries, $strict );

my $ofh = open_file( $options{'output_map'}, 'out' );
foreach my $p ( @pairs ) {
    print $ofh join("\t", @{$p})."\n";
}
close($ofh);


sub map_bam_files {
    my ($bam_list, $queries, $strict) = @_;
    my $bfh = open_file( $bam_list, 'in' );
    chomp( my @bam_files = <$bfh> );
    close($bfh);
    
    my %leftover_queries = %{$queries};

    my @pairs;
    for( @bam_files ) {
	my $basename = basename( $_, qw(.bam) );

	my $flag = 0;
	$flag = 1 if( $basename =~ /VC_18/ );

	my $query_match;
	foreach my $qname ( keys %{$queries} ) {
	    next unless( $basename =~ /$qname[\._]/ );
	    $query_match = $queries->{$qname};
	    delete( $leftover_queries{$qname} ) if( exists( $leftover_queries{$qname} ) );
	    last;
	}
	if( defined( $query_match ) ) {
	    push(@pairs, [$query_match,$_]);
	} else {
	    &_log($WARN, "Couldn't find matching query for bam file $_");
	}

    }

    if( keys %leftover_queries ) {
	&_log($WARN, Dumper( \%leftover_queries ) );
	&_log($WARN, "Couldn't find matching bam files for the previous queries");
    }

    return @pairs;
	
}

sub parse_query_list {
    my ($qlist, $s) = @_;
    my $qfh = open_file( $qlist, 'in' );
    chomp( my @qfiles = <$qfh> );
    close($qfh);

    my %retval;
    for( @qfiles ) {
	next if( /^\s*$/ );
	my $basename = basename( $_, @possible_query_extensions );
	$retval{$basename} = $_;
    }

    return %retval;
}


sub check_options {
   my $opts = shift;

   &_pod if( $opts->{'help'} );
   (open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)")) if( $opts->{'log'} );

   $debug = $opts->{'debug'} if( exists( $opts->{'debug'} ) );
   $strict = $opts->{'strict'} if( exists( $opts->{'strict'} ) );

   exit(0) unless( $opts->{'query_list'} );

   foreach my $req ( qw(bam_list output_map) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }
}

sub _log {
   my ($level, $msg) = @_;
   if( $level <= $debug ) {
       print STDOUT "$msg\n"; 
       print $logfh "$msg\n" if( defined( $logfh ) );
   }
   exit(1) if( $level == $ERROR );
}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
