#!/usr/bin/env perl

=head1 NAME

create_fastq_contigs_map.pl

=head1 SYNOPSIS

 USAGE: create_fastq_contigs_map.pl
       --input_map=/path/to/celera_or_velvet.map
       --fastq_list=/path/to/fastq.list
       --output_map=/path/to/output.map
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1  DESCRIPTION

    This relies on how filter_fastq names the output files it creates. If that
    method changes, this script should be changed.

=head1  CONTACT

    Kevin Galens
    kgalens@gmail.com

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Basename;

use Data::Dumper;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
my $format;
my ($VELVET, $CELERA) = (0,1);
####################################################

my %options;
my $results = GetOptions (\%options,
			  "input_map|i=s",
			  "fastq_list|f=s",
			  "output_map|o=s",
			  "celera_input_iter_list|c=s",
			  "log|l=s",
			  "debug|d=s",
			  "help|h"
                          );

&check_options(\%options);

my %map = &parse_map( $options{'input_map'} );

my @pairs;
if( $format eq $VELVET ) {
    @pairs = &parse_fastq_list( $options{'fastq_list'}, \%map );
} elsif( $format eq $CELERA ) {
    @pairs = &parse_celera_pairs( $options{'fastq_list'}, $options{'celera_input_iter_list'}, \%map );
}

open(OUT, "> $options{'output_map'}") or die("Couldn't open $options{'output_map'} for writing: $!");
print OUT join("\t",@{$_})."\n" for( @pairs );
close(OUT);

sub parse_celera_pairs {
    my ($flist, $ciil, $map) = @_;

    my %iter_map;
    open(IN, "< $ciil") or die("Can't open $ciil: $!");
    my $burn = <IN>;
    while( my $line = <IN> ) {
	chomp($line);
	my @c = split(/\t/, $line);
	die("Could not find gkpstore $c[2] in input_map [$options{'input_map'}]") 
	    unless( exists( $map->{$c[2]} ) );
	$iter_map{$c[0]} = $c[2];
    }
    close(IN);

    my @pairs;
    open(IN, "< $flist") or die("Can't open $flist: $!");
    while( my $fastq = <IN> ) {
	chomp($fastq);
	my $basename = basename( $fastq, qw(.all.fastq) );
	die("Could not find $basename from fastq file in iter list lookup") 
	    unless( exists( $iter_map{$basename} ) );
	my $gkpstore = $iter_map{$basename};
	my $contigs = $map{$gkpstore};
	push(@pairs, [$contigs, $fastq]);
    }
    close(IN);

    return @pairs;
}

sub parse_fastq_list {
    my ($list, $map) = @_;
    my @pairs;

    open(IN, "< $list") or die("Can't open $list: $!");
    while( my $l = <IN> ) {
	chomp($l);

	my $list_basename = basename( $l, qw(.fastq.list) );
	if( exists( $map{ $list_basename } ) ) {
	    push(@pairs, [$map->{$list_basename},$l]);
	} else {
	    print Dumper( $map );
	    die("Couldn't find mapping for fastq list $l [$list_basename]");
	}

    }
# 	open(LIST, "< $l") or die("Can't open $l: $!");
# 	chomp( my @fastqs = <LIST> );
# 	close(LIST);
	
# 	my @bnames;
# 	map { push(@bnames,basename($_)) } @fastqs;

# 	my $lcp = &lcp( @bnames );
# 	if( exists( $map->{$lcp} ) ) {
# 	    push(@pairs, [$map->{$lcp},$l]);
# 	} else {
# 	    die("Couldn't find mapping for fastq list $l");
# 	}
#     }
    close(IN);
    return @pairs;
}

sub parse_map {
    my ($map_file) = @_;
    my %map;

    my %subs = (
		'store_velvet' => sub {
		    my ($map, $contigs, $fastqs, $read_names) = @_;
		    my @f = split(/,/, $fastqs);
		    map { $_ = basename( $_, qw(.txt .fastq.gz .fastq .fq.gz .fq) ) } @f;
		    #my $lcp = &lcp( @f );
		    #$map->{$lcp} = $contigs;
		    $map->{$f[0]} = $contigs;
		},
		'store_celera' => sub { 
		    my ($map, $contigs, $gkpstore) = @_;
		    $map->{$gkpstore} = $contigs;
		});

    open(IN, "< $map_file") or die("Can't open $map_file: $!");
    my $store;
    while( my $line = <IN> ) {
	chomp($line);
	my @c = split(/\t/, $line);

	if( defined( $store ) ) {
	    $store->(\%map, @c);
	} else {
	    if( @c == 3 ) {
		$format = $VELVET;
		$store = $subs{'store_velvet'};
	    } elsif( @c == 2 ) {
		$format = $CELERA;
		$store = $subs{'store_celera'};
	    } else {
		die("Incorrect number of columns in map file. Expected 2 or 3, got ".scalar(@c));
	    }
	    $store->(\%map, @c);
	}

    }
    close(IN);

    return %map;
}

#longest common prefix
sub lcp {
    my ($s1, $s2) = @_;
    return $s1 unless( defined( $s2 ) );
    chop($s1) while( $s2 !~ /^\Q$s1\E/ );
    return $s1;
}

sub check_options {

   my $opts = shift;

   if( $opts->{'help'} ) {
       &_pod;
   }

   if( $opts->{'log'} ) {
       open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
   }

   $debug = $opts->{'debug'} if( exists( $opts->{'debug'} ) );

   foreach my $req ( qw(input_map fastq_list output_map) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }
}

sub _log {
   my ($level, $msg) = @_;
   if( $level <= $debug ) {
      print STDOUT "$msg\n";
   }
   print $logfh "$msg\n" if( defined( $logfh ) );
   exit(1) if( $level == $ERROR );
}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
