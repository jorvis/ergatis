#!/usr/bin/env perl

=head1 NAME

parse_snp_blast.pl - Description

=head1 SYNOPSIS

 USAGE: parse_snp_blast.pl
	--input_file=/path/to/genome_file.m0
	--output_file=/path/to/output.m0.parsed

=head1 OPTIONS

B<--input_file,-i>

B<--output_file,-o>

=head1  DESCRIPTION

 DESCRIPTION
 
=head1  INPUT
    Describe the input

=head1 OUTPUT
    Describe the output

=head1  CONTACT

    Kevin Galens
    kgalens@gmail.com

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Find::Rule;
use File::Basename;
use Pod::Usage;

############# GLOBALS AND CONSTANTS ################
####################################################

my %options;
my $results = GetOptions (\%options,
		"input_dir|i=s",
		"file_extension|e=s",
		"output_dir|o=s",
		"help|h"
		);

my $file_ext;
&check_options(\%options);
my $ext = qr/\.$file_ext$/;  
my @blast_files = File::Find::Rule->file()->name($ext)->in($options{'input_dir'});
print "@blast_files\n";
if(@blast_files > 0) {
	foreach my $file(@blast_files) {
		if (-e $file) {
			&parse_blast($file);		
		}
	}
} else {
	die "No BLAST results (*.$options{'file_extension'} files) found in the input directory $options{'input_directory'}\n";
}

sub parse_blast {
	my ($raw_file) = @_;
	open(IN, "< $raw_file") or die("Could not open file $raw_file for reading.\nReason :$!\n");
	my ($file_base,$file_dir,$file_ext) = fileparse($raw_file,qr/\.[^.]*/);
	my $parsed_file = $options{'output_dir'}."/".$file_base.".parsed";
	open(OUT, "> $parsed_file") or die("Could not open file $parsed_file for writing.\nReason : $!\n");
	my $snpid;
	my $flag = 0;
	my @hit = ();
	my $no_print = 0;
	while( my $line = <IN> ) {
		next if( $line =~ /^\s*$/ );
		if( $line =~ /Query=/ ) {
#			print "$line";
			if( @hit ) {
				unless( $no_print ) {
#					print "printing: $hit[0]";
					map { print OUT $_ } @hit;
				} else {
					$no_print = 0;
				}
				@hit = ();
			}
			push(@hit, $line);
		} elsif( $line =~ /^Query:/ || $line =~ /^[\s\|]+$/ || $line =~ /^Sbjct:/ ) {
			push(@hit, $line);
		} elsif( $line =~ /No hits found/ ) {
			$no_print = 1;
		}
	}
	close(IN);
	unless( $no_print ) {
#		print "lastly, printing: $hit[0]";
		map { print OUT $_ } @hit;
		close(OUT);
	}
}

sub check_options {

	my $opts = shift;

	if( $opts->{'help'} ) {
		&_pod;
	}

	foreach my $req ( qw(input_dir output_dir) ) {
		die("Option $req is required") unless( $opts->{$req} );
	}
	if(exists($opts->{'file_extension'})) {
		$file_ext = $opts->{'file_extension'};
	} else {
		$file_ext = "raw";
	}
}

sub _pod {
	pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
