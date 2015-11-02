#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Storable;

my $options = {};
GetOptions( $options, 'input_file|i=s', 'dataset_info|d=s', 'output_file|o=s' );
my $dataset = retrieve( $$options{ 'dataset_info' } ) or die "Error in retrieving the dataset info, $$options{ 'dataset_info' }, $!\n";

open( OFH, ">$$options{ 'output_file' }" ) or die "Error in opening the file, $$options{ 'output_file' }, $!\n";
my $root = process_input( $$options{ 'input_file' }, $dataset );
foreach my $file( keys %$root ) {
	print_query( $file );
}
close OFH or die "Error in closing the file, $$options{ 'output_file' }, $!\n";
error_check( );
exit $?;


sub error_check {
	while( my( $file, $id ) = each %$root ) {
		foreach( keys %$id ) {
			unless( $$root{ $file }{ $_ }{ 'fetched' } ) {
				print STDERR "ERROR: '$_' in '$file' IS NOT processed\n";
			}
		}
	}	
}

sub process_input {
	my( $file, $info ) = @_;
	my( $root ) = {};
	open( FH, "<$file" ) or die "Error in opening the file, $file, $!\n";
	while( my $line = <FH> ) {
		chomp $line;
		my( $query_name ) = split( "\t", $line );
		my( $dataset ) = $$info{ $query_name };
		$$root{ $dataset }{ $query_name }{ 'fetched' } = 0;
	}
	close FH or die "Error in closing the file, $file, $!\n";
	print STDERR "--------------   Returning files is over -----------------\n";
	return $root;
}

sub print_query {
	my( $file ) = @_;
	open( FH, "<$file" ) or die "Error opening the file, $file, $!\n";
	my $num = keys %{$$root{ $file }};
	print STDERR "PROCESSING FILE: $file \t Should have sequences: ", $num,"\n";
	my $count = 0;
	WHILE: while( my $line = <FH> ) {
		if( $line =~ /^>(\S+)/ ) {
			my( $key ) = $1;
			if( exists $$root{ $file }{ $key } ) {
				print OFH $line;
				$line = print_sequence();
				++$count;
				$$root{ $file }{ $key }{ 'fetched' } = 1;
				redo WHILE if( $line );		
			}
		}
	}
	if( $num != $count ) {
		print STDERR "ERROR: Number of seqs expected ( ", ( $num ),' ) != ( ', $count, " ) seqs processed\n";
	}
	close FH or die "Error closing the file, $file, $!\n";
}

sub print_sequence {
	WHILE: while( my $line = <FH> ) {
		if( $line =~ /^[^>\s]+/ ) {
			print OFH $line;
		} else {
			return $line;
		}
	}
}
