#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;


my $options = parse_options( );
my( $header, $taxa_info ) = gather_info( $$options{ 'taxonomy_file' } );
print_taxa_of_reads( $$options{ 'input_file' }, $$options{ 'output_file' }, $taxa_info, $header );
exit $?;


sub print_taxa_of_reads {
	my( $input_file, $output_file, $info, $header ) = @_;
	open( FH, "<$input_file" ) or die "Error in opening the file, $input_file, $!\n";
	open( OFH, ">$output_file" ) or die "Error in writing to the file, $output_file, $!\n";
	my @head = split( "\t", $header );
	print OFH join( "\t", ( "query", "subject", @head[ 1 .. scalar @head -1 ] ) ),"\n";
	while( my $line = <FH> ) {
		chomp $line;
		my( $query_id, $sub_id ) = split( "\t", $line );
		unless( exists $$info{ $sub_id } ) {
			print STDERR "ERROR: No taxonomy info found for '$sub_id'\n";
		} else {
			print OFH join( "\t", ( $query_id, $sub_id, @{ $$info{ $sub_id } } ) ),"\n";
		}
	}
	close FH or die "Error in closing the file, $input_file, $!\n";
	close OFH or die "Error in closing the file, $output_file, $!\n";
}

sub gather_info {
	my( $file ) = @_;
	my $info = {};
	open( FH, "<$file" ) or die "Error in opening the file, $file, $!\n";
	my $header = <FH>;
	chomp $header;
	while( my $line = <FH> ) {
		chomp $line;
		my @array = split( "\t", $line );
		$$info{ $array[0] } = [ @array[ 1 .. scalar @array - 1 ] ];
	}
	close FH or die "Error in closing the file, $file, $!\n";
	return ( $header, $info );
}

sub parse_options {
	my $options = {};
	GetOptions( $options, 'input_file|i=s', 'taxonomy_file|t=s', 'output_file|o=s', 'help|h' ) || pod2usage( );
	if( $$options{ 'help' } ) {
		pod2usage( { 'exitval' => 1, 'verbose' => 2, 'output' => \*STDOUT } );
	}
	return $options;
}
