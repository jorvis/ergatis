#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Storable;

my $options = get_options( );
my $mapping_info = retrieve( $$options{ 'read_dataset_info' } ) or die "Error in retrieving the info from $$options{ 'read_database_info' }, $!\n";
my ( $info, $header ) = get_info( $$options{ 'input_file' }, $mapping_info );
print_info( $info, $header, $$options{ 'output_file' } );
exit $?;

sub print_info {
	my ( $info, $header, $output_file ) = @_;
	open( OFH, ">$output_file" ) or die "Error in writing to the file, $output_file, $!\n";
	print OFH join( "\t", ( 'dataset', @$header ) ), "\n";
	foreach my $dataset ( sort {$a cmp $b} keys %$info ) {
		print OFH $dataset;
		foreach my $resprot ( @$header ) {
			my $count = ( exists $$info{ $dataset }{ $resprot } ) ? $$info{ $dataset }{ $resprot } : '0';
			print OFH "\t", $count;
		}
		print OFH "\n";
	}
	close OFH;
}

sub get_info {
	my ( $file, $mapping_info ) = @_;
	my ( $info, $header );
	open( FH, "<$file" ) or die "Error in opening the file, $file, $!\n";
	while( my $line = <FH> ) {
		chomp $line;
		my ( @values ) = split( "\t", $line );
		my ( $query_id, $source ) = ( $values[0], $values[ scalar@values - 1 ] );
		my $dataset = $$mapping_info{ $query_id };
		unless( $dataset ) {
			print STDERR "There is no mapping info found for $query_id\n";
		}
#		if( $query_id =~ /^(.+_)(.+)/ ) {
#			$dataset = $2;
#		}
		unless( exists $$header{ $source } ) {
			$$header{ $source }++;
		}
		$$info{ $dataset }{ $source }++;
	}
	close FH;
	my @header_keys = sort { $a cmp $b } keys %$header;
	return ( $info, \@header_keys );
}

sub get_options {
	my %options;
	GetOptions( \%options, 'input_file|i=s', 'output_file|o=s', 'read_dataset_info|r=s', 'help|h' ) || pod2usage( );
	if( $options{ 'help' } ) {
		pod2usage( { -exitval => 0, -verbose => 1, -output => \*STDOUT } );
	}
	return \%options;
}

