#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Storable;

my $options = { };
my $info = { };

GetOptions( $options, 'input_file_list|l:s', 'input_file|i:s', 'output_file|o=s' );

unless( $$options{ 'input_file_list' } || $$options{ 'input_file' } ) {
	die "options: <input_file_list> or <input_file> must be provided\n";	
}

if( $$options{ 'input_file_list' } ) {
	foreach my $file( get_file_list( $$options{ 'input_file_list' } ) ) {
		process_file( $file );		
	}
} elsif( $$options{ 'input_file' } ) {
	process_file( $$options{ 'input_file' } );
}

store_info( $info, $$options{ 'output_file' } );

exit $?;


sub store_info {
	my( $info, $file ) = @_;
	store( $info, $file ) or die "Error in storing the binary info to file, $file, $!\n";
}

sub process_file {
	my( $file ) = @_;
	open( FH, "<$file" ) or die "Error in opening the file, $file, $!\n";
	while( my $line = <FH> ) {
		if( $line =~ /^>(\S+)/ ) {
			$$info{ $1 } = $file;
		}
	}
	close FH or die "Error in closing the file, $file, $!\n";
}

sub get_file_list {
	my( $list_file ) = @_;
	my @files;
	open( FH, "<$list_file" ) or die "Error in opening the file, $list_file, $!\n";
	while( my $file = <FH> ) {
		chomp $file;
		if( -e $file ) {
			push @files, $file;
		} else {
			print STDERR "$file doens't exists\n";
		}
	}
	return @files;
}
