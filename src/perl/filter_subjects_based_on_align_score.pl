#!/usr/bin/perl
=head1 NAME
	some name

=head1 SYNOPSIS
	Usage: $0 --input-file-list|l <file with paths for all files> --output_file|o  [--help|h]

=head1 Descriptions

=head1 AUTHOR

	Mahesh Vangala
	mvangala@som.umaryland.edu
	vangalamaheshh@gmail.com

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my $info = {};

my $options = get_options( );
parse_file_list( $$options{ 'input_file_list' } );
get_difference( );
print_info( $$options{ 'output_file' } );
exit $?;


sub print_info {
	my( $file ) = @_;
	open( OFH, ">$file" ) or die "Error in writing to the file, $file, $!\n";
	while( my( $k, $v ) = each %$info ) {
		my( $match, $score ) = ( undef, 0 );
		while( my( $k2, $v2 ) = each %$v ) {
			if( $score < $$info{ $k }{ $k2 }{ 'total_align_length' } ) {
				$score = $$info{ $k }{ $k2 }{ 'total_align_length' };
				$match = $k2;
			}
		}
		print OFH join( "\t", ( $k, $match, $score ) ), "\n";
	}
	close OFH or die "Error in closing the file, $file, $!\n";
}

sub get_difference {
	while( my( $k, $v ) = each %$info ) {
		while( my( $k2, $v2 ) = each %$v ) {
			my @sorted = sort by_num @{ $$info{ $k } { $k2 } { 'array' } };
			$$info{ $k } { $k2 } { 'total_align_length' } = difference( @sorted );
		}
	}
}

sub difference {
	my( @array ) = @_;
	my( $start_flag, $end_flag, $total, $start_num, $counter ) = ( 1, 0, 0, 0, 0 );
	foreach my $val ( @array ) {
		my( $token, $num ) = ( $val =~ /(\w+)\W+(\d+)/ );
		$token eq 'start' ? $counter++ : $counter--;
		if( $start_flag && $counter == 1 ) {
			$start_num = $num;
			( $start_flag, $end_flag ) = ( 0, 1 );
		} elsif( $end_flag && $counter == 0 ) {
			my $end_num = $num;
			$total += ( $end_num - $start_num );
			( $start_flag, $end_flag, $start_num ) = ( 1, 0, 0 );
		}
	}
	return $total+1;
}

sub by_num {
	my( $start, $end ) = ( $a, $b );
	$start =~ s/\D+//g;
	$end =~ s/\D+//g;
	return $start <=> $end;
}


sub parse_file_list {
	my( $file_list ) = @_; 
	open( FH, "<$file_list" ) or die "Error in opening the file, $file_list, $!\n";
	while( my $file = <FH> ) {
		chomp $file;
		process_file( $file );
	}
	close FH or die "Error closing the file, $file_list, $!\n";
}


sub process_file {
	my( $file ) = @_;
	open( FH2, "<$file" ) or die "Error in opening the file, $file, $!\n";
        while( my $line = <FH2> ) {
                chomp $line;
		my @array = split( "\t", $line);
		my ( $start, $end ) = $array[ 8 ] <= $array[ 9 ] ? ( "start:".$array[ 8 ], "end:".$array[ 9 ] ) : ( "end:".$array[ 8 ], "start:".$array[ 9 ] );
		push @{ $$info{ $array[ 0 ] } { $array[ 1 ] } { 'array' } }, ( $start, $end );
        }
        close FH2 or die "Error closing the file, $file, $!\n";

}

sub get_options {
	my $hash = {};
	GetOptions( $hash, 'input_file_list|l=s', 'output_file|o=s', 'help|h' ) || pod2usage( );
	if( $$hash{ 'help' } ) {
		pod2usage( { -verbose => 1, -exitval => 1, -output => \*STDOUT } );
	}
	return $hash;
}
