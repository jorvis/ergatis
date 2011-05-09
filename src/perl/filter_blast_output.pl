#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my $options = get_options( );
my $files = get_files( $$options{'blast_file_list'} );
my ( $reads_length, $reads_source ) = get_reads_info( $$options{ 'read_length_file' } );
## we need a hash to keep track of sub start and sub end for each query
## we need to hold on to these values through out all files
my $filter = {};
my $data = {};
foreach my $file ( @$files ) {
	parse_file( $file, $reads_length, $reads_source, $$options{ 'identity_percentage' }, $$options{ 'read_length_coverage' } );
}
print_output( $data, $$options{ 'output_file' } );
exit $?;

sub print_output {
	my( $data, $file ) = @_;
	open( OFH, ">$file" ) or die "Error in writing to the file, $file, $!\n";
	while( my( $key, $value ) = each %$data ) {
		print OFH $key, "\t", join( "\t", @$value ), "\n";
	}
	close OFH;
}

sub parse_file {
	my ( $file, $reads_length, $reads_source, $identity_cutoff, $coverage_cutoff ) = @_;
	open( FH, "<$file" ) or die "Error in reading the file, $file, $!\n";
	while( my $line = <FH> ) {
		chomp $line;
		my @values = split( "\t", $line );
		my ( $query_id, $sub_id, $identity, $coverage, $sub_start, $sub_end ) = ( @values[0 .. 3], @values[8 .. 9] );
		if ( $identity >= $identity_cutoff && ( $coverage >= $coverage_cutoff * $$reads_length{ $sub_id } / 100 ) ) {
			unless( exists $$filter{ $query_id } ) {
				$$filter{ $query_id }  = $sub_end - $sub_start;
				$$data{$query_id} = [ $sub_id, $identity, $coverage, $$reads_length{ $sub_id }, $sub_start, $sub_end, $$reads_source{ $sub_id } ];
			} elsif( is_valid( $sub_start, $sub_end, $$filter{ $query_id } ) ) {
					$$filter{ $query_id } = $sub_end - $sub_start;
					$$data{$query_id} = [ $sub_id, $identity, $coverage, $$reads_length{ $sub_id }, $sub_start, $sub_end, $$reads_source{ $sub_id } ];
			}
		}
	}
}

sub is_valid {
	my ( $sub_start, $sub_end, $length ) = @_;
	return (( $sub_end - $sub_start ) > $length ) ? 1 : 0;
}


sub get_options {
	my ( %options );
	GetOptions( \%options, 'read_length_file|r=s', 'blast_file_list|b=s', 'output_file|o=s',
			'identity_percentage|i=i', 'read_length_coverage|c=i', 'help|h' ) || pod2usage( );
	if( $options{'help'} ) {
		pod2usage( {-exitval => 0, -verbose => 1, -output => \*STDOUT } );
	}
	return \%options;
}

sub get_files {
	my ( $file_list ) = @_;
	my @files;
	open( FH, "<$file_list" ) or die "Error in opening the file, $file_list, $!\n";
	while( my $file = <FH> ) {
		chomp $file;
		unless( -f $file ) {
			print STDERR "Cannot find find: $file\n";
		} else {
			push @files, $file;
		}		
	}
	close FH;
	return \@files;
}

sub get_reads_info {
	my ($file) = @_;
	my $read_length_info = {};
	my $read_source_info = {};
	open( FH, "<$file" ) or die "Error in reading the read length file, $file, $!\n";
	while( my $line = <FH> ) {
		chomp $line;
		my ($id, $length, $source) = split( "\t", $line );
		$$read_length_info{ $id } = $length;
		$$read_source_info{ $id } = $source;
	}
	close FH;
	return ( $read_length_info, $read_source_info );
}
