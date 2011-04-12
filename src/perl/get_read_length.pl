#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;

my $options = get_options ( );
open( OFH, ">$$options{'output_file'}" ) or die "Error in opening the file, $$options{'output_file'}, $!\n";
if( $$options{'input_file'} ) {
	get_read_length( $$options{'input_file'} );
}

if( $$options{'input_file_list'} ) {
	open( FH, "<$$options{'input_file_list'}" ) or die "Error in opening the file, $$options{'input_file_list'}, $!\n";
	while( my $file = <FH> ) {
		chomp $file;
		unless( -e $file ) {
			print STDERR "$file doen't exists\n";
			next;
		}
		get_read_length( $file );
	}	
}

close OFH or die "Error closing the outputfile, $$options{'output_file'}, $!\n";
exit $?;


sub get_read_length {
	my ($file) = @_;
	my $seq_in = Bio::SeqIO -> new ( -file => $file, -format => $$options{'format'} );
	while ( my $seq = $seq_in -> next_seq( ) ) {
		my $desc;
		if( $seq -> desc =~ /.+Source:\s*(.+)$/ ) {
			$desc = $1;
		}
		print OFH $seq -> id, "\t", $seq -> length, "\t", $desc, "\n";
	}
}


sub get_options {
	my %options;
	GetOptions ( \%options, 'format|f=s', 'input_file_list|l:s', 'input_file|i:s', 'output_file|o=s', 'help|h' ) || pod2usage ( ) ;
	if ( $options{'help'} ) {
		pod2usage ( { -exitval => 0, -verbose => 1, -output => \*STDOUT } );
	}
	return \%options;
}
