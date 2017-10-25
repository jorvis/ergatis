#!/usr/bin/env perl

=head1 NAME

genbank2gff3.pl - create a gff3 file from a gbk file

=head1 SYNOPSIS

 USAGE: genbank2gff3.pl 
	--gb_file=<genbank file> 
	--gff3_file=<gff3 file> 
	--feat_type=<feature type in genbank file>

=head1  CONTACT

    Umar Farooq, Todd Creasy
    ufarooq@som.umaryland.edu
    tcreasy@som.umaryland.edu

=cut

use strict;
use warnings;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use NGS::GenBank;
use NGS::GFF3;

my %options = ();
my $results = GetOptions (\%options, 
						  'gb_file=s',
						  'gff3_file=s',
						  'feat_type=s',
                          'help|h') || pod2usage();

if( $options{'help'} ) {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

# check for required params
check_parameters(\%options);

########################################

# genbank
my $gb = NGS::GenBank->new( file => $options{gb_file} );

# change genbank seq id value
$gb->format_seq_id( format => "display_id:seq_version" );

# gff3 file
my $gff3 = NGS::GFF3->new( file => ">$options{gff3_file}" );

# write gff3 file with selected features
if( !defined($options{feat_type}) ) {
	$gff3->write( $gb->get_all_features );

} else {
	$gff3->write( $gb->get_features_from_type( feat_type => $options{feat_type} ) );
}


########################################

sub check_parameters {
    my ($options) = @_;
    
    # make sure required arguments were passed
    my @required = qw( gb_file gff3_file );

    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }

	# default feat_type is CDS
	if( !defined( $options{feat_type} ) ) {
		$options{feat_type} = "CDS";
	}
}

