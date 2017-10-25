#!/usr/bin/env perl

eval 'exec /usr/bin/env perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/env perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

=head1 NAME

bam2sam.pl - convert bam file back to sam file

=head1 SYNOPSIS

  USAGE: bam2sam.pl [
	--samtools-exec=path to samtools binary
	--help
	]
=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
						  'samtools-exec=s',
						  'view_options=s',
						  'sort_options=s',
						  'input_file=s',
						  'output_directory=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
check_parameters(\%options);

# sam to bam
convert_bam_to_sam(\%options);

exit(0);

#######################

sub convert_bam_to_sam {
	my ($options) = @_;

	# samtools view
	my $cmd = qq{$options->{'samtools-exec'} view $options->{view_options} -o $options->{output_file} $options->{input_file}};
	
    system( $cmd );
}

sub check_parameters {
    my ($options) = @_;
    
    # make sure required arguments were passed
    my @required = qw( samtools-exec input_file );
	
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }

	my $bam_file = $options->{input_file};
	my ($base_name) = ($bam_file =~ /.*\/(.*)\.bam/);

	# bam output file
	$options->{output_file} = $options->{output_directory} . "/" . $base_name . ".sam";
}

