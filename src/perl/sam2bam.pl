#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

=head1 NAME

sam2bam.pl - run the samtools convert sam files to sorted bam (by position and name) and indexed bam files (by position)

=head1 SYNOPSIS

  USAGE: sam2bam.pl [
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
convert_sam_to_bam(\%options);

# sort the bam file by position
my $sorted_bam_file = sort_bam_by_position(\%options);

# index the bam file
index_bam($sorted_bam_file);

# sort the bam file by name
my $sorted_bam_file_by_name = sort_bam_by_name(\%options, $sorted_bam_file);

exit(0);

#######################

sub convert_sam_to_bam {
	my ($options) = @_;

	# samtools view
	my $cmd = qq{$options->{'samtools-exec'} view $options->{view_options} -o $options->{output_file} $options->{input_file}};
	
    system( $cmd );
}

sub sort_bam_by_position {
	my ($options) = @_;

	# just setting up the prefix file (file.sorted instead of file.bam)
	my $sorted_file = $options->{output_file};
	$sorted_file =~ s/\.bam/\.position_sorted/;

	# this is what the final output file will be from the sort command
	my $sorted_bam_file = $sorted_file . ".bam";
	
	# samtools sort by position
	my $cmd = qq{$options->{'samtools-exec'} sort $options->{sort_options} $options->{output_file} $sorted_file};

    system( $cmd );
	
	return $sorted_bam_file;	
}

sub index_bam {
	my ($sorted_bam_file) = @_;

	# samtools index
	my $cmd = qq{samtools index $sorted_bam_file};

	system( $cmd );
}

sub sort_bam_by_name {
	my ($options, $bam_file) = @_;
	
	my $base_name = $bam_file;
	$base_name =~ s/\.position_sorted\.bam/\.name_sorted/;

	my $cmd = qq{$options->{'samtools-exec'} sort -n $bam_file $base_name};

    system( $cmd );
}

sub check_parameters {
    my ($options) = @_;
    
    # make sure required arguments were passed
    my @required = qw( samtools-exec input_file output_directory);
	
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }

	my $sam_file = $options->{input_file};
	my ($base_name) = ($sam_file =~ /.*\/(.*)\.sam/);

	# bam output file
	$options->{output_file} = $options->{output_directory} . "/" . $base_name . ".bam";
}

