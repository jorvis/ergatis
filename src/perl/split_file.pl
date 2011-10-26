#!/usr/bin/perl

=head1 NAME

split_file.pl - will split a file (or list of files) into chunks of the specified size or into the specified number of chunks.

=head1 SYNOPSIS

 USAGE: split_file.pl
    --input_files=/path/to/list.txt
    --input_file_list=/path/to/list_of_lists.txt
    --output_directory=/path/to/output_dir/
    --lines_per_file=100 || --num_files=100

=head1 OPTIONS

B<--input_file,-f>
    The file to split. Could also be a comma separated list of files to redistribute lines from.

B<--input_file_list,-l>
    A list of files to redistribute into files of the selected size

B<--output_directory,-o>
    The directory to put the split files.

B<--lines_per_file,-l>
    The number of lines to print per file. Cannot be used with num files.

B<--num_files,-f>
    The number of total files to make. Cannot be used with lines_per_file.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

=head1  AUTHOR

    David Riley
    driley@som.umaryland.edu

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use File::Basename;

my ( $input_files, $input_lists, $output_dir, $num_lines, $num_files, $prefix, $help);
my $results = GetOptions ("input_files:s" => \$input_files,
                          "input_file_list:s" => \$input_lists,
                          "lines_per_file:s" => \$num_lines,
                          "num_files:s" => \$num_files,
                          "prefix:s" => \$prefix,
                          "output_directory=s" => \$output_dir,
                          "help|h" => \$help
                          ) || pod2usage();

if( $help ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

&check_parameters();

my @input_files = split( ",", $input_files ) if( $input_files );
my @input_lists = split( ",", $input_lists ) if( $input_lists );

foreach my $list ( @input_lists ) {
    open IN1, "<$list" or die "Unable to open $list\n";
    while(<IN1>) {
        chomp;
        push( @input_files, $_ );
    }
    close(IN1);
}

# Need to count the lines
my $total_lines = 0;
foreach my $f (@input_files) {
    my $c = `wc -l $f`;
    my($count,$fname) = split(/\s/,$c);
    $total_lines += $count;
    print STDERR "Found $total_lines lines\n";
}
if($num_files) {
    my $lpf = $total_lines/$num_files;
    if(!$num_lines) {
        $num_lines = ($lpf == int($lpf)) ? $lpf : int($lpf + 1);
    }
}
else {
    my $nf = $total_lines/$num_lines;
    $num_files = ($nf == int($nf)) ? $nf : int($nf + 1);
}

print STDERR "Going to print $num_lines to $num_files files\n";

my $num_digits = length($num_files);
my $line_count =0;
my $file_num = 0;
my $fname = $prefix."_".(sprintf("%0$num_digits\d",$file_num)).".splitlist";
open my $fh, ">$output_dir/$fname" or die "Unable to open $fname for writing\n";

foreach my $input_file ( @input_files ) {
    open IN, "<$input_file" or die "Unable to open $input_file\n";
    while(<IN>) {
        &process_line($_);
    }
}

sub process_line {
    my $line = shift;
    
    if($line_count >= $num_lines) {
        $line_count = 0;
        $file_num++;
        close $fh;
        $fname = $prefix."_".(sprintf("%0$num_digits\d",$file_num)).".splitlist";
        open $fh, ">$output_dir/$fname" or die "Unable to open $fname for writing\n";
    }
    $line_count++;
    print $fh $line;
}

sub check_parameters {

    if($num_lines && $num_files) {
        die "Can't specify both num_files and lines_per_file\n";
    }
    if(!($num_lines || $num_files)) {
        die "Need to specify --num_files or --lines_per_file\n";
    }
    if(!$prefix) {
        $prefix = 'split_file';
    }

}
