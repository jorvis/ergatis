#!/usr/local/bin/perl

=head1 NAME

concatenate_files.pl - will concatenate files.

=head1 SYNOPSIS

 USAGE: concatenate_files.pl
    --input_files=/path/to/file.1.tab,/path/to/file2.tab
    --input_lists=/path/to/list1,/path/to/list2
    --output=/path/to/new.file.tab
    [ --log=/path/to/file.log
      --debug=4
      --help
    ]

=head1 OPTIONS

B<--input_files,-f>
    A comma separated list of files.  Will cat each file.

B<--input_lists,-l>
    A comma separated list of lists.  Will cat each files in the list.

B<--output,-o>
    Output file

B<--log,-l>
    Logfile.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

 DESCRIPTION
 
=head1  INPUT
    Can be any type of input.  Will simply cat the files and append to same output

=head1 OUTPUT
    All the input files specified in one file

=head1  CONTACT

    Kevin Galens
    kgalens@som.umaryland.edu

=cut

use strict;
use warnings;
use File::OpenFile qw( open_file );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Data::Dumper;

my ( $input_files, $input_lists, $output );
my $results = GetOptions ("input_files=s" => \$input_files,
                          "input_lists=s" => \$input_lists,
                          "output=s" => \$output,
                          );

my $oh = open_file( $output, 'out' );

my @input_files = split( ",", $input_files ) if( $input_files );
my @input_lists = split( ",", $input_lists ) if( $input_lists );

foreach my $list ( @input_lists ) {
    my $lh = open_file( $list, 'in' );
    chomp( my @tmp = <$lh> );
    push( @input_files, @tmp );
    close($lh);
}

foreach my $input_file ( @input_files ) {
    my $in = open_file( $input_file, 'in' );
    print $oh (<$in>);
    close($in);
}

close($oh);
