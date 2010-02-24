#!/usr/bin/perl

=head1 NAME

prepare_for_tbl2asn.pl - Copies files to a directory which tbl2asn expects

=head1 SYNOPSIS

USAGE: prepare_for_tbl2asn.pl
            --tbl_file=/path/to/file.tbl
            --fasta_list=/path/to/fasta.list
            --output_directory=/path/to/dir
            --help
          
=head1 OPTIONS

B<--tbl_file,-t>
    A tbl file

B<--fasta_list,-f>
    Fasta list of files one of which must be named to match the name of the tbl file
    (with a .fsa extension)

B<--output_directory,-o>
    Directory to copy files to

B<--help,-h>
    Print this message

=head1  DESCRIPTION

This script is used to convert the output from a glimmer3 search into BSML.

=head1  CONTACT

    Kevin Galens
    kevingalens@gmail.com

=cut


use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

my $tbl_file;
my $fsa_file;
my $outdir;

&check_options();

system("cp $tbl_file $outdir");
system("cp $fsa_file $outdir");

sub check_options {
    my %options = ();
    my $results = GetOptions (\%options, 
                              'tbl_file|t=s',
                              'fasta_list|f=s',
                              'output_directory|o=s',
                              'help|h') || &_pod;

    foreach my $req_opt( qw(tbl_file fasta_list output_directory) ) {
        die("Option $req_opt is required") unless( $options{$req_opt} );
    }

    $outdir = $options{'output_directory'};
    $tbl_file = $options{'tbl_file'};
    my $basename = $1 if( $tbl_file =~ m|/([^/]+)\.tbl| );
    die("Could not parse basename from tbl file: $tbl_file") unless( $basename );

    open(IN, "< $options{'fasta_list'}");
    while(<IN>) {
        chomp;
        if( /$basename\.fsa/ ) {
            $fsa_file = $_;
            last;
        }
    }
    close(IN);
    die("Could not find a fsa file which matched the basename: $basename") unless( $fsa_file );
}




sub _pod {   
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
