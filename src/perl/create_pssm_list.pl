#!/usr/local/bin/perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case bundling);
use IO::File;
use File::Basename;

my $dir = undef;
my $out = *STDOUT;

parse_options();
create_list();

sub print_usage
{
    my $progname = basename($0);
    die << "END";
usage: $progname --input_directory|-d <input_directory>
        [--output|-o <output_list>] [--help|-h]
END
}

sub parse_options
{
    my %opts = ();
    GetOptions(\%opts, "input_directory|d=s", "output|o=s", "help|h");
    print_usage() if $opts{help};
    if (exists $opts{input_directory}) {
        $dir = $opts{input_directory};
    }
    else {
        print STDERR "No input directory provided\n";
        print_usage();
    }
    $out = new IO::File($opts{output}, "w") or
        die "Error writing output list $opts{output}: $!"
        if $opts{output};
}

sub create_list
{
    print $out "\$;UID\$;\t\$;CHD\$;\t\$;CHK\$;\n";
    while (my $chd = <$dir/*.chd>) {
        my $uid = basename($chd);
        $uid =~ s/\.chd$//;
        my $chk = $chd;
        $chk =~ s/chd$/chk/;
        print $out "$uid\t$chd\t$chk\n" if -e $chk;
    }
}
