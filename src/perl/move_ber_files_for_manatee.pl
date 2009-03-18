#!/usr/bin/perl

=head1 NAME

move_ber_files_for_manatee.pl - will copy a set of files in a given directory to a remote server

=head1 SYNOPSIS

USAGE: move_ber_files_for_manatee.pl 
      --ber_output_directory=/usr/local/projects/aengine/output_repository/ber/1234_default
      --remote_directory=/remove/server/destination/dir
      --local_directory=/local/path/to/replace
      --server=khan
    [ --log 
      --help
    ]

=head1 OPTIONS

B<--ber_output_directory,-b>
    Directory where output files live for ber run.  The directory should also
    contain a list file of all the ber output files in the directory.

B<--remote_directory,-d>
    Path on remote server for destination of copied files

B<--local_directory,-l>
    Path on local machine to remove when copying to remote server.  Example:
    If a file is located in 

    /usr/local/projects/p1/output_repository/ber/1234_default/i1/g1/file.raw

    and the --remote_directory option is /export/www/data
    and the --local_directory option is /usr/local/projects

    The file will be copied to 
    /export/www/data/p1/output_repository/ber/1234_default/i1/g1/file.raw
    on the remote server

B<--server,-s>
    The name of the server to copy the files to

B<--log,-o> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

Takes the content of the directory passed in --ber_output_directory and copies that to --server.
It will remove the part of the absolute path specified in --local_directory and replace with
value of --remote_directory.

=head1  CONTACT

    Kevin Galens
    kgalens@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

my $ber_dir;
my $server;
my $user;
my $pass;
my $local_dir;
my $remote_dir;

my %options = ();
my $results = GetOptions (\%options, 
                          'ber_output_directory|b=s',
                          'server|s=s',
                          'local_directory|l=s',
                          'remote_directory|d=s',
                          'log|o=s',
                          'help|h') || pod2usage();

&check_options(\%options);

my $new_dir = $ber_dir;
$new_dir =~ s/$local_dir//;
$new_dir = $remote_dir."/".$new_dir;

my $mkdir = "ssh $server 'mkdir -p $new_dir'";
my $cmd = "scp -r $ber_dir/* $server:$new_dir";

eval {
    print "$mkdir\n";
    system( $mkdir );
};
if( $@ ) {
    die("Could not create remote directory $new_dir on $server [$mkdir]");
}

eval {
    print "$cmd\n";
    system( $cmd );
};
if( $@ ) {
    die("Could not scp directory ($ber_dir) to remove server ($server) [$cmd]");
}


sub check_options {
    my ($opts) = @_;

    if( $opts->{'help'} ) {
        pod2usage();
        exit(0);
    }

    my @reqs = qw(ber_output_directory server local_directory remote_directory);
    foreach my $req ( @reqs ) {
        die("Option $req is required") unless( exists( $opts->{$req} ) );
    }

    $ber_dir = $opts->{'ber_output_directory'};
    $server = $opts->{'server'};
    $local_dir = $opts->{'local_directory'};
    $remote_dir = $opts->{'remote_directory'};

    if( !defined( $ber_dir ) || !defined( $server ) || !defined( $local_dir ) ||
        !defined( $remote_dir ) ) {
        die("Something's not defined");
    }

    #make sure remote_dir is an absolute path
    if( $remote_dir !~ m|^/| ) {
        die("option --directory should be absolute path (from / )");
    }
}
