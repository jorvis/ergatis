#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib ".";

=head1  NAME 

store_chado_version.pl - 

=head1 SYNOPSIS

USAGE:  store_chado_version

=head1 OPTIONS

=over 8

=item B<--bin_dir> 

The chado version will be derived from the BIN_DIR path

=item B<--database> 

Chado database
 
=item B<--log4perl,-l> Log file

Optional - log file (default is /tmp/store_chado_version.pl.log

=item B<--debug,-d> Debug level

Optional

=item B<--workflow> workflow_repository

Optional - directory to store the update.sql

=item B<--username> Username

username must have permissions to update common..genomes

=item B<--password> Password

password

=item B<--help,-h> This help message

=head1   DESCRIPTION

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

use File::Basename;
BEGIN {
use Ergatis::Logger;
}

umask(0000);

my %options = ();

my ($bindir, $database, $log4perl, $help, $debug, $workflow, $username, $password);

my $results = GetOptions (\%options, 
                          'bin_dir=s'     => \$bindir, 
                          'database=s'    => \$database,
                          'log4perl|l=s'  => \$log4perl, 
			  'debug|d=s'     => \$debug,
			  'workflow=s'    => \$workflow,
			  'username=s'    => \$username,
			  'password=s'    => \$password,
                          'help|h'        => \$help ) || pod2usage();



my $log4perl = $options{'log4perl'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$log4perl,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Ergatis::Logger::get_logger();

if (!defined($bindir)){
    die "bin_dir was not defined";
}
if (!defined($database)){
    die "database was not defined";
}
if (!defined($username)){
    die "username was not defined";
}
if (!defined($password)){
    die "password was not defined";
}
if (!defined($workflow)){
    $workflow = "/tmp";
}

## Extract the version e.g. chado-v1r7b1
my $dirname = dirname($bindir);
my $version = basename($dirname);

my $file = $workflow . "/update.sql";
open (OUTFILE, ">$file") or $logger->logdie("Could not open file '$file' for output: $!");

## Keep a record
print OUTFILE "update common..genomes set type ='$version' where db = '$database'\ngo\n";

my $execstring = "sqsh -U $username -P $password -D $database -i $file";

## Update common..genomes
eval { qx{$execstring} };

if ($@){
    $logger->logdie("Error detected: $!\nexecstring was '$execstring'");
}

exit(0);


