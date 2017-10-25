#!/usr/bin/env perl

eval 'exec /usr/bin/env perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
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

=item B<--password_file, -f>

Optional - Read password from a file rather than use the --password option.  Useful if using --password results in the password being publicly visible.  If both --password_file and --password are provided, the entry from --password_file will take priority
 
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

my ($bindir, $database, $log4perl, $help, $debug, $workflow, $username, $password,
$pass_file, $database_type);

my $results = GetOptions (\%options, 
                          'bin_dir=s'     => \$bindir, 
                          'database=s'    => \$database,
                          'log4perl|l=s'  => \$log4perl, 
			  'debug|d=s'     => \$debug,
			  'workflow=s'    => \$workflow,
			  'username=s'    => \$username,
			  'password=s'    => \$password,
			  'password_file|f=s'	=> \$pass_file,
			  'database_type=s' => \$database_type,
                          'help|h'        => \$help ) || pod2usage();



my $log4perl = $options{'log4perl'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$log4perl,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Ergatis::Logger::get_logger();

if ($database_type eq 'postgresql'){
    $logger->warn("$0 will not store revision information in PostgreSQL databases");
    exit(0);
}

if (!defined($bindir)){
    die "bin_dir was not defined";
}
if (!defined($database)){
    die "database was not defined";
}
if (!defined($username)){
    die "username was not defined";
}
if (!defined($password) && !defined($pass_file)){
    die "Neither password or a password file were defined.  Please use one or the other";
}
if (!defined($workflow)){
    $workflow = "/tmp";
}

#Assign password to be read from the file if it exists.
if (defined ($pass_file)  && -s $pass_file ) {
	open PASS, $pass_file or die ("Cannot open password file $pass_file : $!\n");
	print STDERR ("Password from file will take priority over --password option\n") if (defined $password);
	$password= <PASS>;
	chomp $password;
	close PASS;
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


