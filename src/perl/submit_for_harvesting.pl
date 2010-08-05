#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

submit_for_harvesting.pl - A small script to submit a specific output directory for harvesting by a master node ina CloVR cluster.

=head1 SYNOPSIS

USAGE: ./submit_for_harvesting.pl --harvest_dir=/path/to/harvest/directory

=head1 OPTIONS

B<--harvest_dir, -h>
    The directory that should be harvested by the master node.

B<--log, -l>
    Log file

B<--debug, -d>
    Debug level

=head1 DESCRIPTION

=head1 INPUT

A string containing the output directory that should be harvested by the master node.

=head1 CONTACT

    Cesar Arze  
    carze@som.umaryland.edu

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Ergatis::Logger;

my %options = ();
my $results = GetOptions (\%options,
                          'harvest_dir|h=s',
                          'log|l=s',
                          'debug|d=s' ) || die ("Unprocessable option");

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger( 'LOG_FILE'    =>  $logfile,
                                  'LOG_LEVEL'   =>  $options{'debug'} );
$logger = Ergatis::Logger::get_logger();

## Check if harvesting directory exists.
my $harvest_dir = $options{'harvest_dir'};
if (-d $harvest_dir && -e $harvest_dir && -r $harvest_dir) {

} else {
    $logger->logdie("Directory supplied for harvesting is not readable");
} 

my $hostname = `hostname -f`;
chomp($hostname);

my $cmd = "/opt/sge/bin/lx24-amd64/qsub -o /mnt/scratch -e /mnt/scratch -S /bin/sh -b n -sync y -q harvesting.q /opt/vappio-scripts/harvesting.sh $hostname $harvest_dir";
system($cmd);

###################################################################

sub get_overall_status {
    my ($id_refs, $request) = @_;       
    my $success;         
                     
    foreach my $id (@$id_refs) {
        if ($request->get_status($id) eq "FINISHED") {
            $success = 1;
        } else {
            $success = 0;
            last;
        }
    }

    return $success;
}
