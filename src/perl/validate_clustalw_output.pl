#!/usr/local/bin/perl
# $Id$

=head1 NAME

validate_clustalw_output.pl - Verifies the clustalw output is correct

=head1 SYNOPSIS

USAGE:  validate_clustalw_output.pl [-d debug_level] [-h] -i indir [-l log4perl] [-m]

=head1 OPTIONS

=over 8

=item B<--debug_level,-d>

    Optional: Coati::Logger log4perl logging level.  Default is 0

=item B<--help,-h>

    Print this help

=item B<--indir,-i>

    Jaccard cluster subflow output directory which contains the .clw, .dnd and jkcluster.out files

=item B<--log4perl,-l>

    Optional - log4perl log file.  Default is /tmp/validate_clustalw_output.pl.log

=back

=head1 DESCRIPTION

    validate_clustalw_output.pl - Verifies the clustalw output is correct

    Assumptions:
    1. User has appropriate permissions (to execute script, access chado database, write to output directory).
    2. All software has been properly installed, all required libraries are accessible.

    Sample usage:
    ./validate_clustalw_output.pl -i /usr/local/scratch/nema/jaccard/14225


=cut

use Mail::Mailer;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Data::Dumper;
use Log::Log4perl qw(get_logger);
use Workflow::Logger;


$|=1;

#-------------------------------------------------------------
# Parse command line options
#-------------------------------------------------------------

my ($debug_level, $help, $log4perl, $man, $indir);


my $results = GetOptions (
			  'log4perl|l=s'        => \$log4perl,
			  'debug_level|d=s'     => \$debug_level, 
			  'help|h'              => \$help,
			  'man|m'               => \$man,
			  'indir|i=s'           => \$indir
			  );


&pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if ($man);
&pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT}) if ($help);



&print_usage if(!$indir);

$debug_level = 5;

#
# initialize the logger
#
$log4perl = "/tmp/validate_clustalw_output.pl.log" if (!defined($log4perl));
my $mylogger = new Workflow::Logger('LOG_FILE'=>$log4perl,
				 'LOG_LEVEL'=>$debug_level);

my $logger = Workflow::Logger::get_logger(__PACKAGE__);




#
# Prep exec strings
#
my $clusterexecstring = "grep COG $indir/jkcluster.out | grep -v 'size 1,' | wc -l";
$logger->debug("cluster count execution string '$clusterexecstring'") if $logger->is_debug;

my $fastaexecstring = "find $indir -name \"*.fasta\" -type f |wc -l ";
$logger->debug("fasta count execution string '$fastaexecstring'") if $logger->is_debug;

my $dndexecstring = "find $indir -name \"*.dnd\" -type f |wc -l ";
$logger->debug("dnd count execution string '$dndexecstring'") if $logger->is_debug;


my $clustercount = qx{$clusterexecstring};
my $fastacount   = qx{$fastaexecstring};
my $dndcount     = qx{$dndexecstring};


$clustercount =~ s/\s+//g;
$fastacount =~ s/\s+//g;
$dndcount =~ s/\s+//g;


$logger->debug("clustercount '$clustercount' fastacount '$fastacount' dndcount '$dndcount'") if $logger->is_debug;

my $fatalctr=0;

if ($clustercount != $fastacount){
    $fatalctr++;
}
if ($clustercount != $dndcount){
    $fatalctr++;
} 
if ($dndcount != $fastacount){
    $fatalctr++;
}

if ($fatalctr > 0 ){

    $logger->logdie("Please verify whether clustalw execution failed abruptly.\n".
		    "clustercount '$clustercount'\n".
		    "dndcount '$dndcount'\n".
		    "fastacount '$fastacount'\n".
		    "clusterexecstring '$clusterexecstring'\n".
		    "fastaexecstring '$fastaexecstring'\n".
		    "dndexecstring '$dndexecstring'");
}


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
#                                  END OF MAIN  -- SUBROUTINES FOLLOW
#
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------
# print_usage()
#
#------------------------------------------------------
sub print_usage {

    print STDERR "SAMPLE USAGE:  $0 [-d debug_level] [-f filelist] [-h] -i indir [-l log4perl] [-m]\n".
    "  -d|--debug_level         = Optional - Coati::Logger log4perl logging level.  Default is 0\n".
    "  -h|--help                = Optional - Display pod2usage help screen\n".
    "  -i|--indir               = Jaccard subflow directory containing .fasta, .clw, .dnd and jkcluster.out files\n".
    "  -l|--log4perl            = Optional - Log4perl log file (default: /tmp/validate_clustalw_output.pl.log)\n".
    "  -m|--man                 = Optional - Display pod2usage pages for this utility\n";
    exit 1;

}
