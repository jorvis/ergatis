#!/usr/local/bin/perl
# $Id$

=head1 NAME

validate_clustalw_output.pl - Verifies the clustalw output is correct

=head1 SYNOPSIS

USAGE:  validate_clustalw_output.pl -c clwfile [-d debug_level] [-f fastafile] [-h] [-l log4perl] [-m]

=head1 OPTIONS

=over 8

=item B<--clwfiler,-c>

    Clustalw .clw file

=item B<--debug_level,-d>

    Optional: Coati::Logger log4perl logging level.  Default is 0

=item B<--fastafile,-f>

    Input fasta file

=item B<--help,-h>

    Print this help

=item B<--log4perl,-l>

    Optional - log4perl log file.  Default is /tmp/validate_clustalw_output.pl.log

=back

=head1 DESCRIPTION

    validate_clustalw_output.pl - Verifies the clustalw output is correct

    Assumptions:
    1. User has appropriate permissions (to execute script, access chado database, write to output directory).
    2. All software has been properly installed, all required libraries are accessible.

    Sample usage:
    ./validate_clustalw_output.pl -c /usr/local/scratch/nema/jaccard/14225/jaccard_jaccard_998.10538.clw


=cut

use Mail::Mailer;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Data::Dumper;
use Log::Log4perl qw(get_logger);
BEGIN {
use Workflow::Logger;
}


$|=1;

#-------------------------------------------------------------
# Parse command line options
#-------------------------------------------------------------

my ($debug_level, $help, $log4perl, $man, $fastafile, $clwfile);


my $results = GetOptions (
			  'log4perl|l=s'        => \$log4perl,
			  'debug_level|d=s'     => \$debug_level, 
			  'help|h'              => \$help,
			  'man|m'               => \$man,
			  'fastafile|f=s'       => \$fastafile,
			  'clwfile|c=s'         => \$clwfile
			  );


&pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if ($man);
&pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT}) if ($help);



&print_usage if(!$clwfile);

$debug_level = 5;

#
# initialize the logger
#
$log4perl = "/tmp/validate_clustalw_output.pl.log" if (!defined($log4perl));
my $mylogger = new Workflow::Logger('LOG_FILE'=>$log4perl,
				 'LOG_LEVEL'=>$debug_level);

my $logger = Workflow::Logger::get_logger(__PACKAGE__);



$logger->debug("Checking clustalw file '$clwfile'") if $logger->is_debug;

#
# Prep exec strings
#
if (!-e $clwfile){
    $logger->logdie("clustalw file '$clwfile' does not exist");
}
if (-z $clwfile){
    $logger->logdie("clustalw file '$clwfile' has zero size");
}

if (defined($fastafile)){

    $logger->debug("Checking fasta file '$fastafile'") if $logger->is_debug;

    my $nameexec = "grep \"Name:\" $clwfile | wc -l";
    my $namecount = qx{$nameexec};
    $namecount =~ s/\s+//g;

    my $headexec = "grep \">\" $fastafile | wc -l";
    my $headcount = qx{$headexec};
    $headcount =~ s/\s+//g;

    if ($namecount != $headcount){
	$logger->logdie("Sequence counts were off.\n".
			"namecount '$namecount'\n".
			"headcount '$headcount'\n".
			"nameexec '$nameexec'\n".
			"headexec '$headexec'\n");
    }



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

    print STDERR "SAMPLE USAGE:  $0 -c clwfile [-d debug_level] [-f fastafile] [-h] [-l log4perl] [-m]\n".
    "  -c|--clwfile             = Clustalw output file\n".
    "  -d|--debug_level         = Optional - Coati::Logger log4perl logging level.  Default is 0\n".
    "  -h|--help                = Optional - Display pod2usage help screen\n".
    "  -f|--fastafile           = Optional - fasta file to check header counts\n".
    "  -l|--log4perl            = Optional - Log4perl log file (default: /tmp/validate_clustalw_output.pl.log)\n".
    "  -m|--man                 = Optional - Display pod2usage pages for this utility\n";
    exit 1;

}
