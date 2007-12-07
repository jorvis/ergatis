#!/usr/bin/perl
=head1  NAME 

generate_asmbl_list.pl - Default output is a workflow iterator that
can be used to iterator over a set of asmbl_ids

=head1 SYNOPSIS

USAGE:  generate_asmbl_list

=head1 OPTIONS

=item *

B<--bcp_extension> Optional - BCP tab-delimited file extension

=item *

B<--bcp_dir> Optional - Directory containing BCP files with bcp_extension

=item *

B<--output_file> Name of output file whose file size will be verified

=item *

B<--debug,-d> Debug level.  Use a large number to turn on verbose debugging. 

=item *

B<--logfile> Log4perl log file

=item *

B<--help,-h> This help message

=head1   DESCRIPTION

=cut

#use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);


BEGIN {
    require '/usr/local/devel/ANNOTATION/cas/lib/site_perl/5.8.5/Workflow/Logger.pm';
    import Workflow::Logger;
}
use File::Basename;

umask(0000);

my %options = ();


my $results = GetOptions (\%options, 
                          'bcp_extension=s', 
			  'bcp_dir=s',
			  'logfile=s',
			  'output_file=s',
                          'debug=s', 
                          'help|h' ) || pod2usage();

my $logfile = $options{'logfile'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Workflow::Logger::get_logger();


my $output_file;
if (!defined($options{'output_file'})){
    $logger->logdie("output_file was not defined");
}
else {
    $output_file = $options{'output_file'};
}

if (!-e $output_file){
    $logger->logdie("output_file '$output_file' does not exist");
}
if (!-f $output_file){
    $logger->logdie("output_file '$output_file' is not a file");
}
if (!-r $output_file){
    $logger->logdie("output_file '$output_file' does not have read permissions");
}
if (-z $output_file){

    if (defined($options{'bcp_dir'})){
	
	my $bcp_dir = $options{'bcp_dir'};

	if (-e $bcp_dir){
	    if (-d $bcp_dir){
		if (-r $bcp_dir){
		    if (defined($options{'bcp_extension'})){

			my $bcp_extension = $options{'bcp_extension'};

			my $execstring = "wc -l " . $bcp_dir . "/*" . $bcp_extension;

			my @fileRecordCounts;

			eval { @fileRecordCounts = qx($execstring); };
			if ($@){
			    $logger->logdie("Some error was encountered while executing '$execstring': $!");
			}
			$logger->fatal("@fileRecordCounts");
		    }
		    else {
			$logger->logdie("bcp_dir was specified, but bcp_extension was not defined");
		    }
		}
		else {
		    $logger->logdie("bcp_dir '$bcp_dir' does not have read permissions");
		}
	    }
	    else {
		$logger->logdie("bcp_dir '$bcp_dir' is not a directory");
	    }
	}
	else {
	    $logger->logdie("bcp_dir '$bcp_dir' does not exist");
	}
    }

    $logger->logdie("output_file '$output_file' does not have any content");
}
else {
    print "output_file '$output_file' does have content\n";
}

print "$0 program execution completed\n";
print "Log file is '$logfile'\n";
exit(0);


