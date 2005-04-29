#!/usr/local/bin/perl

=head1 NAME

verify_logfiles.pl - greps all log4perl logfiles and checks for FATAL, ERROR, WARN

=head1 SYNOPSIS

USAGE:  verify_logfiles.pl [-d debug_level] [-f filelist] [-h] [-l log4perl] [-m] [-r repository] [-w workflow_id] -U username

=head1 OPTIONS

=over 8

=item B<--debug_level,-d>

    Optional: Coati::Logger log4perl logging level.  Default is 0

=item B<--filelist,-f>

    Optional - List of log4perl logfiles to be checked for FATAL, ERROR, WARN

=item B<--help,-h>

    Print this help

=item B<--log4perl,-l>

    Optional - log4perl log file.  Default is /tmp/verify_logfiles.pl.log

=item B<--man,-m>

    Display the pod2usage page for this utility

=item B<--workflow_id,-w>

    Optional - workflow_id

=item B<--repository,-r>

    Optional - $WORKFLOW_REPOSITORY$;

=item B<--username,-U>

    Username of person to be notified via email

=back

=head1 DESCRIPTION

    verify_logfiles.pl - greps all log4perl logfiles and checks for FATAL, ERROR, WARN

    Assumptions:
    1. User has appropriate permissions (to execute script, access chado database, write to output directory).
    2. All software has been properly installed, all required libraries are accessible.

    Sample usage:
    ./check_logfiles.pl -f rebuild.log,drop_indexes.log,load.log -l check_logfiles.pl.log -U sundaram


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

my ($debug_level, $help, $log4perl, $man, $filelist, $workflow_id, $username, $repository);


my $results = GetOptions (
			  'log4perl|l=s'     => \$log4perl,
			  'debug_level|d=s'  => \$debug_level, 
			  'help|h'           => \$help,
			  'man|m'            => \$man,
			  'filelist|f=s'     => \$filelist,
			  'workflow_id|w=s'  => \$workflow_id,
			  'username|U=s'     => \$username,
			  'repository|r=s'   => \$repository
			  );


&pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if ($man);
&pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT}) if ($help);



&print_usage if(!$username);

$debug_level = 5;

#
# initialize the logger
#
$log4perl = "/tmp/check_logfiles.pl.log" if (!defined($log4perl));
my $mylogger = new Workflow::Logger('LOG_FILE'=>$log4perl,
				 'LOG_LEVEL'=>$debug_level);

my $logger = Workflow::Logger::get_logger(__PACKAGE__);

$logger->info("Processing the following list of log4perl logfiles: '$filelist'");    



#
# if error is defined then the username better be too
#
if (!defined($username)){
    $logger->logdie("username was not defined");
}
else{
    $username .= "\@tigr.org";
    $logger->debug("username set to '$username'") if $logger->is_debug;
}

#
# Either repository or filelist should be defined
#
if ((!defined($filelist)) and (!defined($repository))){
    $logger->logdie("Either filelist or repository (or both) should be defined");
}


#
# if defined, repository should be a valid directory
#
if (defined($repository)){
    if (!-d $repository){
	$logger->logide("repository '$repository' is not a directory");
    }
}


my $list = &get_file_list($filelist, $repository);
my $filectr=0;

my $fatalmaster=0;
my $errormaster=0;
my $warnmaster=0;

my $filehash = {};

foreach my $file (sort @{$list}){

    $filectr++;

    $logger->info("Processing the following log4perl logfile '$file'");


    if (!-e $file){
	$logger->error("file '$file' does not exist");
	$filectr--;
	next;
    }
    if (!-r $file){
	$logger->error("file '$file' does not have read permissions");
	$filectr--;
	next;
    }
    if (-z $file){
	$logger->debug("file '$file' had zero size and therefore will not be processed");
	$filectr--;
	next;
    }


    open (INFILE, "<$file") or $logger->logdie("Could not open file '$file'");


    my $linectr=0;
    my $fatalctr=0;
    my $errorctr=0;
    my $warnctr=0;

    while ( my $line = <INFILE> ){
	
	chomp $line;

	$linectr++;

	if ($line =~ / FATAL /){
	    $fatalctr++;
	    my ($desc, $message) = split(/\d+>/,$line);
	    $logger->logdie("desc '$desc' and message '$message' were both not defined") if ((!defined($desc)) and (!defined($message)));
	    my ($prgline, $samplemsg) = split(/ - /,$message);
	    if ((defined($prgline)) and (defined($samplemsg))){
		$filehash->{$file}->{$prgline}->{'count'}++;
		$filehash->{$file}->{$prgline}->{'samplemsg'} = $samplemsg;
		$filehash->{$file}->{$prgline}->{'level'} = 'FATAL';
	    }
	    else{
		$logger->logdie("prgline '$prgline' and samplemsg '$samplemsg' were not both defined, could not parse file '$file' at line '$linectr'. Line was '$line'");
	    }
	}
	if ($line =~ / ERROR /){
	    $errorctr++;
	    my ($desc, $message) = split(/\d+>/,$line);
	    $logger->logdie("desc '$desc' and message '$message' were both not defined") if ((!defined($desc)) and (!defined($message)));
	    my ($prgline, $samplemsg) = split(/ - /,$message);
	    if ((defined($prgline)) and (defined($samplemsg))){
		$filehash->{$file}->{$prgline}->{'count'}++;
		$filehash->{$file}->{$prgline}->{'samplemsg'} = $samplemsg;
		$filehash->{$file}->{$prgline}->{'level'} = 'ERROR';
	    }
	    else{
		$logger->logdie("prgline '$prgline' and samplemsg '$samplemsg' were not both defined, could not parse file '$file' at line '$linectr'. Line was '$line'");
	    }
	}
	if ($line =~ / WARN /){
	    $warnctr++;
	    my ($desc, $message) = split(/\d+>/,$line);
	    $logger->logdie("desc '$desc' and message '$message' were both not defined") if ((!defined($desc)) and (!defined($message)));
	    my ($prgline, $samplemsg) = split(/ - /,$message);
	    if ((defined($prgline)) and (defined($samplemsg))){
		$filehash->{$file}->{$prgline}->{'count'}++;
		$filehash->{$file}->{$prgline}->{'samplemsg'} = $samplemsg;
		$filehash->{$file}->{$prgline}->{'level'} = 'WARN';
	    }
	    else{
		$logger->logdie("prgline '$prgline' and samplemsg '$samplemsg' were not both defined, could not parse file '$file' at line '$linectr'. Line was '$line'");
	    }
	}

    }


    close INFILE or $logger->logdie("Could not close file '$file'");

    
    if (($fatalctr < 1) and ($errorctr < 1) and ($warnctr < 1)){
	$logger->info("Found    '$fatalctr' fatals    '$errorctr' errors    '$warnctr' warns    in file '$file' which had '$linectr' lines");
    }
    else{
	$logger->fatal("Found    '$fatalctr' fatals    '$errorctr' errors    '$warnctr' warns    in file '$file' which had '$linectr' lines");
    }


    $fatalmaster += $fatalctr;
    $errormaster += $errorctr;
    $warnmaster += $warnctr;
}


if (($fatalmaster > 0) or ($errormaster > 0) or ($warnmaster > 0)){
    my $subject  = "Ran check_logfiles.pl ";
    if ($fatalmaster > 0){
	$subject .= "FATAL ";
    }
    if ($errormaster > 0){
	$subject .= "ERROR ";
    }
    if ($warnmaster > 0 ){
	$subject .= "WARN ";
    }
    $subject .= "detected";

    my $body;

    $body .= "Workflow link http://xmen:8080/tigr-scripts/papyrus/cgi-bin/show_pipeline.cgi?xmltemplate=" . $workflow_id . "\n\n" if (defined($workflow_id));

    $body  .= "The following log4perl logfiles were scanned:\n\n@$list\n\nTotal fatals '$fatalmaster' total errors '$errormaster' total warns '$warnmaster'.\n\nThe unique log messages reported and their number of occurrences are listed below\n\n";

    foreach my $file (sort keys %{$filehash}){


	$body .= "\nLOGFILE '$file'\n\n";

	foreach my $prgline (sort keys %{$filehash->{$file}}){

	    my $count     = $filehash->{$file}->{$prgline}->{'count'};
	    my $samplemsg = $filehash->{$file}->{$prgline}->{'samplemsg'};
	    my $level     = $filehash->{$file}->{$prgline}->{'level'};
	    $body .= "LEVEL: $level\nPROGRAM LINE: $prgline\nSAMPLE MESSAGE: $samplemsg\nCOUNT: $count\n\n";
	}

    }


    $body .= "Please review logfile '$log4perl'";


    &send_notification($username, $subject, $body);
    
    $logger->fatal("Total fatals '$fatalmaster' total errors '$errormaster' total warns '$warnmaster'. Please review $log4perl");
}

$logger->info("Number of log4perl logfiles processed was '$filectr'. Please review '$log4perl'");


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
#                                  END OF MAIN  -- SUBROUTINES FOLLOW
#
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



#------------------------------------------------------
# write_report()
#
#------------------------------------------------------
sub write_report {

    my ($subject, $body, $rep) = @_;

    my $reportfile = $rep . '/check_logfiles.report.txt';
    if (-e $reportfile){
	rename($reportfile, "$reportfile.bak");
    }
    open (REPORTFILE, ">$reportfile") or $logger->logdie("Could not open report file '$reportfile'");
    print REPORTFILE $subject . "\n\n" . $body . "\n";
    close REPORTFILE;

}

#------------------------------------------------------
# send_notification()
#
#------------------------------------------------------
sub send_notification {

    my ($username, $subject, $body) = @_;

    my $mailer = Mail::Mailer->new ('sendmail');
    $mailer->open({
	             To      => $username,
		     Subject => $subject,
		     From    => $username
		 }) or $logger->logdie("Could not create and send message");
    
    print $mailer $body;
    
    $mailer->close;

    $logger->debug("Notification sent to $username") if $logger->is_debug;


}



#------------------------------------------------------
# get_file_list()
#
#------------------------------------------------------
sub get_file_list {

    my ($list, $rep) = @_;

    #
    # Store command-line specified logfiles
    #
    my @filelist;

    if (defined($list)){
	$logger->info("filelist was defined as '$list'");
	@filelist = split(/,/,$list);
    }

    #
    # Add all log4perl logfiles found under the specified 
    # workflow subdirectory
    #
    if (defined($rep)){

	$logger->info("repository was defined as '$rep'");

	my $string = "find $rep -name \"*.log\" -type f";
	$logger->debug("string '$string'") if $logger->is_debug;
	my @workflowlist = qx{$string};

	chomp @workflowlist;

	foreach my $log4perl (@workflowlist){

	    push (@filelist, $log4perl);
	    
	}
    }
    return \@filelist;
}


#------------------------------------------------------
# show_count()
#
#------------------------------------------------------
sub show_count{
    
    $logger->info("Entered show_count");

    my $string = shift;
    $logger->logdie("string was not defined") if (!defined($string));

    print "\b"x(30);
    printf "%-30s", $string;

}#end sub show_count()


#-------------------------------------------------------------------
# is_file_readable()
#
#-------------------------------------------------------------------
sub is_file_readable {

    my ( $file) = @_;

    $logger->fatal("file was not defined") if (!defined($file));

    my $fatal_flag=0;

    if (!-e $file){
	$logger->fatal("$file does not exist");
	$fatal_flag++;
    }

    else{#if ((-e $file)){
	if ((!-r $file)){
	    $logger->fatal("$file does not have read permissions");
	    $fatal_flag++;
	}
	if ((-z $file)){
	    $logger->fatal("$file has no content");
	    $fatal_flag++;
	}
    }


    return 0 if ($fatal_flag>0);
    return 1;
   

}#end sub is_file_readable()

#------------------------------------------------------
# print_usage()
#
#------------------------------------------------------
sub print_usage {

    print STDERR "SAMPLE USAGE:  $0 [-d debug_level] [-f filelist] [-h] [-l log4perl] [-m] [-r repository] [-w workflow_id] -U username\n";
    print STDERR "  -d|--debug_level         = Optional - Coati::Logger log4perl logging level.  Default is 0\n";    
    print STDERR "  -f|--filelist            = Optional - Comma-separated list of log4perl logfiles to check for FATAL, ERROR, WARN\n";
    print STDERR "  -h|--help                = Optional - Display pod2usage help screen\n";
    print STDERR "  -l|--log4perl            = Optional - Log4perl log file (default: /tmp/check_logfiles.pl.log)\n";
    print STDERR "  -m|--man                 = Optional - Display pod2usage pages for this utility\n";
    print STDERR "  -w|--repository          = Optional - $;WORKFLOW_REPOSITORY$;\n";
    print STDERR "  -w|--workflow_id         = Optional - workflow pipeline XML\n";
    print STDERR "  -U|--username            = Username of person to be notified by email (only if -e=1)\n";
    exit 1;

}
