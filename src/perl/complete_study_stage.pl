#!/usr/local/bin/perl -w

=head1 NAME

 complete_study_stage.pl - Connects to IPD and marks a study stage of type "automatic annotation" with the status "complete"

=head1 SYNOPSIS

# USAGE : perl complete_study_stage.pl --study_stage_id=42


=head1 OPTIONS
 --study_stage_id (-s)	--	ID of the study stage whose status needs to be changed to "complete"

 --log (-l)	--	Path to log file

 --help (-h)	--	Print this information


=head1 DESCRIPTION

 Takes the study stage ID from the options, and verifies that it is of type "automatic annotation".  If it is so, then the status is changed to "complete".
 
 This was intended to be used as an compnent in the Ergatis prokaryotic annotation pipeline during times when the user 
wishes to sync their pipeline progress with a sample and study from the IGS Projects Database (IPD)

=head1 INPUT
 None


=head1 OUTPUT
 None


=head1 CONTACT

	Shaun Adkins
	sadkins@som.umaryland.edu

=cut

use strict;
use warnings;
use IPD::Client;
use IPD::IPDObject::StudyStage;
# Extended processing of command line options
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
# Prints usage message from embedded pod documentation
use Pod::Usage;
# Include packages here 

#############
# CONSTANTS #
#############

my $ss_id;
my $config = "devel";

###########
# GLOBALS #
###########
my %options;
# Log file handle;
my $logfh;
my ($ERROR, $WARN, $DEBUG) = (1, 2, 3);

################
# MAIN PROGRAM #
################
GetOptions(\%options,
	   'study_stage_id|s=s',
	   'config|c=s',
	   'log|l=s', 
	   'help|h'
	  ) or pod2usage();

pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ) if ($options{'help'});

&check_parameters();

exit(0) if $ss_id == -1;	#since no automatic annotation study stage was specified... just exit

IPD::Client->set_client($config);
my $ss = IPD::IPDObject::StudyStage->get_study_stage($ss_id) or die("Cannot create new IPD::IPDObject::StudyStage object... perhaps your study stage ID does not exist $!\n");
$ss->set_status('complete');
$ss->save_study_stage($ss);
print "Study Stage " . $ss->get_id() . " has been marked as \'complete\' in IPD\n";
exit(0);

###############
# SUBROUTINES #
###############

# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing. 
# Parameters    : NA
# Returns       : NA
# Modifications :

sub check_parameters{
	if(exists($options{'log'})) {
		open($logfh, "> $options{'log'}") or die "Could not open $options{'log'} file for writing: $!\n"
	}	

	if (exists($options{'study_stage_id'})){
	 	 $ss_id = $options{'study_stage_id'};
	} else { 
		printLogMsg(1,"No Study Stage ID options present\n") if (exists($options{'log'}));
	}

	if (exists($options{'config'})){
		$config = $options{'config'};
	}
}

####################################################################################################################################################

# Description   : Used to handle logging of messages(errors and warnings) during the execution of the script
# Parameters    : level = can be ERROR, WARNING or INFO
#		  msg   = msg to be printed in the log file or to STDERR
# Returns       : NA
# Modifications : 

sub printLogMsg {
	my ($level, $msg) = @_;
	if( $level <= $DEBUG ) {
		print STDERR "$msg\n";
		print $logfh "$msg\n" if(defined($logfh));
		die "$msg\n" if($level == $ERROR);
	}	
}

__END__

