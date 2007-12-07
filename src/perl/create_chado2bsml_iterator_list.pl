#!/usr/bin/perl
=head1  NAME 

create_chado2bsml_iterator_list.pl will take a new-line separated list of assembly identifiers and create an ergatis workflow iterator list file

=head1 SYNOPSIS

USAGE:  create_chado2bsml_iterator_list.pl --control_file=inputfile --output_file=outputfile

=head1 OPTIONS

=item *

B<--control_file> Input file containing new-line separated list of assembly identifiers

=item *

B<--output_file> Output iterator file that will be processed by the chado2bsml workflow component

=item *

B<--debug,-d> Debug level.  Use a large number to turn on verbose debugging. 

=item *

B<--logfile> Ergatis Workflow::Logger logfile

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
                          'control_file=s', 
			  'output_file=s',
			  'logfile=s',
                          'debug=s', 
                          'help|h' ) || pod2usage();

my $logfile = $options{'logfile'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Workflow::Logger::get_logger();

if (&verify_control_file($options{'control_file'})){

    if (defined($options{'output_file'})){

	my $asmblCtr=0;
	my $lineCtr=0;
	my $asmblLookup={};

	open (INFILE, "<$options{'control_file'}") || $logger->logdie("Could not open control_file '$options{'control_file'}' for input mode: $!");

	while (my $line = <INFILE>){
	    chomp $line;
	    $lineCtr++;

	    if ($line =~ /^\s*$/){
		next; ## skip blank lines
	    }

	    if ($line =~ /^\#+/){
		next;
	    }

	    $line =~ s/^\s+//; ## remove leading white spaces
	    $line =~ s/\s+$//; ## remove trailing white spaces

	    if (! exists $asmblLookup->{$line}){
		$asmblCtr++;
		$asmblLookup->{$line}++;
	    }
	    else {
		$logger->warn("Found duplicate asmbl '$line' at line '$lineCtr' of control_file '$options{'control_file'}'");
	    }		
	}

	close INFILE;

	if ($asmblCtr>0){

	    open (OUTFILE, ">$options{'output_file'}") || $logger->logdie("Could not open output_file '$options{'output_file'}' for output mode: $!");

	    print OUTFILE  '$;ASMBL$;' . "\n";

	    foreach my $asmblId (sort keys %{$asmblLookup} ){
		print OUTFILE "$asmblId\n";
	    }

	    close OUTFILE;
	}
	else {
	    $logger->logdie("Did not find any assembly identifiers in control_file '$options{'control_file'}'");
	}
    }
    else{
	$logger->logdie("output_file was not defined");
    }
}


exit;
						     
#---------------------------------------------------------------------------------------------------------
#
#                           END OF MAIN  --  SUBROUTINES FOLLOW
#
#---------------------------------------------------------------------------------------------------------
 
#-------------------------------------------
# verify_control_file()
#
#-------------------------------------------
sub verify_control_file {
    my $file = shift;

    if ((defined($file)) && (-e $file) && (-r $file) && (-s $file)) {
	# control file was defined, exists, has read permissions
	# and has content
	return 1;
    }
    else {
	if (!defined($file)){
	    $logger->logdie("control file '$file' was not defined");
	}
	if (!-e $file){
	    $logger->logdie("control file '$file' does not exist");
	}
	if (!-r $file){
	    $logger->logdie("control file '$file' does not have read permissions");
	}
	if (!-s $file){
	    $logger->logdie("control file '$file' has zero content");
	}
    }
}
