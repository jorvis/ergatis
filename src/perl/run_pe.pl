#!/usr/local/bin/perl

=head1 NAME

run_pe.pl - Run position effect program

=head1 SYNOPSIS

USAGE:  run_pe.pl -g genefile -p matchfile [-c] [-g gap penalty] [-b pe_binary] [-l log] [-d debug] [-o output]

=head1 OPTIONS

=over 8

=item B<--pegene,-g>
    PE format XML file containing gene lists

=item B<--pematch,-p>
    PE format XML file containing gene matches

=item B<--output,-o>
    Optional. Output file. Default is STDOUT.  This is required when using the -c option.

=item B<--pebin,-p>
    Optional. PE format XML file containing gene matches

=item B<--gap_penalty,-g>
    Optional. Gap penalty. eg -50.

=item B<--condor,-c>
    Optional. Run on condor

=item B<--log,-l>
    Optional. Log file

=item B<--debug,-d>
    Optional. Debug level

=back

=head1 DESCRIPTION

    run_pe.pl runs the postion effect program

=cut
 
use strict;
use Log::Log4perl qw(get_logger :levels :easy);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

my ($pegene,$pematch,$pebin,$debug,$condor,$outputfile,$gap_penalty,$log,$help);
my $results = GetOptions ('pegene|g=s' => \$pegene, 
			  'pematch|m=s' => \$pematch, 
			  'pebin|b=s' => \$pebin,
			  'gap_penalty|g=s' => \$gap_penalty,
			  'condor|c' => \$condor,
			  'output|o=s' => \$outputfile,
			  'debug|D=s' => \$debug,
			  'log|l=s' => \$log,
			  'help|?|h' => \$help);
pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if($help || !$pegene || !$pematch);
Log::Log4perl-> Log::Log4perl::init_and_watch($ENV{LOG4PERL_CONF}) if($ENV{LOG4PERL_CONF});
my $logger = get_logger('papyrus::pe');
$logger->level($INFO);
$logger->more_logging($debug) if($debug);
# Define a file appender or a screen appender
if($log){
    my $file_appender = Log::Log4perl::Appender->new(
						     "Log::Dispatch::File",
						     mode => "append",
						     filename  => $log);
    
    my $layout = 
	Log::Log4perl::Layout::PatternLayout->new(
						  "%d %p> %F{1}:%L %M - %m%n");
    $file_appender->layout($layout);
    $logger->add_appender($file_appender);
}else{
    my $screen_appender = Log::Log4perl::Appender->new(
						       "Log::Dispatch::Screen");	
    
    $logger->add_appender($screen_appender);
}

$pebin = "/usr/local/devel/ANNOTATION/shared/bin/linux/peffect" if($pebin eq "");

$logger->debug("PE binary set to $pebin");
$gap_penalty = "-50" if($gap_penalty eq "");
$logger->debug("PE gap_penalty set to $gap_penalty");

if($condor){
    &run_pe_condor($pebin,$pegene,$pematch);
}
else{
    my $outputref;
    $outputref = &run_pe($pebin,$pegene,$pematch);
    $logger->debug("Outputing file $outputfile");
    open FILE,">$outputfile" or die "Can't open file $outputfile";
    print FILE $$outputref;
    close FILE;
}

#$resultsref = &run_pe($pebin,$pegene,$pematch)
#$pebin = position effect binary
#$pegene = position effect gene XML file
#$pematch = position effect match XML file
#$resultsref = reference to string containing output
sub run_pe{
    my($pebin,$pegene,$pematch) = @_;
    my $logger = get_logger('papyrus::pe');
    $logger->info("Running local execution of '$pebin  -w 10 -g $gap_penalty -r -100 -m 2 -o 3 -f $pegene < $pematch'");
    my $output = qx($pebin /tmp/done -w 10 -g $gap_penalty -r -100 -m 4 -o 3 -f $pegene < $pematch);
    if($?) {
	$logger->fatal("Unable to execute $pebin");
	return;
    }else {
	return \$output;
    }
}

#$resultsref = &run_pe_condor($pebin,$pegene,$pematch)
#$pebin = position effect binary
#$pegene = position effect gene XML file
#$pematch = position effect match XML file
#$resultsref = reference to string containing output
sub run_pe_condor{
    my($pebin,$pegene,$pematch) = @_;
    my $logger = get_logger('papyrus::pe');
    $logger->info("Running condor execution of '$pebin  -w 10 -g $gap_penalty -r -100 -m 4 -o 3 -f $pegene < $pematch'");
    
    open FILE,">/usr/local/scratch/pe$$.condor.config";
    print FILE <<ENDCONFIG;
executable = $pebin
Requirements = ((Arch == "INTEL") && (OpSys == "LINUX"))
arguments = $outputfile.done -w 10 -g $gap_penalty -r -100 -m 4 -o 3 -f $pegene
input = $pematch
log = /usr/local/scratch/pe$$.condor.log
error = /usr/local/scratch/pe$$.condor.error
output = $outputfile
initialdir = /usr/local/scratch/
notification = error
universe = vanilla 
queue 1

ENDCONFIG

;

    close FILE;
    
    $logger->debug("Running 'condor_submit /usr/local/scratch/pe$$.condor.config'");
    unlink "$outputfile.done" if(-e "$outputfile.done");
    print `condor_submit /usr/local/scratch/pe$$.condor.config`;
    my $counter=0;
    while (! (-e "$outputfile.done")){
	if($counter>6){
	    print STDERR ".";
	    $counter=0;
	}
	sleep 5;
	$counter++;
    }
    $logger->debug("Condor job /usr/local/scratch/pe$$.condor.config finished");
    $logger->debug("Unlinking $outputfile.done");
    unlink "$outputfile.done";
}











