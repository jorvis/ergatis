#!/usr/local/bin/perl
#!/usr/local/bin/perl

=head1 NAME

run_pe.pl - Run position effect program

=head1 SYNOPSIS

USAGE:  run_pe.pl -g genefile -p matchfile [-c] [-b pe_binary] [-l log] [-d debug]

=head1 OPTIONS

=over 8

=item B<--pegene,-g>
    PE format XML file containing gene lists

=item B<--pematch,-p>
    PE format XML file containing gene matches

=item B<--pe_binary,-p>
    Optional. PE format XML file containing gene matches

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

my ($pegene,$pematch,$pebin,$debug,$condor,$log,$help);
my $results = GetOptions ('pegene|g=s' => \$pegene, 
			  'pematch|m=s' => \$pematch, 
			  'pebin|b=s' => \$pebin,
			  'condor|c' => \$condor,
			  'debug|D=s' => \$debug,
			  'log|l=s' => \$log,
			  'help|?|h' => \$help);

pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if($help || !$pegene || !$pematch);
Log::Log4perl-> Log::Log4perl::init_and_watch($ENV{LOG4PERL_CONF}) if($ENV{LOG4PERL_CONF});
my $logger = get_logger('papyrus::pe');
$logger->level($INFO);
$logger->more_logging($debug) if($debug);

$pebin = "/usr/local/devel/ANNOTATION/shared/bin/linux/peffect" if($pebin ne "");

my $logger = get_logger();
$logger->debug("PE binary set to $pebin");

my $outputref;
if($condor){
    $outputref = &run_pe_condor($pebin,$pegene,$pematch);
}
else{
    $outputref = &run_pe($pebin,$pegene,$pematch);
}

#$resultsref = &run_pe($pebin,$pegene,$pematch)
#$pebin = position effect binary
#$pegene = position effect gene XML file
#$pematch = position effect match XML file
#$resultsref = reference to string containing output
sub run_pe{
    my($pebin,$pegene,$pematch) = @_;
    $logger->info("Running local execution of '$pebin  -w 10 -g -50 -r -100 -m 4 -o 3 -f $pegene < $pematch 2>&1'");
    my $output = qx($pebin  -w 10 -g -50 -r -100 -m 4 -o 3 -f $pegene < $pematch 2>&1);
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
    $logger->info("Running local execution of '$pebin  -w 10 -g -50 -r -100 -m 4 -o 3 -f $pegene < $pematch 2>&1'");
    my $output = qx($pebin  -w 10 -g -50 -r -100 -m 4 -o 3 -f $pegene < $pematch 2>&1);
    if($?) {
	$logger->fatal("Unable to execute $pebin");
	return;
    }else {
	return \$output;
    }
}










