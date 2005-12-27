#!/usr/local/bin/perl

=head1 NAME

region2bsml - Convert gene_boundaries(region) output to BSML

=head1 SYNOPSIS

USAGE: region2bsml -r region_file -o output_file [--debug debug_level] [--log log_file] [--help]

=head1 DESCRIPTION

    Convert gene_boundaries(region) output to BSML

=cut


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use File::Basename;
BEGIN {
    require '/usr/local/devel/ANNOTATION/ard/chado-v1r5b1/lib/site_perl/5.8.5/BSML/BsmlBuilder.pm';
    import BSML::BsmlBuilder;
}
use Log::Log4perl qw(get_logger :levels :easy);
use Pod::Usage;

my ($region,$output_file,$debug,$log,$help,$class);
my $results = GetOptions (
			  'region|r=s' => \$region,
			  'output|o=s' => \$output_file,
			  'debug|D=s'  => \$debug,
			  'log|l=s'    => \$log,
			  'help|?|h'   => \$help,
			  'class|c=s'  => \$class
			  );

pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if($help || !$region || !$output_file);

###-------------------------------------------------------###
#Init logger
#Use central install log file as default if available
Log::Log4perl-> Log::Log4perl::init_and_watch($ENV{LOG4PERL_CONF}) if($ENV{LOG4PERL_CONF});

my $logger = get_logger('ergatis::regions');
$logger->level($INFO);
$logger->more_logging($debug) if($debug);
# Define a file appender or a screen appender
if($log){
    my $file_appender = Log::Log4perl::Appender->new(
						     "Log::Dispatch::File",
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


###-------------------------------------------------------###


if (!defined($class)){
    $logger->logdie("class was not defined");
}



my $doc = BSML::BsmlBuilder->new();

open (IN, $region) or die "unable to read $region due to $!";
while (my $line = <IN>) {
    if($line !~ /^\#/) {
	chomp($line);
	my @elements = split(/\s+/, $line);
	my ($refseq,$compseq,$refstart,$compstart,$refend,$compend, $derived_class) = @elements;

	if (!defined($derived_class)){
	    $derived_class = $class;
	}

	$logger->debug("Parsing match line $refseq,$compseq,$refstart,$compstart,$refend,$compend");
	if (scalar(@elements) == 6) {
	    my $align = $doc->createAndAddSequencePairAlignment('refseq' => $refseq, 'compseq' => $compseq, 'class' => $derived_class);
	    my $complement= ($compstart > $compend) ? 1 : 0;
	    if($complement){
		($compstart,$compend) = ($compend,$compstart);
	    }
	    if($refstart >= $refend){
		$logger->logdie("Bad reference coordinates $refstart- $refend in output between $refseq $compseq\n");
	    }
	    my $query_align_length = abs($refend-$refstart);
	    my $match_align_length = abs($compend-$compstart);
	    $doc->createAndAddSequencePairRun('alignment_pair' => $align, 
					      'refpos' => $refstart, 
					      'runlength' =>$query_align_length,
					      'complement' => 0,
                                              'comppos'=> $compstart, 
					      'comprunlength'=>$match_align_length,
					      'compcomplement'=>$complement);
        } else {
	    $logger->info("Bad match line: $line");
        }
    }
}

# write the altered file to disk
$doc->write( $output_file );
chmod 0666, $output_file;
$logger->info("Wrote output file $output_file");

