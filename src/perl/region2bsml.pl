#!/usr/local/bin/perl


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
#use BSML::BsmlReader;
#use BSML::BsmlParserTwig;
use File::Basename;
use BSML::BsmlBuilder;


my %options = ();
my $results = GetOptions (\%options, 'pe_region|f=s', 'bsml_dir|b=s', 'output|o=s', 'verbose|v', 'help|h',);


###-------------PROCESSING COMMAND LINE OPTIONS-------------###
my $pe_region   = $options{'pe_region'};
my $output_file = $options{'output'};
my $BSML_dir   = $options{'bsml_dir'};
$BSML_dir =~ s/\/+$//;         #remove terminating '/'s
my $verbose    = $options{'verbose'};

if(!$pe_region or !$output_file  or exists($options{'help'})) {
    #$logger->fatal("Not all of the required options have been defined.  Exiting...");
    &print_usage();
}

###-------------------------------------------------------###

my $doc = BSML::BsmlBuilder->new();

open (IN, $pe_region) or die "unable to read $pe_region due to $!";
while (my $line = <IN>) {
    if($line !~ /^\#/) {
	chomp($line);
	my @elements = split(/\s+/, $line);
	if (@elements == 6) {
	    my $align = $doc->createAndAddSequencePairAlignment('query_name' => $elements[0], 'dbmatch_accession' => $elements[1]);
	    my $query_align_length = $elements[3]-$elements[2]+1;
	    my $match_align_length = $elements[5]-$elements[4]+1;
	    $doc->createAndAddSequencePairRun('alignment_pair' => $align, 'start_query' => $elements[2], 'runlength' =>$query_align_length,
                                              'start_hit'=> $elements[4], 'comprunlength'=>$match_align_length);
        } else {
	    print STDERR "less than 6 elements\n";
        }
    }
}
# write the altered file to disk
$doc->write( $output_file );

