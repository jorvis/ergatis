#!/usr/local/bin/perl


use lib("shared", "/usr/local/annotation/PNEUMO/clu_dir/BSML/bsml/src");
use strict;
#use Log::Log4perl qw(get_logger);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use English;
use BSML::BsmlBuilder;
#use BSML::BsmlReader;
#use BSML::BsmlParserTwig;
use File::Basename;


my %options = ();
my $results = GetOptions (\%options, 'pe_file|p=s', 'bsml_dir|d=s', 'verbose|v', 'output|o=s', 'verbose|v', 'help|h',);

###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $output     = $options{'output'};
my $pe_file     = $options{'pe_file'};
my $BSML_dir   = $options{'bsml_dir'};
$BSML_dir =~ s/\/+$//;         #remove terminating '/'s
my $verbose    = $options{'verbose'};
#Log::Log4perl->init("log.conf");
#my $logger = get_logger();
my $verbose    = $options{'verbose'};

if(!$pe_file or !$output or exists($options{'help'})) {
    #$logger->fatal("Not all of the required options have been defined.  Exiting...");
    &print_usage();
}

###-------------------------------------------------------###


read_pe_output($pe_file);

sub read_pe_output {

    my $pe_out = shift;

    my $doc = BSML::BsmlBuilder->new();

    my ($ClusterGeneCount, $ClusterGapCount, $ClusterScore, $ClusterId);
    open (PE, "$pe_out") or die "Unable to open \"$pe_out\" due to $!";
    while (my $line = <PE> ) {
	chomp($line);
	my @pe = split(/\s+/, $line);
	if(@pe == 5) {
	    #Encountered a cluster definition
	    $ClusterGeneCount = $pe[0];
	    $ClusterGapCount = $pe[1];
	    $ClusterScore = $pe[2];
	    $ClusterId = $pe[4];
	} elsif(@pe == 4) {
	    my $aln = $doc->createAndAddSequencePairAlignment( 'query_name'        => $pe[2],
							       'dbmatch_accession' => $pe[3] 
							       );
	    my $s = $doc->createAndAddSequencePairRun( 'alignment_pair' => $aln,
						       'start_query'    => $pe[0],
						       'runlength'      => 0,
						       'start_hit'      => $pe[1],
						       'runscore'       => $ClusterScore,
						       'PEffect_Cluster_Id' => $ClusterId,
						       'PEffect_Cluster_Gap_Count' => $ClusterGapCount,
						       'PEffect_Cluster_Gene_Count' => $ClusterGeneCount 
						       );
	}
    }
    close PE;
    
    $doc->write($output);
}

