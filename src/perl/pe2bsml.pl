#!/usr/local/bin/perl

=head1  NAME 

pe2bsml.pl  - convert PEffect output files into BSML documents

=head1 SYNOPSIS

USAGE:  pe2bsml.pl -p peffect.out -o pe.bsml

=head1 OPTIONS

=over 4

=item *

B<--pe_file,-p> [REQUIRED] PEffect output file

=item *

B<--output,-o> [REQUIRED] output BSML file

=item *

B<--help,-h> This help message

=back

=head1   DESCRIPTION

pe2bsml.pl is designed to convert PEffect output files into BSML documents.  

Samples:

1. Convert PE file a.pe into a BSML doc a.bsml
   pe2bsml.pl -p a.pe -o a.bsml


NOTE:  

Calling the script name with NO flags/options or --help will display the syntax requirement.


=cut

use strict;
#use Log::Log4perl qw(get_logger);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BSML::BsmlBuilder;
use File::Basename;
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 'pe_file|p=s', 'bsml_dir|b=s', 'verbose|v', 'output|o=s', 
                                     'verbose|v', 'help|h', 'man') || pod2usage();

###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $output     = $options{'output'};
my $pe_file    = $options{'pe_file'};
my $verbose    = $options{'verbose'};
#Log::Log4perl->init("log.conf");
#my $logger = get_logger();

&cmd_check();
###-------------------------------------------------------###


read_pe_output($pe_file);


sub read_pe_output {

    my $pe_out = shift;

    my $doc = BSML::BsmlBuilder->new();

    my $seq_pair = {};

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
	    if(! (exists $seq_pair->{$pe[2]}->{$pe[3]})){
		$seq_pair->{$pe[2]}->{$pe[3]} = [];
	    }
	    push @{$seq_pair->{$pe[2]}->{$pe[3]}}, {
		'start_query'    => $pe[0],
		'runlength'      => 0,
		'start_hit'      => $pe[1],
		'runscore'       => $ClusterScore,
		'PEffect_Cluster_Id' => $ClusterId,
		'PEffect_Cluster_Gap_Count' => $ClusterGapCount,
		'PEffect_Cluster_Gene_Count' => $ClusterGeneCount 
		};
	}
    }

    close PE;

    print STDERR "Created lookup with ",scalar(keys %$seq_pair)," matches\n";

    foreach my $query_name (keys %$seq_pair){
	foreach my $dbmatch_accession (keys %{$seq_pair->{$query_name}}){
	    my $aln = $doc->createAndAddSequencePairAlignment( 'refseq'  => $query_name,
							       'compseq' => $dbmatch_accession 
							     );
	    my $runs = $seq_pair->{$query_name}->{$dbmatch_accession};
	    my $hitnum = 0;
	    foreach my $run (sort {$b->{'runscore'} <=> $a->{'runscore'}} (@$runs)){
		if($hitnum <3){
		    my $s = $doc->createAndAddSequencePairRun( 'alignment_pair' => $aln,
							       'refpos'    => $run->{'start_query'},
							       'runlength'      => $run->{'runlength'},
							       'comppos'      => $run->{'start_hit'},
							       'runscore'       => $run->{'runscore'},
							       );
		    #additional attributes
		    $s->addBsmlAttr( 'PEffect_Cluster_Id',  $run->{'PEffect_Cluster_Id'} );
		    $s->addBsmlAttr( 'PEffect_Cluster_Gap_Count', $run->{'PEffect_Cluster_Gap_Count'} );
		    $s->addBsmlAttr( 'PEffect_Cluster_Gene_Count', $run->{'PEffect_Cluster_Gene_Count'} );

		}
		$hitnum++;
	    }
	}
    
    }

    $doc->createAndAddAnalysis("program" => "peffect", "programversion" => '1.0', 'sourcename' =>$output,
                               "bsml_link_relation" => 'SEQ_PAIR_ALIGNMENTS', 'bsml_link_url' => '#BsmlTables');
    print STDERR "Writing output\n";
    
    $doc->write($output);

    print STDERR "Writing done\n";
    chmod 0666, $output;
}

sub cmd_check {
#quality check

    if( exists($options{'man'})) {
	pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT});
    }   

    if( exists($options{'help'})) {
	pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT});
    }

    if(!$output or !$pe_file) {
	pod2usage({-exitval => 2,  -message => "$0: All the required options are not specified", -verbose => 1, -output => \*STDERR});
    }

}
