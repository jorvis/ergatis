#!/usr/local/bin/perl

=head1  NAME 

dummy.pl - do nothing

=head1 SYNOPSIS

USAGE:  dummy.pl --debug debug_level --log log_file

=head1 OPTIONS

=item *

B<--debug,-d> Debug level.  Use a large number to turn on verbose debugging. 

=item *

B<--log,-l> Log file

=item *

B<--help,-h> This help message

=back

=head1   DESCRIPTION

=cut


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Workflow::Logger;
use BSML::BsmlBuilder;
use BSML::BsmlRepository;
use BSML::BsmlParserSerialSearch;

my %options = ();
my $results = GetOptions (\%options, 
			  'file|f=s',
			  'bsml_repository|b=s',
			  'output|o=s',
			  'num_hits|n=s',
			  'log|l=s',
			  'debug=s',
			  'help|h') || pod2usage();

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Workflow::Logger::get_logger();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDERR} );
}

$options{'num_hits'} = 1 if(!$options{'num_hits'});

&check_parameters(\%options);

my $aalookup = &get_aa_lookup($options{'bsml_repository'});

my $doc = BSML::BsmlBuilder->new();

my $seq_pair = {};

my ($ClusterGeneCount, $ClusterGapCount, $ClusterScore, $ClusterId);

open (PE, "$options{'file'}") or die "Unable to open \"$options{'file'}\" due to $!";
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
	my $length = $aalookup->{$pe[2]}->{'length'} || 0;
	my $comprunlength = $aalookup->{$pe[3]}->{'length'} || 0;
	push @{$seq_pair->{$pe[2]}->{$pe[3]}}, {
	    'start_query'    => 0,
	    'runlength'      => $length,
	    'comprunlength'  => $comprunlength,
	    'start_hit'      => 0,
	    'runscore'       => $ClusterScore,
	    'PEffect_Cluster_Id' => $ClusterId,
	    'PEffect_Cluster_Gap_Count' => $ClusterGapCount,
	    'PEffect_Cluster_Gene_Count' => $ClusterGeneCount 
	    };
    }
}

close PE;

$logger->debug("Created lookup with ",scalar(keys %$seq_pair)," matches") if($logger->is_debug());

foreach my $query_name (keys %$seq_pair){
    if($query_name ne "GAP"){
	foreach my $dbmatch_accession (keys %{$seq_pair->{$query_name}}){
	    if($dbmatch_accession ne "GAP"){
		if( !( $doc->returnBsmlSequenceByIDR($query_name) )){
		    $doc->createAndAddSequence($query_name,$query_name,$aalookup->{$query_name}->{'length'}, 'aa' );
		    my $seq = $doc->returnBsmlSequenceByIDR($query_name);
		    $seq->addBsmlAttr('ASSEMBLY',$aalookup->{$query_name}->{'asmbl'});
		}
		if( !( $doc->returnBsmlSequenceByIDR($dbmatch_accession) )){
		    $doc->createAndAddSequence($dbmatch_accession,$dbmatch_accession,$aalookup->{$dbmatch_accession}->{'length'}, 'aa' );
		    my $seq = $doc->returnBsmlSequenceByIDR($dbmatch_accession);
		    $seq->addBsmlAttr('ASSEMBLY',$aalookup->{$dbmatch_accession}->{'asmbl'});
		}
		my $aln = $doc->createAndAddSequencePairAlignment( 'refseq'  => $query_name,
								   'compseq' => $dbmatch_accession 
								   );
		
		my $runs = $seq_pair->{$query_name}->{$dbmatch_accession};
		my $hitnum = 0;
		foreach my $run (sort {$b->{'runscore'} <=> $a->{'runscore'}} (@$runs)){
		    if($hitnum <$options{'num_hits'}){
			my $s = $doc->createAndAddSequencePairRun( 'alignment_pair' => $aln,
								   'refpos'    => $run->{'start_query'},
								   'runlength'      => $run->{'runlength'},
								   'comprunlength'      => $run->{'comprunlength'},
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
    }
 }

$doc->write($options{'output'});

sub get_aa_lookup{
    my($repository) = @_;
    
    my $bsmlrepo = new BSML::BsmlRepository('BSML_repository'=>$repository);
    my ($files) = $bsmlrepo->list_bsml_files();

    my $lookup = {};
    foreach my $bsml_doc (@$files) {
	my $seqParser = new BSML::BsmlParserSerialSearch( ReadFeatureTables => 0,
							  SequenceCallBack =>sub 
							  {
							      my $seqRef = shift;
							      my $id = $seqRef->returnattr( 'id' );
							      my $type = $seqRef->returnattr( 'molecule' );
							      my $length = $seqRef->returnattr( 'length' );
							      if($type eq "aa"){
								  $logger->debug("Storing sequence $id in lookup with ".$seqRef->returnBsmlAttr('ASSEMBLY'))  if($logger->is_debug());
								  $lookup->{$id}->{'asmbl'} = $seqRef->returnBsmlAttr('ASSEMBLY');
								  $lookup->{$id}->{'length'} = $length;
							      }
							  }); 
	
	
	$seqParser->parse($bsml_doc);
    }
    return $lookup;
}


sub check_parameters{
    my ($options) = @_;
    
    if(0){
	pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});    
    }
}

