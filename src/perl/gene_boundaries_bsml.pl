#!/usr/local/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BSML::BsmlReader;
use BSML::BsmlParserTwig;
use Data::Dumper;
use File::Basename;

my $DEBUG=0;

my %options = ();
my $results = GetOptions (\%options, 'asmbl_id|a=s@', 'bsml_file|f=s', 'bsml_dir|d=s',
                                     'maxgap|m=s', 'minsize|s=s', 'coordinates|c', 'help|h' );
###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $bsml_file = $options{'bsml_file'};
my $maxgap     = $options{'maxgap'} || '10000';
my $minsize    = $options{'minsize'} || '10000';
my $BSML_dir = $options{'bsml_dir'};
$BSML_dir =~ s/\/$//;

if(!$options{'asmbl_id'} or !$BSML_dir or exists($options{'help'})) {
    &print_usage();
}

my @asmbls = @{$options{'asmbl_id'}};


###-------------------------------------------------------###



my($asmbl_id1) = $asmbls[0];

if(scalar(@asmbls) == 1){
    my @files = <$BSML_dir/*.bsml>;
    foreach (@files) {
        my $basename = basename($_);
        if(my ($newasmbl) = ($basename =~ /(.+)\.bsml/)) {
	    push(@asmbls, $newasmbl);
	    print STDERR "Adding $newasmbl to comparison\n";
        }
    }
}

my $parser = new BSML::BsmlParserTwig;
my $reader = new BSML::BsmlReader;
$parser->parse( \$reader, $bsml_file );

my %asmbllookup;
my %proteinlookup;

print "Parsed $bsml_file\n" if($DEBUG);

&fillAnnotation(\%asmbllookup, \%proteinlookup, $asmbl_id1);

for(my $i=1;$i<@asmbls;$i++){
    my $asmbl_id2 = $asmbls[$i];
    &fillAnnotation(\%asmbllookup, \%proteinlookup, $asmbl_id2);
    &populateMatches($reader,$asmbl_id1,$asmbl_id2,\%asmbllookup,\%proteinlookup);
    my($xref,$yref,$matches) = &getGroupings($asmbl_id1,$asmbl_id2,\%asmbllookup);
    &printGroupings($matches,$xref,$yref,$asmbl_id1,$asmbl_id2,$minsize);
}

#foreach my $asmbl_id1 (@asmbls){
#    print STDERR "Filling for $asmbl_id1\n";
#    &fillAnnotation(\%asmbllookup, \%proteinlookup, $asmbl_id1);
#    foreach my $asmbl_id2 (@asmbls){
#	if($asmbl_id1 ne $asmbl_id2){
#	    &fillAnnotation(\%asmbllookup, \%proteinlookup, $asmbl_id2);
#	    print STDERR "Populating matches for $asmbl_id1 $asmbl_id2\n";
#	    &populateMatches($reader,$asmbl_id1,$asmbl_id2,\%asmbllookup,\%proteinlookup);
#	    print STDERR "Getting groupings\n";
#	    my($xref,$yref,$matches) = &getGroupings($asmbl_id1,$asmbl_id2,\%asmbllookup);
#	    print STDERR "Printing groupings\n";
#	    &printGroupings($matches,$xref,$yref,$asmbl_id1,$asmbl_id2,$minsize);
#	    print STDERR "Done\n";
#	}
#    }
#}

sub populateMatches{
    my ($reader, $query_asmbl_id, $match_asmbl_id,$asmbl_lookup, $protein_lookup) = @_;
    $reader->makeCurrentDocument();
    #print STDERR "$query_asmbl_id Proteins($protein_lookup->{$query_asmbl_id}->{'proteins'})",@{$protein_lookup->{$query_asmbl_id}->{'proteins'}},"\n";
    foreach my $seq (@{$protein_lookup->{$query_asmbl_id}->{'proteins'}}) {   #grab all seq obj given query_asmbl_id
	my $seq_id = $seq->returnattr( 'id' );  #return seq_id of an seq obj
	if($seq->returnattr( 'molecule' ) eq 'aa'){
	    my $alnpairs = $reader->fetch_all_alignmentPairs( $seq_id );
	    foreach my $aln (@$alnpairs) { #return all alignment with query as $seq_id
		if(ref($aln)) {
		    my $match_ref = $reader->readSeqPairAlignment($aln);            #return all pair_runs for an alignment_pair
		
		    if(exists $asmbl_lookup->{$match_asmbl_id}->{$match_ref->{'compseq'}}){
			print  "Checking $match_asmbl_id $match_ref->{'compseq'}\n" if($DEBUG);
			print  "Fetching alignment pairs $seq_id ($aln)\n" if($DEBUG);
			my $bestscore=0;
			print  "Process alignment $aln between  $match_ref->{'refseq'} $match_ref->{'compseq'}\n" if($DEBUG);
			foreach my $pair_run(@{ $match_ref->{'seqPairRuns'} }) {
			    $bestscore = $pair_run->{'runscore'} if($pair_run->{'runscore'} > $bestscore);  #store best bit_score
			}
			
			my $id1 = $match_ref->{'refseq'};
			my $id2 = $match_ref->{'compseq'};
			if(exists $asmbl_lookup->{$query_asmbl_id}->{$id1}){
			    if(exists $asmbl_lookup->{$match_asmbl_id}->{$id2}){
				if($asmbl_lookup->{$query_asmbl_id}->{$id1}->{'yhit'} ne ""){
				    if($asmbl_lookup->{$query_asmbl_id}->{$id1}->{'yhit'} ne $id2){
					if($bestscore > $asmbl_lookup->{$query_asmbl_id}->{$id1}->{'yscore'}){
					    $asmbl_lookup->{$query_asmbl_id}->{$id1}->{'yhit'} = $id2;
					    $asmbl_lookup->{$query_asmbl_id}->{$id1}->{'yscore'} = $bestscore;
					    print  "Saving y best score for $id1 as $id2 $bestscore\n" if($DEBUG);
					}
				    }
				    elsif($bestscore > $asmbl_lookup->{$query_asmbl_id}->{$id1}->{'yscore'}){
					print  "Saving y best score for $id1 as $id2 $bestscore\n" if($DEBUG);
					$asmbl_lookup->{$query_asmbl_id}->{$id1}->{'yscore'} = $bestscore;
				    }
				}
				else{
				    print  "Saving y best score for $id1 as $id2 $bestscore\n" if($DEBUG);
				    $asmbl_lookup->{$query_asmbl_id}->{$id1}->{'yscore'} = $bestscore;
				    $asmbl_lookup->{$query_asmbl_id}->{$id1}->{'yhit'} = $id2;
				}
				if($asmbl_lookup->{$match_asmbl_id}->{$id2}->{'xhit'} ne ""){
				    if($asmbl_lookup->{$match_asmbl_id}->{$id2}->{'xhit'} ne $id1){
					if($bestscore > $asmbl_lookup->{$match_asmbl_id}->{$id2}->{'xscore'}){
					    $asmbl_lookup->{$match_asmbl_id}->{$id2}->{'xhit'} = $id1;
					    $asmbl_lookup->{$match_asmbl_id}->{$id2}->{'xscore'} = $bestscore;
					    print  "Saving x best score for $id2 as $id1 $bestscore\n" if($DEBUG);
					}
				    }
				    elsif($bestscore > $asmbl_lookup->{$match_asmbl_id}->{$id2}->{'xscore'}){
					$asmbl_lookup->{$match_asmbl_id}->{$id2}->{'xscore'} = $bestscore;
					print  "Saving x best score for $id2 as $id1 $bestscore\n" if($DEBUG);
				    }
				}
				else{
				    print  "Saving x best score for $id2 as $id1 $bestscore\n" if($DEBUG);
				    $asmbl_lookup->{$match_asmbl_id}->{$id2}->{'xscore'} = $bestscore;
				    $asmbl_lookup->{$match_asmbl_id}->{$id2}->{'xhit'} = $id1;
				}

				print  "Marking $id1 $id2\n" if($DEBUG);
				print Dumper ($asmbl_lookup->{$query_asmbl_id}->{$id1}),"\n" if($DEBUG);
				print Dumper ($asmbl_lookup->{$match_asmbl_id}->{$id2}),"\n" if($DEBUG);

				$asmbl_lookup->{$query_asmbl_id}->{$id1}->{'marked'} = 1;
				$asmbl_lookup->{$match_asmbl_id}->{$id2}->{'marked'} = 1;
			    }
			    else{
				print  "Can't find $id2 in $match_asmbl_id\n" if($DEBUG);
			    }
			}
			else{
			    print  "Can't find $id1 in $query_asmbl_id\n" if($DEBUG);
			}
		    }
		    else{
			#print STDERR "Can't find $match_ref->{'compseq'} in $match_asmbl_id $match_ref->{'refseq'} ($query_asmbl_id)\n";
		    }
		}
	    }
	}
    }
}
		
sub getGroupings{
    my($asmbl_id1,$asmbl_id2,$asmbl_lookup) = @_;
    my $xref = $asmbl_lookup->{$asmbl_id1};
    #$xref->{$feat_name}->{'end5'} = end5 coord of $feat_name in genome on x axis
    #$xref->{$feat_name}->{'end3'} = end3 coord of $feat_name in genome on x axis
    #$xref->{$feat_name}->{'yhit'} =  best matching genome y feat_name 
    my $yref = $asmbl_lookup->{$asmbl_id2};
    #$yref->{$feat_name}->{'end5'} = end5 coord of $feat_name in genome on y axis
    #$yref->{$feat_name}->{'end3'} = end3 coord of $feat_name in genome on y axis

    print  "XREF $asmbl_id1 ",scalar(keys %$xref),"\n" if($DEBUG);
    print  "YREF $asmbl_id2 ",scalar(keys %$yref),"\n" if($DEBUG);
    my $hitsref = createPositionBasedLookup($xref);
    print STDERR "Found ",scalar(keys %$xref)," hits (",scalar (keys %$hitsref),")\n" if($DEBUG);
    #hash to keep track of hit positions to genome on x axis
    #$pos - position(bp) of match gene on x axis genome.  The larger coord of end5,end3 is used.
    #$hitsref->{$pos}->{'match'} = 1 indicates a match 
    #$hitsref->{$pos}->{'xfeat'} = feat_name of x axis genome at position $pos
    #$hitsref->{$pos}->{'yfeat'} = best matching genome y feat_name that matches xfeat
    
    my @sortedhits = sort{$a <=> $b} keys(%$hitsref);
    my($xstartfeat); #current feat_name on genome x
    my($ystartfeat); #current feat_name on genome y
    my($lastxfeat); #prev feat_name on genome x
    my($lastyfeat); #prev feat_name on genome y
    my($xgaplength)=0; #distance between current feat_name and prev feat_name on genome x
    my($ygaplength)=0; #distance between current feat_name and prev feat_name on genome y
    my($step)=$maxgap; #minimum linear separation(bp) between groupings
    my($matches)={};	#hash containing the groupings
    #$i = group number
    #$matches->{$i}->{'xstartfeat'} = feat_name of first gene in matching region $i in genome on x axis
    #$matches->{$i}->{'xstopfeat'} = feat_name of last gene in matching region $i in genome on x axis
    #$matches->{$i}->{'ystartfeat'} = feat_name of first gene in matching region $i in genome on y axis 
    #$matches->{$i}->{'ystopfeat'} = feat_name of last gene in matching region $i in genome on y axis
    #$matches->{$i}->{'xgap'} = size(bp) of gap in x coord space
    #$matches->{$i}->{'ygap'} = size(bp) of gap in y coord space
    #$matches->{$i}->{'xgapfeat'} = feat_name after gap in x coord space
    #$matches->{$i}->{'ygapfeat'} = feat_name after gap in y coord space
    my($matchcount)=0;
    foreach my $pos (@sortedhits){
	$xgaplength = &getDistance($xref->{$hitsref->{$pos}->{'xfeat'}},$xref->{$lastxfeat});
	$ygaplength = &getDistance($yref->{$hitsref->{$pos}->{'yfeat'}},$yref->{$lastyfeat});
	if($xstartfeat eq ""){
	    #first match
	    $xstartfeat = $hitsref->{$pos}->{'xfeat'};
	    $ystartfeat = $hitsref->{$pos}->{'yfeat'};
#	print "Starting match at $xref->{$xstartfeat}->{'end5'} $yref->{$ystartfeat}->{'end5'}\n";
	}
	elsif($xgaplength>$step || $ygaplength>$step){
	    #print "Ending match at $xref->{$lastxfeat}->{'end5'} $yref->{$lastyfeat}->{'end5'}\n";
	    $matches->{$matchcount}->{'xstartfeat'} = $xstartfeat;
	    $matches->{$matchcount}->{'ystartfeat'} = $ystartfeat;
	    $matches->{$matchcount}->{'xstopfeat'} = $lastxfeat;
	    $matches->{$matchcount}->{'ystopfeat'} = $lastyfeat;
	    $matches->{$matchcount}->{'xgap'} = $xgaplength;
	    $matches->{$matchcount}->{'xgapfeat'} = $hitsref->{$pos}->{'xfeat'};
	    $matches->{$matchcount}->{'ygap'} = $ygaplength;	
	    $matches->{$matchcount}->{'ygapfeat'} = $hitsref->{$pos}->{'yfeat'};
	    $matchcount++;
	    #print "Starting match at $xref->{$xstartfeat}->{'end5'} $yref->{$ystartfeat}->{'end5'}\n";
	    $xstartfeat = $hitsref->{$pos}->{'xfeat'};
	    $ystartfeat = $hitsref->{$pos}->{'yfeat'};
	}
	$lastxfeat = $hitsref->{$pos}->{'xfeat'};
	$lastyfeat = $hitsref->{$pos}->{'yfeat'};
	#print "Adding $xref->{$lastxfeat}->{'end5'} $yref->{$lastyfeat}->{'end5'} ($xgaplength) ($ygaplength)\n";
    }
	#Clean up end of last grouping
    if($xstartfeat ne ""){
	#print "Ending match at $xref->{$lastxfeat}->{'end5'} $yref->{$lastyfeat}->{'end5'}\n";
	$matches->{$matchcount}->{'xstartfeat'} = $xstartfeat;
	$matches->{$matchcount}->{'ystartfeat'} = $ystartfeat;
	$matches->{$matchcount}->{'xstopfeat'} = $lastxfeat;
	$matches->{$matchcount}->{'ystopfeat'} = $lastyfeat;
    }

    print  "Matches ",scalar(keys %$matches),"\n" if($DEBUG);

    return ($xref,$yref,$matches);
}

#printGroupings - Prints the groupings
#$matches->{$i}->{'xstartfeat'} = feat_name of first gene in matching region $i in genome on x axis
#$matches->{$i}->{'xstopfeat'} = feat_name of last gene in matching region $i in genome on x axis
#$matches->{$i}->{'ystartfeat'} = feat_name of first gene in matching region $i in genome on y axis
#$matches->{$i}->{'ystopfeat'} = feat_name of last gene in matching region $i in genome on y axis
#$matches->{$i}->{'xgap'} = size(bp) of gap in x coord space
#$matches->{$i}->{'ygap'} = size(bp) of gap in y coord space
#$matches->{$i}->{'xgapfeat'} = feat_name after gap in x coord space
#$matches->{$i}->{'ygapfeat'} = feat_name after gap in y coord space
#$minsize = minimum size(bp) of matching regions to print
#$xref->{$feat_name}->{'end5'} = end5 coord of $feat_name in genome on x axis
#$yref->{$feat_name}->{'end5'} = end5 coord of $feat_name in genome on y axis
#$xref->{$feat_name}->{'end3'} = end3 coord of $feat_name in genome on x axis
#$yref->{$feat_name}->{'end3'} = end3 coord of $feat_name in genome on y axis
#$asbml_idx = asmbl_id of genome x
#$asmbl_idy = asmbl_id of genome y 
sub printGroupings{
    my($matches,$xref,$yref,$asmbl_idx,$asmbl_idy,$minsize) = @_;
    foreach my $match (sort {$matches->{$a} <=> $matches->{$b}} (keys %$matches)){
	    my($xdist) =$xref->{$matches->{$match}->{'xstopfeat'}}->{'end5'} -  $xref->{$matches->{$match}->{'xstartfeat'}}->{'end5'};
	    my($ydist) =$yref->{$matches->{$match}->{'ystopfeat'}}->{'end5'} - $yref->{$matches->{$match}->{'ystartfeat'}}->{'end5'};
	    if(abs($xdist)>= $minsize || abs($ydist)>= $minsize){
		print "#Match $match $xdist $ydist ($matches->{$match}->{'xgap'}:$matches->{$match}->{'xgapfeat'} $matches->{$match}->{'ygap'}:$matches->{$match}->{'ygapfeat'})"; 
		print "#$asmbl_idx:$matches->{$match}->{'xstartfeat'} $asmbl_idy:$matches->{$match}->{'ystartfeat'} $asmbl_idx:$matches->{$match}->{'xstopfeat'} $asmbl_idy:$matches->{$match}->{'ystopfeat'}\n";
		print "$asmbl_idx $asmbl_idy $xref->{$matches->{$match}->{'xstartfeat'}}->{'end5'} $yref->{$matches->{$match}->{'ystartfeat'}}->{'end5'} $xref->{$matches->{$match}->{'xstopfeat'}}->{'end5'} $yref->{$matches->{$match}->{'ystopfeat'}}->{'end5'}\n";
	    }
	}
}
    

#$xref->{$feat_name}->{'end5'} = end5 coord of $feat_name in genome on x axis
#$xref->{$feat_name}->{'end3'} = end3 coord of $feat_name in genome on x axis
#$xref->{$feat_name}->{'marked'} == 1 IAOI there is a match
#$xref->{$feat_name}->{'yhit'} = best scoring $feat_name in genome on y axis

#returns
#$hitsref->{[x feat max postion]} = hash keyed by coordinate postion of xref genes
#$hitsref->{[x feat max postion]}->{'xfeat'} = matching x feat_name at this postion
#$hitsref->{[x feat max postion]}->{'yfeat'} = matching y feat_name at this postion



sub createPositionBasedLookup{
    my($xref) = @_;
    my $hitsref = {};
    foreach my $feat (keys %$xref){
	my $pos = $xref->{$feat}->{'end5'} > $xref->{$feat}->{'end3'} ? $xref->{$feat}->{'end5'} : $xref->{$feat}->{'end3'};
	print  "POS: $pos $feat [ $xref->{$feat}->{'marked'} ]\n" if($DEBUG);
	print Dumper($xref->{$feat}),"\n" if($DEBUG);
	if($xref->{$feat}->{'marked'} == 1){
	    $hitsref->{$pos}->{'match'} = 1;
	    $hitsref->{$pos}->{'xfeat'} = $feat;
	    $hitsref->{$pos}->{'yfeat'} = $xref->{$feat}->{'yhit'};
	}
    }
    return $hitsref;
}


sub getDistance{

    my ($f1, $f2) = @_;
    if($f1 eq ""){
#	print STDERR "At first hit using 0 as dist (f1)\n";
	return 0;
    }
    elsif($f2 eq ""){
#	print STDERR "At first hit using 0 as dist (f2)\n";
	return 0;
    }
    else{
	return abs($f2->{'end5'} - $f1->{'end5'});
    }
}

sub fillAnnotation  {

    my $asmbl_lookup = shift;
    my $protein_lookup = shift;
    my $asmbl_id = shift;

    if(! exists $asmbl_lookup->{$asmbl_id}){
	
	my $reader = new BSML::BsmlReader;  #use local reader for now
	
	my $parser = new BSML::BsmlParserTwig;
	$parser->parse( \$reader, "$BSML_dir/$asmbl_id.bsml" );
	
	print STDERR "Parsing $BSML_dir/$asmbl_id.bsml\n" if($DEBUG);
	
	my $idlookup = $reader->returnAllIdentifiers();
	
	foreach my $gene (keys %{$idlookup->{$asmbl_id}}){
	    foreach my $transcript (keys %{$idlookup->{$asmbl_id}->{$gene}}){
		my $fobj = BSML::BsmlDoc::BsmlReturnDocumentLookup($transcript);
		my $locs = [];
		foreach my $interval_loc (@{$fobj->returnBsmlSiteLocListR()})
		  {
		      push( @{$locs}, $interval_loc->{'sitepos'});
		  }
		my $start_loc = $locs->[0];
		my $end_loc = $locs->[1];
		$asmbl_lookup->{$asmbl_id}->{$idlookup->{$asmbl_id}->{$gene}->{$transcript}->{'proteinId'}}->{'end5'} = $start_loc;
		$asmbl_lookup->{$asmbl_id}->{$idlookup->{$asmbl_id}->{$gene}->{$transcript}->{'proteinId'}}->{'end3'} = $end_loc;
		#print STDERR "Reading $transcript $start_loc $end_loc\n";

	    }
	}
	$protein_lookup->{$asmbl_id}->{'proteins'} = $reader->assemblyIdtoSeqList($asmbl_id);
	print STDERR "Lookup built for $asmbl_id Proteins($protein_lookup->{$asmbl_id}->{'proteins'})\n" if($DEBUG);

    }
    return $asmbl_lookup;
}


sub print_usage {


    print STDERR "SAMPLE USAGE:  gene_boundaries_bsml.pl -d /usr/local/annotation/PNEUMO/BSML_repository/ -a gbs_799_assembly -f /usr/local/annotation/PNEUMO/BSML_repository/PEffect/gbs_799_assembly_vs_all.pe.bsml > pe.region\n";
    print STDERR "  --database   = database name\n";
    print STDERR "  --bsml_dir   = directory containing bsml doc\n";
    print STDERR "  --asmbl_id = asmbl_id \n";
    print STDERR "  *optional: --maxgap = maximum gap (default 10000)\n";
    print STDERR "             --minsize = minimum size (default 10000)\n";
    print STDERR "  --help = This help message.\n";
    exit 1;

}








