#!/usr/local/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BSML::BsmlReader;
use BSML::BsmlParserTwig;


my %options = ();
my $results = GetOptions (\%options, 'database=s', 'asmbl_id|a=s@', 'bsml_file=s', 'bsml_dir|d=s',
                                     'maxgap|m=s', 'minsize|s=s', 'coordinates|c', 'help|h' );
###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $database   = $options{'database'};
my $bsml_file = $options{'bsml_file'};
my $maxgap     = $options{'maxgap'} || '10000';
my $minsize    = $options{'minsize'} || '10000';
my $BSML_dir = $options{'bsml_dir'};
$BSML_dir =~ s/\/$//;

if(!$database or $options{'asmbl_id'} or !$BSML_dir or exists($options{'help'})) {
    &print_usage();
}

my @asmbls = @{$options{'asmbl_id'}};


###-------------------------------------------------------###


my $parser = new BSML::BsmlParserTwig;
my $reader = new BsmlReader;
$parser->parse( \$reader, $bsml_file );

my %asmbllookup;

foreach my $asmbl_id1 (@asmbls){
    &fillAnnotation(\%asmbllookup, $asmbl_id1, $reader);
    foreach my $asmbl_id2 (@asmbls){
	&fillAnnotation(\%asmbllookup, $asmbl_id2, $reader);
	&populateMatches($reader,$asmbl_id1,$asmbl_id2,\%asmbllookup);
	my($xref,$yref,$matches) = &getGroupings($asmbl_id1,$asmbl_id2,\%asmbllookup);
	&printGroupings($matches,$xref,$yref,$asmbl_id1,$asmbl_id2,$minsize);
    }
}

sub populateMatches{
    my ($reader, $query_asmbl_id, $match_asmbl_id,$asmbl_lookup) = @_;
    foreach my $seq (@{$reader->assemblyIdtoSeqList($query_asmbl_id )}) {   #grab all seq obj given query_asmbl_id
	my $seq_id = $seq->returnattr( 'id' );                              #return seq_id of an seq obj
	foreach my $aln (@{$reader->fetch_all_alignmentPairs( $seq_id )}) { #return all alignment with query as $seq_id
	    if(ref($aln)) {
		my $match_ref = $reader->readSeqPairAlignment($aln);            #return all pair_runs for an alignment_pair
		my $bestscore=0;
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
			    }
#			    if($coord_only){
#				print "$asmbl_lookup{$database}->{$query_asmbl_id}->{$id1}->{'end5'} $asmbl_lookup{$database}->{$query_asmbl_id}->{$id1}->{'end3'} $asmbl_lookup{$database}->{$match_asmbl_id}->{$id2}->{'end5'} $asmbl_lookup{$database}->{$match_asmbl_id}->{$id2}->{'end3'}\n";
#			    } 
			}
			elsif($bestscore > $asmbl_lookup->{$query_asmbl_id}->{$id1}->{'yscore'}){
			    $asmbl_lookup->{$query_asmbl_id}->{$id1}->{'yscore'} = $bestscore;
			}
		    }
		    else{
#			if($coord_only){
#			    print "$asmbl_lookup{$database}->{$query_asmbl_id}->{$id1}->{'end5'} $asmbl_lookup{$database}->{$query_asmbl_id}->{$id1}->{'end3'} $asmbl_lookup{$database}->{$match_asmbl_id}->{$id2}->{'end5'} $asmbl_lookup{$database}->{$match_asmbl_id}->{$id2}->{'end3'}\n";
#			} 
			$asmbl_lookup->{$query_asmbl_id}->{$id1}->{'yhit'} = $id2;
		    }
			if($asmbl_lookup->{$match_asmbl_id}->{$id2}->{'xhit'} ne ""){
			    if($asmbl_lookup->{$match_asmbl_id}->{$id2}->{'xhit'} ne $id1){
				if($bestscore > $asmbl_lookup->{$match_asmbl_id}->{$id2}->{'xscore'}){
				    $asmbl_lookup->{$match_asmbl_id}->{$id2}->{'xhit'} = $id1;
				    $asmbl_lookup->{$match_asmbl_id}->{$id2}->{'xscore'} = $bestscore;
				}
			    }
			    elsif($bestscore > $asmbl_lookup->{$match_asmbl_id}->{$id2}->{'xscore'}){
				$asmbl_lookup->{$match_asmbl_id}->{$id2}->{'xscore'} = $bestscore;
			    }
			}
			else{
			    $asmbl_lookup->{$match_asmbl_id}->{$id2}->{'xhit'} = $id1;
			}
			$asmbl_lookup->{$query_asmbl_id}->{$id1}->{'marked'} = 1;
			$asmbl_lookup->{$match_asmbl_id}->{$id2}->{'marked'} = 1;
		    }
		    else{
			print STDERR "Can't find $id2 in $match_asmbl_id\n";
		    }
		}
		else{
		    print STDERR "Can't find $id1 in $query_asmbl_id\n";
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

    my $hitsref = createPositionBasedLookup($xref);
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
	print STDERR "At first hit using 0 as dist (f1)\n";
	return 0;
    }
    elsif($f2 eq ""){
	print STDERR "At first hit using 0 as dist (f2)\n";
	return 0;
    }
    else{
	return abs($f2->{'end5'} - $f1->{'end5'});
    }
}

sub fillAnnotation  {

    my $asmbl_lookup = shift;
    my $asmbl_id = shift;
    my $reader = new BsmlReader;  #use local reader for now

    my $parser = new BSML::BsmlParserTwig;
    $parser->parse( \$reader, "BSML_dir/$asmbl_id.bsml" );

    my $idlookup = $reader->returnAllIdentifiers();
    foreach my $seq (keys %{$idlookup->{$asmbl_id}}){
	foreach my $gene (keys %{$idlookup->{$seq}}){
	    foreach my $transcript (keys %{$idlookup->{$gene}}){
		my $fobj = $reader->id_to_object($transcript);
		$asmbl_lookup->{$seq}->{$transcript->{'protein_id'}}->{'end5'} = $fobj->{'start_loc'};
		$asmbl_lookup->{$seq}->{$transcript->{'protein_id'}}->{'end3'} = $fobj->{'stop_loc'};
		
	    }
	}
    }
    return $asmbl_lookup;
}


sub print_usage {


    print STDERR "SAMPLE USAGE:  gene_boundaries_bsml.pl -d bsml_dir -a 1 -b 5 --database PNEUMO < pe.out > pe.region\n";
    print STDERR "  --database   = database name\n";
    print STDERR "  --bsml_dir   = directory containing bsml doc\n";
    print STDERR "  --asmbl_id = list of asmbl_ids (at least 2 asmbl_ids are required) \n";
    print STDERR "  *optional: --maxgap = maximum gap (default 10000)\n";
    print STDERR "             --minsize = minimum size (default 10000)\n";
    print STDERR "  --help = This help message.\n";
    exit 1;

}








