#!/usr/local/bin/perl

use lib("../..", "/usr/local/annotation/PNEUMO/clu_dir/BSML/ANNOTATION/bsml/src");
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
#use English;
use BsmlReader;
use BsmlParserTwig;


my %options = ();
my $results = GetOptions (\%options, 'database=s', 'asmbl_id1|a=s', 'asmbl_id2|b=s', 'bsml_dir|d=s',
                                     'maxgap|m=s', 'minsize|s=s', 'coordinates|c', 'help|h' );
###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $database   = $options{'database'} || 'pneumo';
my $asmbl_id1  = $options{'asmbl_id1'};
my $asmbl_id2  = $options{'asmbl_id2'};
my $maxgap     = $options{'maxgap'} || '10000';
my $minsize    = $options{'minsize'} || '10000';
my $coord_only = $options{'coordinates'};
my $BSML_dir = $options{'bsml_dir'};
$BSML_dir =~ s/\/$//;

if(!$database or !$asmbl_id1 or !$asmbl_id2 or !$BSML_dir or exists($options{'help'})) {
    &print_usage();
}

###-------------------------------------------------------###
my $parser = new BsmlParserTwig;
my %genes;
my $hitsref = {};

foreach my $id ($asmbl_id1, $asmbl_id2) {
    my $bsml_file = "$BSML_dir/asmbl_${id}.bsml";
    if(!-s $bsml_file) {
	print STDERR "The $bsml_file does not exist!  Aborting...\n";
	exit 5;
    } else {
	my $reader = new BsmlReader;
	$parser->parse( \$reader, $bsml_file );
	fillAnnotation(\%genes, $database, $id, $reader);
    }
}


#fillAnnotation(\%genes, $database, $asmbl_id1);
#fillAnnotation(\%genes, $database, $asmbl_id2);

my $currscore;
while(my $line= <STDIN>){
    if($line =~/:/) {
	my($s, $t) = split(/:/, $line);
	my @scores = split(/\s/, $s);
	$currscore = $scores[2];
    }elsif($line =~ /\-?\d+\s+\-?\d+\s+\w+\s+\w+/) {
    #elsif($line =~ /----/){
	#my($id1,$id2) = split(/----/,$line);
	my @ids = split(/\s+/,$line);   
        my ($id1, $id2) = @ids[2,3];
	$id1 =~ s/\s//g;
	$id2 =~ s/\s//g;
	if($id1 ne "GAP" && $id2 ne "GAP"){
	    if(exists $genes{$database}->{$asmbl_id1}->{$id1}){
		if(exists $genes{$database}->{$asmbl_id2}->{$id2}){
		    if($genes{$database}->{$asmbl_id1}->{$id1}->{'yhit'} ne ""){
			if($genes{$database}->{$asmbl_id1}->{$id1}->{'yhit'} ne $id2){
			    if($currscore > $genes{$database}->{$asmbl_id1}->{$id1}->{'yscore'}){
				$genes{$database}->{$asmbl_id1}->{$id1}->{'yhit'} = $id2;
				$genes{$database}->{$asmbl_id1}->{$id1}->{'yscore'} = $currscore;
			    }
			    if($coord_only){
				print "$genes{$database}->{$asmbl_id1}->{$id1}->{'end5'} $genes{$database}->{$asmbl_id1}->{$id1}->{'end3'} $genes{$database}->{$asmbl_id2}->{$id2}->{'end5'} $genes{$database}->{$asmbl_id2}->{$id2}->{'end3'}\n";
			    } 
			}
			elsif($currscore > $genes{$database}->{$asmbl_id1}->{$id1}->{'yscore'}){
			    $genes{$database}->{$asmbl_id1}->{$id1}->{'yscore'} = $currscore;
			}
		    }
		    else{
			if($coord_only){
			    print "$genes{$database}->{$asmbl_id1}->{$id1}->{'end5'} $genes{$database}->{$asmbl_id1}->{$id1}->{'end3'} $genes{$database}->{$asmbl_id2}->{$id2}->{'end5'} $genes{$database}->{$asmbl_id2}->{$id2}->{'end3'}\n";
			} 
			$genes{$database}->{$asmbl_id1}->{$id1}->{'yhit'} = $id2;
		    }
		    
		    if($genes{$database}->{$asmbl_id2}->{$id2}->{'xhit'} ne ""){
			if($genes{$database}->{$asmbl_id2}->{$id2}->{'xhit'} ne $id1){
			    if($currscore > $genes{$database}->{$asmbl_id2}->{$id2}->{'xscore'}){
				$genes{$database}->{$asmbl_id2}->{$id2}->{'xhit'} = $id1;
				$genes{$database}->{$asmbl_id2}->{$id2}->{'xscore'} = $currscore;
			    }
			}
			elsif($currscore > $genes{$database}->{$asmbl_id2}->{$id2}->{'xscore'}){
			    $genes{$database}->{$asmbl_id2}->{$id2}->{'xscore'} = $currscore;
			}
		    }
                    else{
			$genes{$database}->{$asmbl_id2}->{$id2}->{'xhit'} = $id1;
		    }
                    $genes{$database}->{$asmbl_id1}->{$id1}->{'marked'} = 1;
		    $genes{$database}->{$asmbl_id2}->{$id2}->{'marked'} = 1;
		}
		else{
		    print STDERR "Can't find $id2 in $database $asmbl_id2\n";
		}
            }
	    else{
		print STDERR "Can't find $id1 in $database $asmbl_id1\n";
	    }
	}
    }
}

my $xref = $genes{$database}->{$asmbl_id1};
my $yref = $genes{$database}->{$asmbl_id2};


foreach my $feat (keys %$xref){
    my $pos = $xref->{$feat}->{'end5'} > $xref->{$feat}->{'end3'} ? $xref->{$feat}->{'end5'} : $xref->{$feat}->{'end3'};
    if($xref->{$feat}->{'marked'} == 1){
	$hitsref->{$pos}->{'match'} = 1;
	$hitsref->{$pos}->{'xfeat'} = $feat;
	$hitsref->{$pos}->{'yfeat'} = $xref->{$feat}->{'yhit'};
    }
    #$hitsref->{$x}->{'xfeat'} = $feat;
}

my @sortedhits = sort{$a <=> $b} keys(%$hitsref);
my($xstartfeat);
my($ystartfeat);
my($lastxfeat);
my($lastyfeat);
my($xgaplength)=0;
my($ygaplength)=0;
my($step)=$maxgap;
my($run)=0;
my($matches)={};
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
if($xstartfeat ne ""){
    #print "Ending match at $xref->{$lastxfeat}->{'end5'} $yref->{$lastyfeat}->{'end5'}\n";
    $matches->{$matchcount}->{'xstartfeat'} = $xstartfeat;
    $matches->{$matchcount}->{'ystartfeat'} = $ystartfeat;
    $matches->{$matchcount}->{'xstopfeat'} = $lastxfeat;
    $matches->{$matchcount}->{'ystopfeat'} = $lastyfeat;
}

if(!$coord_only){
    foreach my $match (sort {$matches->{$a} <=> $matches->{$b}} (keys %$matches)){
	my($xdist) =$xref->{$matches->{$match}->{'xstopfeat'}}->{'end5'} -  $xref->{$matches->{$match}->{'xstartfeat'}}->{'end5'};
	my($ydist) =$yref->{$matches->{$match}->{'ystopfeat'}}->{'end5'} - $yref->{$matches->{$match}->{'ystartfeat'}}->{'end5'};
	if(abs($xdist)>= $minsize || abs($ydist)>= $minsize){
	    print "#Match $match $xdist $ydist ($matches->{$match}->{'xgap'}:$matches->{$match}->{'xgapfeat'} $matches->{$match}->{'ygap'}:$matches->{$match}->{'ygapfeat'})"; 
	    print "$matches->{$match}->{'xstartfeat'} $matches->{$match}->{'ystartfeat'} $matches->{$match}->{'xstopfeat'} $matches->{$match}->{'ystopfeat'}\n";
	    print "$xref->{$matches->{$match}->{'xstartfeat'}}->{'end5'} $yref->{$matches->{$match}->{'ystartfeat'}}->{'end5'}\n";
	    print "$xref->{$matches->{$match}->{'xstopfeat'}}->{'end5'} $yref->{$matches->{$match}->{'ystopfeat'}}->{'end5'}\n";
	}
    }
}
elsif($coord_only == 2){
    foreach my $feat (keys %$xref){
        if($xref->{$feat}->{'marked'} == 1){
	    print "$xref->{$feat}->{'end5'} $xref->{$feat}->{'end3'} $yref->{$xref->{$feat}->{'yhit'}}->{'end5'} $yref->{$xref->{$feat}->{'yhit'}}->{'end3'}\n";
	}
    }
} 



sub getDistance{

    my ($f1, $f2) = @_;
    if($f1 eq ""){
	print STDERR "At first hit using 0 as dist (f1)\n";
	#return $f2->{'end5'};
	return 0;
    }
    elsif($f2 eq ""){
	print STDERR "At first hit using 0 as dist (f2)\n";
	#return $f1->{'end5'};
	return 0;
    }
    else{
#	my($o1) = $f1->{'end5'} < $f1->{'end3'} ? 0 : 1;
#	my($o2) = $f2->{'end5'} < $f2->{'end3'} ? 0 : 1;
#	if($o1 != $o2){
#	    return $step + 1;
#	}
#	elsif($o1 == 1){
#	    return $f2->{'end3'} - $f1->{'end5'};
#	}
#	elsif($o1 == 0){
#	    return $f2->{'end5'} - $f1->{'end3'};
#	}
	return abs($f2->{'end5'} - $f1->{'end5'});
    }
}

sub fillAnnotation  {

    my $genes    = shift;
    my $database = shift;
    my $asmbl_id = shift;
    my $reader   = shift;

    my $gene_pos = $reader->fetch_gene_positions("PNEUMO_${asmbl_id}");
    foreach (@$gene_pos) {
	foreach my $gene (keys %$_) {
	    my $end5 = $_->{$gene}->{'startpos'};
	    my $end3 = $_->{$gene}->{'endpos'};
	    $genes->{$database}->{$asmbl_id}->{$gene}->{'end5'} = $end5;
	    $genes->{$database}->{$asmbl_id}->{$gene}->{'end3'} = $end3;
	}
    }

}


sub print_usage {


    print STDERR "SAMPLE USAGE:  gene_boundaries_bsml.pl -d bsml_dir -a 1 -b 5 --database PNEUMO < pe.out > pe.region\n";
    print STDERR "  --database   = database name\n";
    print STDERR "  --bsml_dir   = directory containing bsml doc\n";
    print STDERR "  --asmbl_id1(-a) = reference asmbl_id \n";
    print STDERR "  --asmbl_id2(-b) = query asmbl_id \n";
    print STDERR "  *optional: --maxgap = maximum gap (default 10000)\n";
    print STDERR "             --minsize = minimum size (default 10000)\n";
    print STDERR "  --help = This help message.\n";
    exit 1;

}








