#!/usr/local/bin/perl

=head1 NAME

gene_boundaries_bsml - Clusters gene matches into syntenic regions

=head1 SYNOPSIS

USAGE:  gene_boundaries_bsml -a asmbl_id [-f asmbl_file] -b bsml_file [-c file] -d bsml_directory [-e match_all_asmbls]  [-m min_gap_between_clusters] [-s min_cluster_size] [-l log_file] [--debug debug_level]

=head1 OPTIONS

=over 8

=item B<--asmbl_id,-a>
    
    The asmbl_id to analyze for syntenic regions.  Multiple asmbl_ids can be specified by a comma separated list or through multiple -a parameters.  See --asmbl_file for an alternate forms for input.

=item B<--asmbl_file,-f>
    
    The a file containing a list of asmbl_ids to analyze for syntenic regions.  The file format is one asmbl_id per line.

=item B<--bsml_file,-b>
    
    A BSML file containing a set of matches to analyze.

=item B<--bsml_file_list,-c>
    
    A file containing a list of BSML files.  This is an alternative to -b.

=item B<--bsml_directory,-d>
    
    Directory containing BSML files with the genes referenced in the match file specified with --bsml_file.

=item B<--min_gap,-m>

    Optional. The minimum gap size in bp between clusters.  This is calculated based on gene position mapped to genomic sequence.  The default value is 10000.

=item B<--min_cluster,-s>

    Optional. The minimum size in bp of a syntenic region/cluster.  The size of the cluster is determined by the gene coordinates of the start and end of the cluster.

=item B<--match_all_asmbls,-e>

    Optional. Find regions against all asmbl_ids.  The list of asmbl_ids is retrieved from the BSML directory specified with -d.

=item B<--log,-l>
    
    Optional.  Write debug information to a log file.

=item B<--debug,-d>
    
    Optional.  Debug level

=item B<--help,-h>

    Print this help


=back

=head1 DESCRIPTION

    Clusters gene matches into syntenic regions.

=cut
 
use strict;
use Log::Log4perl qw(get_logger :levels :easy);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BSML::BsmlReader;
use BSML::BsmlParserTwig;
use Data::Dumper;
use File::Basename;
use Pod::Usage;

my $DEBUG=0;

my(@asmbl_ids,$asmbl_file,$bsml_file,$bsml_file_list,$mingap,$mincluster,$BSML_dir,$all_asmbl_flag,$debug,$log,$help);
my $results = GetOptions ('asmbl_id|a=s@' => \@asmbl_ids, 
			  'asmbl_file|f=s' => \$asmbl_file,
			  'bsml_file|b=s' => \$bsml_file,
			  'bsml_file_list|c=s' => \$bsml_file_list,
			  'bsml_dir|d=s' => \$BSML_dir,
			  'min_gap|m=s' => \$mingap,
			  'min_cluster|s=s' => \$mincluster,
			  'match_all_asmbls|e' => \$all_asmbl_flag,
			  'debug|D=s' => \$debug,
			  'log|l=s' => \$log,
			  'help|?|h' => \$help);
###-------------PROCESSING COMMAND LINE OPTIONS-------------###
pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if($help || (!$bsml_file && !$bsml_file_list) || (!@asmbl_ids && !$asmbl_file));

$BSML_dir =~ s/\/$//;
$mingap = 10000 if($mingap eq "");
$mincluster = 3000 if($mincluster eq "");

###-------------------------------------------------------###
#Init logger
#Use central install log file as default if available
Log::Log4perl-> Log::Log4perl::init_and_watch($ENV{LOG4PERL_CONF}) if($ENV{LOG4PERL_CONF});
my $logger = get_logger('papyrus::regions');
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


my %asmbllookup;


########
#Get bsml files containing matches
########
my @allbsmlfiles;

if($bsml_file_list){
    open FILE,$bsml_file_list or die "Can't open $bsml_file_list\n";
    while (my $line = <FILE>){
	chomp $line;
	push @allbsmlfiles,$line;
	$logger->info("Adding matches from $line to analysis");
    }
    close FILE;
}
else{
    $logger->info("Adding matches from $bsml_file to analysis");
    push @allbsmlfiles,$bsml_file;
}

######
#Get all asmbls to use for comparison
######
if($asmbl_file){
    open FILE,$asmbl_file or die "Can't open $asmbl_file\n";
    while (my $line = <FILE>){
	chomp $line;
	push @asmbl_ids,$line;
	$logger->info("Adding asmbl_id $line to analysis");
    }
}
if($all_asmbl_flag){
    $logger->info("Comparing ",join(',',@asmbl_ids), " vs. all asmbl_ids in $BSML_dir");
}
else{
    $logger->info("Comparing ",join(',',@asmbl_ids), " against each other");
}
 
foreach my $curr_bsml_file (@allbsmlfiles){
    my $parser = new BSML::BsmlParserTwig;
    my $reader = new BSML::BsmlReader;
    $parser->parse( \$reader, $curr_bsml_file );

    $logger->info("Parsed match file $curr_bsml_file");
    
    

    foreach my $asmbl_id1 (@asmbl_ids){
	
	&fillAnnotation(\%asmbllookup, "$BSML_dir/${asmbl_id1}.bsml") if(! (exists $asmbllookup{$asmbl_id1}));;;    
	
	my @allasmbls = &get_all_asmbls(\@asmbl_ids,$all_asmbl_flag);
	my $proteinmatches = &getMatchingIds($reader,$asmbllookup{$asmbl_id1});

	foreach my $asmbl_id2 (@allasmbls){ 
	    if($asmbl_id1 ne $asmbl_id2){
		&fillAnnotation(\%asmbllookup, "$BSML_dir/${asmbl_id2}.bsml") if(! (exists $asmbllookup{$asmbl_id2}));;
		my $match_lookup = &populateMatches($reader,$asmbl_id1,$asmbl_id2,\%asmbllookup,$proteinmatches);
		my($xref,$yref,$matches) = &getGroupings($asmbl_id1,$asmbl_id2,\%asmbllookup,$match_lookup);
		&printGroupings($matches,\%asmbllookup,$xref,$yref,$asmbl_id1,$asmbl_id2,$mincluster);
	    }
	}
    }
}

sub getMatchingIds{
    my($reader,$asmblref) = @_;
    my @matchingprots;
    my $logger = get_logger('papyrus::regions');
    my $seqalns = $reader->returnAllSeqPairAlignmentsListR();
    foreach my $seqaln (@$seqalns){
	my $match_ref = $reader->readSeqPairAlignment($seqaln);
	if(ref($match_ref)){
	    $logger->debug("Checking $match_ref->{'refseq'}...");
	    if(exists $asmblref->{$match_ref->{'refseq'}}){
		push @matchingprots,$match_ref;
		$logger->debug("Found seqpair alignment on current assembly $match_ref->{'refseq'}. Match is to $match_ref->{'compseq'}");
	    }
	}
    }
    return \@matchingprots;
}


sub get_all_asmbls{
    my($asmblref,$flag) = @_;
    my @asmbls;
    if($flag){
	my @files = <$BSML_dir/*.bsml>; 
	foreach (@files) {
	    my $basename = basename($_);
	    if(my ($newasmbl) = ($basename =~ /(.+)\.bsml/)) {
		push(@asmbls, $newasmbl);
	    }
	}
    }
    else{
	@asmbls = @$asmblref;
    }

    return @asmbls;
}
sub populateMatches{
    my ($reader, $query_asmbl_id, $match_asmbl_id,$asmbl_lookup, $proteinmatches) = @_;
    my $match_lookup = {};
    $reader->makeCurrentDocument();
    my $logger = get_logger('papyrus::regions');
    $logger->debug ("Searching ",scalar(@$proteinmatches)," proteins on $query_asmbl_id for seqpairs");
    foreach my $match_ref (@$proteinmatches) {   #grab all seq obj given query_asmbl_id
	if($match_ref->{'compseq'} ne $match_ref->{'refseq'}){ #ignore self hit
	    $logger->debug("Checking if $match_ref->{'compseq'} exists in $match_asmbl_id");
	    if(exists $asmbl_lookup->{$match_asmbl_id}->{$match_ref->{'compseq'}}){
		my $bestscore=0;
		$logger->debug("MATCH FOUND: Process alignment between $match_ref->{'refseq'} $match_ref->{'compseq'}");
		foreach my $pair_run(@{ $match_ref->{'seqPairRuns'} }) {
		    $bestscore = $pair_run->{'runscore'} if($pair_run->{'runscore'} > $bestscore);  #store best bit_score
		}
		$logger->debug("Best scoring run $bestscore between $match_ref->{'refseq'} $match_ref->{'compseq'}");
		
		my $id1 = $match_ref->{'refseq'};
		my $id2 = $match_ref->{'compseq'};
		if(exists $asmbl_lookup->{$query_asmbl_id}->{$id1}){
		    if(exists $asmbl_lookup->{$match_asmbl_id}->{$id2}){
			if($match_lookup->{$query_asmbl_id}->{$id1}->{'yhit'} ne ""){
			    if($match_lookup->{$query_asmbl_id}->{$id1}->{'yhit'} ne $id2){
				if($bestscore > $match_lookup->{$query_asmbl_id}->{$id1}->{'yscore'}){
				    $logger->debug("Saving best score $bestscore for $id1(y) as $id2(x).  Replaces $match_lookup->{$query_asmbl_id}->{$id1}->{'yhit'} with $match_lookup->{$query_asmbl_id}->{$id1}->{'yscore'}.");
				    $match_lookup->{$query_asmbl_id}->{$id1}->{'yhit'} = $id2;
				    $match_lookup->{$query_asmbl_id}->{$id1}->{'yscore'} = $bestscore;
				}
				else{
				    $logger->debug("Lower scoring hit $id2(x) found for $id1(y). Ignoring");
				}
			    }
			    elsif($bestscore > $match_lookup->{$query_asmbl_id}->{$id1}->{'yscore'}){
				$logger->debug("Saving best score $bestscore for $id1(y) as $id2(x). Replaces score $match_lookup->{$query_asmbl_id}->{$id1}->{'yscore'}");
				$match_lookup->{$query_asmbl_id}->{$id1}->{'yscore'} = $bestscore;
			    }
			}
			else{
			    $logger->debug("Saving best score $bestscore for $id1(y) as $id2(x).  First hit encountered.");
			    $match_lookup->{$query_asmbl_id}->{$id1}->{'yscore'} = $bestscore;
			    $match_lookup->{$query_asmbl_id}->{$id1}->{'yhit'} = $id2;
			}
			if($match_lookup->{$match_asmbl_id}->{$id2}->{'xhit'} ne ""){
			    if($match_lookup->{$match_asmbl_id}->{$id2}->{'xhit'} ne $id1){
				if($bestscore > $match_lookup->{$match_asmbl_id}->{$id2}->{'xscore'}){
				    $logger->debug("Saving best score $bestscore for $id2(x) as $id1(y).  Replaces $match_lookup->{$query_asmbl_id}->{$id2}->{'xhit'} with $match_lookup->{$query_asmbl_id}->{$id2}->{'xscore'}.");
				    $match_lookup->{$match_asmbl_id}->{$id2}->{'xhit'} = $id1;
				    $match_lookup->{$match_asmbl_id}->{$id2}->{'xscore'} = $bestscore;
				}
				else{
				    $logger->debug("Lower scoring hit $id1(y) found for $id2(x). Ignoring");
				}
			    }
			    elsif($bestscore > $match_lookup->{$match_asmbl_id}->{$id2}->{'xscore'}){
				$logger->debug("Saving best score $bestscore for $id2(x) as $id1(y). Replaces score $match_lookup->{$query_asmbl_id}->{$id2}->{'yscore'}");
				$match_lookup->{$match_asmbl_id}->{$id2}->{'xscore'} = $bestscore;
			    }
			}
			else{
			    $logger->debug("Saving best score $bestscore for $id2(x) as $id1(y).  First hit encountered.");
			    $match_lookup->{$match_asmbl_id}->{$id2}->{'xscore'} = $bestscore;
			    $match_lookup->{$match_asmbl_id}->{$id2}->{'xhit'} = $id1;
			}
			
			$logger->debug("Marking $query_asmbl_id --> $id1\n");
			$logger->debug("Marking $match_asmbl_id --> $id2\n");
			$match_lookup->{$query_asmbl_id}->{$id1}->{'marked'} = 1;
			$match_lookup->{$match_asmbl_id}->{$id2}->{'marked'} = 1;
			
		    }
		    else{
			$logger->debug("Can't find $id2 in $match_asmbl_id. Skipping");
		    }
		}
		else{
		    $logger->debug("Can't find $id1 in $query_asmbl_id. Skipping");
		}
	    }
	    else{
		$logger->debug("Can't find match $match_ref->{'compseq'} in $match_asmbl_id. Skipping");
	    }
	}
	else{
	    $logger->debug("Ignoring self hit");
	}
    }
    return $match_lookup;
}
		
sub getGroupings{
    my($asmbl_id1,$asmbl_id2,$asmbl_lookup,$match_lookup) = @_;

    my $logger = get_logger('papyrus::regions');

    my $xref = $match_lookup->{$asmbl_id1};
    #$xref->{$feat_name}->{'end5'} = end5 coord of $feat_name in genome on x axis
    #$xref->{$feat_name}->{'end3'} = end3 coord of $feat_name in genome on x axis
    #$xref->{$feat_name}->{'yhit'} =  best matching genome y feat_name 
    my $yref = $match_lookup->{$asmbl_id2};
    #$yref->{$feat_name}->{'end5'} = end5 coord of $feat_name in genome on y axis
    #$yref->{$feat_name}->{'end3'} = end3 coord of $feat_name in genome on y axis

    $logger->debug("XREF $asmbl_id1 has ",scalar(keys %$xref)," matches");
    $logger->debug("YREF $asmbl_id2 has ",scalar(keys %$yref)," matches");
    my $hitsref = createPositionBasedLookup($xref,$asmbl_lookup->{$asmbl_id1});
    $logger->debug("Stored ",scalar (keys %$hitsref)," hits");
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
    my($step)=$mingap; #minimum linear separation(bp) between groupings
    $logger->debug("Setting minimum linear separation(bp) between groupings at $step");
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
	$logger->debug("At position $pos checking $hitsref->{$pos}->{'xfeat'} $asmbl_lookup->{$asmbl_id1}->{$hitsref->{$pos}->{'xfeat'}}->{'end5'} against $hitsref->{$pos}->{'yfeat'} $asmbl_lookup->{$asmbl_id2}->{$hitsref->{$pos}->{'yfeat'}}->{'end5'}");
	$xgaplength = &getDistance($asmbl_lookup->{$asmbl_id1}->{$hitsref->{$pos}->{'xfeat'}},$asmbl_lookup->{$asmbl_id1}->{$lastxfeat});
	$ygaplength = &getDistance($asmbl_lookup->{$asmbl_id2}->{$hitsref->{$pos}->{'yfeat'}},$asmbl_lookup->{$asmbl_id2}->{$lastyfeat});
	$logger->debug("Distance between $hitsref->{$pos}->{'xfeat'}-->$lastxfeat $hitsref->{$pos}->{'yfeat'}-->$lastyfeat is $xgaplength(x) and $ygaplength(y)"); 
	if($xstartfeat eq ""){
	    #first match
	    $xstartfeat = $hitsref->{$pos}->{'xfeat'};
	    $ystartfeat = $hitsref->{$pos}->{'yfeat'};
	    $logger->debug("START REGION $xstartfeat $ystartfeat");
	}
	elsif($xgaplength>$step || $ygaplength>$step){
	    $logger->debug("CLOSING REGION $xstartfeat-->$lastxfeat $ystartfeat-->$lastyfeat. $xgaplength or $ygaplength exceeded max separation of $step.");
	    $matches->{$matchcount}->{'xstartfeat'} = $xstartfeat;
	    $matches->{$matchcount}->{'ystartfeat'} = $ystartfeat;
	    $matches->{$matchcount}->{'xstopfeat'} = $lastxfeat;
	    $matches->{$matchcount}->{'ystopfeat'} = $lastyfeat;
	    $matches->{$matchcount}->{'xgap'} = $xgaplength;
	    $matches->{$matchcount}->{'xgapfeat'} = $hitsref->{$pos}->{'xfeat'};
	    $matches->{$matchcount}->{'ygap'} = $ygaplength;	
	    $matches->{$matchcount}->{'ygapfeat'} = $hitsref->{$pos}->{'yfeat'};
	    $matchcount++;
	    $xstartfeat = $hitsref->{$pos}->{'xfeat'};
	    $ystartfeat = $hitsref->{$pos}->{'yfeat'};
	    $logger->debug("START REGION $xstartfeat $ystartfeat");
	}
	$lastxfeat = $hitsref->{$pos}->{'xfeat'};
	$lastyfeat = $hitsref->{$pos}->{'yfeat'};
    }
    #Clean up end of last grouping
    if($xstartfeat ne ""){
	$logger->debug("CLOSING REGION $xstartfeat-->$lastxfeat $ystartfeat-->$lastyfeat. $xgaplength or $ygaplength exceeded max separation of $step.");
	$matches->{$matchcount}->{'xstartfeat'} = $xstartfeat;
	$matches->{$matchcount}->{'ystartfeat'} = $ystartfeat;
	$matches->{$matchcount}->{'xstopfeat'} = $lastxfeat;
	$matches->{$matchcount}->{'ystopfeat'} = $lastyfeat;
    }

    $logger->debug("Stored ",scalar(keys %$matches)," groupings");

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
#$mincluster = minimum size(bp) of matching regions to print
#$xref->{$feat_name}->{'end5'} = end5 coord of $feat_name in genome on x axis
#$yref->{$feat_name}->{'end5'} = end5 coord of $feat_name in genome on y axis
#$xref->{$feat_name}->{'end3'} = end3 coord of $feat_name in genome on x axis
#$yref->{$feat_name}->{'end3'} = end3 coord of $feat_name in genome on y axis
#$asbml_idx = asmbl_id of genome x
#$asmbl_idy = asmbl_id of genome y 
sub printGroupings{
    my($matches,$asmbl_lookup,$xref,$yref,$asmbl_idx,$asmbl_idy,$mincluster) = @_;
    foreach my $match (sort {$matches->{$a} <=> $matches->{$b}} (keys %$matches)){
	    my($xdist) = &getDistance($asmbl_lookup->{$asmbl_idx}->{$matches->{$match}->{'xstopfeat'}},$asmbl_lookup->{$asmbl_idx}->{$matches->{$match}->{'xstartfeat'}});
#$asmbl_lookup->{$asmbl_idx}->{$matches->{$match}->{'xstopfeat'}}->{'end5'} -  $asmbl_lookup->{$asmbl_idx}->{$matches->{$match}->{'xstartfeat'}}->{'end5'};
	    my($ydist) = &getDistance($asmbl_lookup->{$asmbl_idy}->{$matches->{$match}->{'ystopfeat'}},$asmbl_lookup->{$asmbl_idy}->{$matches->{$match}->{'ystartfeat'}});
#$asmbl_lookup->{$asmbl_idy}->{$matches->{$match}->{'ystopfeat'}}->{'end5'} - $asmbl_lookup->{$asmbl_idy}->{$matches->{$match}->{'ystartfeat'}}->{'end5'};
	    if(abs($xdist)>= $mincluster || abs($ydist)>= $mincluster){
		print "#Match $match $xdist $ydist ($matches->{$match}->{'xgap'}:$matches->{$match}->{'xgapfeat'} $matches->{$match}->{'ygap'}:$matches->{$match}->{'ygapfeat'})\n"; 
		print "#$asmbl_idx:$matches->{$match}->{'xstartfeat'} $asmbl_idy:$matches->{$match}->{'ystartfeat'} $asmbl_idx:$matches->{$match}->{'xstopfeat'} $asmbl_idy:$matches->{$match}->{'ystopfeat'}\n";
		print "$asmbl_idx $asmbl_idy $asmbl_lookup->{$asmbl_idx}->{$matches->{$match}->{'xstartfeat'}}->{'end5'} $asmbl_lookup->{$asmbl_idy}->{$matches->{$match}->{'ystartfeat'}}->{'end5'} $asmbl_lookup->{$asmbl_idx}->{$matches->{$match}->{'xstopfeat'}}->{'end5'} $asmbl_lookup->{$asmbl_idy}->{$matches->{$match}->{'ystopfeat'}}->{'end5'}\n";
	    }
	    else{
		$logger->debug("Group $matches->{$match}->{'xstartfeat'}-->$matches->{$match}->{'xstopfeat'}=$xdist $matches->{$match}->{'ystartfeat'}-->$matches->{$match}->{'ystopfeat'}=$ydist is below minimum size of $mincluster");
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
    my($xref,$asmblref) = @_;
    my $hitsref = {};
    my $logger = get_logger('papyrus::regions');
    foreach my $feat (keys %$xref){
	my $pos = $asmblref->{$feat}->{'end5'} > $asmblref->{$feat}->{'end3'} ? $asmblref->{$feat}->{'end5'} : $asmblref->{$feat}->{'end3'};
	$logger->debug("POS: $pos $feat [ $xref->{$feat}->{'marked'} ]");
	if($xref->{$feat}->{'marked'} == 1){
	    $logger->debug("Saving match to $feat as $xref->{$feat}->{'yhit'} in position $pos");
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
    my $asmbl_file = shift;

    my $logger = get_logger('papyrus::regions');

    my $reader = new BSML::BsmlReader;  #use local reader for now
    
    my $parser = new BSML::BsmlParserTwig;
    $parser->parse( \$reader, "$asmbl_file" );
    
    $logger->debug("Parsed asmbl_file $asmbl_file");
	
    my $idlookup = $reader->returnAllIdentifiers();
    
    $logger->debug("Storing ",scalar(keys %$idlookup), " assemblies");
    foreach my $asmbl_id (keys %$idlookup){
	$logger->debug("Storing ",scalar(keys %{$idlookup->{$asmbl_id}}), " genes");
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
		$logger->debug("Storing $idlookup->{$asmbl_id}->{$gene}->{$transcript}->{'proteinId'} with coords $start_loc,$end_loc");
		$asmbl_lookup->{$asmbl_id}->{$idlookup->{$asmbl_id}->{$gene}->{$transcript}->{'proteinId'}}->{'end5'} = $start_loc;
		$asmbl_lookup->{$asmbl_id}->{$idlookup->{$asmbl_id}->{$gene}->{$transcript}->{'proteinId'}}->{'end3'} = $end_loc;
	    }
	}
    }
}









