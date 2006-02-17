#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

# This program takes as input a single BTAB file containing
# 1-way best hits, typically between the proteins in 2 or more 
# distinct (and preferably complete) genomes.  It generates 
# clusters of proteins that are connected by bidirectional 
# best hits.  In the special case of two genomes it will find 
# clusters that correspond exactly to all pairs of bidirectional
# best hits.

# A special 22nd column in the BTAB input may be used to store
# the ids of proteins that have already been clustered with the
# protein named in the first column.  This is commonly used 
# to incorporate a set of within-genome Jaccard cluster analyses
# into the best hit analysis.  When a protein in one of these
# input clusters is placed into an output cluster, all of the
# other proteins in that input cluster will be added to the 
# output also.

require "getopts.pl";
$usage = <<_MESSAGE_;
Usage:
best_hit.pl [options] > output_file
    -D     DEBUG              default = off        show diagnostic messages
    -h     help                                    print this message
    -c     cutoff             default = 0.0001     best hits with an e-value above this cutoff will be ignored
    -i     input file                              input btab file
_MESSAGE_

&Getopts('Dhi:c:');
die "$usage\n" if ($opt_h == 1);

$DEBUG = $opt_D;

$cutoff = (length($opt_c) == 0) ? 0.0001 : $opt_c;

$DEBUG && print ("DEBUG is ON, cutoff is $cutoff\n");

my %jaccard;

## open the input btab file
open(my $ifh, "<$opt_i") || die "can't open input file: $!";

# first load all the best hits data
while ($btab = <$ifh>) {
    chomp $btab;
    @data = split (/\t/,$btab);
    if ($data[20] > $cutoff)  {next;}    # discard if a weak hit by E-value.        
    if ($data[20] < 0)        {next;}    # a little quality control; E-value should always be positive.
    if ($data[5] eq "")       {next;}    # a little quality control
    $data[0] =~ s/\-/:/g;
    ($from_tag, $nothing, $from_gi) = split(/:/,$data[0]);
    $data[5] =~ s/[\|\-]/:/g;            # because the YEAST file has -, not |
    ($to_tag, $nothing, $to_gi) = split(/:/,$data[5]);
    if ($to_tag eq $from_tag) {
	# $DEBUG && print ("Oops!, self-self of $from_tag:$to_tag  $to_gi\n");
        next;
    }
    # found a good one: beat the cutoff, not a self-self comparison.
    # Put in a hash.
    # To avoid duplicate info, store in only one direction...
    if ($to_tag gt $from_tag) {
	$munge = $from_tag."..".$to_tag;
        $all_best_hits{$munge}->{'from_score'} = $data[20]; 
        $all_best_hits{$munge}->{'from_Nterm'} = $data[6];  # gets overwritten if bidirectional.
        $all_best_hits{$munge}->{'from_Cterm'} = $data[7];  # 
        $all_best_hits{$munge}->{'to_Nterm'} = $data[8];    #
        $all_best_hits{$munge}->{'to_Cterm'} = $data[9];    #
        $jaccard{$from_tag} = $data[21] if($data[21] ne "");  
    }  else {
        $munge = $to_tag."..".$from_tag;
        $all_best_hits{$munge}->{'from_Nterm'} = $data[8];  # gets overwritten if bidirectional.
        $all_best_hits{$munge}->{'from_Cterm'} = $data[9];  # 
        $all_best_hits{$munge}->{'to_Nterm'} = $data[6];    #
        $all_best_hits{$munge}->{'to_Cterm'} = $data[7];    #
	$all_best_hits{$munge}->{'to_score'} = $data[20];
        $jaccard{$from_tag} = $data[21] if($data[21] ne ""); 
    }
    # bidirectional best hit if non-zero 'to_score' and 'from_score'
}

# alright, all hits are loaded.  Now loop through all bidirectionals.

local (%cogs) = ();
local ($next_cog) = 1;
local (%prot2cog) = ();

foreach $bidirect (keys(%all_best_hits))        {
  if ($all_best_hits{$bidirect}->{'to_score'} && $all_best_hits{$bidirect}->{'from_score'})     {
      ($front, $back) = split (/\.\./,$bidirect);
      ($found1, $found2) = &search_cogs(\%cogs, \%prot2cog, $front, $back);
    if ($found2 > $found1) {       # make sure $found2 < $found1
        ($found1, $found2) = ($found2, $found1);
        ($front, $back) = ($back, $front);
    }
    $DEBUG && print ("Found#1 = $found1, Found#2 = $found2\n");

    if ($found1 == -1 && $found2 == -1)      {
        $DEBUG && print("Starting new COG, members $front and $back\n");
        &new_cog(\%cogs, \%prot2cog, $next_cog++, $front, $back);
    }    
    elsif ($found2 == -1)                    {
        #add $back to COG that contains $front
        $DEBUG && print("WILL ATTEMPT TO ADD  $back TO COG $found1\n");
        &add_to_cog(\%cogs, \%prot2cog, $found1, $back);
    }
    elsif ($found1 == $found2)               {
        $DEBUG && print("CONGRATULATIONS: a confirmation of cluster $found1\n");
        $cogs{$found1}->{'connectivity'} += 1;
    }
    else                                     {
        $DEBUG && print("MERGING current cogs $found1 and $found2\n");
        &merge_cogs(\%cogs, \%prot2cog, $found1, $found2);
    }
  }
}

$DEBUG && print ("============================\n");

foreach $key (keys %cogs)      {
    $n =  $cogs{$key}->{'size'};
    $perfect = $n * ($n - 1) / 2;
    # JC: note that the reported size does NOT include any Jaccard members that get printed below.
    # Don't want to change this without first checking that no downstream analysis relies on the incorrect counts.
    print("COG = $key, size $n, connections = $cogs{$key}->{'connectivity'}, perfect = $perfect;\n");
    *return_array = $cogs{$key}->{'listPtr'};
    foreach $term (sort @return_array)    {
	if(exists $jaccard{$term}){
	    my @jmembers = split(/\,/,$jaccard{$term});
	    foreach my $member (@jmembers){
		print "\t$member\n";
	    }
	}
	else{
	    print "\t$term\n";
	}
    }

    # HERE, COULD TEST CLOSURE UNDER TRANSITIVITY

    # HERE, COULD [OPTIONALLY] PULL IN UNIDIRECTIONAL BEST HITS
}

# The bidirectional best hit is two proteins.
# Search both for membership in up to one cluster each.
sub search_cogs{
    local ($cogs, $p2c, $front, $back) = @_;
    my $found1 = $p2c->{$front} || -1;
    my $found2 = $p2c->{$back} || -1;
    if ($DEBUG) {
	my $current_cog_count = scalar(keys %$cogs);
	print ("Current COG count = $current_cog_count\n");
	print ("RESEARCH  $front, $back, $found1, $found2 \n");
    }
    return ($found1, $found2);
}

# The bidirectional best hit is two proteins.
# Create a new cluster with exactly those two members
sub new_cog {
    local ($cogs, $p2c, $key, $front, $back) = @_;
    local (@array) = ($front, $back);
    $cogs->{$key}->{'listPtr'} = \@array;
    $cogs->{$key}->{'size'} = 2;
    $cogs->{$key}->{'connectivity'} = 1;     
    $p2c->{$front} = $key;
    $p2c->{$back} = $key;
}

# The bidirectional best hit is two proteins.
# One was already a member of a cluster. Now add in the other.
sub add_to_cog {
    local ($cogs, $p2c, $key, $back) = @_; 
    local (*temp_array);
    *temp_array = $cogs->{$key}->{'listPtr'};
    push (@temp_array, $back);
    $cogs->{$key}->{'listPtr'} = \@temp_array;
    $cogs->{$key}->{'size'} += 1;
    $cogs->{$key}->{'connectivity'} += 1;  
    $p2c->{$back} = $key;
}

# The bidirectional best hit is two proteins.
# Since these were in two different clusters, those clusters are linked.
# Append all the members of the donor cluster to the keep cluster,
# then delete the donor cluster.
sub merge_cogs {
    local ($cogs, $p2c, $key1, $key2) = @_;
    local (*keep_array, *donor_array);
    $DEBUG && print ("MERGE WORKING! $key1, $key2 \n");
    *keep_array = $cogs->{$key1}->{'listPtr'};
    *donor_array = $cogs->{$key2}->{'listPtr'};
    push(@keep_array, @donor_array);
    $cogs->{$key1}->{'listPtr'} = \@keep_array;
    $cogs->{$key1}->{'size'}   += $cogs->{$key2}->{'size'};
    $cogs->{$key1}->{'connectivity'} += (1 +  $cogs->{$key2}->{'connectivity'});
    foreach my $p (@donor_array) {
	$p2c->{$p} = $key1;
    }
    delete $cogs->{$key2}->{'listPtr'};
    delete $cogs->{$key2}->{'size'};
    delete $cogs->{$key2}->{'connectivity'};
    delete $cogs->{$key2};
}









