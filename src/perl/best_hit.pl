#!/local/perl/bin/perl
# this program accepts a stream of skim results of the best
# hits among a number of genomes and generates clusters of
# proteins connected by best hits. The core of each group is
# a set of proteins connected by BIDIRECTIONAL BEST HITS.
# In the special case of two genomes, it finds exactly all
# pairs of bidirectional best hits.

# PROGRAMMING STRATEGY:
# Load all best hits into hashes. Run through a first time and
# calculate "cog pieces" and build a hash of them.  Merge cog pieces
# as necessary.  Optionally, in a final run, attach "friends".


# This version being redone as "structures"

require "getopts.pl";
$usage = <<_MESSAGE_;
Usage:  fish.pl feat_name[s]
cat source_of_hit_results  |  bi-best-hitter.pl [options] > output_file
    -D     DEBUG              default = off        Show diagnostic messages
    -h     help               print this message
    -c     cutoff             default = 0.0001
    -s     secondary cutoff   see -u option
    -u     unidirectional     add unidirectional hits better than secondary cutoff   
    -n     minimum            default = 2      minimum set of transitively related.
_MESSAGE_

&Getopts('Dhuc:s:n:');
die "$usage\n" if ($opt_h == 1);

$DEBUG = $opt_D;

$cutoff = (length($opt_c) == 0) ? 0.0001 : $opt_c;

$quorum = ($opt_n < 2) ? 2 : $opt_n;  

$DEBUG && print ("DEBUG is ON, cog minimum is $quorum, cutoff is $cutoff\n");

# first load all the hits data
while ($btab = <STDIN>) {
    @data = split (/\t/,$btab);
    if ($data[20] > $cutoff)  {next;}    # discard if a weak hit by E-value.        
    if ($data[20] < 0)        {next;}    # a little quality control
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
    }  else {
        $munge = $to_tag."..".$from_tag;
        $all_best_hits{$munge}->{'from_Nterm'} = $data[8];  # gets overwritten if bidirectional.
        $all_best_hits{$munge}->{'from_Cterm'} = $data[9];  # 
        $all_best_hits{$munge}->{'to_Nterm'} = $data[6];    #
        $all_best_hits{$munge}->{'to_Cterm'} = $data[7];    #
        $all_best_hits{$munge}->{'to_score'} = $data[20]; 
    }
    # bidirectional best hit if non-zero 'to_score' and 'from_score'
}

# alright, all hits are loaded.  Now loop through all bidirectionals.

local (%cogs) = ();
local ($next_cog) = 1;


foreach $bidirect (keys(%all_best_hits))        {
  if ($all_best_hits{$bidirect}->{'to_score'} && $all_best_hits{$bidirect}->{'from_score'})     {
      ($front, $back) = split (/\.\./,$bidirect);
      ($found1, $found2) = &search_cogs(\%cogs, $front, $back);
    if ($found2 > $found1) {       # make sure $found2 < $found1
        ($found1, $found2) = ($found2, $found1);
        ($front, $back) = ($back, $front);
    }
    $DEBUG && print ("Found#1 = $found1, Found#2 = $found2\n");

    if ($found1 == -1 && $found2 == -1)      {
        $DEBUG && print("Starting new COG, members $front and $back\n");
        &new_cog(\%cogs, $next_cog++, $front, $back);
    }    
    elsif ($found2 == -1)                    {
        #add $back to COG that contains $front
        $DEBUG && print("WILL ATTEMPT TO ADD  $back TO COG $found1\n");
        &add_to_cog(\%cogs, $found1, $back);
    }
    elsif ($found1 == $found2)               {
        $DEBUG && print("CONGRATULATIONS: a confirmation of cluster $found1\n");
        $cogs{$found1}->{'connectivity'} += 1;

    }
    else                                     {
        $DEBUG && print("MERGING current cogs $found1 and $found2\n");
        &merge_cogs(\%cogs, $found1, $found2);
    }
  }
}

$DEBUG && print ("============================\n");

foreach $key (keys %cogs)      {
    $n =  $cogs{$key}->{'size'};
    $perfect = $n * ($n - 1) / 2;
    print("COG = $key, size $n, connections = $cogs{$key}->{'connectivity'}, perfect = $perfect;\n");
    *return_array = $cogs{$key}->{'listPtr'};
    foreach $term (sort @return_array)    {
        print "\t$term\n";
    }

    # HERE, COULD TEST CLOSURE UNDER TRANSITIVITY

    # HERE, COULD [OPTIONALLY] PULL IN UNIDIRECTIONAL BEST HITS
}

# The bidirectional best hit is two proteins.
# Search both for membership in up to one cluster each.
sub search_cogs{
    local (\%cogs, $front, $back) = @_;
    local ($found1, $found2) = (-1, -1);
    local (*temp_array);
    $current_cog_count = 0;
    foreach $key (keys %cogs)                      {
        $current_cog_count++;
	*temp_array = $cogs{$key}->{'listPtr'};
        for ($k=0; $k < @temp_array; $k++)     {
	    $found1 = $key if ($temp_array[$k] eq $front);
            $found2 = $key if ($temp_array[$k] eq $back);
        }
    }
    $DEBUG && print ("Current COG count = $current_cog_count\n");
    $DEBUG && print ("RESEARCH  $front, $back, $found1, $found2 \n");
    return ($found1, $found2);
}

# The bidirectional best hit is two proteins.
# Create a new cluster with exactly those two members
sub new_cog {
    local (\%cogs, $key, $front, $back) = @_;
    local (@array) = ($front, $back);
    $cogs{$key}->{'listPtr'} = \@array;
    $cogs{$key}->{'size'} = 2;
    $cogs{$key}->{'connectivity'} = 1;     
}


# The bidirectional best hit is two proteins.
# One was already a member of a cluster. Now add in the other.
sub add_to_cog {
    local (\%cogs, $key, $back) = @_; 
    local (*temp_array);
    *temp_array = $cogs{$key}->{'listPtr'};
    push (@temp_array, $back);
    $cogs{$key}->{'listPtr'} = \@temp_array;
    $cogs{$key}->{'size'} += 1;
    $cogs{$key}->{'connectivity'} += 1;  
}


# The bidirectional best hit is two proteins.
# Since these were in two different clusters, those clusters are linked.
# Append all the members of the donor cluster to the keep cluster,
# then delete the donor cluster.
sub merge_cogs {
    local (\%cogs, $key1, $key2) = @_;
    local (*keep_array, *donor_array);
    $DEBUG && print ("MERGE WORKING! $key1, $key2 \n");
    *keep_array = $cogs{$key1}->{'listPtr'};
    *donor_array = $cogs{$key2}->{'listPtr'};
    push(@keep_array, @donor_array);
    $cogs{$key1}->{'listPtr'} = \@keep_array;
    $cogs{$key1}->{'size'}   += $cogs{$key2}->{'size'};
    $cogs{$key1}->{'connectivity'} += (1 +  $cogs{$key2}->{'connectivity'});
    delete $cogs{$key2}->{'listPtr'};
    delete $cogs{$key2}->{'size'};
    delete $cogs{$key2}->{'connectivity'};
    delete $cogs{$key2};
}









