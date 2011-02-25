#!/usr/bin/perl

use strict;
#use warnings;
use Carp;
#use lib ($ENV{EUK_MODULES});
use lib "/home/sdaugherty/bin/legacy_bin";
use Overlap_piler;
use Getopt::Long qw(:config no_ignore_case bundling);
use Data::Dumper;
BEGIN {
        use Ergatis::Logger;
}

my ($mummer_coords_file, $mummer_delta_file, $strain, $molnum, $output_dir, $help);

my $min_match_length = 0;
my $max_overlap_bases = 100; #default
my $max_dist_between_segs = 1000;
my $method = "";


&GetOptions ('mummer_coords_file=s' => \$mummer_coords_file,
             'mummer_delta_file=s' => \$mummer_delta_file,
             'min_match_length=i' => \$min_match_length,
             'max_overlap_bases=i' => \$max_overlap_bases,
             'max_dist_between_segs=i' => \$max_dist_between_segs,
             'method=s' => \$method,
	     'strain=s' => \$strain,
             'pseudonum=i' => \$molnum,
	     'output_dir=s' => \$output_dir,
             'h' => \$help);

my $usage = <<_EOUSAGE_;

############################################################################################
#
# required:
#
# --mummer_coords_file     nucmer coords-formatted file 'use show-coords -T out.delta'
#
# --method                 nucmer | promer
#
# --strain		   name of the strain
#
# --pseudonum		   number of the pseudomolecule to be geberated
# optional:
#  
# --min_match_length       default zero
# --mummer_delta_file      if specified, only the relevant subset will be reported.
# --max_overlap_bases      number of overlap allowed along the reference sequence and matching contigs (default 100)
# --max_dist_between_segs  maximum length between two chained segments along query or reference   (default 1000);
#
#
#############################################################################################


_EOUSAGE_

    ;

## Getting the log file
my $logfile = Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                'LOG_LEVEL'=>4);
$logger = $logger->get_logger();


if  ((! $mummer_coords_file) || $help ) {
    die $usage;
}

if ($mummer_delta_file && ! -s $mummer_delta_file) {
	$logger->logdie("Error, cannot locate $mummer_delta_file");
}

unless ($method =~ /^(nucmer|promer)$/i) {
	$logger->logdie("Error, don't understand method $method");
}

my %QUERY_LENGTHS;

my $alignments = $output_dir."/alignments.txt";
open (my $log_fh, "> $alignments") || $logger->logdie("Could not open $alignments file for writing.$!");

main: {
    
    my @matches = &parse_match_file($mummer_coords_file);
    
    &build_DP_trellis(@matches);

    my $highest_scoring_node = &find_highest_scoring_node(@matches);

    ## trace back from highest scoring node. (node and match are synonymous)
    my $trace_node = $highest_scoring_node;
    
    my @traced_nodes;
    while ($trace_node) {
        push (@traced_nodes, $trace_node);
        $trace_node = $trace_node->{prev_best_node};
    }
    
    # report traced nodes:
    if ($mummer_delta_file) {
        &retrieve_subset_delta_file(@traced_nodes);
    }
    else {
        foreach my $node (reverse @traced_nodes) {
            print &get_match_text($node);
        }
    }
    exit(0);
}




####
sub parse_match_file {
    my $match_file = shift;

    my @matches;

    open (my $fh, $match_file) or $logger->logdie("Error, cannot open file $match_file\n");
    while (<$fh>) {
        unless (/^\d+/) { next; }
        chomp;
        my @x = split (/\t/);
        
        my ($ref_end5, $ref_end3, $query_end5, $query_end3, $ref_match_len, $query_match_len, $per_id, $ref_acc, $query_acc) = @x;

        if ($method =~ /promer/) {
            # output is a bit different...
            # set the ref and query accessions as the last two fields.
            
            $query_acc = pop @x;
            $ref_acc = pop @x;
        }
                

        my $ref_orient = ($ref_end5 < $ref_end3) ? '+' : '-';
        my $query_orient = ($query_end5 < $query_end3) ? '+' : '-';
        
        my ($ref_lend, $ref_rend) = sort {$a<=>$b} ($ref_end5, $ref_end3);
        my ($query_lend, $query_rend) = sort {$a<=>$b} ($query_end5, $query_end3);

        if ( (! exists $QUERY_LENGTHS{$query_acc}) || ($query_rend > $QUERY_LENGTHS{$query_acc}) ) {
            $QUERY_LENGTHS{$query_acc} = $query_rend;
        }
        
        unless ($ref_match_len >= $min_match_length) { next; }

        my $match_struct = {  
            ref_lend => $ref_lend,
            ref_rend => $ref_rend,
            ref_orient => $ref_orient,
            ref_acc => $ref_acc,
            
            query_lend => $query_lend,
            query_rend => $query_rend,
            query_orient => $query_orient,
            query_acc => $query_acc,
            
            per_id => $per_id,

            match_length => $ref_match_len,
            
            ## attributes for DP trace
            sum_novel_query_bp => $query_match_len,
            prev_best_node => undef, # link to previous best node in graph for DP tracing.
            query_matching_regions_thus_far => { $query_acc => [ 
                                                                 [$query_lend, $query_rend] 
                                                                 ] 
                                                             }, # list of all nonoverlapping regions of query matched up to this point in the graph, indexed by query accession.
                                                                 
                                                                 
           };
        
        push (@matches, $match_struct);
    }
    
    return (@matches);
    
}


####
sub build_DP_trellis {
    my @matches = @_;

    # sort matches by reference coordinates
    @matches = sort {$a->{ref_lend} <=> $b->{ref_lend}} @matches;
    
    # compare each node to the previous node
    for (my $i = 1; $i <= $#matches; $i++) {
        
        my $node_i = $matches[$i];
        my $node_i_novel_bp = $node_i->{sum_novel_query_bp};
        my $node_i_query_matching_regions_thus_far = $node_i->{query_matching_regions_thus_far};
        
        
        # node j comes before node i in the graph.
        for (my $j = $i - 1; $j >= 0; $j--) {
                        
            my $node_j = $matches[$j];
            my $node_j_novel_bp = $node_j->{sum_novel_query_bp};
            my $node_j_query_matching_regions_thus_far = $node_j->{query_matching_regions_thus_far};

            unless (&ref_coords_overlap($node_i, $node_j) <= $max_overlap_bases) { 
                # comparison not allowed.  
                next; 
            }
            
            my %join_i_j_novel_coords = &join_i_j_novel($node_i_query_matching_regions_thus_far,
                                                        $node_j_query_matching_regions_thus_far);
            
            my $sum_bases = &sum_regions(%join_i_j_novel_coords);

            if ($sum_bases > $node_i_novel_bp) {
                ## got a new higher score
                $node_i_novel_bp = $node_i->{sum_novel_query_bp} = $sum_bases;
                $node_i->{prev_best_node} = $node_j;
                $node_i->{query_matching_regions_thus_far} = {%join_i_j_novel_coords};
            }
        }
    }
}


####
sub join_i_j_novel {
    my ($regions_B_href, $regions_A_href) = @_;

    # A may have many query_accs, B should have only one providing it's single component alignment
    
    my @query_accs = keys %$regions_B_href;
    
    if (scalar @query_accs != 1) {
        confess "Error, regions_B_href should have only one query alignment.";
    }
    
    my $query_acc = $query_accs[0];
    my $matches_B_aref = $regions_B_href->{$query_acc};

    my $matches_A_aref = $regions_A_href->{$query_acc};
    unless ($matches_A_aref) {
        # region_B contributes completely novel sequence here.
        my %combined = (%$regions_A_href, %$regions_B_href);
        return (%combined);
    }
    
    ## join them into non-redundant overlapping regions
    my @combined_coords = &Overlap_piler::simple_coordsets_collapser(@$matches_A_aref, @$matches_B_aref);
    
    my %combined = %$regions_A_href;
    $combined{$query_acc} = [@combined_coords];
    
    return (%combined);
}

####
sub sum_regions {
    my %regions = @_;
    
    my $sum_length = 0;
    
    foreach my $query_acc (keys %regions) {
        foreach my $match (@{$regions{$query_acc}}) {
            my ($lend, $rend) = @$match;
            my $len = abs ($rend - $lend) + 1;
            $sum_length += $len;
        }
    }
    return ($sum_length);
    
}

####
sub find_highest_scoring_node {
    my @matches = @_;

    my $best_scoring_node = shift @matches;
    foreach my $match (@matches) {
        if ($match->{sum_novel_query_bp} > $best_scoring_node->{sum_novel_query_bp}) {
            $best_scoring_node = $match;
        }
    }

    return ($best_scoring_node);
}

####
sub get_match_text {
    my ($match) = @_;

    my $ref_lend = $match->{ref_lend};
    my $ref_rend = $match->{ref_rend};
    my $ref_orient = $match->{ref_orient};
    my $ref_acc = $match->{ref_acc};

    my $query_lend = $match->{query_lend};
    my $query_rend = $match->{query_rend};
    my $query_orient = $match->{query_orient};
    my $query_acc = $match->{query_acc};

    my $per_id = $match->{per_id};
    my $match_length = $match->{match_length};

    my ($ref_end5, $ref_end3) = ($ref_orient eq '+') ? ($ref_lend, $ref_rend) : ($ref_rend, $ref_lend);
    my ($query_end5, $query_end3) = ($query_orient eq '+') ? ($query_lend, $query_rend) : ($query_rend, $query_lend);

    return ("$ref_acc\t$ref_end5\t$ref_end3\t$query_acc\t$query_end5\t$query_end3\t$per_id\t$match_length\n");
    
}


####
sub retrieve_subset_delta_file {
    my @traced_nodes = @_;

    my %alignment_tokens; #just those alignments to print
    foreach my $traced_node (@traced_nodes) {
        my $match_text = &get_match_text($traced_node);
        chomp $match_text;
        my ($ref_acc, $ref_end5, $ref_end3, $query_acc, $query_end5, $query_end3, $per_id, $match_len) = split(/\t/, $match_text);

        my $token = join ("__", $ref_acc, $ref_end5, $ref_end3, $query_acc, $query_end5, $query_end3);
        $alignment_tokens{$token} = 1;
    }

    ## parse the delta file:
    my $curr_header;
    my $curr_text = "";
    my ($ref_acc, $query_acc);

    open (my $fh, $mummer_delta_file) || $logger->logdie("Error could not open $mummer_delta_file file for reading");
    # report the mummer header info (top 2 lines)
    my $line = <$fh>;
#    print $line;
    $line = <$fh>;
#    print $line;
    
    my @entries; ## store delta alignment text so we can sort it according to reference coordinates:
           
    my $ref_end5_stored;
    my $query_acc_stored;

    while (<$fh>) {
        if (/^>/) {
            if ($curr_text) {
                # print $curr_text;
                push (@entries, { text => $curr_text,
                                  coord => $ref_end5_stored,
                                  query_acc => $query_acc_stored } );
            }
            $curr_text = "";
            $curr_header = $_;
            my @rest;
            ($ref_acc, $query_acc, @rest) = split (/\s+/);
            $ref_acc =~ s/>//;
        }
        else {
            my @x = split (/\s+/);
            if (scalar (@x) > 1) {
                my ($ref_end5, $ref_end3, $query_end5, $query_end3, @rest) = @x;
                
                my $token = join ("__", $ref_acc, $ref_end5, $ref_end3, $query_acc, $query_end5, $query_end3);
                if ($alignment_tokens{$token}) {
                    # got one. record it.
                    unless ($curr_text) {
                        # add the header
                        $curr_text = $curr_header;
                        $ref_end5_stored = $ref_end5;
                        $query_acc_stored = $query_acc;
                    }
                    $curr_text .= $_;
                    while (my $line = <$fh>) {
                        $curr_text .= $line;
                        my $value = $line;
                        chomp $value;
                        if ($value == 0) {
                            last;
                        }
                    }
                }
            }
        }
    }
    if ($curr_text) {
        push (@entries, { text => $curr_text,
                          coord => $ref_end5_stored,
                          query_acc => $query_acc_stored } );
    }
    
    
    # print ordered list of query accessions for prettier nucmer plot.
    my $outfile = $output_dir."/".$strain."_".$molnum.".map"; 
    open ($fh, "> $outfile") or $logger->logdie("Error could not open $outfile file for writing.");
   
    ## sort according to coordinate and print delta text:
    
    # sort based on an average coordinate weighting scheme.
    
    foreach my $entry (@entries) {
        my $text = $entry->{text};
        
        my $sum_lengths = 0;
        my $sum_match_pos = 0;
        my @lines = split (/\n/, $text);
        shift @lines; # rid the header;
        
        my %orientation_voter = ( '+' => 0, '-' => 0);
        
        my @component_matches;

        foreach my $line (@lines) {
            chomp $line;
            my @x = split (/\s+/, $line);
            if (scalar @x > 1) {
                my ($ref_lend, $ref_rend, $query_end5, $query_end3) = ($x[0], $x[1], $x[2], $x[3]);
                my $ref_length = abs ($ref_rend - $ref_lend) + 1;
                my $midpt = ($ref_lend + $ref_rend) / 2;
                $sum_lengths += $ref_length;
                $sum_match_pos += $ref_length * $midpt;
                
                my $match_orient = ($query_end5 < $query_end3) ? '+' : '-';
                $orientation_voter{$match_orient} += $ref_length;
                
                my ($query_lend, $query_rend) = sort {$a<=>$b} ($query_end5, $query_end3);
                
                push (@component_matches, {  ref_lend => $ref_lend,
                                             ref_rend => $ref_rend,
                                             query_lend => $query_lend,
                                             query_rend => $query_rend, 
                                             orient => $match_orient,
                                         } );
                
            }
        }
        
        if ($sum_lengths == 0) { 
            $logger->logdie("$text\nError, no sum lengths.\n");
        }
        
        # use weighted average for coordinate positioning:
        my $coord = $sum_match_pos / $sum_lengths;
        $entry->{coord} = $coord;

        ## set the orientation based on majority vote:
        my $contig_orient = ($orientation_voter{'+'} > $orientation_voter{'-'}) ? '+' : '-';
        
        $entry->{orient} = $contig_orient;
        
        $entry->{component_matches} = [@component_matches];
        
        my $CHAIN_CONTIG_SEGS = 1;
        if ($CHAIN_CONTIG_SEGS) {
            &assign_entry_coordinate_by_longest_chain_midpt($entry);
        }
        

    }
    

  

    ## sort by coordinate
    @entries = sort {$a->{coord}<=>$b->{coord}} @entries;
    foreach my $entry (@entries) {
        
#        print $entry->{text};
        
        my $query_acc = $entry->{query_acc};
        my $query_orient = $entry->{orient};

#        print $fh "$query_acc\t" . $QUERY_LENGTHS{$query_acc} . "\t$query_orient\n";
    	print $fh "$strain.pseudomolecule.$molnum\t" . "$query_acc\t" . "$query_orient\n";
    }
    close $fh;
    
    return;
}

####
sub ref_coords_overlap {
    my ($match_A, $match_B) = @_;

    my ($match_A_lend, $match_A_rend) = ($match_A->{ref_lend}, $match_A->{ref_rend});
    my ($match_B_lend, $match_B_rend) = ($match_B->{ref_lend}, $match_B->{ref_rend});

    # if they do not overlap, then overlap = 0
    
    unless ($match_A_lend <= $match_B_rend && $match_A_rend >= $match_B_lend) {
        return (0);
    }

    return (&nucs_in_common($match_A_lend, $match_A_rend, $match_B_lend, $match_B_rend));
}

    
####
sub nucs_in_common {
    my ($e5, $e3, $g5, $g3) = @_;
    ($e5, $e3) = sort {$a<=>$b} ($e5, $e3);
    ($g5, $g3) = sort {$a<=>$b} ($g5, $g3);
    my $length = abs ($e3 - $e5) + 1;
    my $diff1 = ($e3 - $g3);
    $diff1 = ($diff1 > 0) ? $diff1 : 0;
    my $diff2 = ($g5 - $e5);
    $diff2 = ($diff2 > 0) ? $diff2 : 0;
    my $overlap_length = $length - $diff1 - $diff2;
    return ($overlap_length);
}


####
sub assign_entry_coordinate_by_longest_chain_midpt {
    my ($entry) = @_;

    my $query_acc = $entry->{query_acc};
    
    my @component_matches = @{$entry->{component_matches}};
    
    ## sort by reference coordinate position.
    @component_matches = sort {$a->{ref_lend}<=>$b->{ref_lend}} @component_matches;
    
    # init score and prev for DP
    foreach my $match (@component_matches) {
        my $length = $match->{query_rend} - $match->{query_lend} + 1;
        $match->{score} = $length;
        $match->{sum_score} = $length;
        $match->{prev} = undef;
    }
     
    for (my $i = 1; $i <= $#component_matches; $i++) {
        
        my $match_i = $component_matches[$i];
        my $orient_i = $match_i->{orient};
        my ($i_lend, $i_rend) = ($match_i->{query_lend}, $match_i->{query_rend});
        my ($i_ref_lend, $i_ref_rend) = ($match_i->{ref_lend}, $match_i->{ref_rend});
        
        # score values.
        my $i_score = $match_i->{score};
        my $i_sum_score = $match_i->{sum_score};
        
        
        for (my $j = $i - 1; $j >= 0; $j--) {
            my $match_j = $component_matches[$j];
            my $orient_j = $match_j->{orient};
            my ($j_lend, $j_rend) = ($match_j->{query_lend}, $match_j->{query_rend});
            my ($j_ref_lend, $j_ref_rend) = ($match_j->{ref_lend}, $match_j->{ref_rend});
            
            unless ($orient_i eq $orient_j) { next; } 
            
            if ($orient_j eq '-') {
                ## swap i and j query coordinates  
                ($i_lend, $i_rend, $j_lend, $j_rend) = ($j_lend, $j_rend, $i_lend, $i_rend);
            }
            
            unless ($j_rend < $i_lend + $max_overlap_bases) { next; } # don't allow more than that overlap along query
            unless ($j_ref_rend < $i_ref_lend + $max_overlap_bases) { next; } # ditto along the reference.

            my $dist_between_segs = $i_lend - $j_rend;
            if ($dist_between_segs > $max_dist_between_segs) { next; } # intermatch length too high!
            
            my $ref_dist_between_segs = $i_ref_lend - $j_ref_rend;  # i after j
            if ($ref_dist_between_segs > $max_dist_between_segs) { next; }
            
            my $j_sum_score = $match_j->{sum_score};

            if ($i_score + $j_sum_score > $i_sum_score) {
                ## got new best score
                $i_sum_score = $match_i->{sum_score} = $i_score + $j_sum_score;
                $match_i->{prev} = $match_j;
            }
            
        }
    }
    
    ## get longest chain of matches;
    my $highest_scoring_match = shift @component_matches;
    foreach my $other_match (@component_matches) {
        if ($other_match->{sum_score} > $highest_scoring_match->{sum_score}) {
            $highest_scoring_match = $other_match;
        }
    }

    my @best_chain;
    while ($highest_scoring_match) {
        push (@best_chain, $highest_scoring_match);
        $highest_scoring_match = $highest_scoring_match->{prev};
    }

    ## determine midpt of chain:
    my @coords;
    
    print $log_fh "// Alignment Chain\n";
    foreach my $match (@best_chain) {
        
        print $log_fh "$query_acc\t" . $match->{ref_lend} . "-" . $match->{ref_rend} 
        . "\t" . $match->{query_lend} . "-" . $match->{query_rend} 
        . "\t" . $match->{orient} . "\n";
        
        push (@coords, $match->{ref_lend}, $match->{ref_rend});
    }
    @coords = sort {$a<=>$b} @coords;

    my $lend = shift @coords;
    my $rend = pop @coords;

    my $midpt = ($rend + $lend)/2;

    $entry->{coord} = $midpt;
    
    print $log_fh "$query_acc\tMID: $midpt\n\n";
    
    # print STDERR Dumper (\@best_chain);
    
    return;
}


__END__

