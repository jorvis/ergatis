#!/usr/local/bin/perl


use lib("/usr/local/annotation/PNEUMO/clu_dir/BSML/ANNOTATION/bsml/src");
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Log::Log4perl qw(get_logger);
use BsmlReader;
use BsmlParserTwig;


my %options = ();
my $results = GetOptions (\%options, 'match_asmbl_id|b=s', 'asmbl_id|a=s', 'output_dir|f=s',
                                     'bsml_dir|d=s', 'help|h', 'order|o=s' );

Log::Log4perl->init("/usr/local/annotation/PNEUMO/clu_dir/papyrus/src/perl/log.conf");
my $logger = get_logger();
###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $asmbl_id       = $options{'asmbl_id'};
my $match_asmbl_id = $options{'match_asmbl_id'};
my $BSML_dir       = $options{'bsml_dir'};
$BSML_dir =~ s/\/+$//;       #remove terminating '/'s
my $order = $options{'order'} || '1';


if(!$asmbl_id or !$match_asmbl_id or !$BSML_dir or exists($options{'help'})) {
    &print_usage();
}

###-------------------------------------------------------###


my $bsml_file = "$BSML_dir/asmbl_${asmbl_id}_blastp.bsml";
if(! -d $BSML_dir ) {
    $logger->fatal("The directory \"$BSML_dir\" cannot be found.  Exiting...");
    exit 10;
}elsif(! -s $bsml_file) {
    $logger->fatal("The file \"$bsml_file\" cannot be found. Exiting...");
    exit 5;
}

my $parser = new BsmlParserTwig;
my $reader = BsmlReader->new();

$parser->parse( \$reader, $bsml_file );
my @orders = split(/,/, $order);          #parses the order input, as comma separated argument        
@orders = grep {/\d+/} @orders;           #weed out non digit argument
@orders = sort { $a <=> $b } @orders;     #sort in descending order

$logger->debug("The defined orders are \"@orders\"");

#my $gene_asmbl_id = get_gene_asmbl_id();  #hash_ref stores each gene and the asmbl_id to which it belongs
my $bit_score_comparison = fetch_score_for_top_n_hit($asmbl_id, $match_asmbl_id, $orders[-1]);
$bit_score_comparison = calculate_percent_bit_score($bit_score_comparison);

foreach (@orders) {
    my $output_ref = generate_bit_score_report($bit_score_comparison, $_);
    if($output_ref) {
	$logger->info("bit_score_output defined for order $_");
	my @sorted_lines = sort { $b->[4] <=> $a->[4] } @$output_ref;     #sort based on raw bit score in descending order
	foreach (@sorted_lines) {
	    my $line = join("\t",@$_);
	    print "$line\n";
	}
    }
} 


sub generate_bit_score_report {
#generates a list_ref containing output for a defined order

    my $hash_ref = shift;
    my $order = shift || '1';

    $logger->info("Begin generating bit_score_report for order $order");
    
    my $avg_bit_score = calculate_avg_bit_score($hash_ref, $order);
    if(defined($avg_bit_score)) {
	$avg_bit_score = sprintf('%.3f', $avg_bit_score);
	$logger->debug("avg_percent_bit_score for order $order is $avg_bit_score");
    }
    my $array_ref;
    foreach my $query_gene(keys %{ $hash_ref }) {
	my @list = grep {$_ ne $query_gene} keys(%{$hash_ref->{$query_gene}});  #remove gene vs itself from the presorted list
	@list = sort { $hash_ref->{$query_gene}->{$b}->{'bit_score'} <=> $hash_ref->{$query_gene}->{$a}->{'bit_score'} } @list;
	my $match_gene = $list[$order-1];
	if(defined($match_gene)) {
	    my $bit_score = $hash_ref->{$query_gene}->{$match_gene}->{'bit_score'}; 
	    my $base_score = $hash_ref->{$query_gene}->{$query_gene}->{'bit_score'}; 
            my $percent_score = $hash_ref->{$query_gene}->{$match_gene}->{'percent_bit_score'};
            $query_gene =~ s/\_aa$//;
	    $match_gene =~ s/\_aa$//;	    
	    push(@$array_ref, [$query_gene, $match_gene, $bit_score, $base_score, $percent_score, $order, $avg_bit_score]);
        }				

    }

    $logger->info("Finished generating bit_score_report for order $order");
    
    return $array_ref;

}


sub calculate_avg_bit_score {
#calculates the average percent bit score of a given order

    my $hash_ref = shift;
    my $order = shift || '1';

    $logger->info("Begin calculating avg percent bit scores for order $order");

    my $total_bit_score=0;
    my $number_of_scores=0;
    foreach my $query_gene(keys %{ $hash_ref }) {
	my @list = grep {$_ ne $query_gene} keys(%{$hash_ref->{$query_gene}});  #remove gene vs itself from the presorted list
        @list = sort { $hash_ref->{$query_gene}->{$b}->{'bit_score'} <=> $hash_ref->{$query_gene}->{$a}->{'bit_score'} } @list;	
	my $index = $order - 1;
	if(defined($list[$index])) {
	    $total_bit_score += $hash_ref->{$query_gene}->{$list[$index]}->{'percent_bit_score'};
	    $number_of_scores++;
        }
    }

    $logger->info("Finished calculating avg percent bit scores for order $order");
    if(!$total_bit_score or !$number_of_scores) {
	$logger->info("Unable to calcuate avg_percent_bit_score.  Returning undef");
	return undef;
    } else {
	return $total_bit_score/$number_of_scores;
    }
	
}

     
sub calculate_percent_bit_score {

    my $hash_ref = shift;

    $logger->info("Begin calculating percent bit scores");

    foreach my $query_gene(keys %{ $hash_ref }) {
	if(!exists($hash_ref->{$query_gene}->{$query_gene})) {
	    #print STDERR "NO BASE SCORE for $query_gene\n";
            $logger->error("Can't Find base score for $query_gene!!!  Skipping calculation");
	    next;
        }
	my $base_score  = $hash_ref->{$query_gene}->{$query_gene}->{'bit_score'};       #the bit score of gene vs itself
        foreach my $match_gene (keys %{ $hash_ref->{$query_gene} }) {
	    my $match_score = $hash_ref->{$query_gene}->{$match_gene}->{'bit_score'};   #the bit score of a match
	    $hash_ref->{$query_gene}->{$match_gene}->{'percent_bit_score'}= sprintf('%.3f', $match_score/$base_score); # %bit_score
        }
    }

    $logger->info("Finished calculating percent bit score");
    return $hash_ref;
}



sub fetch_score_for_top_n_hit {

    my $query_asmbl_id = shift;
    my $match_asmbl_id = shift;
    my $order = shift || '1';

    my $bit_score_hash = {};
    my $hash_ref = $reader->fetchAlignmentScoresBetweenAssemblies("PNEUMO_${query_asmbl_id}", "PNEUMO_${match_asmbl_id}");

    foreach my $query_gene(keys %{ $hash_ref }) {
	my @list = grep {$_ ne $query_gene} keys(%{$hash_ref->{$query_gene}});  #remove gene vs itself from the presorted list
	@list = sort { $hash_ref->{$query_gene}->{$b}->{'bit_score'} <=> $hash_ref->{$query_gene}->{$a}->{'bit_score'} } @list;
	my $current_order = 0;
	foreach my $m_gene (@list) {
	    if($current_order < $order) {
		$bit_score_hash->{$query_gene}->{$m_gene} = $hash_ref->{$query_gene}->{$m_gene};
		$current_order++;
            }
        }
	$bit_score_hash->{$query_gene}->{$query_gene} = $hash_ref->{$query_gene}->{$query_gene};
    }

    return $bit_score_hash;

}



sub print_usage {


    print STDERR "SAMPLE USAGE:  bit_score_analysis.pl -d bsml_dir -a 1, -b 2 -o 1,2,3\n";
    print STDERR "  --bsml_dir  = location of BSML documents\n";
    print STDERR "  --asmbl_id  = query genome\n";
    print STDERR "  --match_asmbl_id = match_genome\n";
    print STDERR "  *optional: --order default= 1 The top n hit to calculate\n";
    print STDERR "  --help = This help message.\n";
    exit 1;

}




