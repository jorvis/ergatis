package NGS::BLAT;

use strict;
use warnings;
use Data::Dumper;
use Bio::SearchIO;

sub new {
	my ($class, $filename) = @_;
	
	my $self = bless {}, ref($class) || $class;

	my $psl = Bio::SearchIO->new( -file => $filename, 
								  -format => 'psl' );
	
	$self->{psl} = $psl;
	return $self;
}

sub get_gff3_features {
	my ($self) = @_;
	my %features;

	my $id = 0;
	
	while( my $results = $self->{psl}->next_result ) {


		while( my $hit = $results->next_hit ) {

			#print "NAME: " . $hit->name . ", " . $hit->length . ", " . $hit->raw_score . ", " . $hit->rank . "\n";

			while( my $hsp = $hit->next_hsp ) {
				
				#print Dumper $hsp;
				
				#print "\t" . $hsp->length('total') . ", " . $hsp->percent_identity . ", " . $hsp->frac_identical . ", " . $hsp->strand . "\n";
				#my $blocksizes = $hsp->gap_blocks('query');
				#for(my $i=0; $i<@$blocksizes; $i++) {
				#	print "I: $i --> $blocksizes->[$i][0], $blocksizes->[$i][1]\n";
				#}
				my @features_old;
				push @features_old, { 'match' => $hsp->num_identical, 
									  'mismatches' => $hsp->mismatches, 
									  'rep_matches' => "unkn", 
									  'orientation' => $hsp->strand, 
									  'source' => $results->algorithm,
									  'query_id' => $results->query_name,
									  'query_size' => $results->query_length,
									  'query_start' => $hsp->start('query'),
									  'query_end' => $hsp->end('query'),
									  
									  'hit_id' => $hit->name,
									  'hit_size' => $hit->length, 
									  'hit_start' => $hsp->start('hit'),
									  'hit_end' => $hsp->end('hit'),
									  
									  'blocksizes' => "unkn",
									  'qstarts' => "unkn", 
									  'hstarts' => "unkn"
										  
									  };

				#my $unique_hit_id = $hit->name . "_" . ++$id;
				my $unique_hit_id = $hit->name;
				
			  # example. http://iubio.bio.indiana.edu/gmod/tandy/blat2gff.pl
				push @{ $features{$results->query_name} }, ( 
															 $unique_hit_id, $results->algorithm, "match_part", 
															 $hit->length, $hit->start('hit'), $hit->end('hit'), 
															 $hsp->score, $hsp->strand, ".", $results->query_name, 
															 $results->query_name . " " . $hsp->start('query') . " " . $hsp->end('query') 
															 );
			}
		}
	}

	foreach my $feat (keys %features) {
		print "\nFEAT: $feat\n----------------------------\n";
		foreach my $el (@{ $features{$feat} }) {
			print "\t$el\n";
		}
		exit;
	}
	
	#return \@features;

}

1;
