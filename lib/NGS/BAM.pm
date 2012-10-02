package NGS::BAM;

use strict;
use warnings;
use Data::Dumper;
use Bio::DB::Sam;

my $library_type;

sub new {
	my ($class, %args) = @_;
	
	my $self = bless {}, ref($class) || $class;
	$self->_init(\%args);
	return $self;
}

sub _init {
	my ($self, $args) = @_;
	my $bam = Bio::DB::Sam->new( -bam => $args->{bam}, -expand_flags => $args->{ 'uniq' } );
	$self->{_bam} = $bam;
	$self->{_uniq} = $args->{ 'uniq' };
	if(defined $args->{total_mapped_reads}){
	    $self->{_total_mapped_reads} = $args->{total_mapped_reads};
	}
	if(defined $args->{stranded}){
	    $self->{_stranded} = $args->{stranded};
	}
	else{
	    $self->{_stranded} = 0;
	}
	
	if(defined $args->{library_type}){
	    $library_type = $args->{library_type};
	}
	else{
	    $library_type = "firststrand";
	}

	my $seq_ids = $self->seq_ids;

	my %sizes;
	foreach (@{$seq_ids}){
	    my $segment = $self->{_bam}->segment( -seq_id => $_);
	    $sizes{$_} = $segment->end;
	}
	$self->{sizes} = \%sizes;

    }

sub seq_ids(){
    my ($self) = @_;
    my @seqs = $self->{_bam}->seq_ids;
    return \@seqs;
}

sub sizes(){
    my ($self) = @_;
    return $self->{sizes};
}


sub get_total_mapped_reads {
    my ($self) = @_;

    return $self->{_total_mapped_reads} if(defined $self->{_total_mapped_reads});
	
    my $total_mapped_reads = 0;

    my @seq_ids = $self->{_bam}->seq_ids;
    
    foreach my $seq_id (@seq_ids) {
		my $segment = $self->{_bam}->segment( -seq_id => $seq_id);
		my $iterator = $segment->features(-iterator=>1);
		
		while (my $a = $iterator->next_seq){
			my $matched_regions = $a->get_tag_values( 'NH' );
			$total_mapped_reads += $self->{_uniq} ? ( 1/$matched_regions ) : 1;
		}
		
    }
    $self->{_total_mapped_reads} = $total_mapped_reads;
    return $self->{_total_mapped_reads};
}

#
# Options include strand (+/-), segment (seq_id, start, stop)
# Chromosome is required to be defined, and incoming coords are 1-based
# Spanning reads, ie spliced alignment regions, are counted as being covered
# Ex if a read has a CIGAR score of 31M164N3M , the 164bp spliced region will also show 1x coverage
# Maybe this should be changed?
#
sub get_per_bp_coverage {
    my ($self, $args) = @_;
    
    if( !defined($$args{seq_id}) ) {
		die "Must define a seq_id to retrieve coverage on a bp level!";
    }
    
    my $segment;
    my $iterator;
	
    if( defined($$args{start}) && defined($$args{end}) ) {
		$segment = $self->{_bam}->segment( -seq_id => $$args{seq_id}, 
										   -start => $$args{start}, 
										   -end => $$args{end} );

    } else {
		$segment = $self->{_bam}->segment(-seq_id => $$args{seq_id});
    }
    
    my @bp_coverage = ((0) x (($segment->end - $segment->start)+1)); #because segment stuff is always 1-based
	
    if( !defined ($$args{strand}) || ($$args{strand} eq ".") ) {
		$iterator = $segment->features(-iterator=>1);
		
    } elsif( $$args{strand} eq "+" ) {
		
		# this should get only alignments to the forward strand
		$iterator = $segment->features(-iterator=>1,
									   -filter => sub { my $a = shift;
														return get_alignment_strand($a) > 0;
													});
		
    } elsif( $$args{strand} eq "-" ) {
		
		# this should get only alignments to the reverse strand
		$iterator = $segment->features(-iterator=>1,
									   -filter => sub { my $a = shift;
														return get_alignment_strand($a) < 0;
													});
	    
    }
	
	#---------------------------------------------------------
	#  See if user gave start and end coords for chromosome;
	#  Otherwise, take entire chromosome into consideration.
	#                                      # Mahesh Vangala
	#---------------------------------------------------------
	my $s_start = $$args{ 'start' } || $segment->start; 
   	my $s_end = $$args{ 'end' } || $segment->end;
	#---------------------------------------------------------
	
    while( my $a = $iterator->next_seq ) {

		#-----------------------------------------------------------
        # We need to consider CIGAR string instead
        #                                       # Mahesh Vangala
        #-----------------------------------------------------------
		
		my $read_start = $a->start;
        my $read_end = $a->end;
        my $cigar = $a->cigar_str;

        while( $read_start <= $s_end && $cigar =~ /(\d+)(\w)/g ) {
			my( $num, $token ) = ( $1, $2 );
			if( $token eq 'M' || $token eq 'X' || $token eq '=' ) {
				foreach my $int( 1 .. $num ) {
					if( $read_start + $int >= $s_start && $read_start + $int <= $s_end ) {
						#Remember, segment is 1 based but coverage is 0 based!
						$bp_coverage[ $read_start - $s_start + $int - 1 ]++;
					}
				}
				$read_start += $num;
			}
			elsif( $token eq 'D' || $token eq 'N' ) {
				$read_start += $num;
			}
			elsif( $token eq 'I' ) {
				#
				# just go to next token
				#
			}
        }
		#--------------------------------------------------------------
    }
    
    return \@bp_coverage;
}

#
# counts the number of hits to a region, paired end mates separately
#
sub count_region_hits {
    my ($self, $region) = @_;  

    my $segment = $self->{_bam}->segment( -seq_id => $region->{seq_id}, -start => $region->{start}, -end => $region->{end});
	
    my $total_aligned_reads = 0; 
    my $region_length = $segment->end - $segment->start;
    
    my $region_iterator;

    if( !defined $region->{strand} ) {
	$region_iterator = $segment->features(-iterator=>1,
					      -seq_id => $region->{seq_id});		
    } 
    else {
	if( $region->{strand} eq "+" ) {
	    $region_iterator = $segment->features(-iterator=>1,
						  -seq_id => $region->{seq_id}, 
						      -filter => sub {
							  my $a = shift;
							  return get_alignment_strand($a) > 0;
						      });
	} elsif( $region->{strand} eq "-" ) {
	    $region_iterator = $segment->features(-iterator=>1,
						  -seq_id => $region->{seq_id}, 
						      -filter => sub {
							  my $a = shift;
							  return get_alignment_strand($a) < 0;
						      });
	} else {
	    $region_iterator = $segment->features(-iterator=>1,
						  -seq_id => $region->{seq_id});
	}
    }

    while (my $a = $region_iterator->next_seq){	
	$total_aligned_reads++;
    }

    return { 'name' => $region->{'name'}, 
	     'seq_id' => $region->{seq_id}, 
	     'hits' => $total_aligned_reads };
}

#
# Calculates the coverage percentage, avg depth for a region
# Requires the seq_id, start, and end, and optionally takes the strand
#
sub calculate_region_coverage {
    my ($self, $region) = @_;  
  
    my $bps_covered = 0;
    my $total_depth = 0;
    my $coverage;

    if( !defined $region->{strand} ) {
		$coverage = $self->get_per_bp_coverage({ seq_id => $region->{seq_id}, 
												 start => $region->{start}, 
												 end => $region->{end} });
		
    } else {
		$coverage = $self->get_per_bp_coverage({ seq_id => $region->{seq_id}, 
												 start => $region->{start}, 
												 end => $region->{end}, 
												 strand => $region->{strand} });
    }	

	my $region_length = $region->{end} - $region->{start} + 1; 	
	
    # total bases covered
    map { if($_ > 0) { $total_depth += $_; $bps_covered++;} } @$coverage;
	
    # calculate average depth of areas with coverage
    my $depth = sprintf("%.2f", ($bps_covered > 0 ) ? $total_depth / $bps_covered : 0);
	
    # calculate percent coverage of the region
	my $percent_coverage = sprintf("%.2f", $bps_covered / $region_length * 100);

    return { 'name' => $region->{'name'}, 
			 'id' => $region->{'id'},
			 'seq_id' => $region->{seq_id}, 
			 'coverage' => $percent_coverage, 
			 'bps_covered' => $bps_covered,
			 'depth' => $depth, 
			 'start' => $region->{start},
			 'end' => $region->{end},
			 'bps' => $region_length,
			 'hits' => $total_depth };
}


sub get_alignment_strand {    
    my $a = shift;

    #single-stranded read
    if(!$a->paired){
	return $a->strand;
    }

    if($library_type eq "firststrand"){ #the second (right most) mate dictates the strand

	#properly paired paired-end reads
	if($a->proper_pair){

	    my @second = $a->get_tag_values('SECOND_MATE');
	    if($second[0] == 1){
		return $a->strand;
	    }
	    else{#if first mate, we use the second mate's strand
		return $a->mstrand;
	    }
	}
	#orphaned, a bit complicated.  reverse the strand if first mate
	else{
	    my @first = $a->get_tag_values('FIRST_MATE');
	    if($first[0] == 1){
		return ($a->strand > 0) ? -1 : 1;
	    }
	    else{
		return $a->strand;
	    }
	}
    }

    elsif($library_type eq "secondstrand"){ #the first (left most) mate dictates the strand

	#properly paired paired-end reads
	if($a->proper_pair){

	    my @first = $a->get_tag_values('FIRST_MATE');
	    if($first[0] == 1){
		return $a->strand;
	    }
	    else{#if second mate, we use the first mate's strand
		return $a->mstrand;
	    }
	}
	#orphaned, a bit complicated.  reverse the strand if second mate
	else{
	    my @second = $a->get_tag_values('SECOND_MATE');
	    if($second[0] == 1){
		return ($a->strand > 0) ? -1 : 1;
	    }
	    else{
		return $a->strand;
	    }
	}
    }

    die "ERROR: Could not determine the strand for read:".$a->query->name."!";
}

1;
