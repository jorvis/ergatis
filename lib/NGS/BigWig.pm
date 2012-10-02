package NGS::BigWig;

use strict;
use warnings;
use Data::Dumper;
use Bio::DB::BigWig 'binMean','binStdev';

sub new {
	my ($class, %args) = @_;
	
	my $self = bless {}, ref($class) || $class;
	$self->_init(\%args);
	return $self;
}

sub _init {
	my ($self, $args) = @_;
	
	die "Must define the BigWig file!"  if(!defined $args->{big_wig_file});

	$self->{big_wig_file} = $args->{big_wig_file};
	$self->{_bwig} = Bio::DB::BigWig->new( -bigwig => $self->{big_wig_file} );
	
	my @seq_ids = $self->{_bwig}->seq_ids;
	$self->{seq_ids} = \@seq_ids;

	my %sizes = map { $_ => $self->{_bwig}->length($_) } @seq_ids;
	$self->{sizes} = \%sizes;
}

sub seq_ids(){
    my ($self) = @_;    
    return $self->{seq_ids};
}

sub sizes(){
    my ($self) = @_;
    return $self->{sizes};
}

#this assumes incoming coords are 1-based and outgoing will be in a 0-based array
sub get_per_bp_coverage {
    my ($self, $args) = @_;
    
    die "Must define a chromosome to retrieve coverage on a bp level!" if(! defined($$args{seq_id}));
    
    my $start;
    my $end;

    if(! defined($$args{start}) || ! defined($$args{end}) ){
		$start = 1;
		$end = $self->{sizes}{$$args{seq_id}};

    } else {
		$start = $$args{start};
		$end = $$args{end};
    }
    
    my @data = ((0) x ( ($end - $start) + 1) );
    
    #get_seq_stream expects 1-based coords and returns 1-based coords
    my $iter = $self->{_bwig}->get_seq_stream( -seq_id=>$$args{seq_id}, 
											   -start=>$start, 
											   -end=>$end );

    while( my $b = $iter->next_seq ) {
		for( my $i=$b->start; $i<=$b->end; $i++ ) {
			$data[$i-$start] = abs($b->score); # we convert to a 0-based array
		}
    }

    return \@data;
}

# Calculates the coverage percentage, avg depth for a region, counts pairedend mates separately
# Requires the seq_id, start, and end, and optionally takes the strand
sub calculate_region_coverage {
    my ($self, $region) = @_;  
    
    my $bps_covered = 0;
    my $total_depth = 0;

    my $coverage = $self->get_per_bp_coverage({ seq_id => $region->{seq_id}, 
												start => $region->{start}, 
												end => $region->{end} });

    # total bases covered
    my $region_length = $region->{end} - $region->{start} + 1;

    map { if($_ > 0) { $total_depth += $_; $bps_covered++; } } @$coverage;

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

1;


