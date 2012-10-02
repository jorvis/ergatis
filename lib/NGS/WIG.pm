package NGS::WIG;

use strict;
use warnings;
use Data::Dumper;

sub new {
	my ($class, %args) = @_;
	my $self = bless {}, ref($class) || $class;
	$self->_init(\%args);
	return $self;
}

sub _init {
	my ($self, $args) = @_;

	die "Must define a sizes file, and the WIG file !" if(!defined $args->{wig_file} || !defined $args->{sizes_file});

	$self->{wig_file} = $args->{wig_file};
	$self->{sizes_file} = $args->{sizes_file};

	my %coverage;
	my @seq_ids;
	my %sizes;

	open SIZES, $args->{sizes_file};
	while(<SIZES>) {
	    chomp $_;
	    my ($seq_id, $len) = split("\t", $_);
		
	    $sizes{$seq_id} = $len;
	    push @seq_ids, $seq_id;

	    if(!defined $coverage{$seq_id}){
			$coverage{$seq_id} = {};
	    }
	}
	close SIZES;

	my $seq_id;
	my $type; #either variable or fixed
	my $vstart;

	my $vspan;
	my $vstep;
	my $val;

	my $a_start;

	open WIG, $args->{wig_file};

	while(<WIG>){
	    chomp $_;
	    next if($_ =~ /track\stype/);
	    
	    if($_=~/chrom=(\S+)/){
			$seq_id = $1;
			$a_start = -1;
			
			if($_ =~ /span=(\d+)/){ #init the span
				$vspan = $1;
			} else {
				$vspan = 1;
			}
			
			if( $_ =~ /fixedStep/ ) {
				$type = "fixed";
				
				if($_ =~ /start=(\d+)/){ #for fixedStep, init the start, stepsize
					$vstart = $1;
				}
				if($_ =~ /step=(\d+)/){
					$vstep = $1;
				}
				
			} else {
				$type = "variable";	
				$vstep = 1;
			}
			
	    } else {
			
			die "Could not identify WIG file type (fixed or variable)!" if(!defined $type);
			
			if($type eq "fixed"){
				$val = abs($_);
				if($a_start == -1) { $a_start = $vstart;}
				else{
					$a_start += $vstep;
				}
			} elsif($type eq "variable") {
				my ($pos, $v) = split("\t", $_);
				$vstart = $pos;
				$a_start = $vstart;
				$val = abs($v);
			}
			
			die "The wig file refers to a seq_id $seq_id that is not present in the sizes file!" if( !defined $coverage{$seq_id} );
			
			for (my $i = 0; $i<$vspan; $i++){
				$coverage{$seq_id}{($a_start+$i)} = $val;
			}
	    }	          
	}

	close WIG;
	
	$self->{coverage} = \%coverage; #the per-bp coverage is stored as a 1-based hash
	$self->{seq_ids} = \@seq_ids;
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

# Requires a seq_id
# Optionally takes start, stop coordinates, 1-based
# Returns a 0-based array of coverage depth at each base pair
#
sub get_per_bp_coverage {
    my ($self, $args) = @_;

    my $start;
    my $end;

    die "Must define a seq_id to retrieve coverage on a bp level!" if(! defined($$args{seq_id}));
    
    if(! defined($$args{start}) || ! defined($$args{end}) ){
		$start = 1;
		$end = $self->{sizes}{$$args{seq_id}};
		
    } else {
		$start = $$args{start};
		$end = $$args{end};
    }
	
    my @data;

    foreach( my $i = $start; $i<=$end; $i++ ) {
		
		if( defined( $self->{coverage}{$$args{seq_id}}{$i} ) ) {
			push @data, $self->{coverage}{$$args{seq_id}}{$i};
		} else {
			push @data, 0;
		}
    }
    return \@data;
}

#
# Calculates the coverage percentage, avg depth for a region
# Requires the seq_id, start, and end, in 1-based coordinates
#
sub calculate_region_coverage {
    my ($self, $region) = @_;  
  
    my $bps_covered = 0;
    my $total_depth = 0;
	
    my $coverage;
    my $hits;

    if( !defined $region->{strand} ) {
		$coverage = $self->get_per_bp_coverage( { seq_id => $region->{seq_id}, 
												  start => $region->{start}, 
												  end => $region->{end} });
		
    } else {
		# TODO - stranded
		#($coverage, $hits) = $self->get_per_bp_coverage_and_hits({ seq_id => $region->{seq_id}, 
		#														   start => $region->{start}, 
		#														   end => $region->{end}, 
		#														   strand => $region->{strand} });
    }	
	
    # total bases covered
    map { if($_ > 0) { $total_depth += $_; $bps_covered++;} } @$coverage;
	
    my $region_length = @$coverage;
	
    # calculate average depth of areas with coverage
    my $depth = sprintf("%.2f", ($bps_covered > 0 ) ? $total_depth / $bps_covered : 0);
	
    # calculate percent coverage of the region
    my $percent_coverage = sprintf("%.4f", ($region_length > 0) ? $bps_covered / $region_length * 100 : 0);
	
    return { 'name' => $region->{'name'}, 
			 'seq_id' => $region->{seq_id}, 
			 'coverage' => $percent_coverage, 
			 'bps_covered' => $bps_covered,
			 'depth' => $depth, 
			 'bps' => $region_length,
			 'hits' => $total_depth };
}

1;
