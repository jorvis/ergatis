#A wrapper class around various types of alignment and coverage file formats

package NGS::Alignment;

use strict;
use warnings;
use NGS::BAM;
use NGS::WIG;
use NGS::BigWig;

use Data::Dumper;

#initialization requires file_type argument (BAM, WIG, or BigWig currently supported)
#also requires file, being the path to the alignment/coverage file
#BAM files can be used to retrieve data such as hits by location/strand etc
#Coverage files are flat, so only include depth of coverage, and require seperate instances to track each strand


sub new {
	my ($class, %args) = @_;
	
	my $self = bless {}, ref($class) || $class;
	$self->_init(\%args);
	return $self;
}


sub _init {
    my ($self, $args) = @_;

    die "Must provide file type!" if(!defined $args->{file_type});
    die "Must provide alignment or coverage file!" if(!defined $args->{file});
    
    $args->{stranded} = 0 if(!defined $args->{stranded});

    my $alignment;

    if($args->{file_type} eq "BAM") {
		
		if(defined $args->{total_mapped_reads}){
			$alignment = NGS::BAM->new( bam => $args->{file}, 
										stranded => $args->{stranded}, 
										total_mapped_reads => $args->{total_mapped_reads}, 
										uniq => $args->{uniq} );
		} else{
			$alignment = NGS::BAM->new( bam => $args->{file}, 
										stranded => $args->{stranded},
										uniq => $args->{uniq} );
		}
		
    } elsif( uc($args->{file_type}) eq "WIG") {
		die "WIG files must have a sizes file!" if(!defined $args->{sizes_file});
		$alignment = NGS::WIG->new( wig_file => $args->{file}, 
									sizes_file => $args->{sizes_file} );
		
    } elsif( uc($args->{file_type}) eq "BIGWIG") {
		$alignment = NGS::BigWig->new( big_wig_file => $args->{file} );
    }
    # TODO add other file type support
    else{
		die "Unsupported file type ".$args->{file_type};
    }

    $self->{_alignment} = $alignment;
    $self->{_stranded} = $args->{stranded};
    $self->{_file_type} = $args->{file_type};
    $self->{_file} = $args->{file};
    $self->{_total_mapped_reads} = $args->{total_mapped_reads} if(defined $args->{total_mapped_reads});
}

sub seq_ids(){
    my ($self) = @_;
    return $self->{_alignment}->seq_ids;
}

sub sizes(){
    my ($self) = @_;
    return $self->{_alignment}->sizes;
}

sub get_total_mapped_reads(){

    my ($self) = @_;

    if(defined $self->{_total_mapped_reads}) {
		return $self->{_total_mapped_reads};
		
    } elsif($self->{_file_type} eq "BAM") {
		return $self->{_alignment}->get_total_mapped_reads();
		
    } else {
		die "Cannot calculate total mapped reads for non BAM files!\n It must be passed in.\n";
    }
    
}
 
#incoming coordinates are 1-based, to match GFF standards, and outgoing data is in a 0-based array   
sub get_per_bp_coverage(){
    my ($self, $args) = @_;
    return $self->{_alignment}->get_per_bp_coverage($args);
}

# Incoming coordinates are 1-based, to match GFF standards		
sub calculate_region_coverage(){
    my ($self, $args) = @_;    
    return $self->{_alignment}->calculate_region_coverage($args);
}

#
# Get coverage over the entire genome
#
sub get_genome_coverage {
    my ($self) = @_;
	
    my @coverages;
    my $seq_ids = $self->seq_ids();
    my $sizes = $self->sizes();
	
    foreach my $seq_id (@{$seq_ids}) {
		my $end = $sizes->{$seq_id};
		
		push @coverages, $self->calculate_region_coverage({ name => $seq_id, 
															seq_id => $seq_id, 
															start => 1, 
															end => $end });
    }
	
    return \@coverages;
}
    
    
#
# Calculates coverage over a given feature type (mRNA, transcript etc)
# Requires a gff3 object and feature type
#

sub get_feature_coverages {
    my ($self, $args) = @_;
    
    my $feature = $args->{feature};
    my $gff3 = $args->{gff3};	

    # tracks coverage for each feature (gene, CDS etc)
    my @coverages;

    # get the chromosome ids
    my $seq_ids = $self->seq_ids;

    foreach my $seq_id (sort {$a cmp $b } @{$seq_ids}) {
		
		# get the coords for the feature type we are interested in
		my $features = $gff3->get_features( seq_id => $seq_id, feat_type => $feature );
		
		if( !defined $args->{stranded} ) {
			
			foreach my $feat (@{$features}) {
				
				my $name = $feat->primary_id;
				if( $feat->has_tag('Name') ) {
					my @vals = $feat->get_tag_values('Name');
					$name = $vals[0];
				}

				push @coverages, $self->calculate_region_coverage({ id => $feat->primary_id,
																	name => $name,
																	seq_id => $seq_id, 
																	start => $feat->start, 
																	end => $feat->end });
			}
			
		} else {
			
			foreach my $feat (@{$features}) {
				push @coverages, $self->calculate_region_coverage({ name => $feat->primary_id, 
																	seq_id => $seq_id, 
																	start => $feat->start, 
																	end => $feat->end , 
																	strand => $feat->{strand} });
			}
		}
    }
	
    return \@coverages;
}

1;
