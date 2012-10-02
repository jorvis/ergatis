package NGS::BED;

use strict;
use warnings;
use Data::Dumper;
use Bio::FeatureIO;
use Bio::SeqFeature::Annotated;

sub new {
	my ($class, %args) = @_;
	
	my $self = bless {}, ref($class) || $class;
	$self->_init(\%args);
	return $self;
}

sub _init {
	my ($self, $args) = @_;
	
	my $bed = Bio::FeatureIO->new( -file => $args->{file},
								   -format => "bed",
								   -name => $args->{name});

	$self->{_bed} = $bed;
}

sub write_bed {
	my ($self, $features) = @_;

	foreach my $feat (@$features) {
		# must change object type from SeqFeature::Generic to SeqFeature::Annotated
		my $af = Bio::SeqFeature::Annotated->new();
		$af->from_feature( $feat );

		# write feature to bed file
		$self->{_bed}->write_feature( $af );
	}
}

1;
