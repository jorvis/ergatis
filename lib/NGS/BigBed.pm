package NGS::BigBed;

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
}

1;
