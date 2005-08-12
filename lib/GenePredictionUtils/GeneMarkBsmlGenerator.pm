package GeneMarkBsmlGenerator;
@ISA = qw(GenePredictionUtils::GenePredictionBsmlGenerator);

use strict;
use warnings;

use GenePredictionUtils::GenePredictionBsmlGenerator;

sub new
{
	my $class = shift;
	my $this = {};
	bless $this, $class;
	$this->Init;
	return $this;
}

sub Init
{
	my $this = shift;
	$this->SUPER::Init();
}

1;
