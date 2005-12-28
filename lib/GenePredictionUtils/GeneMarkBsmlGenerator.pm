package GenePredictionUtils::GeneMarkBsmlGenerator;
@ISA = qw(GenePredictionUtils::GenePredictionBsmlGenerator);

use strict;
use warnings;

BEGIN {
    require '/usr/local/devel/ANNOTATION/ard/chado-v1r5b1/lib/site_perl/5.8.5/GenePredictionUtils/GenePredictionBsmlGenerator.pm';
    import GenePredictionUtils::GenePredictionBsmlGenerator;
}

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
