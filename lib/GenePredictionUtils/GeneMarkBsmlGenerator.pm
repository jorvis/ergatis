package GenePredictionUtils::GeneMarkBsmlGenerator;
@ISA = qw(GenePredictionUtils::GenePredictionBsmlGenerator);

use strict;
use warnings;
use GenePredictionUtils::GenePredictionBsmlGenerator;

## keep things properly scoped
{

my %_atts = ( 
    id_generator => undef,
);

sub new
{
	my ($class, %args) = @_;
    
	my $this = bless { %_atts }, $class;
	$this->Init;
    
    ## an id_generator reference is required
    if ( ! $args{id_generator} ) {
        die "you must pass an id_generator object reference when creating a new GenePredictionBsmlGenerator";
    }
    
    $this->{id_generator} = $args{id_generator};
    
	return $this;
}

sub Init
{
	my $this = shift;
	$this->SUPER::Init();
}

}
1;
