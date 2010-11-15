package UnirefAnnotation::Lookup;

use strict;
use warnings;
use PolypeptideAnnotation::Database;
use PolypeptideAnnotation::NestedDB;

use Data::Dumper;

my $OPTIONS = {
    'PolypeptideAnnotation::Database' => 
    {
        "username" => { 'required' => 1 },
        "password" => { 'required' => 1 },
        "test" => { 'required' => 0, 'default' => 1 },
    },
    'PolypeptideAnnotation::NestedDB' =>
    {
        "path" => { 'required' => 1 }
    }
};

## Only takes one arg, a string of options
## key=val key=val ...
## Will make the correct object based on
## those params
sub new {
    my ($class, $opts) = @_;
    
    my $self = {};
    bless( $self, $class );
    $self->init( $opts );
}


sub init {
    my ($self, $kvs) = @_;
    my %opts;
    map { $opts{$1} = $2 if( /([^=]+)\s*=\s*([^=]+)/ ) } split(/\s+/, $kvs);
    my $type;
    my %defaults;
    
  OBJECT:
    foreach my $obj_string ( keys %{$OPTIONS} ) {
      OPTION:
        foreach my $opt( keys %{$OPTIONS->{$obj_string}} ) {
            if( $OPTIONS->{$obj_string}->{$opt}->{'required'} &&
                !exists( $opts{$opt} ) ) {
                %defaults = ();
                next OBJECT;
            }

            $defaults{$opt} = $OPTIONS->{$obj_string}->{$opt}->{'default'} 
            unless( $OPTIONS->{$obj_string}->{$opt}->{'required'} );
            
        }
        $type = $obj_string;
        last OBJECT;
    }

    die("Could not determine object type from parameters passed: $kvs")
        unless( $type );
    map { $defaults{$_} = $opts{$_} } keys %opts;
    my @obj_opts;
    map { push(@obj_opts, "$_ => '$defaults{$_}'") } keys %defaults;
    my $eval_string = "\$obj = new $type( ".join(", ", @obj_opts)." )";
    my $obj;
    eval($eval_string);
    die("Could not create object type: $type [$@]\n") unless( defined( $obj ) );
    $obj;
}
1;
