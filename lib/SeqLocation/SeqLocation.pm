package SeqLocation::SeqLocation;

=head1 NAME

SeqLocation.pm - an object that holds sequence location information
    for use within the adjust_*_coordinates suite of scripts

=head1 SYNOPSIS

my $seqLocation = new SeqLocation::SeqLocation($id, $type, $start, $end);

=head1 DESCRIPTION

=cut


use warnings;
use strict;
use Ergatis::Logger;
use Data::Dumper;


####################################################################
#                      Constructor and Init                        #
####################################################################


sub new {
    my ($class, @args) = @_;
    my $self = {};
    bless $self, $class;

    $self->init(@args);
    return $self;
}

sub init {
    my $self = shift;
    my @args = @_;

    #Get the logger.
    $self->{'logger'} = Ergatis::Logger::get_logger();

    if(@args < 2) {
        $self->{'logger'}->logdie("SeqLocation object requires four arguments".
                                  " in constructor. (id, type)");
    }

    $self->{'id'} = $args[0];
    $self->{'type'} = $args[1];

    #The start hash data structure allows the quick retrieval of
    #seqLocations by the start coordinate.
    $self->{'start'} = {};
    

    #The end hash data structure allows the quick retrieval of
    #seqLocations by the end coordinate.
    #Format: $self->{'end'}->{$end}->{$start}->{('id'|'type'|'properties')}
    $self->{'end'} = {};
}

####################################################################
#                     Instance Methods                             #
####################################################################
sub addSeqLocation {
    my ($self, $id, $type, $start, $end) = @_;
    my $retval;

    
    print STDERR "Start is not defined\n" 
        if(!defined($self->{'overlapStart'}));
    print STDERR "End is not defined\n"
        if(!defined($self->{'overlapEnd'}));
    
    
    if(defined($self->{'overlapStart'}) && defined($self->{'overlapEnd'}) &&
       ($start >= $self->{'overlapStart'} && $start <= $self->{'overlapEnd'}) ||
       ($end <= $self->{'overlapEnd'} && $end >= $self->{'overlapStart'})) {

        $self->{'start'}->{$start}->{$end}->{'id'} = $id;
        $self->{'start'}->{$start}->{$end}->{'type'} = $type;
        
        $self->{'end'}->{$end}->{$start}->{'id'} = $id;
        $self->{'end'}->{$end}->{$start}->{'type'} = $type;
        $retval = 1;
    } else {
        $retval = 0;
    }

    return $retval;

}

sub checkOverlap {
    my ($self, $type, $start, $end, $properties, $seq) = @_;
    my $retval = 0;

    #Check the start coodinate
    if(defined($self->{'start'}->{$start}) && defined($self->{'start'}->{$start}->{$end})) {
        if($seq eq 'next') {
            $retval = 1;
        } else {
            $retval = 0;
        }
    } elsif(defined($self->{'start'}->{$start}) && $seq eq 'next') {
        foreach my $keyEnd(keys %{ $self->{'start'}->{$start} }) {
            last if($keyEnd > $end &&
                    ($retval = $self->compareProperties($self->{'start'}->{$start}->{$keyEnd}->{'properties'},$properties)));
        }
    } 


    #Check the end data structure
    if(defined($self->{'end'}->{$end}) && $seq eq 'prev' && !$retval) {
        foreach my $keyStart(keys %{ $self->{'end'}->{$end} } ) {
            last if($keyStart < $start &&  
                    ($retval = $self->compareProperties($self->{'end'}->{$end}->{$keyStart}->{'properties'},
                                             $properties)));
        }
    }
    
    return $retval;
    
}

sub compareProperties {
    my ($self, $pro1, $pro2) = @_;
    my $retval = 1;

    if(defined($pro1) && defined($pro2)) {

        foreach my $key1(keys %{$pro1}) {
            if(!defined($pro2->{$key1}) || 
               (defined($pro2->{$key1}) && $pro2->{$key1} ne $pro1->{$key1})) {
                $retval = 0;
            }
            last unless($retval);
        }
    } elsif(!defined($pro1) && !defined($pro2)) {
        $retval = 1;
    } else {
        $retval = 0;
    }
    
    return $retval;

}

sub addProperties {
    my ($self, $start, $end, $prop) = @_;

    $self->{'start'}->{$start}->{$end}->{'properties'} = $prop;
    $self->{'end'}->{$end}->{$start}->{'properties'} = $prop;

    die "Start and end hashes have different prop values"
        if($self->{'start'}->{$start}->{$end}->{'properties'} != 
           $self->{'end'}->{$end}->{$start}->{'properties'});
}

sub removeSeqLocation {
    my ($self, @loc) = @_;
    my $retval = 1;

    unless(defined($self->{'start'}->{$loc[2]}->{$loc[3]}->{'id'})) {
        $self->{'logger'}->logdie("Does not have ID");
    }
    unless(defined($self->{'start'}->{$loc[2]}->{$loc[3]}->{'type'})) {
        $self->{'logger'}->logdie("Does not have type");
    }

    if(defined($self->{'start'}->{$loc[2]}->{$loc[3]}) &&
       $self->{'start'}->{$loc[2]}->{$loc[3]}->{'id'} eq $loc[0] &&
       $self->{'start'}->{$loc[2]}->{$loc[3]}->{'type'} eq $loc[1]) {
        delete($self->{'start'}->{$loc[2]}->{$loc[3]})|| die
            "Could not delete @loc from start\n";
    } else {
        $retval = 0;
    }

    if(defined($self->{'end'}->{$loc[3]}->{$loc[2]}) &&
       $self->{'end'}->{$loc[3]}->{$loc[2]}->{'id'} eq $loc[0] &&
       $self->{'end'}->{$loc[3]}->{$loc[2]}->{'type'} eq $loc[1] &&
       $retval) {
        delete( $self->{'end'}->{$loc[3]}->{$loc[2]}) || die
            "Could not delete @loc from end\n";
    } else {
        $retval = 0;
    }

    return $retval;
}

sub setOverlapRange {
    my ($self, $overlapStart, $overlapEnd) = @_;
    $self->{'overlapStart'} = $overlapStart;
    $self->{'overlapEnd'} = $overlapEnd;

}
1;
