package IntervalTree;

=head1 NAME

  IntervalTree.pm - Builds a tree of intervals that can be searched for overlaps to points or other intervals.

=head1 DESCRIPTION

    For a description of the data structure and search algorithm, see wikipedia.

=head1 SYNOPSIS

my $iTree = new IntervalTree;

$iTree->addInterval($id, $start, $stop);
.
.
.
$iTree->buildTree();
my @overlaps = $iTree->searchInterval($left, $right);
my @otherOverlaps = $iTree->searchPoint( $point );

=head1 AUTHOR

Kevin Galens
kgalens@tigr.org

=cut


use strict;
use Data::Dumper;

sub new {
    my $class = shift;

    my $self = { 
        'max' => 0,
        'min' => undef
             };

    bless( $self, $class );
    return $self;
}

=item $intervalTree->addInterval( $id, $start, $stop );

B<Description:> Adds an interval to the tree.

B<Parameters:> $id - the id of the interval
    $start - begining of the interval
    $stop  - end of the interval

B<Returns:> Nothing

=cut

sub addInterval {
    my ($self, $id, $start, $stop) = @_;
    push(@{$self->{'intervals'}}, [$start,$stop,$id]);

    my ($min, $max) = ($self->{'min'}, $self->{'max'});
        
    if($stop > $max) {
        $self->{'max'} = $stop;
    }
    if(!defined($min) || $start < $min) {
        $self->{'min'} = $start;
    }
  
}

=item $intervalTree->buildTree();

B<Description:> Builds the data structure from the added intervals.

B<Parameters:> None

B<Returns:> 1 on success

=cut
sub buildTree {
    my $self = shift;
    $self->{'root'} = 
        &_separateIntervals( $self->{'intervals'}, ($self->{'max'} - $self->{'min'})/2 + $self->{'min'});

    return 1;

}


#Recursive private function.
#Takes a set of intervals and center point and returns a node of the tree
#(which is actually a tree itself).
sub _separateIntervals {    
    my ($intervals, $center) = @_;
    
    #Maximum of those on the right, and those on the left of this center.
    my ($l_max, $l_min, $r_max, $r_min) = (0, undef, 0, undef);

    #Data structure representing the node itself.
    my $node = {
        's_left'   => {},     #Node of intervals which are completely less than center
        's_right'  => {},     #Node of intervals which are completely higher than center
        's_center' => [],     #Intervals that overlap the center. 
        'point'    => $center #The center.
        };

    my ($leftList, $rightList); #holds the list of intervals that will be on the right of this node,
       #and on the left of this node.

    #Cycle through all the intervals and put them to the right or left of this one.
    foreach my $interval ( @{$intervals} ) {

        #Error checking, if the left is bigger than the right, die.
        if($interval->[0] > $interval->[1]) {
            die("Start greater than stop\n");
        }

        #Left and right
        my ($left, $right) = @{$interval};
        
        #If the right most point is less than the center (totally to the left (less than)).
        if($right < $center) {

            #Push it on the left list
            push(@{$leftList}, $interval); 

            #Check for new maximums or new minimums
            if($right > $l_max) {
                $l_max = $right;
            } 
            if( !defined($l_min) || $left < $l_min) {
                $l_min = $left;
            }

            #If the left most point is more than the center (Totally to the right (more than)).
        } elsif( $left > $center ) {

            #Push it on the right list.
            push(@{$rightList}, $interval);

            #Check for new max/min
            if($right > $r_max) {
                $r_max = $right;
            } 
            if( !defined($r_min) || $left < $r_min ) {
                $r_min = $left;
            }
            #Otherwise (must be overlapping the center
        } else {
            #Push it on the current.
            push(@{$node->{'s_center'}}, $interval);
        }
    }

    #Sort those intervals that overlap the center (first by left, then by right)
    my @sLeftTmp = sort { $a->[0] <=> $b->[0] } @{$node->{'s_center'}};
    $node->{'sort_left'} = \@sLeftTmp;
    my @sRightTmp= sort { $a->[1] <=> $b->[1] } @{$node->{'s_center'}};
    $node->{'sort_right'}= \@sRightTmp;


    #These are the recurcsive calls
    if(defined($leftList) && @{$leftList} > 0) {

        #If we have intervals on the left, find the center and make a tree.
        my $leftCenter = ($l_max - $l_min)/2 + $l_min;
        $node->{'s_left'} = &_separateIntervals($leftList, $leftCenter);
    }
    if(defined($rightList) && @{$rightList} > 0) {

        #If we have intervals on the right, find the center and make a tree.
        my $rightCenter = ($r_max-$r_min)/2 + $r_min;
        $node->{'s_right'} = &_separateIntervals($rightList, $rightCenter );
    }

    #Return the node (which is actually a tree).
    

    unless(defined($node->{'s_center'})) {
        print Dumper($node);
        die("Node s_center was not defined");
    }
    return $node;
}

=item $intervalTree->searchInterval( $start, $stop );

B<Description:> Searches intervals overlapping the interval passed in

B<Parameters:> $start - begining of the interval (lowest)
               $stop  - end of the interval (highest)

B<Returns:> Array of array references in no special order

    ( [ 0,526 ], [ 42,698 ] ... [ 664366,674211 ] )

=cut

sub searchInterval {
    my ($self, $s_left, $s_right) = @_;

    my @retval;

    #Start the recursion.
    die("Self root is undef") unless($self->{'root'});
    push(@retval, &_searchNode($self->{'root'}, $s_left, $s_right));
    
    return(@retval);

}

sub _searchNode {
    my ($node, $s_left, $s_right) = @_;
    my @retval;
   
    #If we overlap the center, we need to check both sides for overlaps.
    #Otherwise we only have to check one side.
    if($s_left < $node->{'point'} && $node->{'point'} < $s_right) {
        push(@retval, @{$node->{'s_center'}});
        push(@retval, &_searchNode($node->{'s_left'}, $s_left, $s_right))
            if( scalar( keys %{$node->{'s_left'}} ) );
        push(@retval, &_searchNode($node->{'s_right'}, $s_left, $s_right))
            if( scalar( keys %{$node->{'s_right'}} ) );
    } elsif($s_right < $node->{'point'} ) {

        #Add all the intervals with the left most point less than passed in
        my $counter = -1;
        foreach my $num ( @{$node->{'sort_left'}} ) {
            last if($num->[0] > $s_right);
            $counter++;
        }
        push(@retval, @{$node->{'sort_left'}}[0..$counter]);

        #Search now, to the left.
        push(@retval, &_searchNode($node->{'s_left'}, $s_left, $s_right)) if(scalar(keys %{$node->{'s_left'}}));
    } elsif($s_left => $node->{'point'} ) {

        #Add all the intervals with the left most point less than passed in
        my $counter = -1;
        foreach my $num ( @{$node->{'sort_right'}} ) {
            last if($num->[1] <= $s_left);
            $counter++;
        }
        push(@retval, @{$node->{'sort_right'}}[0..$counter]);

        #Search now, to the right.
        push(@retval, &_searchNode($node->{'s_right'}, $s_left, $s_right)) if(scalar(keys %{$node->{'s_right'}}));
    }
        
    return @retval;
    
}


sub searchPoint {

}

1;
