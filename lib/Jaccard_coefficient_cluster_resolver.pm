package main;

our $SEE;

package Jaccard_coefficient_cluster_resolver;

use strict;

BEGIN {
    require '/usr/local/devel/ANNOTATION/ard/chado-v1r5b1/lib/site_perl/5.8.5/SingleLinkageClusterer.pm';
    import SingleLinkageClusterer;
}

=head1 NAME

package Jaccard_coefficient_cluster_resolver


=cut


=head1 DESCRIPTION

    Module is used to resolve the clusters within a graph by using the Jaccard's similarity coefficient to break edges within loosely coupled clusters.

    The Jaccard's similarity coefficient is defined as:

    given nodes A and B, 
    X = set of nodes connected to A including A
    Y = set of nodes connected to B including B

    
                     #(nodes intersecting X and Y)            (X && Y)
    Jlink (A,B)  =   -----------------------------     =     ----------
                     #(nodes in set X union set Y)            (X || Y)

    

    Jlink =  1 when nodes A and B are identically connected.

    Jlink = 0 when nodes A and B have no connected neighbors in common, and are themselves unconnected.

    For connected nodes A and B which have different neighbors, the Jlink score will be between 0 and 1, providing a similiarity coefficient for the level of similarity between sets of connections.



    Another description of the Jaccard similarity coefficient: (http://www.ergometrika.org/Volume3/mulqueen-rev-2003.htm)


The Jaccard index was originally developed to assess similarity among distributions of flora in different geographic areas (Jaccard, 1912). The procedure results in a matching coefficient for binary variables in which joint absences are excluded from both the denominator and the numerator and equal weight is given to matches and non-matches:
SJ = a/(a+b+c) x 100, where

SJ = Jaccard similarity coefficient,

a = number of elements shared by all groups,

b = number of elements unique to the first group, and

c = number of elements unique to the second group.



Jaccard, P. (1912). The distribution of flora in the alpine zone. The New Phytologist, 11(2), 37-50.


=cut


=over 4

=item new()

B<Description:> Constructor: Instantiate an object of Jaccard_coefficient_cluster_resolver 

B<Parameters:> $linkScore

$linkScore is a real number between 0 and 1 defined as the Jaccard similarity coefficient.


B<Returns:> $Jaccard_coefficient_cluster_resolver_object 

=back

=cut


sub new {
    my $packagename = shift;
    
    my $linkScore = shift;

    unless ($linkScore >= 0 && $linkScore <= 1) {
	die "Invalid link score ($linkScore). \n\n 0 <= link_score <= 1\n\n\n";
    }
    my $self = { 
	linkScore => $linkScore
	};
    bless ($self, $packagename);
    return ($self);

}


=over 4

=item resolve_clusters()

B<Description:> Given a set of paired elements, the edges between pairs are removed if the elements have a link score less than that set in the constructor.

B<Parameters:> @pairs

@pairs is a list of paired elements in the form of array references.  For example:

    @pairs = ( 
	       [ a, b],
	       [b, c],
	       [e, f] 
	       )

B<Returns:> @clusters


@clusters is a list of array references where each array references provides a list of elements belonging to a single cluster.

For example, given the inputted pairs above as input, 
    
    @clusters = ( 
		  [a, b, c],
		  [e, f], 
		  )
    
    
See SingleLinkageClusterer.pm for more details.


=back

=cut


sub resolve_clusters {
    my $self = shift;
    my @pairs = @_;
    
    my $linkScore = $self->{linkScore};

    ## Transform each pair to a data structure which provides sort of a graph, with each element (node) pointing to a list of other elements(nodes), implemented with a hash{node_id} = (node id list)
    
    my %inputGraph;

    foreach my $pair (@pairs) {
	my ($a, $b) = @$pair;
	
	my $a_aref = $inputGraph{$a};
	unless ($a_aref) {
	    $a_aref = $inputGraph{$a} = [];
	}
	
	my $b_aref = $inputGraph{$b};
	unless ($b_aref) {
	    $b_aref = $inputGraph{$b} = [];
	}

	&add_element($a_aref, $b);
	&add_element($b_aref, $a);
    }

    ## Now examine the Jaccard similarity coefficient between each member
    my @resolved_pairs; #all pairs meeting link score restrictions.
    my %seen; #avoid analyzing the same pairs twice.

    my @all_elements = keys %inputGraph;
    
    foreach my $element (@all_elements) {
	my $connected_elements_aref = $inputGraph{$element};
	
	my @listA = ($element, @$connected_elements_aref);
	
	foreach my $connected_element (@$connected_elements_aref) {
	    
	    if ($seen{$element}->{$connected_element}) { next;} #avoid duplicate comparisons.
	    
	    my $neighbors_list_aref = $inputGraph{$connected_element};

	    my @listB = ($connected_element, @$neighbors_list_aref);

	    my $current_link_score = &calculate_Jaccard_coeff(\@listA, \@listB);
	    
	    if ($current_link_score >= $linkScore) { #passed test
		push (@resolved_pairs, [$element, $connected_element]);
	    }

	    $seen{$element}->{$connected_element} = 1;
	    $seen{$connected_element}->{$element} = 1;
	}
    }
    
    my @clusters;

    if (@resolved_pairs) {
	@clusters = &SingleLinkageClusterer::build_clusters(@resolved_pairs);
    }
    
    return (@clusters);

}



#private
sub add_element {
    my ($list_aref, $element) = @_;
    
    my $has_element = 0;
    foreach my $existing_element (@$list_aref) {
	if ($element eq $existing_element) {
	    $has_element = 1;
	    last;
	}
    }
    if (! $has_element) {
	push (@$list_aref, $element);
    }
}

#private
####
sub calculate_Jaccard_coeff {
    my ($a_list_aref, $b_list_aref) = @_;

    if ($SEE) {
	print "\nList a: (@$a_list_aref)\n"
	    . "List b: (@$b_list_aref)\n";
    }
    
    ## Track each element
    my (%all_entries, %a_entries, %b_entries);
    
    foreach my $element (@$a_list_aref) {
	$a_entries{$element} = 1;
	$all_entries{$element} = 1;
    }

    foreach my $element (@$b_list_aref) {
	$b_entries{$element} = 1;
	$all_entries{$element} = 1;
    }

    ## Determine number of unique entries.
    my @unique_entries = keys %all_entries;
    my $num_unique_entries = $#unique_entries + 1;
    
    print "$num_unique_entries unique entries between a and b: (@unique_entries)\n" if $SEE;
	
    ## See how many a and b have in common.
    my $num_in_common = 0;
    foreach my $element (@unique_entries) {
	if ($a_entries{$element} && $b_entries{$element}) {
	    $num_in_common++;
	    print "($element) is common.\n" if $SEE;
	}
    }

    my $link_score = ($num_in_common / $num_unique_entries);
    
    print "link score (a,b) = ($num_in_common / $num_unique_entries) = $link_score\n" if $SEE;
    
    return ($link_score);
}


1; #EOM
