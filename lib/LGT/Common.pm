package LGT::Common;

=head1 NAME

LGT::Common.pm - Collection of common or shared LGT functions

=head1 SYNOPSIS

Need to put something useful here

=head1 DESCRIPTION

A module that has functions for commonly used actions, or were previously duplicated among several LGT packages

=head1 AUTHOR - Shaun Adkins

e-mail: sadkins@som.umaryland.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Exporting functions so we don't have to include package name
use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(find_lca parse_flag);

use strict;
use warnings;


=head2 find_lca

 Title	 : find_lca
 Function: Determine new LCA from currently stored LCA and new lineage information
 Returns : String of he new lowest common ancestor (LCA)
 Args	 : Array reference with 2 elements:
	[1] - Assigned lineage for the current hit
	[2] - Currently assigned LCA for the query hit
=cut

sub find_lca {
    my $lineages = shift;

    # prime it
    my @lca = split( ';', $lineages->[0] );

    # if LCA was never defined, hit lineage will be LCA
    # Each additional time, LCA will only extend to the deepest common ancestor
    foreach my $l (@$lineages) {
        my $newlca = []; 
        my @lineage = split( ';', $l );
        for ( my $i = 0; $i < @lca; $i++ ) { 
            if ( $lca[$i] eq $lineage[$i] ) { 
                push( @$newlca, $lineage[$i] );
            } else {
                last;
            }
        }
        @lca = @$newlca;
    }
    return join( ';', @lca );
}

=head2 parse_flag

 Title   : parse_flag
 Function: Determines properties of a SAM read based on flag bit number 
 Returns : A hash containing 0/1 values for individual properties
 Args    : 
    [1] - A bitwise flag int from the 5th field in a SAM read when using 'samtools view'

=cut

sub parse_flag {
    my $int    = shift;
    my $rawbin = &_dec2bin($int);
    my $rev    = scalar $rawbin;
    my $bin = sprintf( "%012d", $rev );
    my $final_bin = reverse $bin;
    return {
        'paired'        => substr( $final_bin, 0,  1 ), 
        'proper'        => substr( $final_bin, 1,  1 ), 
        'qunmapped'     => substr( $final_bin, 2,  1 ), 
        'munmapped'     => substr( $final_bin, 3,  1 ), 
        'qrev'          => substr( $final_bin, 4,  1 ), 
        'mrev'          => substr( $final_bin, 5,  1 ), 
        'first'         => substr( $final_bin, 6,  1 ), 
        'last'          => substr( $final_bin, 7,  1 ), 
        'secondary'     => substr( $final_bin, 8,  1 ), 
        'failqual'      => substr( $final_bin, 9,  1 ), 
        'pcrdup'        => substr( $final_bin, 10, 1 ), 
        'supplementary' => substr( $final_bin, 11, 1 ), 
    };   
}


=head2 _dec2bin

 Title	 : _dec2bin
 Function: [Private method] Converts decimal to binary format
 Returns : The raw binary number
 Args	 :
	[1] - The decimal-formatted bitflag number

=cut 

sub _dec2bin {
    my $str = unpack( "B32", pack( "N", shift ) ); 
    $str =~ s/^0+(?=\d)//;    # otherwise you'll get leading zeros
    return $str;
}

1;
