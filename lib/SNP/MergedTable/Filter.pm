package SNP::MergedTable::Filter;

use strict;
use warnings;
use IntervalTree;

my %valid_filters = (
    'excluded_regions' => sub { filter_excluded_regions(@_) }
);

sub filter_is_valid {
    my ($filter_name) = @_;
    return exists( $valid_filters{$filter_name} );
}

sub get_valid_filters {
    return keys %valid_filters;
}

sub filter {
    my ($filter, $table, @args) = @_;
    die("A MergedTable is required for filtering") unless( defined( $table ) && ref($table) eq 'SNP::MergedTable' );
    die("Invalid filter: $filter") unless( &filter_is_valid( $filter ) );
    $valid_filters{$filter}->($table, @args);
}


sub filter_excluded_regions {
    my ($table, @regions) = @_;
    die("filter_excluded_regions requires a set of regions and a MergedTable as arguments")
        unless( defined( $table ) && ref( $table ) eq 'SNP::MergedTable' && @regions > 0 );
    
    my %int_trees = ();
    my $fake_id = 1;
    
    foreach my $region ( @regions ) {
        die("Region did not contain 3 elements. Check region input into filter_excluded regions.")
            unless( @{$region} == 3 );
        my ($molecule, $left, $right) = @{$region};
        ($left, $right) = ($right, $left) if( $left > $right );
        
        unless( exists( $int_trees{$molecule} ) ) {
            $int_trees{$molecule} = new IntervalTree();
        }
        $int_trees{$molecule}->addInterval( $fake_id++, $left, $right );
    }
    map { $_->buildTree(); } values %int_trees;
    
    my $removed = 0;
    foreach my $row ( $table->get_rows ) {
       my $mol = $row->molecule();
       next unless( exists( $int_trees{$mol} ) );
       my @overlaps = $int_trees{$mol}->searchInterval( ($row->refpos -1), $row->refpos );
       if( @overlaps > 0 ) {
           my $c = $table->remove_row( $row );
           $removed += $c;
       }
    }
    
    return $removed;
}
