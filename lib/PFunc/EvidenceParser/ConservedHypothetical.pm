package PFunc::EvidenceParser::ConservedHypothetical;

use strict;
use warnings;
use PFunc::Annotation;

##### class vars #####
my $conserved_hypo = {
    'gene_product_name' => "conserved hypothetical protein",
    'GO'                => ['GO:0008150', 'GO:0003674', 'GO:0005575'],
    'TIGR_Role'         => 156,
};
##########


# assigns values set in $conserved_hypo global hash
# or empty string if info is not available
sub _assign_as_conserved_hypothetical {
    my ($self, $annotation, $source, $type ) = @_;
    $annotation->clear_annotation;

    my @fields = PFunc::Annotation::get_valid_fields;

    foreach my $field ( @fields ) {
        my $value;
        if( exists( $conserved_hypo->{$field} ) ) {
            $value = $conserved_hypo->{$field};
        } else {
            $value = "";
        }
        $annotation->set( $field, $value, $source, $type );
    }
}
1;
