package PFunc::AnnotationOrder;

use strict;
use warnings;
use base qw(Exporter);
our @annotation_order;
our @EXPORT_OK  = qw(@annotation_order);


# an array of array refs. Each reference defines an annotation source as well
# as each annotation type (i.e. common name, EC number, etc).  
# $annotation_order[0] is the highest priority annotation type and will replace others.
@annotation_order = 
    ( 
      [ 'HMM::equivalog', ['all'] ],
      [ 'BER::characterized::full::full', ['all'] ],
      [ 'HMM::equivalog_domain', ['all'] ],
      [ 'Priam', ['EC'] ],
      [ 'BER::characterized::partial::full', ['all'] ],
      [ 'HMM::subfamily', ['all'] ],
      [ 'HMM::superfamily', ['all'] ],
      [ 'HMM::subfamily_domain', ['all'] ],
      [ 'HMM::domain', ['all'] ],
      [ 'HMM::pfam', ['all'] ],
      [ 'BER::characterized::full::partial', ['all'] ],
      [ 'TMHMM', ['all'] ],
      [ 'LipoproteinMotif', ['all'] ],
      [ 'HMM::hypoth_equiv', ['all'] ],
      [ 'BER::uncharacterized::full::full', ['all'] ],
      [ 'BER::uncharacterized::partial::full', ['all'] ],
      [ 'BER::uncharacterized::full::partial', ['all'] ],
      [ 'BER::conserved_hypothetical', ['all'] ],
      [ 'hypothetical', ['all'] ],
      );

1;
