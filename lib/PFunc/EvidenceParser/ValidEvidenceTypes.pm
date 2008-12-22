package PFunc::EvidenceParser::ValidEvidenceTypes;

use strict;
use warnings;
use base qw(Exporter);
 
our %valid_evidence_types;
our @EXPORT = qw( %valid_evidence_types );

%valid_evidence_types = 
    ('BER' => 'PFunc::EvidenceParser::BERParser2',
     'HMM' => 'PFunc::EvidenceParser::HMMParser',
     'PriamEC' => 'PFunc::EvidenceParser::PriamECParser',
     'Lipoprotein' => 'PFunc::EvidenceParser::LipoproteinMotifParser',
     'TMHMM' => 'PFunc::EvidenceParser::TMHMMParser',
     'Hypothetical' => 'PFunc::EvidenceParser::Hypothetical',
     );
     
