package PFunc::EvidenceParser::ValidEvidenceTypes;

use strict;
use warnings;
use base qw(Exporter);
 
our %valid_evidence_types;
our @EXPORT = qw( %valid_evidence_types );

%valid_evidence_types = 
    ('BER' => 'PFunc::EvidenceParser::BERParser',
     'HMM' => 'PFunc::EvidenceParser::HMMParser',
     'PriamEC' => 'PFunc::EvidenceParser::PriamECParser',
     'Lipoprotein' => 'PFunc::EvidenceParser::LipoproteinMotifParser',
     'LipoP' => 'PFunc::EvidenceParser::LipoPParser',
     'TMHMM' => 'PFunc::EvidenceParser::TMHMMParser',
     'Hypothetical' => 'PFunc::EvidenceParser::Hypothetical',
     );
     
