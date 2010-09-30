package PFunc::EvidenceParser::HMMParser;

#Required options: hmm_info, tigr_roled_db_dir

use strict;
use warnings;
use XML::Twig;
use File::OpenFile qw( open_file );
use Carp;
use Carp qw( cluck );
use Fcntl qw( O_RDONLY );
use MLDBM "DB_File";
use TIGR::Roles::HMM::PfamToRoleLookup;
use TIGR::Roles::HMM::TIGRFamToRoleLookup;
use PFunc::EvidenceParser::ConservedHypothetical;
use Data::Dumper;

use base qw(PFunc::EvidenceParser);

##### class vars #####
my $annotation_type = "HMM";
my $default_hmm_info = "coding_hmm.lib.db";

sub new {
    my ($class, %args) = @_;
    my $self = $class->SUPER::new( %args );
    bless($self, $class);
    $self->_init_hmm_parser( %args );
    return $self;
}

sub _init_hmm_parser {
    my ($self, %args) = @_;
    
    ## hmm info
    my %tmp;
    if( $args{'hmm_info'} ) {
        tie(%tmp, 'MLDBM', $args{'hmm_info'}, O_RDONLY )
            or die("Could not tie hash to $args{'ber_info'}");
    } else {
        die("Option hmm_info is required to create ".__PACKAGE__);
    }
    $self->{'_hmm_info'} = \%tmp;

    ## tigr_roles dir
    my $tigr_roles_db_dir;
    if( $args{'tigr_roles_db_dir'} ) {
        $tigr_roles_db_dir = $args{'tigr_roles_db_dir'}
    } else {
        die("Option tigr_roles_db_dir is required to create ".__PACKAGE__);
    }

    ## tigrfams dir
    my $tigrfams_dir;
    if( $args{'tigrfams_dir'} ) {
        $tigrfams_dir = $args{'tigrfams_dir'};
    } else {
        die("Option tigrfams_dir is required to create ".__PACKAGE__);
    }

    $self->_init_valid_isotypes();
    $self->_init_name_suffixes();
    $self->_init_tigr_roles_lookup($tigr_roles_db_dir, $tigrfams_dir);
    
}

sub _init_valid_isotypes {
    my ($self) = @_;
    my $a = $annotation_type;
    $self->{'_valid_isotypes'} = {
        $a."::equivalog"        => 1,
        $a."::equivalog_domain" => 2,
        $a."::subfamily"        => 3,
        $a."::superfamily"      => 4,
        $a."::subfamily_domain" => 5,
        $a."::domain"           => 6,
        $a."::pfam"             => 7,
        $a."::hypoth_equivalog" => 8,
    };
}

sub _init_name_suffixes {
    my ($self) = @_;
    my $a = $annotation_type;
    $self->{'_name_suffixes'} = {
        $a."::subfamily"        => {'string' => " family protein",
                                    'action' => "append" },
        $a."::superfamily"      => {'string' => " family protein",
                                    'action' => "append" },
        $a."::pfam"             => {'string' => " family protein",
                                    'action' => "append" },
        $a."::subfamily_domain" => {'string' => " domain protein",
                                    'action' => "append" },
        $a."::domain"           => {'string' => " domain protein",
                                    'action' => 'append' },
        $a."::hypoth_equivalog" => {'string' => "conserved hypothetical protein",
                                    'action' => "replace" },
    };
}

sub _init_tigr_roles_lookup {
    my ($self, $tr_db_dir, $tigrfams_dir) = @_;
    my $pfam_lookup = new TIGR::Roles::HMM::PfamToRoleLookup('roles_db_dir' => $tr_db_dir );
    my $tigr_lookup = new TIGR::Roles::HMM::TIGRFamToRoleLookup( 'roles_db_dir' => $tr_db_dir, 
                                                                 'tigrfams_dir' => $tigrfams_dir );
    $self->_tigr_role_lookup( 'PFAM', $pfam_lookup );
    $self->_tigr_role_lookup( 'TIGRFam', $tigr_lookup );
}

sub _parse {
    my ($self, $fh) = @_;

    my $twig = XML::Twig->new( 'twig_roots' => {
        'Seq-pair-alignment' => sub { $self->_process_hmm_coding_alignment( @_ ) },
    });

    $twig->parse($fh);
}

sub _is_valid_isotype {
    my ($self, $isotype) = @_;
    my $retval = 0;
    $retval = 1 if( exists( $self->{'_valid_isotypes'}->{$isotype} ) );
    return $retval;
}

sub _get_isotype_score {
    my ($self, $isotype) = @_;
    my $retval;
    if( $self->_is_valid_isotype($isotype) ) {
        $retval = $self->{'_valid_isotypes'}->{$isotype};
    }
    return $retval;
}

sub _hmm_info  {
    my ($self, $id) = @_;
    my $retval;
    if( !defined( $id ) ) {
        $retval = $self->{'_hmm_info'};
    } elsif( exists( $self->{'_hmm_info'}->{$id} ) ) {
        $retval = $self->{'_hmm_info'}->{$id};
    }
    return $retval;
}

sub _tigr_role_lookup {
    my ($self, $fam_type, $lookup) = @_;
    my $retval;
    if( defined( $lookup ) ) {
        $retval = $lookup;
        $self->{'_tigr_roles_lookup'}->{$fam_type} = $lookup;
    } else {
        if( exists( $self->{'_tigr_roles_lookup'}->{$fam_type} ) ) {
            $retval = $self->{'_tigr_roles_lookup'}->{$fam_type};
        }
    }
    return $retval;
}

sub _process_hmm_coding_alignment {
    my ($self, $t, $spa) = @_;
    
    # genomic prediction id
    my $ref_id = $spa->att('refseq');
    # hmm id
    my $comp_id = $spa->att('compseq');

    my ($total_score_att) = $spa->find_nodes('Attribute[@name="total_score"]');
    die("Could not find total_score Attribute element") unless( $total_score_att );
    my $total_score = $total_score_att->att('content') if( defined( $total_score_att ) );
    
    ## we need both IDs and the score to continue
    unless ( $ref_id && $comp_id && defined $total_score ) {
        die "failed to get ref_id, comp_id and total_score from Seq-pair-alignment";
    }

    #get the id to annotate on
    my $annotation_feature_id = $self->lookup_feature_id( $ref_id, $self->annotate_on );
    return if( !defined( $annotation_feature_id ) );

    ## get the annotation object
    my $annotation = $self->get_feature_annotation( $annotation_feature_id );
    return if( !defined( $annotation ) );
    
    ## make sure we have info on this HMM
    my $hmm_info = $self->_hmm_info( $comp_id );
    if ( $hmm_info ) {
        
        if ( $total_score >= $hmm_info->{'trusted_cutoff'} ) {
            
            ## what is the priority of this match, based on its type?
            my $isotype = $hmm_info->{'isotype'};
            my $current_type = $annotation_type."::".$isotype;
            if ( $self->_is_valid_isotype( $current_type ) ) {
                
                # get the score the previous annotations gene product name type;
                my $gpn_type_score = $self->_get_isotype_score( $annotation->_get_type( 'gene_product_name' ) );

                ## we can just compare with the product_score here since all scores are currently
                #   set together.  in the future this could change.
                if( !defined( $gpn_type_score ) || 
                    $self->_get_isotype_score( $current_type ) < $gpn_type_score ) {

                    #set the current annotation
                    die("Could not get current_type [$current_type] for $comp_id")
                        if( !defined( $current_type ) );
                    $self->_assign_annotation( $annotation, $hmm_info, $comp_id, $current_type );
                    
                } elsif( defined( $gpn_type_score ) && 
                         $self->_get_isotype_score( $current_type ) == $gpn_type_score ) {

                    #count how many annotation fields annotation currently has
                    my $cur_count = 0;
                    foreach my $field ( PFunc::Annotation::get_valid_fields ) {
                        $cur_count++ if( $annotation->has_annotation($field) );
                    }

                    #count how many
                    my $possible_count = 0;
                    $possible_count++ if( $hmm_info->{'hmm_com_name'} );
                    foreach my $field ( qw(gene_symbol ec_num go) ) {
                        $possible_count++ if( exists( $hmm_info->{ $field } ) &&
                                              defined( $hmm_info->{ $field } ) );
                    }

                    if( $possible_count > $cur_count ) {
                        #set the current annotation
                        $self->_assign_annotation( $annotation, $hmm_info, $comp_id, $current_type );
                    }

                }
                
            }
        }
        
    } else {
        die("Could not find hmm $comp_id in lookup");
    }
}

sub _assign_annotation {
    my ($self, $annotation, $hmm_info, $hmm_acc, $hmm_annot_iso) = @_;
    
    #clear all annotation
    $annotation->clear_annotation;

    #if the isotype of the hmm is 
    #set the new annotation $a."::hypoth_equivalog", should be conserved hypothetical
    die("HMM_ANNOT is not defined") unless( defined( $hmm_annot_iso ) );
    die("annotation_type is not defined") unless( defined( $annotation_type ) );
    if( $hmm_annot_iso eq $annotation_type."::hypoth_equivalog" ) {
        PFunc::EvidenceParser::ConservedHypothetical->_assign_as_conserved_hypothetical( $annotation, 
                                                                                         $hmm_acc,
                                                                                         $annotation_type."::hypoth_equivalog" );
        return;
    }

    #product name
    my $com_name = $self->_append_to_com_name( $hmm_info->{'hmm_com_name'}, $hmm_annot_iso );
    $annotation->set_gene_product_name( $com_name,
                                        $hmm_acc, $hmm_annot_iso );

    #gene symbol
    $annotation->set_gene_symbol( $hmm_info->{'gene_symbol'},
                                  $hmm_acc, $hmm_annot_iso );

    #ec number
    $annotation->set_EC( $hmm_info->{'ec_num'},
                         $hmm_acc, $hmm_annot_iso );

    # go
    $annotation->set_GO( $hmm_info->{'go'},
                         $hmm_acc, $hmm_annot_iso );

    #TIGR roles
    #These are keyed with the .\d\d version numbers. So make a tmp acc and
    #remove this if it's there
    my $tmp = $hmm_acc;
    $tmp =~ s/\.\d+$//;
    $annotation->set_TIGR_Role( $self->_get_tigr_role( $tmp, $com_name ),
                                $hmm_acc, $hmm_annot_iso );
    
}

sub _get_tigr_role {
    my ($self, $hmm_acc, $com_name) = @_;
    my $retval;
    if( $hmm_acc =~ /^PF/ ) {
        my $lookup = $self->_tigr_role_lookup('PFAM');
        $retval = $lookup->pfam2tigr_role( $hmm_acc );
    } elsif( $hmm_acc =~ /TIGR/ ) {
        my $lookup = $self->_tigr_role_lookup('TIGRFam');
        $retval = $lookup->tigrfam2tigr_role( $hmm_acc );
    }

    if( !$retval ) {
        $retval = 157;
    }

    return $retval;
    
}

sub _append_to_com_name {
    my ($self, $com_name, $isotype) = @_;
    
    if( exists( $self->{'_name_suffixes'}->{$isotype} ) ) {
        my $suffix = $self->{'_name_suffixes'}->{$isotype}->{'string'};
        my $action = $self->{'_name_suffixes'}->{$isotype}->{'action'};

        if( $action eq 'append' ) {
            $com_name =~ s/$suffix//;
            $com_name .= $suffix;
        } elsif( $action eq 'replace' ) {
            $com_name = $suffix;
        }
    }

    return $com_name;

}
