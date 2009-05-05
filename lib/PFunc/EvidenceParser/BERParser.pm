package PFunc::EvidenceParser::BERParser;

use strict;
use warnings;
use Carp;
use XML::Twig;
use File::OpenFile qw( open_file );
use Carp qw(cluck);
use Fcntl qw( O_RDONLY );
use MLDBM "DB_File";
use TIGR::Roles::Omnium::OmniumToRoleLookup;
use PFunc::EvidenceParser::ConservedHypothetical;
use PFunc::Annotation;
use Data::Dumper;

use base qw(PFunc::EvidenceParser);

######################### Class Variables ###########################
my $annotation_type = "BER";
my $default_ber_info = "/usr/local/projects/db/tchar/tchar.db";
my $percent_id_cutoff = 30;
my $percent_coverage_cutoff = 80;
my $ber_annot_levels = {
    'BER::characterized::full::full'         => 1,   
    'BER::characterized::full::partial'      => 2,   
    'BER::characterized::partial::full'      => 2,
    'BER::uncharacterized::full::full'       => 3,
    'BER::uncharacterized::full::partial'    => 4,
    'BER::uncharacterized::partial::full'    => 4,   
    'BER::characterized::partial::partial'   => 5,
    'BER::uncharacterized::partial::partial' => 6,
    'BER::conserved_hypothetical'            => 7,
};

my $flag = 0;
my $sprFlag = 0;
######################################################################

sub new {
    my ($class, %args) = @_;
    my $self = $class->SUPER::new( %args );
    bless($self, $class);
    $self->_init_ber_parser( %args );
    return $self;
}

##### private subroutines #####
sub _init_ber_parser {
    my ($self, %args ) = @_;
    
    ## initialize instance vars
    $self->{'_seq_lengths'} = {};
    $self->{'_seq_ids'} = {};
    $self->{'_sequence_titles'} = {};

    ## ber info
    my %tmp;
    if( $args{'ber_info'} ) {
        tie(%tmp, 'MLDBM', $args{'ber_info'}, O_RDONLY )
            or die("Could not tie hash to $args{'ber_info'}");
    } else {
        tie(%tmp, 'MLDBM', $default_ber_info, O_RDONLY )
            or die("Could not tie hash to $args{'ber_info'}");
    }
    $self->{'_ber_info'} = \%tmp;

    ## tigr roles lookup
    $self->_tigr_roles_lookup( new TIGR::Roles::Omnium::OmniumToRoleLookup() );
    
    
}

sub _parse {
    my ($self, $fh) = @_;

    my $twig = new XML::Twig( 'twig_handlers' => {
        'Sequence' => sub { $self->_handle_sequence(@_) },
        'Seq-pair-alignment' => sub { $self->_handle_seq_pair_alignment(@_) },
    });

    $twig->parse($fh);
}

sub _pre_parse {
    my ($self) = @_;
    $self->_parse_sequence_lengths( $self->bsml );
}

sub _parse_sequence_lengths {
    my ($self, $bsmls) = @_;
    my $twig = new XML::Twig( 'twig_roots' => {
        'Sequence' => sub {
            my ($t,$e) = @_;
            my $id = $e->att('id');
            my $len = $e->att('length');
            $self->_seq_length( $id, $len )
                if( defined( $len ) && defined( $id ) );
        },
        'Feature' => sub {
            my($t, $e) = @_;
            my $id = $e->att('id');
            my $intloc = $e->first_child('Interval-loc');
            my ($start,$end) = ($intloc->att('startpos'), $intloc->att('endpos'));
            $self->_seq_length( $id, $end-$start );
        }
    } );

    foreach my $bsml( @{$bsmls} ) {
        my $in = open_file( $bsml, 'in' );
        $twig->parse( $in );
        close($in);
    }
}

sub _tigr_roles_lookup {
    my ($self, $lookup)  = @_;
    my $retval;
    if( $lookup ) {
        $self->{'_tigr_roles_lookup'} = $lookup;
    }
    $retval = $self->{'_tigr_roles_lookup'};
    return $retval;
}

sub _get_tigr_role {
    my ($self, $prot_id) = @_;
    die("No protein id was passed in") unless( defined( $prot_id ) );
    my $lookup = $self->_tigr_roles_lookup();
    my $retval = $lookup->omnium2tigr_role( $prot_id );
    return $retval;
}

sub _seq_length {
    my ($self, $seq_id, $length) = @_;
    my $retval;
    if( defined($length) ) {
        $self->{'_seq_lengths'}->{$seq_id} = $length;
        $retval = $length;
    } elsif( exists( $self->{'_seq_lengths'}->{$seq_id} ) ) {
        $retval = $self->{'_seq_lengths'}->{$seq_id};
    }
    return $retval;
}

sub _formatted_seq_id {
    my ($self, $bsml_id, $real_id) = @_;
    my $retval = $bsml_id;
    if( defined($real_id) ) {
        $self->{'_seq_ids'}->{$bsml_id} = $real_id;
        $retval = $real_id;
    } elsif( exists( $self->{'_seq_ids'}->{$bsml_id} ) ) {
        $retval = $self->{'_seq_ids'}->{$bsml_id};
    }
    return $retval;
}

sub _sequence_title {
    my ($self, $bsml_id, $title) = @_;
    my $retval;
    if( defined( $title ) ) {
        $self->{'_sequence_titles'}->{$bsml_id} = $title;
        $retval = $title;
    } elsif( exists( $self->{'_sequence_titles'}->{$bsml_id} ) ) {
        $retval = $self->{'_sequence_titles'}->{$bsml_id};
    }
    return $retval;
}

sub _ber_info {
    my ($self, $id) = @_;
    my $retval;
    if( !defined( $id ) ) {
        $retval = $self->{'_ber_info'};
    } elsif( exists( $self->{'_ber_info'}->{$id} ) ) {
        $retval = $self->{'_ber_info'}->{$id};
    }
    return $retval;
}

sub _handle_sequence {
    my ($self, $twig, $seq) = @_;
    
    # grab the id
    my $id = $seq->att('id') or 
        die("failed to get ID attribute of Sequence in BER alignment file");

    ## grab the sequence length. If we already have a length for this sequence don't
    ## set it again.  This is because we've parsed other bsml to grab the real lengths
    ## of the input sequences.  Since CDS sequences are extended by 300 nucs, we grab
    ## this from the gene describing bsml that was passed in.  If this wasn't passed
    ## will use lengths from this file.
    unless( $self->_seq_length( $id ) ) {
        $self->_seq_length( $id, $seq->att('length') ) or 
            return;
    }
    
    ## check for input sequences.  They will be a CDS. If it is a cds, return
    return if ( $seq->att('class') eq 'CDS' );
    
    my ($defline_att) = $seq->find_nodes('Attribute[@name="defline"]');
    if( defined( $defline_att ) ) {
        
        ## store the formatted id
        my $formatted_id = $defline_att->att('content');
        $formatted_id = $1 if( $formatted_id =~ /^([A-Z]+\|[A-Z0-9_\.]+)/ );
        $self->_formatted_seq_id( $id, $formatted_id );

    } else {
        # if it didn't have a defline, not sure what it is. Skip it.
        return;
    }
    
    ## store the panda header
    my $title = $seq->att('title');
    $self->_sequence_title( $id, $title );
}

sub _handle_seq_pair_alignment {
    my ($self, $twig, $spa) = @_;

    my $ref_id = $spa->att('refseq');
    my $comp_id = $spa->att('compseq');

    ## get the correct annotation id
    my $annotation_feature_id = $self->lookup_feature_id( $ref_id, $self->annotate_on );
    return unless( $annotation_feature_id );

    ## get the annotation object. If nothing is returned, we shouldn't be annoting
    ## this feature.
    my $annotation = $self->get_feature_annotation( $annotation_feature_id );
    return if( !defined( $annotation ) );

    my ($query_coverage, $subject_coverage) = $self->_calculate_spr_coverage( $ref_id, $comp_id, $spa );
    my $confidence_level = $self->_assign_confidence_level( $query_coverage, $subject_coverage, 
                                                            $self->_ber_info( $self->_formatted_seq_id( $comp_id ) ) );

    ## get the annotation related to the compseq
    my ($comp_annot, $is_char) = $self->_get_compseq_annotation( $comp_id, $confidence_level );

    # don't take annotation from proteins containing the words hypothetical protein
    my $gp_name = $comp_annot->_get_value( 'gene_product_name' )->[0];
    return if( $gp_name =~ /hypothetical\s+protein/ );

    # if we don't have a common name, skip it. 
    return if( !defined($gp_name) || $gp_name eq "" );

    # does it pass cutoffs?
    return unless( &_match_passes_cutoff( $spa ) );

    # if we've made it this far and the name is ambiguous (contains general terms such as
    # putative, probable, etc.) and the match protein is not characterized, we should mark 
    # the protein as conserved hypothetical.
    if( $self->_is_name_ambiguous( $gp_name ) && !$is_char ) {
        $self->_assign_as_conserved_hypothetical( $comp_annot );
    } else {
        # if the current protein hasn't been annotated yet, automatically assign the current annotation
        if( !($annotation->has_annotation()) ) {
            $self->_assign_annotation( $annotation, $comp_annot );
        } else {

            my $comp_type = $comp_annot->_get_type( 'gene_product_name' );
            my $anno_type = $annotation->_get_type( 'gene_product_name' );

            #otherwise check to see if the comp_annot confidence level is better than the current annotation
            if( $ber_annot_levels->{ $comp_annot->_get_type( 'gene_product_name' ) } <
                $ber_annot_levels->{ $annotation->_get_type( 'gene_product_name' ) } ) {
                $self->_assign_annotation( $annotation, $comp_annot );
                #if they are the same
            } elsif( $ber_annot_levels->{ $comp_annot->_get_type( 'gene_product_name' ) } ==
                     $ber_annot_levels->{ $annotation->_get_type( 'gene_product_name' ) } ) {
                
                #does one have more annotation than the other?
                my ($cur_count, $comp_count);
                foreach my $field( PFunc::Annotation::get_valid_fields ) {
                    $cur_count++  if( $annotation->has_annotation( $field ) );
                    $comp_count++ if( $comp_annot->has_annotation( $field ) );
                }

                if( $comp_count > $cur_count ) {
                    $self->_assign_annotation( $annotation, $comp_annot );
                }
            }
        }
    }
}

sub _get_compseq_annotation {
    my ($self, $id, $confidence_level) = @_;
    my $formatted_id = $self->_formatted_seq_id( $id );

    # what infomration can we get from the protein header?
    my @header_info = $self->_clean_panda_title( $id, $self->_sequence_title( $id ) );

    # make the annotation object
    my $source = shift @header_info;
    # just set the feature id to the compseq id.  This will be changed later.
    my $ret_annot = new PFunc::Annotation( 'feature_id' => $source );

    foreach my $field ( qw(gene_product_name EC gene_symbol TIGR_Role GO) ) {
        $ret_annot->set( $field, shift @header_info, $source, $confidence_level );
    }

    if( my $com_name = $ret_annot->_get_value( 'gene_product_name' )->[0] ) {
        if( $com_name !~ /hypothetical/ ) {
            my $tigr_roles = $ret_annot->_get_value( 'TIGR_Role' );
            if( @{$tigr_roles} == 0 ) {
                $ret_annot->set( 'TIGR_Role', [157], $source, $confidence_level );
            }
        }
    }

    return ($ret_annot, shift @header_info );
    
}

sub _clean_panda_title {
    my ($self, $id, $title ) = @_;
    my @retval;

    #append the id back onto the first entry.  This was removed somewhere
    #when converting to bsml from ber btab.
    $title = $self->_formatted_seq_id( $id )." ".$title;

    my @entries = split(/\^\|\^/, $title );

    my $max_count = -1;
    foreach my $entry( @entries ) {
        my @stats = $self->_parse_protein_header_line( $entry );
        my @valid_info = grep {  defined( $_ ) && $_ ne "" && $_ !~ /hypothetical protein/ } @stats;
        if( scalar(@valid_info) > $max_count ) {
            @retval = @stats;
            $max_count = scalar(@valid_info);
        }
    }

    return @retval;
}

sub _parse_protein_header_line {
    my ($self, $line) = @_;
    my ($com_name, $ec, $gene_symbol, $role, $go);
    my $is_char = 0;
    
    my $full_id = $1 if( $line =~ /^(\S+)\s+.+/ );
    my $look_id = "not in db";

    if( $full_id ) {

        $look_id = $1 if( $full_id =~ /^([^\|]+\|[^\|]+)/ );

        my $info = $self->_ber_info( $look_id );

        #try to parse things from the name
        $com_name = $1 if($line =~ /^\S+\s+(.*)/);

        #clean the common name up a little bit
        $com_name =~ s/taxon:.*//;
        $com_name =~ s/^\s+//;
        $com_name =~ s/\s+$//;

        if( defined( $com_name ) ) {
            while( $com_name =~ /\(EC\s+([^\)]+)\)/g ) {
                push( @{$ec}, $1 );
            }
        }

        if( defined( $info ) ) {
            $is_char = 1;
            @{$ec} = split(/[,\s]+/, $info->{'ec_num'}) if( exists( $info->{'ec_num'} ) );
            map { $_ =~ s/^EC\s+// } @{$ec};
            @{$ec} = grep( $_ !~ /^\s*$/, @{$ec} );
            $gene_symbol = $info->{'gene_symbol'} if( exists( $info->{'gene_symbol'} ) && 
                                                      $info->{'gene_symbol'} !~ /^\s*$/ );
            $go = $info->{'go'} if( exists( $info->{'go'} ) );

            my $tmp_com_name = $info->{'com_name'};

            #if it's characterized, we usually don't use the name in the lookup. But if we
            #have a crappy name from the header, we can use it.
            if( $self->_is_name_ambiguous( $com_name ) && !$self->_is_name_ambiguous( $tmp_com_name ) ) {
                $com_name = $tmp_com_name;
            }

        }

        $role = $self->_get_tigr_role( $look_id );
        
    }

    return ($look_id, $com_name, $ec, $gene_symbol, $role, $go, $is_char);
    
}


# will return a 1 if match passes cutoff
# currently, if percent identity is greater than or equal to
# $percent_id_cutoff
sub _match_passes_cutoff {
    my ($spa) = @_;
    my $retval = 0;

    # should only have 1 spr
    my @sprs = $spa->children('Seq-pair-run');
    die("SPA had more than one SPR.  Is this supposed to happen with BER results?")
        unless( @sprs == 1 );
    my $spr = shift @sprs;

    my ($attribute) = $spr->find_nodes('Attribute[@name="percent_identity"]');
    if( $attribute ) {
        my $per_id = $attribute->att('content');
        $retval = 1 if( $per_id >= $percent_id_cutoff );
    } else {
        die("Could not find percent identity from spr");
    }

    return $retval;
    
}

sub _is_name_ambiguous {
    my ($self, $name) = @_;
    my $retval = 0;
    
    my @ambiguous_words = qw( hypothetical probably unknown putative related );
    foreach my $aword ( @ambiguous_words ) {
        if( $name =~ /$aword/ ) {
            $retval = 1;
            last;
        }
    }
    return $retval;
}

sub _assign_confidence_level {
    my ($self, $q_cov, $s_cov, $char) = @_;
    my $retval = "BER";

    if( $char ) {
        $retval .= "::characterized";
    } else {
        $retval .= "::uncharacterized";
    }

    if( $q_cov >= $percent_coverage_cutoff ) {
        $retval .= "::full";
    } else {
        $retval .= "::partial";
    }
    
    if( $s_cov >= $percent_coverage_cutoff ) {
        $retval .= "::full";
    } else {
        $retval .= "::partial";
    }
    
    return $retval;
}

sub _calculate_spr_coverage {
    my ($self, $query_id, $subject_id, $spa) = @_;
    
    my @sprs = $spa->children( 'Seq-pair-run' );

    #I'm pretty sure ber will only ever have one of these
    die("There were multiple Seq-pair-run children for this Seq-pair-alignment.")
        unless( @sprs == 1 );
    my $spr = $sprs[0];

    my $query_length = $self->_seq_length( $query_id );
    my $subject_length = $self->_seq_length( $subject_id );
    
    my $query_hit_length = $spr->att('runlength');
    my $subject_hit_length = $spr->att('comprunlength');

    my $query_per_cov = int(( ($query_hit_length/$query_length)*10000 )+.5)/100;
    my $subject_per_cov = int(( ($subject_hit_length/$subject_length)*10000 )+.5)/100;
    
    return ($query_per_cov, $subject_per_cov);
}

sub _assign_annotation {
    my ($self, $annotation, $new_annotation) = @_;

    my $feature_id = $annotation->get_feature_id;
    $new_annotation->set_feature_id( $feature_id );
    $annotation->transfer_annotation( $new_annotation );

}

# assigns values set in $conserved_hypo global hash
# or empty string if info is not available
sub _assign_as_conserved_hypothetical {
    my ($self, $annotation ) = @_;
    $annotation->clear_annotation;

    my $source = $annotation->get_feature_id;
    my $type = "BER::conserved_hypothetical";

    PFunc::EvidenceParser::ConservedHypothetical->_assign_as_conserved_hypothetical( $annotation, $source, $type );

}
1;
