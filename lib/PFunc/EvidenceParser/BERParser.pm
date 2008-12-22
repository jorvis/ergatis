package PFunc::EvidenceParser::BERParser;

use strict;
use warnings;
use Carp;
use XML::Twig;
use File::OpenFile qw( open_file );
use Carp qw(cluck);
use Fcntl qw( O_RDONLY );
use MLDBM "DB_File";
use Data::Dumper;
use TIGR::Roles::Omnium::OmniumToRoleLookup;

use base qw(PFunc::EvidenceParser);

######################### Class Variables ###########################
my $annotation_type = "BER";
my $default_ber_info = "/usr/local/projects/db/tchar/tchar.db";
my $percent_id_cutoff = 30;
my $percent_coverage_cutoff = 80;
my $conserved_hypo = {
    'gene_product_name' => "conserved hypothetical protein",
    'GO'                => ['GO:0008150', 'GO:0003674', 'GO:0005575'],
    'TIGR_Role'         => 156,
};
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
    die("No protein id was pssed in") unless( defined( $prot_id ) );
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

	my $trans_id = $self->lookup_feature_id( $id, 'transcript' );
    if( defined( $trans_id ) && $trans_id eq 'gnmM04.transcript.167597423.1' ) {
        $flag = 1;
    }

    ## grab the sequence length. If we already have a length for this sequence don't
    ## set it again.  This is because we've parsed other bsml to grab the real lengths
    ## of the input sequences.  Since CDS sequences are extended by 300 nucs, we grab
    ## this from the gene describing bsml that was passed in.  If this wasn't passed
    ## will use lengths from this file.
    unless( $self->_seq_length( $id ) ) {
        $self->_seq_length( $id, $seq->att('length') ) or 
            die("Could not get length for sequence $id in BER evidence file");
    }
    
    ## check for input sequences.  They will be a CDS. If it is a cds, return
    return if ( $seq->att('class') eq 'CDS' );
    
    my ($defline_att) = $seq->find_nodes('Attribute[@name="defline"]');
    if( defined( $defline_att ) ) {
        
        ## store the look up
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

    print "\n*****Evaluating HSP: $ref_id and $comp_id*****\n";

    ## get the correct annotation id
    my $annotation_feature_id = $self->lookup_feature_id( $ref_id, $self->annotate_on );
    print "returning because I couldn't find an annotation feature_id\n" unless( $annotation_feature_id );
    return unless( $annotation_feature_id );

    ## get the annotation object. If nothing is returned, we shouldn't be annoting
    ## this feature.
    my $annotation = $self->get_feature_annotation( $annotation_feature_id );
    print "returning because no annotation object was returned\n" if( !defined( $annotation ) );
    return if( !defined( $annotation ) );

    #find the best annotation in the header
    my ($id, $title) = $self->_clean_panda_title( $comp_id, $self->_sequence_title( $comp_id ) );
    print "The best annotation for match is from $id [$title]\n";

    # is our match characterized?
    my $compseq_info = $self->_ber_info( $self->_formatted_seq_id( $id ) );

    #disregard the match if the name of the protein is hypothetical protein and it's not characterized
    print "returning because title eq 'hypo' and it's not characterized\n" 
        if( $title eq 'hypothetical protein' && !defined( $compseq_info ) );
    return if( $title eq 'hypothetical protein' && !defined( $compseq_info ) );

    # cycle through the seq-pair-runs
  SPR:
    foreach my $spr ( $spa->children('Seq-pair-run') ) {

        print "skipping this spr because it didn't pass the cutoff\n" unless( &_match_passes_cutoff( $spr ) );
        next SPR unless( &_match_passes_cutoff( $spr ) );

        my ($query_coverage, $subject_coverage) = $self->_calculate_spr_coverage( $ref_id, $comp_id, $spr );
        my $confidence_level = $self->_assign_confidence_level( $query_coverage, $subject_coverage, $compseq_info );
        print "confidence level for this spr is $confidence_level\n";

        print "skipping because confidence level was not acceptable\n"
            if( $confidence_level eq 'BER::uncharacterized::partial::partial' );
        next if( $confidence_level eq 'BER::uncharacterized::partial::partial' );

        ## We are using annotation from only one BER hit. For example, if our top hit
        ## only contains a product name and the next top hit contains an EC number, we shouldn't
        ## annotate the product name from one source and the EC number from another.  This is for
        ## consistency reasons (i.e. the EC number and product name might not match).
        my $source = $self->_formatted_seq_id( $id );
        print "setting new possible source to $source\n";

        # get the old annotation type
        my $annot_type = $annotation->_get_type( 'gene_product_name' );
        if( $annot_type ) {
            print "old annotation type is $annot_type\n";
        } else {
            print "has not yet been assigned an annotation\n";
        }

        if( !$annot_type || 
            $ber_annot_levels->{ $confidence_level } < $ber_annot_levels->{ $annot_type } ) {
            print "this hit has a better confidence level so were replacing the annotation\n";
            print "NEW SOURCE: $id\n";
            $self->_assign_annotation( $annotation, $title, $compseq_info, $source, $confidence_level );

        } elsif( $annot_type && 
                 $ber_annot_levels->{ $confidence_level } == 
                 $ber_annot_levels->{ $annot_type } ) {

            print "the confidence levels were the same\n";
            
            #count how many annotation fields annotation currently has
            my $cur_count = 0;
            print Dumper( $annotation );
            foreach my $field ( Annotation::get_valid_fields ) {
                $cur_count++ if( $annotation->has_annotation($field) );
            }

            #count how many
            my $possible_count = 0;
            $possible_count++ if( $title );
            print Dumper( $compseq_info );
            foreach my $field ( qw(gene_symbol ec_num go) ) {
                $possible_count++ if( exists( $compseq_info->{ $field } ) &&
                                      defined( $compseq_info->{ $field } ) );
            }

            print "current annotation has a count of $cur_count. new annot has $possible_count\n";

            if( $possible_count > $cur_count ) {
                print "since new annot has a higher count, replace annotation\n";
                $self->_assign_annotation( $annotation, $title, $compseq_info, $source, $confidence_level );
                print "NEW SOURCE: $id\n";
            }

        } 

    } #foreach my $spr ( $spa->children('Seq-pair-run') ) {

    if( $annotation->has_annotation() && $self->_is_name_ambiguous( $annotation->_get_value('gene_product_name')->[0] ) ) {
        print "setting annotation to conserved hypothetical because the annotation sucked\n";
        $self->_assign_as_conserved_hypothetical( $annotation, $self->_formatted_seq_id( $comp_id ), "BER::conserved_hypothetical" );
    }

}

sub _clean_panda_title {
    my ($self, $id, $title ) = @_;
    my $retval;

    #append the id back onto the first entry.  This was removed somewhere
    #when converting to bsml from ber btab.
    $title = $self->_formatted_seq_id( $id )." ".$title;

    my @entries = split(/\^\|\^/, $title );

    my $max_count = -1;
    foreach my $entry( @entries ) {
        my @stats = $self->_parse_protein_header_line( $entry );
        my @valid_info = grep {  defined( $_ ) && $_ ne "" && $_ =~ "hypothetical protein" } @stats;
        if( scalar(@valid_info) > $max_count ) {
            $retval = $entry;
            $max_count = scalar(@valid_info);
        }
    }

    return ($1,$2) if( $retval =~ /^(\S+)(.*)/ );
}

sub _parse_protein_header_line {
    my ($self, $line) = @_;
    my ($com_name, $ec, $gene_symbol, $role, $go);
    
    my $full_id = $1 if( $line =~ /^(\S+)\s+.+/ );
    my $look_id = "not in db";

    if( $full_id ) {

        $look_id = $1 if( $full_id =~ /^([^\|]+\|[^\|]+)/ );

        my $info = $self->_ber_info( $look_id );

        #try to parse things from the name
        $com_name = $1 if($line =~ /^\S+\s+(.*)/);

        #clean the common name up a little bit
        $com_name =~ s/taxon:.*//;

        if( defined( $com_name ) && $com_name =~ /\(EC\s+([^\)]+)\)/ ) {
            $ec = $1;
        }

        if( defined( $info ) ) {
            $ec = $info->{'ec_num'} if( exists( $info->{'ec_num'} ) );
            $gene_symbol = $info->{'gene_symbol'} if( exists( $info->{'gene_symbol'}  ) );
            $go = $info->{'go'} if( exists( $info->{'go'} ) );
        }

        $role = $self->_get_tigr_role( $look_id );
        
    }

    print "parsed and now returning $com_name\n" if( $com_name );

    return ($look_id, $com_name, $ec, $gene_symbol, $role, $go);
    
}


# will return a 1 if match passes cutoff
# currently, if percent identity is greater than or equal to
# $percent_id_cutoff
sub _match_passes_cutoff {
    my ($spr) = @_;
    my $retval = 0;

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
    my ($self, $query_id, $subject_id, $spr) = @_;

    my $query_length = $self->_seq_length( $query_id );
    my $subject_length = $self->_seq_length( $subject_id );
    
    my $query_hit_length = $spr->att('runlength');
    my $subject_hit_length = $spr->att('comprunlength');

    my $query_per_cov = int(( ($query_hit_length/$query_length)*10000 )+.5)/100;
    my $subject_per_cov = int(( ($subject_hit_length/$subject_length)*10000 )+.5)/100;
    
    return ($query_per_cov, $subject_per_cov);
}

sub _assign_annotation {
    my ($self, $annotation, $title, $compseq_info, $source, $specific_type) = @_;
    
    if( $specific_type eq 'BER::conserved_hypothetical' ) {
        $self->_assign_as_conserved_hypothetical( $annotation, $source, $specific_type );

    } else {
        
        #clear all annotation
        $annotation->clear_annotation;

        #parse the title
        my ($id, $com_name, $ec, $gene_symbol, $role, $go) = 
            $self->_parse_protein_header_line( $source." ".$title );

        return unless( $com_name );

        #set the new annotation

        #product name
        $annotation->set_gene_product_name( $com_name,
                                            $source, $specific_type );

        #gene symbol
        my $final_gene_sym = $compseq_info->{'gene_symbol'};
        $final_gene_sym = $gene_symbol unless( $final_gene_sym );
        $annotation->set_gene_symbol( $final_gene_sym,
                                      $source, $specific_type );

        #ec number
        my $final_ec = $compseq_info->{'ec_num'};
        $final_ec = $ec unless( $final_ec );
        $annotation->set_EC( $final_ec, 
                             $source, $specific_type );

        # go
        my $final_go = $compseq_info->{'go'};
        $final_go = $go unless( $final_go );
        $annotation->set_GO( $final_go,
                             $source, $specific_type );

        # TIGR roles
        my $final_role = $self->_get_tigr_role( $source );
        $final_role = $role unless( $final_role );
        $annotation->set_TIGR_Role( $final_role,
                                    $source, $specific_type );

        print Dumper( $annotation );
    }
}

# assigns values set in $conserved_hypo global hash
# or empty string if info is not available
sub _assign_as_conserved_hypothetical {
    my ($self, $annotation, $compseq_id, $type ) = @_;
    $annotation->clear_annotation;

    my @fields = Annotation::get_valid_fields;

    foreach my $field ( @fields ) {
        my $value;
        if( exists( $conserved_hypo->{$field} ) ) {
            $value = $conserved_hypo->{$field};
        } else {
            $value = "";
        }
        $annotation->set( $field, $value, $compseq_id, $type );
    }
}
1;
