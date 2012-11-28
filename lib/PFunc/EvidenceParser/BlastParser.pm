package PFunc::EvidenceParser::BlastParser;

use strict;
use warnings;
use Carp;
use XML::Twig;
use XML::LibXML;
use File::OpenFile qw( open_file );
use Carp qw(cluck);
use PFunc::EvidenceParser::ConservedHypothetical;
use PFunc::Annotation;
use UnirefAnnotation::Lookup;
use Data::Dumper;

use base qw(PFunc::EvidenceParser);

######################### Class Variables ###########################
my $annotation_type = "BER";
my $percent_id_cutoff = 40;
my $percent_coverage_cutoff = 80;
my $ber_annot_levels = {
    'BER::characterized::full::full'         => 1,   
    'BER::characterized::full::partial'      => 2,   
    'BER::characterized::partial::full'      => 2,
    'BER::uncharacterized::full::full'       => 3,
    'BER::uncharacterized::full::partial'    => 4,
    'BER::uncharacterized::partial::full'    => 4,   
    'BER::conserved_hypothetical'            => 5,
};
my $count = 0;
######################################################################

sub new {
    my ($class, %args) = @_;
    my $self = $class->SUPER::new( %args );
    bless($self, $class);
    $self->_init_ber_parser( %args );
    return $self;
}
#################################################################
#                    'Private' Subroutines                      #
#################################################################


#######################################################
#                    Init Subs                        #
#######################################################
sub _init_ber_parser {
    my ($self, %args ) = @_;
    
    ## initialize instance vars
    $self->{'_seq_lengths'} = {};
	$self->{'_specific_seq_lengths'} = {};
    $self->{'_seq_ids'} = {};
    $self->{'_sequence_titles'} = {};

    ## make the uniref clusters object    
    if( $args{'char_db'} ) {
        $self->{'_uniref_clusters_annot'} = new UnirefAnnotation::Lookup("path=$args{'char_db'}");
    } elsif( $args{'username'} && $args{'password'} ) {
        $self->{'_uniref_clusters_annot'} = new UnirefAnnotation::Lookup("username=$args{'username'} ".
                                                                         "password=$args{'password'}");
   } else {
       die("Must provide path to characterized db (char_db) or ".
           "username and password for mysql db");
   }

    
}

#######################################################
#                 Pre Parse Subs                      #
#######################################################
sub _pre_parse {
    my ($self) = @_;
    $self->_parse_sequence_lengths( $self->bsml );
}

sub _parse_sequence_lengths {
    my ($self, $bsmls) = @_;
    my $twig = new XML::Twig('twig_roots' => {
	'Sequence' => sub {
	    my ($t,$e) = @_;
	    my $id = $e->att('id');
	    my $len = $e->att('length');
	    $self->_seq_length( $id, $len )
		if ( defined( $len ) && defined( $id ) && !$self->_seq_length( $id ) );
	    $t->purge;
	},
	'Feature[@class="CDS"]' => sub { $self->_store_feature_length( @_ ); $_[0]->purge; },
    } );

    foreach my $bsml ( @{$bsmls} ) {
	print "parsing $bsml\n";
	my $in = open_file( $bsml, 'in' );
	$twig->parse( $in );
	close($in);
    }
}

sub _store_feature_length {
  my ($self, $t, $e) = @_;
  my $id = $e->att('id');
  my $intloc = $e->first_child('Interval-loc');
  my ($start,$end) = ($intloc->att('startpos'), $intloc->att('endpos'));
  $id =~ s/\.CDS$//;
  $self->_seq_length( $id, $end-$start ) unless( $self->_seq_length( $id ) );
}

#######################################################
#                 Parse Subs                          #
#######################################################
sub _parse {
    my ($self, $fh) = @_;

    #print "\n".scalar(keys %{$self->{"_seq_lengths"}} )." keys in length hash\n";

    my $parser = new XML::LibXML;
    my $doc = $parser->parse_fh( $fh );

    foreach my $seq ( $doc->findnodes("//Sequence") ) {
	$self->_handle_sequence( $seq );
    }

    foreach my $spa ( $doc->findnodes("//Seq-pair-alignment") ) {
	$self->_handle_seq_pair_alignment( $spa );
    }

    $self->{'_specific_seq_length'} = {};
    
}

sub _handle_sequence {
    my ($self, $seq) = @_;
    
    # grab the id
    my $id = $seq->getAttribute('id') or 
        die("failed to get ID attribute of Sequence in BER alignment file");

    ## grab the sequence length. If we already have a length for this sequence don't
    ## set it again.  This is because we've parsed other bsml to grab the real lengths
    ## of the input sequences.  Since CDS sequences are extended by 300 nucs, we grab
    ## this from the gene describing bsml that was passed in.  If this wasn't passed
    ## will use lengths from this file.
    unless( $self->_seq_length( $id ) ) {
        $self->_specific_seq_length( $id, $seq->getAttribute('length') ) or 
            return;
    }
}

sub _handle_seq_pair_alignment {
    my ($self, $spa) = @_;
    
    my $ref_id = $spa->getAttribute('refseq');
    my $comp_id = $spa->getAttribute('compseq');

#    print "Comparing $ref_id with $comp_id\n";

    ## get the correct annotation id (i.e. if we are annotating on transcript)
    ## we need db1.transcript.123456.1 instead of the CDS id
    my $annotation_feature_id = $self->lookup_feature_id( $ref_id, $self->annotate_on );
    unless( $annotation_feature_id ) {
	print "\tDidn't find annotation feature id, so skipping\n";
	return;
    }
  #  print "\tFound annotation feature id: $annotation_feature_id\n";

    ## get the annotation object. If nothing is returned, we shouldn't be annoting
    ## this feature.
    my $annotation = $self->get_feature_annotation( $annotation_feature_id );
    if( !defined( $annotation ) ) {
	print "\tCouldn't find annotation object, so not annotating this feature\n";
	return;
    }
 #   print "\tFound the annotation object, so continuing to annotate\n";

    ## make sure the match passes cutoff
    unless( &_match_passes_cutoff( $spa ) ) {
#	print "\tMatch did not pass cutoff, skipping\n";
	return;
    }

    ## Grab the cluster id and see if the match is characterized
    my $db = $self->{'_uniref_clusters_annot'};
    my $comp_trusted = $db->is_trusted( $comp_id );
#     print "\tComp_trusted? ";
#     if( $comp_trusted ) {
# 	print "$comp_trusted";
#     } else {
# 	print "undefind";
#     } 
#     print "\n";
    
    ## assign confidence level
    my ($query_coverage, $subject_coverage) = $self->_calculate_spr_coverage( $ref_id, $comp_id, $spa );
    my $confidence_level = $self->_assign_confidence_level( $query_coverage, $subject_coverage, 
                                                            $comp_trusted );

    #print "\tConfidence level: $confidence_level\n";

    ## we don't use annotation from partial::partial matches.
    if( $confidence_level eq 'BER::uncharacterized::partial::partial' ||
	$confidence_level eq 'BER::characterized::partial::partial' ) {
	#print "\tSkipping because partial::partial\n";
	return;
    }

    ## get the annotation related to the compseq
    my $comp_annot = $self->_get_compseq_annotation( $comp_id, $confidence_level );
    
    #print "\tFound annotation related to comp_seq\n";
    #print Dumper( $comp_annot );

    ## -> Skip matches that contain the words 'hypothetical' in the common name
    ## -> Skip matches that don't have a common name
    ## -> Add the string 'domain protein' to non-hypothetical, full::partial and 
    ##    partial::full matches
    ##       -> Make sure there isn't the word protein already in the common name to 
    ##          avoid names like "zinc finger protein domain protein"
    ## -> If the match is not characterized and contains vague words ('putatative', 'probable', etc.)
    ##    then consider the match to be annotated as conserved hypothetical
    ## -> Then decide on the best using the following criteria
    ##       -> Highest confidence level
    ##       -> If confidence levels are the same, take one with more fields
    
    ## Skipping matches that don't have a common name 
    ## or it contains the word hypothetical
    my $gp_name = ($comp_annot->get('gene_product_name'))[0]->[0];
    unless( defined( $gp_name ) && $gp_name ne "" ) {
	#print "\tReturning because gp name is blank [$gp_name]\n";
	return;
    }
    if( $gp_name =~ /hypothetical/i ) {
	#print "\tReturning because gp name of comp seq has the word hypothetical in it\n";
	return;
    }
    
    ## Add the words 'domain protein' to full::partial matches.
    ## Make sure there isn't the word protein already in the common name to 
    ## avoid names like "zinc finger protein domain protein"
    if( $confidence_level eq 'BER::characterized::full::partial' ||
        $confidence_level eq 'BER::uncharacterized::full::partial' ||
        $confidence_level eq 'BER::characterized::partial::full' ||
        $confidence_level eq 'BER::uncharacterized::partial::full') {

        $gp_name =~ s/\s*protein$//;
        $gp_name .= " domain protein";
        $comp_annot->_set_value( 'gene_product_name', $gp_name );
	#print "\tSetting value of gp name to $gp_name\n";
    }

    ## If the match is not characterized and contains vague words ('putative', 'probable', etc.)
    ## then annotation is changed to conserved hypothetical
    if( $self->_is_name_ambiguous( $gp_name ) && !$comp_trusted ) {
	#print "\tThe name is abmiguous and it's not characterized, assigning as conserved hypo\n";
        $self->_assign_as_conserved_hypothetical( $comp_annot );

	## We add the string 'possible' to the beginning of uncharacterized 
	## matches but not to those which were just annotated as a conserved hypothetical
    } elsif( $confidence_level =~ /BER::uncharacterized/ ) {
	$gp_name =~ s/(putative|possible|probable)\s//gi;
	$gp_name = "Putative ".lc(substr($gp_name,0,1)).substr($gp_name,1);
	$comp_annot->_set_value( 'gene_product_name', $gp_name );
	#print "\tAdding possible because uncharacterized match [$gp_name]\n";
    }

    ##If the match doesn't have any annotation, assign this
    if( !($annotation->has_annotation()) ) {
	#print "\tAssigning annotation because there was none\n";
        $self->_assign_annotation( $annotation, $comp_annot );
    } else {

        ## These are the confidence levels
        my $comp_conf_level = $comp_annot->_get_type( 'gene_product_name' );
        my $anno_conf_level = $annotation->_get_type( 'gene_product_name' );
        
        ## otherwise check to see if the comp_annot confidence level 
        ## is better than the current annotation
        if( $ber_annot_levels->{ $comp_conf_level } <
            $ber_annot_levels->{ $anno_conf_level } ) {
	    #print "\tThe comp annot level was better than current annotation, so replacing\n";
            $self->_assign_annotation( $annotation, $comp_annot );

	## if they are the same confidence levels
        } elsif( $ber_annot_levels->{ $comp_conf_level } ==
                 $ber_annot_levels->{ $anno_conf_level } ) {
            
            ## does one have more annotation than the other?
            my ($cur_count, $comp_count);
            foreach my $field( PFunc::Annotation::get_valid_fields ) {
                $cur_count++  if( $annotation->has_annotation( $field ) );
                $comp_count++ if( $comp_annot->has_annotation( $field ) );
            }

            ## if the match has more annotation than current annotation
            ## transfer the annotation
            if( $comp_count > $cur_count ) {
		#print "\tThe comp annotation has same level, but was better anyway, so replacing\n";
                $self->_assign_annotation( $annotation, $comp_annot );
            }
        }
    }

    my $annot = $self->get_feature_annotation( $annotation_feature_id );
    
}

## Checks to see if the name contains one of the following words
## hypothetical probably unknown putative related probable possible conserved
sub _is_name_ambiguous {
    my ($self, $name) = @_;
    my $retval = 0;
    
    my $ambiguous_words_regex = 
        join("|", qw(hypothetical probably unknown putative related probable possible conserved) );
    if( $name =~ /($ambiguous_words_regex)/ ) {
        $retval = 1;
    }

    if( !$retval && length( $name ) < 4 ) {
        $retval = 1;
    }
    return $retval;
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

## Will assign some annotation to some other annotation.
sub _assign_annotation {
    my ($self, $annotation, $new_annotation) = @_;

    my $feature_id = $annotation->get_feature_id;
    $new_annotation->set_feature_id( $feature_id );
    $annotation->transfer_annotation( $new_annotation );

}

## Grabs the annotation from the database and creates
## an annotation object from it.
sub _get_compseq_annotation {
    my ($self, $comp_id, $confidence_level) = @_;

    my $db = $self->{'_uniref_clusters_annot'};
    my $assert_array = $db->get_annotation( $comp_id );
    my %assertions = ();

    ## make the annotation object
    ## just set the feature id to the compseq id.  This will be changed later.
    my $ret_annot = new PFunc::Annotation( 'feature_id' => $comp_id );

    map {
        
        #initialize to empty array ref
        $_->{'type'} = 'TIGR_Role' if( $_->{'type'} eq 'TIGRrole' );
        $assertions{$_->{'type'}} = [] unless( exists $assertions{$_->{'type'}} );

        #push on values.
        push( @{$assertions{$_->{'type'}}}, $_->{'value'} );
        
    } @{$assert_array};
    
    foreach my $type ( keys %assertions ) {
        my $field = $type;
        $field = "gene_product_name" if( $type eq 'common_name' );
        $ret_annot->set( $field, $assertions{$type}, $comp_id, $confidence_level );
    }

    return $ret_annot;
    
}

# will return a 1 if match passes cutoff
# currently, if percent identity is greater than or equal to
# $percent_id_cutoff
sub _match_passes_cutoff {
    my ($spa) = @_;
    my $retval = 0;

    my ($att) = $spa->findnodes('Attribute[@name="percent_identity"]');
    if( $att ) {
	my $per_id = $att->getAttribute('content');
	$retval = 1 if( $per_id >= $percent_id_cutoff );
    } else {
	die("Could not parse percent_identity attribute from spa");
    }

    return $retval;
}

sub _calculate_spr_coverage {
    my ($self, $query_id, $subject_id, $spa) = @_;

    my ($pc_query_att)   = $spa->findnodes('Attribute[@name="percent_coverage_refseq"]');
    my ($pc_subject_att) = $spa->findnodes('Attribute[@name="percent_coverage_compseq"]');

    my $pc_query   = $pc_query_att->getAttribute('content');
    my $pc_subject = $pc_subject_att->getAttribute('content');

    return ($pc_query, $pc_subject);
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

sub _specific_seq_length {
  my ($self, $seq_id, $length) = @_;
  my $retval;
  if( defined($length) ) {
	$self->{'_specific_seq_lengths'}->{$seq_id} = $length;
	$retval = $length;
  } elsif( exists( $self->{'_specific_seq_lengths'}->{$seq_id} ) ) {
	$retval = $self->{'_specific_seq_lengths'}->{$seq_id};
  }
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
1;
