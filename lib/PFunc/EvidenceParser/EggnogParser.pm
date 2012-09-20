package PFunc::EvidenceParser::EggnogParser;

use strict;
use warnings;
use XML::Twig;
use Data::Dumper;
use File::OpenFile qw( open_file );
use PFunc::Annotation;

use base qw( PFunc::EvidenceParser );

######################### Class Variables ###########################
my $annotation_type = "eggNOG";
my $percent_id_cutoff = 40;
my $percent_coverage_cutoff = 80;
my $eggnog_annot_levels = {
    'eggNOG::unambiguous'    => 1,   
    'eggNOG::ambiguous'      => 2,   
};

my $flag = 0;
my $sprFlag = 0;
my ($alias, %prot_alias);
my ($members, $description, $fun);
my $nog = "NOG";
my $cog = "COG";

sub new {
    my ($class, %args) = @_;
    my $self = $class->SUPER::new( %args );
    bless($self, $class);

    $self->_init_eggnog_parser(%args);
    return $self;
}

#initialize eggnog parser variables
sub _init_eggnog_parser {
    my ($self, %args ) = @_;
    
    ## initialize instance vars
    $self->{'_seq_lengths'} = {};
    $self->{'_specific_seq_lengths'} = {};
    $self->{'_seq_ids'} = {};
    $self->{'_sequence_titles'} = {};

    if ($args{'alias'} && $args{'fun'} && $args{'description'} && $args{'members'}) {
	$alias = $args{'alias'};
	$fun = $args{'fun'};
	$description = $args{'description'};
	$members = $args{'members'};

    } else {
	die ("Must provide path to eggNOG protein alias file, NOG members file, and NOG description file");
    }

    $self->_set_up_alias_hash($alias);
}

sub _set_up_alias_hash {
    my ($self, $alias) = @_;
    my ($prot, $info);
    my (@line, %desc);

#print "Populating protein alias hash\n";
    #create hash of protein_ids, with EC, gene symbol, COG/NOG IDs, and their gene product names
    my $a = open_file($alias, 'in');
    while (<$a>) {
	chomp;
	($prot, $info) = split(/[|]/);
	$prot =~ s/\d+\.//;
	if ($info =~ /^EC\s+(([\d\-]+\.){3}[\d\-]+)$/i) {
	    $prot_alias{$prot}{'EC'} = $1;
	}
	if ($info =~ /^([a-z]{3}\w*)$/ && !(defined $prot_alias{$prot}{'GS'}) && $info !~ /$prot/ ) {
	    $prot_alias{$prot}{'GS'} = $1;
	}
    }
    close($a);

#print "NOG members/descriptions\n";
    my $nm = open_file($members, 'in');
    my $nd = open_file($description, 'in');
    while (<$nd>) {
	chomp;
	my ($id, $description) = split(/\t/);
	$desc{$id} = $description if defined($description);	#hash with NOG_IDs and descriptions
    }
    close($nd);
    while (<$nm>) {

	chomp;
	@line = split;
	$line[1] =~ s/\d+\.//;	# protein alias ID
	$prot_alias{$line[1]}{'NOG'} = $line[0];	#NOG_ID
	$prot_alias{$line[1]}{'gene_product'} = $desc{$line[0]} if defined($desc{$line[0]});
    }
    close($nm);

#print "COG members/descriptions\n";
    $members =~ s/$nog/$cog/;	#switching to COG members now
    $description =~ s/$nog/$cog/;

    my $cm = open_file($members, 'in');
    my $cd = open_file($description, 'in');
    while (<$cd>) {
	chomp;
	my ($id, $description) = split(/\t/);
	$desc{$id} = $description if defined($description);
    }
    close($cd);
    while (<$cm>) {

	chomp;
	@line = split;
	$prot = $line[1] =~ s/\d+\.//;
	$info = $prot_alias{$line[1]}{'NOG'} = $line[0];
	$prot_alias{$line[1]}{'gene_product'} = $desc{$line[0]} if defined($desc{$line[0]});
    }
    close($cm);
}

#Pre-parsing steps to idenitify initial attributes and definition names
sub _pre_parse {
    my ($self) = @_;
    $self->_parse_sequence_lengths( $self->bsml );
}

sub _parse_sequence_lengths {
  my ($self, $bsmls) = @_;
  my $twig = new XML::Twig('twig_roots' => {'Sequence' => sub {
	my ($t,$e) = @_;
	my $id = $e->att('id');
	my $len = $e->att('length');
	$self->_seq_length( $id, $len )
	if ( defined( $len ) && defined( $id ) && !$self->_seq_length( $id ) );
  	}, 'Feature[@class="polypeptide"]' => sub { $self->_store_feature_length( @_ ) },
  } );

  foreach my $bsml ( @{$bsmls} ) {
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
  $self->_seq_length( $id, $end-$start ) unless( $self->_seq_length( $id ) );
}


#Parsing subroutines to parse out evidence

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
        die("failed to get ID attribute of Sequence in eggNOG alignment file");

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

    ## get the correct annotation id (i.e. if we are annotating on transcript)
    ## we need db1.transcript.123456.1 instead of the CDS id
    my $annotation_feature_id = $self->lookup_feature_id( $ref_id, $self->annotate_on );
    return unless( $annotation_feature_id );

    ## get the annotation object. If nothing is returned, we shouldn't be annoting
    ## this feature.
    my $annotation = $self->get_feature_annotation( $annotation_feature_id );
    return if( !defined( $annotation ) );

    ## make sure the match passes cutoff
    return unless( &_match_passes_cutoff( $spa ) );
    
    ## get the spread coverage
    my ($query_coverage, $subject_coverage) = $self->_calculate_spr_coverage( $ref_id, $comp_id, $spa );

    ## get the annotation related to the compseq and calculate the confidence level
    my $comp_annot = $self->_get_compseq_annotation( $comp_id, $ref_id);
    
    my $gp_name = ($comp_annot->get('gene_product_name'))[0]->[0];
    return unless( defined( $gp_name ) && $gp_name ne "" );

    ##If the match doesn't have any annotation, assign this
    if( !($annotation->has_annotation()) ) {
        $self->_assign_annotation( $annotation, $comp_annot );
    } else {
        
        ## These are the confidence levels
        my $comp_conf_level = $comp_annot->_get_type( 'gene_product_name' );
        my $anno_conf_level = $annotation->_get_type( 'gene_product_name' );
        
        ## otherwise check to see if the comp_annot confidence level 
        ## is better than the current annotation
        if( $eggnog_annot_levels->{ $comp_conf_level } <
            $eggnog_annot_levels->{ $anno_conf_level } ) {
            $self->_assign_annotation( $annotation, $comp_annot );

        ## if they are the same confidence levels
        } elsif( $eggnog_annot_levels->{ $comp_conf_level } ==
                 $eggnog_annot_levels->{ $anno_conf_level } ) {
            
            ## does one have more annotation than the other?
            my ($cur_count, $comp_count);
            foreach my $field( PFunc::Annotation::get_valid_fields ) {
                $cur_count++  if( $annotation->has_annotation( $field ) );
                $comp_count++ if( $comp_annot->has_annotation( $field ) );
            }

            ## if the match has more annotation than current annotation
            ## transfer the annotation
            if( $comp_count > $cur_count ) {
                $self->_assign_annotation( $annotation, $comp_annot );
            }
        }
    }
}

## Checks to see if the name contains one of the following words
## hypothetical probably unknown putative related probable possible conserved
sub _is_name_ambiguous {
    my ($self, $name) = @_;
    my $retval = 0;
    
    my $ambiguous_words_regex = 
        join("|", qw(hypothetical probably unknown putative unrelated) );
    if( $name =~ /($ambiguous_words_regex)/i ) {
        $retval = 1;
    }
    return $retval;
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
## also calculates confidence level
sub _get_compseq_annotation {
    my ($self, $comp_id, $ref_id) = @_;

#print $comp_id, "\n";
    my $gp_name;
# Instantiate each evidence type
    my %assertions = ();
    map {
        $assertions{$_} = '';
    } PFunc::Annotation::get_valid_fields;

#Retrieves the EC, gene symbol, and gene product from the protein alias hash
    $assertions{'EC'} = $prot_alias{$comp_id}{'EC'};
#print "FOUND EC\t $prot_alias{$comp_id}{'EC'}\n" if defined $assertions{'EC'};

    $assertions{'gene_symbol'} = $prot_alias{$comp_id}{'GS'};
#print "FOUND GENE SYMBOL $prot_alias{$comp_id}{'GS'}\n" if defined $assertions{'gene_symbol'};

    $assertions{'gene_product_name'} = $gp_name = $prot_alias{$comp_id}{'gene_product'};
#print "FOUND GENE PRODUCT NAME $gp_name\n" if defined $assertions{'gene_product_name'};

    ## If the match is not characterized and contains vague words ('putative', 'probable', etc.)
    ## then confidence level is set to eggNOG::ambiguous
    my $confidence_level;
    if( !defined($gp_name) || $self->_is_name_ambiguous( $gp_name ) ) {
	$confidence_level = 'eggNOG::ambiguous';
    } else {
	$confidence_level = 'eggNOG::unambiguous';
    }

    ## make the annotation object
    ## just set the feature id to the compseq id.  This will be changed later.
    my $ret_annot = new PFunc::Annotation( 'feature_id' => $comp_id );
    
    foreach my $type ( keys %assertions ) {
        my $field = $type;
#print $field, "\t", $assertions{$type}, "\t", $comp_id, "\t", $confidence_level, "\n";
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


    my ($attribute) = $spa->findnodes('Attribute[@name="percent_identity"]');
    if( $attribute ) {
        my $per_id = $attribute->getAttribute('content');
        $retval = 1 if( $per_id >= $percent_id_cutoff );
    } else {
        die("Could not find percent identity from spr");
    }

    return $retval;
    
}

sub _calculate_spr_coverage {
    my ($self, $query_id, $subject_id, $spa) = @_;
    
    my @sprs = $spa->findnodes( 'Seq-pair-run' );

#    die("There were multiple Seq-pair-run children for this Seq-pair-alignment.")
#        unless( @sprs == 1 );
#    my $spr = $sprs[0];

    my $query_length = $self->_seq_length( $query_id );
    my $subject_length = $self->_specific_seq_length( $subject_id );
	
	die("Could not get length for $query_id") unless( $query_length );	
	die("Could not get length for $subject_id") unless( $subject_length );

    my $query_hit_length = 0;
    my $subject_hit_length = 0;
    foreach my $spr ( @sprs ) {
        $query_hit_length += $spr->getAttribute('runlength');
        $subject_hit_length += $spr->getAttribute('comprunlength');
    }

    my $query_per_cov = int(( ($query_hit_length/$query_length)*10000 )+.5)/100;
    my $subject_per_cov = int(( ($subject_hit_length/$subject_length)*10000 )+.5)/100;
    
    return ($query_per_cov, $subject_per_cov);
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
