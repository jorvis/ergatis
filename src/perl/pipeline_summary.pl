#!/usr/bin/perl

=head1 NAME

pipeline_summary.pl - Creates one BSML document describing annotation from a pipeline 
    from other bsml files.  Combines all the features.

=head1 SYNOPSIS

USAGE: pipeline_summary.pl
            --input_bsml=/path/to/file.bsml
            --other_bsml_lists=/path/to/some.list,/path/to/another.list
            --output=/path/to/someDir/
          [ --locus_prefix=ABO553
            --organism="Acidolophus borneisi"
            --translation_table=11
            --cog_search_bsml=/path/to/wu-blastp.bsml.list
            --cog_lookup=/path/to/COGS_file.info
            --cds_fasta=/path/to/CDS.fasta|CDS.fsa.list
            --polypeptide_fasta=/path/to/polypeptide.fasta|polypeptide.fsa.list
            --log=/path/to/file.log
            --debug=4
            --help
          ]


=head1 OPTIONS

B<--input_bsml,-i>
    [Required] Curated gene set bsml document (auto annotate output) (either a file or a list of files).
    

B<--other_bsml_lists,-b>
    [Optional] Comma separated list of bsml lists.  Any other list of 
    bsml documents that define features.  Will search for the correct
    file based on sequence ids.  ie:

    mc74.assembly.7.glimmer3.bsml
                 AND
    mc74.assembly.7.tRNAscan-SE.bsml

B<--output_dir,-o>
    [Required] Output directory where bsml will go.

B<--locus_prefix,-u>
    [Optional] The prefix to be used in Feature/Attribute@name="locus" elements.

B<--organism,-r>
    [Optional] Must have at least two words separated by spaces.  The first word will be
     the genus and the rest of the string will be entered as the species.

B<--translation_table,-t>
    [Optional] Default: 11.  The translation table used for gene prediction.

B<--cog_search_bsml,-c>
    [Optional] Bsml file containing a blast against NCBI's COGs database.

B<--cog_lookup,-g>
    [Optional] NCBI Cogs file

B<--log,-l>
    [Optional] Logfile.

B<--debug,-d>
    [Optional] Higher number is more verbose.

B<--help,-h>
    [Optional] Print this message

=head1  DESCRIPTION

=head1  INPUT


=head1 OUTPUT

    A file will be created in the output directory.  The name of this file will be based on the auto_annotate
    input.


=head1  CONTACT

    Kevin Galens
    kgalens@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Ergatis::Logger;
use Ergatis::IdGenerator;
use BSML::BsmlBuilder;
use XML::Twig;
use File::OpenFile qw(open_file);
use Data::Dumper;
$|++;

####### GLOBALS AND CONSTANTS ###########
my $debug = 0;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);

my @input_files;
my @other_files;
my $output_directory;
my $locus_prefix;
my $locus_database = 'TIGR_moore';
my $organism;
my $trans_table = 11;
my $logger;

#cog related info
my $include_cog = 0;
my %cogId2desc;
my %top_cogs;

#for parsed fasta
my %seq_data_import_info;

#parsed data
my %input_sequences;

#convenience hash
my %feature_index;

#for the output bsml analysis element
my $analysis_name = "pipeline_summary_analysis";
my $sourcename;

########################################

my %options = ();
my $results = GetOptions(\%options, 
                          'input_bsml|i=s',
                          'other_bsml_lists|b=s',
                          'output|o=s',
                          'locus_prefix|u=s',
                          'locus_database|U=s',
                          'organism|r=s',
                          'translation_table|t=s',
                          'cog_search_bsml|s=s',
                          'cog_lookup|g=s',
                          'cds_fasta|f=s',
                          'polypeptide_fasta|p=s',
                          'analysis_name=s',
                          'sourcename=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || &_pod;

# Check the options.
&check_parameters(\%options);

# We start by processing the input files
foreach my $input_file( @input_files ) {
    &parse_input_bsml_file( $input_file );
}

#Next we parse the other bsml files
my $flag = 0;
foreach my $other_file ( @other_files ) {
    $flag = 1;
    &parse_input_bsml_file( $other_file );
}

#Store cogs if we need to
if( $include_cog ) {
    &add_cogs;
}

#If we included fasta 
if( scalar( keys %seq_data_import_info ) > 0 ) {
    &add_seq_data_import;
}

#if we included a locus prefix, assign loci
if( $locus_prefix ) {
    &assign_loci;
}

&to_bsml;

######################## SUB ROUTINES #######################################
sub parse_input_bsml_file {
    my ($input_file) = @_;
    &_log("Parsing $input_file\n", $DEBUG );
    
    my @sequences;

    my $twig = new XML::Twig( 'twig_roots' => {
        'Sequence' => \&input_sequence_handler,
    });

    my $oh = open_file( $input_file, "in" );
    $twig->parse( $oh );
}

sub input_sequence_handler {
    my ($twig, $seq_elem) = @_;
    
    #get the sequence id
    my $seq_id = $seq_elem->att('id');
    my $features = &get_sequence_features( $seq_elem );
    my $atts = &get_attributes( $seq_elem );
    my $bsml_atts = &get_bsml_attributes( $seq_elem );
    my $links = &get_bsml_links( $seq_elem );
    my $sdi = &get_seq_data_import( $seq_elem );
    my $title = &get_seq_title( $seq_id, $sdi );
    $atts->{'title'} = $title if( defined( $title ) );

    if( exists( $input_sequences{ $seq_id } ) ) {
        &append_data( $seq_id, { 
            'features' => $features, 
            'attributes' => $atts,
            'bsml_attributes', => $bsml_atts,
            'links' => $links,
            'seq-data-import' => $sdi
            });
    } else {
        $input_sequences{$seq_id} = {
            'features' => $features,
            'attributes' => $atts,
            'bsml_attributes' => $bsml_atts,
            'links' => $links,
            'seq-data-import' => $sdi
            };
    }
}

sub add_cogs {
    &_log("Adding COG information\n", $DEBUG );

    foreach my $seq_id ( keys %input_sequences ) {
        foreach my $group_set ( keys %{$input_sequences{$seq_id}->{'features'}} ) {
            my $fg = $input_sequences{$seq_id}->{'features'}->{$group_set};
            
            #if we don't have a polypeptide object for this group-set, this means
            #that it's a feature which doesn't have a translation (i.e. ncRNA) and
            #we wouldn't have performed a COG analysis on it.
            next unless( exists( $fg->{'polypeptide'} ) );

            my $cog_desc = &get_cog_desc( $fg->{'polypeptide'}->{'attributes'}->{'id'} );
            $fg->{'CDS'}->{'bsml_attributes'}->{'top_cog_hit'} = $cog_desc if( $cog_desc );
        }
    }
}

sub add_seq_data_import { 
    #this subroutine verifies that we only keep sequences that were 
    #parsed from the bsml.
    &_log("Adding sequence data from input fasta sequences\n", $DEBUG);
    
    foreach my $seq_id ( keys %seq_data_import_info ) {
        next unless( exists( $feature_index{ $seq_id } ) );
        my $href = $seq_data_import_info{$seq_id};
        $input_sequences{ $seq_id } = {
            'seq-data-import' => {
                'format' => $href->{'format'},
                'source' => $href->{'source'},
                'identifier' => $href->{'identifier'}
            },
            'attributes' => {
                'id' => $seq_id."_seq",
                'title' => $seq_id,
            }
        };
    }
}

sub assign_loci {
    #we want to order the sequences based on length.  And we only want to order the sequences
    #which contain features.
    my @seqs_with_features = grep { exists( $input_sequences{$_}->{'features'} ) } keys %input_sequences;
    my @ordered_seq_ids = sort by_molecule_length @seqs_with_features;

    my $locus_num = 1;
    foreach my $seq_id ( @ordered_seq_ids ) {
        foreach my $group_set ( sort by_end3 values %{$input_sequences{$seq_id}->{'features'}} ) {
            my $gene = $group_set->{'gene'};
            my $tmp = {
                'database' => $locus_database,
                'identifier' => $locus_prefix."_".$locus_num++,
                'identifier-type' => 'locus',
            };
            push( @{$gene->{'cross-references'}}, $tmp );
        }
    }
}

sub by_end3 {
    my $end3_a = ( $a->{'gene'}->{'interval_loc'}->{'complement'} ) ?
        $a->{'gene'}->{'interval_loc'}->{'startpos'} :
        $a->{'gene'}->{'interval_loc'}->{'endpos'};
    
    my $end3_b = ( $b->{'gene'}->{'interval_loc'}->{'complement'} ) ?
        $b->{'gene'}->{'interval_loc'}->{'startpos'} :
        $b->{'gene'}->{'interval_loc'}->{'endpos'};

    return $end3_a <=> $end3_b;
}

sub by_molecule_length {
    if( !exists( $input_sequences{$a}->{'attributes'}->{'length'} ) ) {
        &_log("Could not get length for sequence $a", $ERROR );
    }
    if( !exists( $input_sequences{$b}->{'attributes'}->{'length'} ) ) {
        &_log("Could not get length for sequence $b", $ERROR );
    }
    $input_sequences{$b}->{'attributes'}->{'length'} <=> $input_sequences{$a}->{'attributes'}->{'length'};
}

#### Parsing Subroutines #################
sub get_seq_data_import {
    my ($elem) = @_;
    my $atts = $elem->first_child('Seq-data-import')->atts();
    delete( $atts->{'id'} );
    return $atts;
}

sub get_seq_title {
    my ($seq_id, $sdi) = @_;
    my $retval;
    
    #find and open the fasta file
    my $fh = open_file( $sdi->{'source'}, 'in' );
    my $identifier = $sdi->{'identifier'};

    #find the sequence in the file
    while( my $line = <$fh> ) {
        if( $line =~ /^>(\S+)\s+(\S+)/ ) {
            if( $identifier eq $1 ) {
                $retval = $2;
                last;
            }
        }
    }
    
    return $retval;
}

#returns a hashref of feature groups keyed on gene id
#  feature_id -> {
#     'attributes' => { class => 'gene'...},
#     'bsml_attributes' => ( { name => 'whatever', content => 'whatever' }, ... ),
#     'links' => ( { sequence='polypeptide.1_seq' } ),
#     'attribute_lists => ( bsml_attributes )
#     'cross-reference' => ( { 'database' => 'TIGR_moore', 'identifier' => "a124",
#                              'identifier_type' => 'locus' } )
sub get_sequence_features {
    my ($elem) = @_;
    my %retval;

    my $seq_id = $elem->att('id');

    my $feature_tables = $elem->first_child( 'Feature-tables' );
    return {} if( !defined( $feature_tables ) );

    my @features = $feature_tables->first_child('Feature-table')->children('Feature');
    my $total_features = scalar(@features);
    my $count = 0;
    &_log("There are $total_features features on sequence $seq_id\n", $DEBUG );
    return () if( $total_features == 0 );
    
    foreach my $feature ( @features ) {
        my $atts = &get_attributes( $feature );
        my $bsml_atts = &get_bsml_attributes( $feature );
        my $bsml_links = &get_bsml_links( $feature );
        my $interval_loc = &get_interval_loc( $feature );
        my @attribute_lists = &get_attribute_lists( $feature );
        my @cross_references = &get_cross_references( $feature );

        $feature_index{ $feature->att('id') } = {
            'attributes' => $atts,
            'bsml_attributes' => $bsml_atts,
            'links' => $bsml_links,
            'attribute_lists' => \@attribute_lists,
            'interval_loc' => $interval_loc,
            'cross-references' => \@cross_references
        };

        $count++;
        &_log("$count/$total_features\r", $DEBUG );
    }

    my @feature_groups = $elem->find_nodes('Feature-tables/Feature-group');
    my $total_fgs = scalar(@feature_groups);
    $count = 0;
    &_log("\nThere are $total_fgs feature groups\n", $DEBUG );

    foreach my $feature_group ( @feature_groups ) {
        my $tmp = {};
        my $group_set = $feature_group->att('group-set');
        foreach my $fgm ( $feature_group->children( 'Feature-group-member' ) ) {
            my $id = $fgm->att('featref');
            my $class = $fgm->att('feature-type');
            $tmp->{$class} = $feature_index{ $id };
        }
        $count++;
        &_log("$count/$total_fgs\r", $DEBUG );
        $retval{$group_set} = $tmp;
    }
    &_log("\n\n", $DEBUG );
    return \%retval;
}

sub get_cross_references {
    my ($elem) = @_;
    my @retval;
    map { my $atts = $_->atts; delete( $atts->{'id'} ); push(@retval, $atts); } 
    $elem->children('Cross-reference');
    return @retval;
}

sub get_attribute_lists {
    my ($elem) = @_;
    my @retval;
    map { push(@retval, &get_bsml_attributes( $_ ) ) } 
    $elem->children('Attribute-list');
    return @retval;
}

sub get_interval_loc {
    my ($elem) = @_;
    my $retval = {};
    my $int_loc = $elem->first_child('Interval-loc');
    $retval = $int_loc->atts() if( $int_loc );
    return $retval;
}

sub get_attributes {
    my ($elem) = @_;
    return $elem->atts();
}

sub get_bsml_attributes {
    my ($elem) = @_;
    my %retval;

    map {
        $retval{ $_->att('name') } = $_->att('content');
    } $elem->children( 'Attribute' );
    return \%retval;
}

#Does not return analysis links.
sub get_bsml_links {
    my ($elem) = @_;
    my %retval;
    foreach my $link ( $elem->children( 'Link' ) ) {
        next if( $link->att('rel') eq 'analysis' );
        my $tmp = &get_attributes( $link );
        $retval{$tmp->{'rel'}} = $tmp;
    }
    return \%retval;
}

#### Writing to Bsml Related Subroutines #################
sub to_bsml {
    #each input sequence with features gets it's own file
    my @seqs_w_features = grep { exists( $input_sequences{$_}->{'features'} ) } keys %input_sequences;
    
    foreach my $seq_id ( @seqs_w_features ) {
        my $doc = new BSML::BsmlBuilder;
        my $seq = &add_bsml_sequence( $doc, $input_sequences{ $seq_id } );
        my $ft;
        my @feature_group_ids = keys %{$input_sequences{$seq_id}->{'features'}};
        &_log("Adding features to bsml document for $seq_id\n", $DEBUG);
        my $total_fgi = scalar( @feature_group_ids );
        my $count = 0;
        foreach my $group_set ( @feature_group_ids ) {
            $ft = $doc->createAndAddFeatureTable( $seq ) unless( defined( $ft ) );
            my $feat = &add_bsml_feature_set( $doc, $seq, $ft, 
                                              $input_sequences{$seq_id}->{'features'}->{$group_set} );
            $count++;
            &_log("$count/$total_fgi\r", $DEBUG );
        }
        
         &_log("\n\n", $DEBUG );

        #write the analysis
        $doc->createAndAddAnalysis( 'id' => $analysis_name, 
                                    'sourcename' => $sourcename,
                                    'algorithm' => 'pipeline_summary',
                                    'program' => 'pipeline_summary' );

        #add the genome
        my $genome = &add_bsml_genome( $doc );

        #Add sequence link to the genome
        my $genome_id = $genome->returnattr( 'id' );
        &add_bsml_link( $doc, $seq, { 'rel' => 'genome', 'href' => "#$genome_id" } );

        my $output_file = $output_directory."/$seq_id.bsml";
        &_log("Writing document\n", $DEBUG );
        $doc->write( $output_file );
        &_log("Wrote $output_file\n", $DEBUG);
    }
}

sub add_bsml_genome {
    my ($doc) = @_;
    my $genome = $doc->createAndAddGenome( );
    my $organism = $doc->createAndAddOrganism( 'genome' => $genome, 
                                               'genus' => $organism->[0],
                                               'species' => $organism->[1] );

    $doc->createAndAddBsmlAttr( $organism, 'abbreviation', $locus_prefix );
    $doc->createAndAddBsmlAttr( $organism, 'translation_table', $trans_table );
    return $genome;
}

sub add_bsml_feature {
    my ($doc, $ft, $feat) = @_;
    my $atts = $feat->{'attributes'};
    my $id = $atts->{'id'} || &_log("Could not get id for feature", $ERROR );
    my $title = $atts->{'title'} || $id;
    my $class = $atts->{'class'};
    if( !defined( $class ) ) {
        $class = $1 if( $id =~ /^[^\.]+\.([^\.]+)\./ );
    }
    &_log("Could not parse class out for feature $id [$title]", $ERROR) 
        unless( defined( $class ) );
    
    my $bsml_feat = $doc->createAndAddFeature( $ft, $id, $title, $class );

    if( exists( $feat->{'cross-references'} ) ) {
        foreach my $cr ( @{$feat->{'cross-references'}} ) {
            &add_bsml_cross_reference( $doc, $bsml_feat, $cr );
        }
    }

    if( exists( $feat->{'interval_loc'} ) ) {
        my $il = $feat->{'interval_loc'};
        $doc->createAndAddIntervalLoc( $bsml_feat, $il->{'startpos'},
                                       $il->{'endpos'}, $il->{'complement'} );
    } else {
        &_log("Could not find interval loc in feature: $id", $ERROR );
    }

    if( exists( $feat->{'attribute_lists'} ) ) {
        foreach my $list ( @{$feat->{'attribute_lists'}} ) {
            &add_bsml_attribute_list( $doc, $bsml_feat, $list );
        }
    }

    if( exists( $feat->{'bsml_attributes'} ) ) {
        foreach my $name ( keys %{$feat->{'bsml_attributes'}} ) {
            &add_bsml_attribute( $doc, $bsml_feat, $name, $feat->{'bsml_attributes'}->{$name} );
        }
    }

    &add_bsml_analysis_link( $doc, $bsml_feat );

    return $bsml_feat;
}

sub add_bsml_attribute_list {
    my ($doc, $feat, $list) = @_;
    #for some reason, these api calls in bsml are commented out. So hack it instead of fix it
    #( I know that's right )
    
    if( !defined( $feat->{'BsmlAttributeList'} ) ) {
        $feat->{'BsmlAttributeList'} = [];
    }
    
    &_log("Can only deal with ECO value IEA currently", $ERROR ) unless( exists( $list->{'IEA'} ) );

    my ($eco, $source) = ('IEA', $list->{'IEA'} );
    delete( $list->{'IEA'} );

    my $tmp = [];
    foreach my $name ( keys %{$list} ) {
        push( @{$tmp}, { 'name' => $name, 'content' => $list->{$name} } );
    }

    push( @{$tmp}, { 'name' => 'IEA', 'content' => $list->{'IEA'} } );
    push( @{$feat->{'BsmlAttributeList'}}, $tmp );
}

sub add_bsml_cross_reference {
    my ($doc, $feat, $cr) = @_;
    return $doc->createAndAddCrossReference( 'parent' => $feat, 'database' => $cr->{'database'}, 
                                             'identifier' => $cr->{'identifier'},
                                             'identifier-type' => $cr->{'identifier-type'} );
    
}

sub add_bsml_feature_set {
    my ($doc, $seq, $ft, $feature_set) = @_;
    
    #make a feature group
    my $fg = $doc->createAndAddFeatureGroup( $seq, undef, $feature_set->{'gene'}->{'attributes'}->{'id'} );

    foreach my $class ( keys %{$feature_set} ) {
        my $id = $feature_set->{$class}->{'attributes'}->{'id'};
        my $feat = &add_bsml_feature( $doc, $ft, $feature_set->{$class} );
        $doc->createAndAddFeatureGroupMember( $fg, $id, $class );
        if( exists( $input_sequences{$id} ) ) {
            &add_bsml_sequence( $doc, $input_sequences{$id}, $class );
            &add_bsml_link( $doc, $feat, { 'rel' => 'sequence', 'href' => "#${id}_seq" } );
        }
    }
}

sub add_bsml_sequence {
    my ($doc, $seq, $class) = @_;

    if( !$seq->{'attributes'}->{'id'} ) {
        print Dumper( $seq );
    }

    my $id = $seq->{'attributes'}->{'id'} || &_log("Could not find id for sequence", $ERROR);
    my $title = $seq->{'attributes'}->{'title'} || $id;
    my $len = $seq->{'attributes'}->{'length'} || undef;
    my $mol = $seq->{'attributes'}->{'molecule'} || undef;
    $class = $seq->{'attributes'}->{'class'} unless( !exists( $seq->{'attributes'}->{'class'} ) );

    if( !defined( $class ) ) {
        #try to parse it from the id
        $class = $1 if( $id =~ /^[^\.]+\.([^\.]+)\./ );
        &_log("Could not find a class for sequence $id", $ERROR)
            unless( defined( $class ) );
    }
    
    my $bsml_seq = $doc->createAndAddSequence( $id, $title, $len, $mol, $class );

    if( exists( $seq->{'seq-data-import'} ) ) {
        &add_bsml_sdi( $doc, $bsml_seq, $seq->{'seq-data-import'} );
    }

    if( exists( $seq->{'links'} ) ) {
        foreach my $rel ( keys %{$seq->{'links'}} ) {
            &add_bsml_link( $doc, $bsml_seq, $seq->{'links'}->{$rel} );
        }
    }
    
    #add the analysis link
    &add_bsml_analysis_link( $doc, $bsml_seq );

    if( exists( $seq->{'bsml_attributes'} ) ) {
        foreach my $name ( keys %{$seq->{'bsml_attributes'}} ) {
            &add_bsml_attribute( $doc, $bsml_seq, $name, $seq->{'bsml_attributes'}->{$name} );
        }
    }

    return $bsml_seq;
}

sub add_bsml_analysis_link {
    my ($doc, $elem) = @_;
    return &add_bsml_link( $doc, $elem, { 'rel' => 'analysis', 
                                         'href' => "#$analysis_name", 
                                         'role' => 'input_of' } );
}

sub add_bsml_link {
    my ($doc, $elem, $link) = @_;
    unless( exists( $link->{'rel'} ) && exists( $link->{'href'} ) ) {
        &_log("Link information provided for element was not complete", $ERROR);
    }

    if( !exists( $link->{'role'} ) ) {
        $link->{'role'} = undef;
    }
    return $doc->createAndAddLink( $elem, $link->{'rel'}, $link->{'href'}, $link->{'role'} );
}

sub add_bsml_attribute {
    my ($doc, $elem, $name, $content) = @_;
    return $doc->createAndAddBsmlAttr( $elem, $name, $content );
}

sub add_bsml_sdi {
    my ($doc, $seq, $sdi) = @_;
    return $doc->createAndAddSeqDataImport( $seq, $sdi->{'format'},
                                            $sdi->{'source'}, undef,
                                            $sdi->{'identifier'} );
}

### other #################################
sub get_cog_desc {
    my ($id) = @_;
    &_log("ID was not passed &get_cog_desc", $ERROR) if( !defined( $id ) );
    my $retval;
    $retval = $top_cogs{ $id } if( exists( $top_cogs{ $id } ) );
    return $retval;
}

sub append_data {
    my ($id, $args) = @_;
    my $base = $input_sequences{ $id };

    if( exists( $args->{'features'} ) ) {
        foreach my $fg_id ( keys %{$args->{'features'}} ) {
            if( exists( $base->{'features'}->{$fg_id} ) ) {
                $base->{'features'}->{$fg_id} = 
                    &append_features( $base->{'features'}->{$fg_id}, $args->{'features'}->{$fg_id} );
            } else {
                $base->{'features'}->{$fg_id} = $args->{'features'}->{$fg_id};
            }
        }

    }

    if( exists( $args->{'attributes'} ) ) {
        foreach my $key ( keys %{$args->{'attributes'}} ) {
            if( exists( $base->{'attributes'}->{$key} ) ) {
                if( $base->{'attributes'}->{$key} ne $args->{'attributes'}->{$key} ) {
                    &_log("Attributes didn't match: [".$base->{'attributes'}->{$key}.
                          ", ".$args->{'attributes'}->{$key}."]", $ERROR );
                } 
            } else {
                $base->{'attributes'}->{$key} = $args->{'attributes'}->{$key};
            }
        }
    }

    if( exists( $args->{'bsml_attributes'} ) ) {
        foreach my $key ( keys %{$args->{'bsml_attributes'}} ) {
            if( exists( $base->{'bsml_attributes'}->{$key} ) ) {
                if( $base->{'bsml_attributes'}->{$key} ne $args->{'bsml_attributes'}->{$key} ) {
                    &_log("BSML Attribute values didn't match: [".$base->{'bsml_attributes'}->{$key}.
                          ", ".$args->{'bsml_attributes'}->{$key}."]", $ERROR );
                } 
            } else {
                $base->{'bsml_attributes'}->{$key} = $args->{'bsml_attributes'}->{$key}
            }
        }
    }

    if( exists( $args->{'links'} ) ) {
        foreach my $key ( keys %{$args->{'links'}} ) {
            if( !exists( $base->{'links'}->{$key} ) ) {
                $base->{'links'}->{$key} = $args->{'links'}->{$key};
            }
        }
    }

    #Seq-data-import: if there are two, use current. So do nothing.
    if( exists( $args->{'seq_data_import'} ) && 
        ( !exists( $base->{'seq_data_import'} ) || scalar( keys %{$base->{'seq_data_import'}} ) > 0 )) {
        $base->{'seq_data_import'} = $args->{'seq_data_import'};
    }
}

sub append_features {
    my ($feature_a, $feature_b) = @_;
    return $feature_a;
}

sub get_input_files {
    my $inputStr = shift;
    my @files;
    
    #comma unseparate
    my @tokens = split( /,/, $inputStr);

    foreach my $token ( @tokens ) {
        my $tfh = open_file( $token, "in" );

        my $isList = 0;

        map {
            my $line = $_;
            next if( $line =~ /^\s*$/ );
            if($line =~ /</ && !$isList) {
                push(@files, $token);
                last;
            #It's a list.
            } elsif($line =~ m|\.bsml|) {
                chomp $line;
                push(@files, $line);
                $isList = 1;
            }
        } <$tfh>;
            
    }

    return @files;
    
}

sub check_parameters {
    my $opts = shift;

    #Print perldoc if user wants help
    &_pod if($opts->{'help'});

    #debug
    $debug = 1 if( $opts->{'debug'} );

    #logger
    my $level = ( $debug ) ? 4 : 1;
    my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
    $logger = new Ergatis::Logger('LOG_FILE'=>$logfile, 'LOG_LEVEL'=> $level );
    $logger = $logger->get_logger();

    #The required options
    my @reqs = qw(input_bsml output locus_prefix organism sourcename);
    foreach my $req ( @reqs ) {
        &_log( "Option $req is required", $ERROR ) 
            unless( exists( $opts->{$req} ) );
    }

    #make an array of the input files, if there are none then die
    @input_files = &get_input_files( $opts->{'input_bsml'} );

    #output directory and locus prefix
    $output_directory = $opts->{'output'};
    $locus_prefix = $opts->{'locus_prefix'};
    $locus_database = $opts->{'locus_database'} if( $opts->{'locus_database'} );

    #make sure we can parse the species and genus from organism
    $organism = [ $1, $2 ] if( $opts->{'organism'} =~ /(\S+)\s+(.*)/ );
    &_log( "Could not parse genus and species from organism options $opts->{'organism'}")
        unless( @{$organism} == 2 );

    #other bsml lists
    @other_files = &get_input_files( $opts->{'other_bsml_lists'} )
        if( $opts->{'other_bsml_lists'} );

    #set the translation table (default: 11);
    $trans_table = $opts->{'translation_table'} if( $opts->{'translation_table'} );

    $sourcename = $opts->{'sourcename'};
    $analysis_name = $opts->{'analysis_name'} if( $opts->{'analysis_name'} );
    
    #if they included cogs
    if($opts->{'cog_search_bsml'}) {
        my $cogs = open_file( $opts->{'cog_search_bsml'}, "in" );
        chomp( my @cog_bsml_files = <$cogs> );
        close($cogs);
        
        &_log( "Option --cog_lookup is required when using --cog_search_bsml", $ERROR )
            unless($opts->{'cog_lookup'});
        my $cog_lookup = $opts->{'cog_lookup'};

        &_log("Preparsing COG Lookup file\n", $DEBUG );
        &pre_parse_cog_lookup( $cog_lookup );

        &_log("Preparsing COG results\n", $DEBUG );
        &pre_parse_cog_results( @cog_bsml_files );

        $include_cog = 1;
    }

    #if CDS Fasta was included
    if( $opts->{'cds_fasta'} ) {
        my @files = split(/,\s+/, $opts->{'cds_fasta'} );
        &_log("Processing CDS fasta file[s]\n", $DEBUG);
        map { &process_fasta_file_or_list( $_ ) } @files;
    }

    #If polypeptide fsa was included
    if( $opts->{'polypeptide_fasta'} ) {
        my @files = split(/,\s+/, $opts->{'polypeptide_fasta'} );
        &_log( "Processing polypeptide fasta file[s]\n", $DEBUG );
        map { &process_fasta_file_or_list( $_ ) } @files;
        
    }
    
}

sub pre_parse_cog_results {
    my (@bsml_files) = @_;
    
    my $total = scalar(@bsml_files);
    my $count = 0;
    &_log("$count/$total\r", $DEBUG );
    map {
        my $boh = open_file( $_, 'in' );
        my ($cog,$ref,$pval, $best);
        $best = 1;
        while( my $line = <$boh> ) {
            chomp( $line );            
            if( $line =~ /\<Seq-pair-alignment.*refseq=\"([^\"]+)\"/ ) {
                if( defined( $ref ) ) {
                    if( $pval < $best ) {
                        $best = $pval;
                        my $desc = $cogId2desc{$cog};
                        if( defined( $desc ) ) {
                            $top_cogs{$ref} = $desc;
                        }
                    }
                }
                
                $pval = 1;
                $ref = $1;
                $cog = $1 if( $line =~ /compseq=\"([^\"]+)\"/ );
                
            } elsif( $line =~ /\<Seq-pair-run.*runprob=\"([^\"]+)\"/ ) {
                if( $1 < $pval ) {
                    $pval = $1;
                }
            }
        }
        close($boh);

        #and the last one
        if( defined( $ref ) ) {
            if( $pval < $best ) {
                $best = $pval;
                my $desc = $cogId2desc{$cog};
                if( defined( $desc ) ) {
                    $top_cogs{$ref} = $desc;
                }
            }
        }
        $count++;
        &_log("$count/$total\r", $DEBUG );
        
    } @bsml_files;
}

sub pre_parse_cog_lookup {
    my $cog_lookup = shift;
    my $in = open_file( $cog_lookup, 'in' );
    my $desc;
    while( my $line = <$in>) {
        next if( $line =~ /^\s*$/ || $line =~ /^[_\s]+$/ );
		if($line =~ /COG/) {
			$desc = $line;
		} elsif( $line =~ /^\s*(\w+):\s*(.*)\s*/ ) {
            my @ids = split( /\s+/, $2 );
            map { $cogId2desc{ $_ } = $desc } @ids;
		} elsif( $line =~ /^\s+(.*)\s+/ ) {
            my @ids = split( /\s+/, $1 );
            map { $cogId2desc{ $_ } = $desc } @ids;
        }
	}	
    close($in);
}

sub process_fasta_file_or_list {
    my ($file) = @_;
    my @fsa_files;

    #Determine if it's a list or fasta
    if( $file =~ /list(.gz)?$/ ) {
        my $in = &open_file( $file, 'in' );
        chomp( @fsa_files = <$in> );
        close($in);
    } elsif( $file =~ /fsa(.gz)?$/ ) {
        push(@fsa_files, $file);
    } else {
        my $in = &open_file( $file, 'in' );
        chomp( my $test = <$in> );
        if( $test =~ /^>/ ) {
            push(@fsa_files, $file);
        } else {
            my $tmp;
            eval {
                $tmp = &open_file( $test, 'in' );
            };
            if( $@ ) {
                die("Could not determine the format of file $file. Required either ".
                      "fasta file or list of fasta files");
            }
            close($tmp);

            chomp( @fsa_files = <$in> );
            push(@fsa_files, $test );
        }
        close($in);
    }

    #index for seq data import usage
    my $total = scalar(@fsa_files);
    my $count = 0;
    foreach my $fsa_file ( @fsa_files ) {
        my $fsa = &open_file( $fsa_file, 'in' );
        
        while( <$fsa> ) {
            if( /^>(\S*)/ ) {
                $seq_data_import_info{$1} = {
                    'source' => $fsa_file,
                    'format' => 'fasta',
                    'identifier' => $1,
                }
            }
        }
        close($fsa);
        $count++;
        print "$count/$total\r";
    }
    print "\n";
    
}

sub _log {
    my ($msg, $level) = @_;
    if( $level == $DEBUG && $debug ) {
        print $msg;
        $logger->debug( $msg );
    } elsif( $level == $WARN ) {
        print STDERR $msg."\n";
        $logger->warn( $msg );
    } elsif( $level == $ERROR ) {
        $logger->fatal( $msg );
        die( $msg );
    }

}

sub _pod {   
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
