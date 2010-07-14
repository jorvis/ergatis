#!/usr/bin/perl

=head1 NAME

annotation_collection2bsml.pl - Uses output from AnnotationCollection module (using Annotation::to_string 
    subroutine) and writes information to bsml

=head1 SYNOPSIS

 USAGE: parse_evidence.pl
       --input_file=/path/to/file.txt
       --output=/path/to/output.bsml
       --source_bsml_file=/path/to/start_site_curation.all.bsml.list
     [ --sourcename=/path/to/output/files
       --log=/path/to/file.log
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>
    The input files to parse.  Can be a comma separated list of input files.

B<--input_list,-n>
    input list.  Can be a comma separated list of evidence lists. 

<--output,-o>
    Output directory.  Will name files based on parent sequence name.  Organizes all features annotated
    from input file as indicated by the feature_relationship_bsml files.  If there are features in the
    input annotation file that are not a child feature of a sequence described in a feature_relationship_bsml
    file, it will not be printed.  Check the log for the names of the sequences (but not yet, because I'm
    not logging yet)

B<--source_bsml_file,-f>
    Input file used for parent genomic sequence information as well as Seq-data information ( location of 
    sequence data files )

B<--log,-l>
    Logfile.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

    Will parse input annotation file (created by Annotation::to_string subroutine).  Uses AnnotationCollection
    module to parse the information and return base Annotation objects, which will be used to create bsml file. 
 
=head1  INPUT
    
    The input is an annotation file output from Annotation::to_string method.  See this module for specific format
    details.

=head1 OUTPUT
    
    The output

=head1  CONTACT

    Kevin Galens
    kgalens@som.umaryland.edu

=cut

use strict;
use warnings;
use PFunc::AnnotationCollection;
use File::OpenFile qw( open_file );
use XML::Twig;
use BSML::FeatureRelationshipLookup;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use BSML::BsmlBuilder;
use Data::Dumper;

my $input;
my @bsml_fr_files = ();
my $output;
my $sequences;
my %parent_seq_lookup;
my $fr_lookup;
my $printed_features;
my %seq_data_lookup;
my %feat_locs;
my $sourcename;
my $analysis_name = "p_func_analysis";
my %docs;

&check_options;

# parse the feature relationship bsml files
print "Creating the lookup\n";
$fr_lookup = new BSML::FeatureRelationshipLookup( 'bsml' => \@bsml_fr_files );

#create the feature to parent lookup from feature relationship bsml
print "Creating the parent seq lookup\n";
$sequences = &create_parent_seq_lookup( @bsml_fr_files );

# print the features
my $annot_col = new PFunc::AnnotationCollection;

$annot_col->annotations_from_file( $input );

foreach my $seq ( keys %{$sequences} ) {
    $docs{$seq}->{'doc'} = new BSML::BsmlBuilder;
    my ( $bsml_seq, $ft ) = &add_sequence_to_doc( $docs{$seq}->{'doc'}, $seq );
    $docs{$seq}->{'bsml_seq'} = $bsml_seq;
    $docs{$seq}->{'feature_table'} = $ft;
}

while( my $annot = $annot_col->next_annotation_from_file ) {
    #find the parent seq and make sure it's added to the document
    die("Could not find the parent sequence for ".$annot->get_feature_id) 
        unless( exists( $parent_seq_lookup{$annot->get_feature_id} ) );
    my $parent_seq = $parent_seq_lookup{$annot->get_feature_id};

    if( !exists( $docs{$parent_seq} ) ) {
        die("Could not find document for $parent_seq");
    }
    
    &add_gene_to_feature_table( $docs{$parent_seq}->{'doc'}, $docs{$parent_seq}->{'bsml_seq'}, 
                                $docs{$parent_seq}->{'feature_table'}, $annot );

}

#Add analyses and write
foreach my $pseq ( keys %docs ) {
    &add_analysis_to_doc( $docs{ $pseq }->{'doc'}, $analysis_name );
    my $outfile = $output."/$pseq.bsml";
    $docs{$pseq}->{'doc'}->write( $outfile );
    print "Wrote $outfile\n";
}

sub add_analysis_to_doc {
    my ($doc, $a_name) = @_;
    $doc->createAndAddAnalysis( id => $a_name,
                                sourcename => $sourcename,
                                algorithm => "p_func",
                                program => "p_func");
}

sub add_gene_to_feature_table {
    my ($doc, $seq, $ft, $annot) = @_;

    my $loc = $feat_locs{ $annot->get_feature_id };

    my $fg_info = [];
    my $gene_id;

    foreach my $class ( qw( gene polypeptide transcript CDS exon ) ) {
        my $feature_id = $fr_lookup->lookup( $annot->get_feature_id, $class );
        
        #add the feature
        die("Feature table was not defined") unless(defined($ft));
        my $feat =  $doc->createAndAddFeature( $ft, $feature_id, 
                                               $feature_id, $class );
        #add the interval loc
        my $intloc = $doc->createAndAddIntervalLoc( $feat, $loc->[0], $loc->[1], $loc->[2] );

        #add the link
        die("Feature $feature_id is undefined") unless( defined($feat) );
        my $link = $doc->createAndAddLink( $feat, 'analysis', "#$analysis_name", "input_of");
        
        #indirect way of indicating where our annotation is
        if( $feature_id eq $annot->get_feature_id ) {
            &add_annotation_to_feature( $doc, $feat, $annot );
        }

        if( $class eq 'gene' ) {
            $gene_id = $feature_id;
        }
        push(@$fg_info, [$feature_id, $class]);
    }

    my $fg = $doc->createAndAddFeatureGroup( $seq, '', $gene_id );
    foreach my $fgm ( @{$fg_info} ) {
        $fg->addBsmlFeatureGroupMember( $fgm->[0], $fgm->[1] );
    }
}

sub add_annotation_to_feature {
    my ($doc, $feat, $annot) = @_;
    
    #gene_product_name
    if( $annot->has_annotation( 'gene_product_name' ) ) {
        my ($name, $source) = $annot->get_gene_product_name;
        $doc->createAndAddBsmlAttribute( $feat, 'gene_product_name', $name->[0] );
        $doc->createAndAddBsmlAttribute( $feat, 'gene_product_name_source', $source );
    }

    #gene_symbol
    if( $annot->has_annotation( 'gene_symbol' ) ) {
        my ($sym, $source) = $annot->get_gene_symbol;
        $doc->createAndAddBsmlAttribute( $feat, 'gene_symbol', $sym->[0] );
        $doc->createAndAddBsmlAttribute( $feat, 'gene_symbol_source', $source );
    }

    #EC
    if( $annot->has_annotation( 'EC' ) ) {
        my ($ec, $source) = $annot->get_EC;
        foreach my $entry ( @{$ec} ) {
            my $attList = [ { 'name' => "EC", 'content' => 
                                  $entry },
                            { 'name' => "IEA", 'content' => 
                                  $source } ];
            
            push( @{$feat->{'BsmlAttributeList'}}, $attList );
        }
    }

    #GO
    if( $annot->has_annotation( 'GO' ) ) {
        my ($gos, $source) = $annot->get_GO;
        foreach my $go ( @{$gos} ) {

            my $attList = [ { 'name' => "GO", 'content' => 
                                  $go },
                            { 'name' => "IEA", 'content' => 
                                  $source } ];
            
            push( @{$feat->{'BsmlAttributeList'}}, $attList );
        }
    }

    #TIGR_Role
    if( $annot->has_annotation( 'TIGR_Role' ) ) {
        my ($roles, $source) = $annot->get_TIGR_Role;
        foreach my $role ( @{$roles} ) {
            my $attList = [ { 'name' => "TIGR_role", 'content' => 
                                  $role },
                            { 'name' => "IEA", 'content' => 
                                  $source } ];
            
            push( @{$feat->{'BsmlAttributeList'}}, $attList );
        }
    }
}

sub add_sequence_to_doc {
    my ($doc, $parent_seq) = @_;

    my $seq;
    $seq = $doc->returnBsmlSequenceByIDR( $parent_seq );
    if( !defined( $seq ) ) {
        $seq = $doc->createAndAddSequence( $parent_seq, $sequences->{$parent_seq}->{'title'},
                                           $sequences->{$parent_seq}->{'length'},
                                           $sequences->{$parent_seq}->{'molecule'},
                                           $sequences->{$parent_seq}->{'class'} );
    }

    #create a feature table
    my $ft;
    my $ft_list = $seq->returnBsmlFeatureTableListR;
    if( !defined($ft_list) || @{$ft_list} == 0 ) {
        $ft = $doc->createAndAddFeatureTable( $seq );
        die("But it wasn't defined! wtf") if( !defined( $ft ) );
    } else {
        $ft = $ft_list->[0];
    }

    #add the seq data import
    my $sdi;
    $sdi = $seq->returnBsmlSeqDataImport();
    if( !defined( $sdi ) || scalar(keys %$sdi) == 0 ) {
        $sdi = $doc->createAndAddSeqDataImport( $seq, $seq_data_lookup{$parent_seq}->{'format'},
                                                $seq_data_lookup{$parent_seq}->{'source'}, '',
                                                $seq_data_lookup{$parent_seq}->{'identifier'} );
    }

    return ($seq, $ft);
}


sub create_parent_seq_lookup {
    my @bsmls = @_;

    my %retval;
    
    my $twig = new XML::Twig( 'twig_roots' => {
        'Sequence' => sub {
            my ($t, $el) = @_;
            my $seq_id = $el->att('id');
            $retval{$seq_id} = {};
            map { push(@{$retval{$seq_id}->{'features'}}, $_->att('id') );
                  $parent_seq_lookup{ $_->att('id') } = $seq_id;
              } $el->find_nodes('Feature-tables/Feature-table/Feature');
            map { $feat_locs{ $_->att('id') } = &get_feature_locs( $_ ) } $el->find_nodes('Feature-tables/Feature-table/Feature[@class="transcript"]');

            my ($sdi) = $el->find_nodes( 'Seq-data-import' );
            die("Could nto find sdi for $seq_id") unless( defined( $sdi ) );
            my ($source, $identifier, $format) = ($sdi->att('source'), $sdi->att('identifier'),
                                                  $sdi->att('format') );

            $seq_data_lookup{$seq_id} = {
                'source'     => $source,
                'identifier' => $identifier,
                'format'     => $format,
            };

            $retval{$seq_id}->{'length'} = $el->att('length') if( $el->att('length') );
            $retval{$seq_id}->{'class'} = $el->att('class') if( $el->att('class') );
            $retval{$seq_id}->{'molecule'} = $el->att('molecule') if( $el->att('molecule') );
            $retval{$seq_id}->{'title'} = $el->att('title') if( $el->att('title') );

        },
    } );

    foreach my $bsml ( @bsmls ) {       
        my $in = open_file( $bsml, 'in' );
        $twig->parse( $in );
        close($in);
    }

    return \%retval;
}

sub get_feature_locs {
    my ($feat) = @_;
    my $intloc = $feat->first_child('Interval-loc');
    my @coords = ($intloc->att('startpos'), $intloc->att('endpos'), 
                  $intloc->att('complement') );
    return \@coords;
}

sub check_options {
    my %options;
    my $results = GetOptions (\%options,
                              'input_file|i=s',
                              'output|o=s',
                              'source_bsml_file|f=s',
                              'sourcename|s=s',
                              'log|l=s',
                              'help|h',
                              );
    
    if( $options{'help'} ) {
        &_pod;
    }

    if( $options{'input_file'} ) {
        $input = $options{'input_file'};
    } else {
        die("Option --input_file is required");
    }

    if( $options{'source_bsml_file'} ) {
        my $in = open_file( $options{'source_bsml_file'}, 'in' );
        chomp( my @tmp = <$in> );
        close($in);
        push(@bsml_fr_files, @tmp);
    }

    unless( @bsml_fr_files  > 0 ) {
        die("Option --source_bsml_file is required");
    }

    if( $options{'output'} ) { 
        die("Value for option --output was not a directory [$options{'output'}]") 
            if( ! -d $options{'output'} );
        $output = $options{'output'};
    }

    if( $options{'sourcename'} ) {
        $sourcename = $options{'sourcename'};
    } else {
        use Cwd;
        $sourcename = getcwd;
    }

}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2} );
}
