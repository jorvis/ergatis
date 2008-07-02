#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";


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
    kgalens@tigr.org

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Ergatis::Logger;
use Ergatis::IdGenerator;
use BSML::BsmlBuilder;
use XML::Twig;
use PerlIO::gzip;
use File::OpenFile qw(open_file);
use Data::Dumper;

####### GLOBALS AND CONSTANTS ###########
my @inputFiles;               #Holds input (gene prediction bsml) files
my @otherFiles;               #Holds other bsml files.
my $output;                   #Output directory

my $outputFile;               #Output file (includes directory).
my $debug;                    #The debug variable
my $c;
my $analysisId = 'pipeline_summary';
my $sourcename;
my $features;
my $featureGroups;
my $topCogHits;
my @cogBsmls;
my $cogLookup;
my $locusPrefix;
my $prefixCount = 1;
my $featElemsBsml = {};
my $organism = [];
my $transTable;
my $attListId = 0;
my $crid = 0;               #Cross reference id counter
my %seq_data_import_info;
########################################
$|++;

my %options = ();
my $results = GetOptions (\%options, 
                          'input_bsml|i=s',
                          'other_bsml_lists|b=s',
                          'output|o=s',
                          'correct_mistakes|c',
                          'locus_prefix|u=s',
                          'organism|r=s',
                          'translation_table|t=s',
                          'cog_search_bsml|s=s',
                          'cog_lookup|g=s',
                          'cds_fasta|f=s',
                          'polypeptide_fasta|p=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || &_pod;

#Setup the logger
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

my $genomeId;

# Check the options.
&check_parameters(\%options);

foreach my $file (@inputFiles) {
    print STDOUT "using $file\n";

    $sourcename = $1 if($file =~ m|(.*\d+_[^/]+)/|);
    &_die("Could not parse the sourcename directory for analysis section from $file") 
        unless($sourcename);

    #Make the output file name.
    my $base = $1 if($file =~ m|.*/(([^/\.]+\.){3})|);
    $base =~ s/\.$//;
    my $fileName = $1 if($file =~ m|.*/([^/]+)$|);
    if( -d $output ) {
        $outputFile = "$output/$fileName";
    } else {
        $logger->logdie("Option --output must be a directory");
    }

    my $bsml = new BSML::BsmlBuilder;

    #Add genome/organism
    my $genome = $bsml->createAndAddGenome( );
    $genomeId = $genome->returnattr( 'id' );
    my $organism = $bsml->createAndAddOrganism( 'genome' => $genome, 
                                                'genus' => $organism->[0],
                                                'species' => $organism->[1] );
    $logger->logdie("Couldn't get the genome id") unless($genomeId);

    $bsml->createAndAddBsmlAttr( $organism, 'abbreviation', $locusPrefix );
    $bsml->createAndAddBsmlAttr( $organism, 'translation_table', $transTable );

    #Parse the bsml file
    my $twig = new XML::Twig( twig_handlers => {
        'Sequence' => sub { &printSequenceFeatures( $bsml, @_ ) },
        });

    my $openFile;
    if( -e $file ) {
        open($openFile, "< $file") or 
            $logger->logdie("Unable to open $file ($!)");
    } elsif( -e $file.".gz" ) {
        open($openFile, "<:gzip", $file.".gz") or
            $logger->logdie("Can't open $file.gz ($!)");
    } else {
        $logger->logdie("Couldn't find $file");
    }
    

    $twig->parse($openFile);
    print STDOUT "Done parsing the input file $file\n";

    my @addThese;
    foreach my $otherFile (@otherFiles) {
        push(@addThese, $otherFile) if($otherFile =~ /$base\./);
    }

    foreach my $addThis ( @addThese ) {
        print STDOUT "adding $addThis\n";
        my $newTwig = new XML::Twig( twig_handlers => {
            'Sequence' => sub { &printSequenceFeatures( $bsml, @_ ) }
            });

        my $otherOpenFile;

        if( -e $addThis ) {
            open($otherOpenFile, "< $addThis") or 
                $logger->logdie("Unable to open $addThis ($!)");
        } elsif( -e $addThis.".gz" ) {
            open($otherOpenFile, "<:gzip", $addThis.".gz") or
                $logger->logdie("Can't open $addThis.gz ($!)");
        } else {
            $logger->logdie("Couldn't find $addThis");
        }

        $newTwig->parse($otherOpenFile);
    }

    #Add missing sequence elements for features
    &add_feature_sequence_elements( $bsml );

    print STDOUT "Writing to $outputFile\n";
    $bsml->write($outputFile);

    print STDOUT "Fixing attributes\n";
    my $tmpOutputFile = $outputFile.".part";
    my $toh;
    open($toh, "> $tmpOutputFile") or $logger->logdie("Unable to open tmp output $tmpOutputFile");

    #Now we have to go back and fix the attribute lists
    my $attTwig = new XML::Twig( 'twig_roots' => { 'Feature' => \&fixAttributes },
                                 'twig_print_outside_roots' => $toh,
                                 'pretty_print' => 'indented' );

    $attTwig->parsefile( $outputFile );

    close($toh);

    print STDOUT "mv $tmpOutputFile $outputFile\n";
    system("mv $tmpOutputFile $outputFile")

}


######################## SUB ROUTINES #######################################
sub add_feature_sequence_elements {
    my ($bsml) = @_;
    
    #Grab all the Sequence elements to grab all the Feature elements
    my $sequences = $bsml->returnBsmlSequenceListR;

    foreach my $seq ( @{$sequences} ) {
        my $feature_tables = $seq->returnBsmlFeatureTableListR;
        foreach my $ft ( @{$feature_tables} ) {
            my $features = $ft->returnBsmlFeatureListR;
            foreach my $feat ( @{$features} ) {
                my $feat_id = $feat->returnattr('id');
                my $class = $feat->returnattr('class');
                $logger->logdie("Could not parse class from feature $feat_id") unless( $class );

                my $molecule;
                
                if( $class eq 'CDS' ) {
                    $molecule = 'dna';
                } elsif( $class eq 'polypeptide' ) {
                    $molecule = 'aa';
                } else {
                    next;
                }
                
                #do we have a sequence element for it already?
                my $feature_sequence = $bsml->returnBsmlSequenceByIDR( $feat_id );
                undef $feature_sequence if( ref($feature_sequence) ne 'BSML::BsmLSequence' );
                $feature_sequence = $bsml->returnBsmlSequenceByIDR( $feat_id."_seq" )
                    unless( $feature_sequence );


                my $seq_data_import_info = &lookup_seq_data_import_info( $feat_id );

                #if we don't have the info, this probably means it's a tRNA CDS
                next unless( $seq_data_import_info );

                my $add_sdi_flag = 0;
                if( $feature_sequence ) {
                    #make sure the seq data import is up to date
                    my $seq_data_import = $feature_sequence->returnBsmlSeqDataImport;

                    unless( $seq_data_import->returnattr( 'source' ) eq $seq_data_import_info->{'source'} &&
                            $seq_data_import->returnattr( 'identifier' ) eq $seq_data_import_info->{'identifier'} ) {
                        $feature_sequence->dropBsmlSeqDataImport;

                        $add_sdi_flag = 1;

                    }
                    

                } else {
                    #add a new sequence element
                    #  my ( $id, $title, $length, $molecule, $class ) = @_;
                    $feature_sequence = $bsml->createAndAddSequence( $feat_id."_seq", $feat_id,
                                                                         '', $molecule, $class );

                    $add_sdi_flag = 1;

                    #my ($elem, $rel, $href, $role) = @_;
                    my $link = $bsml->createAndAddLink( $feature_sequence, 'analysis', "#".$analysisId."_analysis",
                                                        'input_of' );

                }
                
                if( $add_sdi_flag ) {
                    #my ( $seq, $format, $source, $id, $identifier ) = @_;
                    my $sdi = $bsml->createAndAddSeqDataImport( $feature_sequence, $seq_data_import_info->{'format'},
                                                                $seq_data_import_info->{'source'}, '', 
                                                                $seq_data_import_info->{'identifier'} );
                }
            } #Feature
        } #Feature-table
    } #Sequence
}
sub fixAttributes {
    my ($twig, $featElem) = @_;
    my $featId = $featElem->att('id');
    $logger->logdie("Couldn't parse feature id from feature") unless($featId);
    
    my @elemOrder = ( 'Attribute', 'Interval-loc', 'Attribute-list', 'Cross-reference', 'Link' );
    my %finalElts;

    #Grab each child and cut it out and save it for later.
    my %elts;
    
    foreach my $child ( $featElem->children() ) {
        my $tag = $child->gi();
        push( @{$elts{$tag}}, $child );
        $child->cut;
    }


    #where did the annotation come from?
    my $annotationFrom;
    foreach my $att ( @{$elts{'Attribute'}} ) {
        $annotationFrom = $att->att('content') if( $att->att('name') eq 'gene_product_name_source' );
    }

    foreach my $att ( @{$elts{'Attribute'}} ) {
        
        if($att->att('name') eq 'ec_number') {
            #Where did we get the ec number from?
            my $ecFrom = $att->att('ec_number_source') if($att->att('ec_number_source'));
            $ecFrom = $annotationFrom if( (!$ecFrom) && $annotationFrom );
            $logger->logdie("We should not have an ec_number if we don't know where the annotation came from") unless($ecFrom);
            
            my $attList = new XML::Twig::Elt( 'Attribute-list', { 'id' => "attList_".$attListId++ } );

            my $firstAtt = new XML::Twig::Elt( 'Attribute', { 'name' => 'EC',
                                                              'content' => $att->att('content') } );
            my $secAtt = new XML::Twig::Elt( 'Attribute', { 'name' => 'IEA',
                                                            'content' => $ecFrom } );

            $firstAtt->paste( $attList );
            $secAtt->paste( 'after' => $firstAtt );
          
            push(@{$elts{'Attribute-list'}}, $attList);
            
        } elsif($att->att('name') eq 'role_id') {
            #Where did we get the TIGR role from?
            my $roleFrom = $att->att('role_from') if($att->att('role_id_from'));
            $roleFrom = $annotationFrom if( ( (!$roleFrom) || ($roleFrom eq 'guess_role') ) && $annotationFrom );
            $logger->logdie("We should not have an ec_number if we don't know where the annotation came from") unless($roleFrom);

            my $attList = new XML::Twig::Elt( 'Attribute-list', { 'id' => "attList_".$attListId } );
            $attListId++;

            my $firstAtt = new XML::Twig::Elt( 'Attribute', { 'name' => 'TIGR_role',
                                                              'content' => $att->att('content') } );
            my $secAtt = new XML::Twig::Elt( 'Attribute', { 'name' => 'IEA',
                                                            'content' => $roleFrom } );

            $firstAtt->paste( $attList );
            $secAtt->paste( 'after' => $firstAtt );
       
            push(@{$elts{'Attribute-list'}}, $attList);

        } elsif($att->att('name') eq 'go_id') {
            #Where did we get the GO id from?
            my $goFrom = $att->att('go_id_from') if($att->att('go_id_from'));
            $goFrom = $annotationFrom if( (!$goFrom) && $annotationFrom );
            $logger->logdie("We should not have a go_id if we don't know where the annotation came from") unless($goFrom);

            my $attList = new XML::Twig::Elt( 'Attribute-list', { 'id' => "attList_".$attListId } );
            $attListId++;

            my $firstAtt = new XML::Twig::Elt( 'Attribute', { 'name' => 'GO',
                                                              'content' => $att->att('content') } );
            my $secAtt = new XML::Twig::Elt( 'Attribute', { 'name' => 'IEA',
                                                            'content' => $goFrom } );

            $firstAtt->paste( $attList );
            $secAtt->paste( 'after' => $firstAtt );
           
            push(@{$elts{'Attribute-list'}}, $attList);
            

        } elsif( $att->att('name') eq 'locus_tag' ) {
            my $crossRef = new XML::Twig::Elt( 'Cross-reference', { 
                'database' => 'TIGR_moore',
                'identifier' => $att->att('content'),
                'id' => 'CrossReference_'.$crid++,
                'identifier-type' => 'locus' } );
            
            push(@{$elts{'Cross-reference'}}, $crossRef);

        } elsif( $att->att('name') eq 'role_from' ) {
            
        } elsif( $att->att('name') eq 'top_cog_hit' ) {
            my $str = $att->att('content');
            $str =~ s/\s+$//g;
            $att->set_att('content', $str);
            push(@{$finalElts{'Attribute'}}, $att);
        } else {
            push(@{$finalElts{'Attribute'}}, $att);
        }
    }
    
    foreach my $key ( keys %elts ) {
        next if($key eq 'Attribute');
        $finalElts{$key} = $elts{$key};
    }
 
    if($featId =~ /(polypeptide|CDS)/) {
            
        #Add a link element
        my $newLink = new XML::Twig::Elt( 'Link', {
            'rel' => 'sequence',
            'href' => "#${featId}_seq" } );

        push( @{$finalElts{'Link'}}, $newLink);
        
    }

    #Now paste them in the correct order
    foreach my $tagName ( reverse @elemOrder ) {
        foreach my $featChild( @{$finalElts{$tagName}} ) {
            $featChild->paste( $featElem );
        }
    }

    $featElem->print;
    
}

sub getTopCogHit {
	my $bsmlFiles = shift;
    my $polyId = shift;
 	my $cogId;
	my $maxPval = 0;
	
    my $parseThisFile = "";
	foreach my $bsmlFile ( @{$bsmlFiles} ) {
        #Find the file the contains this polypeptide and parse it
        if($bsmlFile =~ /$polyId\./) {
            $parseThisFile = $bsmlFile;
            last;
        }
    }

    $logger->warn("couldn't find cog bsml file for $polyId") unless($parseThisFile);
    return "" unless($parseThisFile);

    my $cogTwig = new XML::Twig( twigHandlers => {
		'Seq-pair-alignment' => sub { &cogSpaHandler( @_, \$cogId, \$maxPval ) } }
                                 );

    my $fh;

    #Check for gzip extension
    $parseThisFile =~ s/\.gz$//;
    if( -e $parseThisFile ) {
        open($fh, "< $parseThisFile") or 
            $logger->logdie("Unable to open $parseThisFile ($!)");
    } elsif( -e $parseThisFile.".gz" ) {
        open($fh, "<:gzip", $parseThisFile.".gz") or
            $logger->logdie("Can't open $parseThisFile.gz ($!)");
    }
    
    $cogTwig->parse($fh);	
    close($fh);

    return $cogId;
}

sub cogSpaHandler {
    my ($twig, $spaElem, $retval, $maxPval) = @_;

    #Get the seq-pair-run children
    my @sprs = $spaElem->children('Seq-pair-run');
    $logger->logdie("Couldn't find any Seq-pair-run elements in Seq-pair-aligment")
        unless( @sprs > 0 );

    my $queryId = $spaElem->att('refseq');
    my $matchId = $spaElem->att('compseq');
    my $cogId = &getCog($matchId);

    foreach my $spr ( @sprs ) {
        my $curPval = 'noPval';

        #Get the attribute elements
        my @attributes = $spr->children('Attribute');
        $logger->logdie("Unable to find Attribute elements in Seq-pair-run")
            unless( @attributes > 0 );

        foreach my $attribute ( @attributes ) {
            next unless($attribute->att('name') eq 'p_value');
            $curPval = $attribute->att('content');
        }
        
        $logger->logdie("Seq-pair-run didn't have p_value Attribute") 
            if( $curPval eq 'noPval' );  

        #Check the pvalue to see if it's higher
        if($$maxPval < $curPval) {
            $$retval = $cogId;
            $$maxPval = $curPval;
        }
    }

}

sub getCog {
	my $id = shift;
	my $retval = 'nocog';

    my $in;
    $cogLookup =~ s/\.gz$//;
	if( -e $cogLookup ) {
        open($in, "< $cogLookup") or
            $logger->logdie("Unable to open cog lookup file $cogLookup");
    } elsif( -e $cogLookup.".gz" ) {
        open($in, "<:gzip", $cogLookup.".gz") or
            $logger->logdie("Unable to open $cogLookup.gz");
    }

	while( my $line = <$in>) {
		if($line =~ /COG/) {
			$retval = $line;
		} elsif($line =~ /$id/) {
			last;
		}
	}	
    close($in);
    
	return $retval;
}

sub printSequenceFeatures {
    my ($bsml, $twig, $seqElem) = @_;
    
    my $id = $seqElem->att('id');
    &_die("Can't parse id from Sequence in bsml file") unless($id);

    my $bsmlSequence = &addSequence($bsml, $seqElem);

    &_die("bsmlSequence ".$seqElem->att('id')." does not have a class") unless($seqElem->att('class'));

    return unless($seqElem->att('class') eq 'assembly');

    my $asmblId = $id;

    my $fTables = $seqElem->first_child('Feature-tables');
    my $fT = $fTables->first_child('Feature-table') if($fTables);
    my @features = $fT->children('Feature') if($fT);

    my $newFeatureTable;
    #Create new feature table
    if(@{$bsmlSequence->{'BsmlFeatureTables'}}) {
        $newFeatureTable = $bsmlSequence->{'BsmlFeatureTables'}->[0];
    } else {
        $newFeatureTable = $bsml->createAndAddFeatureTable( $bsmlSequence );
    }
    &_die("Couldn't make the new feature table for $asmblId") unless($newFeatureTable);

    print STDOUT "There are a total of ".scalar(@features)."\n";
    my $count = 1;
    foreach my $featElem ( @features ) {
        print STDOUT "\r$count";
        $count++;
        my $bsmlFeature = &addFeature( $bsml, $newFeatureTable, $featElem );
        $featElemsBsml->{$featElem->att('id')} = $bsmlFeature;
        &_die("Unable to create new bsml feature for sequence $asmblId") unless($bsmlFeature);
    }
    print STDOUT "\n";

    my @fgs = $fTables->children('Feature-group') if($fTables);
    
    print STDOUT "Starting feature groups\n";
    print STDOUT "There are a total of ".scalar(@fgs)."\n";
    $count = 1;
    foreach my $fg ( @fgs ) {
        print STDOUT "\r$count";
        $count++;
        my $bsmlFg = 1;#&addFeatureGroup( $bsml, $bsmlSequence, $fg );
        &_die("Could not create new Feature-group") unless($bsmlFg);
    }
    print STDOUT "\nFinished feature groups\n";

    my $analysis;
    if($bsml->returnBsmlAnalysisR(0)) {
        $analysis = $bsml->returnBsmlAnalysisR(0);
    } else {
        $analysis = $bsml->createAndAddAnalysis( 'id' => $analysisId."_analysis",
                                                 'program' => 'pipeline_summary',
                                                 'algorithm' => 'pipeline_summary', 
                                                 'sourcename' => $sourcename );
    }
    print STDOUT "Added analysis\n";
    &_die("Could not find or make an analysis") unless($analysis);
    
    
}

sub addFeatureGroup {
    my ($bsml, $bsmlSequence, $fgElem) = @_;

    my $fgId = $fgElem->att('id');
    &_die("Can't parse id from Feature-group") unless($fgId);

    
    my $groupSet = $fgElem->att('group-set');
    &_die("Could not parse group-set attribute from Feature-group $fgId") unless($groupSet);

    return $featureGroups->{$groupSet} if(exists($featureGroups->{$groupSet}));

    my @fgms = $fgElem->children('Feature-group-member');
    &_die("There were no Feature-group-members in Feature-group id ($fgId), group-set ($groupSet)")
        unless(@fgms > 0);


    my ($geneId, $topCogHit, $cdsId);
    foreach my $fgm ( @fgms ) {
        if ($fgm->att('feature-type') eq 'gene') {
            $geneId = $fgm->att('featref');
            &_die("Can't parse featref from Feature-group-member in $groupSet (Feature-group)") unless($geneId);
        } elsif($fgm->att('feature-type') eq 'polypeptide') {
            $topCogHit = &getTopCogHit(\@cogBsmls, $fgm->att('featref'));
        } elsif($fgm->att('feature-type') eq 'CDS') {
            $cdsId = $fgm->att('featref');
        }
    }

    if($topCogHit) {

        my $featureElem = $featElemsBsml->{$cdsId};
        $logger->logdie("Wasn't able to get the feature from featElemsBsml for $cdsId\n") 
            unless($featureElem);
        $bsml->createAndAddBsmlAttribute( $featElemsBsml->{$cdsId}, 'top_cog_hit', $topCogHit);
    }

    &_die("There was no gene Feature-group-member in the Feature-group $groupSet") unless($geneId);
    
    $groupSet = $geneId;

    my $newFeatureGroup = $bsml->createAndAddFeatureGroup( $bsmlSequence, '', $geneId );

    foreach my $fgm (@fgms) {
        
        my $featref = $fgm->att('featref');
        &_die("Can't parse featref from Feature-group-member in Feature-group $groupSet") unless($featref);

        my $feattype = $fgm->att('feature-type');
        unless($feattype) {
            &_die("Can't parse feattype from Feature-group-member $featref in Feature-group $groupSet")
                unless($c);
            $feattype = $1 if($featref =~ /^[^\.]+\.([^\.]+)\./);
            &_die("Could not parse feattype from $featref") unless($feattype);
        }
        &_die("Could not get a feattype for Feature-group-member ($featref)") unless($feattype);

        $bsml->createAndAddFeatureGroupMember( $newFeatureGroup, $featref, $feattype );
    }

    $featureGroups->{$groupSet} = $newFeatureGroup;
    return $newFeatureGroup;
    
    
}

sub addFeature { 
    my ($bsml, $fTable, $feature) = @_;

    my $featId = $feature->att('id');
    &_die("Can't parse feature id from a feature") unless($featId);

    my $class = $feature->att('class');
    
    unless($class) {
        &_die("Can't get Feature class from feature $featId located") unless($c);
        $class = $1 if($featId =~ /^[^\.]+\.([^\.]+)\./);
    }
    &_die("Can't parse class from $featId") unless($class);

    my $title = $feature->att('title');
    unless($title) {
        &_die("Unable to parse title from feature $featId") unless($c);
        $title = $featId;
    }

    my $bsmlFeature;
    if(exists($features->{$featId})) {
        $bsmlFeature = $features->{$featId};
    } else {
        $bsmlFeature = $bsml->createAndAddFeature( $fTable, $featId, $title, $class );
    }

    #We only want to add this once.
    my $intLoc;

    unless($bsmlFeature && $bsmlFeature->{'BsmlInterval-loc'} && @{$bsmlFeature->{'BsmlInterval-loc'}}) {

        #Now create the interval loc.
        $intLoc = $feature->first_child('Interval-loc');

        if($intLoc) {

            #addBsmlIntervalLoc( $startpos, $endpos, $complement, $startopen, $endopen )
            my ($startpos, $endpos, $complement) = ($intLoc->att('startpos'), $intLoc->att('endpos'), $intLoc->att('complement'));
            
            &_die("Feature $featId does not contain a complete Interval-loc element")
                unless(defined($startpos) && defined($endpos) && defined($complement)); 
            
            if($startpos > $endpos) {
                $logger->logdie("Startpos ($startpos) is greater than endpos ($endpos).  Feature $featId");
                ($startpos, $endpos) = ($endpos, $startpos);
            }
            
            $bsmlFeature->addBsmlIntervalLoc( $startpos, $endpos, $complement );
        }

    } else {
        $intLoc = $bsmlFeature->{'BsmlInterval-loc'}->[0];
    }

    my $link;

    if( $bsmlFeature->{'BsmlLink'} && $bsmlFeature->{'BsmlLink'}->[0] ) {
        $link = $bsmlFeature->{'BsmlLink'}->[0];
    } else {
        $bsml->createAndAddLink( $bsmlFeature, 'analysis', "#${analysisId}_analysis", "input_of" );
        $link = $bsmlFeature->{'BsmlLink'}->[0];
    }
    &_die("Couldn't add or find link element associated with Feature $featId") unless($link);


    #Add attributes if any
    my @attributes = $feature->children('Attribute');
    &addAttributes( $bsml, $bsmlFeature, \@attributes ) unless( @attributes < 1 );

    unless($class ne 'gene' || exists($bsmlFeature->{'BsmlAttr'}->{'locus_tag'})) {
        $bsmlFeature->addBsmlAttr( 'locus_tag', $locusPrefix."_".$prefixCount++);
    }

    #Add attribute lists if any
    my @att_lists = $feature->children( 'Attribute-list' );
    &addAttributeLists( $bsml, $bsmlFeature, \@att_lists ) unless( @att_lists < 1 );
    


    $features->{$featId} = $bsmlFeature;
    return $bsmlFeature;

}

sub addAttributeLists {
    my ($bsml, $bsmlFeature, $att_lists) = @_;

    foreach my $att_list ( @{$att_lists} ) {
        my ($name, $content);
        my $bsml_att_list = [];
        
        foreach my $attribute ( $att_list->children('Attribute') ) {
            ($name, $content) = ($attribute->att('name'), $attribute->att('content') );
            my $att = { 'name' => $name, 'content' => $content };
            push( @{$bsml_att_list}, $att );
        }

        push( @{$bsmlFeature->{'BsmlAttributeList'}}, $bsml_att_list );

    } 
}

sub addAttributes {
    my ($bsml, $bsmlFeature, $atts) = @_;

    foreach my $att ( @{$atts} ) {
        my ($name, $value) = ($att->att('name'), $att->att('content'));
        $name =~ s/gene_sym$/gene_symbol/;
        $name =~ s/gene_symbolbol/gene_symbol/;
        $name =~ s/ec_num$/ec_number/;
        $name =~ s/annotation_from/gene_product_name_source/;
        $name =~ s/com_name/gene_product_name/;
        my $bsmlAtt = $bsml->createAndAddBsmlAttribute( $bsmlFeature, $name, $value);
    }

    my $hash = $bsmlFeature->returnattrHashR();
    my $tmpClass = $hash->{'class'};

    unless(exists($bsmlFeature->{'BsmlAttr'}->{'gene_product_name'}) || $tmpClass eq 'tRNA' ) {
        my $bsmlAtt = $bsml->createAndAddBsmlAttribute( $bsmlFeature, 'gene_product_name', 'hypothetical protein');
    }
    
}

sub addSequence {
    my ($bsml, $seq) = @_;

    my $seqId = $seq->att('id');
    &_die("Can't parse sequence id from Sequence element") unless($seqId);

    my $newSeq = $bsml->returnBsmlSequenceByIDR( $seqId );
    return $newSeq if($newSeq);

    my $title = $seq->att('title');
    unless($title) {
        &_die("Can't parse sequence title from Sequence element $seqId") unless($c);
        $title = $seqId if($c);
    }

    my $length;
    $length = $seq->att('length');

    my $sdi = $seq->first_child('Seq-data-import');
    &_die("Sequence $seqId does not have Seq-data-import element") unless($sdi);

    my $fsa = $sdi->att('source');
    &_die("Seq-data-import (for sequence $seqId) does not have a source attribute") unless($fsa);

    my $ident = $sdi->att('identifier');
    &_die("Seq-data-import (for sequence $seqId) does not have an identifier attribute") unless($ident);
    
    unless($length) {
        $length = &determineLength($fsa, $ident);
    }
    
    my $class;
    $class = $seq->att('class');
    unless($class) {
        &_die("Can't parse attribute class from Sequence $seqId") unless($c);
        $class = $1 if($seqId =~ /^[^\.]+\.([^\.]+)\./);
        $seq->set_att('class', $class);
    }
    &_die("Sequence $seqId does not have class attribute") unless($class);

    my $molecule = $seq->att('molecule');
    unless($molecule) {
        $molecule = 'dna';
        $molecule = 'aa' if($class eq 'polypeptide');
    }
    &_die("Could not assign either aa or dna for molecule attribute") unless($molecule);

    my $bsmlSequence = $bsml->createAndAddSequence( $seqId, $title, $length, $molecule, $class );

    #Add the seq data import  ( $seq, $format, $source, $id, $identifier ) = @_;
    my $newSeqDataImport = $bsml->createAndAddSeqDataImport( $bsmlSequence, 'fasta', $fsa, '', $ident );
    &_die("Can't make a new seq data import") unless($newSeqDataImport);
    
    # ($elem, $rel, $href, $role) = @_
    my $link = $bsml->createAndAddLink( $bsmlSequence, 'analysis', "#${analysisId}_analysis", 'input_of');
    
    my $newLink;
    if($class eq 'assembly') {
        $newLink = $bsml->createAndAddLink( $bsmlSequence, 'genome', "#$genomeId" );
    }
    

    return $bsmlSequence;

    
}

sub determineLength {
    my ($fsaFile, $identifier) = @_;
    my $length;

    
}

sub findFilesToCombine {
    my ($id, $otherFiles) = @_;
    my @retval;

    foreach my $otherFile ( @{$otherFiles} ) {
        push(@retval, $otherFile) if($otherFile =~ /$id\./);
    }

    return @retval;
}

sub getInputFiles {
    my $inputStr = shift;
    my @files;
    
    #comma unseparate
    my @tokens = split( /,/, $inputStr);

    foreach my $token ( @tokens ) {
        my $tfh;
        $token =~ s/\.gz$//;
        if( -e $token ) {
            open($tfh, "< $token") or die("Can't open $token ($!)");
        } elsif( -e $token.".gz" ) {
            my $compFile = $token.".gz";
            open($tfh, "<:gzip", "$compFile") or die("Can't open $compFile ($!)");
        } else {
            die("Can't find file $token");
        }
        

        my $isList = 0;

        while( my $line = <$tfh> ) {
            next if($line =~ /^\s+$/);
            
            #It's a bsml file.
            if($line =~ /</ && !$isList) {
                push(@files, $token);
                last;
            #It's a list.
            } elsif($line =~ m|\.bsml|) {
                chomp $line;
                &_die("$line doesn't exist (from list $token) ($!)")
                    unless(-e $line || -e $line.".gz");
                push(@files, $line);
                $isList = 1;
            }
            
        }

    }

    return @files;
    
}

sub process_fasta_file_or_list {
    my ($file, $class) = @_;
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
                croak("Could not determine the format of file $file. Required either ".
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

sub lookup_seq_data_import_info {
    my ($id) = @_;
    my $info;
    if( exists( $seq_data_import_info{$id} ) ) {
        $info = $seq_data_import_info{$id};
    }
    return $info;
}

sub check_parameters {
    my $opts = shift;

    my $error = "";

    &_pod if($opts->{'help'});

    my $input;
    if($opts->{'input_bsml'}) {
        $input = $opts->{'input_bsml'};
    } else {
        $error = "Option input_bsml is required\n";
    }

    @inputFiles = &getInputFiles($input) if( $input );

    if($opts->{'other_bsml_lists'}) {
        $input = $opts->{'other_bsml_lists'};
    }
    @otherFiles = &getInputFiles($input) if( $input );

    unless($opts->{'output'}) {
        $error .= "Option output is required\n";
    } else {
        $output = $opts->{'output'};
    }
    
    if($opts->{'debug'}) {
        $debug = $opts->{'debug'};
    }

    if($opts->{'correct_mistakes'}) {
        $c = 1;
    }

    if($opts->{'cog_search_bsml'}) {
        my $cogs;

        my $file = $opts->{'cog_search_bsml'};
        $file =~ s/\.gz$//;
        if( -e $file  ) {
            open($cogs, "< $file") or $logger->logdie("Can't open cog_search_bsml $file ($!)");
        } elsif( -e $file.".gz" ) {
            open($cogs, "<:gzip", "$file.gz") or $logger->logdie("Can't open cog_search_bsml $file.gz ($!)");
        } else {
            $logger->logdie("Can't find $file (or a gzipped version)");
        }
        

        while (my $line = <$cogs>) {
            chomp $line;
            push(@cogBsmls, $line);
        }
        close($cogs);

        $cogLookup = $opts->{'cog_lookup'} if($opts->{'cog_lookup'});
        $error .= "Option --cog_lookup is required when using --cog_search_bsml\n" unless($opts->{'cog_lookup'});
    }

    if($opts->{'locus_prefix'}) {
        $locusPrefix = $opts->{'locus_prefix'};
    }

    if( $opts->{'organism'} ) {
        $organism = [ $1, $2 ] if($opts->{'organism'} =~ /(\S+)\s+(.*)/);
        $error .= "Coult not parse genus and species out of organism: ".$opts->{'organism'}."\n" 
            if( @{$organism} != 2 );
    }

    if( $opts->{'translation_table'} ) {
        $transTable = $opts->{'translation_table'};
    } else {
        $transTable = 11;
    }

    
    if( $opts->{'cds_fasta'} ) {
        my $option = $opts->{'cds_fasta'};
        my @files = split(/,\s+/, $option );
        print "Processing CDS fasta file[s]\n";
        foreach my $file ( @files ) {
            &process_fasta_file_or_list( $file, 'CDS' );
        }
    }

    if( $opts->{'polypeptide_fasta'} ) {
        my $option = $opts->{'polypeptide_fasta'};
        my @files = split(/,\s+/, $option );
        print "Processing polypeptide fasta file[s]\n";
        foreach my $file ( @files ) {
            &process_fasta_file_or_list( $file, 'polypeptide' );
        }
        
    }
    
    unless($error eq "") {
        &_die($error);
    }
    
}

sub _pod {   
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

sub _die {
    my $msg = shift;
    $logger->logdie($msg);
}
