#!/usr/local/bin/perl
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
            --output_dir=/path/to/someDir/
          [ --log=/path/to/file.log
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
    [Required] Output bsml file.

B<--log,-l>
    [Optional] Logfile.

B<--debug,-d>
    [Optional] Higher number is more verbose.

B<--help,-h>
    [Optional] Print this message

=head1  DESCRIPTION

=head1  INPUT


=head1 OUTPUT


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
########################################
$|++;

my %options = ();
my $results = GetOptions (\%options, 
                          'input_bsml|i=s',
                          'other_bsml_lists|b=s',
                          'output_dir|o=s',
                          'correct_mistakes|c',
                          'log|l=s',
                          'debug=s',
                          'help|h') || &_pod;

#Setup the logger
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

# Check the options.
&check_parameters(\%options);

foreach my $file (@inputFiles) {

    $sourcename = $1 if($file =~ m|(.*\d+_[^/]+)/|);
    &_die("Could not parse the sourcename directory for analysis section from $file") 
        unless($sourcename);
    print STDOUT $sourcename."\n";

    #Make the output file name.
    my $base = $1 if($file =~ m|.*/([^/\.]+\.[^/\.]+\.[^/\.]+)\.[^/]+$|);
    $outputFile = "$output/$base.summary.bsml";
    print "input: $file\n";

    my $bsml = new BSML::BsmlBuilder;

    #Parse the bsml file
    my $twig = new XML::Twig( twig_handlers => {
        'Sequence' => sub { &printSequenceFeatures( $bsml, @_ ) },
        });

    print STDOUT "Parsing $file\n";
    $twig->parsefile($file);

    my @addThese;
    foreach my $otherFile (@otherFiles) {
        push(@addThese, $otherFile) if($otherFile =~ /$base\./);
    }

    foreach my $addThis ( @addThese ) {
        my $newTwig = new XML::Twig( twig_handlers => {
            'Sequence' => sub { &printSequenceFeatures( $bsml, @_ ) }
            });

        print STDOUT "adding $addThis to bsml file\n";
        $newTwig->parsefile($addThis);
    }

    

    print STDOUT "Writing to $outputFile\n";
    $bsml->write($outputFile);

}


######################## SUB ROUTINES #######################################
sub printSequenceFeatures {
    my ($bsml, $twig, $seqElem) = @_;
    
    my $id = $seqElem->att('id');
    &_die("Can't parse id from Sequence in bsml file") unless($id);

    my $bsmlSequence = &addSequence($bsml, $seqElem);

    &_die("bsmlSequence ".$seqElem->att('id')." does not have a class") unless($seqElem->att('class'));

    return unless($seqElem->att('class') eq 'assembly');

    my $asmblId = $id;
    print STDOUT "parsing: $asmblId\n";

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

    foreach my $featElem ( @features ) {
        
        my $bsmlFeature = &addFeature( $bsml, $newFeatureTable, $featElem );
        &_die("Unable to create new bsml feature for sequence $asmblId") unless($bsmlFeature);
    }

    my @fgs = $fTables->children('Feature-group') if($fTables);
    
    foreach my $fg ( @fgs ) {
        my $bsmlFg = &addFeatureGroup( $bsml, $bsmlSequence, $fg );
        &_die("Could not create new Feature-group") unless($bsmlFg);
    }

    my $analysis;
    if($bsml->returnBsmlAnalysisR(0)) {
        $analysis = $bsml->returnBsmlAnalysisR(0);
    } else {
        $analysis = $bsml->createAndAddAnalysis( 'id' => $analysisId."_analysis",
                                                 'program' => 'pipeline_summary',
                                                 'algorithm' => 'pipeline_summary', 
                                                 'sourcename' => $sourcename );
    }
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


    my $geneId;
    foreach my $fgm ( @fgms ) {
        next unless($fgm->att('feature-type') eq 'gene');
        $geneId = $fgm->att('featref');
        &_die("Can't parse featref from Feature-group-member in $groupSet (Feature-group)") unless($geneId);
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

    unless($bsmlFeature->{'BsmlLink'}) {
        $link = $bsmlFeature->{'BsmlLink'}->[0];
    } else {
        $bsml->createAndAddLink( $bsmlFeature, 'analysis', "#${analysisId}_analysis", "input_of" );
        $link = $bsmlFeature->{'BsmlLink'}->[0];
    }
    &_die("Couldn't add or find link element associated with Feature $featId") unless($link);

    #Add attributes if any
    my @attributes = $feature->children('Attribute');

    foreach my $att ( @attributes ) {
        my $bsmlAtt = $bsml->createAndAddBsmlAttribute( $bsmlFeature, $att->att('name'), $att->att('content') );
    }

    $features->{$featId} = $bsmlFeature;
    return $bsmlFeature;

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
        open($tfh, "< $token") or die("Can't open $token ($!)");

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
                    unless(-e $line);
                push(@files, $line);
                $isList = 1;
            }
            
        }

    }

    return @files;
    
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

    @inputFiles = &getInputFiles($input);

    if($opts->{'other_bsml_lists'}) {
        $input = $opts->{'other_bsml_lists'};
    }
    @otherFiles = &getInputFiles($input);

    unless($opts->{'output_dir'}) {
        $error .= "Option output_dir is required\n";
    } else {
        $output = $opts->{'output_dir'};
    }
    
    if($opts->{'debug'}) {
        $debug = $opts->{'debug'};
    }

    if($opts->{'correct_mistakes'}) {
        $c = 1;
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
