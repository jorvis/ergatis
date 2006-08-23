#!/local/packages/perl-5.8.8/bin/perl

BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

    parse_for_ncRNAs.pl - parse and combine the output for blast searches
    (using 16s and 23s rRNA), tRNAscan-SE and infernal (covariance model
    search) of Rfam CMs.

=head1 CONTACT

    Kevin Galens (kgalens@tigr.org)

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use BSML::BsmlBuilder;
use Workflow::Logger;
use Workflow::IdGenerator;
use Class::Struct;
use XML::Twig;
use Data::Dumper;


################GLOBALS AND CONSTANTS###################
struct( RNAFEATURE => { id => '$', 
                        start => '$', 
                        end => '$',
                        title => '$',
                        class => '$'} );
my %rnaFeats = ();
my $output;
my $bsmlDoc = new BSML::BsmlBuilder();
my @bsmlLists;
my @inputFiles;
my @bsmlToCombine;
my $idGenerator;
my $project;
my $asmblId;
########################################################

my %options = ();
my $results = GetOptions (\%options, 
                          'input_list|i=s',
                          'input_file|f=s',
                          'bsml_lists|b=s',
                          'output|o=s',
                          'project|p=s',
                          'id_repository|r=s',
                          'log|l=s',
                          'debug|d=s',
                          'help|h') || &_pod;

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

# display documentation
&_pod if( $options{'help'} );

# make sure all passed options are peachy
&check_parameters(\%options);


foreach my $inFile (@inputFiles) {
    @bsmlToCombine = &findBsmlFiles($inFile, \@bsmlLists);
    
    foreach my $bsmlFile (@bsmlToCombine) {
        &parseBsml($bsmlFile, \%rnaFeats);
    }
}

#Write the bsml file
&writeRnaBsmlFile($output, \%rnaFeats);


##############################SUB ROUTINES##############################################

#Uses the bsml_lists to find all the bsml files that are predictions made from the fsa
#input_file sequence
sub findBsmlFiles {
    my ($inputFile, $bsmlLists) =@_;
    my $retval;

    my $baseName = $1 if($inputFile =~ m|([^/]+)\.fsa$|);
    &_die("Could not parse base name from $inputFile") unless($baseName);
    $asmblId = $baseName;
    my $escapedName = $baseName;
    $escapedName =~ s/\./\\./g;

    foreach my $list (@{$bsmlLists}) {
        my $fileFound;

        open(IN, "< $list");
        my @files = <IN>;
        close(IN);
        
        foreach my $file (@files) {
            chomp $file;
            $fileFound = $file if($file =~ /$escapedName/);
        }
             
        &_die("Could not find a match for $baseName in $list") unless($fileFound);
        push(@{$retval}, $fileFound);

    }

    return @{$retval};
}

sub parseBsml {
    my ($bsmlFile, $rnaFeats) = @_;
    my $onlyFeatures = 0;
    

    my $twig = XML::Twig->new(twig_roots => {
        'Feature'  => sub { &parseFeat(@_, $rnaFeats, \$onlyFeatures) },
        'Seq-pair-alignment' => sub { &parseSpa(@_, $rnaFeats) unless($onlyFeatures); }});

    $twig->parsefile($bsmlFile);
   
}

sub parseFeat {
    my ($twig, $elt, $rnaFeats) = @_;
    my $tmpRna;
    
    if($elt->{'att'}->{'class'} =~ /(ncRNA|tRNA)/ ) {
        $tmpRna = new RNAFEATURE;
        $tmpRna->id($idGenerator->next_id( project => $project, type => $elt->{'att'}->{'class'}));        
        $tmpRna->class($elt->{'att'}->{'class'});
        $tmpRna->title($elt->{'att'}->{'title'}) if(exists($elt->{'att'}->{'title'}));

        #Retrieve the assocatiated interval loc
        my $intLoc = $elt->first_child('Interval-loc');
        &_die("Feature element $elt->{'att'}->{'id'} does not have a child Interval-loc element")
            unless($intLoc);
        $tmpRna->start($intLoc->{'att'}->{'startpos'}) if($intLoc->{'att'}->{'startpos'});
        $tmpRna->end($intLoc->{'att'}->{'endpos'}) if($intLoc->{'att'}->{'endpos'});

        &addRnaFeat($tmpRna, $rnaFeats);

    }

    

}

sub parseSpa {
    my ($twig, $spa, $rnaFeats) = @_;
    my $tmpRna = new RNAFEATURE;
    my ($max, $highSpr) = (0,0);
    
    #We don't really know which type
    $tmpRna->id($idGenerator->next_id( project => $project, type => 'ncRNA' ));
    $tmpRna->class('ncRNA');
    $tmpRna->title($spa->{'att'}->{'compxref'});

    
    
    foreach my $spr($spa->children('Seq-pair-run')) {
        unless($spr->{'att'}->{'runscore'}) {
            $highSpr = $spr;
            last;
        }

        #Find the highest scoring SPR. (For the blast output)
        if($spr->{'att'}->{'runscore'} > $max) {
            $max = $spr->{'att'}->{'runscore'};
            $highSpr = $spr;
        }        
    }

    $tmpRna->start($highSpr->{'att'}->{'refpos'});
    $tmpRna->end($tmpRna->start + $highSpr->{'att'}->{'runlength'});

    &addRnaFeat($tmpRna, $rnaFeats);
   
}

#Checks for overlaps and fixes it.
sub addRnaFeat {
    my ($addRna, $rnaFeats) = @_;
    my $add = 1;

    foreach my $rnaId (keys %{$rnaFeats}) {
        my $tmpRna = $rnaFeats->{$rnaId};
        
        #Is there an overlap (will produce a positive number if so).
        my @sortedBound = sort($tmpRna->start, $tmpRna->end, $addRna->start, $addRna->end);
        my $ol = abs($tmpRna->start - $tmpRna->end) + abs($addRna->start + $addRna->end)
            - abs($sortedBound[0] - $sortedBound[3]);
        
        #If there is an overlap and also if they are on the same strand.
        if($ol > 0 && (($tmpRna->start - $tmpRna->end) * ($addRna->start - $addRna->end) > 0)) {
            my ($s, $e) = $tmpRna->start < $tmpRna->end ? ($sortedBound[0], $sortedBound[3]) :
                ($sortedBound[3], $sortedBound[0]);
            $tmpRna->start($s);
            $tmpRna->end($e);
            $tmpRna->class('ncRNA');
            $tmpRna->title('overlap');
            $add = 0;
            last;
        }

       
    }

    $rnaFeats->{$addRna->id} = $addRna if($add);

}

sub writeRnaBsmlFile {
    my ($outFile, $rnaFeats) = @_;

    #Create the sequence element
    my $class = 'assembly';
    $class = $1 if($asmblId =~ /^.*?\.(\w+)/);
    my $seq = $bsmlDoc->createAndAddSequence( $asmblId, $asmblId, '', 'dna', $class);
    my $ft = $bsmlDoc->createAndAddFeatureTable( $seq );

    while(my ($id, $rna) = each(%{$rnaFeats})) {
        my $feat = $bsmlDoc->createAndAddFeature($ft, $rna->id, $rna->title, $rna->class);
        $feat->addBsmlLink('analysis', '#parse_for_ncRNA_analysis', 'computed_by');
        my $strand = $rna->start < $rna->end ? 0 : 1;
        $feat->addBsmlIntervalLoc($rna->start, $rna->end, $strand);
    }

    $bsmlDoc->createAndAddAnalysis( id => "parse_for_ncRNA_analysis",
                                    sourcename => $options{'input_file'}
                                    );

    $bsmlDoc->write($outFile);
}

sub check_parameters {

    if($options{'input_file'}) {
        #&_die("$options{input_file} does not exist") unless(-e $options{'input_file'});
        push(@inputFiles, $options{'input_file'});
    } elsif($options{'input_list'}) {
        &_die("$options{'input_list'} does not exist") unless(-e $options{'input_list'});
        open(IN, "< $options{'input_list'}") or &_die("Unable to open $options{'input_list'}".
                                                      " ($!)");
        @inputFiles = <IN>;
        close(IN);
    } else {
        &_die("Either input_list or input_file option must exist");
    }


    if($options{'bsml_lists'}) {
        @bsmlLists = split(/\s/,$options{'bsml_lists'});
    } else {
        &_die("A space separated string of bsml lists must be provided (option bsml_lists)");
    }

    unless($options{'id_repository'}) {
        &_die("Option id_repository is required for id creation.  See Workflow::IdGenerator ".
              "for id_repository information.");
    }
    $idGenerator = new Workflow::IdGenerator( id_repository => "$options{'id_repository'}" );

    unless($options{'project'}) {
        &_die("Option project is required for id creation.");
    }
    $project = $options{'project'};

    unless($options{'output'}) {
        &_die("Option output is required (bsml output file name)");
    }
    $output = $options{'output'};
    

}

sub _pod {
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

sub _die {
    my $msg = shift;
    $logger->logdie($msg);
}






