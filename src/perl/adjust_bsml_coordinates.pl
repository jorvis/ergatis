#!/usr/local/packages/perl-5.8.5/bin/perl

use lib (@INC, "/usr/local/devel/ANNOTATION/ard/current/lib/5.8.8/");

# eval 'exec /usr/local/packages/perl-5.8.5/bin/perl  -S $0 ${1+"$@"}'
#     if 0; # not running under some shell

# BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
# use lib (@INC,$ENV{"PERL_MOD_DIR"});
# no lib "$ENV{PERL_MOD_DIR}/i686-linux";
# no lib ".";

=head1  NAME 

adjust_bsml_coordinates.pl - adjust the positional coordinates within a BSML file

=head1 SYNOPSIS

USAGE:  adjust_bsml_coordinates.pl
            --map_file=/path/to/somemapfile.bsml
            --map_dir=/or/path/to/somedir
            --list_file=/path/to/somefile.list
            --input_file=/paht/to/somefile.bsml
            --output_dir=/path/to/somedir
            --removed_log=/path/to/some/removed.log
          [ --filter_ends=1|0
            --debug=4
            --log=/path/to/somefile.log
          ]

=head1 OPTIONS

B<--map_file,-m> 
    Input BSML map file.

B<--map_dir,-d> 
    Directory of input BSML map files.

B<--list_file,-i> 
    Path to a list of analysis BSML output files.

B<--input_file,-n>
    Path to the input file to be adjusted and searched for overlap

B<--output_list,-u>
    Optional.  If passed, will create an output list with the full paths to each of the 
    BSML files created by this script.

B<--output_dir,-o> 
    Directory where output BSML files will be written.

B<--removed_log,-r>
    The log file where duplicate elements will be printed to after
    removal.

B<--filter_ends,-f> 
    NOT YET IMPLEMENTED

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--help,-h> 
    This help message

=head1   DESCRIPTION

Because some sequences are too large to be analyzed directly by some analysis tools, 
it is sometimes necessary to first split those sequences into pieces, perform the 
analysis on each fragment, and then stitch them back together. This last step causes 
problems because any coordinates given in the result file will represent the 
coordinates of the query in the fragment, and not be representative of the coordinates 
of the original sequence.

The purpose of this script is to read through the BSML results of such an analysis and 
adjust the coordinates to reflect those in the original input file.  This is done using
a set of BSML output files along with a mapping file.  Read the INPUT section for more
information.

Additionally, this script will search neighboring fragments of the input_file fragment
and search for overlapping location specific predictions.  If a prediction is matched
it is removed from the BSML file.  Also, if a location is contained completely within another,
the smaller of the locations is removed.

To use this component, you would most likely have used a splitting component previously
in the pipeline, such as split_fasta, which generates the input maps.

The list of attributes that are adjusted are defined within the script, and are currently
refstart, refend, refpos, start, startpos, endpos, sitepos.

=head1 INPUT

Three types of input files are required to for this to work.  The first is a list of the
BSML result files from an analysis, defined using --list_file.  The list file should
look something like this:

    /some/path/cpa1.assem.1.0.aat_aa.bsml
    /some/path/cpa1.assem.1.1.aat_aa.bsml
    /some/path/cpa1.assem.1.10.aat_aa.bsml
    /some/path/cpa1.assem.1.11.aat_aa.bsml

In this case, an aat_aa analysis was run on 4 sequences, probably named like cpa1.assem.1.0

The BSML mapping file, passed with the --map_file option, gives information about how the
coordinates in these files relate to the original non-split file, here named cpa1.assem.1

If mapping information for the files in your list come from multiple BSML maps, you can use
the --map_dir option.  This will load all maps in the specified directory whose names end
in '.map.bsml'

Also note this script handles compressed input.  If the input BSML files have been gzipped
and end in either the '.gz' or '.gzip' extension, they will be decompressed, modified and
recompressed on the fly.

=head1 OUTPUT

The output BSML files will be written to the directory specified by --output_dir.  You may
specify the output_dir to be the same location as the original files to overwrite them.

As each BSML file is analyzed, a few changes are made to reflect the coordinates back onto
the non-fragmented sequence.  The following elements are searched and the attributes listed
are modified:

    <Aligned-sequence> - start
    <Interval-loc> - startpos, endpos
    <Seq-pair-alignment> - refseq, refxref, refstart, refend
    <Seq-pair-run> - refpos
    <Site-loc> - sitepos

In addition to these, any references to the fragment ID are changed back to the original
sequence.  E.g.

    cpa1.assem.1.0  -->  cpa1.assem.1

Also, an Analysis element is added to reflect that this script modified the file.

Note that the file names, such as cpa1.assem.1.0.aat_aa.bsml, are not changed, though
they no longer contain information directly concerning the fragment.  This way these
sequences can be loaded as they are - else we would have to concatenate each of the
fragmented files into one large BSML file before loading, which is not desirable.

NOTE:  The BSML output currently has a return character between elements that are
analyzed.  This is an artifact of XML::Twig's printing methods that I have yet to 
resolve.  

=head1 CONTACT

    Kevin Galens
    kgalens@tigr.org

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
BEGIN {
use Workflow::Logger;
}
use XML::Twig;
use File::Basename;
use SeqLocation::SeqLocation;
use Data::Dumper;

#######
## ubiquitous options parsing and logger creation
my %options = ();
my $results = GetOptions (\%options, 
                            'map_file|m=s',
                            'map_dir|d=s',
                            'list_file|i=s',
                            'input_file|n=s',
                            'removed_log|r=s',
                            'output_dir|o=s',
                            'filter_ends|f=s',
                            'debug=s',
                            'log|l=s',
                            'help|h') || pod2usage();

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = Workflow::Logger::get_logger();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

&check_parameters(\%options);

## how many maps do we have?
my @maps;
if (defined $options{map_file}) {
    $logger->debug("adding $options{map_file} to map list");
    push @maps, $options{map_file};
}

if (defined $options{map_dir}) {
    for (glob "$options{map_dir}/*.map.bsml") {
        $logger->debug("adding $_ to map list");
        push @maps, $_;
    }
}

my $output_dir = $options{output_dir};

## load the analysis info from the map(s)
my %analysis;
my %sequence_map;
my ($map, $mtwig);
my $analysisId;

for $map (@maps) {
    $mtwig = XML::Twig->new(
                            twig_roots => {
                                            'Sequence' => \&process_map_sequence,
                                            'Analysis' => \&process_map_analysis
                                          }
                           );
    $analysisId=$map;
    $mtwig->parsefile($map);
    $mtwig->purge;  ## free memory
}

## we want to loop over a list of bsml result files, each of which should have a
##  record of the same id in the mapping file.
## first get the list of files
my @bsml_files;
open(my $lfh, "<$options{list_file}") || $logger->logdie("can't open list file: $!");
for (<$lfh>) {
    chomp;
    push @bsml_files, $_;
}

## open an output list file if the user requested it
my $olfh;

my $spaCounter = 0;

## loop through each BSML file and adjust coordinates, adding an Analysis element to each
my $adjustment;
my $ofh;
my $analyses_found = 0;
my $feat_id;  ## like cpa1.assem.2.1

## some analysis types (such as AAT) have chain numbers which need to be refactored
##  along with the coordinate adjustment.  this hash stores the next available numbers
##  for each feat_name
my %next_chain_num;

## the first sequence in each of these bsml file should be a sequence stub for the query sequence.
##  it's ID (after changing the reference back to the original) goes in here.
my $root_seq_stub;  ## like cpa1.assem.2

#If we are parsing an adjacent fragment, this flag is one.
my $adjFlag = 0;

#Some constants.  Pseudo enumeration
my ($CUR, $PREV, $NEXT) = (-1, 0, 1);

#Open the removed.log file handle.
my $rFH;
open($rFH, "> $options{removed_log}") || 
    $logger->logdie("Unable to open $options{removed_log} ($!)");

my @adjSeqLoc;

my $bf = $options{input_file};

$adjSeqLoc[$PREV] = new SeqLocation::SeqLocation('prev', 'adj');
$adjSeqLoc[$NEXT] = new SeqLocation::SeqLocation('next', 'adj');
&findAdjacent($bf, \@bsml_files);

## get the filename
my $fname = basename($bf);
$root_seq_stub = '';

## it is assumed that files are named like featid.analysiscomponent.bsml (eg. cpa1.assem.2.1.aat_aa.bsml)
##  we can pull the source ID from this file then using a regex (here, cpa1.assem.2.1)
$fname =~ /(.+\d)\..+\.bsml/ || $logger->logdie("$fname does not match naming convention.  can't extract id");
$feat_id = $1;

$logger->debug("using $feat_id as feat_id in file $bf") if ($logger->is_debug);

## get the adjustment
$adjustment = $sequence_map{$feat_id}{offset};

$logger->debug("setting adjustment as $adjustment for $feat_id in file $bf") if ($logger->is_debug);

## open the input and output files.  how we do this depends on whether the input was zipped or not
my $ifh;
if ($bf =~ /\.(gz|gzip)$/) {
    open ($ifh, "<:gzip", $bf)                       || $logger->logdie("can't read zipped input file '$bf': $!");
    open ($ofh, ">:gzip", "$output_dir/$fname.part") || $logger->logdie("can't create output file: $!");
} else {
    open ($ifh, "<$bf")                     || $logger->logdie("can't read input file $bf: $!");
    open ($ofh, ">$output_dir/$fname.part") || $logger->logdie("can't create output file: $!");
}

my $twig = XML::Twig->new(
                          twig_roots               => {
                              'Seq-pair-alignment' => \&processSeqPairAlignment,
                              'Aligned-sequence'   => \&processAlignedSequence,
                              'Feature-tables'     => \&processFeatureTables,
                              #'Interval-loc'       => \&processIntervalLoc,
                              #'Site-loc'           => \&processSiteLoc,
                              'Analyses'           => \&processAnalyses,
                              'Sequence'           => \&processSequence,
                          },
                          twig_print_outside_roots => $ofh,
                          pretty_print => 'indented',
                          );

## do the parse
#$twig->parsefile($bf);
$adjFlag = $CUR;
$twig->parse($ifh);

## error if we didn't find an analysis
if ($analyses_found) {
    $analyses_found = 0;
} else {
    $logger->error("Analysis element not found in $bf") if ($logger->is_error);
}

## mv the temp file over the target (can't read and write to same file)
system("mv $output_dir/$fname.part $output_dir/$fname");

exit(0);

sub check_parameters {
    my ($options) = @_;

    ## map_dir or map_file must have been passed
    unless ( $options{map_dir} || $options{map_file} ) {
        $logger->logdie("map_dir or map_file must be passed");
    }

    ## check the map dir, if passed
    if ( $options{map_dir} && ! -e $options{map_dir} ) {
        $logger->logdie("map_dir was passed but does not exist");
    }

    ## check map file
    if ( $options{map_file} && ! -e $options{map_file} ) {
        $logger->logdie("map_file either not passed or does not exist");
    }
    
    ## check list file
    unless ( $options{list_file} && -e $options{list_file} ) {
        $logger->logdie("list_file either not passed or does not exist");
    }    

    ## output_dir is required
    if (! $options{output_dir} ) {
        $logger->logdie("output_dir option is required");
    }
    
    ## check input file
    if( !($options{input_file}) || !(-e $options{input_file}) ) {
        if(!(-e $options{input_file}) && -e $options{input_file}.".gz") {
            $options{input_file}.=".gz";
        } else {
            $logger->logdie("input_file either not passed or does not exist");
        }
        
    } 

    unless ($options{removed_log}) {
        $logger->logdie("removed_log is required");
    }
    
    ## set some defaults
    $options{output_list} = 0 unless($options{output_list});

    if(0){
        pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});    
    }
}

sub process_map_analysis {
    my ($mtwig, $Analysis) = @_;
    
    for my $Attribute ( $Analysis->children('Attribute') ) {
        ## save each of these
        $analysis{$analysisId}{$Attribute->{att}->{name}} = $Attribute->{att}->{content};
    }
}

sub process_map_sequence {
    my ($mtwig, $Sequence)  = @_;

    ## get the id of this sequence
    my $id = $Sequence->{att}->{id} || $logger->logdie("Sequence element missing id!");

    ## we're only considering tiling path fragments
    if ($Sequence->{att}->{class} ne 'tiling_path_fragment') {
        $logger->info("skipping sequence $id - not a tiling_path_fragment") if ($logger->is_info);
        return 0;
    }

    $sequence_map{$id}{length} = $Sequence->{att}->{length};

    ## there should only be one Numbering, get the refnum (offset) from it
    my $offset = $Sequence->first_child('Numbering')->{att}->{refnum};
    if (! defined $offset) { $logger->logdie("sequence $id missing Numbering->refnum") }
    $sequence_map{$id}{offset} = $offset;
}


sub processAlignedSequence { 
    my ($twig, $element) = @_;

    $logger->logdie("This script does not handle the Aligned-sequence tag yet");
    
    ## need to replace any start attributes
    if (defined $element->{att}->{start}) {
        $element->{att}->{start} += $adjustment;
    }
    
    $element->print($ofh);
}

sub processAnalyses {
    my ($twig, $element) = @_;
    
    $analyses_found++;
    
    ## create an Analysis element
    my $analysis = XML::Twig::Elt->new('Analysis');
    $analysis->set_att('id', 'adjust_bsml_coordinates_analysis');
    
    ## add each option
    for my $option ( keys %options ) {
        my $attribute = XML::Twig::Elt->new('Attribute');
        $attribute->set_att('name', $option);
        $attribute->set_att('content', $options{$option});
        $attribute->paste('last_child', $analysis);
    }
    
    ## add this Analysis to the Analyses group
    $analysis->paste('last_child', $element);
    
    $element->print($ofh);
}

# sub processIntervalLoc { 
#     my ($twig, $element) = @_;
    
#     ## need to replace any startpos or endpos attributes
#     for my $attribute ( qw(startpos endpos) ) {
#         if (defined $element->{att}->{$attribute}) {
#             $element->{att}->{$attribute} += $adjustment;
#         }
#     }
    
#     ## don't print if we are within a Feature.  Since Features
#     #   are within a sequence, we'll just print twice.
#     if (($twig->context)[-1] eq 'Feature') {
#         $logger->debug("not printing Interval-loc since it is within a Feature") if ($logger->is_debug);
#     } else {
#         $element->print($ofh);
#     }
# }

sub processFeatureTables {
    my ($twig, $Featuretables) = @_;
    my (@featStats);

    my @featStats = ('','Feature','','',{});
    my %numFeats;

    #When a feature is removed, we must keep record of it for printing and
    #also to check for its presence in a feature group.
    my %removedFeats;

    for my $Featuretable($Featuretables->children('Feature-table')) {
        my $featTableID;
        if(defined($Featuretable->{att}->{id})) {
            $featTableID = $Featuretable->{att}->{id};
            $numFeats{$featTableID} = [0, $Featuretable];
        } else { 
            $logger->logdie("No id attribute for Feature-table element");
        }
    
        for my $Feature($Featuretable->children('Feature')) {
            if(defined($Feature->{att}->{id})) {
                $featStats[0] = $Feature->{att}->{id};
                $numFeats{$featTableID}->[0]++;
            } else {
                $logger->logdie("No id for Feature element");
            }


            if(defined($Feature->{att}->{class})) {
                $featStats[4]->{class} = $Feature->{att}->{class};
            } else {
                $logger->logdie("Feature is missing class attribute");
            }

            my $Intervalloc = $Feature->first_child('Interval-loc');

            $logger->logdie("Could not retrieve Interval-loc") 
                if(!$Intervalloc);

            foreach my $att(qw/startpos endpos/) {
                if(defined($Intervalloc->{att}->{$att})) {
                    $Intervalloc->{att}->{$att} += $adjustment;
                } else {
                    $logger->logdie("Interval-loc: $att was not defined");
                }
                
                $featStats[2] = $Intervalloc->{att}->{$att} if($att eq 'startpos');
                $featStats[3] = $Intervalloc->{att}->{$att} if($att eq 'endpos');
                
            }

            $featStats[4]->{complement} = $Intervalloc->{att}->{complement};

            if(!($adjFlag == $CUR)) {
                $adjSeqLoc[$adjFlag]->addSeqLocation($featStats[0], $featStats[1],$featStats[2],$featStats[3] );
            } else {
                if($adjSeqLoc[$PREV] && $adjSeqLoc[$PREV]->checkOverlap(@featStats[1..4],'prev')) {
                    $removedFeats{$featStats[0]} = [\@featStats, $Feature];
                    $adjSeqLoc[$PREV]->removeSeqLocation(@featStats[0..3]);
                } elsif($adjSeqLoc[$NEXT] && $adjSeqLoc[$NEXT]->checkOverlap(@featStats[1..4],'next')) {
                    $removedFeats{$featStats[0]} = [\@featStats, $Feature];
                    $adjSeqLoc[$NEXT]->removeSeqLocation(@featStats[0..3]);
                } 
            }
        }

    }

    if($adjFlag == $CUR) {
        my $numFeaturemems = 0;
        my $delFeatMems = 0;
        my @deletedFM;

        foreach my $Featuregroup($Featuretables->children('Feature-group')) {
            $numFeaturemems = scalar $Featuregroup->children('Feature-group-member');
            foreach my $Featuregroupmember($Featuregroup->children('Feature-group-member')) {
                if(defined($Featuregroupmember->{att}->{featref}) &&
                   defined($removedFeats{$Featuregroupmember->{att}->{featref}}) ) {
                    push(@deletedFM, $Featuregroupmember);
                    $delFeatMems++;
                } else {
                  #   print STDERR "featxref not defined\n" 
#                         if(!defined($Featuregroupmember->{att}->{featref}));
#                     print STDERR "didn't find featxref in removed\n"
#                         if(defined($removedFeats{$Featuregroupmember->{att}->{featref}}) );
                }
            }
            if($delFeatMems == $numFeaturemems) {
                $Featuregroup->print($rFH);
                $Featuregroup->delete();
            } else {
                foreach my $delThis(@deletedFM) {
                    $delThis->print($rFH);
                    $delThis->delete();
                }
            }
            $delFeatMems = 0;
            $numFeaturemems = 0;
            @deletedFM = ();
                
        }
        foreach my $featID(keys %removedFeats) {
            $removedFeats{$featID}->[1]->print($rFH);
            $removedFeats{$featID}->[1]->delete();
        }

    }


}

sub processSeqPairAlignment {
    my ($twig, $spa) = @_;
    ## spa = SeqPairAligment
    my %properties = {};
    my $sprsDeleted = 0;
    my @removedSprs;
    my $totSprs = 0;
    $spaCounter++;
    
    ## need to switch the IDs in refseq and refxref with the root stub
    for my $attribute ( qw(refseq refxref) ) {
        if (defined $spa->{att}->{$attribute}) {
            $spa->{att}->{$attribute} =~ s/$feat_id/$root_seq_stub/;
        }
    }
    
    ## needs to check for refstart refend
    ## (compstart and compend refer to the subject sequence and are skipped)
    for my $attribute ( qw(refstart refend) ) {
        if (defined $spa->{att}->{$attribute}) {
            $spa->{att}->{$attribute} += $adjustment;
        }
    }
    if(defined $spa->{att}->{compseq}) {
        $properties{compseq} = $spa->{att}->{compseq};
    }
    
    my $sprID = 0;
    my @sprStats = ();

    ## reserve next chain number for this Seq-pair-alignment
    #   this may only be used by some components, such as AAT
    my $chain_num = ++$next_chain_num{$root_seq_stub};
    $logger->debug("chain_num $chain_num reserved for root_seq_stub $root_seq_stub (feat_id $feat_id)") if ($logger->is_debug);
    
    ##  then check each child Seq-pair-run (spr) for refpos
    ##  (comppos refers to the subject sequence and is skipped)
    for my $spr ( $spa->children('Seq-pair-run') ) {
        $totSprs++;
        @sprStats = ($sprID, 'Seq-pair-run');

        if (defined $spr->{att}->{refpos}) {
            $spr->{att}->{refpos} += $adjustment;
            $sprStats[2] = $spr->{att}->{refpos};
            $sprStats[3] = $sprStats[2] +  $spr->{att}->{runlength};
        }

        
        foreach(qw( refcomplement comppos comprunlength compcomplement)) {
            if(defined($spr->{att}->{$_})) {
                $properties{$_} = $spr->{att}->{$_};
            }
        }

        $sprStats[4] = \%properties;

        #print STDERR "\nFlag :: $adjFlag\nCur :: $CUR\n";

        if(!($adjFlag == $CUR)) {
            my $sprStart = $spr->{att}->{refpos};
            my $sprEnd = $sprStart+$spr->{att}->{runlength};
            $adjSeqLoc[$adjFlag]->addSeqLocation(@sprStats[0..3]);
            $sprID++;
        } else {
          #   print STDERR "In the else at least\n";
#             print STDERR "next isn't defined\n" if(!$adjSeqLoc[$NEXT]);
#             print STDERR "prev isn't defined\n" if(!$adjSeqLoc[$PREV]);
#             print STDERR "@sprStats\n";
#             print STDERR Dumper($adjSeqLoc[$PREV]);
#             exit(0);
            if($adjSeqLoc[$NEXT] && $adjSeqLoc[$NEXT]->checkOverlap(@sprStats[1..4], 'next')) {
                $adjSeqLoc[$NEXT]->removeSeqLocation(@sprStats[0..3]);
                push @removedSprs, $spr; 
                #$spr->print($rFH);
                $sprsDeleted++;
            } elsif($adjSeqLoc[$PREV] && $adjSeqLoc[$PREV]->checkOverlap(@sprStats[1..4], 'prev')) { 
                $adjSeqLoc[$PREV]->removeSeqLocation(@sprStats[0..3]);
                push @removedSprs, $spr; 
                #$spr->print($rFH);
                $sprsDeleted++;
            }
        }
        
        ## check for a chain_number attribute
        for my $att ( $spr->children('Attribute') ) {
            if ($att->{att}->{name} && $att->{att}->{name} eq 'chain_number') {
                $logger->debug("reassigned chain_number from " . $att->{att}->{content} . " to $chain_num") if ($logger->is_debug);
                $att->{att}->{content} = $chain_num;
            }
        }
    }

    if($adjFlag == $CUR) {

        if ($sprsDeleted == $totSprs) {
            $spa->print($rFH); 
        } else {
            foreach (@removedSprs) {
                $_->print($rFH);
                $_->delete();
            }
            $spa->print($ofh);
        }

    }
}

sub processSequence { 
    my ($twig, $element) = @_;
    
    ## if this is the first Sequence, replace the ID.
    ## the id should be $feat_id
    unless ( $root_seq_stub ) {
        if ($element->{att}->{id} eq $feat_id) {
            
            ## strip the partial id out of the title (if there is one).  it should be the first word
            $element->{att}->{title} =~ s/^$feat_id //;

            ## the root_seq_stub should be the id up to and excluding the last .N
            #   capturing like (aa1.assembly.17619).0
            $feat_id =~ /^(.+)\.[0-9]+$/;
            
            $root_seq_stub = $1;
            
            ## get the new id (now first word of title)
            #$root_seq_stub = (split('\s', $element->{att}->{title}))[0];

            ## set the id
            $logger->debug("changing id from $element->{att}->{id} to $root_seq_stub") if ($logger->is_debug);
            $element->{att}->{id} = $root_seq_stub;
            
            
        } else {
            $logger->logdie($element->{att}->{id} . " doesn't match $feat_id");
        }
    }
    
    $element->print($ofh);
}

# sub processSiteLoc { 
#     my ($twig, $element) = @_;
    
#     ## need to replace any sitepos
#     if (defined $element->{att}->{sitepos}) {
#         $element->{att}->{sitepos} += $adjustment;
#     }
    
#     $element->print($ofh);
# }

sub findAdjacent {
    my ($infile, $bsmlList) = @_;

    my @prevAndNext;

    #Parse out this number
    ($infile =~ /.*\.(\d+)\.(\d+)\.[\w-.]+\.bsml/) 
        || $logger->logdie("Could not extract fragment number from $infile");
    my $fragNo = $2;
    my $assembl = $1;
    foreach my $bsmlFile(@{$bsmlList}) {
        if(($bsmlFile =~ /.*\.(\d+)\.(\d+)\.[\w-.]+\.bsml/) && $1 == $assembl && ($2-$fragNo)**2 == 1) {
            $prevAndNext[(($2-$fragNo+1)/2)] = $bsmlFile;
        }
    }
    

    #print STDERR "@prevAndNext\n";
    #exit(0);

    if(defined($prevAndNext[$PREV]) ) {
        my @over = &overlap($prevAndNext[$PREV], 'prev');
        $adjSeqLoc[$PREV]->setOverlapRange(&overlap($prevAndNext[$PREV], 'prev'));
        $adjFlag = $PREV;
        &processAdjacent($prevAndNext[$PREV]);
    } else {
        $adjSeqLoc[$PREV] = 0;
    }

    if(defined($prevAndNext[$NEXT]) ) {
        $adjSeqLoc[$NEXT]->setOverlapRange(&overlap($prevAndNext[$NEXT], 'next'));
        $adjFlag = $NEXT;
        &processAdjacent($prevAndNext[$NEXT]);
    } else {
        $adjSeqLoc[$NEXT] = 0;
    }

}

sub processAdjacent {
    my $adjacentFile = shift;
    my $aFH;
    my $null;
    open($null, "> /dev/null");

    #Were only doing SeqPairAlignments for now
    my $twig = XML::Twig->new(
                              twig_roots               => {
                                  'Seq-pair-alignment' => \&processSeqPairAlignment,
                                  'Feature-tables'     => \&processFeatureTables,
                              },
                              twig_print_outside_roots => $null,
                              );

    $adjacentFile.='.gz' unless(-e $adjacentFile);

    $logger->logdie("$adjacentFile doesn't exist") unless(-e $adjacentFile);

    if ($adjacentFile =~ /\.(gz|gzip)$/) {
        open($aFH, "<:gzip", $adjacentFile) || $logger->logdie("can't read zipped input file '$adjacentFile': $!");
    } else {
        open($aFH, "$adjacentFile") || $logger->logdie("can't read zipped input file '$adjacentFile': $!");
    }
    
    $twig->parse($aFH);
}

sub overlap {
    my ($file, $adj) = @_;
    ($file =~ /.*\/(.*)(\d+)\.[\w-.]+\.bsml/) 
        || $logger->logdie("$file doesn't match naming scheme.  Can't use name to find map file\n");
    my $mapAnID = "$options{map_dir}/$1map.bsml";
    my $length = $sequence_map{"$1$2"}{length};
    my $overlapLen = $analysis{$mapAnID}{overlap_length};
    $adjustment = $sequence_map{"$1$2"}{offset};
    my ($overStart, $overEnd);
    if($adj eq 'prev') {
        $overStart = $adjustment+($length-$overlapLen);
        $overEnd = $adjustment+$length;
    } elsif($adj eq 'next') {
        $overStart = $adjustment;
        $overEnd = $adjustment+$overlapLen;
    } else {
        $logger->logdie("overlap requires either prev or next as second argument.  $adj was passed");
    }

    return($overStart,$overEnd);
    
}
