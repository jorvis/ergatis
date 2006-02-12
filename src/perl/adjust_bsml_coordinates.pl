#!/usr/local/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

adjust_bsml_coordinates.pl - adjust the positional coordinates within a BSML file

=head1 SYNOPSIS

USAGE:  adjust_bsml_coordinates.pl
            --map_file=/path/to/somemapfile.bsml
            --map_dir=/or/path/to/somedir
            --list_file=/path/to/somefile.list
            --output_dir=/path/to/somedir
          [ --output_list=/path/to/somefile.list
            --output_subdir_size=1000
            --output_subdir_prefix=fasta
            --list_file_glob='*.aat_aa.bsml'
            --filter_ends=1|0
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

B<--list_file_glob,-g> 
    Use to filter the BSML files considered in --list_file.

B<--output_list,-u>
    Optional.  If passed, will create an output list with the full paths to each of the 
    BSML files created by this script.

B<--output_dir,-o> 
    Directory where output BSML files will be written.

B<--output_subdir_size,-z>
    If defined, this script will create numbered subdirectories in the output directory, each
    containing this many sequences files.  Once this limit is reached, another subdirectory
    is created.

B<--output_subdir_prefix,-x>
    To be used along with --output_subdir_size, this allows more control of the names of the
    subdirectories created.  Rather than just incrementing numbers (like 10), each subdirectory 
    will be named with this prefix (like prefix10).

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

To use this component, you would most likely have used a splitting component previously
in the pipeline, such as split_fasta, which generates the input maps.

The list of attributes that are adjusted are defined within the script, and are currently
refstart, refend, refpos, start, startpos, endpos, sitepos.

=head1 INPUT

Two types of input files are required to for this to work.  The first is a list of the
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

#######
## ubiquitous options parsing and logger creation
my %options = ();
my $results = GetOptions (\%options, 
                            'map_file|m=s',
                            'map_dir|d=s',
                            'list_file|i=s',
                            'list_file_glob|g=s',
                            'output_list|u=s',
                            'output_dir|o=s',
                            'output_subdir_size|z=s',
                            'output_subdir_prefix|x=s',
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

## load the analysis info from the map(s)
my %analysis;
my %sequence_map;
my ($map, $mtwig);

for $map (@maps) {
    $mtwig = XML::Twig->new(
                            twig_roots => {
                                            'Sequence' => \&process_map_sequence,
                                            #'Analysis' => \&process_map_analysis
                                          }
                           );

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
    ## are we filtering entries in the list file?
    if (defined $options{list_file_glob}) {
        next unless (/$options{list_file_glob}/);
    }
    
    push @bsml_files, $_;
}


## open an output list file if the user requested it
my $olfh;
if ($options{output_list}) {
    open($olfh, ">$options{output_list}") || $logger->logdie("can't create $options{output_list} : $!");
}

## these are used if we are grouping output into subdirectories
my $sub_dir   = 1;
my $files_in_dir = 0;
my $output_dir = $options{output_dir};

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

for my $bf (@bsml_files) {
    ## get the filename
    my $fname = basename($bf);
    $root_seq_stub = '';
    
    ## it is assumed that files are named like featid.analysiscomponent.bsml (eg. cpa1.assem.2.1.aat_aa.bsml)
    ##  we can pull the source ID from this file then using a regex (here, cpa1.assem.2.1)
    $fname =~ /(.+)\..+\.bsml/ || $logger->logdie("$fname does not match naming convention.  can't extract id");
    $feat_id = $1;

    $logger->debug("using $feat_id as feat_id in file $bf") if ($logger->is_debug);

    ## get the adjustment
    $adjustment = $sequence_map{$feat_id}{offset};

    $logger->debug("setting adjustment as $adjustment for $feat_id in file $bf") if ($logger->is_debug);
    
    ## are we grouping the output files?
    if ($options{output_subdir_size}) {
        $output_dir = "$options{output_dir}/$options{output_subdir_prefix}$sub_dir";
        
        ## if the output directory doesn't exist, create it.
        mkdir($output_dir) unless (-e $output_dir);
        
        ## increment the sub_dir label if we've hit our sequence limit
        $files_in_dir++;
        $sub_dir++ if ( $files_in_dir == $options{output_subdir_size} );
    }
    
    ## open the input and output files.  how we do this depends on whether the input was zipped or not
    my $ifh;
    if ($fname =~ /\.(gz|gzip)$/) {
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
                                                             'Interval-loc'       => \&processIntervalLoc,
                                                             'Site-loc'           => \&processSiteLoc,
                                                             'Analyses'           => \&processAnalyses,
                                                             'Sequence'           => \&processSequence,
                                                           },
                               twig_print_outside_roots => $ofh,
                               pretty_print => 'indented',
                             );
    
    ## do the parse
    #$twig->parsefile($bf);
    $twig->parse($ifh);
    
    ## error if we didn't find an analysis
    if ($analyses_found) {
        $analyses_found = 0;
    } else {
        $logger->error("Analysis element not found in $bf") if ($logger->is_error);
    }
    
    ## mv the temp file over the target (can't read and write to same file)
    system("mv $output_dir/$fname.part $output_dir/$fname");
    
    ## write to the list file, if requested
    if ($options{output_list}) {
        print $olfh "$output_dir/$fname\n";
    }
}


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

    ## set some defaults
    $options{filter_ends} = 0 unless ($options{filter_ends});
    $options{output_list} = 0 unless($options{output_list});
    $options{output_subdir_size}   = 0  unless ($options{output_subdir_size});
    $options{output_subdir_prefix} = '' unless ($options{output_subdir_prefix});

    if(0){
        pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});    
    }
}

sub process_map_analysis {
    my ($mtwig, $Analysis) = @_;
    
    for my $Attribute ( $Analysis->children('Attribute') ) {
        ## save each of these
        $analysis{$map}{$Attribute->{att}->{name}} = $Attribute->{att}->{content};
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

    ## there should only be one Numbering, get the refnum (offset) from it
    my $offset = $Sequence->first_child('Numbering')->{att}->{refnum};
    if (! defined $offset) { $logger->logdie("sequence $id missing Numbering->refnum") }
    $sequence_map{$id}{offset} = $offset;
}


sub processAlignedSequence { 
    my ($twig, $element) = @_;
    
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

sub processIntervalLoc { 
    my ($twig, $element) = @_;
    
    ## need to replace any startpos or endpos attributes
    for my $attribute ( qw(startpos endpos) ) {
        if (defined $element->{att}->{$attribute}) {
            $element->{att}->{$attribute} += $adjustment;
        }
    }
    
    ## don't print if we are within a Feature.  Since Features
    #   are within a sequence, we'll just print twice.
    if (($twig->context)[-1] eq 'Feature') {
        $logger->debug("not printing Interval-loc since it is within a Feature") if ($logger->is_debug);
    } else {
        $element->print($ofh);
    }
}

sub processSeqPairAlignment {
    my ($twig, $spa) = @_;
    ## spa = SeqPairAligment
    
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
    
    ## reserve next chain number for this Seq-pair-alignment
    #   this may only be used by some components, such as AAT
    my $chain_num = ++$next_chain_num{$root_seq_stub};
    $logger->debug("chain_num $chain_num reserved for root_seq_stub $root_seq_stub (feat_id $feat_id)") if ($logger->is_debug);
    
    ##  then check each child Seq-pair-run (spr) for refpos
    ##  (comppos refers to the subject sequence and is skipped)
    for my $spr ( $spa->children('Seq-pair-run') ) {
    
        if (defined $spr->{att}->{refpos}) {
            $spr->{att}->{refpos} += $adjustment;
        }
        
        ## check for a chain_number attribute
        for my $att ( $spr->children('Attribute') ) {
            if ($att->{att}->{name} && $att->{att}->{name} eq 'chain_number') {
                $logger->debug("reassigned chain_number from " . $att->{att}->{content} . " to $chain_num") if ($logger->is_debug);
                $att->{att}->{content} = $chain_num;
            }
        }
    }
    
    $spa->print($ofh);
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

sub processSiteLoc { 
    my ($twig, $element) = @_;
    
    ## need to replace any sitepos
    if (defined $element->{att}->{sitepos}) {
        $element->{att}->{sitepos} += $adjustment;
    }
    
    $element->print($ofh);
}

