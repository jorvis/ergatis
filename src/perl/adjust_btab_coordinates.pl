#!/usr/local/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

adjust_btab_coordinates.pl - adjust the positional coordinates within a nap
output alignment file

=head1 SYNOPSIS

USAGE:  adjust_btab_coordinates.pl
            --map_file=/path/to/somemapfile.bsml
            --map_dir=/or/path/to/somedir
            --list_file=/path/to/somefile.list
            --output_dir=/path/to/somedir
          [ --list_file_glob='.*.nap.btab'
            --output_list=/path/to/some.list
            --output_subdir_size=1000
            --output_subdir_prefix=fasta
            --debug=4
            --log=/path/to/somefile.log
          ]

=head1 OPTIONS

B<--map_file,-m> 
    Input BSML map file.

B<--map_dir,-d> 
    Directory of input BSML map files.

B<--list_file,-i> 
    Path to a list of btab files.

B<--list_file_glob,-g> 
    Use to filter the btab files considered in --list_file.  This should probably
    be given as in the example above.

B<--output_list,-u>
    Optional.  If passed, will create an output list with the full paths to each of the 
    BSML files created by this script.

B<--output_dir,-o> 
    Directory where output btab files will be written.

B<--output_subdir_size,-z>
    If defined, this script will create numbered subdirectories in the output directory, each
    containing this many sequences files.  Once this limit is reached, another subdirectory
    is created.

B<--output_subdir_prefix,-x>
    To be used along with --output_subdir_size, this allows more control of the names of the
    subdirectories created.  Rather than just incrementing numbers (like 10), each subdirectory 
    will be named with this prefix (like prefix10).

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

The purpose of this script is to read through the results of such an analysis and 
adjust the coordinates to reflect those in the original input file.  This is done using
a set of btab result files along with a mapping file.  Read the INPUT section for more
information.

To use this component, you would most likely have used a splitting component previously
in the pipeline, such as split_fasta, which generates the input maps.

=head1 INPUT

Two types of input files are required to for this to work.  The first is a list of the
btab result files from an analysis, defined using --list_file.  The list file should
look something like this:

    /some/path/cpa1.assem.1.0.nap.btab
    /some/path/cpa1.assem.1.1.nap.btab
    /some/path/cpa1.assem.1.10.nap.btab
    /some/path/cpa1.assem.1.11.nap.btab

In this case, a nap analysis was run on 4 sequences, probably named like cpa1.assem.1.0

If the list file contains more lines than you actually want to consider, you can filter which
will be included by using the --list_file_glob option.

The BSML mapping file, passed with the --map_file option, gives information about how the
coordinates in these files relate to the original non-split file, here named cpa1.assem.1

If mapping information for the files in your list come from multiple BSML maps, you can use
the --map_dir option.  This will load all maps in the specified directory whose names end
in '.map.bsml'

It is important to note that currently the files must match the naming convention like the
examples above.  This convention is:

    featid.component.btab

where featid must match the id of some Sequence element in the mapping file.  The "component"
portion of the name is currently ignored by this script.  This is how the relationship 
between the fragment and the map is defined.

Also note this script handles compressed input.  If the input BSML files have been gzipped
and end in either the '.gz' or '.gzip' extension, they will be decompressed, modified and
recompressed on the fly.

=head1 OUTPUT

The output alignment files will be written to the directory specified by --output_dir.  You may
specify the output_dir to be the same location as the original files to overwrite them.

As each alignment file is analyzed, a few changes are made to reflect the coordinates back onto
the non-fragmented sequence.  Specifically, the 7th and 8th columns of the tab-delimited
file are adjusted.  These define the query start and stop coordinates, respectively.

Note that the file names, such as cpa1.assem.1.0.nap.btab, are not changed, though
they no longer contain information directly concerning the fragment.  This way these
sequences can be loaded as they are - else we would have to concatenate each of the
fragmented files into one large file before loading, which is not desirable.

The first column of the output file, however, is stripped of the fragment sequence ID,
since the coordinates no longer map to that fragmented sequence.  What remains should be
the ID of the original, larger sequence.

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
                            'output_dir|o=s',
                            'output_list|u=s',
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
                                          }
                           );

    $mtwig->parsefile($map);
    $mtwig->purge;  ## free memory
}

## we want to loop over a list of btab files, each of which should have a
##  record of the same id in the mapping file.
## first get the list of files
my @btab_files;
open(my $lfh, "<$options{list_file}") || $logger->logdie("can't open list file: $!");
for (<$lfh>) {
    chomp;
    ## are we filtering entries in the list file?
    if (defined $options{list_file_glob}) {
        next unless (/$options{list_file_glob}/);
    }
    
    push @btab_files, $_;
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

## loop through each btab file and adjust coordinates
my $adjustment;
my $feat_id;

for my $bf (@btab_files) {
    ## get the filename
    my $fname = basename($bf);
    
    ## it is assumed that files are named like featid.analysiscomponent.bsml (eg. cpa1.assem.2.1.nap.raw)
    ##  we can pull the source ID from this file then using a regex (here, cpa1.assem.2.1)
    $fname =~ /(.+)\..+\.btab/ || $logger->logdie("$fname does not match naming convention.  can't extract id");
    $feat_id = $1;

    ## get the adjustment
    if (defined $sequence_map{$feat_id}{offset}) {
        $adjustment = $sequence_map{$feat_id}{offset};
    } else {
        $logger->logdie("adjustment value not found for feature $feat_id in file $bf");
    }

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
    my ($ifh, $ofh);
    if ($fname =~ /\.(gz|gzip)$/) {
        open ($ifh, "<:gzip", $bf)                       || $logger->logdie("can't read zipped input file '$bf': $!");
        open ($ofh, ">:gzip", "$output_dir/$fname.part") || $logger->logdie("can't create output file: $!");
    } else {
        open ($ifh, "<$bf")                     || $logger->logdie("can't read input file $bf: $!");
        open ($ofh, ">$output_dir/$fname.part") || $logger->logdie("can't create output file: $!");
    }

    ## do the adjustment (same as brian's previous code)
    ##  the numbers to adjust are in columns 7 and 8
    while (<$ifh>) {
        my @cols = split(/\t/);
        
        ## brian reports that the 11th column will have a negative value on rows
        ##  that don't represent real results.  the btab binary writes such a row
        ##  if there were no results.
        if ($cols[10] > 0) {

            ## remove the portion ID from the first column, leaving the original that
            ##  should be after it
            $cols[0] =~ s/^$feat_id\s*//;
            $cols[6] += $adjustment;
            $cols[7] += $adjustment;
        }

        ## error if there is nothing left in the first column
        $logger->error("id completely removed from btab row") unless ($cols[0]);
        
        ## write out this line
        print $ofh join("\t", @cols);
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

    ## some defaults
    $options{output_list} = 0 unless ($options{output_list});
    $options{output_subdir_size}   = 0  unless ($options{output_subdir_size});
    $options{output_subdir_prefix} = '' unless ($options{output_subdir_prefix});

    if(0){
        pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});    
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


