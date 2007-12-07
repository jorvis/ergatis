#!/usr/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

split_fasta.pl - split a single-sequence FASTA file into separate, optionally overlapping, portions.

=head1 SYNOPSIS

    USAGE: split_fasta.pl 
                --input_file=/path/to/some_file.fsa 
                --output_dir=/path/to/somedir
                --fragment_length=50000
              [ --overlap_length=1000 
                --file_numbering=incremental|positional
                --file_name_root=string
                --file_name_suffix=string
                --output_list=/path/to/somefile.list
                --bsml_map=/path/to/somefile.map.bsml
              ]

=head1 OPTIONS

B<--input_file,-i>
    The input multi-fasta file to split.

B<--output_dir,-o>
    The directory to which the output files will be written.

B<--fragment_length,-f>
    Given in base pairs, the input sequence will be split into portions of this size (plus optional overhanging sequence)

B<--overlap_length,-p>
    Given in base pairs, a portion of sequence of this size will be included from the next fragment.

B<--file_numbering,-n>
    Either 'incremental' or 'positional'.  See the OUTPUT section for more information.

B<--file_name_root,-r>
    Forms the root part of the file name created.  See the OUTPUT section for more information.

B<--file_name_suffix,-u>
    Forms the suffix part of each file name created.  If omitted, defaults to 'fsa'. See the OUTPUT section for more information.

B<--output_list,-s>
    Write a list file containing the paths of each of the regular output files.  This may be useful
    for later scripts that can accept a list as input.

B<--bsml_map,-b>
    Creates a bsml file containing information which maps each of the fasta files created back
    onto the source file.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

This script takes a single-sequence FASTA file and breaks it into chunks for analysis. This 
is useful for performing analyses on large input files, such as large assemblies, using 
computational tools that do not handle or scale well with large input sequences (such as AAT).

Aside from simply splitting the input sequence into user-defined chunks, each chunk can 
optionally include a defined-length overhang from the next sequence. This helps to ensure 
that no features will be missed over the break-point when performing an analysis. Gene 
finders, for example, will miss a gene call if the gene spans the region where the sequence 
was fragmented. 

=head1  INPUT

The input is defined with --input_file and should be a FASTA file with a single sequence.  File 
extensions are ignored.  For example:

    >gi53791237 Tragulus javanicus p97bcnt gene for p97Bcnt
    ACAGGAGAAGAGACTGAAGAGACACGTTCAGGAGAAGAGCAAGAGAAGCCTAAAGAAATGCAAGAAGTTA
    AACTCACCAAATCACTTGTTGAAGAAGTCAGGTAACATGACATTCACAAACTTCAAAACTAGTTCTTTAA
    AAAGGAACATCTCTCTTTTAATATGTATGCATTATTAATTTATTTACTCATTGGCGTGGAGGAGGAAATG

Whitespace is ignored within the input file.  See the OUTPUT section for more on creation of 
output files.

=head1  OUTPUT

As the input file is fragmented, multiple output files (one for each fragment) need to be written.  
How these files are named depends on the value of the optional --file_numbering and --file_name_root 
options.  If --file_name_root is not given, the first part of the file name will simply be the 
original file name.  The --file_numbering option will control the part of each output file name
after the root portion.  The default value is 'incremental'.  Examples are in order:

    original file name: somefile.fsa
    --file_numbering='incremental'
    --file_name_root was NOT PASSED
    --file_name_suffix='fsa'
    --fragment_length=50000

    somefile.fsa.1.fsa
    somefile.fsa.2.fsa
    somefile.fsa.3.fsa
    ...

But that's probably not what you want, since .fsa is now in the file name twice.  It would be better
to do it like this:

    original file name: somefile.fsa
    --file_numbering='incremental'
    --file_name_root='somefile'
    --file_name_suffix='fsa'
    --fragment_length=50000

    somefile.1.fsa
    somefile.2.fsa
    somefile.3.fsa
    ...

This may be fine for some applications, but these file names give no positional information
about how the contents of the file relate to the original, larger file.  If this is needed,
use --file_numbering=positional .  More examples:

If, for example, your large file is called 'somefile.fsa' and
your --fragment_length is set to 50000, the files generated will be named like:

    original file name: somefile.fsa
    --file_numbering='incremental'
    --file_name_root='blastres'
    --file_name_suffix='fsa'
    --fragment_length=50000

    blast_res.0.fsa
    blast_res.50000.fsa
    blast_res.100000.fsa
    ...

The FASTA headers for each of the fragment files created will have a slightly modified header.
For each, the ID of the file (such as blast_res.0 above) will be inserted as the first element
of the header, followed by a space.


=head1  CONTACT

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Basename;
use Pod::Usage;
BEGIN {
use Ergatis::Logger;
use BSML::BsmlBuilder;
}

my %options = ();
my $results = GetOptions (\%options, 
                          'input_file|i=s',
                          'fragment_length|f=s',
                          'overlap_length|p=s',
                          'file_numbering|n=s',
                          'file_name_root|r=s',
                          'file_name_suffix|u=s',
                          'output_list|s=s',
                          'bsml_map|b=s',
                          'output_dir|o=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

## open the list file for appending if one was passed
my $listfh;
if (defined $options{output_list}) {
    open($listfh, ">>$options{output_list}") || $logger->logdie("couldn't create $options{output_list} list file");
}

## open the input file
open(my $ifh, "<$options{input_file}") || $logger->logdie("the input file $options{input_file} could not be read");

## parse until you find the header
my $header;
while (<$ifh>) {
    if (/^\>(.*)/) {
        $header = $1;
        last;
    }
}

## make sure we found one
unless (defined $header) { $logger->logdie("couldn't parse a header from the input sequence.  cowardly refusing to go on.") }

## get the name portion of the input file (without the preceeding path)
my $fname = fileparse($options{input_file});

## start a BSML mapping document if requested
my ($doc, $refseq, $analysis);
if (defined $options{bsml_map}) {
    $doc = new BSML::BsmlBuilder();
    
    ## add the reference sequence.  do the root name first, if passed, else
    ##  just use the fname.
    $refseq = $doc->createAndAddSequence( ($options{file_name_root} || $fname), undef, undef, undef, 'assembly' );
    $refseq->addBsmlLink('analysis', '#split_fasta_analysis');
    $refseq->addBsmlAttr('sourceuri', $options{input_file});
    
    ## add this analysis
    $analysis = $doc->createAndAddAnalysis(
                                            id => 'split_fasta_analysis',
                                            name => 'split_fasta',
                                            program => $0,
                                            sourceuri => $options{input_file},
                                          );
    $analysis->addBsmlAttr('fragment_length', $options{fragment_length});
    $analysis->addBsmlAttr('overlap_length', $options{overlap_length});
    $analysis->addBsmlAttr('input_file', $options{input_file});
    $analysis->addBsmlAttr('file_numbering', ($options{file_numbering} || ''));
    $analysis->addBsmlAttr('file_name_root', ($options{file_name_root} || ''));
    $analysis->addBsmlAttr('output_list', ($options{output_list} || 'X'));
    $analysis->addBsmlAttr('bsml_map', $options{bsml_map});
    $analysis->addBsmlAttr('output_dir', $options{output_dir});
}

my $flength = $options{fragment_length};
my $olength = $options{overlap_length};
my $tlength = $flength + $olength;
my $chain = '';
my $seqnum = 0;
my $offset = 0;

## we're going to build a sequence chain, adding the sequences on each line
##  as we get to it.  Once it is large enough (fragment + overlap), write out
##  that sequence and continue.
while (<$ifh>) {
    ## die if you hit another header
    if (/\>/) { $logger->logdie("found more than one header in input file!  single-sequence fasta files only") }

    ## remove whitespace from this line
    s/\s//g;

    ## add this sequence to the chain
    $chain .= $_;
    
    ## if the individual sequence lines were especially long, each line could be longer
    ##  than our fragment_length.  so we need a loop here
    for (;;) {
        if (length $chain > $tlength) {
            &writeSequence($seqnum, \$header, $tlength, \substr($chain, 0, $tlength));

            ## purge the part we just wrote (- the overlap)
            $chain = substr($chain, $flength);
            
            $offset += $flength;
            
            if ($options{file_numbering} eq 'positional') {
                $seqnum = $offset;
            } else {
                ++$seqnum;
            }
        } else {
            last;
        }     
    }
}

## don't forget to purge the chain
&writeSequence($seqnum, \$header, length($chain), \$chain);

## write the bsml map doc
if ($doc) {
    $doc->write($options{bsml_map});
}

exit;

sub check_parameters {
    my $options = shift;
    
    ## make sure input_file and output_dir were passed
    unless ( $options{input_file} && $options{output_dir} && $options{fragment_length} ) {
        #$logger->logdie("Required options are:\n\t--input_file\n\t--output_dir\n\t--fragment_");
        pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});    
    }
    
    ## make sure input_file exists
    if (! -e "$options{input_file}") {
        $logger->logdie("the input file passed ($options{input_file}) cannot be read or does not exist");
    }
    
    ## make sure the output_dir exists
    if (! -e "$options{output_dir}") {
        $logger->logdie("the output directory passed could not be read or does not exist");
    }
    
    ## set any defaults
    $options{overlap_length} = 0 if (! defined $options{overlap_length});
    $options{file_numbering} = 'incremental' if (! defined $options{overlap_length});
    
    ## if the file name suffix was defined, keep it, but take off the dot if passed.
    if (defined $options{file_name_suffix}) {
        if ( $options{file_name_suffix} =~ /^\.(.+)/ ) {
            $options{file_name_suffix} = $1;
        }
    ## else assign the default
    } else {
        $options{file_name_suffix} = 'fsa';
    }
}

sub writeSequence {
    my ($id, $header, $len, $seq) = @_;
    
    ## what's the root of the filename
    my $froot;
    if (defined $options{file_name_root} && $options{file_name_root} ne '') {
        $froot = $options{file_name_root};
    } else {
        $froot = $fname;
    }
    
    my $filepath = "$options{output_dir}/$froot.$id.$options{file_name_suffix}";
    
    ## write the sequence
    $logger->debug("Writing sequence to $filepath");
    open (my $ofh, ">$filepath") || $logger->logdie("can't create $filepath");
    print $ofh ">$froot.$id $$header\n$$seq\n";
    
    ## should we add to the list file?
    if (defined $listfh) {
        print $listfh "$filepath\n";
    }
    
    ## should we add to a bsml map doc?
    if ($doc) {
        ## add the sequence
        my $sequence = $doc->createAndAddSequence( "$froot.$id", undef, $len, undef, 'tiling_path_fragment' );
        
        ## add a ref to the file
        $sequence->addBsmlAttr('sourceuri', $filepath);
        
        ## add a link to the analysis
        $sequence->addBsmlLink('analysis', '#split_fasta_analysis');
        
        ## add a Numbering element that maps this on the original
        $doc->createAndAddNumbering( 
                                        seq => $sequence, 
                                        seqref => ($options{file_name_root} || $fname),
                                        refnum => $offset,
                                        ascending => 1,
                                   );
    }
    
}
