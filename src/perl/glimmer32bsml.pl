#!/usr/local/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";


=head1 NAME

glimmer32bsml.pl - Creates a bsml document from glimmer3 raw
    output

=head1 SYNOPSIS

USAGE: glimmer32bsml.pl
            --input_list=/path/to/some/glimmer3.raw.list
            --input_file=/path/to/some/glimmer3.raw
            --output=/path/to/transterm.bsml
            --project=aa1
            --id_repository=/path/to/id_repository
            --fasta_input=/path/to/glimmer3/input.fsa
          [ --log=/path/to/file.log
            --debug=4
            --help
          ]

=head1 OPTIONS

B<--input_list,-i>
    Input list of glimmer3 raw output files (.predict)

B<--input_file,-f>
    Input glimmer3 raw file (.predict)

B<--output,-o>
    The output bsml file.

B<--project,-p>
    The project (used for id generation).

B<--id_repository,-r>
    Id repository for use by Workflow::IdGenerator.pm

B<--fasta_input,-a>
    The input file that was used as input for the glimmer3 run

B<--log,-l>
    Logfile.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

This script is used to convert the output from a glimmer3 search into BSML.

=head1  INPUT

You define the input file using the --input option.  This file does not need any
special file extension.  The regular output of glimmer3 looks like this:


>zma1.assembly.5808
orf00001   319296      174  +1     3.95
orf00002      430      969  +1     3.91
orf00003     1197     1051  -1    12.79
orf00004     1258     1368  +1     6.95
orf00005     1541     1681  +2     8.58
orf00007     1889     1755  -3    11.14
orf00008     2014     1874  -2     0.37
orf00009     2153     2359  +2    48.95

Where the columns from left to right contain:
orf_id   start_position   end_position   reading_frame   "raw"_score

A few things to note:
- glimmer3 will find genes spanning the end/beginning of the given sequence, as
if dealing with a circular molecule (unless given the -l or --linear option).  These
are dealt with in bsml by converting the start_pos into a coordinate to the left of
the origin (a negative number).
- The raw_score in the .predict file may be different than that found in the 
.details file, due to adjustments made for the PWM and start codon frequency.
This value from the .predict file can therefore be used to make direct comp-
arisons between predictions.
- The ranges given by glimmer3 differ from those given by glimmer2 because they now
include the stop codon.

=head1 OUTPUT


After parsing the input file, a file specified by the --output option is created.  This script
will fail if it already exists.  The file is created, and temporary IDs are created for
each result element.  They are only unique to the document, and will need to be replaced
before any database insertion.

Base positions from the input file are renumbered so that positions start at zero.

=head1  CONTACT

    Kevin Galens
    kgalens@tigr.org

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Gene;
use BSML::GenePredictionBsml;
use Workflow::IdGenerator;
use Ergatis;:Logger;

####### GLOBALS AND CONSTANTS ###########
my @inputFiles;               #Holds input files
my $project;                  #The project (ex aa1)
my $output;                   #Output file
my $idMaker;                  #The Workflow::IdGenerator
my $bsml;                     #BSML::BsmlBuilder object object.
my $data;                     #Holds parsed glimmer3 information
my $inputFsa;                 #The fasta file input to glimmer3
my $debug;                    #The debug variable
########################################

my %options = ();
my $results = GetOptions (\%options, 
                          'input_list|i=s',
                          'input_file|f=s',
                          'output|o=s',
                          'project|p=s',
                          'id_repository|r=s',
                          'fasta_input|a=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || &_pod;

#Setup the logger
my $logfile = $options{'log'} || Ergatis;:Logger::get_default_logfilename();
my $logger = new Ergatis;:Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

# Check the options.
&check_parameters(\%options);

foreach my $file (@inputFiles) {
    $data = &parseGlimmer3Data($file);
    $bsml = &generateBsml($data);
}

$bsml->writeBsml($output);

exit(0);

######################## SUB ROUTINES #######################################
sub parseGlimmer3Data {
    my ($file) = @_;
    my $genes;
    my $foundId = 0;

    open(IN, "<$file") or &_die("Unable to open $file");
    
    while(<IN>) {
        if(/^>(.*?)\s/) {

            $foundId = $1;

        } elsif($foundId) {

            my @cols = split(/\s+/,$_);
            ($cols[1],$cols[2]) = ($cols[2],$cols[1]) if($cols[1] > $cols[2]);
                
             
            #Create some genes and push them ontot he $genes array
            my $tmp = new Gene( $idMaker->next_id( 'type' => 'gene',
                                                   'project' => $project),
                                $cols[1]-1, $cols[2]-1, ($cols[3] > 0) ? 0 : 1,
                                $foundId);

            foreach my $type(qw(exon CDS transcript polypeptide)) {
                my $typeid =$idMaker->next_id( 'type' => $type,
                                               'project' => $project);
                $tmp->addFeature($typeid, $cols[1]-1, $cols[2]-1, ($cols[3] > 0) ? 0 : 1,
                                $type);
            }

            my $count = $tmp->addToGroup($tmp->getId, { 'all' => 1 });
            &_die("Nothing was added to group") unless($count);

            push(@{$genes}, $tmp);

        } else {
            &_die("Didn't find the id");
        }
    }

    close(IN);

    return $genes;

}

sub generateBsml {
    my $data = shift;

    #Create the document
    my $doc = new GenePredictionBsml( 'glimmer3', $inputFsa );
    
    foreach my $gene(@{$data}) {
        $doc->addGene($gene);
    }
    
    my $seqId;
    open(IN, "< $inputFsa") or &_die("Unable to open $inputFsa");
    while(<IN>) {
        #assume it's a single fasta file
        if(/^>([^\s+]+)/) {
            $seqId = $1;
            last;
        }
    }
    close(IN);

    my $addedTo = $doc->setFasta($seqId, $inputFsa);
    &_die("$seqId was not a sequence associated with the gene") unless($addedTo);

    return $doc;
    
    
}

sub check_parameters {
    my $options = shift;

    my $error = "";

    &_pod if($options{'help'});

    if($options{'input_list'}) {
        $error .= "Option input_list ($options{'input_list'}) does not exist\n" 
            unless(-e $options{'input_list'});
        open(IN, "< $options{'input_list'}") || &_die("Unable to open $options{'input_list'} ($!)");
        @inputFiles = <IN>;
        close(IN);
    }
    if($options{'input_file'}) {
        $error .= "Option input_file ($options{'input_file'}) does not exist\n" unless(-e $options{'input_file'});
        push(@inputFiles, $options{'input_file'});
    }
    unless($options{'input_list'} || $options{'input_file'}) {
        $error .= "Either option input_list or input_file is required\n";
    }

    unless($options{'output'}) {
        $error .= "Option output is required\n";
    } else {
        $output = $options{'output'};
    }

    unless($options{'project'}) {
        $error .= "Option project is required\n";
    } else {
        $project = $options{'project'};
    }

    unless($options{'id_repository'}) {
        $error .= "Option id_repository is required.  Please see Workflow::IdGenerator ".
            "for details.\n";
    } else {
        $idMaker = new Workflow::IdGenerator( 'id_repository' => $options{'id_repository'} );
        $idMaker->set_pool_size( 'exon'        => 20,
                                 'transcript'  => 20,
                                 'gene'        => 20,
                                 'polypeptide' => 20,
                                 'CDS'         => 20 );
                                 
    }

    unless($options{'fasta_input'}) {
        $error .= "Option fasta_input is required\n";
    } else {
        $error .= "$options{'fasta_input'} (fasta_input) does not exist\n" 
            unless(-e $options{'fasta_input'});
        $inputFsa = $options{'fasta_input'};
    }
    
    if($options{'debug'}) {
        $debug = $options{'debug'};
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
