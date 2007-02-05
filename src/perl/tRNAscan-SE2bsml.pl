#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

tRNAscan-SE2bsml.pl - convert tRNAcan-SE output to BSML

=head1 SYNOPSIS

USAGE: tRNAscan-SE2bsml.pl 
        --input=/path/to/tRNAscanfile 
        --output=/path/to/output.bsml
        --fasta_input=/path/to/fastafile
        --id_repository=/path/to/id_repository
        --project=aa1 
        [ --gzip_output=1       
          --log=/path/to/logfile 
          --debug=3
        ]

=head1 OPTIONS

B<--input,-i> 
    Input file file from a tRNAscan-SE search.

B<--output,-o> 
    Output BSML file (will be created, must not exist)

B<--fasta_input,-a>
    The input file that was used as input for the tRNAscan-SE run

B<--id_repository,-r>
    Path to the project's id repository

B<--gzip_output,-g>
    Optional. A non-zero value will make compressed output.
    
B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--project,-p> 
    [DEPRECATED] Project ID.  Used in creating feature ids.  Defaults to 'unknown' if
    not passed.

    *Now parses project from input sequence id parsed from tRNAscan-SE raw file.
    *ex:  project.type.num     (adg.assembly.1)

B<--log,-l> 
    Log file

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used to convert the output from a tRNAscan-SE search into BSML.

=head1 INPUT

tRNAscan-SE can be run using multiple input sequences simultaneously, and this
script supports parsing single or multiple input result sets.  Usual tRNAscan-SE
output looks like:

    Sequence                tRNA            Bounds          tRNA    Anti    Intron Bounds   Cove
    Name            tRNA #  Begin           End             Type    Codon   Begin   End     Score
    --------        ------  ----            ------          ----    -----   -----   ----    ------
    51595           1       101064          101137          Asn     GTT     0       0       82.29
    51595           2       705796          705868          Ala     AGC     0       0       68.12
    ...
    51595           13      3488675         3488743         Pseudo  GTA     0       0       22.49
    51595           14      3493468         3493555         Ser     GCT     3493506 3493513 40.19

You can elimate the headers in the original tRNAscan-SE output file by running 
tRNAscan-SE using the -b option.  If they are present, they should be ignored by 
this script.

You define the input file using the --input option.  This file does not need any
special file extension.

=head1 OUTPUT

After parsing the input file, a file specified by the --output option is created.  This script
will fail if it already exists.  The file is created, and temporary IDs are created for
each result element.  They are only unique to the document, and will need to be replaced
before any database insertion.

Base positions from the input file are renumbered so that positions start at zero.  The
current output elements from tRNAscan-SE that are not represented in the BSML file are:

    tRNA type
    Anti Codon
    Cove Score

These need to be included later.

=head1 CONTACT

    Jason Inman
    jinman@tigr.org

=cut

use warnings;
use strict;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Ergatis::Logger;
use Ergatis::IdGenerator;
use Chado::Gene;
use BSML::GenePredictionBsml;
use BSML::BsmlBuilder;


my $input;      # Input file name
my $output;     # Output file name
my $fasta_input;   # fasta input to tRNAscan-SE
my $project;    # project id
my $idcreator;  # Ergatis::IdGenerator object
my $bsml;       # BSML::BsmlBuilder object
my $data = [];  # parsed tRNAscan-SE data
my $debug = 4;  # debug value.  defaults to 4 (info)

my %options = ();
my $results = GetOptions (\%options, 
                          'input|i=s',
                          'output|o=s',
                          'fasta_input|a=s',
                          'project|p=s',
                          'log|l=s',
                          'id_repository|r=s',
                          'gzip_output|g=s',
                          'debug=s',
              'help|h') || &_pod;


# Set up the logger
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

## make sure all passed options are peachy
&check_parameters(\%options);

# Use the Gene.pm module if there exists a non-empty input
# otherwise just dump an 'empty' bsml file.
$data = &parse_tRNAscanSE_input( $input );
$bsml = &generateBsml( $data );
$bsml->writeBsml( $output, '', $options{'gzip_output'});
exit (0);

################ SUBS ######################

sub parse_tRNAscanSE_input {
# Open the input file and parse the data into an array of gene models.
    my $infile = shift;
    my $genes;

    ## open the input file for parsing
    open (IN,"< $infile") || $logger->logdie("can't open input file for reading");

    ## Strip out uninteresting lines of data from the file
    my %rawdata;
    while (<IN>) {

        my @cols = split;
 
        #check whitespace, no warn
        next if ( /^\s*$/ );
    
        ## make sure we don't parse the tRNAscan-SE output header lines
        next if ( /^sequence.*bounds.*cove/i ||
                  /^name.*end.*score/i ||
                  /\-\-.*\-\-/);

        ## there should be 9 elements in cols, unless we have an unrecognized format.
        unless (scalar @cols == 9) {
            $logger->error("the following tRNAscan-SE line was not recognized and could not be parsed:\n$_\n") 
                if ($logger->is_error);
            next;
        }
    
        ## add this data row to this sequence
        push( @{$rawdata{shift @cols}}, \@cols );

    }

    close IN;

    ## loop through each of the matches that we found and make a entries in the genes array.
    for my $seqid (keys %rawdata) {

        #Parse the project from the sequence id
        unless($project) {
            $project = $1 if($seqid =~ m|^([^/\.])\.|);
            $logger->logdie("Could not parse project from sequence id $seqid")
                unless($project);
        }

        ## loop through each array reference of this key, adding to the data array as necessary.
        foreach my $arr ( @{$rawdata{$seqid}} ) {
        
            ## 1 is subtracted from each position to give interbase numbering
            $$arr[1]--;     ## tRNA begin
            $$arr[2]--;     ## Bounds End

            ## Determine strandedness
            my $complement = ($$arr[1] > $$arr[2]) ? 1 : 0;

            ## First, create the gene model object
            my $currGene = new Chado::Gene ( $idcreator->next_id( 'type' => 'gene',
                                                           'project' => $project ),
                                      ($complement) ? $$arr[2] : $$arr[1],
                                      ($complement) ? $$arr[1] : $$arr[2],
                                      $complement,
                                      $seqid
                                    );

            ## Next, with the same coords, create the tRNA feature
            $currGene->addFeature( $idcreator->next_id ( 'type' => 'tRNA',
                                                         'project' => $project ),
                                   ($complement) ? $$arr[2] : $$arr[1],
                                   ($complement) ? $$arr[1] : $$arr[2],
                                   $complement,
                                   'tRNA'
                                 );

            ## an exon needs to be added.
            ##  the following will eval as true if tRNAscan-SE reported an
            ## intron in the sequence, else the entire range is added as an exon.
            ##this seems to limit tRNAscan-SE to only report tRNAs with 0 or 1 intron.
            if ($$arr[5] && $$arr[6]) {

                ## 1 is subtracted from each position to give interbase numbering
                $$arr[5]--;     ## Intron Begin
                $$arr[6]--;     ## Bounds End
            
                ## exon1
                &add_exon_and_cds($currGene, $$arr[1], $$arr[5], $complement);

                ## exon2
                &add_exon_and_cds($currGene, $$arr[6], $$arr[2], $complement);
            
            } else {
                ## just add the whole thing as an exon
                &add_exon_and_cds($currGene, $$arr[1], $$arr[2], $complement);
            }

            # Handle Group now:
            my $count = $currGene->addToGroup( $currGene->getId, { 'all' => 1} );
            &_die("Nothing was added to group") unless ($count);

            push (@{$genes}, $currGene);

        }
    }

    return $genes;

}

sub generateBsml {

    my $data = shift;

    #Create the document
    my $gene_pred = new BSML::GenePredictionBsml( 'tRNAscan-SE' );
    my $doc = $gene_pred->{'doc'};

    unless($data && @{$data} > 0) {
        $gene_pred->addSequence( '', $fasta_input );
    }

    foreach my $gene( @{$data} ) {
        $gene_pred->addGene($gene);
        my $addedTo = $gene_pred->addSequence($gene->{'seq'}, $fasta_input);
        $logger->logdie("Could not find identifier ".$gene->{'seq'}." in $fasta_input")
            unless($addedTo || !(-e $fasta_input));
        $logger->logdie("fasta_input file: '$fasta_input' does not exist") unless(-e $fasta_input);
        
    }
    
    return $gene_pred;

}

sub add_exon_and_cds {
# Add an exon and cds to the gene object.  These features should share 
# start and stop coordinates.
     my ($gm, $start, $stop, $complement) = @_;

     foreach my $type (qw(exon CDS)) {
         
         # Add the exon or CDS to the gene model object
         $gm->addFeature( $idcreator->next_id( 'type' => $type,
                                               'project' => $project ),
                          ($complement) ? $stop : $start,
                          ($complement) ? $start : $stop,
                          $complement,
                          $type
                        );

     }

}

sub check_parameters {
# Parse options hash, make sure certain ones exist and have meaningful values.
# Also take the time to setup some globally used values.

    my $error = '';

    &_pod if $options{'help'}; 
 
    ## make sure input file was given and exists
    if ($options{'input'}) {

        if (! -e $options{'input'}) {
            $logger->logdie("input file $options{'input'} does not exist")
        }

        $input = $options{'input'};

    } else {
        $error .= "--input_file is a required opttion\n";
    }
    
    ## make sure output file was given and doesn't exist yet
    if ($options{'output'}) {

        $output = $options{'output'};

    } else {
        $error .= "--output_file is a required option\n";
    }

    ## Make sure we're given the input fasta
    if ($options{'fasta_input'}) {
        $fasta_input = $options{'fasta_input'};
    } else {
        $error .= "--fasta_input is a required option\n";
    } 
    
    ## Now set up the id generator stuff
    if ($options{'id_repository'}) {

        # we're going to generate ids
        $idcreator = new Ergatis::IdGenerator('id_repository' => $options{'id_repository'});

        # Set the pool size
        $idcreator->set_pool_size('gene'=>30,'tRNA'=>30,'exon'=>30,'CDS'=>30);

    } else {
        $error .= "--id_repository is a required option\n";
    }

    # The debug option is not required... but let's set it up if it's been given
    if ($options{'debug'}) {
        $debug = $options{'debug'};
    }

    &_die($error) if $error;

    return 1;

} # END OF check_parameters 

sub _pod {
# Used to display the perldoc help.
    pod2usage( {-exitval => 0, -verbose => 2} );
}

sub _die {
# Kinda like a hitman, this is called when something needs to DIE.
    my $msg = shift;
    $logger->logdie($msg);
}
