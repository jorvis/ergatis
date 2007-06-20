#!/usr/local/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";


=head1 NAME

augustus2bsml.pl - Creates a bsml document from augustus raw
    output

=head1 SYNOPSIS

USAGE: augustus2bsml.pl
            --input_file=/path/to/some/augustus.raw
            --output=/path/to/augustus.bsml
            --project=aa1
            --id_repository=/path/to/id_repository
            --fasta_input=/path/to/glimmer3/input.fsa
          [ --log=/path/to/file.log
            --debug=4
            --help
          ]

=head1 OPTIONS

B<--input_file,-i>
    Raw augustus output.

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

This script is used to convert the output from an augustus search into BSML.

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
use Chado::Gene;
use BSML::GenePredictionBsml;
use Ergatis::IdGenerator;
use Ergatis::Logger;
use Data::Dumper;

####### GLOBALS AND CONSTANTS ###########
my $inputFile;               #Holds input files
my $project;                  #The project (ex aa1)
my $output;                   #Output file
my $idMaker;                  #The Ergatis::IdGenerator
my $bsml;                     #BSML::BsmlBuilder object object.
my $data;                     #Holds parsed augustus information
my $inputFsa;                 #The fasta file input to augustus
my $sourcename;               #Used for the analysis section of the bsml
my $debug;                    #The debug variable
########################################

my %options = ();
my $results = GetOptions (\%options, 
                          'input_file|i=s',
                          'output|o=s',
                          'project|p=s',
                          'id_repository|r=s',
                          'fasta_input|a=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || &_pod;

#Setup the logger
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE' =>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

# Check the options.
&check_parameters(\%options);

$data = &parseAugustusData($inputFile);
$bsml = &generateBsml($data);

$bsml->writeBsml($output);

exit(0);

######################## SUB ROUTINES #######################################
sub parseAugustusData {
    my ($file) = @_;
    my @genes; #Array of gene objects.
    my ($tmp, $curGene, $curTransc) = ({}, 0, 0);
    
    open(IN, "< $file") or &_die("Unable to open $file (input_file) ($!)");
    
    while(<IN>) {
        next if(/^\#/);
        chomp;
        
        my @cols = split(/\t/);

        #The first column should contain the input sequence id.  We use it for
        #the project name if it was not passed in or if the keyword parse
        #was passed in.  
        if( ! $project || $project eq 'parse' ) {
            $project = $1 if( $cols[0] =~ /^([^\.]+)/ );
        }
        die("Could not parse project name from $cols[0]") unless( $project  && $project ne 'parse' );
        
        #Parse out the augustus assigned gene id.
        #Column two is the feature column.
        #The id will be something like g1 (for gene 1) or t1 (for transcript 1).
        my ($newGeneId, $newTransId); 
        if($cols[2] eq 'gene') {
            $newGeneId = $cols[8];
        } elsif($cols[2] eq 'transcript' && $cols[8] =~ /^(.+?)\.(.+)/) {
            $newGeneId = $1;
            $newTransId = $2;
        } elsif($cols[8] =~ /transcript_id\s\"(.+?)\.(.+?)\"/) {
            $newGeneId = $1;
            $newTransId = $2;
        } else {
            &_die("could not parse out fake gene id from column 8 ($cols[8]) on line $_");
        }
        
        #Make sure the start is always less than the stop.  And parse the strand out.
        my ($start, $end) = ($cols[3] > $cols[4]) ? ($cols[4], $cols[3]) : ($cols[3], $cols[4]);
        my $strand = ($cols[6] eq '+') ? 0 : 1;

        #If we have come across a new gene id (a fake one that augustus assigns), then 
        #create a new gene.
        unless($newGeneId eq $curGene) {
            $tmp->{$newGeneId} = new Chado::Gene( $idMaker->next_id( 'type'    => 'gene',
                                                              'project' => $project),
                                           $start, $end, $strand, $cols[0]);  
            $curGene = $newGeneId;
        }

        #Deal with the lines according to column 2 (the type column).
        if($cols[2] =~ /^(initial|internal|terminal|single)$/) {
            my $exonId = $idMaker->next_id( 'type' => 'exon', 'project' => $project );
            $tmp->{$newGeneId}->addExon($exonId, $start, $end, $strand);
            $tmp->{$newGeneId}->addToGroup($newTransId, { 'id' => $exonId } );

            my $score = $cols[5];
            $tmp->{$newGeneId}->addFeatureScore($exonId, 'p-value', $score) if($score);
        } elsif($cols[2] eq 'transcript') {
            foreach my $type(qw(CDS transcript polypeptide) ) {
                my $id = $idMaker->next_id( 'type' => $type, 'project' => $project );
                $tmp->{$newGeneId}->addFeature( $id, $start, $end, $strand, $type );
                $tmp->{$newGeneId}->addToGroup( $newTransId, { 'id' => $id } );
                
                my $score = $cols[5];
                $tmp->{$newGeneId}->addFeatureScore($id, 'p-value', $score) if($score);
            }
        }        
        
    }

    @genes = values %{$tmp};
    return \@genes;

}

sub generateBsml {
    my $data = shift;

    #Create the document
    my $doc = new BSML::GenePredictionBsml( 'augustus', $sourcename );
    
    #Add all the genes to the document
    foreach my $gene(@{$data}) {
        $doc->addGene($gene);
    }
    
    #Get the sequence id from the input fasta file.
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

    #Set the fasta input (Returns zero if not succesful)
    my $addedTo = $doc->setFasta($seqId, $inputFsa);
    &_die("$seqId was not a sequence associated with the gene") unless($addedTo);

    #Return the document.  
    return $doc;
    
}

sub check_parameters {
    my $options = shift;

    my $error = "";

    &_pod if($options{'help'});

    if($options{'input_file'}) {
        $error .= "Option input_file ($options{'input_file'}) does not exist\n" 
            unless(-e $options{'input_file'});
        $inputFile = $options{'input_file'};
    } else {
        $error = "Option input_file is required\n";
    }

    #Parse the sourcename from the inputfile
    if( $inputFile ) {
        $sourcename = $1 if( $inputFile =~ m|^(.*)/[^/]+$| );
    }
    
    #If the sourcename contains the iterator and group directories
    #( ie if run in ergatis ), remove them.
    $sourcename =~ s|/.*/g\d+$|| if( $sourcename );

    unless($options{'output'}) {
        $error .= "Option output is required\n";
    } else {
        $output = $options{'output'};
    }

    unless($options{'project'}) {
        $project = "parse";
    } else {
        $project = $options{'project'};
    }

    unless($options{'id_repository'}) {
        $error .= "Option id_repository is required.  Please see Ergatis::IdGenerator ".
            "for details.\n";
    } else {
        $idMaker = new Ergatis::IdGenerator( 'id_repository' => $options{'id_repository'} );
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
