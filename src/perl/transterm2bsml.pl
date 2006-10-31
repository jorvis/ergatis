#!/usr/local/bin/perl

BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";


=head1 NAME

transterm2bsml.pl - turns tranterm raw output to bsml

=head1 SYNOPSIS

USAGE: transterm2bsml.pl 
            --input_list=/path/to/some/transterm.raw.list
            --input_file=/path/to/some/transterm.raw
            --output=/path/to/transterm.bsml
            --fasta_input=/path/to/transter.input.fsa
            --compress_bsml_output=1
            --id_repository=/path/to/id_repository/
            --project=aa1
          [ --log=/path/to/file.log
            --debug=4
            --help
          ]

=head1 OPTIONS
B<--input_list,-i>
    List of input transterm raw files.

B<--input_file,-f>
    Input raw transterm file.

B<--output,-o>
    Output bsml file. If compress_bsml_output is a non zero value, a .gz will be appended to the end
    of the filename if it does not already exist.

B<--fasta_input,-a>
    The input fasta file used as input to transterm.

B<--compress_bsml_output,-c>
    A non zero value will produce gzipped output.

B<--id_repository,-d>
    A valid id repository.  See Workflow::IdGenerator for details.

=head1  DESCRIPTION

    Creates a bsml document from raw transterm output.  See format for input and of output below.

=head1  INPUT

    Raw output from transterm run.  Looks something like this:

    ######################################################################
    TransTerm v2.01 BETA (Aug  3 2006)
    
    #HEADER, yada yada

    SEQUENCE gms.contig.1 gms:3447|cmr:632 chromosome 1 Mycobacterium smegmatis MC2 {Mycobacterium smegmatis MC2} (length 6988209)

    prok.gene.66911.1      414 - 1976     + | 
    prok.gene.66913.1     4509 - 6206     + | 
    prok.gene.66912.1     6177 - 6617     + | 
    prok.gene.66914.1     6647 - 9175     + | 
    prok.gene.66915.1    10065 - 10184    + | 
    prok.gene.66916.1    10368 - 12041    + | 
    ...
    prok.gene.69489.1  6969437 - 6971797  + | 
    prok.gene.69490.1  6971794 - 6975444  + | 
    prok.gene.69491.1  6981746 - 6982972  + | 
    prok.gene.69492.1  6984030 - 6985334  + | 
    prok.gene.69493.1  6986276 - 6986407  + | 
    prok.gene.69494.1  6986599 - 6988113  + | 
    #####################################################################

=head1  OUTPUT

    Looks something like this:

    ######################################################################
    <?xml version="1.0"?>

    <Bsml>
      <Definitions>
        <Sequences>
          <Sequence class="contig" title="gms.contig.1" id="gms.contig.1" molecule="dna">
            <Feature-tables>
              <Feature-table id="Bsml1">
                <Feature class="terminator" title="prok.terminator.1.1" id="prok.terminator.1.1">
                  <Attribute name="confidence_score" content="99"></Attribute>
                  <Attribute name="hairpin_score" content="-27.6"></Attribute>
                  <Attribute name="tail_score" content="-4.11692"></Attribute>
                  <Interval-loc complement="0" startpos="6250506" endpos="6250544"></Interval-loc>
                  <Link rel="analysis" href="#transterm_analysis" role="computed_by"></Link>
                </Feature>

                ...

              </Feature-table>
            </Feature-tables>
            <Seq-data-import source="/usr/local/annotation/scratch/kgalens/PROK/output_repository/transterm/5139_default/0/gms.contig.1.transterm.raw" identifier="gms.contig.1" format="fasta" id="Bsml0"></Seq-data-import>
            <Link rel="analysis" href="#transterm_analysis" role="input_of"></Link>
          </Sequence>
         </Sequences>
      </Definitions>
      <Research>
        <Analyses>
          <Analysis id="transterm_analysis">

           ...

          </Analysis>
        </Analyses>
      </Research>
    </Bsml>
    ##############################################################################################


=head1  CONTACT

    Kevin Galens
    kgalens@tigr.org

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use BSML::BsmlBuilder;
use Workflow::Logger;
use Workflow::IdGenerator;
use Data::Dumper;

####### GLOBALS AND CONSTANTS ###########
my @inputFiles;
my $project;
my $output;
my $debug;
my $idMaker;
my $inputFasta;
my $compressOutput;
my $defline;
########################################

my %options = ();
my $results = GetOptions (\%options, 
                          'input_list|i=s',
                          'input_file|f=s',
                          'fasta_input|a=s',
                          'compress_bsml_output|c=s',
                          'output|o=s',
                          'project|p=s',
                          'id_repository|d=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || &_pod;

#Setup the logger
my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

# Check the options.
&check_parameters(\%options);

foreach my $file (@inputFiles) {
    my $tmp = &parseTranstermData($file);
    my $data;
    $data->{$file} = $tmp;
    my $bsml = &generateBsml($data);
    $bsml->write($output,'',$compressOutput);
}

exit(0);

######################## SUB ROUTINES #######################################
sub parseTranstermData {
    my $file = shift;
    my $retval;
    my $seq = "";

    
    my $str = "< $file";
    if(-e $file && $file =~ /\.gz$/) {
        $str = "<:gzip $file";
    } elsif(! -e $file && -e $file.".gz") {
        $str = "<:gzip $file";
    }
    open(IN, $str) ||
        &_die("Unable to open the transterm raw file $file ($!)");
    
    while(<IN>) {
        if(/^SEQUENCE\s(.*?)\s.*length\s(\d+)/) {
            $seq = $1;
            $retval->{$seq}->{'length'} = $2;
        } elsif(/^\s+TERM/) {
            my @tmp = split(/\s+/);

            my $strand = 0; #0 means forward strand.

            my ($start, $stop) = ($tmp[3], $tmp[5]);
            if($start > $stop) {
                my $tmpBound = $start;
                $start = $stop;
                $stop = $tmpBound;
                $strand = 1;
            }

            $retval->{$seq}->{$tmp[2]}->{'start'} = $start;
            $retval->{$seq}->{$tmp[2]}->{'stop'} = $stop;
            $retval->{$seq}->{$tmp[2]}->{'strand'} = $strand;
            $retval->{$seq}->{$tmp[2]}->{'conf'} = $tmp[8];
            $retval->{$seq}->{$tmp[2]}->{'hp'} = $tmp[9];
            $retval->{$seq}->{$tmp[2]}->{'tail'} = $tmp[10];

        }

    }

    close(IN);
    
    return $retval;

}

sub generateBsml {
    my $data = shift;
    my $doc = new BSML::BsmlBuilder();

    foreach my $file(keys %{$data}) {
        foreach my $seq(keys %{$data->{$file}}) {
    
            #Try to figure out the class
            my $class = 'sequence';
            $class = $1 if($seq =~ /\w+\.(\w+)\./);
            
            my $seqObj = $doc->createAndAddSequence( $seq, $seq, $data->{$seq}->{'length'}, 'dna', $class );
            $doc->createAndAddSeqDataImport( $seqObj, 'fasta', $inputFasta, $seq, $seq );
            $doc->createAndAddLink( $seqObj, 'analysis', '#transterm_analysis', 'input_of' );
            $doc->createAndAddBsmlAttribute( $seqObj, 'defline', $defline );
            my $featTable = $doc->createAndAddFeatureTable($seqObj);

            #Start adding features to the sequence
            my ($key, $term);
            while(($key, $term) = each(%{$data->{$file}->{$seq}})) {
                next if($key eq 'length');                
                
                my $id = $idMaker->next_id( 'type' => 'terminator', 'project' => $project );
                my $feat = $doc->createAndAddFeature( $featTable, $seq, $seq, 'terminator');
                $feat->addBsmlLink('analysis', '#transterm_analysis', 'computed_by');
                $feat->addBsmlIntervalLoc($data->{$file}->{$seq}->{$key}->{'start'},
                                          $data->{$file}->{$seq}->{$key}->{'stop'},
                                          $data->{$file}->{$seq}->{$key}->{'strand'});
                $doc->createAndAddBsmlAttribute( $feat, 'confidence_score', 
                                                 $data->{$file}->{$seq}->{$key}->{'conf'} );
                $doc->createAndAddBsmlAttribute( $feat, 'hairpin_score',
                                                 $data->{$file}->{$seq}->{$key}->{'hp'} );
                $doc->createAndAddBsmlAttribute( $feat, 'tail_score',
                                                 $data->{$file}->{$seq}->{$key}->{'tail'} );
            }
        }

    $doc->createAndAddAnalysis( 'id' => 'transterm_analysis',
                                'sourcename' => $file );
    }
    return $doc;
}

sub check_parameters {
    my $options = shift;

    my $error = "";

    &_pod if($options{'help'});

    if($options{'input_list'}) {
        $error .= "Option input_list ($options{'input_list'}) does not exist\n" unless(-e $options{'input_list'});
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

    if($options{'compress_bsml_output'}) {
        $compressOutput = 1;
    } else {
        $compressOutput = 0;
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
        $error .= "Option id_repository is required\n";
    } else {
        $idMaker = new Workflow::IdGenerator( 'id_repository' => $options{'id_repository'} );
        $idMaker->set_pool_size('terminator' => 25 );
    }

    if($options{'debug'}) {
        $debug = $options{'debug'};
    }

    if($options{'fasta_input'}) {
        $inputFasta = $options{'fasta_input'};
        unless(-e $inputFasta) {
            $error .= "Option fasta_input ($inputFasta) does not exist\n";
        } else {
            open(FSA, "< $inputFasta") 
                or &_die("Unable to open $inputFasta");
        }

        open(FSA, "< $inputFasta") or
            &_die("Unable to open $inputFasta ($!)");

        while(<FSA>) {
            next unless(/^>(.+)/);
            $defline = $1;
            last;
        }
        close(FSA);
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
