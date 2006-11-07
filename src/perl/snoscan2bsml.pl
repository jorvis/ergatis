#!/usr/local/bin/perl

use lib $ENV{'PERL_MOD_DIR'};

=head1 NAME

snoscan2bsml.pl - converts snoscan output to BSML

=head1 SYNOPSIS

    USAGE: snoscan2bsml.pl 
            --input=/path/to/snoscanfile
            --output=/path/to/output.bsml
         [  --project=projectName 
            --class=class of sequence
            --degug=integer(level of debug)
            --fasta_file=input file to snoscan
         ]

=head1 OPTIONS

B<--input_file,-i>
    Input file from the snoscan search.

B<--output_file,-o>
    The file will be created

B<--class,-c>
    This is used by the workflow (for Id purporses)

B<--project,-p>
    Project ID. Used in creating feature IDs.  Defaults to unknown 
    if not passed.

B<--command_id>
    The command passed in from the workflow

B<--debug,-d>
    Debug level.  Use a large number to turn on verbose debugging.

B<--log,-l>
    Log file

B<--help,-h>
    This help message.

=head1  DESCRIPTION

snoscan2bsml.pl will take an output file from snoscan, sorted or raw,
and convert the format to BSML.  

=head1  INPUT

The input is defined with --input_file and should be an output file from the program snoscan.  For an example, see 

http://lowelab.ucsc.edu/snoscan/snoscanReadme.html

=head1  OUTPUT

After parsing the input file, a file specified by the --output option will 
be generated (this file should not exist already).  Temporary ID's will be
created for each result element and are only unique to this document.

Base positions are renumbered to start at zero and represent interbase
numbering.

=head1  CONTACT

    Kevin Galens
    kgalens@tigr.org

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Ergatis;:Logger;
use BSML::BsmlRepository;
use Papyrus::TempIdCreator;
use BSML::BsmlBuilder;
use BSML::BsmlParserTwig;

my %options = ();
my $results = GetOptions (\%options, 
                          'input|i=s',
                          'output|o=s',
                          'class|c=s',
                          'project|p=s',
                          'command_id=s', #Input from workflow
                          'log|l=s',
                          'debug|d=s',
                          'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis;:Logger::get_default_logfilename();
my $logger = new Ergatis;:Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});

$logger = $logger->get_logger();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

unless($options{'compress_bsml_output'}) {
    $options{'compress_bsml_output'} = 0;
}

#Create a new BSML document
my $doc = new BSML::BsmlBuilder();

#Create a new id creator
my $idcreator = new Papyrus::TempIdCreator();

#Open the input file (.sorted)
my $inputFile;

#The fasta file header line (Bug 3781)
my $defline;

#Should be able to handle gzip input without hassle (Bug 2591)
my $mode = "<";
$mode = "<:gzip" if($options{'input'} =~ /\.gz$/);

open (IFH, $mode, "$options{input}")
    or $logger->logdie("can't open input file for reading");
    #or die ("Unable to open output file ($!)");

#If there was no sorted file made (which is indicated by the
#words 'No output generated') then switch the sorted to raw
#and open the file again.  This means that the user did not
#choose to sort the snoscan output and therefore has no
#data in the .sorted file and the .raw file should be used 
#instead.
if(<IFH> =~ /^No output generated/) {
    $inputFile = $options{input};
    if(!($inputFile =~ s/sorted/raw/)) {
        $logger->logdie ("Input file: $inputFile\n");
    }
    close(IFH);
    print "Opening $inputFile\n";
    open(IFH, "< $inputFile") 
        or die "Unable to open $inputFile ($!)";
}

#%data will hold the summary line of output for each snoRNA predicted.
my %data;

#Holds start and end positions of snoRNAs already predicted.
#Key is the location interval and value is the sequence id
my %locs = {};    

#Go through the input file
while(<IFH>) {
    my @cols = split;
    if(/^>>/) {
        unless(@cols == 14 || @cols == 15) {
            warn ("Cannot recognize input line\n(@cols)");
            next;
        }
    
        #Unique ID = seqId (column1) and loc int (col3)
        #Value is the score of the prediction.  Only keeping
        #high score.
        if(!(defined($data{"$cols[1]".";;$cols[3]"})) 
           || $data{"$cols[1]".";;$cols[3]"} < $cols[2] ){
           $data{"$cols[1]".";;$cols[3]"} =  $cols[2];
       } 
    }
}

#Done with IFH
close(IFH);

foreach my $unique (keys %data) {
    my ($seqID, $intLoc) = split(';;', $unique);

    #Create a sequence object
    my $seq = $doc->createAndAddSequence($seqID, undef, undef, 'na',
                                         $options{'class'});

    #Link it to the analysis as the input of
    $seq->addBsmlLink('analysis', '#snoscan_analysis', 'input_of');

    #Create a feature group
    my $ft = $doc->createAndAddFeatureTable($seq);

    #If an input file was listed add the SeqDataImport
    if($options{'fasta_file'}) {
        $doc->createAndAddSeqDataImport( $seq, 'fasta', $options{fasta_file}, '', $seqID );
        $doc->createAndAddBsmlAttribute( $seq, 'defline', $defline );
    }

    #Create some variables used in following foreach loop
    my $startNuc;      #Start of predicted snoRNA (interbase)
    my $endNuc;        #End of predicted snoRNA (interbase)
    my $gene;          #Feature Object
    my $snoRNA;        #Also feature ojbect
    my $exon;          #Again, feature object
    my $fg;            #The feature group
            
    #Parse the start and end nucleotide positions
    if( $intLoc =~ /\((\d+)-(\d+)\)/) {
        $startNuc = $1-1; #Minus one to translate to interbase
        $endNuc = $2;
        
        #Find the snoRNA gene predicted and create a feature
        #object.
        $gene = $doc->
            createAndAddFeature($ft,
                                $idcreator->
                                new_id(db => $options{project},
                                       so_type=>'gene',
                                       prefix=>$options{command_id}
                                       ), '', 'gene');
        
        #Link the gene to the analysis
        $gene->addBsmlLink('analysis', '#snoscan_analysis', 'computed_by');
        
        #Add an interval location element, to define start and stop
        #positions.
        &add_interval_loc($gene, $startNuc, $endNuc);
       
         
        #Create a snoRNA feature
        $snoRNA = $doc->createAndAddFeature($ft,
                                            $idcreator->
                                            new_id(db => $options{project},
                                                   so_type=>'snoRNA',
                                                   prefix=>$options{command_id}
                                                   ),
                                            '', 'snoRNA');
        
        #Link the snoRNA feature to the analysis
        $snoRNA->addBsmlLink('analysis', '#snoscan_analysis', 
                             'computed_by');

        #Add an Attribute element to the feature group in order
        #to document the bit score of the prediction.
        $snoRNA->addBsmlAttr("bit_score", $data{$unique});
        
        #Add an interval location element.
        &add_interval_loc($snoRNA, $startNuc, $endNuc);
        
        #Create an exon feature
        $exon = $doc->
            createAndAddFeature($ft,
                                $idcreator->
                                new_id(db => $options{project},
                                       so_type=>'exon',
                                       prefix=>$options{command_id}
                                       ),
                                '', 'exon');
        
        #Link the snoRNA feature to the analysis
        $exon->addBsmlLink('analysis', '#snoscan_analysis', 
                           'computed_by');
        
        #Add an interval location element.
        &add_interval_loc($exon, $startNuc, $endNuc);
        
        #Now create and add a feature group containing these three members
        #and an Attribute element to hold the bit_score.
        $fg = $doc->createAndAddFeatureGroup( $seq, '', 
                                              $gene->returnattr('id') );
        $fg->addBsmlFeatureGroupMember( $gene->returnattr('id'), 
                                        $gene->returnattr('class') );
        $fg->addBsmlFeatureGroupMember( $snoRNA->returnattr('id'), 
                                        $snoRNA->returnattr('class') );
        $fg->addBsmlFeatureGroupMember( $exon->returnattr('id'), 
                                        $exon->returnattr('class') );
        
        
    }
    
}

#Add the analysis documentation
my $analysis = $doc->createAndAddAnalysis( id=> 'snoscan_analysis',
                                           sourcename=> $options{'output'}
                                           );

#Write everything to the BSML document
$doc->write($options{'output'},, $options{'compress_bsml_output'});

exit;

##########################################################################
#                               SUB-ROUTINES                             #
##########################################################################    

sub check_parameters {
    my $options = shift;
    
    ## make sure input_file and output_dir were passed
    unless ( $options{input} && $options{output}) {

        #$logger->logdie("Required options are:\n\t--input_file\n".
	#"\t--output_dir\n\t--fragment_");
        if(!($options{input})) {
            print "no input\n";
        }
        if(!($options{output})) {
            print "no output\n";
        }
        pod2usage({-exitval => 2,  -message => "error message", 
		   -verbose => 1, -output => \*STDERR}); 
	
    }

    #If the files input doesn't exist, check for the gzipped version.
    #Bug 2591
    unless(-e $options{'input'}) {
        if(-e $options{'input'}.".gz") {
            $options{'input'} .= ".gz";
        }
    }
    
    ## make sure input_file exists
    if (! -e "$options{input}") {
        $logger->logdie("the input file passed ($options{input_file})".
			"cannot be read or does not exist");
        #die ("The input file does not exist ($options{input})");
    }

    ##Make sure the output file doesn't exist
    if(-e "$options{output}") {
        $logger->logdie("The output file already exists");
        #die ("The output file already exists");
    }

    if(!defined($options{'class'})) {
        $logger->logdie("The class was not defined");
        #die ("The class was not defined");
    } 

    #Parse the defline out of the fasta file.  And save it for later.
    if($options{'fasta_input'}) {
        my $mod = "<";
        if($options{'fasta_input'} =~ /\.gz$/) {
            $mod = "<:gzip";
        } elsif(!-e $options{'fasta_input'} && -e $options{'fasta_input'}.".gz") {
            $mod = "<:gzip";
        }

        open(IN, $mod, $options{'fasta_input'}) or
            $logger->logdie("Unable to open $options{'fasta_input'} ($!)");
        while(<IN>) {
            chomp;
            if(/^>(.*)/) {
                $defline = $1;
                #Assume it's a single fasta file
                last;
            }
        }
        close(IN);
    } else {
        $logger->logdie("Option fasta_input is required");
    }
        
}

sub add_interval_loc {
    my ($feat, $n, $m) = @_;
    
    ## was it found on the forward or reverse strand?
    if ($n <= $m) {
        $feat->addBsmlIntervalLoc($n, $m, 0);
    } else {
        $feat->addBsmlIntervalLoc($m, $n, 1);
    }
}
