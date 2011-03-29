#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
if 0; # not running under some shell

BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

###########################################################
# POD DOCUMENTATION                                       #
###########################################################
=head1 NAME

pmarks2bsml.pl - Creates a bsml document for pmarks location on a pseudomolecule

=head1 SYNOPSIS

USAGE: pmarks2bsml.pl
            --input_file=/path/to/some/pseudomol.pmarks
            --output=/path/to/pmarks.bsml
            --id_repository=/path/to/id_repository
          [ --fasta_input=/path/to/pseudomolecule.fsa
            --log=/path/to/file.log
            --debug=4
            --help
          ]

=head1 OPTIONS

B<--input_file,-f>
    Input pmarks positions file

B<--output,-o>
    The output bsml file.

B<--id_repository,-r>
    Id repository for use by Workflow::IdGenerator.pm

B<--fasta_input,-a>
    The input file that was used to locate pmarks location in a pseudomolecule

B<--log,-l>
    Logfile.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

This script is used to convert the output from a create_pseudomolecule script specifying the pmarks location into BSML.

=head1  INPUT

You define the input file using the --input option.  This file does not need any
special file extension.  The regular output of create_pseudomolecule looks like this:


>zma1.assembly.5808
1	30
174	204 
430	460  
1258    1288 
1541    1571 
2153    2183 

Where the columns from left to right contain:
start_position   end_position 

=head1 OUTPUT

After parsing the input file, a file specified by the --output option is created. 
The file is created, and temporary IDs are created for each result element.  
They are only unique to the document, and will need to be replaced
before any database insertion.

=head1  CONTACT

    Sonia Agrawal
    sagrawal@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Basename;
use BSML::BsmlBuilder;
use Ergatis::IdGenerator;
use Ergatis::Logger;

####### GLOBALS AND CONSTANTS ###########
my $inputFile;               #Holds input files
my $project;                  #The project (ex aa1)
my $output;                   #Output file
my $idMaker;                  #The Workflow::IdGenerator
my $doc;                      #BSML::BsmlBuilder object object.
my $inputFsa;                 #The fasta file input (pseudomolecule)
my $debug;                    #The debug variable
my $length;
my $ft;
my $seq;
my $thing;
########################################

my %options = ();
my $results = GetOptions (\%options, 
		'input_file|f=s',
		'output|o=s',
		'project|p=s',
		'id_repository|r=s',
		'fasta_input|a=s',
		'log|l=s',
		'debug=s',
		'help|h') || &_pod;

#Setup the logger
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

# Check the options.
&check_parameters(\%options);

$doc = new BSML::BsmlBuilder();

&parsePmarksData($inputFile);

$doc->createAndAddAnalysis(
				id => 'pseudomolecule_analysis',
				sourcename => $output,
				version => 'current',
				algorithm => 'pseudomolecule_generation',
				program => 'pmark_spacer'
			  );

$doc->write($output);

exit(0);

######################## SUB ROUTINES #######################################
sub parsePmarksData {
	my ($file) = @_;
	my $genes;
	my $foundId = 0;
	open(IN, "<$file") or $logger->logdie("Unable to open $file file");
	while(<IN>) {
		if(/^>(.*?)\s/) {
			$foundId = $1;
			print "$foundId is foundId\n";
			if( $project eq 'parse' ) {
				$project = $1 if($foundId =~ /^(\w+?)\./);
		                $logger->logdie("Could not parse project name out of id $foundId.") unless($project);
			}
			$seq = $doc->createAndAddSequence($foundId, undef, $length->{$foundId}, 'dna', 'assembly');
			$doc->createAndAddBsmlAttribute($seq,'defline',$foundId);
			$doc->createAndAddSeqDataImport($seq,'fasta',$inputFsa,'',$foundId);
			$seq->addBsmlLink('analysis', '#pseudomolecule_analysis', 'input_of');
			$ft = $doc->createAndAddFeatureTable($seq);

		} elsif($foundId) {
			print "$_\n";
			my @cols = split(/\s+/,$_);
			my $id = $idMaker->next_id( type => 'pmark_spacer', project => $project);
			$thing = $doc->createAndAddFeature( $ft, $id, '', 'pmark_spacer' );
			$thing->addBsmlLink('analysis', '#pseudomolecule_analysis', 'computed_by');
			$thing->addBsmlIntervalLoc($cols[0], $cols[1], 0);	    
		} else {
			$logger->logdie("Did not find the id in the file");
		}
	}
	close(IN);
}

sub check_parameters {
	my $options = shift;
	if($options{'help'}) {
		pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
	}
	
	if($options{'output'}) {
		$output = $options{'output'};
	} else {
		$output = $inputFile.".bsml"; 
	}

	unless($options{'id_repository'}) {
	        $logger->logdie("Option id_repository is required.  Please see Ergatis::IdGenerator for details.\n");
	} else {
		$idMaker = new Ergatis::IdGenerator( 'id_repository' => $options{'id_repository'} );
		$idMaker->set_pool_size( 'pmark_spacer'      => 20);
	}

	if(defined($options{'fasta_input'}) && (-e $options{'fasta_input'})) {
		$inputFsa = $options{'fasta_input'};
#Find the length of the sequence
		open(IN, "<$inputFsa") or $logger->logdie("Cannot open $inputFsa file for reading\n$!\n");
		my $seq = "";
		my $curId;
		while(<IN>) {
			chomp;
			if(/^>(\S+)/) {
				if(defined($curId)) {
					$length->{$curId} = length($seq);
					$seq = "";
				}
				$curId = $1;
			} else {
				$seq.= $_;
			}
		}
		$length->{$curId} = length($seq);
		close(IN);
	} else {
#		my ($file_base,$file_dir,$file_ext) = fileparse($inputFile,qr/\.[^.]*/);	
#		$inputFsa = $file_base;
		$logger->logdie("--fasta_input is mandatory and $options{'fasta_input'} file should exist.");
	}
	
	if($options{'input_file'}) {
		$inputFile = $options{'input_file'};
	} else {
		$inputFile = $inputFsa.".pmarks";
		$logger->logdie("$inputFile file does not exist. Either pass it as a parameter or execute pseudomolecule creation script\n") unless(-e $inputFile);
	}

	if( $options{'project'} ) {
		$project = $options->{'project'};
	} else {
		$project = "parse";
	}

	if($options{'debug'}) {
		$debug = $options{'debug'};
	}
}
