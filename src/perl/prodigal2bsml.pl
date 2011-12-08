#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

=head1 NAME

prodigal2bsml.pl - Creates a bsml document from prodigal raw output

=head1 SYNOPSIS

USAGE: prodigal2bsml.pl
            --input_list=/path/to/some/prodigal.raw.list
            --input_file=/path/to/some/prodigal.raw
            --output=/path/to/transterm.bsml
            --id_repository=/path/to/id_repository
            --fasta_input=/path/to/prodigal/input.fsa
          [ --log=/path/to/file.log
            --debug=4
            --help
          ]

=head1 OPTIONS

B<--input_list,-i>
    Input list of prodigal raw output files (.predict)

B<--input_file,-f>
    Input prodigal raw file (.predict)

B<--output,-o>
    The output bsml file.

B<--id_repository,-r>
    Id repository for use by Workflow::IdGenerator.pm

B<--fasta_input,-a>
    The input file that was used as input for the prodigal run

B<--log,-l>
    Logfile.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

This script is used to convert the output from a prodigal search into BSML.

=head1  INPUT

You define the input file using the --input option.  This file does not need any
special file extension.  The regular output of prodigal looks like this:


# Sequence Data: seqnum=1;seqlen=5127136;seqhdr="gec3273.pseudomolecule.11"
# Model Data: version=Prodigal.v2.50;run_type=Single;model="Ab initio";gc_cont=50.70;transl_table=4;uses_sd=1
 >17_20759_21505_+
 >18_21508_22068_-
 >19_22103_22444_-
 >20_22579_22905_+

Where the columns from left to right separated by '_' contain:
prediction_id   start_position   end_position   strand 
start_position is always less than the end_position

=head1 OUTPUT

After parsing the input file, a file specified by the --output option is created.  This script
will fail if it already exists.  The file is created, and temporary IDs are created for
each result element.  They are only unique to the document, and will need to be replaced
before any database insertion.

Base positions from the input file are renumbered so that positions start at zero.

NOTE:
This script is created by modifying Kevin's glimmer32bsml.pl for prodigal

=head1  CONTACT

    Sonia Agrawal
    sagrawal@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Chado::Gene;
use BSML::GenePredictionBsml;
use Ergatis::IdGenerator;
use Ergatis::Logger;

####### GLOBALS AND CONSTANTS ###########
my @inputFiles;               #Holds input files
my $project;                  #The project (ex aa1)
my $output;                   #Output file
my $idMaker;                  #The Workflow::IdGenerator
my $bsml;                     #BSML::BsmlBuilder object object.
my $data;                     #Holds parsed prodigal information
my $inputFsa;                 #The fasta file input to prodigal
my $debug;                    #The debug variable
my $length;
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
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
		'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

# Check the options.
&check_parameters(\%options);

foreach my $file (@inputFiles) {
	$data = &parseProdigalData($file);
	$bsml = &generateBsml($data);
}

$bsml->writeBsml($output);

exit(0);

######################## SUB ROUTINES #######################################
sub parseProdigalData {
	my ($file) = @_;
	my $genes;
	my $foundId = 0;

	open(IN, "<$file") or &_die("Unable to open $file");

	while(<IN>) {
		chomp($_);	
		# return header to first space to be consistent with glimmer3
                if(/seqhdr\=\"(\S+)[^\"]*\"/g) {
			$foundId = $1;
			if( $project eq 'parse' ) {
				$project = $1 if($foundId =~ /^(\w+?)\./);
				$logger->logdie("Could not parse project name out of id $foundId.") unless($project);
			}

		} elsif($foundId) {
			next if(/^#/);
			my @cols = split(/_/,$_);
			my $strand = ($cols[3] eq '+') ? 0 : 1;
# Correcting for interbase numbering
			$cols[1]--;
#Create some genes and push them onto the $genes array
			my $tmp = new Chado::Gene( $idMaker->next_id( 'type' => 'gene',
					'project' => $project),
				$cols[1], $cols[2], $strand,
				$foundId);

			foreach my $type(qw(exon CDS transcript polypeptide)) {
				my $typeid =$idMaker->next_id( 'type' => $type,
					'project' => $project);
				$tmp->addFeature($typeid, $cols[1], $cols[2], $strand,
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
	my $doc = new BSML::GenePredictionBsml( 'prodigal', $inputFsa );

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

	unless($options{'output'}) {
		$error .= "Option output is required\n";
	} else {
		$output = $options{'output'};
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
		$error .= "$options{'fasta_input'} (fasta_input) does not exist\n" unless(-e $options{'fasta_input'});
		$inputFsa = $options{'fasta_input'};

#Find the length of the sequence
		open(IN, "<$options{'fasta_input'}") or $logger->logdie("Can't open $options->{'fasta_input'}");
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
				$curId =~ s/\s+$//;
			} else {
				$seq.= $_;
			}

		}

		$length->{$curId} = length($seq);
		close(IN);
	}

	if( $options{'project'} ) {
		$project = $options->{'project'};
	} else {
		$project = "parse";
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
