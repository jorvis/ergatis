#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

predotar2bsml.pl - convert predotar output to BSML

=head1 SYNOPSIS

USAGE: predotar2bsml.pl --input=/path/to/predotar_file --output=/path/to/output.bsml

=head1 OPTIONS

B<--input,-i> 
    Input file file from a predotar run.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--output,-o> 
    Output BSML file (will be created, must not exist)

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used to convert the output from predotar into BSML.

=head1 INPUT

predotar can be run using multiple input sequences simultaneously, and this
script supports parsing single or multiple input result sets.  predotar
output using 'p' option for plant sequence looks like:

Predotar v. 1.03 using plant predictors
initializing predictors...
initializing sequence formats...
reading temp_multi
format assumed to be FastA
Seq     Mit     Plast   ER      None    Prediction
NCAPP1  0.01    0.01    0.98    0.02    ER
two     0.59    0.05    0.00    0.39    mitochondrial
three   Discarding three: no terminal Met
finished reading temp_multi
predicted 2 sequences out of 3

or using 'a' for animal/fungi:

Predotar v. 1.03 using animal/fungal predictors
initializing predictors...
initializing sequence formats...
reading temp_multi
format assumed to be FastA
Seq     Mit     ER      None    Prediction
NCAPP1  0.02    0.99    0.01    ER
two     0.63    0.01    0.36    mitochondrial
three   Discarding three: no terminal Met
finished reading temp_multi
predicted 2 sequences out of 3

Output from both analysis types are supported. Discarded sequences are ignored.

You define the input file using the --input option.  This file does not need any
special file extension.

Weak predictions are qualified with the term 'possibly' (eg: 'possibly ER').
By default, weak predictions and strong predictions will be treated as equally
acceptable by the parser. To modify this behaviour and reject weak predictions,
use the option --strict.

=head1 OUTPUT

After parsing the input file, a file specified by the --output option is created.  This script
will fail if it already exists.

=head1 CONTACT

Brett Whitty
bwhitty@tigr.org

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
BEGIN {
use Workflow::Logger;
use BSML::BsmlRepository;
use Papyrus::TempIdCreator;
use BSML::BsmlBuilder;
use BSML::BsmlParserTwig;
}

my %options = ();
my $results = GetOptions (\%options, 
	      'input|i=s',
              'output|o=s',
              'debug|d=s',
              'command_id=s',       ## passed by workflow
              'logconf=s',          ## passed by workflow (not used)
              'project|p=s',
              'log|l=s',
	      'help|h',
	      'strict',
                         ) || pod2usage();

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

## make sure all passed options are peachy
&check_parameters(\%options);

## we want to create ids unique to this document, which will be replaced later.  they must
##  contain the prefix that will be used to look up a real id, such as ath1.gen.15
my $next_id = 1;

## we want a new doc
my $doc = new BSML::BsmlBuilder();

## we're going to generate ids
my $idcreator = new Papyrus::TempIdCreator();

## open the input file for parsing
open (my $ifh, $options{'input'}) || $logger->logdie("can't open input file for reading");

my @results;
my @result_fields;
my @field_ids;
my $finished_flag = 0;
while (<$ifh>) {
    chomp;

    #check whitespace, no warn
    next if ( /^\s*$/ );

    ##recognize header fields
    if (/^Seq\t/) {
	##save the field ids for verifying data later
	@field_ids = split("\t");

	while (my $result_line = <$ifh>) {
		chomp $result_line;
		if ($result_line =~ /^finished reading/) {
			$finished_flag = 1;
			last;
		} else {
			my @result_fields = split("\t", $result_line);
			if ($result_line =~ /Discarding/) {
				$logger->error(
					"input sequence ".$result_fields[0]." was discarded by predotar."
				              ) if ($logger->is_error);
				next;
			}
			push (@results, \@result_fields); 
		}

	}
    }
}

## if we reached EOF without hitting finished reading line
## we need to die because something's probably gone very wrong
if (! $finished_flag) {
    $logger->logdie("Reached end of input file without seeing 'finished reading' line.\nRaw file is corrupted or in an unrecognized format.\n");
}

my @clean_results = @results; ## removed skipping of discarded sequences
#my @clean_results;
#foreach my $result_fields_ref(@results) {
	## Check if any sequences were discarded 
#	if ((scalar @$result_fields_ref) == 2 && $result_fields_ref->[1] =~ /^Discarding ([^:]+): (.*)$/) {
#        	$logger->error("The sequence with id \"$1\" was discarded by predotar.\nThe reason given was: $2\n") if ($logger->is_error);
		
#	} elsif (scalar @$result_fields_ref != scalar @field_ids) {
#        	$logger->error("The number of fields in the following line did not match the results header:\n".join("\t", @$result_fields_ref)."\n") if ($logger->is_error);
#	} else {
#		push (@clean_results, $result_fields_ref);
#	} 
#}
@results = ();

my $result_count = scalar @clean_results;

## GO mappings tables for output
my %go_term;
my %go_id;
if (defined($options{strict})) {
	%go_term = 	( 
		  	  'mitochondrial'   		=> 'mitochondrion',
			  'possibly mitochondrial'   	=> 'cellular component unknown',
		  	  'plastid' 			=> 'plastid',
			  'possibly plastid' 		=> 'cellular component unknown',
		  	  'ER'	  			=> 'endoplasmic reticulum',
			  'possibly ER'	  		=> 'cellular component unknown',
		  	  'none'  			=> 'cellular component unknown',
	  		);
	%go_id = 	( 
	          	  'mitochondrial'   		=> 'GO:0005739',
			  'possibly mitochondrial'   	=> 'GO:0008372',
		  	  'plastid' 			=> 'GO:0009536',
			  'possibly plastid' 		=> 'GO:0008372',
		  	  'ER'	  			=> 'GO:0005783',
			  'possibly ER'	  		=> 'GO:0008372',
		  	  'none'  			=> 'GO:0008372',
		  	);
} else {
	%go_term = 	( 
			  'mitochondrial'   		=> 'mitochondrion',
			  'possibly mitochondrial'   	=> 'mitochondrion',
			  'plastid' 			=> 'plastid',
			  'possibly plastid' 		=> 'plastid',
			  'ER'	  			=> 'endoplasmic reticulum',
			  'possibly ER'	  		=> 'endoplasmic reticulum',
		  	  'none'  			=> 'cellular component unknown',
	  		);
	%go_id =	( 
	          	  'mitochondrial'   		=> 'GO:0005739',
		          'possibly mitochondrial'   	=> 'GO:0005739',
		  	  'plastid' 			=> 'GO:0009536',
			  'possibly plastid' 		=> 'GO:0009536',
		  	  'ER'	  			=> 'GO:0005783',
			  'possibly ER'	  		=> 'GO:0005783',
		  	  'none'  			=> 'GO:0008372',
		  	);
}

## loop through each of the matches that we found
for (my $i = 0; $i < $result_count; $i++) {
	my @result_fields = @{$clean_results[$i]};
	my @result_field_ids = @field_ids;
	
	## pull the sequence id off the front of the array
	my $sequence_id = shift @result_fields;
	## and discard the label
	shift @result_field_ids;
	## pull the prediction off the end of the array
	my $prediction = pop @result_fields;
	## and again discard the label
	pop @result_field_ids;
	
	my $seq = $doc->createAndAddSequence($sequence_id, undef, '', 'aa', 'polypeptide');
        $seq->addBsmlLink('analysis', '#predotar_analysis', 'input_of');
	unless ($prediction =~ /^Discarding ([^:]+): (.*)$/) {
    	my $ft  = $doc->createAndAddFeatureTable($seq);
	my $new_id = $idcreator->new_id( db      => $options{project},
                                         so_type => 'transit_peptide',                                           	 prefix  => $options{command_id},
                                       );
        my $feature = $doc->createAndAddFeature($ft, $new_id, '', 'transit_peptide');
	$feature->addBsmlLink('analysis', '#predotar_analysis', 'computed_by');
	
	unless (defined($go_id{$prediction})) {
		$logger->logdie("predotar prediction '$prediction' was not an expected result");
	}
	
	my $attribute_array_ref;
	push( @{$attribute_array_ref}, { name    => 'GO',
		          		 content => $go_id{$prediction}}
	    );
	push( @{$attribute_array_ref}, { name    => 'IEA',
		          		 content => 'predotar prediction'}
	    );
        $feature->addBsmlAttributeList($attribute_array_ref);

	my $attribute_count = scalar @result_field_ids;
	for (my $j = 0; $j < $attribute_count; $j++) {
		$doc->createAndAddBsmlAttribute($feature, $result_field_ids[$j], $result_fields[$j]);
	}
	}
}

## add the analysis element
my $analysis = $doc->createAndAddAnalysis(
                            id => 'predotar_analysis',
                            sourcename => $options{'output'},
                          );

$doc->createAndAddBsmlAttribute( $analysis, 'version', 'current' );
$doc->createAndAddBsmlAttribute( $analysis, 'algorithm', 'predotar' );

## now write the doc
$doc->write($options{'output'});

exit;

sub check_parameters {
    
    ## make sure input file exists
    if (! -e $options{'input'}) { $logger->logdie("input file $options{'input'} does not exist") }
    
    ## make sure output file doesn't exist yet
    if (-e $options{'output'}) { $logger->logdie("can't create $options{'output'} because it already exists") }
    
    $options{'project'}    = 'unknown' unless ($options{'project'});
    $options{'command_id'} = '0' unless ($options{'command_id'});
    
    return 1;
}

