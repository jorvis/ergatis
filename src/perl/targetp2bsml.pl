#!/usr/local/packages/perl-5.8.5/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

targetp2bsml.pl - convert TargetP output to BSML

=head1 SYNOPSIS

USAGE: targetp2bsml.pl --input=/path/to/targetp_file --output=/path/to/output.bsml

=head1 OPTIONS

B<--input,-i> 
    Input file file from a targetp run.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--output,-o> 
    Output BSML file (will be created, must not exist)

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used to convert the output from TargetP into BSML.

=head1 INPUT

TargetP can be run using multiple input sequences simultaneously, and this
script supports parsing single or multiple input result sets.
 
You define the input file using the --input option.  This file does not need any
special file extension.

=head1 OUTPUT

After parsing the input file, a file specified by the --output option is created.  This script
will fail if it already exists.

=head1 CONTACT

Brett Whitty
bwhitty@tigr.org

=cut

use strict;
use Log::Log4perl qw(get_logger);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BSML::BsmlBuilder;
use BSML::BsmlReader;
use BSML::BsmlParserTwig;
use BSML::BsmlRepository;
use Papyrus::TempIdCreator;
use Pod::Usage;
use Workflow::Logger;

my %options = ();
my $results = GetOptions (\%options, 
			  'input|i=s',
              'output|o=s',
			  'sequence_id|s=s',
              'debug|d=s',
              'command_id=s',       ## passed by workflow
              'logconf=s',          ## passed by workflow (not used)
              'project|p=s',
              'log|l=s',
			  'help|h') || pod2usage();

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

my $plant_flag = 0;

my @result_line_ref;
while (<$ifh>) {
    chomp;
	if (/using plant networks/i) {
		$plant_flag = 1;
	}
    if ($plant_flag && /^cutoff\s+([^ ]+)\s+([^ ]+)\s+([^ ]+)\s+([^ ]+)/) {
    	## the cutoff values used will be added to the analysis block
		## so we don't really need to parse them here
		my $ctp_cutoff   = $1;
		my $mtp_cutoff   = $2;
		my $sp_cutoff    = $3;
		my $other_cutoff = $4;
	} elsif (!$plant_flag && /^cutoff\s+([^ ]+)\s+([^ ]+)\s+([^ ]+)/) {
    	## the cutoff values used will be added to the analysis block
		## so we don't really need to parse them here
		my $mtp_cutoff   = $1;
		my $sp_cutoff    = $2;
		my $other_cutoff = $3;
	}

	## skip any lines we're not looking for specifically
	next if ( ! /^-{70}$/ );

    ##recognize start of results for an individual sequence
    if (/^-{70}/) {
		$logger->debug("start of results table found") if($logger->is_debug());  
		while (<$ifh>) {
			chomp;
			
			if (/^-{70}/) {
				last;
			}
			if ($plant_flag && /^(.{20})\s+(\d+)\s+([^ ]+)\s+([^ ]+)\s+([^ ]+)\s+([^ ]+)\s+([^ ]+)\s+([^ ]+)\s+([^ ]+)?$/) {
			
				my $result_line = {
									'name'  => $1,
									'len'   => $2,
									'ctp'   => $3,
									'mtp'   => $4,
									'sp'    => $5,
									'other' => $6,
									'loc'   => $7,
									'rc'    => $8,
									'tplen' => $9,
					  			  };
				push (@result_line_ref, $result_line);
			} elsif (!$plant_flag && /^(.{20})\s+(\d+)\s+([^ ]+)\s+([^ ]+)\s+([^ ]+)\s+([^ ]+)\s+([^ ]+)\s+([^ ]+)?$/) {
			
				my $result_line = {
									'name'  => $1,
									'len'   => $2,
									'mtp'   => $3,
									'sp'    => $4,
									'other' => $5,
									'loc'   => $6,
									'rc'    => $7,
									'tplen' => $8,
					  			  };
				push (@result_line_ref, $result_line);
			} else {
				if ($plant_flag) {print STDERR "died in plant mode\n";} else {print STDERR "died in nonplant mode\n";}
				$logger->logdie("failed parsing result line:\n$_\n");
			}
		}		
	}
}

## GO mappings tables for output
my %go_term = 	( 'C' => 'plastid',
		  'M' => 'mitochondrion',
		  'S' => 'extracellular space',
		  '_' => 'cytoplasm',
		  '*' => 'cellular component unknown',
		  '?' => 'cellular component unknown',
	  	);
my %go_id =     ( 'C' => 'GO:0009536',
		  'M' => 'GO:0005739',
		  'S' => 'GO:0005615',
		  '_' => 'GO:0005737',
		  '*' => 'GO:0008372',
		  '?' => 'GO:0008372',
	  	);

if ($options{'sequence_id'} && scalar(@result_line_ref) > 1) {
	$logger->logdie("Sequence ID was provided as an argument, but input file has multiple result lines. Sequence ID will override IDs read from the raw file.");  
}
		
    foreach my $line_ref(@result_line_ref) {
		$line_ref->{'name'} =~ s/\s+$//;
		my $seq_id = '';
		if ($options{'sequence_id'}) {
			$seq_id = $options{'sequence_id'};
		} else {
			$seq_id = $line_ref->{'name'};
		}
    	my $seq = $doc->createAndAddSequence(
			$seq_id,
		   	undef, 
			'', 
			'aa',
		   	'polypeptide'
			);
	   	$seq->addBsmlLink('analysis', '#targetp_analysis', 'input_of');
	
		if ($line_ref->{'loc'} ne '*' && $line_ref->{'loc'} ne '?') {
        	my $ft  = $doc->createAndAddFeatureTable($seq);
			my $new_id = $idcreator->new_id( 
				db      => $options{project},
                so_type => 'transit_peptide',
				prefix  => $options{command_id},
                                                        );
        		my $feature = $doc->createAndAddFeature($ft, $new_id, '', 'transit_peptide');
	   	   	$feature->addBsmlLink('analysis', '#targetp_analysis', 'computed_by');
			my $attribute_array_ref;
			push( @{$attribute_array_ref}, { 
				name    => 'GO',
		          	content => $go_id{$line_ref->{'loc'}}
			       }
	    		    );
			push( @{$attribute_array_ref}, { 
				name    => 'IEA',
		          	content => 'targetp prediction'
			                               }
	    		    );
        $feature->addBsmlAttributeList($attribute_array_ref);


			
			foreach my $k(keys %{$line_ref}) {
				$doc->createAndAddBsmlAttribute($feature, $k, $line_ref->{$k});
			}
		}
	}
   
my $analysis = $doc->createAndAddAnalysis(
	                            id => 'targetp_analysis',
                            sourcename => $options{'output'},
                          );

## now write the doc
$doc->write($options{'output'});

exit();

sub check_parameters {
    
    ## make sure input file exists
    if (! -e $options{'input'}) { $logger->logdie("input file $options{'input'} does not exist") }
    
    ## make sure output file doesn't exist yet
    if (-e $options{'output'}) { $logger->logdie("can't create $options{'output'} because it already exists") }
    
    $options{'project'}    = 'unknown' unless ($options{'project'});
    $options{'command_id'} = '0' unless ($options{'command_id'});
    
    return 1;
}

## parse_position parses the position field to give an array containing
## start [and stop] sites
sub parse_position {
	my ($position) = @_;
	if ($position =~ /^\s*([^\s]+)\s*$/) {
		return split("-", $1);
	} else {
		return ();
	}
}

