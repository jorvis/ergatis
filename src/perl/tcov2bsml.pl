use lib "/export/CVS/bsml/src/";
use BSML::BsmlBuilder;
use BSML::BsmlReader;
use BSML::BsmlParserTwig;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

# Utility to convert coverage data stored in the TIGR tcov format from
# getCoverage into BSML

# BSML encoding to support externally linked BSML sequences with sequences
# containing consensus quality scores, read depth, and most commonly occuring
# alternative base mapped to them via the numbering element.

my %options = ();
my $results = GetOptions (\%options, 
			  'tcov|t=s',
			  'bsml_repository_path|b=s',
			  'output|o=s',  
			  'log|l=s',
			  'debug=s',
			  'help|h') || pod2usage();

my $tcovFilePath = $options{'tcov'};
my $repositoryPath = $options{'bsml_repository_path'};
my $outputFile = $options{'output'};

# filename is expected to encode the assembly id
# ie /usr/local/annotation/MYCOTB/chauser/bmt_coverage_data/443.tcov

$tcovFilePath =~ /[\S]*\/([\d]*).tcov/;
my $assemblyId = $1;

open( TCOV, $tcovFilePath ) or die "Could not open $tcovFilePath\n";

# an array containing the consensus quality score at each position
my @consensus_quality_score;

# an array containing the read depth at each position
my @read_depth;

# an array containing the most commonly occuring alternative base at 
# each position
my @most_common_alternate_base;

while( my $line = <TCOV> )
{
    # tcov input is space delimited with each line representing the following
    # base position | consensus base | consensus quality value | base calls covering the position 
    # | quality values of the bases | Indices of the corresponding bases and quality values
    # see http://intranet/software_docs/getCoverage.html for a more thorough description

    my ($pos, $consensus_base, $consensus_quality_value, $base_calls, $quality_values, $indices ) =
	split( ' ', $line );

    # does not support ambiguous base calls 
    if( $consensus_base =~ /[AGCT]/ )
    {
	$consensus_quality_score[$pos-1] = $consensus_quality_value;
	$read_depth[$pos-1] = length( $base_calls );

	my $baseCount = { 'A' => 0, 'G' => 0, 'C' => 0, 'T' => 0 };

	for( my $i=0; $i<length( $base_calls ); $i++ )
	{
	    my $base = substr( $base_calls, $i, 1 );
	    $baseCount->{$base} += 1;
	}

	# set the consensus base to -1 to get the second most common base

	$baseCount->{$consensus_base} = -1;

	# if the only base encountered is the consensus base and therefore the largest count is
	# an unencountered base, the consensus base will be used

	if( ($baseCount->{'A'}) >= $baseCount->{'G'} && ($baseCount->{'A'} >= $baseCount->{'C'}) &&  ($baseCount->{'A'} >= $baseCount->{'T'}) )
	{
	    if( $baseCount->{'A'} == 0 )
	    {
		$most_common_alternate_base[$pos-1] = $consensus_base;
	    }
	    else
	    {
		$most_common_alternate_base[$pos-1] = 'A';
	    }
	}

	if( ($baseCount->{'G'}) >= $baseCount->{'A'} && ($baseCount->{'G'} >= $baseCount->{'C'}) &&  ($baseCount->{'G'} >= $baseCount->{'T'}) )
	{
	    if( $baseCount->{'G'} == 0 )
	    {
		$most_common_alternate_base[$pos-1] = $consensus_base;
	    }
	    else
	    {
		$most_common_alternate_base[$pos-1] = 'G';
	    }
	}

	if( ($baseCount->{'C'}) >= $baseCount->{'G'} && ($baseCount->{'C'} >= $baseCount->{'A'}) &&  ($baseCount->{'C'} >= $baseCount->{'T'}) )
	{
	    if( $baseCount->{'C'} == 0 )
	    {
		$most_common_alternate_base[$pos-1] = $consensus_base;
	    }
	    else
	    {
		$most_common_alternate_base[$pos-1] = 'C';
	    }
	}

	if( ($baseCount->{'T'}) >= $baseCount->{'G'} && ($baseCount->{'T'} >= $baseCount->{'C'}) &&  ($baseCount->{'T'} >= $baseCount->{'A'}) )
	{
	    if( $baseCount->{'T'} == 0 )
	    {
		$most_common_alternate_base[$pos-1] = $consensus_base;
	    }
	    else
	    {
		$most_common_alternate_base[$pos-1] = 'T';
	    }
	}
    }
}

# bsml builder class to construct the coverage document
my $bsmlDoc = new BSML::BsmlBuilder;

my $asmbl_length = @consensus_quality_score;


# add the reference sequence and it's link in the repository

my $seq = $bsmlDoc->createAndAddExtendedSequenceN( id => "bmt_".$assemblyId."_assembly",
						   class => "CONTIG",
						   length => $asmbl_length );

$bsmlDoc->createAndAddSeqDataImportN( seq => $seq,
				      format => "BSML",
                                      source => $repositoryPath."/bmt_".$assemblyId."_assembly.bsml",
                                      id => "_bmt_".$assemblyId."_assembly");


# add the consensus quality score sequence

my $tmpstr = '';

for( my $i=0; $i<$asmbl_length; $i++ )
{
    $tmpstr .= ":$consensus_quality_score[$i]";
}

my $consensus_quality_score_seq = $bsmlDoc->createAndAddExtendedSequenceN( id => "bmt_".$assemblyId."_assembly_consensus_quality_score",
									   class => "quality_score",
									   length => $asmbl_length );

$bsmlDoc->createAndAddSeqData( $consensus_quality_score_seq, $tmpstr );

# map it to the reference sequence

my $numbering = $bsmlDoc->createAndAddNumbering( seq => $consensus_quality_score_seq,
						 seqref => "bmt_".$assemblyId."_assembly",
						 refnum => 0,
						 ascending => 1 );

# add the sequence depth score sequence

$tmpstr = '';

for( my $i=0; $i<$asmbl_length; $i++ )
{
    $tmpstr .= ":$read_depth[$i]";
}

my $read_depth_score_seq = $bsmlDoc->createAndAddExtendedSequenceN( id => "bmt_".$assemblyId."_assembly_read_depth",
									   class => "depth_score",
									   length => $asmbl_length );

$bsmlDoc->createAndAddSeqData( $read_depth_score_seq, $tmpstr );

# map it to the reference sequence

my $numbering = $bsmlDoc->createAndAddNumbering( seq => $read_depth_score_seq,
						 seqref => "bmt_".$assemblyId."_assembly",
						 refnum => 0,
						 ascending => 1 );

# add the alternate base sequence

$tmpstr = '';

for( my $i=0; $i<$asmbl_length; $i++ )
{
    $tmpstr .= ":$most_common_alternate_base[$i]";
}

my $alt_base_score_seq = $bsmlDoc->createAndAddExtendedSequenceN( id => "bmt_".$assemblyId."_assembly_alt_base",
									   class => "alt_base",
									   length => $asmbl_length );

$bsmlDoc->createAndAddSeqData( $alt_base_score_seq, $tmpstr );

# map it to the reference sequence

my $numbering = $bsmlDoc->createAndAddNumbering( seq => $alt_base_score_seq,
						 seqref => "bmt_".$assemblyId."_assembly",
						 refnum => 0,
						 ascending => 1 );

$bsmlDoc->write( $outputFile ); 





