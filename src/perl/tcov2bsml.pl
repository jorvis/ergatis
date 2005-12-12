#!/usr/local/bin/perl

# this is temporary for development
# use lib "/export/CVS/bsml/src/";

=head1  NAME 

tcov2bsml.pl  - convert tcov files into BSML documents

=head1 SYNOPSIS

USAGE:  tcov2bsml.pl -t tcov_file -b bsml_repository_path -o output_file_bsml

=head1 OPTIONS

=over 4

=item *

B<--tcov,-t> [REQUIRED] file containing tcov data

=item *

B<--output,-o> [REQUIRED] output BSML file containing coverage information

=item * 

B<--bsml_repository,-b> [REQUIRED] filepath to the bsml repository 

=item *

B<--chunk_size,-c> [OPTIONAL] optional chunk size (default is 1000 bases) 

=item *

B<--help,-h> This help message

=back

=head1   DESCRIPTION

tcov2bsml.pl is designed to convert information in tcov files (coverage data)
into BSML documents.

=cut

BEGIN {
    require '/usr/local/devel/ANNOTATION/cas/lib/site_perl/5.8.5/BSML/BsmlReader.pm';
    import BSML::BsmlReader;
    require '/usr/local/devel/ANNOTATION/cas/lib/site_perl/5.8.5/BSML/BsmlBuilder.pm';
    import BSML::BsmlBuilder;
    require '/usr/local/devel/ANNOTATION/cas/lib/site_perl/5.8.5/BSML/BsmlParserTwig.pm';
    import BSML::BsmlParserTwig;
}

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
			  'database|d=s',
			  'chunk_size|c=s',
			  'output|o=s',  
			  'log|l=s',
			  'debug=s',
			  'help|h') || pod2usage();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

if( !( $options{'bsml_repository_path' } ) )
{
    pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
}

if( !( $options{'tcov' } ) )
{
    pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
}

if( !( $options{'database' } ) )
{
    pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
}

my $tcovFilePath = $options{'tcov'};
my $repositoryPath = $options{'bsml_repository_path'};
my $outputFile = $options{'output'};
my $chunk_count;
my $database = $options{'database'};

if( $options{'chunk_size'} )
    {
	$chunk_count = $options{'chunk_size'};
    }
else
{
    $chunk_count = 1000;
}


# filename is expected to encode the assembly id
# ie /usr/local/annotation/MYCOTB/chauser/bmt_coverage_data/443.tcov

$tcovFilePath =~ /[\S]*_([\d]*).tcov/;
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

my $seq = $bsmlDoc->createAndAddExtendedSequenceN( id => $database."_".$assemblyId."_assembly",
						   class => "contig",
						   length => $asmbl_length );

$seq->addattr( 'class', "assembly" );

$bsmlDoc->createAndAddSeqDataImportN( seq => $seq,
				      format => "BSML",
                                      source => $repositoryPath."/".$database."_".$assemblyId."_assembly.bsml",
                                      id => "_".$database."_".$assemblyId."_assembly");


# add the consensus quality score sequence

my $tmpstr = '';
my $tmpstr_count = 0;
my $offset = 0;
my $count = 0;


for( my $i=0; $i<$asmbl_length; $i++ )
{
    $tmpstr .= ":$consensus_quality_score[$i]";
    $tmpstr_count++;

    if( ((($i+1) % $chunk_count) == 0))
    {
	my $consensus_quality_score_seq = $bsmlDoc->createAndAddExtendedSequenceN( id => $database."_".$assemblyId."_assembly_consensus_quality_score_".$count,
										   class => "consensus_quality_value",
										   length => $tmpstr_count );
	
	$consensus_quality_score_seq->addattr( 'class', "consensus_quality_value" );
    
	$bsmlDoc->createAndAddSeqData( $consensus_quality_score_seq, $tmpstr );
	
	my $numbering = $bsmlDoc->createAndAddNumbering( seq => $consensus_quality_score_seq,
						 seqref => $database."_".$assemblyId."_assembly",
						 refnum => $offset,
						 ascending => 1 );

	$offset += $chunk_count;
	$count++;

	$tmpstr = '';
	$tmpstr_count = 0;
    }
}

my $consensus_quality_score_seq = $bsmlDoc->createAndAddExtendedSequenceN( id => $database."_".$assemblyId."_assembly_consensus_quality_score_".$count,
									   class => "consensus_quality_value",
									   length => $tmpstr_count );

$consensus_quality_score_seq->addattr( 'class', "consensus_quality_value" );

$bsmlDoc->createAndAddSeqData( $consensus_quality_score_seq, $tmpstr );

# map it to the reference sequence

my $numbering = $bsmlDoc->createAndAddNumbering( seq => $consensus_quality_score_seq,
						 seqref => $database."_".$assemblyId."_assembly",
						 refnum => $offset,
						 ascending => 1 );

# add the sequence depth score sequence

$tmpstr = '';
$tmpstr_count = 0;
$offset = 0;
$count = 0;


for( my $i=0; $i<$asmbl_length; $i++ )
{
    $tmpstr .= ":$read_depth[$i]";
    $tmpstr_count++;

    if( ((($i+1) % $chunk_count) == 0))
    {
	my $read_depth_score_seq = $bsmlDoc->createAndAddExtendedSequenceN( id => $database."_".$assemblyId."_assembly_read_depth_".$count,
										   class => "read_depth",
										   length => $tmpstr_count );
	    
	$read_depth_score_seq->addattr( 'class', "read_depth" );

	$bsmlDoc->createAndAddSeqData( $read_depth_score_seq, $tmpstr );
	
	my $numbering = $bsmlDoc->createAndAddNumbering( seq => $read_depth_score_seq,
						 seqref => $database."_".$assemblyId."_assembly",
						 refnum => $offset,
						 ascending => 1 );

	$offset += $chunk_count;
	$count++;

	$tmpstr = '';
	$tmpstr_count = 0;
    }
}

my $read_depth_score_seq = $bsmlDoc->createAndAddExtendedSequenceN( id => $database."_".$assemblyId."_assembly_read_depth_".$count,
									   class => "read_depth",
									   length => $tmpstr_count );

$read_depth_score_seq->addattr( 'class', "read_depth" );

$bsmlDoc->createAndAddSeqData( $read_depth_score_seq, $tmpstr );

# map it to the reference sequence

my $numbering = $bsmlDoc->createAndAddNumbering( seq => $read_depth_score_seq,
						 seqref => $database."_".$assemblyId."_assembly",
						 refnum => $offset,
						 ascending => 1 );

# add the alternate base sequence

$tmpstr = '';
$tmpstr_count = 0;
$offset = 0;
$count = 0;

for( my $i=0; $i<$asmbl_length; $i++ )
{
    $tmpstr .= "$most_common_alternate_base[$i]";
    $tmpstr_count++;

      if( ((($i+1) % $chunk_count) == 0))
      {
	my $alt_base_score_seq = $bsmlDoc->createAndAddExtendedSequenceN( id => $database."_".$assemblyId."_assembly_alt_base_".$count,
										   class => "alternative_base",
										   length => $tmpstr_count );
	    
	$alt_base_score_seq->addattr( 'class', "alternative_base" );


	$bsmlDoc->createAndAddSeqData( $alt_base_score_seq, $tmpstr );
	
	my $numbering = $bsmlDoc->createAndAddNumbering( seq => $alt_base_score_seq,
						 seqref => $database."_".$assemblyId."_assembly",
						 refnum => $offset,
						 ascending => 1 );

	$offset += $chunk_count;
	$count++;

	$tmpstr = '';
	$tmpstr_count = 0;
    }


}

my $alt_base_score_seq = $bsmlDoc->createAndAddExtendedSequenceN( id => $database."_".$assemblyId."_assembly_alt_base_".$count,
									   class => "alternative_base",
									   length => $tmpstr_count );

$alt_base_score_seq->addattr( 'class', "alternate_base" );

$bsmlDoc->createAndAddSeqData( $alt_base_score_seq, $tmpstr );

# map it to the reference sequence

my $numbering = $bsmlDoc->createAndAddNumbering( seq => $alt_base_score_seq,
						 seqref => $database."_".$assemblyId."_assembly",
						 refnum => $offset,
						 ascending => 1 );

$bsmlDoc->write( $outputFile ); 





