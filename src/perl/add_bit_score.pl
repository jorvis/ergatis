#!/usr/local/bin/perl


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BSML::BsmlReader;
use BSML::BsmlParserTwig;
use File::Basename;
use BSML::BsmlBuilder;


my %options = ();
my $results = GetOptions (\%options, 'bsml_file|f=s', 'bsml_dir|b=s', 'output|o=s', 'verbose|v', 'help|h',);


###-------------PROCESSING COMMAND LINE OPTIONS-------------###
my $bsml_file   = $options{'bsml_file'};
my $output_file = $options{'output'} || $bsml_file;
my $BSML_dir   = $options{'bsml_dir'};
$BSML_dir =~ s/\/+$//;         #remove terminating '/'s
my $verbose    = $options{'verbose'};

if(!$bsml_file or !$output_file  or exists($options{'help'})) {
    #$logger->fatal("Not all of the required options have been defined.  Exiting...");
    &print_usage();
}

###-------------------------------------------------------###

my $reader = BSML::BsmlReader->new();
my $parser = BSML::BsmlParserTwig->new();

$parser->parse( \$reader, $bsml_file );

foreach my $seqobj (@{$reader->returnAllSequences()}) {
    my $seqId = $seqobj->returnattr('id');
    # return the alignment object corresponding to the query sequences verses itself.
    my $alignment_list = $reader->fetch_all_alignmentPairs( $seqId, $seqId );
    if( $alignment_list->[0] ) {
	my $base_alignment = $reader->readSeqPairAlignment( $alignment_list->[0] );
	my $base_score = $base_alignment->{'seqPairRuns'}->[0]->{'runscore'};
	# get all the alignment objects having $seqId as the query sequence
	$alignment_list = $reader->fetch_all_alignmentPairs( $seqId );
	# foreach alignment, calculate the percent bit score for each seq pair run
	foreach my $align (@{$alignment_list}) {
	    my $alignmentDat = $reader->readSeqPairAlignment( $align );
	    foreach my $run ( @{$alignmentDat->{'seqPairRuns'}} ) {
		my $score = $run->{'runscore'};
		my $PercentBitScore = sprintf('%.3f', $score/$base_score);
		# add a bsml attribute to the document for percent bit score
		print STDERR "adding percent bit score $PercentBitScore\n" if($verbose);
		$run->{'bsmlRef'}->addBsmlAttr( 'percent_bit_score', $PercentBitScore );
	    }
	}
    } else {
	print STDERR "NO BASE MATCH found\n" if($verbose);
    }
}

# write the altered file to disk
$reader->write( $output_file );
chmod 0666, $output_file;

sub print_usage {

    print STDERR "SAMPLE USAGE:  add_bit_score.pl -f bsml_file -o output_file\n";
    print STDERR "  --bsml_file    = BSML doc containing blastp info\n";
    print STDERR "  --output       = dir to save output to (if omitted, overwrite input file)\n";
    print STDERR "  --help = This help message.\n";
    exit 1;

}
