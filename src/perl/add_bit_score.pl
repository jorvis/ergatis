#!/usr/local/bin/perl

=head1  NAME 

add_bit_score.pl  - adds bit_score percentage to existing BLASTP Bsml documents

=head1 SYNOPSIS

USAGE:  add_bit_score.pl -f blastp.bsml -o blastp.bsml

=head1 OPTIONS

=over 4

=item *

B<--bsml_file,-f> [REQUIRED] The Bsml file containing pairwise alignments (i.e. BLASTP)

=item *

B<--output,-o> output BSML file containing bit_score information (default: input file)

=item *

B<--help,-h> This help message

=back

=head1   DESCRIPTION

add_bit_score.pl is designed to add bit_score percentage information to an 
existing BSML document containing pairwise alignment (i.e. BLASTP).  Bit_score
percentage is defined as the bit_score of gene A vs gene B divided by the
bit_score of gene A vs itself.

Samples:

1. add bit_score percentage to file "blastp.bsml"

   add_bit_score.pl -f blastp.bsml -o blastp.bsml 

   Note*  if -o is not specified, the output file will be same as input file.

NOTE:  

Calling the script name with NO flags/options or --help will display the syntax requirement.


=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BSML::BsmlReader;
use BSML::BsmlParserTwig;
use File::Basename;
use BSML::BsmlBuilder;
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 'bsml_file|f=s', 'bsml_dir|b=s', 'output|o=s', 'verbose|v', 'help|h', 'man') || pod2usage();


###-------------PROCESSING COMMAND LINE OPTIONS-------------###
my $bsml_file   = $options{'bsml_file'};
my $output_file = $options{'output'} || $bsml_file;
my $BSML_dir   = $options{'bsml_dir'};
$BSML_dir =~ s/\/+$//;         #remove terminating '/'s
my $verbose    = $options{'verbose'};

&cmd_check();
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

sub cmd_check {
#quality check

    if( exists($options{'man'})) {
	pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT});
    }   
    
    if( exists($options{'help'})) {
	pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT});
    }
    
    
    if(!$bsml_file or !$output_file) {
	pod2usage({-exitval => 2,  -message => "$0: All the required options are not specified", -verbose => 1, -output => \*STDERR});    
    }

    if(! -e $bsml_file) {
	print STDERR "$bsml_file CANNOT be found... Aborting!!!\n";
	exit 3;
    }
}
