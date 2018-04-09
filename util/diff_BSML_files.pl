#!/usr/local/env perl

# diff_BSML_files.pl - perform a diff on two BSML documents, ignoring attributes
# that are generated with the Ergatis ID generator

# Shaun Adkins - sadkins@som.umaryland.edu

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Array::Compare;
use Data::Dumper;
use List::Util qw(first);

# For BsmlParserTwig and BsmlDoc
use lib '/local/projects/ergatis/package-nightly/lib/perl5/';
use BSML::BsmlDoc;
use BSML::BsmlParserTwig;

## Command-line options
my ($debug, $infile, $infile2, $help, $logfile, $man, $asmbl_id);

my $results = GetOptions (
    'debug=s'    => \$debug,
    'help|h'     => \$help,
    'logfile=s'  => \$logfile,
    'man|m=s'    => \$man,
    'infile=s'   => \$infile,
    'infile2=s'  => \$infile2,
    'asmbl_id=s' => \$asmbl_id
    );

print_usage() if ($help);
pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if($man);

my $fatalCtr=0;

if (!$infile || !$infile2){
    print STDERR "--infile and --infile2 must be defined\n";
    $fatalCtr++;
} elsif (checkFile($infile) == 0 || checkFile($infile2) == 0){
    $fatalCtr++;
}

if ($fatalCtr>0){
    print_usage();
}

if (!defined($logfile)){
    $logfile = '/tmp/diff_BSML.log';
}

open(LOGFILE, ">$logfile") || die "Could not open log file '$logfile':$!";

my $masterProfileLookup={};
my $bsmlObj1 = profileBSMLFile($infile);
my $bsmlObj2 = profileBSMLFile($infile2);
compareFiles($bsmlObj1, $bsmlObj2);
close LOGFILE;

print "BSML files comparison ended in success.  Exiting.\n";
exit(0);

sub profileBSMLFile {
    my ($file) = @_;
    my $doc = BSML::BsmlDoc->new();
    my $parser = BSML::BsmlParserTwig->new();
    $parser->parse( \$doc, $file );
    # print Dumper($doc->{'BsmlSequences'});
    return $doc;
}

sub compareFiles {
    my ($bsmlObj1, $bsmlObj2) = @_;

    # Ensure all outermost elements are present
    # Convert blessed object into hash for comparisons;
    my $bsml_h1 = { %$bsmlObj1 };
    my $bsml_h2 = { %$bsmlObj2 };
    my $bsml_pass = are_keys_common($bsml_h1, $bsml_h2);
    die "BSML document elements are different between both files. Check logfile at $logfile\n" unless $bsml_pass;

    # BSMLAnalyses - SKIP
    # BsmlTableId - SKIP
    # GZIP - SKIP
    # UID_COUNTER - SKIP
    # BsmlMultipleAlignmentTables - SKIP
    # BsmlGenomes - SKIP
    # BsmlSegmentSets - SKIP
    # EXT_DTD - SKIP
    # BsmlAttr - SKIP
    # 'attr' - SKIP

    # BsmlSequences - SKIP if BsmlAttr->'defline' defined (source sequence)
    my $seq_pass compare_bsml_sequences($bsmlObj1->{'BsmlSequences'}, $bsmlObj2->{'BsmlSequences'});
    print "Some BSML sequences were different. Check logfile at $logfile\n" unless $seq_pass;

    # BsmlSeqPairAlignments
    my $spa_pass += compare_bsml_seq_pair_alignments($bsmlObj1->{'BsmlSeqPairAlignment'}, $bsmlObj2->{'BsmlSeqPairAlignment'});
    print "Some BSML seq-pair-alignments were different. Check logfile at $logfile\n" unless $spa_pass;

    die "BSML document elements are different between both files. Check logfile at $logfile\n" unless $seq_pass && $spa_pass;
}

sub compare_bsml_sequences {
    my ($seq_arr1, $seq_arr2) = @_;

    # Number of elements should match between both BSML files
    if (scalar @$seq_arr1 != scalar @$seq_arr2) {
      print LOGFILE "BSML <Sequences> tag has differing number of <Sequence> inner tags.\n";
      print LOGFILE "--infile1 : " . scalar @$seq_arr1 ."\n";
      print LOGFILE "--infile2 : " . scalar @$seq_arr2 ."\n";
      print LOGFILE "--- BSML <Sequences> from infile1 ---\n";
      print LOGFILE Dumper($seq_arr1);
      print LOGFILE "--- BSML <Sequences> from infile2 ---\n";
      print LOGFILE Dumper($seq_arr2);
      return 0;
    }

    # PASS if BsmlAttr->'defline' defined (source sequence)
    # BsmlSeqDataImport
    # BsmlAttributeList - SKIP
    # BsmlSeqData - SKIP
    # BsmlNumbering - SKIP
    # BsmlLink - SKIP
    # BsmlCrossReference
    # 'attr'
    # BsmlAttr
    # BsmlFeatureGroups
    # BsmlFeatureTables

    # Now that we know the array lengths are equal, let's see if each Sequence
    # pairs 1-to-1 with a Sequence on the other BSML file

    my $seq_tests = 1;

    my @seen_seq_objs = ();
    # First push other seq_obj with the defline
    foreach my $seq_obj2 (@$seq_arr2) {
      if exists $seq_obj2->{'BsmlAttr'}->{'defline'} {
        push @seen_seq_obj, $seq_obj2;
        last;
      }
    }

    foreach my $seq_obj1 (@$seq_arr1){
      next if exists $seq_obj1->{'BsmlAttr'}->{'defline'};
      foreach my $seq_obj2 (@$seq_arr2){
        # Skip if IDs don't match or
        next if first {$_ eq $seq_obj2} @seen_seq_objs;
        next if $seq_obj1->{'attr'}->{'id'} ne $seq_obj2->{'attr'}->{'id'};

        # We have found the right ID.  I'm going to trust the source database put the correct length and title in.
        if ($seq_obj1->{'BsmlSeqDataImport'}->{'identifier'} ne $seq_obj2->{'BsmlSeqDataImport'}->{'identifier'} ) {
          print LOGFILE "BSML <Sequences> have same <attr> tags but different <Seq-data-import> tags\n";
          print LOGFILE "--- BSML <Sequence> from infile1 ---\n";
          print LOGFILE Dumper($seq_obj1);
          print LOGFILE "--- BSML <Sequence> from infile2 ---\n";
          print LOGFILE Dumper($seq_obj2);
          $seq_tests = 0;
        }
        # Push if object is a 1-to-1 match with the first object
        push @seen_seq_obj $seq_obj2;
      }
    }

    # Verify all of the other Sequence objects have been seen.
    if (scalar @seen_seq_obj != scalar @seq_arr2) {
      print LOGFILE "There was not a 1-to-1 mapping of BSML <Sequence> tags.  Some tags in --infile1 and --infile2 did not map to each other.\n";
      print LOGFILE "Matches found in " . scalar @seen_seq_obj . " of " . scalar @seq_arr2 . " tags.\n";
      return 0;
    }

    return $seq_tests;

}

sub compare_bsml_seq_pair_alignments {
    my ($spa_arr1, $spa_arr2) = @_;

    # Number of elements should match between both BSML files
    if (scalar @$spa_arr1 != scalar @$spa_arr2) {
      print LOGFILE "BSML <Seq-pair-alignments> tag has differing number of <Seq-pair-alignment> inner tags.\n";
      print LOGFILE "--infile1 : " . scalar @$spa_arr1 ."\n";
      print LOGFILE "--infile2 : " . scalar @$spa_arr2 ."\n";
      print LOGFILE "--- BSML <Seq-pair-alignments> from infile1 ---\n";
      print LOGFILE Dumper($spa_arr1);
      print LOGFILE "--- BSML <Seq-pair-alignments> from infile2 ---\n";
      print LOGFILE Dumper($spa_arr2);
      return 0;
    }

    # BsmlLink - SKIP
    # BsmlSeqPairRuns
    # 'attr'
    # BsmlAttr
    # BsmlCrossReference
    #my $spa_obj1;
    #my $spa_obj2;
    #my $spr_test = compare_bsml_seq_pair_run($spa_obj1->{'BsmlSeqPairRun'}, $spa_obj2->{'BsmlSeqPairRun'});
    return 1
}

sub compare_bsml_seq_pair_run {
    my ($spr_arr1, $spr_arr2) = @_;

    # Number of elements should match between both BSML files
    if (scalar @$spr_arr1 != scalar @$spr_arr2) {
      print LOGFILE "BSML <Seq-pair-runs> tag has differing number of <Seq-pair-run> inner tags.\n";
      print LOGFILE "--infile1 : " . scalar @$spr_arr1 ."\n";
      print LOGFILE "--infile2 : " . scalar @$spr_arr2 ."\n";
      print LOGFILE "--- BSML <Seq-pair-runs> from infile1 ---\n";
      print LOGFILE Dumper($spr_arr1);
      print LOGFILE "--- BSML <Seq-pair-runs> from infile2 ---\n";
      print LOGFILE Dumper($spr_arr2);
      return 0;
    }
    # BsmlLink - SKIP
    # BsmlAttr
    # 'attr'
    return 1
}

# Check to see if two hashes have the exact same keys
sub are_keys_common {
    my ($bsml_h1, $bsml_h2) = @_;

    my $comp = Array::Compare->new();
    my @bsml_arr1 = keys %{$bsml_h1};
    my @bsml_arr2 = keys %{$bsml_h2};
    return 1 if $comp->perm(\@bsml_arr1, \@bsml_arr2);

    # Hashes are not the same, so fail
    print LOGFILE "Both sets of BSML have differing keys\n";
    print LOGFILE "--- BSML from infile1 ---\n";
    print LOGFILE Dumper($bsml_h1);
    print LOGFILE "--- BSML from infile2 ---\n";
    print LOGFILE Dumper($bsml_h2);
    return 0
}

sub checkFile {

    my ($file) = @_;
    my $status=1;
    if (!-e $file){
        print "input BSML file '$file' does not exist\n";
        $status=0;
    }
    if (!-f $file){
        print "input BSML file '$file' is not a file\n";
        $status=0;
    }
    if (!-r $file){
        print "input BSML file '$file' does not have read permissions\n";
        $status=0;
    }
    if (!-s $file){
        print "input BSML file '$file' has no content\n";
        $status=0;
    }
    return $status;
}

sub print_usage {
    print "Sample usage: perl $0 --infile [--infile2] [--logfile]\n";
    exit(1);
}
