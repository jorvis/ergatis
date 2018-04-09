#!/usr/local/env perl

# diff_BSML_files.pl - perform a diff on two BSML documents, ignoring attributes
# that are generated with the Ergatis ID generator

# Shaun Adkins - sadkins@som.umaryland.edu

use strict;
use warnings;

# For BsmlParserTwig and BsmlDoc
use lib '/local/projects/ergatis/package-nightly/lib/perl5/';

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Data::Dumper;
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
profileBSMLFile($infile);
profileBSMLFile($infile2);
compareFiles($infile, $infile2);
close LOGFILE;

print "$0 program execution has completed\n";
print "Please review log file '$logfile'\n";
exit(0);

sub profileBSMLFile {
    my ($file) = @_;
    my $doc = BSML::BsmlDoc->new();
    my $parser = BSML::BsmlParserTwig->new();
    $parser->parse( \$doc, $file );
    print Dumper($doc);
    exit(0);
}

sub compareFiles {

    my ($file1, $file2) = @_;

    print LOGFILE "Comparing '$file1' versus '$file2'\n\n";


    my $profileHeadings = { 'organismCtr' => "<Organism>",
      'crossReferenceCtr' => "<Cross-reference>",
      'sequenceCtr' => "<Sequence> (loadable)",
      'nonLoadedSequenceCtr' => "<Sequence> (non-loadable)",
      'featureCtr' => "<Feature>",
      'attributeCtr' => "<Attribute>",
      'attributeListCtr' => "<Attribute-list>",
      'featureGroupCtr' => "<Feature-group>",
      'featureGroupMemberCtr' => "<Feature-group-member>",
      'uniqueFeatureGroupMemberCtr' => "<Feature-group-member> (unique)",
      'intervalLocCtr' => "<Interval-loc>",
      'analysisCtr' => "<Analysis>",
      'seqPairAlignmentCtr' => "<Seq-pair-alignment>",
      'seqPairRunCtr' => "<Seq-pair-run>",
      'linkCtr' => "<Link>",
      'multipleAlignmentTableCtr' => "<Multiple-alignment-table>"
  };

    my $profileLookupHeadings = { 'organismTypeCtr' => "<Organism> by genus-species",
    'sequenceClassCtr' => "<Sequence> (loadable) by class",
    'nonLoadedSequenceClassCtr' => "<Sequence> (non-loadable) by class",
    'featureClassCtr' => "<Feature> by class",
    'featureGroupMemberByClassCtr' => "<Feature-group-member> by class",
    'attributeTypeCtr' => "<Attribute> by name",
    'linkRelCtr' => "<Link> by rel",
    'linkHrefCtr' => "<Link> by href",
    'linkRoleCtr' => "<Link> by role"
        };

    foreach my $counterType (keys %{$profileHeadings} ){
      my $count1=$masterProfileLookup->{$file1}->{$counterType};
      my $count2=$masterProfileLookup->{$file2}->{$counterType};
      if ($count1 != $count2){
          print LOGFILE "$profileHeadings->{$counterType}\n".
            "$file1: $count1\n".
            "$file2: $count2\n\n";
      }
    }

    foreach my $lookupType (keys %{$profileLookupHeadings} ){
      ## Collect all of the different types of keys present in either file
      my $uniqueKeysLookup = {};
     	foreach my $key (keys %{$masterProfileLookup->{$file1}->{$lookupType}} ){
          $uniqueKeysLookup->{$key}++;
      }
      foreach my $key (keys %{$masterProfileLookup->{$file2}->{$lookupType}} ){
          $uniqueKeysLookup->{$key}++;
      }
      foreach my $key (keys %{$uniqueKeysLookup}) {
        my $count1=$masterProfileLookup->{$file1}->{$lookupType}->{$key};
        my $count2=$masterProfileLookup->{$file2}->{$lookupType}->{$key};
        if ($count1 != $count2){
        print LOGFILE "$profileLookupHeadings->{$lookupType} '$key'\n".
          "$file1: $count1\n".
          "$file2: $count2\n\n";
        }
      }
    }
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
