#!/usr/bin/perl
# Input a list of BSML files
# perform a series of QC checks on the files
# output summary results
use strict;
use warnings;
use File::Basename;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use XML::Twig;
use Data::Dumper;

my %opts = &parse_options();
my $input_list = $opts{input_list};
my $odir = $opts{output_dir};

# create list of input files
my @input_files = compile_input_files();


my $LOG_BSML; # for purposes of the logger, the current BSML file
my $OUT_BASE; # because we want a single output file for each execution of the script (not for each input file)
              # but when run as a component it might utilize the same output directory,
              # so use the name of the input_file or input_list to individually identify the output data

my %qc; # generic hash to store all qc count info
my $oname; # name of currently processing organism

my %file_qc;  # like above, but counts on a per file/assembly basis
my $db;
my $asmbl_id;

# create twig
my $twig= new XML::Twig( 
			twig_handlers =>                 
			{ 
 			 'Organism' => \&process_Organism,
 			 'Feature-group' => \&process_Featuregroup,
 			 'Feature' => \&process_Feature,
 			 'Seq-data-import' => \&process_Seqdataimport,
 			 'Sequence' => \&process_Sequence,
 			 'Genome' => \&process_Genome_Crossreference,


#			  "Genome/Cross-reference/Attribute[name=source_database]" => \&process_source_db,
#			 'Genome/Cross-reference/Attribute' => \&process_source_db,
#			 'Genome/Cross-reference' => \&process_Genome_Crossreference,

			}
                       );


open(my $LOG, ">$odir/bsml_qc.$OUT_BASE.log") || die "Unable to write to log ($odir/bsml_qc.log):$!";
open(my $AC,  ">$odir/bsml_qc.$OUT_BASE.annotation_counts.txt") || die "Unable to write to log ($odir/bsml_qc.annotation_counts.txt):$!";
# process each bsml file
foreach my $bsmlfile (@input_files) {
  $LOG_BSML = $bsmlfile;
  # TODO: what if the file is .gz?

  # Deal w/ name junk
  $db = undef;
  $asmbl_id = undef;
  my $basename = basename($bsmlfile);
  ($db, $asmbl_id, my $other) = split("_", $basename);

  # if neither were defined then assume we're opering on a genbank2bsml output file
  unless ( defined($db) && defined($asmbl_id) ) {
    if ( $basename =~ /^([^\.]+)\./ ) {
      $db = $1;
      $asmbl_id = $db;
    }
    else { warn "going to die"; }
  }
  die "No db in basename $basename" unless (defined $db);
  die "No asmbl_id in basename $basename of $bsmlfile" unless (defined $asmbl_id);

  die "Dup db/asmbl_id ($db / $asmbl_id)" if (exists($file_qc{$db}->{$asmbl_id}));

  # reset the organism name to 'unknown' in case it's missing
  $oname = 'UNKNOWN';
  print "Parsing $bsmlfile\n";
  
  my $IFH;
  if ($bsmlfile =~ /\.(gz|gzip)$/) {
    open ($IFH, "<:gzip", $bsmlfile) || die "couldn't open '$bsmlfile' for reading: $!";
  } else {
    open ($IFH, "<".$bsmlfile) || die "couldn't open '$bsmlfile' for reading: $!";
  }
  $twig->parse($IFH);

  # print out per-file counts to match Jay's counts
  unless (exists $file_qc{$db}->{$asmbl_id}->{transcript}) { $file_qc{$db}->{$asmbl_id}->{transcript} = 0; }
  unless (exists $file_qc{$db}->{$asmbl_id}->{polypeptide}) { $file_qc{$db}->{$asmbl_id}->{polypeptide} = 0; }
  print {$AC} "$db\t$asmbl_id\t".$file_qc{$db}->{$asmbl_id}->{transcript}."\t".$file_qc{$db}->{$asmbl_id}->{polypeptide}."\n";

  $twig->purge();
  close($IFH);

}
#close($FIN);
close($LOG);
close($AC);

open(my $SOUT, ">$odir/summary.out") || die "Unable to write to $odir/summary.out: $!";
# output count info for each organism
print "number of organisms: ".scalar(keys %qc)."\n";
foreach (keys %qc) {
  print_counts($_);
}
close($SOUT);
print "Done w/ script!\n";


#
# SUBROUTINES
#

sub log_say {
  my $message = shift;

  print {$LOG} "$LOG_BSML: $message\n";

}

sub check_att {
  my ($tag, $att) = @_;

  my $name = $tag->name;

  my $value = $tag->{'att'}->{$att};
  unless (defined($value)) {
    log_say("Missing attribute ($att) in element ($name)");
    return undef;
  }

  if ($value eq '') {
    log_say("Empty string for attribute ($att) in element ($name)");
    return undef;
  }

  return $value;
}

# return all input files as an array
sub compile_input_files {

  my @input_files = ();
  
  if ( defined($opts{input_list}) ) {
    open (my $FIN, $opts{input_list}) || die "Unable to open input_list ($opts{input_list}):$!";
    while (my $bsmlfile = <$FIN>) {
      chomp($bsmlfile);
      push (@input_files, check_file($bsmlfile));
      $OUT_BASE = basename($opts{input_list});
    }
  }
  elsif ( defined($opts{input_file}) ) {
    push (@input_files, check_file($opts{input_file}));
    $OUT_BASE = basename($opts{input_file});
  } 

  return @input_files;
}

# check if file exists, or as a zip
sub check_file {
  my $file = shift;

  return $file if (-r $file);

  $file = $file.".gz";

  return $file if (-r $file);

  die "Unable to read file ($file)";

}



sub print_counts {
  my $oname = shift;

  print {$SOUT} "Organism\t$oname\n";

  print {$SOUT} "Files\ttotal\t".$qc{$oname}->{file_count}->{total}."\n";
  
  print {$SOUT} "Feature-groups\ttotal\t" .$qc{$oname}->{tag_count}->{'Feature-group'}."\n";
  print {$SOUT} "Feature-groups\tvalid\t" .$qc{$oname}->{valid_count}->{'Feature-group'}."\n";
  print {$SOUT} "Feature-groups\tpartial\t" .$qc{$oname}->{partial_count}->{'Feature-group'}."\n";
  print {$SOUT} "Feature-groups\tdup polypeptide\t" .$qc{$oname}->{dup_polypeptide_count}->{'Feature-group'}."\n";

  foreach my $class (sort keys %{$qc{$oname}->{Featureclass}}) {
    print {$SOUT} "Feature class\t$class\t".$qc{$oname}->{Featureclass}->{$class}."\n";
  }

  foreach my $key (sort keys  %{$qc{$oname}->{transcript_info_count}}) {
    print {$SOUT} "Transcript info\t$key\t".$qc{$oname}->{transcript_info_count}->{$key}."\n";
  }

  print {$SOUT} "Feature/Cross-reference\ttotal\t" .$qc{$oname}->{tag_count}->{'Feature/Cross-reference'}."\n";
  print {$SOUT} "Feature/Cross-reference\tvalid\t" .$qc{$oname}->{valid_count}->{'Feature/Cross-reference'}."\n";

  print {$SOUT} "Seq-data-imports\ttotal\t" .$qc{$oname}->{tag_count}->{'Seq-data-import'}."\n";
  print {$SOUT} "Seq-data-imports\texist\t" .$qc{$oname}->{valid_count}->{'Seq-data-import'}."\n";

  print {$SOUT} "Sequence\ttotal\t" .$qc{$oname}->{tag_count}->{'Sequence'}."\n";
  print {$SOUT} "Sequence\tvalid\t" .$qc{$oname}->{valid_count}->{'Sequence'}."\n";
  foreach my $key (sort keys %{$qc{$oname}->{Sequence_info}}) {
    print {$SOUT} "Sequence\tw/ $key\t" .$qc{$oname}->{Sequence_info}->{$key}."\n";
  }
  print {$SOUT} "source_database\t" .$qc{$oname}->{valid_count}->{'Genome/Cross-reference/Attribute'}."\n";

}


# <Attribute name="source_database" content="bcl"></Attribute> in
# <Cross-reference database="TIGR_Bcl" identifier="c_perfringens atcc13124" id="_1" identifier-type="legacy_annotation_database">
# in <Genome> 
sub process_Genome_Crossreference {
  my( $twig, $tag)= @_;
  my $source_database = 0;

  foreach my $attribute ($tag->get_xpath('//Cross-reference/Attribute')) {
    (my $att_name = check_att($attribute, 'name')) || next;
  
    # content could == 0
    my $att_content = check_att($attribute, 'content');
    next unless defined($att_content);

    if ($att_name eq 'source_database') {
      ++$source_database;
    }
  }

  if ($source_database) {
    ++$qc{$oname}->{valid_count}->{'Genome/Cross-reference/Attribute'}
  }
  else {
    log_say("No Genome/Cross-reference/Attribute for source_database");
  }
}

# note: only care about assemblies!
# a valid sequence has a topology, molecule_name, and molecule_type
sub process_Sequence {
  my( $twig, $tag)= @_;
  # if they're new, they'll be > 1
  # and if a required was missing it'll still = 0
  my %sequence_info = (
		       topology => 0, # <Sequence length="14" topology="linear">
		       molecule_name => 0, # <Attribute name="molecule_name" >
		       molecule_type => 0  # <Attribute name="molecule_type">
		      );
  
  (my $class = check_att($tag, 'class')) || return;
  return unless ($class eq 'assembly');

  # count it as a sequence we saw (really only counting assembly sequences)
  ++$qc{$oname}->{tag_count}->{'Sequence'};

  (my $id = check_att($tag, 'id')) || return;

  # even if topology is absent keep going
  (my $topology = check_att($tag, 'topology')) && ++$sequence_info{topology};

  foreach my $attribute ($tag->get_xpath('Attribute')) {
    (my $att_name = check_att($attribute, 'name')) || next;
    
    # content could == 0
    my $att_content = check_att($attribute, 'content');
    next unless defined($att_content);

    ++$sequence_info{$att_name};
  }

  my $valid_sequence = 1;
  foreach my $key (keys %sequence_info) {
    if ($sequence_info{$key} == 0) {
      #++$qc{$oname}->{sequence_info_count}->{$key};
      log_say("Sequence ($id) missing $key");
#      return;
      $valid_sequence = 0;
    }
    else {
      ++$qc{$oname}->{'Sequence_info'}->{$key};
    }
  }
  if ($valid_sequence) {
    ++$qc{$oname}->{valid_count}->{'Sequence'};
  }
}

# count each Feature class
# note that the same Feature can appear in > 1 Feature-groups!
sub process_Feature {
  my( $twig, $tag)= @_;
  
  (my $class = check_att($tag, 'class')) || return;
  ++$qc{$oname}->{Featureclass}->{$class};
  ++$file_qc{$db}->{$asmbl_id}->{$class};

  # All Feature/Cross-reference elements need to have a database and identifier values
  foreach my $xref ($tag->get_xpath('Cross-reference')) {
    ++$qc{$oname}->{tag_count}->{'Feature/Cross-reference'};
    my $xref_database = check_att($xref, 'database');
    my $xref_identifier = check_att($xref, 'identifier');
    if (defined($xref_database) && defined($xref_identifier)) {
      # valid Feature/Cross-reference
      ++$qc{$oname}->{valid_count}->{'Feature/Cross-reference'};
    }
    else {
#      die "Missing something";
    }
  }

  if ($class eq 'transcript') {
    process_Feature_transcript($tag);
  }
}

# transcript has a locus and gene_product_name 
sub process_Feature_transcript {
  my $tag = shift;

  my %transcript_info = (
			  locus => 0, #  <Cross-reference  identifier-type="locus">
			  gene_product_name => 0  #  <Attribute name="gene_product_name">
			 );

  my $id = $tag->{'att'}->{'id'};

  foreach my $attribute ($tag->get_xpath('Attribute')) {
      (my $att_name = check_att($attribute, 'name')) || next;

      # content could == 0
      my $att_content = check_att($attribute, 'content');
      next unless defined($att_content);

      if ($att_name eq 'gene_product_name') {
	++$transcript_info{gene_product_name};
      }
  }

  foreach my $xref ($tag->get_xpath('Cross-reference')) {
    # this is optional atm so no log
#    (my $xref_type = check_att($xref, 'identifier-type')) || next;
    my $xref_type = $xref->{'att'}->{'identifier-type'};
    next unless defined($xref_type);

    (my $xref_identifier = check_att($xref, 'identifier')) || next;

    if ($xref_type eq 'locus') {
      ++$transcript_info{locus};
    }
  }
  
  foreach my $key (keys %transcript_info) {
    if ($transcript_info{$key} > 0) {
      ++$qc{$oname}->{transcript_info_count}->{$key};
    }
    else {
      log_say("transcript ($id) missing $key");
    }
  }
}

# Handle Feature-group:
# - gene/transcript/CDS/exon/polypeptide present?
sub process_Featuregroup {
    my( $twig, $tag)= @_;
    # feature types to look for in Feature-group
    # if they're new, they'll be > 1
    # and if a required was missing it'll still = 0
    my %featuretypes = (
			gene => 0,
			transcript => 0,
			CDS => 0,
			exon => 0,
			polypeptide => 0
		       );

    ++$qc{$oname}->{tag_count}->{'Feature-group'};

    my $group_set = $tag->{'att'}->{'group-set'};

    foreach my $fgm ($tag->get_xpath('Feature-group-member')) {
      (my $featuretype = check_att($fgm, 'feature-type')) || next;

      ++$featuretypes{$featuretype};
      ++$qc{$oname}->{tag_count}->{'Feature-group-member'};
    }

    # if all the feature model counts are zero, don't do anything
    my $count_sum = 0;
    foreach ('gene', 'transcript', 'CDS', 'exon', 'polypeptide') {
      $count_sum += $featuretypes{$_};
    }
    if ($count_sum == 0) {
      return;
    }

    # otherwise if one is zero, then it's a partial
    foreach (keys %featuretypes) {
      if ($featuretypes{$_} == 0) {
	log_say("Partial Feature-group ($group_set)");
	++$qc{$oname}->{partial_count}->{'Feature-group'};
	return;
      }
    }

    # otherwise if there are more polypeptides than transcripts it's extra
    if ( $featuretypes{polypeptide} > 1 ) {
      log_say("Duplicate polypeptides ($featuretypes{polypeptide}) in ($group_set)");
      	++$qc{$oname}->{dup_polypeptide_count}->{'Feature-group'};
	return;
      }

    # otherwise, it's valid
    ++$qc{$oname}->{valid_count}->{'Feature-group'};
}

# Seq-data-import files actually exist 
sub process_Seqdataimport {
  my( $twig, $tag)= @_;
  ++$qc{$oname}->{tag_count}->{'Seq-data-import'};

  (my $source = check_att($tag, 'source')) || next;
  (my $identifier = check_att($tag, 'identifier')) || next;

  if ( (-e $source) || (-e "$source.gz")) {
    ++$qc{$oname}->{valid_count}->{'Seq-data-import'};
  }
  else {
    log_say("Nonexistant source ($source) in Seq-data-import ($identifier)");
    return; # counting valid so do nothing!
  }
}

# store the genus+species as the current organism name
sub process_Organism {
  my( $twig, $tag)= @_;

  (my $genus = check_att($tag, 'genus')) || die "No Organism genus";
#  my $genus = $tag->{'att'}->{'genus'};
#  die "$bsmlfile: Missing genus attribute for Organism element" unless defined($genus);
  
  (my $species = check_att($tag, 'species')) || die "No Organism species";
#  my $species = $tag->{'att'}->{'species'};
#  die "$bsmlfile: Missing species attribute for Organism element" unless defined($species);
  
  $oname = "$genus $species";

  unless (defined($qc{$oname})) {
    new_oname($oname);
  }

  ++$qc{$oname}->{file_count}->{total};
}

# zero counts (note: not zeroing everying we count, just stuff that throws concat warns)
sub new_oname {
  my $new_name = shift;
  
  $qc{$oname} = {
		 file_count => {
				total => 0,
				invalid => 0
			       },
		 valid_count => {
				 'Feature-group' => 0,
				 'Seq-data-import' => 0,
				 'Sequence' => 0,
				 'Genome/Cross-reference/Attribute' => 0,
				 'Feature/Cross-reference' => 0,
				 },
		 tag_count => {
			       'Feature-group' => 0,
			       'Seq-data-import' => 0,
			       'Sequence' => 0,
			       'Feature/Cross-reference' => 0,
				 },
		 partial_count => {
				 'Feature-group' => 0,
				 },
		 dup_polypeptide_count => {
				 'Feature-group' => 0,
				 },
		 };

}


sub parse_options {
    my %options = ();

    GetOptions( \%options,
                'input_list|i=s',
                'input_file|i=s',
                'output_dir|o=s',
		'help|h',
                ) || print_instructions("Unprocessable option");

    if ($options{help}) {
      print_instructions('');
    }

    unless ( (defined($options{input_list})) || (defined($options{input_file})) ) {
        print_instructions( "input_list ($options{input_list}) not readable");
      }

    (-d $options{output_dir}) 
        || print_instructions( "output_dir ($options{output_dir}) not a directory");

#    print "Executing $0 with options\n";
#    foreach (sort keys %options) { print "  $_: $options{$_}\n";}
    
    return %options;
}


sub print_instructions {
  my $message = shift;

  print qq{$message
BSML QC: Counts and Format
Executed as:
bsml_qc.pl --input_list|i <list file> --output_dir|o <dir>

Output directory will have
bsml_qc.log - lists all failed Quality Control Checks encountered on a file by file basis
summary.out - summary counts for each different Organism encountered in all files

Quality Control Checks:
- Checks that each assemble id is unique to a database.
- Checks that Organism name and species are provided
- Logs missing Cross-reference/Attribute for Genome
- Logs that each Sequence assembly has a topology, molecule_name, and molecule_type Attribute
- Logs that each Feature transcript has a locus Cross-reference and gene_product_name Attribute
- Logs that each Feature-group has a set of gene, transcript, CDS, exon, and polypeptide Features
- Confirm files referenced in Seq-data-import elements are present on filesystem and that they contain a source and identifier
- Log that Genome Cross-references, Sequence Attributes, and Feature transcript Cross-references and Attributes have names and content
- Log that each Sequence has a class 

Organism Counts:
Feature groups (valid/partial/duplicate polypeptides/total)
Number of Features of each class
Transcript info
Seq-data-imports (valid/total)
Source files
};

exit;

}
