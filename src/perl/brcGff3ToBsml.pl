#!/usr/local/bin/perl -w
#-----------------------------------------------------------------------
# program:   brcGff3ToBsml.pl
# author:    Jay Sundaram
# date:      2005-11-25
# 
# purpose:   Parses IOWG GFF and produces BSML document
#
#
#-------------------------------------------------------------------------


=head1 NAME

brcGff3ToBsml.pl - Parse IOWG GFF file and produces BSML document

=head1 SYNOPSIS

USAGE:  brcGff3ToBsml.pl [-d debug_level] [-f fastadir] -g gff_file [-h] [--id_generator_mapping_file] [--id_repository] [-l log4perl] [-m] [--no_id_generator] [--no_report_stats] [-o outdir]

=head1 OPTIONS

=over 8

=item B<--debug_level,-d>

 Optional: Ergatis::Logger log4perl logging level.  Default is 0

=item B<--fastadir,-f>

 Optional: Directory to which multi-fasta file containing all sequences will be written.  Default is current working directory.

=item B<--gff_file,-g>

 IOWG GFF3 file to be parsed

=item B<--id_generator_mapping_file>

 Optional: File from which old-identifier to new-identifier mappings will be read and written to.  The default file will be {$outdir}/brcGff3ToBsml.pl.{$gff_file}.mapping_file.txt

=item B<--id_repository>

 Optional: IdGenerator compliant directory (must contain valid_id_repository file).  Default is current working directory

=item B<--man,-m>

 Display the pod2usage page for this utility

=item B<--no_id_generator>

 Optional: If specified ID values present in the GFF file will be propagated into the resulting BSML file.  The default behaviour is to
           store identifiers produced by IdGenerator in the BSML file.

=item B<--no_report_stats>

 Optional: If specified will disable the default behaviour (will not produce stats file)  The stats file will contain a count of GFF3 records by type as well as some BSML element counts

=item B<--outdir,-o>

 Optional: Directory to which .bsml file will be written.  Default is current directory

=item B<--project>

 Optional: project argument for IdGenerator->next_id().  Default is 'brc'.

=item B<--help,-h>

 Print this help

=back

=head1 DESCRIPTION

brcGff3ToBsml.pl - Parse IOWG GFF file and produces BSML document

 Assumptions:
1. The BSML pairwise alignment encoding should validate against the XML schema:.
2. User has appropriate permissions (to execute script, access chado database, write to output directory).
3. Target chado database already contains all reference features (necessary to build feature and organism lookups) Review and execute db2bsml.pl if required.
4. Target chado database contains the necessary controlled vocabulary terms: "match" etc.
5. All software has been properly installed, all required libraries are accessible.

Sample usage:
./brcGff3ToBsml.pl -g /usr/local/annotation/OMNIUM/brc-central/GFF_BACKUP/TIGR/bacillus_anthracis/b_anthracis_a0071_western_north_america.gff


=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Data::Dumper;
use Digest::MD5 qw(md5);
use BSML::BsmlBuilder;
use Ergatis::Logger;
use Ergatis::IdGenerator;
use Config::IniFiles;
use File::Basename;
use URI::Escape;

$|=1;

#-------------------------------------------------------------
# Parse command line options
#-------------------------------------------------------------

my ($gff_file, $debug_level, $help, $log4perl, $man, $outdir, $fastadir, $no_report_stats, $no_id_generator, $id_generator_mapping_file, $id_repository, $project);


my $results = GetOptions (
			  'log4perl|l=s'        => \$log4perl,
			  'debug_level|d=s'     => \$debug_level, 
			  'help|h'              => \$help,
			  'man|m'               => \$man,
			  'outdir|o=s'          => \$outdir,
			  'gff_file|g=s'        => \$gff_file,
			  'fastadir|f=s'        => \$fastadir,
			  'no_report_stats=s'   => \$no_report_stats,
			  'no_id_generator=s'   => \$no_id_generator,
			  'id_generator_mapping_file=s' => \$id_generator_mapping_file,
			  'id_repository=s'             => \$id_repository,
			  'project=s'                   => \$project
			  );


&pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if ($man);
&pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT}) if ($help);



if (!$gff_file){
    print STDERR ("gff_file was not defined\n");
    &print_usage();
} 


#
# Get the Log4erl logger
#
my $logger = &get_logger($log4perl, $debug_level);


#
# verify and set the output directory
#
$outdir = &verify_and_set_outdir($outdir);


#
# verify and set the fasta output directory
#
$fastadir = &verify_and_set_fastadir($fastadir);

#
# Check permissions of the bsml document
#
&is_file_readable($gff_file);

$logger->info("Processing GFF file '$gff_file'");

#
# File handle for the output fasta file
#
my $fastafh;

my $basename = File::Basename::basename($gff_file);

$project = 'brc' if (!defined($project));

my $id_mapping_lookup = {};

my $idcreator;

if ((defined($no_id_generator)) && ($no_id_generator == 1)) {
    # user specified not to use IdGenerator
    $logger->warn("user specified not to use IdGenerator");
}
else {
    
    if (!defined($id_repository)){
	$id_repository = $outdir;
	$logger->warn("id_repository was not defined and was therefore set to outdir '$outdir'");
    }

    $idcreator = Ergatis::IdGenerator->new('id_repository' => $id_repository);

    $id_generator_mapping_file = &get_id_mapping_lookup($gff_file, $outdir, $id_generator_mapping_file, $id_mapping_lookup);
}


#
# Hash contain known gff3 attribute types
#
my $known_attributes =  { ID                => 1,
			  Dbxref            => 1,
			  Name              => 1,
			  description       => 1,
			  gene_symbol       => 1,
			  locus             => 1,
			  Parent            => 1,
			  localization      => 1,
			  topology          => 1,
			  molecule_type     => 1,
			  organism_name     => 1,
			  strain            => 1,
			  translation_table => 1,
			  ec_number         => 1,
			  Ontology_term     => 1,
			  Alias             => 1,
			  size              => 1,
			  stable_id         => 1,
			  translation_id    => 1
		      };

## Keep track of all contigs and supercontigs
my $primarySequences;

## The prefix is necessary for Cross-reference object creation
my $prefix;

## gffrecords => all GFF records
##
my ($gffrecords, $parent_lookup, $cdsSegmentLookup, $organism_attributes_lookup) = &parseGffFile();

## To ensure ID uniqueness
my $id_lookup = {};

## Lookup for all stats related information
my $report_stats = {};

## This will allow us to link the Sequence objects to the Genome objects
my $seqid_to_genome_lookup = {};

&mergeCDSSegments($gffrecords, $cdsSegmentLookup, $parent_lookup);

&derivePolypeptideFeatureSequenceForCds($gffrecords, $parent_lookup);

&createTranscriptFeatureStubs($gffrecords, $parent_lookup);


## Instantiate BsmlBuilder object
my $bsmlBuilder = new BSML::BsmlBuilder();

## Set the output BSML document name
$bsmlBuilder->{'doc_name'} = &create_out_file($outdir, $basename, 'bsml');

## Increment the Cross-reference 'id'
$bsmlBuilder->{'xrefctr'}++;

my $organismDataProcessed=0;

## This will contain the organism name to genome_id key-value pairs
my $genome_id_lookup = {};

## Maintain references for <Sequence> element objects
my $sequenceIdToBsmlSequenceLookup = {};

## Maintain references for <Feature-table> element objects
my $sequenceIdToBsmlFeatureTableLookup = {};

## Create the output fasta file name and open in write mode
my $fastafile = &create_out_file($fastadir, $basename, 'fsa');

open ($fastafh, ">$fastafile") || $logger->logdie("Could not open fastafile '$fastafile': $!");

my $translation_table = &storePrimarySequences($primarySequences,
					       $gffrecords,
					       $bsmlBuilder,
					       $genome_id_lookup,
					       $seqid_to_genome_lookup,
					       $sequenceIdToBsmlSequenceLookup,
					       $sequenceIdToBsmlFeatureTableLookup,
					       $organism_attributes_lookup);
#
# Create Feature BSML element objects for all Features.
# Also create all related BSML sub-elements (i.e. Interval-loc, Attribute,
# Attribute-list, Cross-reference, etc.)
#
&storeGffRecordsInBsml($primarySequences,
		       $gffrecords, 
		       $bsmlBuilder,
		       $genome_id_lookup,
		       $sequenceIdToBsmlFeatureTableLookup,
		       $sequenceIdToBsmlSequenceLookup,
		       $parent_lookup);

&storeFeatureAsFeatureGroupMember($parent_lookup,
				  $bsmlBuilder, 
				  $sequenceIdToBsmlSequenceLookup,
				  $gffrecords);

#
# Write the BSML document tree to the named BSML document
#
&write_out_bsml_doc($outdir, $bsmlBuilder, $gff_file, $no_report_stats);		       

$logger->info("BSML document created.");

if ((defined($no_id_generator)) && ($no_id_generator == 1)) {
    ## User specified not to use IdGenerator.  Since no
    ## identifiers were generated, there is no need to 
    ## write some mapping file.
}
else {
    &write_mapping_file($id_generator_mapping_file, $id_mapping_lookup);
}

print "End of program $0\nSee log file '$log4perl'\n";
exit(0);

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
#                                  END OF MAIN  -- SUBROUTINES FOLLOW
#
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------
# write_fasta_to_file()
#
#---------------------------------------------------------------
sub write_fasta_to_file {

    my ($seqid, $sequence, $fastafh) = @_;
    
    my $fastaout = &fasta_out($seqid, $sequence);
    
    print $fastafh "$fastaout";
    

}#end sub create_multifasta()


#-------------------------------------------------------------------------
# fasta_out()
#
#-------------------------------------------------------------------------
sub fasta_out {

    #This subroutine takes a sequence name and its sequence and
    #outputs a correctly formatted single fasta entry (including newlines).

    my ($seq_name, $seq) = @_;
    
    my $fasta=">"."$seq_name"."\n";

    $seq =~ s/\s+//g;

    for(my $i=0; $i < length($seq); $i+=60){

	my $seq_fragment = substr($seq, $i, 60);

	$fasta .= "$seq_fragment"."\n";
    }

    return $fasta;

}


#-------------------------------------------------------------------------------------------------------------------------------------------------------------
# create_genome_component()
#
# Sample <Genomes> component:
# <Genomes>
#   <Genome>
#      <Cross-reference id="1" database="TIGR_euk:tba1" identifier="t_brucei" identifier-type="current"></Cross-reference>
#      <Organism species="brucei" genus="Trypanosoma"></Organism>
#   </Genome>
# </Genomes>
#
#
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
sub create_genome_component {    
    
    my ($bsmlBuilder, $organism_name, $translation_table, $attributes, $organism_attributes_lookup, $source) = @_;

    if (!defined($source)){
	$logger->logdie("source was not defined");
    }
    
    my $genus;
    my @straininfo; # capture the species and strain

    ($genus, @straininfo) = split(/\s+/, $organism_name);
    
    my $species = join(" ", @straininfo);
    
    my $identifier = lc(substr($genus,0,1)) . '_' . lc($species);
    

    #
    # The strain value is incorporated into the //Cross-reference/@identifier
    #
    if (( exists $organism_attributes_lookup->{'strain'})  && (defined($organism_attributes_lookup->{'strain'}))) {

	my $strain = $organism_attributes_lookup->{'strain'};
    
	$identifier .= "_" . lc($strain);
	
	#
	# The stain value is incorporated into the //Genome/@species
	#
	$species .= " $strain";
	
	## undefine the strain hashkey-value here
	delete $organism_attributes_lookup->{'strain'};
    }

    my $genome_elem = $bsmlBuilder->createAndAddGenome();
    
    $logger->logdie("Could not create <Genome> element object reference for GFF3 file '$gff_file'") if (!defined($genome_elem));

    my $xref_elem = $bsmlBuilder->createAndAddCrossReference(
						     'parent'          => $genome_elem,
						     'id'              => $bsmlBuilder->{'xrefctr'}++,
						     'database'        => $source,
						     'identifier'      => $identifier,
						     'identifier-type' => 'brc_name'
						     );


    $logger->logdie("Could not create <Cross-reference> element object reference for GFF3 file '$gff_file'") if (!defined($xref_elem));
    

    my $xref_identifier_type=1;

    &store_genome_cross_reference($organism_attributes_lookup,
				  $genome_elem,
				  $bsmlBuilder,
				  $organism_name);

    &store_genome_attributes($organism_attributes_lookup,
			     $genome_elem,
			     $bsmlBuilder,
			     $organism_name);
    

    my $organism_elem = $bsmlBuilder->createAndAddOrganism( 
						    'genome'  => $genome_elem,
						    'genus'   => $genus,  
						    'species' => $species,
						    );
    if (defined($organism_elem)){		
	my $attrElem = $bsmlBuilder->createAndAddBsmlAttribute( $organism_elem,
							'translation_table',
							$translation_table );
	
	if (!defined($attrElem)){
	    $logger->logdie("Could not create <Attribute> for name 'translation_table' ".
			    "content '$translation_table' (for organism_name '$organism_name' ".
			    "GFF3 file '$gff_file')");
	}
    }
    else {
	$logger->logdie("Could not create <Organism> element object reference ".
			"for GFF3 file '$gff_file'");
    }

    #
    # editor:    sundaram@tigr.org
    # date:      2005-08-17
    # bgzcase:   2051
    # URL:       http://serval.tigr.org:8080/bugzilla/show_bug.cgi?id=2051
    # comment:   The <Sequence> will now be explicitly linked with the <Genome>
    #
    if (( exists $genome_elem->{'attr'}->{'id'} ) && ( defined ( $genome_elem->{'attr'}->{'id'} ) )  ){
	return ($genome_elem->{'attr'}->{'id'});
    }
    else {
	$logger->logdie("Genome id was not defined  for GFF3 file '$gff_file'");
    }
}


#-------------------------------------------------------------------
# is_file_readable()
#
#-------------------------------------------------------------------
sub is_file_readable {

    my ( $file) = @_;

    $logger->logdie("file was not defined") if (!defined($file));

    if (!-e $file){
	$logger->logdie("$file does not exist");
    }
    if (!-r $file){
	$logger->logdie("$file does not have read permissions");
    }
    if (-z $file){
	$logger->logdie("$file has no content");
    }

    return 1;
    

}#end sub is_file_readable()

#--------------------------------------------------------
# verify_and_set_outdir()
#
#
#--------------------------------------------------------
sub verify_and_set_outdir {

    my ( $outdir) = @_;

    #
    # strip trailing forward slashes
    #
    $outdir =~ s/\/+$//;
    
    #
    # set to current directory if not defined
    #
    if (!defined($outdir)){
	if (!defined($ENV{'OUTPUT_DIR'})){
	    $outdir = "." 
	}
	else{
	    $outdir = $ENV{'OUTPUT_DIR'};
	}
    }

    $outdir .= '/';

    #
    # verify whether outdir is in fact a directory
    #
    $logger->logdie("$outdir is not a directory") if (!-d $outdir);

    #
    # verify whether outdir has write permissions
    #
    $logger->logdie("$outdir does not have write permissions") if ((-e $outdir) and (!-w $outdir));

    #
    # store the outdir in the environment variable
    #
    $ENV{OUTPUT_DIR} = $outdir;
    
    return $outdir;

}#end sub verify_and_set_outdir()


#-------------------------------------------------------------------------------------------------
# Write bsml document to outdir
#
#-------------------------------------------------------------------------------------------------
sub write_out_bsml_doc {

    my ($outir, $bsmlBuilder, $gff_file, $no_report_stats) = @_;

    my $bsmldocument = $bsmlBuilder->{'doc_name'};

    #
    # If bsml document exists, copy it to .bak
    #
    if (-e $bsmldocument){

	my $bsmlbak = $bsmldocument . '.bak';
	rename ($bsmldocument, $bsmlbak);
	
	chmod (0666, $bsmlbak);
	
	$logger->info("Saving '$bsmldocument' as '$bsmlbak'");
    }

    print "Writing BSML document '$bsmldocument'\n" if $logger->is_debug;

    $logger->info("Writing BSML document '$bsmldocument'");
    
    $bsmlBuilder->write("$bsmldocument");
    
    if(! -e "$bsmldocument"){
	$logger->error("File not created '$bsmldocument'");
    }

    print "Changing permissions on  BSML document '$bsmldocument'\n" if $logger->is_debug;

    chmod (0777, "$bsmldocument");

    if ((defined($no_report_stats)) && ($no_report_stats == 1 )){
	#
	$logger->info("User requested that reporting be suppressed");
    }
    else {
	my $reportfile = $bsmldocument . ".stats";
	&produce_report($reportfile, $gff_file, $bsmlBuilder);
    }

    print "Created BSML file '$bsmldocument'\n";
}



#-------------------------------------------------------------------
# convert_strand()
#
# input:
#
# output:
#
# return:
#
# comment: gff3 records contain a strand column
#
# Lincoln Stein's notes excerpted from:
# http://song.sourceforge.net/gff3.shtml
#
#
#  __BEGIN__
#
# Column 7: "strand"
#
# The strand of the feature.  + for positive strand (relative to the
# landmark), - for minus strand, and . for features that are not
# stranded.  In addition, ? can be used for features whose strandedness
# is relevant, but unknown.
#
# __END__
#
#
# Note that acceptable values for the BSML attribute 'strand' according to the DTD are:
#
# 	strand (std-not-set | ss | ds | mixed | std-other) #IMPLIED
#
#
#-------------------------------------------------------------------
sub convert_strand {

    my ($strand) = @_;


    if ($strand eq '-'){
	$strand = 0;
    }
    else {
	$strand = 1;
    }
    
    return $strand;
}

#---------------------------------------------------------
# store_locus_as_cross_reference()
#
# input:
#
# output:
#
# return:
#
# comment:
#
#
#---------------------------------------------------------
sub store_locus_as_cross_reference {

    my ($bsmlElement, $attributes, $bsmlBuilder, $id, $newId) = @_;

    my $identifier_type = 'locus';

    if (( exists $attributes->{$identifier_type}->[0]) &&
	( defined($attributes->{$identifier_type}->[0]))) {
	
	my $identifier = $attributes->{$identifier_type}->[0];

	#
	# Create <Cross-reference> element object for the feat_name
	#
	my $xref_elem = $bsmlBuilder->createAndAddCrossReference( 'parent'          => $bsmlElement,
							  'id'              => $bsmlBuilder->{'xrefctr'}++,
							  'database'        => $prefix,
							  'identifier'      => $identifier,
							  'identifier-type' => $identifier_type
							 );
	if (!defined($xref_elem)){
	    $logger->logdie("Could not create <Cross-reference> for id '$id' ".
			    "newId '$newId' database '$prefix' identifier '$identifier' ".
			    "identifier-type '$identifier_type' in GFF3 file '$gff_file'");
	}							
    }
}


#------------------------------------------------------
# store_gff_attributes_as_attribute_list()
#
# input:
#
# output:
#
# return:
#
# comment:
#------------------------------------------------------
sub store_gff_attributes_as_attribute_list {
    
    my ($element, $attributes, $type) = @_;

    if (( exists $attributes->{$type}) &&
	(defined($attributes->{$type})) && 
	(scalar(@{$attributes->{$type}}) > 0 )) {

	my $name;

	if ($type eq 'ec_number'){
	    $name = 'EC';
	}
	elsif ($type eq 'Ontology_term'){
	    $name = 'GO';
	}
	else {
	    # Currently, no other types require translation
	}
	
	foreach my $ec (sort @{$attributes->{$type}} ) {
	    
	    my $list = [];
	    
	    push( @{$list},  { name => $name,
			       content => $ec });
	    
	    $element->addBsmlAttributeList($list);
	}	
    }
}


#--------------------------------------------------------
# verify_and_set_fastadir()
#
#--------------------------------------------------------
sub verify_and_set_fastadir {

    my ( $fastadir) = @_;

    #
    # strip trailing forward slashes
    #
    $fastadir =~ s/\/+$//;
    
    #
    # set to current directory if not defined
    #
    if (!defined($fastadir)){
	$fastadir = "." 
    }

    $outdir .= '/';

    #
    # verify whether outdir is in fact a directory
    #
    if (!-e $fastadir){
	$logger->warn("fastadir '$fastadir' does not exist.  Attempting to create now.");
	mkdir($fastadir,0777);
    }
    else {
	if (!-d $fastadir){
	    $logger->logdie("$fastadir is not a directory") 
	}
	else {
	    if (!-w $fastadir){
		$logger->logdie("fastadir '$fastadir' does not have write permissions");
	    }
	}
    }

    return $fastadir;

}#end sub verify_and_set_fastadir()


#------------------------------------------------------
# print_usage()
#
#------------------------------------------------------
sub print_usage {

    print STDERR "SAMPLE USAGE:  $0 [-d debug_level] [-f fastadir] -g gff_file [-h] [--id_generator_mapping_file] [--id_repository] [-l log4perl] [-m] [--no_id_generator] [--no_report_stats] [-o outdir]\n".
    "  -d|--debug_level            = Optional - Ergatis::Logger log4perl logging level.  Default is 0\n".
    "  -f|--fastadir               = Optional - Directory to which the multi-fasta file containing all fasta sequences will be written (Default is current working directory)\n".
    "  -g|--gff_file               = GFF3 file to be converted into BSML encoding\n".
    "  -h|--help                   = Optional - Display pod2usage help screen\n".
    "  --id_generator_mapping_file = Optional - old ID to new ID mapping file.   See man for details\n".
    "  --id_repository             = Optional - IdGenerator compliant directory.   See man for details\n".
    "  -h|--help                   = Optional - Display pod2usage help screen\n".
    "  -l|--log4perl               = Optional - Log4perl log file (default: /tmp/brcGff3ToBsml.pl.log)\n".
    "  -m|--man                    = Optional - Display pod2usage pages for this utility\n".
    "  --no_id_generator           = Optional - disables IdGenerator support\n".
    "  --no_report_stats           = Optional - disables reporting of Feature, Sequence types counted as well as BSML element objects created\n".
    "  -o|--outdir                 = Optional - Directory to which the .bsml document will be written (Default is current working directory)\n";
    exit 1;

}


#----------------------------------------------------------------------
# create_out_file()
#
#----------------------------------------------------------------------
sub create_out_file {

    my ($dir, $basename, $ext) = @_;

    my $file;

    $file = $dir . "/" . $basename . '.' . $ext;

    return $file;
}

#----------------------------------------------------------------------
# fetch_gff_record()
#
#----------------------------------------------------------------------
sub fetch_gff_record {

    my ($gffrecords, $line, $linectr, $parent_lookup,
	$organism_attributes_lookup, $cdsSegmentLookup) = @_;

    #
    # Per the GFF3 spec, the core 9 columns are tab-separated and NOT
    # space-separated.
    #
    my @columns = split(/\t+/, $line);

    my ($seq_id, $source, $type, $start, $end, $score,
	$strand, $phase, $attributes) = &verify_gff_columns(\@columns, $line, $linectr);
    
    $report_stats->{'gff'}->{$type}++;

    my ($id, $attribute_hash) = &fetch_gff_attributes($attributes,
						      $parent_lookup, 
						      $known_attributes,
						      $linectr,
						      $start,
						      $end,
						      $organism_attributes_lookup,
						      $cdsSegmentLookup,
						      $type,
						      $gffrecords);
 
    if($id){
	
	## Keep track of all contig and supercontig sequence types
	if (($type eq 'contig') ||
	    ($type eq 'supercontig')){
	    push(@{$primarySequences->{'list'}}, $id);
	    $primarySequences->{'lookup'}->{$id}++;
	}

	&store_gff_record_reference($gffrecords, $id, $start, $end, 
				    $linectr, $seq_id, $source, $type, 
				    $score, $strand, $phase, $attribute_hash);

    }
    else {
	$logger->warn("ID was not defined for start '$start' end '$end' ".
		      "seq_id '$seq_id' source '$source' type '$type' ".
		      "score '$score' strand '$strand' phase '$phase' ".
		      "in GFF file '$gff_file'");
    }
    
}

#----------------------------------------------------------------------
# verify_gff_columns()
#
#----------------------------------------------------------------------
sub verify_gff_columns {

    my ($columns, $line, $linectr) = @_;

    for (my $i=0; $i < 8; $i++){

	if (!defined($columns->[$i])){
	    $logger->logdie("column '$i' was not defined at line '$linectr'" . Dumper $columns);
	}
    }

    return (@{$columns});
}

#-------------------------------------------------------------------
# fetch_gff_attributes()
#
#-------------------------------------------------------------------
sub fetch_gff_attributes {

    my ($attributes, $parent_lookup, $known_attributes, $linectr,
	$start, $end, $organism_attributes_lookup, $cdsSegmentLookup, $type, $gffrecords) = @_;

    my $attribute_hash = {};

    #
    # If this GFF3 annotation record contains organism_name, then any Dbxref attribute
    # should be associated with the Genome element
    #
    my $organism_flag=0;

    if ($attributes =~ /organism_name/){
	$organism_flag=1;
    }

    my @atts = split(/;/, $attributes);
    
    my $key2attr = {};

    my $id;

    my $Dbxref_values;

    my $cdsid;

    #
    # Need to process the ID attribute first
    #
    foreach my $group ( @atts ){
	
	# remove leading whitespace
	$group =~ s/^\s*//;
	
	# remove trailing whitespace
	$group =~ s/\s*$//;
	
	if ($group =~ /^\s*$/){
	    ## Encountered trailing semicolon and whitespace
	    next;
	}

	my ($key, $values) = split(/=/, $group);

	if ($key eq 'ID'){
	    
	    #
	    # Capture the ID attribute
	    #
	    $id = $values;
	    
	    if ($id =~ /,/){
		$logger->logdie("$0 found ID attribute with possibly multiple values for GFF3 file '$gff_file'");
	    }

	    $id = &clean_id($id);
		
	    if (uc ($type) eq  'CDS'){

		## e.g. id = CDS.17922
		$cdsid = $id;
		## e.g. cdsid = CDS.17922

		$id = &ensure_uniqueness($id, $start, $end);
		## e.g. id = CDS.17922.v2 (on 2nd encounter)
		## e.g. id = CDS.17922 (on 1st encounter)

		$cdsSegmentLookup->{'cds'}->{$id} = $cdsid;
		push ( @{$cdsSegmentLookup->{'cdsorig'}->{$cdsid}}, $id);
	    }
	    else {
		$id = &ensure_uniqueness($id, $start, $end);
	    }


	}
	else {
	    #
	    # Review and qualify the types of attributes being stored
	    # in the GFF3 attribute section
	    #
	    if ( ! exists $known_attributes->{$key} ) {
		$logger->warn("$0 encountered new type of attribute '$key' for ID '$id' at line '$linectr' of GFF file '$gff_file'");
	    }
	    else {

		##
		## Hacks required for non-compliant GFF encoding with respect to Dbxref attributes
		##
		if (($key eq 'Dbxref') && ($organism_flag == 1)){
		    $Dbxref_values .= "$values,";
		}
		elsif (($key eq 'Dbxref') && ($organism_flag == 0)){

		    my @dbxrefValues = split(/,/,$values);

		    foreach my $dbxrefValue ( @dbxrefValues ){

			$dbxrefValue = URI::Escape::uri_unescape(URI::Escape::uri_unescape($dbxrefValue));

			push( @{ $attribute_hash->{'Dbxref'}},  $dbxrefValue);
		    }
		}
		else {
		    $key2attr->{$key} = $values;
		}
	    }
	}
    }

    ##
    ## Hacks required for non-compliant GFF encoding with respect to Dbxref attributes
    ##
    if (($organism_flag == 1) && (defined($Dbxref_values)) && (length($Dbxref_values) > 0 )) {
	
	chop $Dbxref_values;

	push( @{ $organism_attributes_lookup->{'Dbxref'}},  $Dbxref_values);
	
    }

    foreach my $key (keys %{$key2attr} ) {

	my $values = $key2attr->{$key};

	if ($organism_flag == 1){
	    
	    if (($key eq 'organism_name') || ($key eq 'Dbxref') || ($key eq 'strain') || ($key eq 'translation_table')){
		#
		# All of these attributes should be stored as part of the Genome/Organism
		#
		if ($key eq 'Dbxref'){

		    push( @{ $organism_attributes_lookup->{$key} },  $values);
		}
		else {
		    $organism_attributes_lookup->{$key} = $values;
		}
	    }
	    
	    next;
	}
	else {
	    my @vals = split(/,/, $values);
	    
	    foreach my $val (@vals){

		$val = URI::Escape::uri_unescape(URI::Escape::uri_unescape($val));
		
		# remove leading double quotation marks
		$val =~ s/^\"//;
		
		# remove trailing double quotation marks
		$val =~ s/\"$//;
		
		if ($key eq 'Parent'){
		
		    #
		    # This will be important for creating feature_relationships between this subject feature and the object (parent) feature
		    #
		    $val = &clean_id($val);
			

		    if (!defined($id)){
			## create some temporary identifier based on
			## the parent identifier, start, end and feature type
			$id = $val . '_'. $start . '_' . $end . '_' . $type;
			
		    }

		    push (@{$parent_lookup->{$id}}, $val);

		    if (defined($cdsid)){
			## dealing with some CDS
			push ( @{$cdsSegmentLookup->{'parent'}->{$val}}, $id);
		    }
		}
		else {
		    #
		    # All other attributes are stored
		    #
		    push ( @{$attribute_hash->{$key}}, $val)
		}
	    }
	}
    }

    if (!defined($id)){
	#
	# Per the GFF3 spec, the ID attribute is mandatory for features with
	# children subfeatures.
	#
	$logger->warn("ID attribute was not defined for GFF record number '$linectr' for GFF3 file '$gff_file'");
    }


    return ($id, $attribute_hash);
}


#-------------------------------------------------------------------
# store_gff_record_reference()
#
#----------------------------------------------------------------------
sub store_gff_record_reference {

    my ($gffrecords, $id, $start, $end, $linectr, $seq_id, $source, $type, $score, $strand, $phase, $attribute_hash) = @_;


    if ((exists $gffrecords->{$id}) && (defined($gffrecords->{$id} ))) {

	if ($type eq 'exon'){
	    $logger->info("Encountered duplicate ID '$id' start '$start' end '$end' at line '$linectr' for type '$type' in GFF3 file '$gff_file'");
	    #
	    # Store the Parent attribute value(s)
	    #
	    push ( @{$gffrecords->{$id}->{'Parent'}}, $attribute_hash->{'Parent'});
	}
	else {
	    $logger->logdie("Encountered duplicate ID '$id' start '$start' end '$end' at line '$linectr' for type '$type' in GFF3 file '$gff_file'");
	}
	
    }
    else {
	
	$gffrecords->{$id} = { 'seq_id' => $seq_id,
			      'source' => $source,
			      'type'   => $type,
			      'score'  => $score,
			      'strand' => $strand,
			      'phase'  => $phase,
			      'start'  => $start,
			      'end'    => $end };
	
	foreach my $attributetype (sort keys %{$attribute_hash} ) {
	    
	    #
	    # Load the attributes onto the gffrecord hash
	    # 
	    $gffrecords->{$id}->{'attributes'}->{$attributetype} = $attribute_hash->{$attributetype};
	}
    }
}


#----------------------------------------------------------------------
# fetch_gff_comment()
#
#----------------------------------------------------------------------
sub fetch_gff_comment {

    my ($gffcomments, $line, $linectr) = @_;

    $gffcomments->{$line}++;
    
}



#----------------------------------------------------------------------
# storePrimarySequences()
#
#----------------------------------------------------------------------
sub storePrimarySequences {

    my ($primarySequences, $gffrecords, $bsmlBuilder, $genome_id_lookup, 
	$seqid_to_genome_lookup, $sequenceIdToBsmlSequenceLookup,
	$sequenceIdToBsmlFeatureTableLookup, $organism_attributes_lookup) = @_;

    print "Storing primary sequences\n";

    my $translation_table;

    my $genomeId;

    foreach my $id ( @{$primarySequences->{'list'}} ){
	
	my $attributes = $gffrecords->{$id}->{'attributes'};
	
	if (!$organismDataProcessed){
	    ## Only store the organism information once
	    $organismDataProcessed=1;

	    ($genomeId, $translation_table) = &store_organism_data($attributes,
								   $genome_id_lookup,
								   $bsmlBuilder,
								   $organism_attributes_lookup, 
								   $gffrecords->{$id}->{'source'});
	    
	    if (defined($genomeId)){
		if (exists $gffrecords->{$id}->{'seq_id'}){
		    $seqid_to_genome_lookup->{$gffrecords->{$id}->{'seq_id'}} = $genomeId;
		}
		else {
		    $logger->logdie("seq_id was not defined for primary sequence with ".
				    "id '$id' in GFF3 file '$gff_file'");
		}
	    }
	    else {
		$logger->logdie("genomeId was not defined for primary sequence with ".
				"id '$id' in GFF3 file '$gff_file'");
	    }
	}


	&createBsmlForPrimarySequence($id,
				      $gffrecords,
				      $sequenceIdToBsmlSequenceLookup,
				      $sequenceIdToBsmlFeatureTableLookup,
				      $genomeId,
				      $bsmlBuilder);	
    }
    
    return $translation_table;

}

#----------------------------------------------------------------------
# store_organism_data()
#
#----------------------------------------------------------------------
sub store_organism_data {

    my ($attributes, $genome_id_lookup, $bsmlBuilder, $organism_attributes_lookup, $source) = @_;
    
      if (  exists $organism_attributes_lookup->{'organism_name'} && ( defined( $organism_attributes_lookup->{'organism_name'} ))) {

	  my $organism_name = $organism_attributes_lookup->{'organism_name'};

	  ## undefine the organism_name hashkey-value here
	  delete $organism_attributes_lookup->{'organism_name'};

	  my $genome_id;
	  my $translation_table;

	  if (( exists $organism_attributes_lookup->{'translation_table'}) &&
	      ( defined($organism_attributes_lookup->{'translation_table'})) ){ 

	      $translation_table = $organism_attributes_lookup->{'translation_table'};

	      ## undefine the translation_table hashkey-value here
	      delete $organism_attributes_lookup->{'translation_table'};


	  }
	  else {
	      #
	      # translation_table attribute is mandatory
	      #
	      $logger->logdie("translation_table attribute was not defined for organism_name '$organism_name' in GFF3 file '$gff_file'");
	  }


	  if (( exists $genome_id_lookup->{$organism_name}) && 
	      ( defined($genome_id_lookup->{$organism_name}))) {
	      
	      $genome_id = $genome_id_lookup->{$organism_name};
	      
	      # This means that only the attributes associated with the first
	      # organism record will be stored in the BSML //Genome/ component.
	  }
	  else {
	      
	      $genome_id = &create_genome_component( $bsmlBuilder,
						     $organism_name,
						     $translation_table,
						     $attributes,
						     $organism_attributes_lookup,
						     $source);
	      
	      $genome_id_lookup->{$organism_name} = $genome_id;
	  }
	  
	  return ($genome_id, $translation_table);
	  
      }
    else {
	
	$logger->logdie("organism_name was not defined in organism_attributes_lookup");
    }
}

#----------------------------------------------------------------------
# get_complement()
#
#----------------------------------------------------------------------
sub get_complement {

    my ($strand, $id) = @_;

    my $complement;

    if ($strand eq '+'){
	$complement = 0;
    }
    elsif ($strand eq '-'){
	$complement = 1;
    }
    elsif ($strand eq '.'){
	$complement = undef;
    }
    elsif ($strand eq '?'){
	$complement = undef;
    }
    else {
	$logger->logdie("Unrecognized value for strand '$strand' ID '$id' for GFF3 file '$gff_file'");
    }
    
    return $complement;
}

#----------------------------------------------------------------------
# get_logger()
#
#----------------------------------------------------------------------
sub get_logger {

    my ($log4perl, $debug_level) = @_;

    #
    # initialize the logger
    #
    $log4perl = "/tmp/brcGff3ToBsml.pl.log" if (!defined($log4perl));
    my $mylogger = new Ergatis::Logger('LOG_FILE'=>$log4perl,
				     'LOG_LEVEL'=>$debug_level);
    
    my $logger = Ergatis::Logger::get_logger(__PACKAGE__);
    
    return $logger;
}


#----------------------------------------------------------------------
# get_filehandle()
#
#----------------------------------------------------------------------
sub get_filehandle {

    my ($gff_file) = @_;

    my $fh;
    
    # This is not a file handle- is an actual file.
    
    if ($gff_file =~ /\.gff3l\.(gz|gzip)$/) {
	# This file has been zipped
	open ($fh, "<:gzip", $gff_file) || $logger->logdie("Could not open zipped file '$gff_file': $!");
	
    }
    else {
	# Regular unzipped file
	open ($fh, "<$gff_file") || $logger->logdie("Could not open file '$gff_file': $!");
    }

    return $fh;
}

#-------------------------------------------------------------------
# clean_id()
#
#-------------------------------------------------------------------
sub clean_id {
    
    my ($id) = @_;
    
    # remove leading double quotation marks
    $id =~ s/^\"//;
    
    # remove trailing double quotation marks
    $id =~ s/\"$//;
    
    # all pipe symbols are transformed into underscores
    $id =~ s/\|/_/;
    
    return $id;
}

#----------------------------------------------------------------------
# ensure_uniqueness()
#
#----------------------------------------------------------------------
sub ensure_uniqueness {

    my ($id, $start, $end) = @_;

    #
    # This is to ensure that identifiers propagated into BSML and chado are unique
    #
    if (! exists $id_lookup->{$id}) {
	$id_lookup->{$id}++;
    }
    else {
	
	## Found duplicate ID for two different records
	##

	$id_lookup->{$id}++;
	
	my $newid = $id . '.v' . $id_lookup->{$id};
	    
	$id = $newid;

	## Will propagate isoform identifier
	## 
    }
    
    return $id;
}


#-------------------------------------------------------------------
# store_genome_cross_reference()
#
#-------------------------------------------------------------------
sub store_genome_cross_reference {

    my ($organism_attributes_lookup, $genome_elem, $bsmlBuilder, $organism_name) = @_;

    if (( exists $organism_attributes_lookup->{'Dbxref'} ) &&
	  (defined ($organism_attributes_lookup->{'Dbxref'} ))){

	my $processedOrganismAttributes = {};

	foreach my $xref ( @{$organism_attributes_lookup->{'Dbxref'}} ) {

	    my @xrefs = split(/,/,$xref);
	    
	    foreach my $xref (@xrefs) {

		$xref = URI::Escape::uri_unescape(URI::Escape::uri_unescape($xref));
	
		if (!defined($xref)){
		    $logger->logdie("xref was not defined for GFF3 file '$gff_file'");
		}
		
		if (! exists $processedOrganismAttributes->{$xref} ) {
		    ## Only process this Dbxref one time

		    if ($xref =~ /:/){
			
			my ($database, $identifier) = split(/:/, $xref);
	
			if (!defined($database)){
			    $logger->logdie("database was not defined for xref '$xref' while processing ".
					    "Sequence with organism name '$organism_name' in GFF3 file '$gff_file'");
			}
			
			## Create <Cross-reference> element object for the feat_name
			my $xref_elem = $bsmlBuilder->createAndAddCrossReference(
										 'parent'          => $genome_elem,
										 'id'              => $bsmlBuilder->{'xrefctr'}++,
										 'database'        => $database,
										 'identifier'      => $identifier
										 );
			if (!defined($xref_elem)){
			    $logger->logdie("Could not create <Cross-reference> for Sequence with organism ".
					    "name '$organism_name' ".
					    "for database '$database' identifier '$identifier' in GFF3 file ".
					    "'$gff_file'");
			}							
		    }
		    else {
			$logger->warn("Could not derive database, identifier while processing Sequence with organism name ".
				      "'$organism_name' in GFF3 file '$gff_file'");
		    }
		    ## Remember that we processed this Dbxref
		    $processedOrganismAttributes->{$xref}++;
		}
	    }
	}

	## remove all Dbxref haskey-values here
	delete $organism_attributes_lookup->{'Dbxref'};
    }
}

#-------------------------------------------------------------------
# store_genome_attributes()
#
#-------------------------------------------------------------------
sub store_genome_attributes {

    my ($organism_attributes_lookup, $genome_elem, $bsmlBuilder, $organism_name) = @_;

    foreach my $attributetype (sort keys %{$organism_attributes_lookup}){

	my $attributestring = $organism_attributes_lookup->{$attributetype};
	
	my @attributes = split(/,/,$attributestring);
	
	foreach my $attribute (@attributes){
	    
	    $attribute = URI::Escape::uri_unescape(URI::Escape::uri_unescape($attribute));

	    my $attribute_elem = $bsmlBuilder->createAndAddBsmlAttribute( $genome_elem,
								  $attributetype,
								  $attribute );
	    
	    if (!defined($attribute_elem)){
		$logger->logdie("Could not create <Attribute> for name '$attributetype' content '$attribute' (for organism_name '$organism_name'  GFF3 file '$gff_file')");
	    }
	}
	
	## undef the attribute hashkey-values here
	delete $organism_attributes_lookup->{$attributetype};
    }
}

#------------------------------------------------------------------
# mergeCDSSegments()
#
#------------------------------------------------------------------
sub mergeCDSSegments {

    my ($gffrecords, $cdsSegmentLookup, $parent_lookup) = @_;

    print "Merging CDS segments\n";

    ## 	 
    ## Merging CDS segments based on common CDS identfiers. 	 
    ##
    foreach my $master_cds_id ( keys %{$cdsSegmentLookup->{'cdsorig'}} ) { 	 
	## e.g. master_cds_id = CDS.17922
	
	my $member_cds_array = $cdsSegmentLookup->{'cdsorig'}->{$master_cds_id}; 	 
	
#	$logger->fatal("member_cds_array for master_cds_id '$master_cds_id':". Dumper $member_cds_array);
	
	# Get the smallest start and the largest end 	 
	# 	 
	my $minfmincds=0; 	 
	my $maxfmaxcds=0; 	 
	my $strand; 	 
	my $cdsctr = 0; 	 
	
	my $masterCdsGffAttributes = $gffrecords->{$master_cds_id};
	my $lastcdsrecordsparents;

	foreach my $cds_segment_id ( @{$member_cds_array} ) { 	 

#	    $logger->fatal("Processing cds_segment_id '$cds_segment_id' for master_cds_id '$master_cds_id'");
	    
	    my $cdsrecord = $gffrecords->{$cds_segment_id}; 	 
	    
	    my $start  = $cdsrecord->{'start'}; 	 
	    my $end    = $cdsrecord->{'end'}; 	 
	    my $strand = $cdsrecord->{'strand'}; 	 
	    
	    if ($cdsctr == 0 ) { 	 
		$minfmincds = $start; 	 
		$maxfmaxcds = $end; 	 
	    } 	 
	    else { 	 
		if ($start < $minfmincds ){ 	 
		    $minfmincds = $start; 	 
		} 	 
		if ($end > $maxfmaxcds){ 	 
		    $maxfmaxcds = $end; 	 
		} 	 
	    } 	 
	    $cdsctr++; 	 
	     
	     $lastcdsrecordsparents = $parent_lookup->{$cds_segment_id};
	    
	    ## eliminate the individual segments 	
	    if ($master_cds_id ne $cds_segment_id){
		delete $gffrecords->{$cds_segment_id}; 	
		delete $parent_lookup->{$cds_segment_id};
#		$logger->fatal("For master_cds_id '$master_cds_id' deleting cds_segment_id ' $cds_segment_id'");
	    } 	
# 	    else {
# 		$logger->fatal("master_cds_id '$master_cds_id' == cds_segment_id ' $cds_segment_id'");
# 	    }
	}
	

	$gffrecords->{$master_cds_id} = $masterCdsGffAttributes;
	$gffrecords->{$master_cds_id}->{'start'} = $minfmincds; 	 
	$gffrecords->{$master_cds_id}->{'end'} = $maxfmaxcds; 	 

    } 	 

}


#------------------------------------------------------------------
# createTranscriptFeatureStubs()
#
#------------------------------------------------------------------
sub createTranscriptFeatureStubs {

    my ($gffrecords, $parent_lookup) = @_;

    print "Creating transcript feature stubs\n";

#     my $child_lookup = {};

#     my ($parentarray, $child);

#     #
#     # Create a child-to-parent lookup
#     #
#     while (($child, $parentarray) = each %{$parent_lookup}){

# 	foreach my $parent ( @{$parentarray} ) {
	    
# 	    push ( @{$child_lookup->{$parent}}, $child);
# 	}
#     }

    foreach my $featureId ( keys %{$gffrecords}){

	if ( exists $gffrecords->{$featureId}->{'type'}){

	    if (lc($gffrecords->{$featureId}->{'type'}) eq 'cds'){

		## Found a CDS feature
		my $cdsHasParentTranscriptFeature=0;

		if (exists $parent_lookup->{$featureId}){

		    foreach my $parentId ( @{$parent_lookup->{$featureId}} ){

			if ( exists $gffrecords->{$parentId} ){
			    if ( exists $gffrecords->{$parentId}->{'type'} ){
				if ((  lc($gffrecords->{$parentId}->{'type'}) eq 'transcript') || 
				    (  lc($gffrecords->{$parentId}->{'type'}) eq 'mrna')) {
				    ## Found a parent that is a transcript or mRNA
				    ## feature type.
				    $cdsHasParentTranscriptFeature = 1;
				    
				    ## Retrieve all of the CDS feature's attribute
				    my $cdsAttributes = $gffrecords->{$featureId}->{'attributes'};

				    ## Assign all of the attributes to the parent transcript
				    ## or mRNA feature
				    foreach my $attribute_type ( keys %{$cdsAttributes} ) {
					
					foreach my $val ( @{$cdsAttributes->{$attribute_type}} ) {
					
					    push ( @{$gffrecords->{$parentId}->{$attribute_type}}, $val ) ;
					
					}
				    }
				}
				else {
				    $logger->info("CDS with id '$featureId' had a parent ".
						  "feature with type '$gffrecords->{$parentId}->{'type'}' ".
						  "in GFF file '$gff_file'");
				}
			    }
			}
		    }
		}
		 
		if ($cdsHasParentTranscriptFeature == 0) {
		    $logger->info("CDS with id '$featureId' did not have a ".
				  "parent transcript/mRNA feature in ".
				  "GFF file '$gff_file'");
		    
		    ## no mRNA/transcript feature was found for this CDS
		    ## Create some transcript Feature stub and associate all of the CDS feature's
		    ## attributes with the stub.
		    ##
		    my $transcriptStubId = $featureId . '_transcript_stub';
		    
		    $transcriptStubId = &get_id_from_id_generator($transcriptStubId, 'transcript');

		    push ( @{$parent_lookup->{$featureId}}, $transcriptStubId);
		    
		    my $cds_record = $gffrecords->{$featureId};

		    foreach my $attribute_type ( keys %{$cds_record} ) {
			
			if ($attribute_type ne 'attributes'){

			    $gffrecords->{$transcriptStubId}->{$attribute_type} = $cds_record->{$attribute_type};

			}
			else {

			    foreach my $aux_attrs (keys %{$cds_record->{'attributes'}} ){ 
				
				foreach my $val ( @{$cds_record->{'attributes'}->{$aux_attrs}} ) {
				    
				    push ( @{$gffrecords->{$transcriptStubId}->{'attributes'}->{$aux_attrs}}, $val ) ;

				}
			    }
			    
			    delete $cds_record->{'attributes'};
			}

			$gffrecords->{$transcriptStubId}->{'type'} = 'transcript';

		    }
		    
		    $report_stats->{'new'}->{'transcript'}++;

		}
	    }
	}
    }
}


#------------------------------------------------------------
# derivePolypeptideFeatureSequenceForCds()
#
#------------------------------------------------------------
sub derivePolypeptideFeatureSequenceForCds {

    ## Here we'll derive the nucleotide sequence for the CDS features based on the
    ## CDS features' coordinates (start,end).
    ## We'll also create a polypeptide Features and Sequences.  The original CDS 
    ## FASTA sequence (amino acid sequence) will be associated with the polypeptide.

    print "Deriving polypeptide Features for each CDS\n";

    my ($gffrecords, $parent_lookup) = @_;

    foreach my $cdsId ( keys %{$gffrecords} ) {
	
	if ( exists $gffrecords->{$cdsId}->{'type'}) {

	    my $type = $gffrecords->{$cdsId}->{'type'};
	    
	    if ( uc($type) eq 'CDS'){
		## dealing with some CDS feature	

		if ( exists $gffrecords->{$cdsId}->{'fasta'}){

		    ## create some polypeptide feature (at this point I assume that none of the BRCs including TIGR
		    ## are storing polypeptide/protein features in their GFF3 files.
		    
		    my $polypeptide_id = $cdsId . '_polypeptide';
		    
		    $polypeptide_id = &getIdFromIdGenerator($polypeptide_id, 'polypeptide');
		    
		    my $cdsGffRecord = $gffrecords->{$cdsId};
		    
		    ## assign all CDS' attributes to the new polypeptide feature
		    foreach my $attribute_type ( keys %{$cdsGffRecord} ){
			$gffrecords->{$polypeptide_id}->{$attribute_type} = $cdsGffRecord->{$attribute_type};
		    }
		    
		    ## Assign the feature type
		    $gffrecords->{$polypeptide_id}->{'type'} = 'polypeptide';
		    
		    ## assign the CDS' amino acid sequence to the polypeptide feature
		    my $tempfasta = $gffrecords->{$cdsId}->{'fasta'};

		    $gffrecords->{$polypeptide_id}->{'fasta'} = $tempfasta;
		
		    $report_stats->{'new'}->{'polypeptide'}++;
		
		    ## the CDS becomes the parent of the new polypeptide feature
		    push( @{$parent_lookup->{$polypeptide_id}}, $cdsId);

		}
		else {
		    $logger->error("No FASTA sequence for CDS '$cdsId' in GFF file '$gff_file' ".
				   "with GFF record for CDS '$cdsId':" . Dumper $gffrecords->{$cdsId});

		    ## There was no FASTA sequence associated with CDS '$cdsId' (old ID '$oldid').  
		    ## This CDS feature's nucleotide sequence has been derived from the associated contig's FASTA sequence.  
		    ## Note that the derived polypeptide feature will not have any amino acid sequence associated with it.");
		}
	    }
	}
    }    
}

#--------------------------------------------------------
# produce_report()
#
#--------------------------------------------------------
sub produce_report  {

    my ($reportfile, $gff_file, $bsmlBuilder) = @_;

    open (REPORTFILE, ">$reportfile") || $logger->logdie("Could not open file '$reportfile' for output: $!");


    print REPORTFILE "$0 processed the following GFF3 file: '$gff_file'\n";


    print REPORTFILE "Counted the following number of GFF3 record types\n";

    foreach my $type (keys %{$report_stats->{'gff'}} ) {

	print REPORTFILE "$type $report_stats->{'gff'}->{$type}\n";
    }



    print REPORTFILE "Created the following number of record types\n";

    foreach my $type (keys %{$report_stats->{'new'}} ) {

	print REPORTFILE "$type $report_stats->{'new'}->{$type}\n";
    }


    print REPORTFILE "BsmlBuilder created the following BSML elements\n";

    my $bsmlElementCounts = $bsmlBuilder->getBsmlElementCounts();

    foreach my $type (keys %{$bsmlElementCounts}){
	print REPORTFILE "$type $bsmlElementCounts->{$type}\n";
    }
}

#--------------------------------------------------------
# get_id_mapping_lookup()
#
#--------------------------------------------------------
sub get_id_mapping_lookup {

    my ($gff_file, $outdir, $id_generator_mapping_file, $id_mapping_lookup) = @_;

    if (!defined($id_generator_mapping_file)){
	
	## strip trailing forward slashes
	$outdir =~ s/\/+$//;

	$gff_file = basename($gff_file);
	
	$id_generator_mapping_file = $outdir . '/brcGff3ToBsml.pl.' . $gff_file . '.mapping_file.txt';
    }
    
    if (-e $id_generator_mapping_file){
	
	if (-r $id_generator_mapping_file){
	    
	    if (-w $id_generator_mapping_file){
		
		if (-f $id_generator_mapping_file){
		    
		    open (MAPPINGFILE, "<$id_generator_mapping_file") || $logger->logdie("Could not open file '$id_generator_mapping_file' for input: $!");
		    
		    while (my $line = <MAPPINGFILE>){
			chomp $line;
			
			my ($oldid, $newid) = split(/\s+/, $line);
			
			$id_mapping_lookup->{'old2new'}->{$oldid} = $newid;
			$id_mapping_lookup->{'new2old'}->{$newid} = $oldid;
		    }
		}
		else {
		    $logger->logdie("'$id_generator_mapping_file' was not a file");
		}
	    }
	    else {
		$logger->logdie("'$id_generator_mapping_file' did not have write permissions");
	    }
	}
	else {
	    $logger->logdie("'$id_generator_mapping_file' did not have read permissions");
	}
    }
    else {
	qx{touch $id_generator_mapping_file};
    }

    return $id_generator_mapping_file;

}

#--------------------------------------------------------------
# get_id_from_id_generator()
#
#--------------------------------------------------------------
sub get_id_from_id_generator {

    my ($id, $type) = @_;

    my $returnid = $id;

    if ((defined($no_id_generator)) && ($no_id_generator == 1)){
	# not using IdGenerator this round!
	$id_mapping_lookup->{'old2new'}->{$id} = $id;
	$id_mapping_lookup->{'new2old'}->{$id} = $id;
	
    }
    else {

	if ( exists $id_mapping_lookup->{'old2new'}->{$id}) {
	 
	    $returnid = $id_mapping_lookup->{'old2new'}->{$id};
	}
	else {
	    
	    $type = lc($type);

	    $returnid = $idcreator->next_id( project => $project,
					     type    => $type );
	    
	    $id_mapping_lookup->{'old2new'}->{$id} = $returnid;
	    $id_mapping_lookup->{'new2old'}->{$returnid} = $id;
	    
	}
    }
    
    return $returnid;

}

#--------------------------------------------------------------
# write_mapping_file()
#
#--------------------------------------------------------------
sub write_mapping_file {

    my ($id_generator_mapping_file, $id_mapping_lookup) = @_;

    if (!defined($id_generator_mapping_file)){
	$logger->logdie("id_generator_mapping_file was not defined");
    }

    if (-e $id_generator_mapping_file){
	my $bakfile = $id_generator_mapping_file . '.bak';
	rename($id_generator_mapping_file, $bakfile);
    }

    open (MAPPINGFILE, ">$id_generator_mapping_file") || $logger->logdie("Could not open file '$id_generator_mapping_file' for output: $!");

    foreach my $oldid (keys %{$id_mapping_lookup->{'old2new'}} ) {

	my $newid = $id_mapping_lookup->{'old2new'}->{$oldid};

	print MAPPINGFILE "$oldid\t$newid\n";

    }
}

##---------------------------------------------------------------
## parseGffFile()
##
##---------------------------------------------------------------
sub parseGffFile {

    my $fasta_header;
    my $fasta_sequence;
    my $fasta_directive_detected = 0;
    my $organism_attributes_lookup = {};

    my $gffComment = {};
    my $gffrecords = {};
    my $cdsSegmentLookup = {};
    my $parent_lookup = {};

    my $linectr=0;
    
    print "Parsing GFF file '$gff_file'\n";

    # Get gff file handle
    my $fh = &get_filehandle($gff_file);

    while (my $line = <$fh>){
	
	chomp $line;
	
	$linectr++;
	
	if (! ($fasta_directive_detected)){
	    
	    if ($line =~ /^\s*$/){
		# skip blank lines
		next;
	    }
	    
	    if ( $line =~ /^\#\#FASTA/) {
		$fasta_directive_detected = 1;
	    }
	    elsif ( $line =~ /^\#/) {
		next;
		&fetch_gff_comment($gffComment, $line, $linectr);
	    }
	    else {
		&fetch_gff_record($gffrecords, 
				  $line, 
				  $linectr, 
				  $parent_lookup,
				  $organism_attributes_lookup,
				  $cdsSegmentLookup);
	    }
	}
	else  {
	    if ($line =~ /^>(\S+)/) {
		if($fasta_sequence && $fasta_header){
		    if ( exists $gffrecords->{$fasta_header}){
			$gffrecords->{$fasta_header}->{'fasta'} = $fasta_sequence;
		    }
		    else {
			$logger->warn("FASTA header '$fasta_header' does not exist in gffrecords lookup for GFF file '$gff_file'");
		    }
		}

		$fasta_header = &clean_id($1);
		$fasta_sequence = '';

	    }
	    else{
		chomp $line;
		if($line =~ /\w+/){
		    $fasta_sequence .= $line;
		}
	    }
	}
    }

    if($fasta_sequence && $fasta_header){
	#
	#  The last line
	#
	if ( exists $gffrecords->{$fasta_header}){
	    $gffrecords->{$fasta_header}->{'fasta'} = $fasta_sequence;
	}
	else {
	    $logger->warn("FASTA header '$fasta_header' does not exist in gffrecords lookup for GFF file '$gff_file'");
	}
	
    }


    return ($gffrecords, $parent_lookup, $cdsSegmentLookup, $organism_attributes_lookup);
}



#----------------------------------------------------------------------
# createBsmlForPrimarySequence()
#
#----------------------------------------------------------------------
sub createBsmlForPrimarySequence {

    my ($id, $gffrecords, $sequenceIdToBsmlSequenceLookup, $sequenceIdToBsmlFeatureTableLookup,
	$genomeId, $bsmlBuilder) = @_;


    #
    # This should be the first time we attempt to process this Sequence record
    #
    if ( exists $sequenceIdToBsmlSequenceLookup->{$id}){
	$logger->logdie("ID '$id' already exists on sequenceIdToBsmlSequenceLookup- ".
			"meaning this Sequence has already been processed.  Please ".
			"verify that there are not any duplicated records in ".
			"GFF3 file '$gff_file'");
    }
    
    if ( exists $sequenceIdToBsmlFeatureTableLookup->{$id}){
	$logger->logdie("ID '$id' already exists on sequenceIdToBsmlFeatureTableLookup- ".
			"meaning this Sequence has already been processed.  Please ".
			"verify that there are not any duplicated records in ".
			"GFF3 file '$gff_file'");
    }

    my $start  = $gffrecords->{$id}->{'start'};
    my $end    = $gffrecords->{$id}->{'end'};
    my $type   = $gffrecords->{$id}->{'type'};
    my $source = $gffrecords->{$id}->{'source'};
    my $seq_id = $gffrecords->{$id}->{'seq_id'};
    
    my $length = abs($end - $start);
    my $strand = &convert_strand($gffrecords->{$id}->{'strand'});
    
    my $attributes = $gffrecords->{$id}->{'attributes'};

    my $molecule_type = $attributes->{'molecule_type'}->[0];
    my $topology      = $attributes->{'topology'}->[0];
    my $dbsource      = shift (@{$attributes->{'Dbxref'}});
    
    ## Retrieve the new identifier from IdGenerator
    my $newId = &getIdFromIdGenerator($id, $type);

    ## Store the new IdGenerator identifier
    $gffrecords->{$id}->{'nid'} = $newId;

    ## Create a <Sequence> element object
    my $bsmlSequenceElem = $bsmlBuilder->createAndAddExtendedSequenceN(
								       id       => $newId, 
								       title    => $newId,
								       length   => $length,
								       molecule => $molecule_type, 
								       locus    => undef,
								       dbsource => $dbsource,
								       icAcckey => undef,
								       topology => $topology,
								       strand   => $strand,
								       class    => $type
								       );
    
    if (!defined($bsmlSequenceElem)) {
	$logger->logdie("Could not create a <Sequence> element for primary sequence ".
			"with id '$id' newId '$newId' for GFF file '$gff_file'") 
    }
    else {

	## Store reference to this <Sequence> element object
	$sequenceIdToBsmlSequenceLookup->{$id} = $bsmlSequenceElem;

	## Create a <Seq-data-import> element object
	if ( exists ($gffrecords->{$id}->{'fasta'})){
	    my $bsmlSeqDataImportElem = $bsmlBuilder->createAndAddSeqDataImport(
										$bsmlSequenceElem, # Sequence element object reference
										'fasta',           # format
										$fastafile,        # source
										undef,             # id
										$newId                # identifier
										);
	    if (!defined($bsmlSeqDataImportElem)){
		$logger->logdie("Could not create <Seq-data-import> for primary sequence ".
				"with id '$id' newId '$newId' in GFF3 file '$gff_file'");
	    }
	    else{
		&write_fasta_to_file($newId,
				     $gffrecords->{$id}->{'fasta'},
				     $fastafh);
	    }
	}
	else{
	    $logger->logdie("There was not FASTA sequence for the primary ".
			    "sequence with id '$id' newId '$newId' in ".
			    "GFF3 file '$gff_file'");
	}

	
	if (defined($source)){

	    ## Create a <Cross-reference> element object
	    my $bsmlCrossReferenceElem = $bsmlBuilder->createAndAddCrossReference(
										  'parent'    => $bsmlSequenceElem,
										  'id'        => $bsmlBuilder->{'xrefctr'}++,
										  'database'  => $source,
										  'identifier'=> $id,
										  'identifier-type' => 'current'
										  );
	    if (!defined($bsmlCrossReferenceElem)){
		$logger->logdie("Could not create <Cross-reference> for primary Sequence with ".
				"id '$id' newId '$newId' database '$source' identifier '$id' ".
				"identifier-type 'current' for GFF3 file '$gff_file'");
	    }
	}
	
	## Create <Feature-table> element object
	my $bsmlFeatureTableElem = $bsmlBuilder->createAndAddFeatureTable($bsmlSequenceElem);
	
	if (!defined($bsmlFeatureTableElem)){
	    $logger->logdie("Could not create <Feature-table> for primary Sequence with ".
			    "id '$id' newId '$newId' for GFF3 file '$gff_file'");
	}

	## Store reference to this <Feature-table> element object
	$sequenceIdToBsmlFeatureTableLookup->{$id} = $bsmlFeatureTableElem;


	## Link this primary sequence to the Genome
	my $bsmlLinkElem = $bsmlBuilder->createAndAddLink(
							  $bsmlSequenceElem,
							  'genome',       # rel
							  "#$genomeId"    # href
							  );
	
	if (!defined($bsmlLinkElem)){
	    $logger->logdie("Could not create a 'genome' <Link> for primary sequence with ".
			    "id '$id' newId '$newId' for genomeId '$genomeId' - ".
			    "GFF3 file '$gff_file'");
	}

	if (defined($dbsource)){
	    my ($database, $identifier) = &getDatabaseAndIdentifier($dbsource);
	    
	    if (!defined($database)){
		$logger->logdie("database was not defined for dbsource '$dbsource' ".
				"id '$id' newId '$newId' in GFF3 file '$gff_file'");
	    }
	    
	    ## Create <Cross-reference> element object
	    my $bsmlCrossReferenceElem = $bsmlBuilder->createAndAddCrossReference(
										  'parent'          => $bsmlSequenceElem,
										  'id'              => $bsmlBuilder->{'xrefctr'}++,
										  'database'        => $database,
										  'identifier'      => $identifier
										  );
	    if (!defined($bsmlCrossReferenceElem)){
		$logger->logdie("Could not create <Cross-reference> for id '$id' newId '$newId' ".
				"database '$database' identifier '$identifier' ".
				"in GFF3 file '$gff_file'");
	    }							
	}


	if (exists $attributes->{'Dbxref'}){
	    &storeDbxrefAttributesAsBsmlCrossReferences($id,
							$newId,
							$bsmlSequenceElem,
							$bsmlBuilder,
							$attributes->{'Dbxref'});

	    ## Do we really need to ensure that some other code
	    ## will attempt to process the same primary sequence's
	    ## Dbxref attributes?  Oh well, just delete it to be 
	    ## safe.
	    delete $attributes->{'Dbxref'};

	}


	my $bsmlAttributeName = 'molecule_name';

	if ($type ne 'contig'){
	    $bsmlAttributeName = 'gene_name';
	}
	
	&storeAttributeByTypeAsBsmlAttribute($bsmlBuilder,
					     $bsmlSequenceElem,
					     $attributes,
					     $type,
					     'Name',
					     $bsmlAttributeName,
					     $id,
					     $newId);
	

	&store_locus_as_cross_reference($bsmlSequenceElem,
					$attributes,
					$bsmlBuilder,
					$id,
					$newId);


	&storeAttributeByTypeAsBsmlAttribute(  $bsmlBuilder,
					       $bsmlSequenceElem,
					       $attributes,
					       'description',
					       'description',
					       $id,
					       $newId);

	&storeAttributeByTypeAsBsmlAttribute( $bsmlBuilder ,
					      $bsmlSequenceElem,
					      $attributes,
					      'gene_symbol',
					      'gene_symbol',
					      $id,
					      $newId);
	
	
    }
}



#--------------------------------------------------------------
# getIdFromIdGenerator()
#
#--------------------------------------------------------------
sub getIdFromIdGenerator {

    my ($id, $type) = @_;

    my $returnid = $id;

    if ((defined($no_id_generator)) && ($no_id_generator == 1)){
	# not using IdGenerator this round!
	$id_mapping_lookup->{'old2new'}->{$id} = $id;
	$id_mapping_lookup->{'new2old'}->{$id} = $id;
	
    }
    else {

	if ( exists $id_mapping_lookup->{'old2new'}->{$id}) {
	 
	    $returnid = $id_mapping_lookup->{'old2new'}->{$id};
	}
	else {
	    
	    $type = lc($type);

	    $returnid = $idcreator->next_id( project => $project,
					     type    => $type );
	    
	    $id_mapping_lookup->{'old2new'}->{$id} = $returnid;
	    $id_mapping_lookup->{'new2old'}->{$returnid} = $id;
	    
	}
    }
    
    return $returnid;

}

##-----------------------------------------------------------------
## storeDbxrefAttributesAsBsmlCrossReferences()
##
##-----------------------------------------------------------------
sub storeDbxrefAttributesAsBsmlCrossReferences {

    my ($id, $newId, $bsmlSequenceElem, $bsmlBuilder, $dbxrefAttributes) = @_;


    if (scalar(@{$dbxrefAttributes}) > 0 ) {

	foreach my $dbxref ( @{$dbxrefAttributes} ) { 
		
	    my ($database, $identifier) = &getDatabaseAndIdentifer($dbxref);
	    
	    if (!defined($database)){
		$logger->warn("database was not defined for dbxref '$dbxref' ".
			      "id '$id' newId '$newId' in GFF3 file '$gff_file'");
		next;
	    }

	    ## Create <Cross-reference> element object
	    my $bsmlCrossReferenceElem = $bsmlBuilder->createAndAddCrossReference(
										  'parent'          => $bsmlSequenceElem,
										  'id'              => $bsmlBuilder->{'xrefctr'}++,
										  'database'        => $database,
										  'identifier'      => $identifier
									  );
	    if (!defined($bsmlCrossReferenceElem)){
		$logger->logdie("Could not create <Cross-reference> for id '$id' newId '$newId' ".
				"database '$database' identifier '$identifier' ".
				"in GFF3 file '$gff_file'");
	    }							
	}
    }
}	    

#--------------------------------------------------------------------------
# storeAttributeByTypeAsBsmlAttribute()
#
#--------------------------------------------------------------------------
sub storeAttributeByTypeAsBsmlAttribute {

    my ($bsmlBuilder, $bsmlElement, $attributes, $gffAttributeName,
	$bsmlAttributeName, $id, $newId) = @_;
    
    ## We currently do not support multiple values for any
    ## gff attributes - including the following:
    ## 1) Name
    ## 2) locus
    ## 3) description
    ## 4) gene_symbol

    if ( exists $attributes->{$gffAttributeName}){

	my $value = $attributes->{$gffAttributeName}->[0];
	
	if ((defined($value)) &&
	    (length($value) > 0)) {

	    my $bsmlAttributeElem = $bsmlBuilder->createAndAddBsmlAttribute(
									    $bsmlElement,
									    $bsmlAttributeName,
									    $value
									    );
	    
	    if (!defined($bsmlAttributeElem)){
		$logger->logdie("Could not create <Attribute> for name '$bsmlAttributeName' ".
				"content '$value' for id '$id' newId '$newId' in ".
				"GFF3 file '$gff_file'");
	    }
	}
	else {
	    $logger->warn("No value for GFF attribute '$gffAttributeName ".
			  "for id '$id' newId '$newId' in GFF3 file '$gff_file'");
	}
	
	## Remove the gff record attribute from the attributes hash
	delete $attributes->{$gffAttributeName};
    }
}

##--------------------------------------------------
## getDatabaseAndIdentifier()
##
##--------------------------------------------------
sub getDatabaseAndIdentifier {

    my ($data) = @_;

    my $database;
    my $identifier;

    if (defined($data)){
	if ($data =~ /:/){
	    ($database, $identifier) = split(/:/, $data);
	}
    }

    return ($database, $identifier);
}



#-------------------------------------------------------------------
# storeGffRecordsInBsml()
#
#-------------------------------------------------------------------
sub storeGffRecordsInBsml {
    
    my ($primarySequences, 
	$gffrecords,
	$bsmlBuilder,
	$organism_lookup, 
	$sequenceIdToBsmlFeatureTableLookup,
	$sequenceIdToBsmlSequenceLookup, 
	$parent_lookup) = @_;

    print "Processing GFF records\n";
    
    foreach my $id (sort keys %{$gffrecords}) {

	if (( exists $primarySequences->{'lookup'}->{$id}) &&
	    ( defined ($primarySequences->{'lookup'}->{$id})) ) {
	    $logger->info("This primary sequence id '$id' newId '$gffrecords->{$id}->{'nid'}' ".
			  "was already processed");
	    next;
	}

	my $seq_id = $gffrecords->{$id}->{'seq_id'};

	if (!defined($seq_id)){
	    $logger->logdie("seq_id was not defined for this GFF record id '$id' ".
			    "in GFF3 file '$gff_file'");
	}
	
	if ( ! exists $sequenceIdToBsmlFeatureTableLookup->{$seq_id} ){
	    $logger->logdie("No <Feature-table> exists for Feature with id '$id' ".
			    "seq_id '$seq_id' in GFF3 file '$gff_file'");
	}

	my $start  = $gffrecords->{$id}->{'start'};
	my $end    = $gffrecords->{$id}->{'end'};
	my $strand = $gffrecords->{$id}->{'strand'};
	my $type   = $gffrecords->{$id}->{'type'};
	
	my $attributes = $gffrecords->{$id}->{'attributes'};

	my $dbsource = shift (@{$attributes->{'Dbxref'}});

	my $bsmlFeatureTableElem = $sequenceIdToBsmlFeatureTableLookup->{$seq_id};
	
	my $complement = &get_complement($strand, $id);
	
	my $fmin = $start - 1;
	my $fmax = $end;

	## Retrieve the new identifier from IdGenerator
	my $newId = &getIdFromIdGenerator($id, $type);
	
	## Store the new IdGenerator identifier
	$gffrecords->{$id}->{'nid'} = $newId;

	#
	# Create <Feature> element object
	#
	my $bsmlFeatureElem = $bsmlBuilder->createAndAddFeatureWithLoc(
								       $bsmlFeatureTableElem,  # <Feature-table> element object reference
								       $newId,               # id
								       undef,                # title
								       $type,                # class
								       undef,                # comment
								       undef,                # displayAuto
								       $fmin,                # start
								       $fmax,                # stop
								       $complement           # complement
								       );
	
	if (!defined($bsmlFeatureElem)){
	    $logger->logdie("Could not create <Feature> for id '$id' newId '$newId' ".
			    "in GFF3 file '$gff_file'"); 
	}
	else {

	    if (lc($type) eq 'cds'){

		
		my $contigFasta = $gffrecords->{$seq_id}->{'fasta'};
		
		if (!defined($contigFasta)){
		    $logger->logdie("FASTA sequence was not defined for contig '$seq_id' ".
				    "while processing id '$id' in GFF file '$gff_file'");
		}
		
		my $length = $end - $start + 1;

		my $derivedCdsFastaSequence = substr($contigFasta, $start, $length);

		## We will attempt to create the following BSML elements:
		## <Sequence>
		## <Seq-data-import>
		## <Link>
		&createSequenceForFeatureWithFasta($id,
						   $newId,
						   $bsmlFeatureElem,
						   $derivedCdsFastaSequence,
						   $dbsource,
						   $gffrecords->{$id}->{'topology'},
						   $strand,
						   $fmax - $fmin,
						   $gffrecords->{$id}->{'attributes'}->{'molecule'},
						   $bsmlBuilder,
						   $type
						   );
	    }
	    elsif ( exists $gffrecords->{$id}->{'fasta'} ){

		## We will attempt to create the following BSML elements:
		## <Sequence>
		## <Seq-data-import>
		## <Link>
		&createSequenceForFeatureWithFasta($id,
						   $newId,
						   $bsmlFeatureElem,
						   $gffrecords->{$id}->{'fasta'},
						   $dbsource,
						   $gffrecords->{$id}->{'topology'},
						   $strand,
						   $fmax - $fmin,
						   $gffrecords->{$id}->{'attributes'}->{'molecule'},
						   $bsmlBuilder,
						   $type
						   );
	    }


	    #-----------------------------------------------------------------------------------------------------
	    # Create 'current' <Cross-reference> object
	    #
	    #-----------------------------------------------------------------------------------------------------
	    my $xref_elem = $bsmlBuilder->createAndAddCrossReference( 'parent'          => $bsmlFeatureElem,
								      'id'              => $bsmlBuilder->{'xrefctr'}++,
								      'database'        => $gffrecords->{$id}->{'source'},
								      'identifier'      => $id,
								      'identifier-type' => 'current'
								      );
	    if (!defined($xref_elem)){
		$logger->logdie("Could not create <Cross-reference> for Feature with id '$id'  ".
				"newId '$newId' database '$gffrecords->{$id}->{'source'}' ".
				"identifier '$id' identifier-type 'current' in GFF3 file '$gff_file'");
	    }
	    
	    #-----------------------------------------------------------------------------------------------------
	    # Create 'stable_id' <Cross-reference> objects
	    #
	    #-----------------------------------------------------------------------------------------------------
	    if ( exists $attributes->{'stable_id'} ) {

		foreach my $stable_id ( sort @{$attributes->{'stable_id'}} ) {
		    
		    my $xref = $bsmlBuilder->createAndAddCrossReference( 'parent'          => $bsmlFeatureElem,
								 'id'              => $bsmlBuilder->{'xrefctr'}++,
								 'database'        => $gffrecords->{$id}->{'source'},
								 'identifier'      => $stable_id,
								 'identifier-type' => 'stable_id'
								 );
		    if (!defined($xref)){
			$logger->logdie("Could not create <Cross-reference> Feature with id '$id' ".
					"newId '$newId' database '$gffrecords->{$id}->{'source'}' ".
					"identifier '$stable_id' identifier-type 'stable_id' for ".
					"GFF3 file '$gff_file'");
		    }
		}
	    }

	    #-----------------------------------------------------------------------------------------------------
	    # Create 'primary' <Cross-reference> object
	    #
	    #-----------------------------------------------------------------------------------------------------
	    if (defined($dbsource)) {			    

		if ($dbsource =~ /:/){
		    
		    my ($database, $identifier) = split(/:/, $dbsource);
		    
		    my $xref = $bsmlBuilder->createAndAddCrossReference( 'parent'          => $bsmlFeatureElem,
								 'id'              => $bsmlBuilder->{'xrefctr'}++,
								 'database'        => $database,
								 'identifier'      => $identifier
								 );
		    if (!defined($xref)){
			$logger->logdie("Could not create <Cross-reference> for Feature with id '$id' ".
					"newId '$newId' database '$database' identifier '$identifier' ".
					"identifier-type 'primary' in GFF3 file '$gff_file'");
		    }
		}
		else {
		    $logger->logdie("Could not derive database, identifier from '$dbsource' for ID '$id' ".
				    "newId '$newId' in GFF3 file '$gff_file'");
		}
	    }
	    
	    #-----------------------------------------------------------------------------------------------------
	    # Process the 'Alias' attribute
	    #
	    #-----------------------------------------------------------------------------------------------------
	    if (exists $attributes->{'Alias'}) {

		foreach my $alias ( @{$attributes->{'Alias'}} ) { 
		    
		    if ($alias =~ /:/){
			
			## Store the Alias as a <Cross-reference>

			my ($database, $identifier) = split(/:/, $alias);
			
			my $xref = $bsmlBuilder->createAndAddCrossReference( 'parent'          => $bsmlFeatureElem,
									     'id'              => $bsmlBuilder->{'xrefctr'}++,
									     'database'        => $database,
									     'identifier'      => $identifier,
									     'identifier-type' => 'alias'
								     );
			if (!defined($xref)){
			    $logger->logdie("Could not create <Cross-reference> for Feature with id '$id' ".
					    "newId '$newId' database '$database' identifier '$identifier' ".
					    "identifier-type 'Alias' for GFF3 file '$gff_file'");
			}
		    }
		    else {
			## Store the Alias as an <Attribute>
			my $attribute_elem = $bsmlBuilder->createAndAddBsmlAttribute(
									     $bsmlFeatureElem,
									     'alias',
									     $alias
									     );
			
			if (!defined($attribute_elem)){
			    $logger->logdie("Could not create <Attribute> for name 'alias' content '$alias' ".
					    "for Feature with id '$id' newId '$newId' in GFF3 file '$gff_file'");
			}

		    }
		}
		
		delete $attributes->{'Alias'};
	    }

	    #-----------------------------------------------------------------------------------------------------
	    # Create 'Dbxref' <Cross-reference> objects
	    #
	    #-----------------------------------------------------------------------------------------------------
	    if (exists $attributes->{'Dbxref'} ){

		foreach my $dbxref ( @{$attributes->{'Dbxref'}} ) { 
		    
		    if ($dbxref =~ /:/){
			
			my ($database, $identifier) = split(/:/, $dbxref);
			
			my $xref = $bsmlBuilder->createAndAddCrossReference( 'parent'          => $bsmlFeatureElem,
								     'id'              => $bsmlBuilder->{'xrefctr'}++,
								     'database'        => $database,
								     'identifier'      => $identifier
								     );
			if (!defined($xref)){
			    $logger->logdie("Could not create <Cross-reference> element object for Feature ".
					    "element ID '$id' database '$database' identifier '$identifier' ".
					    "identifier-type 'Dbxref' for GFF3 file '$gff_file'");
			}
		    }
		    else {
			$logger->logdie("Could not derive database, identifier from Dbxref '$dbxref' for ID '$id'");
		    }
		}
		
		delete $attributes->{'Dbxref'};
	    }


	    #-----------------------------------------------------------------------------------------------------
	    # Create 'locus' <Cross-reference> objects
	    #
	    #-----------------------------------------------------------------------------------------------------
	    if (exists $attributes->{'locus'} ){
		
		foreach my $locus (@{$attributes->{'locus'}} ){
		    
		    my $xref = $bsmlBuilder->createAndAddCrossReference( 'parent'          => $bsmlFeatureElem,
								 'id'              => $bsmlBuilder->{'xrefctr'}++,
								 'database'        => $gffrecords->{$id}->{'source'},
								 'identifier'      => $locus,
								 'identifier-type' => 'locus'
								 );
		    if (!defined($xref)){
			$logger->logdie("Could not create <Cross-reference> for Feature with id '$id' ".
					"newId '$newId' database '$gffrecords->{$id}->{'source'}' ".
					"identifier '$locus' identifier-type 'locus' in GFF3 file '$gff_file'");
		    }
		}
	    }
	    
	    &storeAttributeByTypeAsBsmlAttribute($bsmlBuilder,
						 $bsmlFeatureElem,
						 $attributes,
						 $type,
						 'Name',
						 'gene_name',
						 $id,
						 $newId);
	    
	    &storeAttributeByTypeAsBsmlAttribute($bsmlBuilder,
						 $bsmlFeatureElem,
						 $attributes,
						 $type,
						 'description',
						 'gene_product_name',
						 $id,
						 $newId);
	    
	    &storeAttributeByTypeAsBsmlAttribute($bsmlBuilder,
						 $bsmlFeatureElem,
						 $attributes,
						 $type,
						 'gene_symbol',
						 'gene_symbol',
						 $id,
						 $newId);


	    #
	    # Store the EC number in an <Attribute-list> element
	    #
	    &store_gff_attributes_as_attribute_list( $bsmlFeatureElem,
						     $attributes,
						     'ec_number');
	    
	    #
	    # Store the GO assignments in an <Attribute-list> element
	    #
	    &store_gff_attributes_as_attribute_list( $bsmlFeatureElem,
						     $attributes,
						     'Ontology_term');
	    

	}
    }
}


##--------------------------------------------------------
## createSequenceForFeatureWithFasta()
##
##--------------------------------------------------------
sub createSequenceForFeatureWithFasta {

    my ($id, $newId, $bsmlFeatureElem, $fasta,$dbsource, 
	$topology, $strand, $length, $molecule, $bsmlBuilder,
	$type) = @_;

    my $seqId = $newId . '_seq';

    ## Create a BSML <Sequence> element object
    my $bsmlSequenceElem = $bsmlBuilder->createAndAddExtendedSequenceN(
								       id       => $seqId, 
								       title    => $seqId,
								       length   => $length,
								       molecule => $molecule, 
								       locus    => undef,
								       dbsource => $dbsource,
								       icAcckey => undef,
								       topology => $topology,
								       strand   => $strand,
								       class    => $type
								       );
    
    if (!defined($bsmlSequenceElem)) {
	$logger->logdie("Could not create a <Sequence> for Feature with id '$id' ".
			"newId '$newId' in GFF3 file '$gff_file'") 
    }

    ## Create <Seq-data-import> elemen object
    my $bsmlSeqDataImportElem = $bsmlBuilder->createAndAddSeqDataImport(
									$bsmlSequenceElem, # Sequence element object reference
									'fasta',           # format
									$fastafile,        # source
									undef,             # id
									$newId                # identifier
									);
    if (!defined($bsmlSeqDataImportElem)){
	$logger->logdie("Could not create <Seq-data-import> for primary sequence ".
			"with id '$id' newId '$newId' in GFF3 file '$gff_file'");
    }
    else{
	&write_fasta_to_file($newId,
			     $fasta,
			     $fastafh);
    }



    ## Create a BSML <Link> element object
    my $bsmlLinkElem = $bsmlBuilder->createAndAddLink(
						      $bsmlFeatureElem, # element object reference
						      'sequence',       # rel
						      "#$seqId"         # href
						      );
	
    if (!defined($bsmlLinkElem)){
	$logger->logdie("Could not create a sequence <Link> for Feature with id '$id' ".
			"newId '$newId' in GFF3 file '$gff_file'");
    }
}


##------------------------------------------------------------------------------------------------
## storeFeatureAsFeatureGroupMember()
##
## Currently, the code can only create a Feature-group-member
## element objects for gff3 records which have a "Parent"
## attribute in their attribute lists.
##
##------------------------------------------------------------------------------------------------
sub storeFeatureAsFeatureGroupMember {

    my ($parent_lookup,
	$bsmlBuilder, 
	$sequenceIdToBsmlSequenceLookup,
	$gffrecords) = @_;

    print "Storing Feature-group data\n";

    ## Retain reference for all Feature-group elements
    my $featureGroupLookup = {};
    
    ## Because we cannot guarantee when we are going to process
    ## a given Feature
    my $featureGroupMemberLookup = {};

    my $featureGroupMemberByParentLookup = {};

    foreach my $id ( keys %{$gffrecords}){
	
	my $type = $gffrecords->{$id}->{'type'};
	my $newId = $gffrecords->{$id}->{'nid'};

	if ($logger->is_debug()){
	    $logger->debug("id '$id' newId '$newId' type '$type'");
	}
	
	## Does this feature reference some parent Feature/Sequence?
	if (( exists $parent_lookup->{$id}) &&
	    ( defined($parent_lookup->{$id})) &&
	    ( scalar(@{$parent_lookup->{$id}}) > 0 )) {
	    
	    foreach my $parent_id ( @{$parent_lookup->{$id}} ) {
		
		if (!defined($parent_id)){
		    next;
		}
		
		if ($logger->is_debug()){
		    $logger->debug("id '$id' newId '$newId' type '$type' parent_id '$parent_id'");
		}

		my $featureGroupMemberRepresentative;
		
		my $bsmlSequenceElemForPrimarySequence = $sequenceIdToBsmlSequenceLookup->{$gffrecords->{$id}->{'seq_id'}};
		
		
		if (!defined($bsmlSequenceElemForPrimarySequence)){
		    $logger->logdie("Sequence does not exist for seq_id '$gffrecords->{$id}->{'seq_id'}' while processing ".
				    "id '$id' newId '$newId' parent_id '$parent_id' in GFF file '$gff_file'");
		}
		## This Feature's parent is not a primary Sequence
		## so we should Feature-group both together
		
		my $newParentId = $gffrecords->{$parent_id}->{'nid'};
		
		if (! exists $featureGroupLookup->{$parent_id}){
		    ## Create a <Feature-group> for the parent_id
		    if ($logger->is_debug()){
			$logger->debug("<Feature-group> does not exist for parent_id '$parent_id'");
		    }

		    
		    $featureGroupMemberRepresentative = $bsmlBuilder->createAndAddFeatureGroup($bsmlSequenceElemForPrimarySequence,
											       undef,
											       $newParentId);
		    if (!defined($featureGroupMemberRepresentative)){
			$logger->logdie("Could not create <Feature-group> for id '$parent_id' ".
					"newParentId '$newParentId' while processing Feature ".
					"with id '$id' newId '$newId' in ".
					"GFF file '$gff_file'");
		    }
		    
		    if ($logger->is_debug()){
			$logger->debug("Created <Feature-group> for parent_id '$parent_id' with newParentId '$newParentId'");
		    }
		    
		    ##
		    $featureGroupLookup->{$parent_id} = $featureGroupMemberRepresentative;
		    

		    if ($logger->is_debug()){
			$logger->debug("Created Feature-group-member for id '$id' with parent_id '$parent_id' ".
				       "having just had its Feature-group created");
		    }

		    &storeFeatureGroupMember($featureGroupMemberRepresentative,
					     $id, 
					     $newId,
					     $type,
					     $featureGroupMemberByParentLookup,
					     $bsmlBuilder,
					     $parent_id);

		}
		else {
		    ## <Feature-group-member> already exists for parent_id
		    if ($logger->is_debug()){
			$logger->debug("Feature-group-member already exists for parent_id '$parent_id' newParentId ".
				       "'$newParentId' while processing id '$id' newId '$newId'");
		    }
		    $featureGroupMemberRepresentative = $featureGroupLookup->{$parent_id};
		}


		if (!exists $featureGroupMemberByParentLookup->{$parent_id}->{$parent_id}){
		    ## Create a <Feature-group-member> for the parent_id
		    my $parentType = $gffrecords->{$parent_id}->{'type'};
		    
		    my $featureGroupMemberElem = $bsmlBuilder->createAndAddFeatureGroupMember($featureGroupMemberRepresentative, # <Feature-group>
											      $newParentId, # featref
											      $parentType,  # feattype
											      undef,      # grouptype
											      undef       # cdata
											      );
		    if (!defined($featureGroupMemberElem)){
			$logger->logdie("Could not create <Feature-group-member> for Feature parent_id '$parent_id' ".
					"newParentId '$newParentId' while processing id '$id' newId '$newId' ".
					"in GFF file '$gff_file'");
		    }
		    
		    if ($logger->is_debug()){
			$logger->debug("Created Feature-group-member for parent_id '$parent_id'");
		    }

		    ## Store reference for this <Feature-group-member>
		    $featureGroupMemberByParentLookup->{$parent_id}->{$parent_id} = $featureGroupMemberElem;
		}
		
		if (! exists $featureGroupMemberByParentLookup->{$parent_id}->{$id}){

		    if ($logger->is_debug()){
			$logger->debug("<Feature-group-member> does not exist for id '$id' newId '$newId'");
		    }

		    if (defined($featureGroupMemberRepresentative)){

			if ($logger->is_debug()){
			    $logger->debug("Feature-group-member for representative of id '$id' newId '$newId' was defined");
			}

			&storeFeatureGroupMember($featureGroupMemberRepresentative,
						 $id, 
						 $newId,
						 $type,
						 $featureGroupMemberByParentLookup,
						 $bsmlBuilder,
						 $parent_id);

		    }
		    else {
			$logger->logdie("<Feature-group> was not defined for parent_id '$parent_id' while ".
					"processing Feature with id '$id' in GFF file '$gff_file'");
		    }
		}
		else {
		    if ($logger->is_debug()){
			$logger->debug("Feature with id '$id' newId '$newId' ".
				       "already has a <Feature-group-member>");
		    }
		}
	    }
	}
    }
}

##----------------------------------------------------------------------------------
## storeFeatureGroupMember()
##
##----------------------------------------------------------------------------------
sub storeFeatureGroupMember {

    my ($featureGroupMemberRepresentative,
	$id,
	$newId,
	$type,
	$featureGroupMemberByParentLookup,
	$bsmlBuilder,
	$parent_id) = @_;

    my $featureGroupMemberElem = $bsmlBuilder->createAndAddFeatureGroupMember($featureGroupMemberRepresentative, # <Feature-group>
									      $newId,        # featref
									      $type,      # feattype
									      undef,      # grouptype
									      undef       # cdata
									      );
    if (!defined($featureGroupMemberElem)){
	$logger->logdie("Could not create <Feature-group-member> for Feature id '$id' ".
			"newId '$newId' with parent_id '$parent_id' in GFF file '$gff_file'");
    }
    
    if ($logger->is_debug()){
	$logger->debug("Created <Feature-group-member> for id '$id' newId '$newId' type '$type' ".
		       "parent_id '$parent_id'");
    }
    
    
    ## Store reference for this <Feature-group-member>
    $featureGroupMemberByParentLookup->{$parent_id}->{$id} = $featureGroupMemberElem;
}
