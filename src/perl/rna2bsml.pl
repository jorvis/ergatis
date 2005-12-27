#!/usr/local/bin/perl

=head1 NAME

rna2bsml - Produces RNA/repeat BSML documents

=head1 SYNOPSIS

USAGE:  rna2bsml -D database -G genus -S species -f infile [-d debug_level] [-h] [-l log4perl] [-o outdir] [-s schema] [-t dtd] -y orgtype -z ftype

=head1 OPTIONS

=over 8

=item B<--database,-D>
    
    Legacy organism database name e.g. "tba1"

=item B<--genus,-G>
    
    Legacy organism database name e.g. "tba1"

=item B<--infile,-f> 

    Tab-delimited file containing the feature data

=item B<--debug_level,-d> 

    Optional -- Coati::Logger log4perl logging level.  Default debug_level is 0

=item B<--help,-h>

    Print this help

=item B<--log4perl,-l>
    
    Optional -- Coati::Logger log4perl log file.  Default is /tmp/rna2bsml.pl.log

=item B<--outdir,-o> 

    Optional -- Directory to write BSML gene model documents.  Default directory is current working directory

=item B<--schema,-s>

    Optional -- Schema to validate against. If blank validate to inline schema

=item B<--dtd,-t>

    Optional -- DTD validation via the DTD validator

=item B<--orgtype,-y>

    organism type e.g. euk, prok or ntprok

=item B<--ftype,-z>

    feature type either 'rna', 'te', or 'repeat'

=back

=head1 DESCRIPTION

    rna2bsml creates RNA/repeat BSML documents
    
=cut

use lib "shared";
use strict;
use Log::Log4perl qw(get_logger :levels :easy);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Data::Dumper;
BEGIN {
    require '/usr/local/devel/ANNOTATION/ard/chado-v1r5b1/lib/site_perl/5.8.5/Workflow/Logger.pm';
    import Workflow::Logger;
    require '/usr/local/devel/ANNOTATION/ard/chado-v1r5b1/lib/site_perl/5.8.5/BSML/BsmlBuilder.pm';
    import BSML::BsmlBuilder;
}

#
# Process command line options
#
my ($database, $species, $genus, $infile, $debug_level, $log4perl ,$help, $outdir, $schema, $dtd, $man, $orgtype, $ftype);

my $results = GetOptions (
			  'database|D=s'     => \$database, 
			  'outdir|o=s'       => \$outdir, 
			  'species|S=s'      => \$species,
			  'genus|G=s'        => \$genus,
			  'schema|s=s'       => \$schema,
			  'dtd|t=s'          => \$dtd,
			  'debug_level|d=s'  => \$debug_level,
			  'log4perl|l=s'     => \$log4perl,
			  'help|h'           => \$help,
			  'infile|f=s'       => \$infile,
			  'man|m'            => \$man,
			  'orgtype|y=s'      => \$orgtype,
			  'ftype|z=s'        => \$ftype
			  );

print STDERR "database was not defined\n" if (!defined($database));
print STDERR "species was not defined\n"  if (!defined($species));
print STDERR "genus was not defined\n"    if (!defined($genus));
print STDERR "infile was not defined\n"   if (!defined($infile));
print STDERR "orgtype was not defined\n"  if (!defined($orgtype));
print STDERR "ftype was not defined\n"    if (!defined($ftype));


&print_usage() if (!defined($database) || !defined($species) || !defined($genus) || !defined($infile) || !defined($orgtype));

&print_usage() if ($help);
&pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if($man);


#
# initialize the logger
#
$log4perl = "/tmp/rna2bsml.pl.log" if (!defined($log4perl));
my $mylogger = new Coati::Logger('LOG_FILE'=>$log4perl,
				 'LOG_LEVEL'=>$debug_level);

my $logger = Coati::Logger::get_logger(__PACKAGE__);


if (($ftype eq 'rna') or ($ftype eq 'te') or ($ftype eq 'repeat')){
    $logger->info("Processing features of type '$ftype'");
}
else{
    $logger->logdie("Unrecognized feature type '$ftype'");
}


$outdir = &verify_and_set_outdir($outdir);

$logger->info("database set to '$database'");

&construct_documents(
		     file     => $infile,
		     database => $database,
		     species  => $species,
		     genus    => $genus,
		     ftype    => $ftype     
		     );


$logger->info("Please verify log4perl log file: $log4perl");






#------------------------------------------------------------------------------------------------------------------------
#
#   END OF MAIN -- SUBROUTINES FOLLOW
#
#------------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------
# construct_documents()
#
# This section uses BSML API writer functions to create a 
# BSML gene model document
#
#------------------------------------------------------------------
sub construct_documents {


    $logger->debug("Entered construct_documents") if $logger->is_debug;

    
    my (%args) = @_;
    my $file     = $args{'file'};
    my $database = $args{'database'};
    my $genus    = $args{'genus'};
    my $species  = $args{'species'};
    my $ftype    = $args{'ftype'};

    
    if (!-e $file){
	$logger->logdie("file '$file' does not exist");
    }
    if (!-r $file){
	$logger->logdie("file '$file' does not have read permissions");
    }

    open (INFILE, "<$file") or $logger->logdie("Could not open file '$file'");
    my @contents = <INFILE>;
    chomp @contents;


    my $linectr=0;
    my $asmblhash = {};


    my $regexp;
    if ($ftype eq 'rna'){
	$regexp = '^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S{1,3}RNA)\s+([\S\s]+)';
    }
    elsif ($ftype eq 'te'){
	$regexp = '^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(TE)\s+([\S\s]+)'
    }
    elsif ($ftype eq 'repeat'){
	$regexp = '^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(repeat)\s+([\S\s]+)'
    }
    else{
	$logger->logdie("unrecognized feature type '$ftype'");
    }



    foreach my $line (@contents){
	
	$linectr++;
	
	my ($locus, $feat_name, $asmbl_id, $end5, $end3, $feat_type, $product_name);

	#
	# For RNAs
	#
#	if ($line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S{1,3}RNA)\s+([\S\s]+)/){

	#
	# For tca1 transposable elements
#	if ($line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(TE)\s+([\S\s]+)/){



	if ($line =~ /$regexp/){


	    $locus     = $1;
	    $feat_name = $2;
	    $asmbl_id  = $3;
	    $end5      = $4;
	    $end3      = $5;
	    $feat_type = $6;
	    $product_name = $7;
	


	    #----------------------------------------------------------------
	    # Transform invalid terms:
	    #
	    # TIGR terms        ==>  Sequence Ontology qualified term
	    # slRNA             ==>  splice_leader_RNA
	    # repeated sequence ==> repeat_region
	    #
	    #
	    #----------------------------------------------------------------
	    if ($feat_type eq 'slRNA'){
		$feat_type = 'spliced_leader_RNA';
	    }
	    if ($feat_type eq 'repeated sequence'){
		$feat_type = 'repeat_region';
	    }





	    my $hash;
	    $hash->{'feat_name'} = $feat_name;
	    $hash->{'locus'}     = $locus;
	    $hash->{'end5'}      = $end5;
	    $hash->{'end3'}      = $end3;
	    $hash->{'feat_type'} = $feat_type;
	    $hash->{'product_name'} = $product_name;


	    my $asmbl_id = $database . '_'.  $asmbl_id . '_assembly';

	    push (@{$asmblhash->{$asmbl_id}}, $hash);

	}
	else{
	    $logger->logdie("Regular expression '$regexp' could not parse line '$linectr': $line");
	}
    }


    my $crossrefctr=0;
    my $xref_database = 'TIGR_' . $orgtype . ':'. $database;

    my $asmblctr ={};

    #
    # Create a new BSML document for each assembly
    #
    foreach my $asmbl_id (sort keys %{$asmblhash}){


	$asmblctr->{'asmblcount'}++;


	my $doc = new BSML::BsmlBuilder(); 
	$logger->logdie("doc was not defined") if (!defined($doc));
	


	#
	# Add <Genomes> element
	#
	
	my $genome = $doc->createAndAddGenome();
	
	
	my $identifier = lc(substr($genus,0,1)) . '_' . lc($species);
	
	
	my $xref = $doc->createAndAddCrossReference(
						    'parent'          => $genome,
						    'id'              => '_' . ++$crossrefctr,
						    'database'        => $xref_database,
						    'identifier'      => $identifier,
						    'identifier-type' => 'current'
						    );
	
	
	my $organism = $doc->createAndAddOrganism( 
						   'genome'  => $genome,
						   'genus'   => $genus,  
						   'species' => $species,
						   );
	
	
	#
	# Add the <Sequence> elements for the assemblies
	#
	
	
	my $asmseq = $doc->createAndAddExtendedSequenceN( 
							  'id'       => $asmbl_id, 
							  'title'    => $asmbl_id, 
							  'molecule' => 'dna', 
							  );
	$asmseq->addattr("class", "assembly");
	
	
	#
	# Creating a feature-table object which will eventually result in the
	# creation of a <Feature-tables> BSML element for the assembly
	#
	my $asmFTable = $doc->createAndAddFeatureTable($asmseq);
	$logger->logdie("asmFTable was not defined") if (!defined($asmFTable));

	#
	# For this assembly, add all rna features as <Feature> elements
	#

	my $FGroup;

	foreach my $feathash (sort @{$asmblhash->{$asmbl_id}}){
	    

	    $asmblctr->{$asmbl_id}->{'featctr'}++;

	    my $end5         = $feathash->{'end5'};
	    my $end3         = $feathash->{'end3'};
	    my $locus        = $feathash->{'locus'};
	    my $feat_type    = $feathash->{'feat_type'};
	    my $feat_name    = $feathash->{'feat_name'};
	    my $product_name = $feathash->{'product_name'};


	    if ($feat_type eq 'TE'){
		$feat_type = 'transposable_element';
	    }


	    #
	    # Chado stores stand as 1 = forward strand
	    #                      -1 = reverse strand
	    #
	    # BSML stores as "is complement" 0 = false
	    #                                1 = true
	    #
	    my $complement;
	    if ($end3 > $end5){
		$complement = 0;
	    }
	    elsif ($end5 > $end3){
		$complement = 1;
		($end5, $end3) = ($end3, $end5);
	    }
	    else{
		$logger->logdie("end3 '$end3' = end5 '$end5'");
	    }

	    #
	    # Convert to space-based coordinate
	    #
	    $end5--;

	    
	    #
	    # Create feature group
	    #
	    $FGroup = $doc->createAndAddFeatureGroup(
						     $asmseq,
						     '',
						     "$feat_name"
						     );  

	    my $Feat = $doc->createAndAddFeatureWithLoc(
							$asmFTable,
							"$feat_name",
							'',
							'$feat_type',
							'',
							'',
							$end5,
							$end3,
							$complement
							);


	    my $xref = $doc->createAndAddCrossReference(
							'parent'          => $Feat,
							'id'              => ++$crossrefctr,
							'database'        => $xref_database,
							'identifier'      => $feat_name,
							'identifier-type' => 'current'
							);


	    $logger->logdie("Feat was not defined") if (!defined($Feat));
	    

	    $Feat->addattr('class', "$feat_type");

	    if ($locus ne 'NULL'){
		$doc->createAndAddBsmlAttribute(
						$Feat, 
						'locus', 
						"$locus"
						);
	    }

	    if ($product_name ne 'NULL'){

		$doc->createAndAddBsmlAttribute(
						$Feat, 
						'gene product name', 
						"$product_name"
						);
	    }

	}

	$logger->info("Writing BSML document $outdir/$asmbl_id.bsml");
	
	$doc->write("$outdir/$asmbl_id.bsml");
	
	if(! -e "$outdir/$asmbl_id.bsml"){
	    $logger->error("File not created $outdir/$asmbl_id.bsml");
	}
	
	if($dtd){
	    my $dtdvalid = `$ENV{DTDVALID} $outdir/$asmbl_id.bsml --dtd $dtd`;
	    if($dtdvalid ne '0'){
		$logger->warn("DTD validation failed for $outdir/$asmbl_id.bsml");
	    }
	    elsif($dtdvalid eq '0') {
		$logger->info("DTD validation passed for $outdir/$asmbl_id}.bsml");
	    }
	}
	if($schema){
	    my $schemavalid = `$ENV{SCHEMAVALID} $outdir/$asmbl_id.bsml --schema $schema`;
	    if ($schemavalid ne '0'){
		$logger->warn("XML Schema validation failed for $outdir/$asmbl_id.bsml");
	    }
	    elsif ($schemavalid eq '0'){
		$logger->info("XML schema validation passed for $outdir/$asmbl_id.bsml");
	    }
	}
	chmod 0777, "$outdir/$asmbl_id.bsml";
    }


    $logger->info("Assemblies counted: $asmblctr->{'asmblcount'}");
    $logger->info("Number of features counted per assembly\n");
    foreach my $assembly (keys %{$asmblctr}){
	next if ($assembly eq 'asmblcount');
	$logger->info("assembly '$assembly' feature count '$asmblctr->{$assembly}->{'featctr'}'");
    }



}

#--------------------------------------------------------
# verify_and_set_outdir()
#
#
#--------------------------------------------------------
sub verify_and_set_outdir {

    my ( $outdir, $option) = @_;

    $logger->debug("Verifying and setting output directory") if ($logger->is_debug());

    #
    # strip trailing forward slashes
    #
    $outdir =~ s/\/+$//;
    
    #
    # set to current directory if not defined
    #
     if (!defined($outdir)){
	if (!defined($option)){
	    $outdir = "." 
	}
	else{
	    $outdir = $option;
	}
    }

    #
    # verify whether outdir is in fact a directory
    #
    $logger->fatal("$outdir is not a directory") if (!-d $outdir);

    #
    # verify whether outdir has write permissions
    #
    $logger->fatal("$outdir does not have write permissions") if ((-e $outdir) and (!-w $outdir));


    $logger->debug("outdir is set to:$outdir") if ($logger->is_debug());

    #
    # store the outdir in the environment variable
    #
    
    return $outdir;
    
}#end sub verify_and_set_outdir()


#----------------------------------------------------------------
# print_usage()
#
#----------------------------------------------------------------
sub print_usage {

    print STDERR "\nSAMPLE USAGE:  $0 -D database -G genus -S species -f infile [-d debug_level] [-h] [-l log4perl] [-o outdir] [-s schema] [-t dtd] -y orgtype -z ftype\n";
    print STDERR "  -D|--database         = Legacy organism database name\n";
    print STDERR "  -G|--genus            = genus\n";
    print STDERR "  -S|--species          = species\n";
    print STDERR "  -f|--infile           = Tab-delimited file containing feature data\n";
    print STDERR "  -d|--debug_level      = Optional - Coati::Logger log4perl logging level.  Default debug level is 0\n";
    print STDERR "  -h|--help             = Print this usage help statement\n";
    print STDERR "  -l|--log4perl         = Optional - Coati::Logger log4perl log file.  Default is /tmp/db2bsml.pl.log\n";
    print STDERR "  -o|--outdir           = Optional - Directory to write the BSML gene model documents to. Default directory is current working directory\n";
    print STDERR "  -s|--schema           = Optional - schema validation via XML schema validator\n";
    print STDERR "  -t|--dtd              = Optional - DTD validation via DTD validator\n";
    print STDERR "  -y|--orgtype          = organism type e.g. euk, prok or ntprok\n";
    print STDERR "  -z|--ftype            = feature type either 'rna' or 'te' or 'repeat'\n";
    exit 1;

}
