#!/local/perl/bin/perl
#-------------------------------------------------------------------------------------------
#
#
#
#
# modification:   2003-12-08 Jay Sundaram sundaram@tigr.org 301-610-5983
#                 1) Facilitation of storage of configuration file key-value pairs in the
#                 <Analysis> component, destined for the chado.analysisprop table
#                 2) Introduced log4perl logging
#
# modification:   2003-12-09 Sundaram
#                 1) -l clusterDir parameter.  If the script finds a file "no_cluster.txt"
#                 in the clusterDir, will produce an empty BSML multiple alignment document.
#                 If clusterDir is not defined, script behaves as previous version.
#
#                 2003-12-30 sundaram
#                 1) print_usage() added
#
#
#-------------------------------------------------------------------------------------------

=head1  NAME 

MSF2Bsml.pl  - convert MSF format alignment files into BSML documents

=head1 SYNOPSIS

USAGE:  MSF2Bsml.pl -a ali_dir -o output [-c configfile] [-l log4perl] [-p program] [-s suffix] [-h help] [-k clusterDir]

=head1 OPTIONS

=over 4

=item *

B<--ali_dir,-a>   [REQUIRED] Dir containing MSF format alignment files that ends in *.msf

=item *

B<--output,-o> [REQUIRED] output BSML file

=item *

B<--suffix,-s> suffix of the multiple sequence alignment file. Default(msf)

=item *

B<--program,-p> projects this analysis is for. Default(clustal)

=item *

B<--configfile,-f> Configuration file which will contain key=value pairs for storage in the <Analysis> component of the BSML document

=item *

B<--log4perl,-l> Log4perl log file (default is msf2bsml.log)

=item *

B<--help,-h> This help message

=item *

B<--clusterDir,-k> Directory that will contain the "no_cluster.txt" file if no clusters were generated

=back

=head1   DESCRIPTION

MSF2Bsml.pl is designed to convert multiple sequence alignments in MSF format
into BSML documents.  As an option, the user can specify the suffix of the 
alignment files that the script should parse in the directory specified by 
--ali_dir.  The default suffix is "msf" i.e. family.msf.  In addition, the
user can specify what program this script is for.  For example if the
multiple sequence alignment represents COGS, the user can specify that
with the --program flag.  This will be reflected in the analysis object
of the BSML document. 

Users can specify a configuration file which will contain key-value pairs.  These
values will be written to the <Analysis> component of the BSML document.  Loader
pipelines will parse and load these values into the chado.analysisprop table.

Samples:

1. Convert all MFS formatted files that end in "aln" in ali_dir
   MSF2Bsml.pl -a ali_dir -o alignment.bsml -s aln

2. Convert all MFS formatted files that end in "aln" in ali_dir
   MSF2Bsml.pl -a ali_dir -o alignment.bsml -s aln -c config.txt -l my.log -k /tmp/cluster_dir


NOTE:  

Calling the script name with NO flags/options or --help will display the syntax requirement.


=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BSML::BsmlBuilder;
use File::Basename;
use Pod::Usage;
use Log::Log4perl qw(get_logger);


my %options = ();
my $results = GetOptions (\%options, 
			  'ali_dir|a=s', 
			  'output|o=s',
			  'verbose|v', 
			  'program|p=s',
			  'suffix|s=s',
			  'help|h',
			  'configfile|c=s',
			  'log4perl|l=s',
			  'clusterDir|k=s',
			  'man') || pod2usage();

###-------------PROCESSING COMMAND LINE OPTIONS-------------###




&print_usage() if ((exists $options{'help'}) and (defined($options{'help'}))); 


my $msf_suffix = $options{'suffix'} if ((exists $options{'suffix'}) and (defined($options{'suffix'}))); 
$msf_suffix = "msf" if (!defined($msf_suffix));

my $program = $options{'program'} if ((exists $options{'program'}) and (defined($options{'program'}))); 
$program = "clustal" if (!defined($program));

my $configfile = $options{'configfile'} if ((exists $options{'configfile'}) and (defined($options{'configfile'}))); 
print STDERR "Configuration file was not defined\n" if (!defined($configfile));

my $output = $options{'output'} if ((exists $options{'output'}) and (defined($options{'output'}))); 
print STDERR "output was not defined\n" if (!defined($output));

my $ali_dir    = $options{'ali_dir'} if ((exists $options{'ali_dir'}) and (defined($options{'ali_dir'})));
print STDERR "ali_dir was not defined\n" if (!defined($ali_dir));
$ali_dir =~ s/\/+$// if($ali_dir);

my $log4perl = $options{'log4perl'} if ((exists $options{'log4perl'}) and (defined($options{'log4perl'})));
if (!defined($log4perl)){
    $log4perl = "msf2bsml.log";
}

my $clusterDir = $options{'clusterDir'} if ((exists $options{'clusterDir'}) and (defined($options{'clusterDir'})));


&cmd_check();

###-------------------------------------------------------###

#
# Initialize log4perl
#
Log::Log4perl->init(
		    \ qq{
			log4perl.logger                       = WARN, A1, Screen
			log4perl.appender.A1                  = Log::Dispatch::File
			log4perl.appender.A1.filename         = $log4perl
			log4perl.appender.A1.mode             = write
			log4perl.appender.A1.Threshold        = WARN
			log4perl.appender.A1.layout           = Log::Log4perl::Layout::PatternLayout
			log4perl.appender.A1.layout.ConversionPattern = %d %p> %F{1}:%L %M - %m%n 
			log4perl.appender.Screen              = Log::Dispatch::Screen
			log4perl.appender.Screen.layout       = Log::Log4perl::Layout::SimpleLayout
                        #log4perl.appender.Screen.layout.ConversionPattern =%d %p> %F{1}:%L %M - %m%n 
			log4perl.appender.Screen.Threshold    = WARN
			Log::Log4perl::SimpleLayout
		    }
		    );

my $logger = get_logger("msf2bsml");


#
# Create and retrieve the key-value pairs configuration hash
#
my $config = &retrieve_config_hash($configfile) if (defined($configfile));


#
# Instantiate the BSML builder object
#
$logger->info("Instantiating the BSML builder object");
my $builder = new BSML::BsmlBuilder;
$logger->logdie("builder was not defined") if (!defined($builder));


#
# Create <Analysis> component of the BSML document
#
&create_analysis_component($builder, $output, $program, $config);

#
# If clusterDir was defined/specified by user, this script checks whether
# a file "no_cluster.txt" is present in clusterDir.  
#
if (defined($clusterDir)){
    if (&no_clusters(\$clusterDir)){
	#
	# "no_cluster.txt" was found and this script will therefore
	# produce an empty BSML document
	#
	$builder->write( $output );
	chmod 0777, $output;
	exit(0);
    }
}

#
# Construct the main portion of the BSML document
#
&build_bsml_document_core($builder, $ali_dir, $output);





#-------------------------------------------------------------------------------------------------------------
#
#           END OF MAIN --- SUBROUTINES FOLLOW
#
#-------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------
# no_clusters()
#
#-------------------------------------------------------------------
sub no_clusters {

    $logger->debug("Entered no_clusters") if $logger->is_debug;
    my $dir = shift;


    opendir(DIR, $$dir) or logger->logdie("Unable to access $$dir due to $!");
    while( my $file = readdir(DIR)) {
	next if ($file =~ /^\.{1,2}$/);  #skip  "." ,  ".."

	if ($file =~ /^no_clusters\.txt$/o) {
	    $logger->warn("\"no_clusters.txt\" flag file was detected in $$dir, this means there were no clusters generated upstream of this step.  This will yield an empty BSML multiple alignment document");
	    return 1;
	}
    }
    return 0;

}#end sub any_clusters()


#------------------------------------------------------------
# create_analysis_component()
#
#------------------------------------------------------------
sub create_analysis_component {

    $logger->debug("Entered create_analysis_component") if $logger->is_debug;

    my ($builder, $filename, $program, $configfile) = @_;

    $logger->logdie("builder was not defined")  if (!defined($builder));
    $logger->logdie("filename was not defined") if (!defined($filename));
    $logger->logdie("program was not defined")  if (!defined($program));
#    $logger->logdie("config was not defined")   if (!defined($config));

    my $analysis = $builder->createAndAddAnalysis( 
						   'sourcename'         => $filename,
						   'programversion'     => '1.0',
						   'program'            => $program,
						   'bsml_link_url'      => 'BsmlTables',
						   'bsml_link_relation' => 'MULTIPLE_ALIGNMENTS',
						   );

    foreach my $key (sort keys %$config){
	$logger->logdie("key was not defined") if (!defined($key));
	my $value = $config->{$key};
	$logger->logdie("key was not defined") if (!defined($key));

	# This appends custom Bsml Attributes to it
	$analysis->addBsmlAttr( $key, $value ); 	
    }


}#end sub create_analysis_component()

#---------------------------------------------------------------------
# build_bsml_document_core()
#
#---------------------------------------------------------------------
sub build_bsml_document_core {

    $logger->debug("Entered build_bsml_document_core") if $logger->is_debug;

    my ($builder, $ali_dir, $output) = @_;
    $logger->logdie("builder was not defined")  if (!defined($builder));
    $logger->logdie("ali_dir was not defined")  if (!defined($ali_dir));
    $logger->logdie("output was not defined")   if (!defined($output));


    my $file_found = 0;
    opendir(DIR, $ali_dir) or die "Unable to access $ali_dir due to $!";
    while( my $file = readdir(DIR)) {
	next if ($file =~ /^\.{1,2}$/);  #skip  "." ,  ".."
	if($file =~ /^(.+)\.$msf_suffix$/o) {
	    $file_found++;
	    my $fam_name = $1;
	    my $MSF_alignments = process_MSF_file("$ali_dir/$file");
	    next if(keys %$MSF_alignments < 1);   #skip empty msf files

	    my $table = $builder->createAndAddMultipleAlignmentTable('molecule-type' => $MSF_alignments->{'mol_type'});
	    $logger->logdie("table was not defined") if (!defined($table));

	    my $summary = $builder->createAndAddAlignmentSummary( 
								  'multipleAlignmentTable' => $table,
								  'seq-type'               => $MSF_alignments->{'mol_type'},
								  'seq-format'             => 'msf' 
								  );
	    $logger->logdie("summary was not defined") if (!defined($summary));
	    
	    
	    my $aln = $builder->createAndAddSequenceAlignment( 'multipleAlignmentTable' => $table );
	    $logger->logdie("aln was not defined") if (!defined($aln));

	    my $seqnum=0;
	    my $sequences_tag;

	    foreach my $seq (keys %{ $MSF_alignments->{'proteins'} }) {
		$logger->logdie("seq was not defined") if (!defined($seq));

		$seqnum++;

		my $alignment = join ('', @{ $MSF_alignments->{'proteins'}->{$seq}->{'alignment'} });
		$logger->logdie("alignment was not defined") if (!defined($alignment));


		my $align_length = $MSF_alignments->{'proteins'}->{$seq}->{'length'} if ((exists $MSF_alignments->{'proteins'}->{$seq}->{'length'}) and (defined($MSF_alignments->{'proteins'}->{$seq}->{'length'})));

		$logger->logdie("align_length was not defined") if (!defined($align_length));

		#IMPORTANT!!!!
		#In order to ensure that each seq in a multiple sequence alignment is truly
		#unique, the seq-name and name will be in the form "protein_accession:seqnum"
		#i.e. (ana1.10005.m00234_protein:1). 
		
		$builder->createAndAddAlignedSequence(
						      'alignmentSummary' => $summary,
						      'seqnum'           => $seqnum,
						      'length'           => $align_length,
						      'name'             => "$seq:$seqnum"
						      );

		$builder->createAndAddSequenceData(
						   'sequenceAlignment' => $aln,
						   'seq-name'          => "$seq:$seqnum",
						   'seq-data'          => $alignment
						   ); 
		$sequences_tag .= "$seqnum:";
	    }
	    $aln->addattr( 'sequences', $sequences_tag );
	    
	}
	
    }
    
    if($file_found) {
	$builder->write( $output );
	chmod 0777, $output;
    } else {
	print STDERR "No files that ends in \"$msf_suffix\" are found in $ali_dir\n";
    }


}#end sub build_bsml_document_core()


#-------------------------------------------------------------
# process_MSF_file()
#
#-------------------------------------------------------------
sub process_MSF_file {

    $logger->debug("Entered process_MSF_file") if $logger->is_debug;

    my $file = shift;
    $logger->logdie("file was not defined") if (!defined($file));

    my $MSF_alignments ={};
    open(MSF, "$file") or die "Unable to open $file due to $!";
    my $line;
    my $msf_type;
    while(defined($line = <MSF>) and $line !~ /^\/\//) {
	if( $line =~ /MSF:\s*([\S]+)\s*Type:\s*([\S]+)\s*Check/) {
	    my $msf_length = $1;
	    return undef if($msf_length == 0);   #abort if align_len = 0

	    if($2 eq 'P') {
		$msf_type = 'protein';
	    }elsif($2 eq 'N') {
		$msf_type = 'nucleotide';
	    }else {
		$msf_type = 'protein';
	    }
	    $MSF_alignments->{'mol_type'} = $msf_type;
	}
	#if($line =~ /Name:\s+([\S]+)\s+[o]{2}\s+Len:\s+([\S]+)\s+Check:\s+([\S]+)\s+Weight:\s+([\S]+)/) {

	if($line =~ /Name:\s*([\S]+)\s*[o]{2}\s*Len:\s*([\S]+)\s*Check:\s*([\S]+)\s*Weight:\s*([\S]+)/) {
	    my $name    = $1;
	    my $ali_len = $2;
	    my $check   = $3;
	    my $weight  = $4;
	    
	    $MSF_alignments->{'proteins'}->{$name}->{'length'} = $ali_len;
	    $MSF_alignments->{'proteins'}->{$name}->{'check'}  = $check;
	    $MSF_alignments->{'proteins'}->{$name}->{'weight'} = $weight;
	    $MSF_alignments->{'proteins'}->{$name}->{'alignment'} = [];
	}
    }

    my $replacements;
    my $spaces;
    while($line = <MSF>) {
	if($line =~ /^([\S]+)/) {
	    my $name = $1;
	    if(exists($MSF_alignments->{'proteins'}->{$name})) {
		push( @{ $MSF_alignments->{'proteins'}->{$name}->{'alignment'} }, $line );
            } else {
		print STDERR "ERROR, $name is not valid protein name for $file\n";
		exit;
            }
	}
    }

    return $MSF_alignments;

}#end sub process_MSF_file()


#------------------------------------------------------------
# retrieve_config_hash()
#
#------------------------------------------------------------
sub retrieve_config_hash {

    $logger->debug("Entered retrieve_config_hash") if $logger->is_debug;

    my $file = shift;

    my $contents = &get_file_contents(\$file);
    $logger->logdie("contents was not defined") if (!defined($contents));

    my $hash = {};

    foreach my $line (@$contents){

	my ($key, $value) = split(/=/, $line);
	$logger->logdie("key was not defined")   if (!defined($key));
	$logger->logdie("value was not defined") if (!defined($value));
	
	$hash->{$key} = $value;
    }

    return $hash;

}#end sub retrieve_config_hash()

#-------------------------------------------------------------------
# get_file_contents()
#
#-------------------------------------------------------------------
sub get_file_contents {

    $logger->debug("Entered get_file_contents") if $logger->is_debug;

    my $file = shift;

    $logger->logdie("file was not defined") if (!defined($file));

    if (&is_file_readable($file)){

	open (IN_FILE, "<$$file") || $logger->logdie("Could not open file: $$file for input");
	my @contents = <IN_FILE>;
	chomp @contents;
	
	return \@contents;

    }
    else{
	$logger->logdie("file $$file does not have appropriate permissions");
    }
    
}#end sub get_contents()


#-------------------------------------------------------------------
# is_file_readable()
#
#-------------------------------------------------------------------
sub is_file_readable {

    my $file = shift;
    $logger->logdie("file was not defined") if (!defined($file));
      
    my $fatal_flag=0;

    if (!-e $$file){
	$logger->fatal("$$file does not exist");
	$fatal_flag++;
    }
    if ((-e $$file) and (!-r $$file)){
	$logger->fatal("$$file does not have read permissions");
	$fatal_flag++;
    }

    if ($fatal_flag>0){
	return 0;
    }
    else{
	return 1;
    }

}#end sub is_file_readable()


#-------------------------------------------------------------
# cmd_check()
#
#-------------------------------------------------------------
sub cmd_check {
#quality check

    if( exists($options{'man'})) {
	pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT});
    }   

    if( exists($options{'help'})) {
	pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT});
    }

    if(!$output or !$ali_dir) {
	pod2usage({-exitval => 2,  -message => "$0: All the required options are not specified", -verbose => 1, -output => \*STDERR});
    }
    
    if(! -d $ali_dir) {
	print STDERR "$ali_dir directory NOT found.  Aborting...\n";
	exit 5;
    }

}#end sub cmd_check()	    




#--------------------------------------------------------------
# print_usage()
#
#
#--------------------------------------------------------------
sub print_usage {

    print STDERR "USAGE:  $0 -a ali_dir -o output [-c configfile] [-l log4perl] [-p program] [-s suffix] [-h help] [-k clusterDir]\n\n";
    print STDERR " -a|--ali_dir      = alignment directory containing input alignment files\n";
    print STDERR " -o|--output       = BSML output filename\n";
    print STDERR " -c|--configfile   = input configuration file (for workflow related key-value pair data to be stored in the <Analysis> component\n";
    print STDERR " -l|--log4perl     = log4perl log file\n";
    print STDERR " -p|--program      = multiple alignment program (to be stored as the BSML attribute \"program\" in the <Analysis> component -default is \"clustal\"\n";
    print STDERR " -s|--suffix       = alignment filename extension (default is \".msf\")\n";
    print STDERR " -h|--help         = display this print_usage message\n";
    print STDERR " -k|--clusterDir   = cluster directory\n\n\n";

    exit(0);


}#end sub print_usage()
