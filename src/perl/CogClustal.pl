#! /local/perl/bin/perl
#---------------------------------------------------------------------------------
#
#
#
# purpose:           Runs clustalw
#
# modification 2003-11-24 Jay Sundaram
#              1) Added print usage details (improved interface)
#              
#              2003-12-09 Jay Sundaram sundaram@tigr.org
#              1) If detects "no_clusters.txt" in cluster directory,
#                 no clustalw alignment will be produced
#              2) Introduced log4perl logging
#           
#              2003-12-24 sundaram
#              1) ls was replaced by find
#    
#
#              2003-12-30 sundaram
#              1) using /usr/local/devel/ANNOTATION/cas/tester/clustalw
#              2) clustalw no longer re-directed to /dev/null
#
#
#
#
#---------------------------------------------------------------------------------
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Log::Log4perl qw(get_logger);


my %options = ();
my $results = GetOptions( 
			  \%options,
			  'fastaDir|f=s',
			  'outputDir|o=s',
			  'help|h',
			  'log4perl|l=s',
			  'clusterDir|c=s',
			  );

my $fastaDir      = $options{'fastaDir'}    if ((exists $options{'fastaDir'}) and (defined($options{'fastaDir'})));
my $outputDir     = $options{'outputDir'}   if ((exists $options{'outputDir'}) and (defined($options{'outputDir'})));
my $help          = $options{'help'}        if ((exists $options{'help'}) and (defined($options{'help'})));
my $log4perl      = $options{'log4perl'}    if ((exists $options{'log4perl'}) and (defined($options{'log4perl'})));
my $clusterDir    = $options{'clusterDir'}  if ((exists $options{'clusterDir'}) and (defined($options{'clusterDir'})));


#
# If help, display the usage info...
#
&print_usage() if ($help);

print STDERR "\nfastaDir was not defined\n\n"  if (!defined($fastaDir));

#
# fasta directory must be specified by user
#
if (!$fastaDir){
    &print_usage();
}

$log4perl = "./cogclustal.log" if (!defined($log4perl));

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

my $logger = get_logger("cogclustal");


#
# Determine whether there were any cluster files produced.
#
&any_clusters(\$clusterDir) if (defined($clusterDir));


#
# If user does not specify the output directory, set to current directory
#
if (!$outputDir){
    $outputDir = ".";
    $logger->info("outputDir not defined, setting to $outputDir");
}

#
# Strip the trailing "/" and verify whether a) the directory exists b)the directory is in fact a directory
#
$outputDir =~ s/\/+$//;
$logger->logdie("$outputDir does not exist") if (!-e $outputDir);
$logger->logdie("$outputDir is not a valid directory") if (!-d $outputDir);


$fastaDir =~ s/\/+$//; 
$logger->logdie("$fastaDir does not exist") if (!-e $fastaDir);
$logger->logdie("$fastaDir is not a valid directory") if (!-d $fastaDir);



my $fastafiles = &retrieve_fasta_files(\$fastaDir);
$logger->logdie("fastafiles was not defined") if (!defined($fastafiles));



#
# For each fasta file in the fasta directory, run clustalw
#

foreach my $fastaFile (@$fastafiles){
    my $clustalFile = $fastaFile;
    $clustalFile =~ s/$fastaDir/$outputDir/;
    $clustalFile =~ s/fasta/clustal/;

    # clustalw spits out a lot of stuff to STDOUT this redirects it to /dev/null
    #my $status = system( "clustalw -output=gcg -infile=$fastaFile -outfile=$clustalFile > /dev/null" );
    my $status = system( "clustalw -output=gcg -infile=$fastaFile -outfile=$clustalFile" );

    my $exit_value = $status >> 8;
    my $signal_num = $status & 127;
    my $dumped_core = $status & 128;

    if( !($exit_value == 0) )
    {
	exit( $exit_value );
    }

    $logger->info("Deleting the temporary .dnd files from $fastaDir\n");
    unlink "$fastaDir/*.dnd";
}

#----------------------------------------------------
# retrieve_fasta_files()
#
#----------------------------------------------------
sub retrieve_fasta_files {

    my $dir = shift;


   my @bsml_files = qx[find $$dir -name "*.fasta" -type f];

    chomp @bsml_files;


    $logger->info("No .fasta files were found in directory $$dir") if (scalar(@bsml_files) < 1);


    return \@bsml_files;



}#end sub retrieve_fasta_files()

#-------------------------------------------------------------------
# any_clusters()
#
#-------------------------------------------------------------------
sub any_clusters {

    $logger->debug("Entered any_clusters") if $logger->is_debug;
    my $dir = shift;

    opendir(DIR, $$dir) or logger->logdie("Unable to access $$dir due to $!");
    while( my $file = readdir(DIR)) {
	next if ($file =~ /^\.{1,2}$/);  #skip  "." ,  ".."

	if ($file =~ /^no_clusters\.txt$/o) {
	    $logger->warn("\"no_clusters.txt\" flag file was detected in $$dir, thus no multifasta files will be created");
	    exit(0);
	}
    }


}#end sub any_clusters()


#----------------------------------------------------
# print_usage()
#
#----------------------------------------------------
sub print_usage {
    
    print STDERR "\nUsage: $0 -f fastaDir -o outputDir [-h|--help]\n";
    print STDERR " -f|--fastaDir    = valid directory containing the multi-fasta input files\n";
    print STDERR " -o|--outputDir   = valid directory to write the clustalw output files to\n";
    print STDERR " -h|--help        = displays this message\n";
    print STDERR "\n";
    exit(0);


}
