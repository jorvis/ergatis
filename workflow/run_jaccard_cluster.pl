#!/usr/local/bin/perl
#---------------------------------------------------------------------
# script name: run_jaccard_cluster.pl
# date:        2003.12.01
# author:      Jay Sundaram 301-610-5983 sundaram@tigr.org
#
#
#
#---------------------------------------------------------------------

=head1 NAME

run_jaccard_cluster.pl - Creates and executes the bsml to chado migration workflow

=head1 SYNOPSIS

USAGE:  run_jaccard_cluster.pl -U username -P password -D database [-g asm_group_name] [-a asmbl_list|ALL] [-F asmbl_file] [-b compute_type] [-S server] [-i pidentity] [-r pvalue] [-l linkscore_list] [-c config] [-h] [-m] [-p] [-v] [-w]

=head1 OPTIONS

=over 8

=item B<--username,-U>
    
    Database username

=item B<--password,-P>
    
    Database password

=item B<--database,-D>
    
    Source database name

=item B<--server,-S>
    
    target server

=item B<--asmbl_file,-F>
    
    File containing list of assembly identifiers

=item B<--asm_group_name,-g>

    group name for the list of assembly identifiers (to be incorporated into the BSML document name)

=item B<--asmbl_list,-a>

    Comma separated list of assembly identifiers or ALL

=item B<--compute_type,-b>

    Source pairwise alignment compute type either blastp OR ber

=item B<--pidentity,-i>

    percent identity threshold. default 75

=item B<--pvalue,-r>

    p_value threshold. default 1e-15

=item B<--linkscore_list,-l>

    List of Jaccard coefficient link score cut-offs

=item B<--config,-c>

    Configuration file containing all substitution keys parameters

=item B<--help,-h>

    Prints out usage

=item B<--printconf,-p>

    Prints out the default configuration key-value pairs

=item B<--man,-m>

    Displays this pod2usage man pages for this script

=item B<--verbose,-v>

    Increases log4perl reporting level to screen from WARN to INFO

=item B<--workflow_monitor,-w>

    Automatically launch the workflow monitor


=back

=head1 DESCRIPTION

    run_jaccard_cluster.pl - Creates and executes the bsml to chado migration workflow

=cut




use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Data::Dumper;
use Log::Log4perl qw(get_logger);
use File::Basename;
use Pod::Usage;

my $verstr = qw($Name$);
my ($verstrparse) = ($verstr =~ /Name:\s+(.*)\s+/);
use constant VERSION_STRING => $verstrparse;
 


my $WorkflowConfigFileDir = $ENV{'WORKFLOW_DOCS_DIR'};
my $WorkflowExecDir       = $ENV{'WORKFLOW_WRAPPERS_DIR'};

#
# The default bsml2chado configuration
#
my $conf = {
    ';ASM_GROUP_NAME;'               => '',
    ';ASMBL_LIST;'                   => '',
    ';ASMBL_FILE;'                   => '',
    ';JACCARD_CONF;'                 => $WorkflowConfigFileDir ."/jaccard_cluster.ini", 
    ';JACCARD_MASTER_CONF;'          => $WorkflowConfigFileDir ."/jaccard_cluster-master.ini", 
    ';JACCARD_MASTER_TEMPLATE;'      => $WorkflowConfigFileDir ."/jaccard_cluster-master_template.xml",
    ';JACCARD_TEMPLATE;'             => $WorkflowConfigFileDir ."/jaccard_cluster_template.xml",
    ';COMPUTE_TYPE;'                 => 'blastp',
    ';CONFIG;'                       => '',
    ';INSTALLATION;'                 => $WorkflowExecDir,
    ';LINKSCORE_LIST;'               => '0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9',
    ';PASSWORD;'                     => '',               
    ';PIDENTITY;'                    => '75',
    ';REPOSITORY_ROOT;'              => '/usr/local/annotation',
    ';RUN_WF;'                       => $WorkflowExecDir ."/run_wf.sh",
    ';SET_RUNTIME;'                  => $WorkflowExecDir ."/set_runtime_vars.pl",
    ';PVALUE;'                       => '1e-15',
    ';SIMILARITY_TYPE;'              => 'jaccard',
    ';DATABASE;'                     => '',               
    ';SERVER;'                       => 'SYBTIGR',
    ';TMP_DIR;'                      => '/usr/local/scratch/jaccard',
    ';USERNAME;'                     => '',
    ';WFNAME;'                       => 'jaccard',	
    ';WORKFLOW_DIR;'                 => $WorkflowExecDir,
    ';VERSION_STRING;'               => VERSION_STRING,
    
};


my %options = ();
my ($username, $password, $database, $server, $pidentity, $pvalue, $linkscore_list, $config, $help, $man, $printconf, $compute_type, $workflow_monitor, $verbose, $asmbl_list, $asmbl_file, $asm_group_name);
my $results = GetOptions (
			  \%options,
			  'username|U=s'            => \$username,
			  'password|P=s'            => \$password,
			  'database|D=s'            => \$database,
			  'compute_type|b=s'        => \$compute_type,
			  'server|S=s'              => \$server,
			  'pidentity|i=s'           => \$pidentity,
			  'pvalue|r=s'              => \$pvalue,
			  'linkscore_list|l=s'      => \$linkscore_list,
			  'config|c=s'              => \$config,
			  'help|h'                  => \$help,
			  'man|m'                   => \$man,
			  'printconf|p'             => \$printconf,
			  'verbose|v'               => \$verbose,
			  'workflow_monitor|w'      => \$workflow_monitor,
			  'asmbl_list|a=s'          => \$asmbl_list,
			  'asmbl_file|F=s'          => \$asmbl_file,
			  'asm_group_name|g=s'      => \$asm_group_name,
			  );

&print_configuration($conf) if ($printconf);
&pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if ($man);
&print_usage() if ($help);


my $log4perl = "jaccard_wf.log";

#
# Increase log4perl reporting level to screen
#
my $screen_threshold = 'ERROR';
$screen_threshold = 'INFO' if ($verbose);

#
# Initialize the log4perl logger
#
Log::Log4perl->init(
		    \ qq{
			log4perl.logger                       = INFO, A1, Screen
			log4perl.appender.A1                  = Log::Dispatch::File
			log4perl.appender.A1.filename         = $log4perl
			log4perl.appender.A1.mode             = write
			log4perl.appender.A1.Threshold        = INFO
			log4perl.appender.A1.layout           = Log::Log4perl::Layout::PatternLayout
			log4perl.appender.A1.layout.ConversionPattern = %d %p> %F{1}:%L %M - %m%n 
			log4perl.appender.Screen              = Log::Dispatch::Screen
			log4perl.appender.Screen.layout       = Log::Log4perl::Layout::SimpleLayout
                        #log4perl.appender.Screen.layout.ConversionPattern =%d %p> %F{1}:%L %M - %m%n 
			log4perl.appender.Screen.Threshold    = $screen_threshold
			Log::Log4perl::SimpleLayout
		    }
		    );

#
# Instantiate logger object
#
my $logger = get_logger();


#
# Load the conf hash with override values defaults if the user has specified a config file
# 
&override_default_configuration(\$config, $conf) if ($config);

#
# So the order of precedence is:
# if the user doesn't specify the particular parameter on the command-line,
# then check for the parameter in the config file.
# if the user didn't specify the particular parameter in the config file,
# then use the default configuration
# if the default configuration does not contain the particular
# parameter, fail.
#
if (!defined($username)){
    $username = $conf->{';USERNAME;'} if ((exists $conf->{';USERNAME;'}) and (defined($conf->{';USERNAME;'})));
}
if (!defined($password)){
    $password = $conf->{';PASSWORD;'} if ((exists $conf->{';PASSWORD;'}) and (defined($conf->{';PASSWORD;'})));
}
if (!defined($database)){
    $database = $conf->{';DATABASE;'} if ((exists $conf->{';DATABASE;'}) and (defined($conf->{';DATABASE;'})));
}
if (!defined($server)){
    $server = $conf->{';SERVER;'} if ((exists $conf->{';SERVER;'}) and (defined($conf->{';SERVER;'})));
}
if (!defined($compute_type)){
    $compute_type = $conf->{';COMPUTE_TYPE;'} if ((exists $conf->{';COMPUTE_TYPE;'}) and (defined($conf->{';COMPUTE_TYPE;'})));
}
if (!defined($pidentity)){
    $pidentity = $conf->{';PIDENTITY;'} if ((exists $conf->{';PIDENTITY;'}) and (defined($conf->{';PIDENTITY;'})));
}
if (!defined($pvalue)){
    $pvalue = $conf->{';PVALUE;'} if ((exists $conf->{';PVALUE;'}) and (defined($conf->{';PVALUE;'})));
}
if (!defined($linkscore_list)){
    $linkscore_list = $conf->{';LINKSCORE_LIST;'} if ((exists $conf->{';LINKSCORE_LIST;'}) and (defined($conf->{';LINKSCORE_LIST;'})));
}
if (!defined($asmbl_list)){
    $asmbl_list = $conf->{';ASMBL_LIST;'} if ((exists $conf->{';ASMBL_LIST;'}) and (defined($conf->{';ASMBL_LIST;'})));
}
if (!defined($asmbl_file)){
    $asmbl_file = $conf->{';ASMBL_FILE;'} if ((exists $conf->{';ASMBL_FILE;'}) and (defined($conf->{';ASMBL_FILE;'})));
}
if (!defined($asm_group_name)){
    $asm_group_name = $conf->{';ASM_GROUP_NAME;'} if ((exists $conf->{';ASM_GROUP_NAME;'}) and (defined($conf->{';ASM_GROUP_NAME;'})));
}


#
# These parameters may be set by the user in the configuration file, but not via the command-line interface
#
my $workflow_dir = $conf->{';WORKFLOW_DIR;'} if ((exists $conf->{';WORKFLOW_DIR;'}) and (defined($conf->{';WORKFLOW_DIR;'})));
$logger->logdie("workflow_dir was not defined") if (!defined($workflow_dir));

my $set_runtime = $conf->{';SET_RUNTIME;'} if ((exists $conf->{';SET_RUNTIME;'}) and (defined($conf->{';SET_RUNTIME;'})));
$logger->logdie("set_runtime was not defined") if (!defined($set_runtime));

my $run_wf = $conf->{';RUN_WF;'} if ((exists $conf->{';RUN_WF;'}) and (defined($conf->{';RUN_WF;'})));
$logger->logdie("run_wf was not defined") if (!defined($run_wf));


#
# Now verify whether the correct combination of parameters have been provided by user at command-line OR in configuration file OR in default configuration
#
print STDERR "\nusername was not defined\n\n"          if (!defined($username));
print STDERR "\npassword was not defined\n\n"          if (!defined($password));
print STDERR "\ndatabase was not defined\n\n"          if (!defined($database));
print STDERR "\nserver was not defined\n\n"            if (!defined($server));
print STDERR "\ncompute_type was not defined\n\n"      if (!defined($compute_type));
print STDERR "\npidentity was not defined\n\n"         if (!defined($pidentity));
print STDERR "\npvalue was not defined\n\n"            if (!defined($pvalue));
print STDERR "\nlinkscore_list was not defined\n\n"    if (!defined($linkscore_list));

&print_usage() if (!$username or !$password or !$database or !$compute_type or !$server or !$pidentity, or !$pvalue or !$linkscore_list);

#
# Determine whether server was defined and whether is either SYBIL or SYBTIGR
#
if ($server !~ /^SYBIL|SYBTIGR$/){
    print STDERR "Invalid Sybase server: $server, must be either SYBIL or SYBTIGR";
    &print_usage();
}

$compute_type = lc($compute_type);
#
# Compute_type must be either 'blastp' or 'ber'
#
if ($compute_type !~/^blastp|ber$/){
    print STDERR "\ncompute_type must be either \"blastp\" or \"ber\"\n\n";
    &print_usage();
}


#
# If asmbl_list defined, but asm_group_name not, set default value for asm_group_name
#
if (defined($asmbl_list) and ((!defined($asm_group_name)) or ($asm_group_name =~ /^\s*$/))){
    $asm_group_name = "jaccard_" . $$;
    $logger->info("The asmbl_list you provided will be assigned the following asm_group_name:$asm_group_name since you did not supply one");
}

#
# Retieve linkscores
#
my $linkscores = &retrieve_linkscores($linkscore_list);
$logger->logdie("linkscore_list was not defined") if (!defined($linkscores));

#
# Retrieve list of assembly identifiers - note that only the unique set will be processed
#
my $asm_hash = &retrieve_asm_hash($asmbl_list, $asmbl_file, $asm_group_name);
$logger->logdie("asm_hash was not defined") if (!defined($asm_hash));


#
# Set the final configuration values
#
$conf->{';PAIRWISE_COMPUTE_TYPE_UC;'}    = uc($compute_type)  if (defined($compute_type));
$conf->{';PAIRWISE_COMPUTE_TYPE_LC;'}    = lc($compute_type)  if (defined($compute_type));
$conf->{';PIDENTITY;'}                   = $pidentity         if (defined($pidentity));
$conf->{';PVALUE;'}                      = $pvalue            if (defined($pvalue));
$conf->{';LINKSCORE_LIST;'}              = $linkscore_list    if (defined($linkscore_list));
$conf->{';USERNAME;'}                    = $username          if (defined($username));
$conf->{';PASSWORD;'}                    = $password          if (defined($password));
$conf->{';SERVER;'}                      = $server            if (defined($server));
$conf->{';DATABASE_LC;'}                 = lc($database)      if (defined($database));
$conf->{';DATABASE_UC;'}                 = uc($database)      if (defined($database));
$conf->{';WFID;'}                        = $$;
$conf->{';LINKSCORE_LIST;'}              = $linkscore_list;

#
# Verify and set the BSML_repository
#
$logger->logdie("repository root was not defined") if ((!exists $conf->{';REPOSITORY_ROOT;'}) or (!defined($conf->{';REPOSITORY_ROOT;'})));
my $bsml_repository = $conf->{';REPOSITORY_ROOT;'} . "/" . $conf->{';DATABASE_UC;'} . "/BSML_repository";
$logger->logdie("$bsml_repository is not a valid directory") if (!-d $bsml_repository);
$conf->{';BSML_REPOSITORY;'} = $bsml_repository;

#
# Verify and set the pairwise compute repository
#
$logger->logdie("pairwise compute_type was not defined") if ((!exists $conf->{';PAIRWISE_COMPUTE_TYPE_LC;'}) or (!defined($conf->{';PAIRWISE_COMPUTE_TYPE_LC;'})));
my $pairwise_repository = $bsml_repository . "/" . $conf->{';PAIRWISE_COMPUTE_TYPE_LC;'};
$logger->logdie("$pairwise_repository is not a valid directory") if (!-d $pairwise_repository);
$conf->{';PAIRWISE_COMPUTE_BSML_DIRECTORY;'} = $pairwise_repository;

#
# The output of the MSF2Bsml.pl script is a BSML search encoding document.  This should be stored in the BSML_repository/jaccard repository
#
my $jaccard_repository = $bsml_repository . "/jaccard";
$conf->{';JACCARD_REPOSITORY;'} = $jaccard_repository;

#
# Verify that the temp directory is valid
#
my $tmpdir = $conf->{';TMP_DIR;'} if ((exists $conf->{';TMP_DIR;'}) and (defined($conf->{';TMP_DIR;'})));
$logger->logdie("tmpdir was not defined") if (!defined($tmpdir));
if (!-d $tmpdir){
    #
    # /usr/local/scratch/jaccard does not exist, therefore create now.
    #
    my $execution_string = "mkdir -p -m 777 $tmpdir";
    &execute(\$execution_string);
}

my $database_tmpdir = $tmpdir . "/" . $conf->{';DATABASE_LC;'};
# e.g. /usr/local/scratch/jaccard/chado_tryp
$conf->{';DATABASE_TMPDIR;'} = $database_tmpdir;


#
# Create the Workflow instance subfolder
#
my $workflow_instance_dir = &create_workflow_instance_subfolder();
$logger->logdie("workflow_instance_dir was not defined") if (!defined($workflow_instance_dir));

#
# Generate and prepare the master subs file.  This will store all of the configuration values and then perform sed replacement on all configuration values
#
&prepare_master_subsfile($linkscores, $workflow_instance_dir, $set_runtime, $conf, $database_tmpdir, $asm_hash);

#
# Generate and prepare the subflows
#
&prepare_subflows($linkscores, $workflow_instance_dir, $set_runtime, $conf, $database_tmpdir, $asm_hash);

#
# Run the jaccard workflow
#
&run_workflow($run_wf, $workflow_instance_dir);


print STDERR "\n$0 execution complete\nPlease review log4perl log file:$log4perl\n\n";




#----------------------------------------------------------------------------------------------------------------------------
#
#                          END OF MAIN -- SUBROUTINES FOLLOW
#
#----------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
# prepare_subflows()
#
# For each linkscore, generate a subflow
#
#------------------------------------------------------------------------------------------
sub prepare_subflows {


    $logger->debug("Entered prepare_subflows") if $logger->is_debug;

    my ($linkscores, $wid, $set_runtime, $conf, $database_tmpdir, $asmb_hash) = @_;
    $logger->logdie("linkscores was not defined")       if (!defined($linkscores));
    $logger->logdie("wid was not defined")              if (!defined($wid));
    $logger->logdie("set_runtime was not defined")      if (!defined($set_runtime));
    $logger->logdie("conf was not defined")             if (!defined($conf));
    $logger->logdie("database_tmpdir was not defined")  if (!defined($database_tmpdir));
    $logger->logdie("asmb_hash was not defined")        if (!defined($asmb_hash));

    #
    # General steps in each subflow pipeline are:
    # 1) Create clusters based on linkscore, pvalue and pidentity thresholds
    # 2) Generate multi-fasta files
    # 3) Produce clustalw alignments
    # 4) convert the MSF alignment document into a BSML document
    #
    #

    my $execution_string;

    
    foreach my $score (sort @$linkscores){
	$logger->logdie("score was not defined") if (!defined($score));
	
	#
	# linkscores should be values between 0 and 1
	#
	if (($score < 0) or ($score >= 1)){
	    $logger->logdie("score $score must be 0 < score < 1");
	    next;
	}
	
	#
	# For each assembly group we will process all of the assembly identifiers 
	# associated to this assembly group
	#
	foreach my $group (sort keys %$asmb_hash){
	    $logger->logdie("group was not defined") if (!defined($group));
	    
	    my $assembly_string_list;
	    #
	    # Grab each assembly from the group hash and build a string list
	    #
	    foreach my $asm (sort keys %{$asmb_hash->{$group}}){
		$logger->logdie("asm was not defined for group:$group") if (!defined($asm));
		
		$assembly_string_list .= $asm . ",";
	    }
    
	    #
	    # Strip the trailing comma from assembly_string_list
	    #
	    chop $assembly_string_list;

	    my $jaccard = $group . "_" . $score;
	    

	    #
	    # Store jaccard group name and list of assembly identifiers associated to this group
	    #
	    $conf->{';JACCARD_GROUP;'}               = $jaccard;
	    $conf->{';ASSEMBLY_IDENTIFIER_LIST;'}    = $assembly_string_list;
	    
	    #
	    # Create a subdirectory in the workflow instance directory for this jaccard
	    # We will store all the associated goodies in this bsml specific directory
	    #
	    $execution_string = "mkdir -p -m 777 ". $wid;
	    &execute(\$execution_string);
	    $execution_string = "mkdir -p -m 777 ". $wid  . "/" . $jaccard;
	    &execute(\$execution_string);
	    
	    #
	    # foreach linkscore, we will create a set of scratch subfolders
	    #
	    my $workflow_scratch_space = $database_tmpdir ."/". $conf->{';WFID;'} ."/". $jaccard;
	    #
	    # e.g. /usr/local/scratch/jaccard/chado_tryp/6144/jaccard_0.1
	    #
	    $conf->{';WORKFLOW_SCRATCH_SPACE;'} = $workflow_scratch_space;
	    
	    my $cluster_output_dir    = $workflow_scratch_space . "/" . "Clustering_dir";
	    #
	    # e.g. /usr/local/scratch/jaccard/chado_tryp/6144/jaccard_0.1/Clustering_dir
	    #
	    my $multifasta_output_dir = $workflow_scratch_space . "/" . "Fasta_dir";
	    #
	    # e.g. /usr/local/scratch/jaccard/chado_tryp/6144/jaccard_0.1/Fasta_dir
	    #
	    my $clustalw_output_dir   = $workflow_scratch_space . "/" . "Clustalw_dir";
	    #
	    # e.g. /usr/local/scratch/jaccard/chado_tryp/6144/jaccard_0.1/Clustalw_dir
	    #
	    my $msf_output            = $conf->{';JACCARD_REPOSITORY;'} . "/" . $jaccard . ".bsml";
	    
	    #
	    # temp scratch space output directories
	    #
	    $conf->{';CLUSTERING_OUTPUT_DIR;'}           = $cluster_output_dir;
	    $conf->{';FASTA_OUTPUT_DIR;'}                = $multifasta_output_dir;
	    $conf->{';CLUSTALW_OUTPUT_DIR;'}             = $clustalw_output_dir;
	    $conf->{';MSF2BSML_OUTPUT;'}                 = $msf_output;
	    
	    #
	    # All of the *.log and *.stats files necessary for the migration
	    #
	    $conf->{';CLUSTERING_LOG4PERL;'}             = $wid ."/". $jaccard . "/clustering.log";
	    $conf->{';CLUSTERING_LINKSCORE;'}            = $score;
	    $conf->{';MULTIFASTA_LOG4PERL;'}             = $wid . "/" . $jaccard . "/multifasta.log";
	    $conf->{';CLUSTALW;'}                        = $wid ."/". $jaccard . "/clustalw.txt";
	    $conf->{';CLUSTALW_LOG4PERL;'}               = $wid ."/". $jaccard . "/clustalw.log";

	    $conf->{';MSF2BSML_PROGRAM_NAME;'}           = "Jaccard-clustering " . $score . " Clustalw";
	    $conf->{';MSF2BSML_LOG4PERL;'}           = $wid ."/". $jaccard . "/msf2bsml.log";
	    $conf->{';MSF2BSML_CONFIG_FILE;'}        = $wid ."/". $jaccard . "/msf2bsml_config.txt";
	    

	    #
	    # Create the intermediary file which will contain all of the workflow parameters
	    # that will be store in the <Analysis> component and then chado.analysisprop table
	    #
	    &write_replacement_vars_to_msf2bsml_configfile($conf);

	    #
	    # Store the configuration for this specific jaccard
	    # and set global permissions
	    #
	    my $jaccard_subflow = $wid ."/". $jaccard ."/". $jaccard .".sub";
	    &store_conf_to_file(\$jaccard_subflow, $conf);
	    
	    $execution_string = "chmod 666 $jaccard_subflow";
	    &execute(\$execution_string);
	    
	    #
	    # Perform sed replacement on the specific bsml subflow configuration .ini file 
	    #
	    $execution_string = $set_runtime ." -c ". $jaccard_subflow ." < ". $conf->{';JACCARD_CONF;'} ." > ". $wid ."/". $jaccard .".ini";
	    &execute(\$execution_string);
	
	    #
	    # Give the subflow instance configuration file global permissions
	    #
	    $execution_string = "chmod 666 " . $wid ."/". $jaccard .".ini";
	    &execute(\$execution_string);
	    
	}
    }

    #
    # Copy the jaccard_cluster_template.xml        to jaccard_cluster-instance_template.xml AND
    #          jaccard_cluster-master_template.xml to jaccard_cluster-master-instance_template.xml
    #
    # and then set global permissions
    #
    $execution_string = "cp ". $conf->{';JACCARD_TEMPLATE;'} ." ". $wid ."/jaccard_cluster-instance_template.xml";
    &execute(\$execution_string);
    $execution_string = "cp ". $conf->{';JACCARD_MASTER_TEMPLATE;'} ." ". $wid ."/jaccard_cluster-master-instance_template.xml";
    &execute(\$execution_string);
    
    $execution_string = "chmod 666 ". $wid . "/jaccard_cluster-instance_template.xml";
    &execute(\$execution_string);
    $execution_string = "chmod 666 ". $wid . "/jaccard_cluster-master-instance_template.xml";
    &execute(\$execution_string);


}#end sub prepare_subflows()




#-------------------------------------------------------------------------
# write_replacement_vars_to_msf2bsml_configfile()
#
#-------------------------------------------------------------------------
sub write_replacement_vars_to_msf2bsml_configfile {

    $logger->debug("Entered write_replacement_vars_to_msf2bsml_configfile") if $logger->is_debug;
    
    my ($conf, $all) = @_;

    $logger->logdie("conf was not defined") if (!defined($conf));


    my $msf2bsml_configfile = $conf->{';MSF2BSML_CONFIG_FILE;'} if ((exists $conf->{';MSF2BSML_CONFIG_FILE;'}) and (defined($conf->{';MSF2BSML_CONFIG_FILE;'})));
    $logger->logdie("msf2bsml_configfile was not defined") if (!defined($msf2bsml_configfile));


    open (OUTFILE, ">$msf2bsml_configfile") or $logger->logdie("Could not open $msf2bsml_configfile in write mode");

    if (defined($all)){
	foreach my $key (sort keys %$conf){
	    $logger->logdie("key was not defined") if (!defined($key));
	    
	    my $value = $conf->{$key};
	    $logger->logdie("value was not defined") if (!defined($value));


	    #
	    # Strip the enclosing "';;'"
	    #
	    $key =~ s/^\;//;
	    $key =~ s/\;$//;
	
	
	    print OUTFILE $key ."=". $value . "\n";
	    
	}
    }
    else{

	# Threshold values
	my $val;

	$val = $conf->{';PIDENTITY;'} if (exists $conf->{';PIDENTITY;'} and (defined($conf->{';PIDENTITY;'})));
	print OUTFILE "PIDENTITY" ."=". $val . "\n" if (defined($val));

	$val = $conf->{';PVALUE;'} if (exists $conf->{';PVALUE;'} and (defined($conf->{';PVALUE;'})));
	print OUTFILE "PVALUE" ."=". $val . "\n" if (defined($val));

	$val = $conf->{';CLUSTERING_LINKSCORE;'} if (exists $conf->{';CLUSTERING_LINKSCORE;'} and (defined($conf->{';CLUSTERING_LINKSCORE;'})));
	print OUTFILE "CLUSTERING_LINKSCORE" ."=". $val . "\n" if (defined($val));

	$val = $conf->{';ASSEMBLY_IDENTIFIER_LIST;'} if (exists $conf->{';ASSEMBLY_IDENTIFIER_LIST;'} and (defined($conf->{';ASSEMBLY_IDENTIFIER_LIST;'})));
	print OUTFILE "ASSEMBLY_IDENTIFIER_LIST" ."=". $val . "\n" if (defined($val));


	# Configuration values:

	$val = $conf->{';WORKFLOW_DIR;'} if (exists $conf->{';WORKFLOW_DIR;'} and (defined($conf->{';WORKFLOW_DIR;'})));
	print OUTFILE "WORKFLOW_DIR" ."=". $val . "\n" if (defined($val));

	$val = $conf->{';VERSION_STRING;'} if (exists $conf->{';VERSION_STRING;'} and (defined($conf->{';VERSION_STRING;'})));
	print OUTFILE "VERSION_STRING" ."=". $val . "\n" if (defined($val));


	$val = $conf->{';COMPUTE_TYPE;'} if (exists $conf->{';COMPUTE_TYPE;'} and (defined($conf->{';COMPUTE_TYPE;'})));
	print OUTFILE "COMPUTE_TYPE" ."=". $val . "\n" if (defined($val));

	$val = $conf->{';SIMLIARITY_TYPE;'} if (exists $conf->{';SIMLIARITY_TYPE;'} and (defined($conf->{';SIMLIARITY_TYPE;'})));
	print OUTFILE "SIMLIARITY_TYPE" ."=". $val . "\n" if (defined($val));

	$val = $conf->{';WFID;'} if (exists $conf->{';WFID;'} and (defined($conf->{';WFID;'})));
	print OUTFILE "WFID" ."=". $val . "\n" if (defined($val));


	# Directories containing the intermediaries:  (self explanatory)

	$val = $conf->{';CLUSTERING_OUTPUT_DIR;'} if (exists $conf->{';CLUSTERING_OUTPUT_DIR;'} and (defined($conf->{';CLUSTERING_OUTPUT_DIR;'})));
	print OUTFILE "CLUSTERING_OUTPUT_DIR" ."=". $val . "\n" if (defined($val));

	$val = $conf->{';FASTA_OUTPUT_DIR;'} if (exists $conf->{';FASTA_OUTPUT_DIR;'} and (defined($conf->{';FASTA_OUTPUT_DIR;'})));
	print OUTFILE "FASTA_OUTPUT_DIR" ."=". $val . "\n" if (defined($val));

	$val = $conf->{';CLUSTALW_OUTPUT_DIR;'} if (exists $conf->{';CLUSTALW_OUTPUT_DIR;'} and (defined($conf->{';CLUSTALW_OUTPUT_DIR;'})));
	print OUTFILE "CLUSTALW_OUTPUT_DIR" ."=". $val . "\n" if (defined($val));


    }

}#end sub write_replacement_vars_to_msf2bsml_configfile()




#-------------------------------------------------------------------------
# prepare_master_subsfile()
#
#-------------------------------------------------------------------------
sub prepare_master_subsfile {

    $logger->debug("Entered prepare_master_subsfile") if $logger->is_debug;

    my ($linkscores, $workflow_instance_dir, $set_runtime, $conf, $database_tmpdir, $asmb_hash) = @_;

    $logger->logdie("linkscores was not defined")            if (!defined($linkscores));
    $logger->logdie("workflow_instance_dir was not defined") if (!defined($workflow_instance_dir));
    $logger->logdie("set_runtime was not defined")           if (!defined($set_runtime));
    $logger->logdie("conf was not defined")                  if (!defined($conf));
    $logger->logdie("database_tmpdir was not defined")       if (!defined($database_tmpdir));
    $logger->logdie("asmb_hash was not defined")             if (!defined($asmb_hash));


    foreach my $score (sort @$linkscores){
	$logger->logdie("score was not defined") if (!defined($score));
	
	#
	# linkscores should be values between 0 and 1
	#
	if (($score < 0) or ($score >= 1)){
	    $logger->logdie("score $score must be 0 < score < 1");
	    next;
	}
	
	#
	# For each assembly group we will process all of the assembly identifiers 
	# associated to this assembly group
	#
	foreach my $group (sort keys %$asmb_hash){
	    $logger->logdie("group was not defined") if (!defined($group));
	    
	    my $jaccard = $group . "_" . $score;
	
	    $conf->{';CONFIG_LIST;'}   .= $workflow_instance_dir ."/" . $jaccard .".ini,";
	    $conf->{';INSTANCE_LIST;'} .= $workflow_instance_dir ."/" . $jaccard .".xml,";
	
	}    
    }
    #
    # Get rid of trailing comma
    #
    chop $conf->{';CONFIG_LIST;'};
    chop $conf->{';INSTANCE_LIST;'};
    
    $conf->{';WORKFLOW_INSTANCE_DIR;'} = $workflow_instance_dir;

    #
    # Store the master configuration for the bsml2chado workflow instances to follow
    # and set global permissions
    #
    my $master_wf_fullpath = $workflow_instance_dir ."/master.subs";
    &store_conf_to_file(\$master_wf_fullpath, $conf);

    my $execution_string = "chmod 666 $master_wf_fullpath";
    &execute(\$execution_string);

    #
    # Perform the substitution on the master bsml2chado workflow
    # and set global permissions
    # 
    $execution_string = $set_runtime ." -c ". $master_wf_fullpath ." < ". $conf->{';JACCARD_MASTER_CONF;'} ." > ". $workflow_instance_dir ."/jaccard_cluster-master-instance.ini";
    &execute(\$execution_string);

    $execution_string = "chmod 666 " . $workflow_instance_dir ."/jaccard_cluster-master-instance.ini";
    &execute(\$execution_string);

}


#-------------------------------------------------------------------
# get_file_contents()
#
#-------------------------------------------------------------------
sub get_file_contents {

    my $file = shift;
    if (!defined($file)){
	$logger->logdie("file was not defined");
    }

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
# is_file_status_ok()
#
#-------------------------------------------------------------------
sub is_file_status_ok {

    my $file = shift;
    my $fatal_flag=0;
    if (!defined($file)){
	$logger->fatal("file was not defined");
	$fatal_flag++;
    }
    if (!-e $$file){
	$logger->fatal("$$file does not exist");
	$fatal_flag++;
    }
    if (!-r $$file){
	$logger->fatal("$$file does not have read permissions");
	$fatal_flag++;
    }

    if ($fatal_flag>0){
	return 0;
    }
    else{
	return 1;
    }

}#end sub is_file_status_ok()

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

#-----------------------------------------------------------------
# retrieve_directory_contents()
#
#-----------------------------------------------------------------
sub retrieve_directory_contents{

    my $dir = shift;
    $logger->logdie("dir was not defined") if (!defined($$dir));
    $logger->logdie("dir:$dir is not a directory") if (!-d $$dir);

    my @bsml_files = qx"ls -1 $$dir/*.bsml";
    chomp @bsml_files;
    
    return \@bsml_files;

}#end sub retrieve_directory_contents()

#-----------------------------------------------------------------
# execute()
#
#-----------------------------------------------------------------
sub execute {
    
    my $string = shift;
    $logger->logdie("string was not defined") if (!defined($string));
    
    $logger->info("$$string");
    system $$string;
    
    
}#end sub execute()
    
#-----------------------------------------------------------------
# print_configuration()
#
#-----------------------------------------------------------------
sub print_configuration {

    my $conf = shift;
    if (!defined($conf)){
	$logger->logdie("conf was not defined");
    }

    foreach my $key (sort keys %$conf){
	if (!defined($key)){
	    $logger->logdie("key was not defined");
	}
	print ("$key=$conf->{$key}\n");
    }

    exit;

}#end sub print_configuration()



#------------------------------------------------------------------
# store_conf_to_file()
#
#------------------------------------------------------------------
sub store_conf_to_file {

    my ($file, $conf) = @_;
    $logger->logdie("file was not defined") if (!defined($file));
    $logger->logdie("conf was not defined") if (!defined($conf));

    open( SUBFILE, ">$$file" ) or $logger->logdie("Could not open $$file");

    foreach my $key (sort keys( %$conf )){
	$logger->logdie("key was not defined") if (!defined($key));

	print SUBFILE "$key=$conf->{$key}\n";
    }

    close( SUBFILE );
}

#-----------------------------------------------------------
# override_default_configuration()
#
#-----------------------------------------------------------
sub override_default_configuration {

    my $file = shift;
    my $confhash = shift;

    $logger->logdie("file was not defined") if (!defined($file));

    my $contents = &get_file_contents($file);
    $logger->logdie("The contents of $config were not retrieved, contents was not defined") if (!defined($contents));

    foreach my $line (@$contents){
	$logger->logdie("line was not defined") if (!defined($line));

	my ($key, $value) = split( '=', $line );
	$logger->logdie("key was not defined")   if (!defined($key));
	$logger->logdie("value was not defined") if (!defined($value));

	$confhash->{$key} = $value;
    }
}

#---------------------------------------------------------------------------------------------
# create_workflow_instance_subfolder()
# 
# Create the workflow archive directory within the BSML 
# repository e.g. mkdir -m 777 -p /usr/local/annotation/CHADO_TEST/Workflow/bsml2chado/3065
#
#---------------------------------------------------------------------------------------------
sub create_workflow_instance_subfolder {

    $logger->logdie("$conf->{';REPOSITORY_ROOT;'} is not a directory") if (!-d $conf->{';REPOSITORY_ROOT;'});

    my $execution_string;
    
    $execution_string = "mkdir -m 777 -p ". $conf->{';REPOSITORY_ROOT;'} ."/". $conf->{';DATABASE_UC;'} ."/Workflow";
    &execute(\$execution_string);
    
    $execution_string = "mkdir -m 777 -p ". $conf->{';REPOSITORY_ROOT;'} ."/". $conf->{';DATABASE_UC;'} ."/Workflow/" .$conf->{';WFNAME;'};
    &execute(\$execution_string);
    
    $execution_string = "mkdir -m 777 -p ". $conf->{';REPOSITORY_ROOT;'} ."/". $conf->{';DATABASE_UC;'} ."/Workflow/" .$conf->{';WFNAME;'} ."/". $conf->{';WFID;'};
    &execute(\$execution_string);
    
    my $workflow_instance_dir = $conf->{';REPOSITORY_ROOT;'} ."/". $conf->{';DATABASE_UC;'} ."/Workflow/" .$conf->{';WFNAME;'} ."/". $conf->{';WFID;'};
    

    return $workflow_instance_dir;
}

#--------------------------------------------------------------------
# retrieve_linkscores()
#
#
#--------------------------------------------------------------------
sub retrieve_linkscores {

    $logger->debug("Entered retrieve_linkscores") if $logger->is_debug;
    
    my $linkscore_list = shift;

    my @list = split(/,/, $linkscore_list);


    my @goodlist;

    foreach my $score (sort @list){
	$score =~ s/\s+//g;
	push (@goodlist, $score);
    }

    return \@goodlist;



}#end sub retrieve_linkscores()


#--------------------------------------------------------------------------
# run_workflow()
#
#--------------------------------------------------------------------------
sub run_workflow {

    $logger->debug("Entered run_workflow") if $logger->is_debug;

    my ($run_wf, $workflow_instance_dir) = @_;
    $logger->logdie("run_wf was not defined")                if (!defined($run_wf));
    $logger->logdie("workflow_instance_dir was not defined") if (!defined($workflow_instance_dir));

    my $execution_string = $run_wf ." -d ". $workflow_instance_dir ." -c ". $workflow_instance_dir ."/jaccard_cluster-master-instance.ini" . " -t ". $workflow_instance_dir ."/jaccard_cluster-master-instance_template.xml" . " -i ". $workflow_instance_dir ."/jaccard_cluster.xml" ." -l ". $workflow_instance_dir ."/log.txt";
    &execute(\$execution_string);
    
    $execution_string = "chmod 666 ". $workflow_instance_dir . "/*.xml";
    &execute(\$execution_string);
    
}#end sub run_workflow()

#----------------------------------------------------------------------
# retrieve_asm_hash()
#
#
#----------------------------------------------------------------------
sub retrieve_asm_hash {

    $logger->debug("Entered retrieve_asm_hash") if $logger->is_debug;

    my ($list, $file, $group) = @_;
    

    my $asmbl_hash = {};


    if ($list){
	$logger->logdie("group was not defined, but list was") if (!defined($group));

	#
	# asmbl_list was defined
	#
	my @tmplist = split(/,/,$list);
	foreach my $asmbl (@tmplist){

	    $logger->logdie("asmbl was not defined") if (!defined($asmbl));
	    
	    $asmbl_hash->{$group}->{$asmbl} = $asmbl if (!exists $asmbl_hash->{$group}->{$asmbl});
	}
    }

    if ($file){

	#
	# asmbl_file was defined
	#
	my $contents = &get_file_contents(\$file);
	$logger->logdie("contents was not defined") if (!defined($contents));

	my $filegroup;

	foreach my $line (@$contents){
	    $logger->logdie("line was not defined") if (!defined($line));
	 

	    if ($line =~ /^$/){
		next;
	    }

	    if ($line =~ /^group:(.+)$/){
		$filegroup = $1;
	    }
	    elsif ($line =~ /_assembly\s*$/){
		if (defined($filegroup)){
		    $asmbl_hash->{$filegroup}->{$line} = $line if (!exists $asmbl_hash->{$filegroup}->{$line});
		}
		else{
		    $logger->logdie("filegroup was not defined while processing line:$line");
		}
	    }
	    else{
		$logger->logdie("un-recognized line:$line");
	    }
	}
    }






#    my $asm_string;

#    foreach my $assembly (sort keys %$asmbl_hash){
#	$logger->logdie("assembly was not defined") if (!defined($assembly));
#	$asm_string .= $assembly .",";
#    }
#    chop $asm_string;
#   return $asm_string;


    return $asmbl_hash;


}#end sub retrieve_asm_hash()


#-----------------------------------------------------------------
# print_usage()
#
#-----------------------------------------------------------------
sub print_usage {



    print STDERR "SAMPLE USAGE:  $0 -U username -P password -D database [-a asmbl_list|ALL] [-F asmbl_file] [-b compute_type] [-S server] [-i pidentity] [-r pvalue] [-l linkscore_list] [-c config] [-h] [-m] [-p] [-v] [-w]\n";
    print STDERR "  -U|--username          = database login username\n";
    print STDERR "  -P|--password          = database login password\n";
    print STDERR "  -D|--database          = name of chado database to load analysis into\n";
    print STDERR "  -a|--asmbl_list        = comma-separated list of assembly identifiers or \"ALL\"\n";
    print STDERR "  -F|--asmbl_file        = file containing new-line separated list of assembly identifiers\n";
    print STDERR "  -b|--compute_type      = compute_type of the source pairwise alignment either \"blastp\" OR \"ber\"\n";
    print STDERR "  -S|--server            = Sybase server\n";
    print STDERR "  -i|--pidentity         = percent identity threshold\n";
    print STDERR "  -r|--pvalue            = pvalue threshold\n";
    print STDERR "  -l|--linkscore_list    = comma-separated list of Jaccard coefficient link score cut-off values\n";
    print STDERR "  -c|--config_file       = configuration file\n";
    print STDERR "  -h|--help              = displays this usage message\n";
    print STDERR "  -m|--man               = displays the pod2usage page for this utility\n";
    print STDERR "  -p|--printconf         = print out the key=value configuration pairs\n";
    print STDERR "  -v|--verbose           = Increases log4perl reporting level to screen from WARN to INFO\n";
    print STDERR "  -w|--workflow_monitor  = automatically launch the workflow monitor\n";
    exit 1;

}#end sub print_usage()

