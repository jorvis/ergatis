#!/usr/local/bin/perl
#---------------------------------------------------------------------
# script name: run_bsml2chado.pl
# date:        2003.10.28
#
#
#
#
#---------------------------------------------------------------------

=head1 NAME

run_bsml2chado.pl - Creates and executes the bsml to chado migration workflow

=head1 SYNOPSIS

USAGE:  run_bsml2chado.pl -U username -P password -D target_database [-t target_server] [-f organism_file|-o organism_list] [-b bsml_type] [-f bsml_batch_file] [-c config] [-h] [-m] [-p] [-q] [-v]

=head1 OPTIONS

=over 8

=item B<--username,-U>
    
    Database username

=item B<--password,-P>
    
    Database password

=item B<--target_database,-D>
    
    Source database name

=item B<--target_server,-t>
    
    target server

=item B<--bsml_batch_file,-f>
    
    File containing list of fullpath to bsml documents.  Use the -q switch to view a sample bsml_batch_file.

=item B<--bsml_list,-o>

    Comma separated list of fullpath's of bsml documents

=item B<--bsml_type,-b>

    One of three: 1=pairwise alignment 2=multiple alignment 3=SNP 4=custom searches 5=sub-assembly

=item B<--config,-c>

    Configuration file containing all substitution keys parameters

=item B<--help,-h>

    Prints out usage

=item B<--printconf,-p>

    Prints out the default configuration key-value pairs

=item B<--bsml_format,-q>

    Prints out sample bsml_batch_file

=item B<--man,-m>

    Displays this pod2usage man pages for this script


=item B<--verbose,-v>

    Increases log4perl reporting level to screen from WARN to INFO


=back

=head1 DESCRIPTION

    run_bsml2chado.pl - Creates and executes the bsml to chado migration workflow

=cut




use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Data::Dumper;
use Log::Log4perl qw(get_logger);
use File::Basename;
use Pod::Usage;


my $WorkflowConfigFileDir = $ENV{'WORKFLOW_DOCS_DIR'};
my $WorkflowExecDir       = $ENV{'WORKFLOW_WRAPPERS_DIR'};

#
# The default bsml2chado configuration
#
my $conf = {
    ';INSTALLATION;'                 => $WorkflowExecDir,
    ';TMP_DIR;'                      => '/usr/local/scratch/bsml2chado',
    ';REPOSITORY_ROOT;'              => '/usr/local/annotation',
    ';WORKFLOW_DIR;'                 => $WorkflowExecDir,
    ';SET_RUNTIME;'                  => $WorkflowExecDir ."/set_runtime_vars.pl",
    ';BSML2CHADO_CONF;'              => $WorkflowConfigFileDir ."/bsml2chado.ini", 
    ';BSML2CHADO_MASTER_CONF;'       => $WorkflowConfigFileDir ."/bsml2chado-master.ini", 
    ';BSML2CHADO_TEMPLATE;'          => $WorkflowConfigFileDir ."/bsml2chado_template.xml",
    ';BSML2CHADO_MASTER_TEMPLATE;'   => $WorkflowConfigFileDir ."/bsml2chado-master_template.xml",
    ';RUN_WF;'                       => $WorkflowExecDir ."/run_wf.sh",
    ';WFNAME;'                       => 'bsml2chado',	
    ';USERNAME;'                     => '',
    ';PASSWORD;'                     => '',               
    ';TARGET_DATABASE;'              => '',               
    ';TARGET_SERVER;'                => 'SYBTIGR',
    ';BSML_BATCH_FILE;'              => '',
    ';BSML_TYPE;'                    => '',
    ';BSML_DIRECTORY;'               => '',
};


my %options = ();
my ($bsml_batch_file, $target_database, $config, $verbose, $username, $password, $target_server, $bsml_type, $bsml_directory, $printconf, $help, $man, $workflow_monitor, $bsml_format);
my $results = GetOptions (
			  \%options,
			  'username|U=s'            => \$username,
			  'password|P=s'            => \$password,
			  'target_database|D=s'     => \$target_database,
			  'target_server|t=s'       => \$target_server,
			  'bsml_batch_file|f=s'     => \$bsml_batch_file,
			  'bsml_type|b=s'           => \$bsml_type,
			  'bsml_directory|d=s'      => \$bsml_directory,
			  'config|c=s'              => \$config,
			  'help|h'                  => \$help,
			  'man|m'                   => \$man,
			  'printconf|p'             => \$printconf,
			  'bsml_format|q'           => \$bsml_format,
			  'workflow_monitor|w'      => \$workflow_monitor,
			  'verbose|v'               => \$verbose,
			  );

&print_configuration($conf) if ($printconf);
&pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if ($man);
&print_usage() if ($help);
&sample_bsml_batch_file() if ($bsml_format);

my $log4perl = "bsml2chado_wf.log";

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
if (!defined($target_database)){
    $target_database = $conf->{';TARGET_DATABASE;'} if ((exists $conf->{';TARGET_DATABASE;'}) and (defined($conf->{';TARGET_DATABASE;'})));
}
if (!defined($target_server)){
    $target_server = $conf->{';TARGET_SERVER;'} if ((exists $conf->{';TARGET_SERVER;'}) and (defined($conf->{';TARGET_SERVER;'})));
}
if (!defined($bsml_batch_file)){
    $bsml_batch_file = $conf->{';BSML_BATCH_FILE;'} if ((exists $conf->{';BSML_BATCH_FILE;'}) and (defined($conf->{';BSML_BATCH_FILE;'})));
}
if (!defined($bsml_type)){
    $bsml_type = $conf->{';BSML_TYPE;'} if ((exists $conf->{';BSML_TYPE;'}) and (defined($conf->{';BSML_TYPE;'})));
}
if (!defined($bsml_directory)){
    $bsml_directory = $conf->{';BSML_DIRECTORY;'} if ((exists $conf->{';BSML_DIRECTORY;'}) and (defined($conf->{';BSML_DIRECTORY;'})));
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
print STDERR "\nusername was not defined\n\n" if (!defined($username));
print STDERR "\npassword was not defined\n\n" if (!defined($password));
print STDERR "\ntarget_database was not defined\n\n" if (!defined($target_database));

&print_usage() if (!$username or !$password or !$target_database);

#
# Either bsml_directory or bsml_batch_file must be specified
#
if (!$bsml_directory && !$bsml_batch_file){
       print STDERR "\nYou must EITHER specify the bsml_batch_file\nOR\nthe bsml_directory AND the bsml_type\nwhere the bsml types are:\n1=pairwise alignment\n2=multiple alignment\n3=SNP\n4=custom searches\n5=sub-assembly\n\nFor an explanation on how to format the bsml_batch_file invoke this script as:\n$0 -q\n\n";
    &print_usage();
}

#
# If user specifies the bsml_directory, must also specify the bsml_type
#
if (($bsml_directory) and (!$bsml_type)){
    print STDERR "\nYou must specify the bsml_type with the bsml_directory\n\n";
    &print_usage();
}

#
# User must specify the bsml type if they don't specify the bsml_batch_file
#
if (!$bsml_type and !$bsml_batch_file){
    print STDERR "\nYou must specify the bsml type \n1=pairwise alignment\n2=multiple alignment\n3=SNP\n4=custom searches\n5=sub-assembly\n\n";
    &print_usage();
}

#
# In case user wishes to specify a bsml_batch_file but accidentally uses -d instead of -f
#
if (($bsml_directory) and (!-d $bsml_directory)){
    print STDERR "\nVerify whether $bsml_directory is a valid directory\nPlease note that \"-d\" is used to specify a bsml directory\nand \"-f\" is used to specify a filename which contains a list of bsml documents.\n\n";
    &print_usage();
}

#
# Determine whether target_server was defined and whether is either SYBIL or SYBTIGR
#
if (!defined($target_server)){
    print STDERR "target_server was not defined";
    &print_usage();
}
if ($target_server !~ /^SYBIL|SYBTIGR$/){
    print STDERR "Invalid Sybase target_server: $target_server, must be either SYBIL or SYBTIGR";
    &print_usage();
}


my $bsml2chado_program = &determine_bsml2chado_program($bsml_type) if ($bsml_type);


#
# Get rid of trailing forward slashes
#
if ($bsml_directory){
    $bsml_directory =~ s/\/+$//;
#$bsml_directory = "." if (!defined($bsml_directory));
    $logger->logdie("$bsml_directory is not a directory") if (!-d $bsml_directory);
    $logger->info("directory is set to:$bsml_directory");
}


#
# Set additional configuration values
#
$conf->{';USERNAME;'}           = $username           if (defined($username));
$conf->{';PASSWORD;'}           = $password           if (defined($password));
$conf->{';BSML2CHADO_PROGRAM;'} = $bsml2chado_program if (defined($bsml2chado_program));
$conf->{';TARGET_SERVER;'}      = $target_server             if (defined($target_server));
$conf->{';TARGET_DATABASE_LC;'} = lc($target_database)       if (defined($target_database));
$conf->{';TARGET_DATABASE_UC;'} = uc($target_database)       if (defined($target_database));
$conf->{';WFID;'}               = $$;


#---------------------------------------------------------------------------------------------
# Create the workflow archive directory within the BSML 
# repository e.g. mkdir -m 777 -p /usr/local/annotation/CHADO_TEST/Workflow/bsml2chado/3065
#
#---------------------------------------------------------------------------------------------
$logger->logdie("$conf->{';REPOSITORY_ROOT;'} is not a directory") if (!-d $conf->{';REPOSITORY_ROOT;'});

my $execution_string;

$execution_string = "mkdir -m 777 -p ". $conf->{';REPOSITORY_ROOT;'} ."/". $conf->{';TARGET_DATABASE_UC;'} ."/Workflow";
&execute(\$execution_string);

$execution_string = "mkdir -m 777 -p ". $conf->{';REPOSITORY_ROOT;'} ."/". $conf->{';TARGET_DATABASE_UC;'} ."/Workflow/" .$conf->{';WFNAME;'};
&execute(\$execution_string);

$execution_string = "mkdir -m 777 -p ". $conf->{';REPOSITORY_ROOT;'} ."/". $conf->{';TARGET_DATABASE_UC;'} ."/Workflow/" .$conf->{';WFNAME;'} ."/". $conf->{';WFID;'};
&execute(\$execution_string);

my $workflow_instance_dir = $conf->{';REPOSITORY_ROOT;'} ."/". $conf->{';TARGET_DATABASE_UC;'} ."/Workflow/" .$conf->{';WFNAME;'} ."/". $conf->{';WFID;'};

#---------------------------------------------------------------------------------------------
# Retrieve the bsml_document hash
#
#---------------------------------------------------------------------------------------------
my $bsml_hash = &retrieve_bsml_document_hash($bsml_batch_file, $bsml_directory, $bsml_type);
$logger->logdie("bsml_document_hash was not defined") if (!defined($bsml_hash));

#---------------------------------------------------------------------------------------------
# For the master workflow, append each bsml to the config and instance list
#
#---------------------------------------------------------------------------------------------
foreach my $bsml_item (sort keys %$bsml_hash){
    
    $logger->logdie("bsml_item was not defined") if (!defined($bsml_item));

    my $bsml_basename = $bsml_hash->{$bsml_item}->{'bsml_basename'} if (exists $bsml_hash->{$bsml_item}->{'bsml_basename'});
    $logger->logdie("bsml_basename was not defined") if (!defined($bsml_basename));
    
    $conf->{';CONFIG_LIST;'}   .= $workflow_instance_dir ."/" . $bsml_basename .".ini,";
    $conf->{';INSTANCE_LIST;'} .= $workflow_instance_dir ."/" . $bsml_basename .".xml,";

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

$execution_string = "chmod 666 $master_wf_fullpath";
&execute(\$execution_string);

#
# Perform the substitution on the master bsml2chado workflow
# and set global permissions
#
$execution_string = $set_runtime ." -c ". $master_wf_fullpath ." < ". $conf->{';BSML2CHADO_MASTER_CONF;'} ." > ". $workflow_instance_dir ."/bsml2chado-master-instance.ini";
&execute(\$execution_string);

$execution_string = "chmod 666 " . $workflow_instance_dir ."/bsml2chado-master-instance.ini";
&execute(\$execution_string);


#------------------------------------------------------------------------------------------
# For each bsml stored in the bsml_document_hash, generate a subflow
#
#------------------------------------------------------------------------------------------
foreach my $bsml_item (sort keys %$bsml_hash){

    $logger->logdie("bsml_item was not defined") if (!defined($bsml_item));
	
    #
    # Extract the basename i.e. bsp_3839_assembly.blastp.bsml
    #
    my $bsml_basename = $bsml_hash->{$bsml_item}->{'bsml_basename'} if (exists $bsml_hash->{$bsml_item}->{'bsml_basename'});
    $logger->logdie("bsml_basename was not defined") if (!defined($bsml_basename));

    $conf->{';BSML_DOCUMENT;'} = $bsml_basename;

    #
    # Determine the BSML program based on BSML type by extracting the bsml_type if available
    #
    if (exists $bsml_hash->{$bsml_item}->{'bsml_type'}){
	$bsml2chado_program = &determine_bsml2chado_program($bsml_hash->{$bsml_item}->{'bsml_type'});
	$conf->{';ANALYSIS_TYPE;'} = $bsml_hash->{$bsml_item}->{'bsml_type'};
    }

    $conf->{';BSML2CHADO_PROGRAM;'} = $bsml2chado_program;

    #
    # Extract the fullpath for the BSML document
    #
    my $bsml_fullpath = $bsml_hash->{$bsml_item}->{'bsml_fullpath'} if (exists $bsml_hash->{$bsml_item}->{'bsml_fullpath'});
    $logger->logdie("bsml_fullpath was not defined for $bsml_hash->{$bsml_item}->{'bsml_basename'}") if (!defined($bsml_fullpath));
    
    $conf->{';BSML_DOCUMENT_FULLPATH;'} = $bsml_fullpath;


    #
    # Create a subdirectory in the scratch space for this bsml 
    #
#    my $execution_string = "mkdir -p -m 777 " . $conf->{';TMP_DIR;'};
#    &execute(\$execution_string);
#    $execution_string = "mkdir -p -m 777 " . $conf->{';TMP_DIR;'} ."/". $conf->{';TARGET_DATABASE_LC;'};
#    &execute(\$execution_string);
#    $execution_string = "mkdir -p -m 777 " . $conf->{';TMP_DIR;'} ."/". $conf->{';TARGET_DATABASE_LC;'} . "/" . $bsml_basename;
#    &execute(\$execution_string);


    #
    # Create a subdirectory in the workflow instance directory for this bsml
    # We will store all the associated goodies in this bsml specific directory
    #
    $execution_string = "mkdir -p -m 777 ". $workflow_instance_dir;
    &execute(\$execution_string);
    $execution_string = "mkdir -p -m 777 ". $workflow_instance_dir  . "/" . $bsml_basename;
    &execute(\$execution_string);

    #
    # All of the *.log and *.stats files necessary for the migration
    #
    $conf->{';BSML2CHADO_LOGFILE;'}    = $workflow_instance_dir ."/". $bsml_basename ."/bsml2chado.log";
    $conf->{';LOAD_LOGFILE;'}          = $workflow_instance_dir ."/". $bsml_basename ."/load.log";
    $conf->{';LOAD_STATS;'}            = $workflow_instance_dir ."/". $bsml_basename ."/load.stats";
    $conf->{';DBSPACE_LOGFILE;'}       = $workflow_instance_dir ."/". $bsml_basename ."/dbspace.log";
    $conf->{';VALIDATION_LOGFILE;'}    = $workflow_instance_dir ."/". $bsml_basename ."/validation.log";
    $conf->{';VALIDATION_OUTDIR;'}     = $workflow_instance_dir ."/". $bsml_basename;
    $conf->{';CONFIRM_LOAD_LOGFILE;'}  = $workflow_instance_dir ."/". $bsml_basename ."/confirm_load.log";
    $conf->{';CONFIRM_LOAD_OUTDIR;'}   = $workflow_instance_dir ."/". $bsml_basename;
    $conf->{';LOG_DIRECTORY;'}         = $workflow_instance_dir ."/". $bsml_basename;

    
    #
    # Store the configuration for this specific bsml
    # and set global permissions
    #
    my $bsml_subflow = $workflow_instance_dir ."/". $bsml_basename ."/". $bsml_basename .".sub";
    &store_conf_to_file(\$bsml_subflow, $conf);
    
    $execution_string = "chmod 666 $bsml_subflow";
    &execute(\$execution_string);

    #
    # Perform sed replacement on the specific bsml subflow configuration .ini file 
    #
    $execution_string = $set_runtime ." -c ". $bsml_subflow ." < ". $conf->{';BSML2CHADO_CONF;'} ." > ". $workflow_instance_dir ."/". $bsml_basename .".ini";
    &execute(\$execution_string);

    #
    # Give the subflow instance configuration file global permissions
    #
    $execution_string = "chmod 666 " . $workflow_instance_dir ."/". $bsml_basename .".ini";
    &execute(\$execution_string);
    
}

#
# Copy the bsml2chado_template.xml        to bsml2chado-instance_template.xml AND
#          bsml2chado-master_template.xml to bsml2chado-master-instance_template.xml
#
# and then set global permissions
#
$execution_string = "cp ". $conf->{';BSML2CHADO_TEMPLATE;'} ." ". $workflow_instance_dir ."/bsml2chado-instance_template.xml";
&execute(\$execution_string);
$execution_string = "cp ". $conf->{';BSML2CHADO_MASTER_TEMPLATE;'} ." ". $workflow_instance_dir ."/bsml2chado-master-instance_template.xml";
&execute(\$execution_string);

$execution_string = "chmod 666 ". $workflow_instance_dir . "/bsml2chado-instance_template.xml";
&execute(\$execution_string);
$execution_string = "chmod 666 ". $workflow_instance_dir . "/bsml2chado-master-instance_template.xml";
&execute(\$execution_string);
    

#
# Run the bsml2chado workflow
#
$execution_string = $run_wf ." -d ". $workflow_instance_dir ." -c ". $workflow_instance_dir ."/bsml2chado-master-instance.ini" . " -t ". $workflow_instance_dir ."/bsml2chado-master-instance_template.xml" . " -i ". $workflow_instance_dir ."/bsml2chado.xml" ." -l ". $workflow_instance_dir ."/log.txt";
&execute(\$execution_string);

$execution_string = "chmod 666 ". $workflow_instance_dir . "/*.xml";
&execute(\$execution_string);


print STDERR "\n$0 execution complete\nPlease review log4perl log file:$log4perl\n\n";
#----------------------------------------------------------------------------------------------------------------------------
#
#                          END OF MAIN -- SUBROUTINES FOLLOW
#
#----------------------------------------------------------------------------------------------------------------------------



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


#------------------------------------------------------------------------
# determine_bsml2chado_program()
#
#------------------------------------------------------------------------
sub determine_bsml2chado_program {

    my $script = shift;
    $logger->logdie("script was not defined") if (!defined($script));
    $logger->logdie("Invalid bsml_type option: $script\nPlease see the man pages for the utility") if ($script !~ /^[1|2|3|4|5]$/);
	
    my $program;

    if ($script eq '1'){
	$program = 'bsmlPairwiseAlignment2Chado';
    }
    elsif ($script eq '2'){
	$program = 'bsmlMultipleAlignment2Chado';
    }
    elsif ($script eq '3'){
	$program = 'bsmlSNP2Chado';
    }
    elsif ($script eq '4'){
	$program  = 'bsmlCustomSearches2Chado';
    }
    elsif ($script eq '5'){
	$program = 'bsmlSubAssembly2Chado';
    }
    else{
	$logger->logdie("Unrecognized option: $script");
    }

    $logger->info("bsml2chado program set to:$program");
    return $program;
}


#---------------------------------------------------------------------------
# retrieve_bsml_document_hash()
#
#---------------------------------------------------------------------------
sub retrieve_bsml_document_hash {
    
    my ($bsml_batch_file, $bsml_directory, $bsml_type) = @_;
    my $bsml_hash={};
    
    #
    # If a bsml_directory was specified, retrieve all BSML documents contained in the
    # specified directory and add unique bsml document names tot he bsml_document_hash
    #
    if ($bsml_directory){
	if ($bsml_type){

	    my $bsml_listref = &retrieve_directory_contents(\$bsml_directory);
	    #
	    # E.g.
	    # $VAR1 = [
	    #           '/usr/local/annotation/CHADO_TEST/BSML_repository/blastp//bsp_3839_assembly.blastp.bsml',
	    #           '/usr/local/annotation/CHADO_TEST/BSML_repository/blastp//gbs18rs21_2501_assembly.blastp.bsml',
	    #           '/usr/local/annotation/CHADO_TEST/BSML_repository/blastp//gbs_799_assembly.blastp.bsml',
	    #           '/usr/local/annotation/CHADO_TEST/BSML_repository/blastp//gmi_1459_assembly.blastp.bsml',
	    #           '/usr/local/annotation/CHADO_TEST/BSML_repository/blastp//gpa_1200_assembly.blastp.bsml'
	    #         ];
	    #
	    $logger->logdie("bsml_listref was not defined") if (!defined($bsml_listref));
	    
	    foreach my $bsml_document (@$bsml_listref){
		#
		# E.g.
		# /usr/local/annotation/CHADO_TEST/BSML_repository/blastp//bsp_3839_assembly.blastp.bsml
		#
		$logger->logdie("bsml_document was not defined") if (!defined($bsml_document));
		
		my $bsml_fullpath = $bsml_document;
		my $bsml_basename = basename($bsml_document);
		#
		# E.g. bsml_basename = bsp_3839_assembly.blastp.bsml
		#
		if (!exists $bsml_hash->{$bsml_fullpath}){
		    if (&is_file_readable(\$bsml_fullpath)){
			#
			# BSML document exists and has read permissions, therefore load into bsml_hash
			#
			$bsml_hash->{$bsml_fullpath}->{'bsml_fullpath'} = $bsml_fullpath;
			$bsml_hash->{$bsml_fullpath}->{'bsml_basename'} = $bsml_basename;
			$bsml_hash->{$bsml_fullpath}->{'bsml_type'}     = $bsml_type;
		    }
		    else{
			#
			# Either the BSML document does not exist or exists but does not have read permissions
			# and therefore shall not be processed.
			#
			$logger->fatal("BSML document:$bsml_fullpath will not be processed");
		    }
		}
		else{
		    #
		    # Specified BSML document name was cited more than one time, though will only be processed once.
		    # Kindly inform the user
		    #
		    $logger->warn("BSML document:$bsml_fullpath has already been flagged for processing\n");
		}
	    }
	}
	else{
	    #
	    # User specified the directory containing BSML document, but did not specify the type of documents.  
	    # Kindly inform the user of their mistake...
	    #
	    $logger->logdie("The bsml_type was not specified.  Unable to determine what type of BSML documents are stored in directory:$bsml_directory\nPlease see the man page for this utility ($0 -m).\n\n");
	}
    }
    
    
    #
    # If bsml_batch_file was specified, retrieve all of the bsml documents listed in the bsml file
    # and add unique bsml document names to the bsml_document_hash
    #
    if ($bsml_batch_file){
	
	my $contents = &get_file_contents(\$bsml_batch_file);
	$logger->logdie("contents was not defined") if (!defined($contents));
	my $linectr=0;
	
	foreach my $line (@$contents){
	    $logger->logdie("line was not defined") if (!defined($line));
	    $linectr++;
	    if ($line =~ /^bsml_doc:(\S+\.bsml)\s*(bsml_type:\d)?$/){
		#
		# E.g.  line = bsml_doc:/usr/local/annotation/CHADO_TEST/BSML_repository/blastp//bsp_3839_assembly.blastp.bsml
		#
		
		my $bsml_fullpath = $1;
		my $bsml_basename = basename($bsml_fullpath);
		my $proposed_bsml_type = $2;
		my $specific_bsml_type;


		if (defined($proposed_bsml_type)){
		    if ($proposed_bsml_type =~ /^bsml_type:(\d)$/){
			$logger->info("\nSetting BSML document type to $1");
			$specific_bsml_type = $1;
		    }
		    else{
			$logger->warn("Could not extract BSML document type from line:$line");
			if (defined($bsml_type)){
			    $logger->info("\nSetting BSML document type to $bsml_type");
			    $specific_bsml_type = $bsml_type;
			}
			else{
			    $logger->logdie("BSML document type variable bsml_type is not defined");
			}
		    }
		}
		else{
		    if ($bsml_type){
			$logger->info("\nBSML document type was not stored in the bsml_batch_file, therefore setting BSML document type to $bsml_type");
			$specific_bsml_type = $bsml_type;
		    }
		    else{
			$logger->error("BSML document type was not specified in the bsml_batch_file, nor on the command-line --bsml_type\nPlease take the time to understand the usage.\n");
			&sample_bsml_batch_file();
		    }
		}
		#
		# E.g. bsml_basename = bsp_3839_assembly.blastp.bsml
		#
		if (!exists $bsml_hash->{$bsml_fullpath}){
		    #
		    # Only process unique bsml documents
		    #
		    if (&is_file_readable(\$bsml_fullpath)){
			#
			# BSML document exists and has read permissions, therefore load into bsml_hash
			#
			$bsml_hash->{$bsml_fullpath}->{'bsml_fullpath'} = $bsml_fullpath;
			$bsml_hash->{$bsml_fullpath}->{'bsml_basename'} = $bsml_basename;
			$bsml_hash->{$bsml_fullpath}->{'bsml_type'}     = $specific_bsml_type;
			
		    }
		    else{
			#
			# Either the BSML document does not exist or exists but does not have read permissions
			# and therefore shall not be processed.
			#
			$logger->fatal("BSML document:$bsml_fullpath will not be processed");
		    }
		}
		else{
		    #
		    # Specified BSML document name was cited more than one time, though will only be processed once.
		    # Kindly inform the user
		    #
		    $logger->warn("BSML document:$bsml_basename has already been flagged for processing\nThe full path to the BSML document was $bsml_hash->{$bsml_fullpath}->{'bsml_fullpath'}\nThe BSML type was set to $bsml_hash->{$bsml_fullpath}->{'bsml_type'}\n");
		}
	    }
	    else{
		#
		# line could not be parsed.  Indicate to user proper formatting for bsml_batch_file
		#
		$logger->fatal("Could not parse line $linectr in file $bsml_batch_file.  Line was:\n$line");
		&sample_bsml_batch_file();
	    }
	}
    }#end if ($bsml_batch_file)
	
	return $bsml_hash;

}#end sub retrieve_bsml_document_hash()


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

#-----------------------------------------------------------------
# sample_bsml_batch_file()
#
#-----------------------------------------------------------------
sub sample_bsml_batch_file {

    print STDERR "\n";
    print STDERR "Each line of the bsml_batch_file should contain the fullpath to a valid BSML document to be loaded AND the BSML document type.\n";
    print STDERR "Note that if the BSML document type is not specified on each line in the file, then the --bsml_type argument on the command-line should be\n";
    print STDERR "specified by the user.  Thus all of the BSML documents specified in the --bsml_batch_file will be treated as the same --bsml_type\n";
    print STDERR "\nHere is a listing of the different bsml types:\n";
    print STDERR "a)\t1=pairwise alignment\n";
    print STDERR "b)\t2=multiple alignment\n";
    print STDERR "c)\t3=SNP\n";
    print STDERR "d)\t4=customer searches\n";
    print STDERR "e)\t5=sub-assembly\n";
    print STDERR "\nUse-case 1: bsml_type specified in the bsml_batch_file\n";
    print STDERR "This is an example of how the bsml_batch_file would be formatted:\n";
    print STDERR "-----------------------------------------------------------------------------------------------\n";
    print STDERR "bsml_doc:/usr/local/annotation/PNEUMO/BSML_repository/all_vs_all/gmi_1459_assembly.allvsall.bsml  bsml_type:1\n";
    print STDERR "bsml_doc:/usr/local/annotation/PNEUMO/BSML_repository/all_vs_all/ntsp01_1_assembly.allvsall.bsml  bsml_type:1\n";
    print STDERR "bsml_doc:/usr/local/annotation/PNEUMO/BSML_repository/all_vs_all/ntsp02_2_assembly.allvsall.bsml  bsml_type:1\n";
    print STDERR "bsml_doc:/usr/local/annotation/PNEUMO/BSML_repository/all_vs_all/gbs_799_assembly.allvsall.bsml  bsml_type:1\n";
    print STDERR "bsml_doc:/usr/local/annotation/PNEUMO/BSML_repository/all_vs_all/gpa_1200_assembly.allvsall.bsml  bsml_type:1\n";
    print STDERR "bsml_doc:/usr/local/annotation/PNEUMO/BSML_repository/all_vs_all/bsp_3839_assembly.allvsall.bsml  bsml_type:1\n";
    print STDERR "-----------------------------------------------------------------------------------------------\n";
    print STDERR "\nUse-case 2: bsml_type not specified in the bsml_batch_file, therefore user must specify the command-line argument --bsml_type\n";
    print STDERR "This is an example of how the bsml_batch_file would be formatted:\n";
    print STDERR "-----------------------------------------------------------------------------------------------\n";
    print STDERR "bsml_doc:/usr/local/annotation/PNEUMO/BSML_repository/all_vs_all/gmi_1459_assembly.allvsall.bsml\n";
    print STDERR "bsml_doc:/usr/local/annotation/PNEUMO/BSML_repository/all_vs_all/ntsp01_1_assembly.allvsall.bsml\n";
    print STDERR "bsml_doc:/usr/local/annotation/PNEUMO/BSML_repository/all_vs_all/ntsp02_2_assembly.allvsall.bsml\n";
    print STDERR "bsml_doc:/usr/local/annotation/PNEUMO/BSML_repository/all_vs_all/gbs_799_assembly.allvsall.bsml\n";
    print STDERR "bsml_doc:/usr/local/annotation/PNEUMO/BSML_repository/all_vs_all/gpa_1200_assembly.allvsall.bsml\n";
    print STDERR "bsml_doc:/usr/local/annotation/PNEUMO/BSML_repository/all_vs_all/bsp_3839_assembly.allvsall.bsml\n";
    print STDERR "-----------------------------------------------------------------------------------------------\n";
    print STDERR "\n";
    print STDERR "\n";
    exit;
}		      




#-----------------------------------------------------------------
# print_usage()
#
#-----------------------------------------------------------------
sub print_usage {

    print STDERR "SAMPLE USAGE:  $0  -U username -P password -t target_server -D target_database [-b bsml_type] [-d bsml_directory|-f bsml_batch_file] [-c config_file] [-m] [-p] [-q] [-v]\n";
    print STDERR "  -U|--username          = target_database login username\n";
    print STDERR "  -P|--password          = target_database login password\n";
    print STDERR "  -t|--target_server     = target server\n";
    print STDERR "  -D|--target_database   = name of target_database to load analysis into\n";
    print STDERR "  -l|--log4perl          = log4perl log file\n";
    print STDERR "  -b|--bsml_type         = 1=pairwise alignment, 2=multiple alignment, 3=SNP, 4=custom searches, 5=sub-assembly\n";
    print STDERR "  -f|--bsml_batch_file   = file containing list of similar type BSML documents to be processed\n";
    print STDERR "  -d|--bsml_directory    = directory containing similar type BSML documents to be processed\n";
    print STDERR "  -c|--config_file       = configuration file\n";
    print STDERR "  -p|--printconf         = print out the key=value configuration pairs\n";
    print STDERR "  -q|--bsml_format       = print out format rules for the bsml_batch_file\n";
    print STDERR "  -v|--verbose           = Increases log4perl reporting level to screen from WARN to INFO\n";
    exit 1;

}#end sub print_usage()

