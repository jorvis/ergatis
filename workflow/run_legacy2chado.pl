#!/usr/local/bin/perl
#---------------------------------------------------------------------
# script name: run_legacy2chado.pl
# date:        2003.10.28
#
#
#
#
#---------------------------------------------------------------------


=head1 NAME

run_legacy2chado.pl - Creates and executes the legacy to chado migration workflow

=head1 SYNOPSIS

USAGE:  run_legacy2chado.pl -U username -P password -D target_database [-S server] [-t target_server] [-f organism_file|-o organism_list] [-b organism_type] [-c config] [-h] [-m] [-p] [-q] [-v]

=head1 OPTIONS

=over 8

=item B<--username,-U>
    
    Database username

=item B<--password,-P>
    
    Database password

=item B<--target_database,-D>
    
    Source database name

=item B<--server,-S>
    
    Source server

=item B<--target_server,-t>
    
    target server

=item B<--organism_file,-f>
    
    File containing list of organisms, organism type and asmbl_ids.  Use the -h switch to view a sample organism_file.

=item B<--organism_list,-o>

    Comma separated list of organisms' legacy database name

=item B<--organism_type,-b>

    One of three: euk,prok or ntprok

=item B<--config,-c>

    Configuration file containing all substitution keys parameters

=item B<--help,-h>

    Prints out usage

=item B<--printconf,-p>

    Prints out the default configuration key-value pairs

=item B<--organism_format,-q>

    Prints out sample organism file

=item B<--man,-m>

    Displays this pod2usage man pages for this script

=item B<--verbose,-v>

    Increases log4perl reporting level to screen from WARN to INFO

=back

=head1 DESCRIPTION

    run_legacy2chado.pl - Creates and executes the legacy to chado migration workflow

=cut







use Pod::Usage;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Data::Dumper;
use Log::Log4perl qw(get_logger);
use File::Basename;



my $WorkflowConfigFileDir = $ENV{'WORKFLOW_DOCS_DIR'};
my $WorkflowExecDir       = $ENV{'WORKFLOW_WRAPPERS_DIR'};

#
# The default bsml2chado configuration
#
my $conf = {
    ';INSTALLATION;'                 => $WorkflowExecDir,
    ';TMP_DIR;'                      => '/usr/local/scratch/legacy2chado',
    ';REPOSITORY_ROOT;'              => '/usr/local/annotation',
    ';WORKFLOW_DIR;'                 => $WorkflowExecDir,
    ';SET_RUNTIME;'                  => $WorkflowExecDir ."/set_runtime_vars.pl",
    ';LEGACY2CHADO_MASTER_CONF;'     => $WorkflowConfigFileDir . "/legacy2chado-master.ini", 
    ';LEGACY2CHADO_CONF;'            => $WorkflowConfigFileDir . "/legacy2chado.ini", 
    ';LEGACY2CHADO_MASTER_TEMPLATE;' => $WorkflowConfigFileDir . "/legacy2chado-master_template.xml",
    ';LEGACY2CHADO_TEMPLATE;'        => $WorkflowConfigFileDir . "/legacy2chado_template.xml",
    ';RUN_WF;'                       => $WorkflowExecDir . "/run_wf.sh",
    ';SERVER;'                       => 'SYBTIGR',
    ';TARGET_SERVER;'                => 'SYBTIGR',
    ';WFNAME;'                       => 'legacy2chado',
    ';USERNAME;'                     => '',
    ';PASSWORD;'                     => '',
    ';TARGET_DATABASE;'              => '',
    ';ORGANISM_FILE;'                => '',
    ';ORGANISM_TYPE;'                => '',
};


my %options = ();
my ($organism_file, $organism_list, $target_database, $config, $man, $verbose, $username, $password, $server, $target_server, $organism_type, $printconf, $organism_format, $help);
my $results = GetOptions (
			  \%options,
			  'username|U=s'            => \$username,
			  'password|P=s'            => \$password,
			  'database|D=s'            => \$target_database,
			  'server|S=s'              => \$server,
			  'server_target|t=s'       => \$target_server,
			  'organism_file|f=s'       => \$organism_file,
			  'organism_list|o=s'       => \$organism_list,
			  'organism_type|b=s'       => \$organism_type,
			  'config|c=s'              => \$config,
			  'help|h'                  => \$help,
			  'man|m'                   => \$man,
			  'verbose|v'               => \$verbose,
			  'printconf|p'             => \$printconf,
			  'organism_format|q'       => \$organism_format,
			  );

&print_configuration($conf) if ($printconf);
&pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if ($man);
&print_usage() if ($help);
&sample_organism_file() if ($organism_format);

#my $log4perl = "/usr/local/annotation/" . $conf->{';TARGET_DATABASE_UC;'} . "/Workflow/" . $conf->{';WFNAME;'} . "/" . $conf->{';WFID;'} . "/legacy2chado_wf.log";
my $log4perl = "legacy2chado_wf.log";

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
if (!defined($server)){
    $server = $conf->{';SERVER;'} if ((exists $conf->{';SERVER;'}) and (defined($conf->{';SERVER;'})));
}
if (!defined($target_server)){
    $target_server = $conf->{';TARGET_SERVER;'} if ((exists $conf->{';TARGET_SERVER;'}) and (defined($conf->{';TARGET_SERVER;'})));
}
if (!defined($organism_file)){
    $organism_file = $conf->{';ORGANISM_FILE;'} if ((exists $conf->{';ORGANISM_FILE;'}) and (defined($conf->{';ORGANISM_FILE;'})));
}
if (!defined($organism_list)){
    $organism_list = $conf->{';ORGANISM_LIST;'} if ((exists $conf->{';ORGANISM_LIST;'}) and (defined($conf->{';ORGANISM_LIST;'})));
}
if (!defined($organism_type)){
    $organism_type = $conf->{';ORGANISM_TYPE;'} if ((exists $conf->{';ORGANISM_TYPE;'}) and (defined($conf->{';ORGANISM_TYPE;'})));
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
# Now determine whether the correct combination of input parameters have been specified/provided by the user OR user's configuration file OR default configuration
#

#
# User must specify the username, password and database
#
print STDERR "\nYou must specify username\n\n" if (!defined($username));
print STDERR "\nYou must specify password\n\n" if (!defined($password));
print STDERR "\nYou must specify database\n\n" if (!defined($target_database));
&print_usage() if (!$username or !$password or !$target_database);

#
# User must either specify the ((organism_list) OR (organism_file AND organism_type))
#
if (!$organism_list and !$organism_file){
    print STDERR "\nYou must either specify the organism_file\nOR\nthe organism_list AND the organism_type\nwhere the organism_type types are:\neuk\nprok\nntprok\n\nFor an explanation on how to format the organism_file invoke this script as:\n$0 -q\n\n";
    &print_usage();
}

#
# User must specify the organism_type if they don't specify the organism_file
#
if (!$organism_type and !$organism_file){
    print STDERR "\nYou must specify the type of organism_type\n1=euk\n2=prok\n3=ntprok\n\n";
    &print_usage();
}

#
# In case user wishes to specify an organism_file but accidentally uses -o instead of -f
#
if (($organism_list) and (-e $organism_list)){
    print STDERR "\nUse \"-o\" to specify a comma separated list of organism names.\nUse \"-f\" to specify a filename which contains a list of organisms.\n\n";
    &print_usage();
}

#
# Verify that valid organism_type was specified
#
if ((defined($organism_type)) and ($organism_type !~ /[euk|prok|ntprok]/) and (!defined($organism_file))){
    print STDERR ("\norganism_type:$organism_type invalid.  Must be either euk, prok or ntprok\n\n");
    &print_usage();
}


if ((!defined($organism_type)) and (!defined($organism_file))){
    print STDERR ("\nYou must specify organism_type on the command-line or in the organism_file\n\n");
    &print_usage();
}




#
# Verify whether servers are either SYBIL or SYBTIGR and nothing else
#
if (!defined($server)){
    print STDERR "server was not defined";
    &print_usage();
}
if (!defined($target_server)){
    print STDERR "target_server was not defined";
    &print_usage();
}
if ($server !~ /^SYBIL|SYBTIGR$/){
    print STDERR "Invalid Sybase server: $server, must be either SYBIL or SYBTIGR";
    &print_usage();
}
if ($target_server !~ /^SYBIL|SYBTIGR$/){
    print STDERR "Invalid Sybase target server: $target_server, must be either SYBIL or SYBTIGR";
    &print_usage();
}

#
# Determine the appropriate migration script
#
my $legacy2chado_program = &determine_legacy2chado_program($organism_type) if ($organism_type);


#---------------------------------------------------------------------------------------------
# Set some configuration values
#
#---------------------------------------------------------------------------------------------
$conf->{';WFID;'}                 = $$;
$conf->{';TARGET_DATABASE_UC;'}   = uc($target_database)  if (defined($target_database));
$conf->{';TARGET_DATABASE_LC;'}   = lc($target_database)  if (defined($target_database));
$conf->{';USERNAME;'}             = $username             if (defined($username));
$conf->{';PASSWORD;'}             = $password             if (defined($password));
$conf->{';LEGACY2CHADO_PROGRAM;'} = $legacy2chado_program if (defined($legacy2chado_program));
$conf->{';TARGET_SERVER;'}        = $target_server        if (defined($target_server));
$conf->{';SERVER;'}               = $server               if (defined($server));
$conf->{';ORGANISM_FILE;'}        = $organism_file        if (defined($organism_file));
$conf->{';ORGANISM_LIST;'}        = $organism_list        if (defined($organism_list));
$conf->{';ORGANISM_TYPE;'}        = $organism_type        if (defined($organism_type));
$conf->{';ORGANISM_LIST;'}        = $organism_list        if (defined($organism_list));
$conf->{';WORKFLOW_DIR;'}         = $workflow_dir         if (defined($workflow_dir));
$conf->{';SET_RUNTIME;'}          = $set_runtime          if (defined($set_runtime));
$conf->{';RUN_WF;'}               = $run_wf               if (defined($run_wf));





#---------------------------------------------------------------------------------------------
# Create the workflow archive directory
# repository e.g. mkdir -m 777 -p /usr/local/annotation/CHADO_TEST/Workflow/legacy2chado/3065
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
# Build the bsml_document_list
#
#---------------------------------------------------------------------------------------------
my $organism_hash;

#
# If organism_file was specified, load unique organism/database names into the organism_hash
#
if (defined($organism_file)){

    $organism_hash = &get_organism_hash(\$organism_file);
    $logger->logdie("organism_hash was not defined") if (!defined($organism_hash));
}
#
# If a organism_list was specified, parse and load only unique organism/database names into the organism_hash
# This assumes all asmbl_ids
#
elsif ($organism_list){
    my @organism_array = split(/,/,$organism_list);
    foreach my $org (@organism_array){
	$logger->logdie("org was not defined") if (!defined($org));

	if (!exists $organism_hash->{$org}){
	    $organism_hash->{$org}->{'org_name'}  = $org;
	    $organism_hash->{$org}->{'org_array'} = 'ALL';
	}
    }
}
    
    
#------------------------------------------------------------------------------------------
# For the master workflow, process each organism stored in the organism_hash
#
#------------------------------------------------------------------------------------------
foreach my $org_ref (sort keys %$organism_hash){
    
    $logger->logdie("org_ref was not defined") if (!defined($org_ref));
    
    #
    # Extract the organism database name from the organism_hash
    #
    my $organism_database;
    if (exists $organism_hash->{$org_ref}->{'org_name'}){
	$organism_database = $organism_hash->{$org_ref}->{'org_name'};
    }
    else{
	$logger->logdie("org_ref: $org_ref is not valid");
    }
    
    $conf->{';CONFIG_LIST;'}   .= $workflow_instance_dir ."/". $organism_database ."-instance.ini, ";
    $conf->{';INSTANCE_LIST;'} .= $workflow_instance_dir ."/". $organism_database .".xml, ";
}

#
# Get rid of trailing space and comma
#
chop $conf->{';CONFIG_LIST;'};
chop $conf->{';CONFIG_LIST;'};
chop $conf->{';INSTANCE_LIST;'};
chop $conf->{';INSTANCE_LIST;'};

$conf->{';LOCALIZATION_LOGFILE;'} = $workflow_instance_dir . "/protein_localization.log";
$conf->{';PROTEIN_LOAD_LOGFILE;'} = $workflow_instance_dir . "/protein.load.log";
$conf->{';PROTEIN_LOAD_STATS;'}   = $workflow_instance_dir . "/protein.load.stats";
$conf->{';WORKFLOW_INSTANCE_DIR;'} = $workflow_instance_dir;

#
# Store the master configuration for the legacy2chado workflow instances to follow
# and set global permissions
#

my $master_wf_fullpath = $workflow_instance_dir ."/master.subs";

&store_conf_to_file(\$master_wf_fullpath, $conf);

$execution_string = "chmod 666 $master_wf_fullpath";
&execute(\$execution_string);

#
# Perform the substitution on the master legacy2chado workflow
# and set global permissions
#
$execution_string = $set_runtime ." -c ". $master_wf_fullpath ." < ". $conf->{';LEGACY2CHADO_MASTER_CONF;'} ." > ". $workflow_instance_dir ."/legacy2chado-master-instance.ini";
&execute(\$execution_string);

$execution_string = "chmod 666 " . $workflow_instance_dir ."/legacy2chado-master-instance.ini";
&execute(\$execution_string);

print STDERR "\n$0 execution complete\nPlease review log4perl log file:$log4perl\n\n";

#------------------------------------------------------------------------------------------
# For the organism/database stored in the organism_hash - generate a subflow
#
#------------------------------------------------------------------------------------------
foreach my $org_ref (sort keys %$organism_hash){
	
    $logger->logdie("org_ref was not defined") if (!defined($org_ref));
    
    #
    # Extract the organism database name from the organism_hash
    #
    my $organism_database;
    if (exists $organism_hash->{$org_ref}->{'org_name'}){
	$organism_database = $organism_hash->{$org_ref}->{'org_name'};
    }
    else{
	$logger->logdie("org_ref: $org_ref is not valid");
    }
    
    
    $conf->{';DATABASE_UC;'} = uc ($organism_database);
    $conf->{';DATABASE_LC;'} = lc ($organism_database);
    
    $logger->info("Processing organism:$organism_database");


    $execution_string = "mkdir -p -m 777 " . $conf->{';TMP_DIR;'};
    &execute(\$execution_string);
    $execution_string = "mkdir -p -m 777 " . $conf->{';TMP_DIR;'} ."/". $conf->{';TARGET_DATABASE_LC;'};
    &execute(\$execution_string);
    $execution_string = "mkdir -p -m 777 " . $conf->{';TMP_DIR;'} ."/". $conf->{';TARGET_DATABASE_LC;'} . "/" .$organism_database;
    &execute(\$execution_string);

    #
    # Extract the legacy2chado_program
    #
    if (exists $organism_hash->{$org_ref}->{'org_type'}){
	$legacy2chado_program = &determine_legacy2chado_program($organism_hash->{$org_ref}->{'org_type'});
	$logger->info("legacy2chado_program for organism:$organism_database is set to:$legacy2chado_program");
    }


    #
    # Extract the array of asmbl_ids if exists in organism_hash
    #
    my $organism_array;
    if (exists $organism_hash->{$org_ref}->{'org_array'}){
	$organism_array = $organism_hash->{$org_ref}->{'org_array'};
    }
    
    my $asmbl_id_list;
    if ($organism_array ne 'ALL'){
	$asmbl_id_list = join(',', @$organism_array);
	$logger->info("assemblies:$asmbl_id_list");
    }
    else{
	$asmbl_id_list = "ALL";
    }
    
    #
    # Create a subdirectory in the workflow instance directory for each organism/database
    # We will store all the associated goodies in these organism/database specific directories
    #
    my $execution_string = "mkdir -p -m 777 ". $workflow_instance_dir ."/". $conf->{';DATABASE_LC;'};
    &execute(\$execution_string);
    
    #
    # reminder- workflow_instance_dir is like: /usr/local/annotation/CHADO_TEST/Workflow/legacy2chado/6409
    #
    $conf->{';ASMBL_ID_LIST;'}        = $asmbl_id_list;
    $conf->{';LEGACY2CHADO_LOGFILE;'} = $workflow_instance_dir ."/". $conf->{';DATABASE_LC;'} ."/". $conf->{';DATABASE_LC;'} .".mig.log";
    $conf->{';LOAD_LOGFILE;'}         = $workflow_instance_dir ."/". $conf->{';DATABASE_LC;'} ."/". $conf->{';DATABASE_LC;'} .".load.log";
    $conf->{';LOAD_STATS;'}           = $workflow_instance_dir ."/". $conf->{';DATABASE_LC;'} ."/". $conf->{';DATABASE_LC;'} .".load.stats";
    $conf->{';DBSPACE_LOGFILE;'}      = $workflow_instance_dir ."/". $conf->{';DATABASE_LC;'} ."/". $conf->{';DATABASE_LC;'} .".dbspace.log";
    $conf->{';VALIDATION_LOGFILE;'}   = $workflow_instance_dir ."/". $conf->{';DATABASE_LC;'} ."/". $conf->{';DATABASE_LC;'} .".validation.log";
    $conf->{';VALIDATION_OUTDIR;'}    = $workflow_instance_dir ."/". $conf->{';DATABASE_LC;'};
    $conf->{';LOCALIZATION_LOGFILE;'} = $workflow_instance_dir ."/". $conf->{';DATABASE_LC;'} ."/". $conf->{';DATABASE_LC;'} .".protein.local.opt.log";
    $conf->{';CONFIRM_LOAD_LOGFILE;'} = $workflow_instance_dir ."/". $conf->{';DATABASE_LC;'} ."/". $conf->{';DATABASE_LC;'} .".confirm.load.log";
    $conf->{';CONFIRM_LOAD_OUTDIR;'}  = $workflow_instance_dir ."/". $conf->{';DATABASE_LC;'};


    #
    # Store the configuration for this specific organism/database
    # and set global permissions
    #
    my $organism_subflow = $workflow_instance_dir ."/". $conf->{';DATABASE_LC;'} ."/". $organism_database .".sub";
    &store_conf_to_file(\$organism_subflow, $conf);
    
    $execution_string = "chmod 666 $organism_subflow";
    &execute(\$execution_string);
    

    #
    # Perform sed replacement on the specific organism/database subflow configuration .ini file 
    #
#    $execution_string = $set_runtime ." -c ". $organism_subflow ." < ". $conf->{';LEGACY2CHADO_CONF;'} ." > ". $workflow_instance_dir ."/". $conf->{';DATABASE_LC;'} ."/legacy2chado_instance.ini";
    $execution_string = $set_runtime ." -c ". $organism_subflow ." < ". $conf->{';LEGACY2CHADO_CONF;'} ." > ". $workflow_instance_dir ."/". $organism_database ."-instance.ini";
    &execute(\$execution_string);

    #
    # Give the subflow instance configuration file global permissions
    #
    $execution_string = "chmod 666 " . $workflow_instance_dir ."/". $organism_database ."-instance.ini";
    &execute(\$execution_string);


}


#
# Copy the legacy2chado_template.xml        to legacy2chado-instance_template.xml AND
#          legacy2chado-master_template.xml to legacy2chado-master-instance_template.xml
#
# and then set global permissions
#
$execution_string = "cp ". $conf->{';LEGACY2CHADO_TEMPLATE;'} ." ". $workflow_instance_dir ."/legacy2chado-instance_template.xml";
&execute(\$execution_string);
$execution_string = "cp ". $conf->{';LEGACY2CHADO_MASTER_TEMPLATE;'} ." ". $workflow_instance_dir ."/legacy2chado-master-instance_template.xml";
&execute(\$execution_string);

$execution_string = "chmod 666 ". $workflow_instance_dir . "/legacy2chado-instance_template.xml";
&execute(\$execution_string);
$execution_string = "chmod 666 ". $workflow_instance_dir . "/legacy2chado-master-instance_template.xml";
&execute(\$execution_string);
    

#
# Run the legacy2chado workflow
#
$execution_string = $run_wf ." -d ". $workflow_instance_dir ." -c ". $workflow_instance_dir ."/legacy2chado-master-instance.ini" . " -t ". $workflow_instance_dir ."/legacy2chado-master-instance_template.xml" ." -i ". $workflow_instance_dir ."/legacy2chado.xml" . " -l ". $workflow_instance_dir ."/log.txt";
&execute(\$execution_string);


print STDERR "\n$0 execution complete\nPlease review log4perl log file:$log4perl\n\n";




#----------------------------------------------------------------------------------------------------------------------------------------
#
#                                          END OF MAIN -- SUBROUTINES FOLLOW
#
#----------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------
# get_file_contents()
#
#-------------------------------------------------------------------
sub get_file_contents {

    my $file = shift;
    $logger->logdie("file was not defined") if (!defined($file));

    if (&is_file_status_ok($file)){

	open (IN_FILE, "<$$file") || $logger->logdie("Could not open file: $$file for input");
	my @contents = <IN_FILE>;
	chomp @contents;
	
	return \@contents;

    }
    else{
	exit;#$logger->logdie("file $$file does not have appropriate permissions");
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

#-----------------------------------------------------------------
# get_directory_contents()
#
#-----------------------------------------------------------------
sub get_directory_contents{

    my $dir = shift;
    $logger->logdie("dir was not defined") if (!defined($$dir));
    $logger->logdie("dir:$dir is not a directory") if (!-d $$dir);
	
    my @bsml_files = qx"ls -1 $$dir/*.bsml";
    chomp @bsml_files;

    return \@bsml_files;

}#end sub get_directory_contents()

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
    $logger->logdie("conf was not defined") if (!defined($conf));

    foreach my $key (sort keys %$conf){
	$logger->logdie("key was not defined") if (!defined($key));
	print  ("$key=$conf->{$key}\n");
    }

    exit;

}#end sub print_configuration()


#-----------------------------------------------------------------
# get_organism_hash()
#
#-----------------------------------------------------------------
sub get_organism_hash {

    my $file = shift;

    $logger->debug("Entered get_organism_hash") if $logger->is_debug();

    $logger->logdie("file was not defined") if (!defined($file));

    my $contents = &get_file_contents($file);
    $logger->logdie("contents was not defined") if (!defined($file));
	
    my $org_obj = {};
    my $organism;
    my $asmbl_id;
    my $org_ctr=0;
    my $type;


    foreach my $inline (@$contents){
	$logger->logdie("inline was not defined") if (!defined($inline));

	
	if (($inline =~ /^organism/) or ($inline =~ /^type/) or ($inline =~ /^asmbl_id/)){

	    #
	    # parse organism data
	    #
	    if ($inline =~ /^organism:(\S+)\s*/){
		$organism = $1;
		$logger->logdie("organism was not defined") if (!defined($organism));
		$org_obj->{$organism}->{'org_name'} = $organism;
		$org_ctr++;
	    }

	    #
	    # parse type data
	    #
	    if ($inline =~ /^type:(\S+)\s*$/){
		$type = $1;
		$logger->logdie("type was not defined") if (!defined($type));
		
		if ($type =~ /[euk|prok|ntprok]/){
		    $org_obj->{$organism}->{'org_type'} = $type;
		    # good
		}
		else{
		    $logger->logdie("invalid type:$type");
		}
		
	    }

	    #
	    # parse asmbl_id data
	    #
	    if ($inline =~ /^asmbl_id:(\S+)\s*$/){
		my $string = $1;
		my $array = &value_array($string);
		$logger->logdie("array was not defined") if (!defined($array));

		if ($array eq 'ALL'){
		    $org_obj->{$organism}->{'org_array'} = "ALL";
		}
		else{
		    foreach my $lk (sort @$array){
			$logger->logdie("lk was not defined") if (!defined($lk));
			push (@{$org_obj->{$organism}->{'org_array'}}, $lk);
		    }
		}
	    }
	
	}
	else{
	    $logger->fatal("The organism_file was poorly formatted");
	    &sample_organism_file();
	}

    }

    return $org_obj;

}#end sub get_organism_hash()

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

#-------------------------------------------------------------------
# determine_legacy2chado_program()
#
#-------------------------------------------------------------------
sub determine_legacy2chado_program {

    my $script = shift;

    $logger->logdie("script was not defined") if (!defined($script));

    my $program;

#    if ($script !~ /^[euk|prok|ntprok]$/){
#	$logger->logdie("You must specify one of three organism types: euk or prok or ntprok\nThis type: $script is not acceptable\n");
#    }
    
    if ($script eq 'euk'){	
	$program = "euktigr2chado";
    }
    elsif ($script eq 'prok'){
	$program = "proktigr2chado";
    }
    elsif ($script eq 'ntprok'){
	$program = "ntproktigr2chado";
    }
    else{
	$logger->logdie("Unrecognized organism type: $script");
    }


    $conf->{';LEGACY2CHADO_PROGRAM;'} = $program;
    return $program;
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


#---------------------------------------------------------------------------
# value_array
#
#---------------------------------------------------------------------------
sub value_array{
    my $string = shift;
    $logger->logdie("string was not defined") if (!defined($string));
    
    #
    # string could be:
    # 
    # a) "1"
    # b) "1,2,3"
    # c) "1-9"
    # d) "1,2,5-9,11,13,21-25"
    # e) "ALL"
    #
    # Note that if "ALL" is found anywhere in the string, all
    # other values will be ignored and "ALL" will be returned.


    #
    # Condition e) string = "ALL"
    #
    if ($string =~ /ALL/){
	$logger->info("\"ALL\" was detected, all other specified values will be ignored.  Returning \"ALL\"");
	return "ALL";
    }

    my $unique = {};

    #
    # Condition a) e.g. string = "1"
    #

    if ($string =~ /^\s*(\d+)\s*$/){
	my $val = $1;
	$logger->info("Only one value was detected.");
	$unique->{$val} = $val if (!exists $unique->{$val});
    }

    #
    # Condition b) e.g. string = "1,2,3"
    #
    if (($string =~ /\,/) and ($string !~ /\-/)){
	my @values = split(/,/,$string);
	foreach my $val (@values){
	    $logger->logdie("val was not defined") if (!defined($val));
	    $unique->{$val} = $val if ((!exists $unique->{$val}) and ($val =~/^\d+$/));
	}
    }

    #
    # Condition c) e.g. string = "1-9"
    #
    if ($string =~ /^\s*(\d+\-\d+)\s*$/){
    
	my $range = $1;
	$unique = &resolve_range($unique, $range);

    }

    #
    # Condition d) e.g. string = "1,2,5-9,11,13,21-25"
    #
    if (($string =~ /\,/) and ($string =~ /\-/)){
	
	my @commaseparated = split(/,/,$string);
	foreach my $val (@commaseparated){
	    $logger->logdie("val was not defined") if (!defined($val));
	    
	    if ((!exists $unique->{$val}) and ($val =~/^\d+$/)){
		$unique->{$val} = $val;
	    }
	    elsif ($val =~ /^\s*(\d+\-\d+)\s*$/){
		my $range = $1;
		$unique = &resolve_range($unique, $range);

	    }	

	}
    }


    #
    # All parsing of string complete, now process unique hash and
    # return reference to array containing unique values
    #
    my @array;
    foreach my $key (sort keys %$unique){
	$logger->logdie("key was not defined") if (!defined($key));

	push (@array, $key);
    }

    $logger->debug("Here is a complete list of unique values:\n@array") if $logger->is_debug();
    return \@array;
}


#-------------------------------------------------------
# resolve_range()
#
#-------------------------------------------------------
sub resolve_range {
    
    my ($unique, $range) = @_;

    $logger->logdie("unique was not defined") if (!defined($unique));
    $logger->logdie("range was not defined")  if (!defined($range));


    if ($range =~ /^\s*(\d+)\-(\d+)\s*$/){

	my $start = $1;
	my $stop = $2;
	foreach my $val ($start..$stop){
	    $unique->{$val} = $val if (!exists $unique->{$val});
	}
    }
    else{
	$logger->logdie("Could not parse $range");
    }

    return $unique;



}#end sub resolve_range()





#-----------------------------------------------------------------
# sample_organism_file()
#
#-----------------------------------------------------------------
sub sample_organism_file {

    print STDERR "The organism file should contain three types of lines.\n";
    print STDERR "1)\"organism\" line contains the name of the organism's legacy database name e.g. gbs or lma2 or ntsp01\n";
    print STDERR "2)\"type\" line contains the type of legacy organism database e.g. prok or euk or ntprok\n";
    print STDERR "3)\"asmbl_id\" line is EITHER \"ALL\" signifying that all assemblies are to be processed,\n";
    print STDERR "OR a positive integer thus indicating the actual asmbl_id to be processed\n";
    print STDERR "Example:\n";
    print STDERR "-------------------------------------------------------------\n";
    print "organism:gbs\n";
    print "type:prok\n";
    print "asmbl_id:1,2,3-9,11,13,21-33\n";
    print "organism:gbs18rs21\n";
    print "type:prok\n";
    print "asmbl_id:1,2,3\n";
    print "asmbl_id:8\n";
    print "asmbl_id:21\n";
    print "organism:lma2\n";
    print "type:euk\n";
    print "asmbl_id:3-12,16\n";
    print "asmbl_id:18\n";
    print "organism:gma\n";
    print "type:prok\n";
    print "asmbl_id:ALL\n";
    print STDERR "-------------------------------------------------------------\n";
    print STDERR "\n";
    print STDERR "\n";
    print STDERR "\n";
    exit;

}#end sub sample_organism_file()

#-----------------------------------------------------------------
# print_usage()
#
#-----------------------------------------------------------------
sub print_usage {

    print STDERR "SAMPLE USAGE:  $0  -U username -P password -D database [-b organism_type] [-f organism_file|-o organism_list] [-S server] [-c config_file] [-h] [-m] [-p] [-q] [-t target_server] [-v]\n";
    print STDERR "  -U|--username            = Database login username\n";
    print STDERR "  -P|--password            = Database login password\n";
    print STDERR "  -S|--server              = Server\n";
    print STDERR "  -D|--database            = Target database to migrate legacy data into\n";
    print STDERR "  -b|--organism_type       = euk OR prok OR ntprok\n";
    print STDERR "  -f|--organism_file       = file containing list of organisms and asmbl_ids to be processed\n";
    print STDERR "  -o|--organism_list       = comma separated list of similar type (euk, prok or ntprok) organisms to be migrated\n";
    print STDERR "  -c|--config_file         = configuration file\n";
    print STDERR "  -h|--help                = print out this help info\n";
    print STDERR "  -m|--man                 = print out pod2usage man page for this utility\n";
    print STDERR "  -p|--printconf           = print out the key=value configuration pairs\n";
    print STDERR "  -q|--organism_format     = print out formatting rules for the organism_file\n";
    print STDERR "  -t|--target_server       = target server\n";
    print STDERR "  -v|--verbose             = Increases log4perl reporting level to screen from WARN to INFO\n";
    exit 1;

}#end sub print_usage()

