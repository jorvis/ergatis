#!/usr/local/bin/perl
#---------------------------------------------------------------------
# script name: run_delete_analysis.pl
# date:        2004.04.20
#
#
#
#
#---------------------------------------------------------------------

=head1 NAME

run_delete_analysis.pl - Creates and executes the bsml to chado migration workflow

=head1 SYNOPSIS

USAGE:  run_delete_analysis.pl -U username -P password -D target_database [-t target_server] -r algorithm -i analysis_id [-c config] [-h] [-m] [-p] [-q]

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

=item B<--algorithm,-r>
    
    Comma separated list of algorithm types to delete from the database

=item B<--analysis_id,-i>

    Comma separated list of analysis_ids to delete from the database


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


=back

=head1 DESCRIPTION

    run_delete_analysis.pl - Deletes analyses from the computational analysis module tables

=cut




use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Data::Dumper;
use Log::Log4perl qw(get_logger);
use File::Basename;
use Pod::Usage;


my $WorkflowConfigFileDir = $ENV{'WORKFLOW_DOCS_DIR'} || ".";
my $WorkflowExecDir       = $ENV{'WORKFLOW_WRAPPERS_DIR'} || ".";

#
# The default delete_analysis configuration
#
my $conf = {
    ';INSTALLATION;'                 => $WorkflowExecDir,
    ';TMP_DIR;'                      => '/usr/local/scratch/delete_analysis',
    ';REPOSITORY_ROOT;'              => '/usr/local/annotation',
    ';WORKFLOW_DIR;'                 => $WorkflowExecDir,
    ';SET_RUNTIME;'                  => $WorkflowExecDir ."/set_runtime_vars.pl",
    ';DELETE_ANALYSIS_MASTER_CONF;'       => $WorkflowConfigFileDir ."/delete_analysis-master.ini", 
    ';DELETE_ANALYSIS_MASTER_TEMPLATE;'   => $WorkflowConfigFileDir ."/delete_analysis-master_template.xml",
    ';RUN_WF;'                       => $WorkflowExecDir ."/run_wf.sh",
    ';WFNAME;'                       => 'delete_analysis',	
    ';USERNAME;'                     => '',
    ';PASSWORD;'                     => '',               
    ';TARGET_DATABASE;'              => '',               
    ';TARGET_SERVER;'                => 'SYBTIGR',
    ';ALGORITHM;'                    => '',
    ';ANALYSIS_ID;'                  => '',
};


my %options = ();
my ($algorithm, $analysis_id, $target_database, $config, $username, $password, $target_server, $printconf, $help, $man, $workflow_monitor, $bsml_format, $log4perl);
my $results = GetOptions (
			  \%options,
			  'username|U=s'            => \$username,
			  'password|P=s'            => \$password,
			  'target_database|D=s'     => \$target_database,
			  'target_server|t=s'       => \$target_server,
			  'algorithm|r=s'           => \$algorithm,
			  'analysis_id|i=s'         => \$analysis_id,
			  'config|c=s'              => \$config,
			  'help|h'                  => \$help,
			  'man|m'                   => \$man,
			  'printconf|p'             => \$printconf,
			  'bsml_format|q'           => \$bsml_format,
			  'log4perl|l'              => \$log4perl,
			  'workflow_monitor|w'      => \$workflow_monitor,
			  );

&print_configuration($conf) if ($printconf);
&pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if ($man);
&print_usage() if ($help);





#
# initialize the logger
#
$log4perl = "/tmp/run_delete_analysis.pl.log" if (!defined($log4perl));;
my $logger = &retrieve_logger($log4perl);





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
if (!defined($algorithm)){
    $algorithm = $conf->{';ALGORITHM;'} if ((exists $conf->{';ALGORITHM;'}) and (defined($conf->{';ALGORITHM;'})));
}
if (!defined($analysis_id)){
    $analysis_id = $conf->{';ANALYSIS_ID;'} if ((exists $conf->{';ANALYSIS_ID;'}) and (defined($conf->{';ANALYSIS_ID;'})));
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
# Either algorithm or analysis_id must be specified
#
if (!$algorithm and !$analysis_id){
       print STDERR "\nYou must EITHER specify the algorithm or analysis_id\n";
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


#
# Set additional configuration values
#
$conf->{';USERNAME;'}           = $username                  if (defined($username));
$conf->{';PASSWORD;'}           = $password                  if (defined($password));
$conf->{';TARGET_SERVER;'}      = $target_server             if (defined($target_server));
$conf->{';TARGET_DATABASE_LC;'} = lc($target_database)       if (defined($target_database));
$conf->{';TARGET_DATABASE_UC;'} = uc($target_database)       if (defined($target_database));
$conf->{';WFID;'}               = $$;


#---------------------------------------------------------------------------------------------
# Create the workflow archive directory within the BSML 
# repository e.g. mkdir -m 777 -p /usr/local/annotation/CHADO_TEST/Workflow/delete_analysis/3065
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

$conf->{';WORKFLOW_INSTANCE_DIR;'} = $workflow_instance_dir;

#
# Store the master configuration for the delete_analysis workflow instances to follow
# and set global permissions
#
my $master_wf_fullpath = $workflow_instance_dir ."/master.subs";
&store_conf_to_file(\$master_wf_fullpath, $conf);

$execution_string = "chmod 666 $master_wf_fullpath";
&execute(\$execution_string);

#
# Perform the substitution on the master delete_analysis workflow
# and set global permissions
#
$execution_string = $set_runtime ." -c ". $master_wf_fullpath ." < ". $conf->{';DELETE_ANALYSIS_MASTER_CONF;'} ." > ". $workflow_instance_dir ."/delete_analysis-master-instance.ini";
&execute(\$execution_string);

$execution_string = "chmod 666 " . $workflow_instance_dir ."/delete_analysis-master-instance.ini";
&execute(\$execution_string);


#
# Copy the delete_analysis-master_template.xml to delete_analysis-master-instance_template.xml
#
$execution_string = "cp ". $conf->{';DELETE_ANALYSIS_MASTER_TEMPLATE;'} ." ". $workflow_instance_dir ."/delete_analysis-master-instance_template.xml";
&execute(\$execution_string);

$execution_string = "chmod 666 ". $workflow_instance_dir . "/delete_analysis-master-instance_template.xml";
&execute(\$execution_string);
    

#
# Run the delete_analysis workflow
#
$execution_string = $run_wf ." -d ". $workflow_instance_dir ." -c ". $workflow_instance_dir ."/delete_analysis-master-instance.ini" . " -t ". $workflow_instance_dir ."/delete_analysis-master-instance_template.xml" . " -i ". $workflow_instance_dir ."/delete_analysis.xml" ." -l ". $workflow_instance_dir ."/log.txt";
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

    my @bsml_files = qx"find $$dir -name '*.bsml'";
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


#-----------------------------------------------------------------
# retrieve_logger()
#
#-----------------------------------------------------------------
sub retrieve_logger {

    my ($log4perl, $verbose) = @_;

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

    return get_logger;
}


#-----------------------------------------------------------------
# print_usage()
#
#-----------------------------------------------------------------
sub print_usage {

    print STDERR "SAMPLE USAGE:  $0  -U username -P password -t target_server -D target_database -r algorithm -i analysis_id [-c config_file] [-d debug_level] [-m] [-p] [-q]\n";
    print STDERR "  -U|--username          = target_database login username\n";
    print STDERR "  -P|--password          = target_database login password\n";
    print STDERR "  -t|--target_server     = target server\n";
    print STDERR "  -D|--target_database   = name of target_database to load analysis into\n";
    print STDERR "  -l|--log4perl          = Optional - Log4perl logfile.  Default is /tmp/run_delete_analysis.pl.log\n";
    print STDERR "  -d|--debug_level       = Optional - Coati::Logger log4perl logging level.  Default is 0\n";
    print STDERR "  -r|--algorithm         = Optional - Algorithm type to delete from the database\n";
    print STDERR "  -i|--analysis_id       = Optional - Analysis identifier to delete from the database\n";
    print STDERR "  -c|--config_file       = Workflow configuration file\n";
    print STDERR "  -p|--printconf         = Optional - print out the key=value configuration pairs\n";
    print STDERR "  -q|--bsml_format       = Optional - print out format rules for the bsml_batch_file\n";
    exit 1;

}#end sub print_usage()

