#!/usr/local/bin/perl
#---------------------------------------------------------------------
# script name: run_chado2bsml.pl
# date:        2003.10.23
#
#
#
#
#---------------------------------------------------------------------

=head1 NAME

run_chado2bsml.pl - Creates and executes the chado to bsml migration workflow

=head1 SYNOPSIS

USAGE:  run_chado2bsml.pl -U username -P password -D database [-S server] [-f asmbl_file|-l asmbl_list] [-c config] [-h] [-m] [-p] [-q] [-v]

=head1 OPTIONS

=over 8

=item B<--username,-U>
    
    Database username

=item B<--password,-P>
    
    Database password

=item B<--database,-D>
    
    Source database name

=item B<--server,-S>
    
    Source server

=item B<--asmbl_file,-f>
    
    File containing list of chado assembly identifiers to be processed.  Note that the list should be newline separated.  Use the -q switch to view a sample asmbl_file.

=item B<--asmbl_list,-l>

    Command-line comma separated list of chado assembly identifiers to be processed


=item B<--config,-c>

    Configuration file containing all substitution keys parameters

=item B<--help,-h>

    Prints out usage

=item B<--man,-m>

    Prints out this pod2usage man page

=item B<--printconf,-p>

    Prints out the default configuration key-value pairs

=item B<--bsml_format,-q>

    Prints out formatting rules for the asmbl_file

=item B<--verbose,-v>

    Increases log4perl reporting level to screen from WARN to INFO


=back

=head1 DESCRIPTION

    run_chado2bsml.pl - Creates and executes the bsml to chado migration workflow

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
# The default chado2bsml configuration
#
my $conf = {
    ';INSTALLATION;'        => $WorkflowExecDir,
    ';TMP_DIR;'             => '/usr/local/scratch/chado2bsml',
    ';REPOSITORY_ROOT;'     => '/usr/local/annotation',
    ';WORKFLOW_DIR;'        => $WorkflowExecDir,
    ';SET_RUNTIME;'         => $WorkflowExecDir ."/set_runtime_vars.pl",
    ';CHADO2BSML_CONF;'     => $WorkflowConfigFileDir ."/chado2bsml.ini", 
    ';CHADO2BSML_TEMPLATE;' => $WorkflowConfigFileDir ."/chado2bsml_template.xml",
    ';RUN_WF;'              => $WorkflowExecDir ."/run_wf.sh",
    ';WFNAME;'              => 'chado2bsml',
    ';USERNAME;'            => '',
    ';PASSWORD;'            => '',
    ';DATABASE;'            => '',
    ';SERVER;'              => 'SYBTIGR',
    ';ASMBL_LIST;'          => '',
    ';ASMBL_FILE;'          => '',
    };



my %options = ();
my ($username, $password, $database, $server, $asmbl_file, $asmbl_list, $asmbl_format, $config, $help, $man, $printconf, $verbose);
my $results = GetOptions (
			  \%options,
			  'username|U=s'        => \$username,
			  'password|P=s'        => \$password,
			  'database|D=s'        => \$database,
			  'server|S=s'          => \$server,
			  'asmbl_file|f=s'      => \$asmbl_file,
			  'asmbl_list|l=s'      => \$asmbl_list,
			  'asmbl_format|q'      => \$asmbl_format,
			  'config|c=s'          => \$config,
			  'help|h'              => \$help,
			  'man|m'               => \$man,
			  'verbose|v'           => \$verbose,
			  'printconf|p'         => \$printconf,
			  );

&print_configuration($conf) if ($printconf);
&pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if ($man);
&print_usage() if ($help);
&sample_asmbl_file() if ($asmbl_format);

#
# Increase log4perl reporting level to screen
#
my $screen_threshold = 'ERROR';
$screen_threshold = 'INFO' if ($verbose);

my $log4perl = "chado2bsml_wf.log";
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
# Load the conf2 hash with override values defaults if the user has specified a config file
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
if (!defined($asmbl_file)){
    $asmbl_file = $conf->{';ASMBL_FILE;'} if ((exists $conf->{';ASMBL_FILE;'}) and (defined($conf->{';ASMBL_FILE;'})));
}
if (!defined($asmbl_list)){
    $asmbl_list = $conf->{';ASMBL_LIST;'} if ((exists $conf->{';ASMBL_LIST;'}) and (defined($conf->{';ASMBL_LIST;'})));
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
# Verify whether the correct combination of parameters have been specified by user on command-line OR configuration file OR default configuration
#


print STDERR "\nYou must specify username\n\n" if (!defined($username));
print STDERR "\nYou must specify password\n\n" if (!defined($password));
print STDERR "\nYou must specify database\n\n" if (!defined($database));

&print_usage() if (!$username or !$password or !$database);

$logger->info("All assemblies in database: $database will be processed") if (!$asmbl_file and !$asmbl_list);


if ($server !~ /^SYBIL|SYBTIGR$/){
    print STDERR "Invalid Sybase server: $server, must be either SYBIL or SYBTIGR";
    &print_usage();
}
if (!defined($server)){
    print STDERR "You must specify the server";
    &print_usage();
}
    

#
# Store the final configuration in conf
#
$conf->{';DATABASE_UC;'} = uc($database)   if (defined($database));
$conf->{';DATABASE_LC;'} = lc($database)   if (defined($database));
$conf->{';USERNAME;'}    = $username       if (defined($username));
$conf->{';PASSWORD;'}    = $password       if (defined($password));
$conf->{';ASMBL_FILE;'}  = $asmbl_file     if (defined($asmbl_file));
$conf->{';ASMBL_LIST;'}  = $asmbl_list     if (defined($asmbl_list));
$conf->{';WFID;'}        = $$;



#---------------------------------------------------------------------------------------------
# Create the workflow archive directory within the BSML 
# repository e.g. mkdir -m 777 -p /usr/local/annotation/CHADO_TEST/Workflow/chado2bsml/3065
#
#---------------------------------------------------------------------------------------------
$logger->logdie("$conf->{';REPOSITORY_ROOT;'} is not a directory") if (!-d $conf->{';REPOSITORY_ROOT;'});

my $execution_string;

$execution_string = "mkdir -m 777 -p ". $conf->{';REPOSITORY_ROOT;'} ."/". $conf->{';DATABASE_UC;'} ."/Workflow";
&execute(\$execution_string);

$execution_string = "mkdir -m 777 -p ". $conf->{';REPOSITORY_ROOT;'} ."/". $conf->{';DATABASE_UC;'} ."/Workflow/" .$conf->{';WFNAME;'};
&execute(\$execution_string);

$execution_string = "mkdir -m 777 -p ". $conf->{';REPOSITORY_ROOT;'} ."/". $conf->{';DATABASE_UC;'} ."/Workflow/" .$conf->{';WFNAME;'} ."/". $conf->{';WFID;'};
&execute(\$execution_string);

my $workflow_instance_dir = $conf->{';REPOSITORY_ROOT;'} ."/". $conf->{';DATABASE_UC;'} ."/Workflow/" .$conf->{';WFNAME;'} ."/". $conf->{';WFID;'};


my $db2bsml_logfile = $workflow_instance_dir . "/db2bsml.log";
$conf->{';DB2BSML_LOGFILE;'} = $db2bsml_logfile;


my $asmbl_id_full_list;
if (($asmbl_file eq "") and ($asmbl_list eq "")){
    $asmbl_id_full_list = "ALL";
}

elsif ((defined($asmbl_file)) or (defined($asmbl_list))){
    my $asmbl_id_hash = &retrieve_asmbl_id_hash($asmbl_file, $asmbl_list);
    $logger->logdie("asmbl_id_hash was not defined") if (!defined($asmbl_id_hash));

    #
    # Build the asmbl_id list for db2bsml and generate_genomic_seq if asmbl_file was specified
    #

    foreach my $asmbl (keys %$asmbl_id_hash){
	
	if (!defined($asmbl)){
	    $logger->logdie("asmbl was not defined");
	}
	
	my $asm = $asmbl_id_hash->{$asmbl}->{'seq'} if (exists $asmbl_id_hash->{$asmbl}->{'seq'});
	$logger->logdie("asm was not defined") if (!defined($asm));
	
	$asmbl_id_full_list .= "$asm,";

    }
    chop $asmbl_id_full_list;
}



$conf->{';ASMBL_ID;'} = $asmbl_id_full_list;

$conf->{';WORKFLOW_INSTANCE_DIR;'} = $workflow_instance_dir;

#
# Store the configuration for this chado2bsml workflow instance
#
my $master_wf_fullpath = $workflow_instance_dir ."/master.subs";
&store_conf_to_file(\$master_wf_fullpath, $conf);


#
# Perform the sed replacements on the chado2bsml-instance.ini
#
    
$execution_string = $set_runtime ." -c ". $workflow_instance_dir ."/master.subs < ". $conf->{';CHADO2BSML_CONF;'} ." > " . $workflow_instance_dir ."/chado2bsml-instance.ini";
&execute(\$execution_string);

#
# Copy the chado2bsml_template.xml to chado2bsml-instance_template.xml 
#
$execution_string = "cp ". $conf->{';CHADO2BSML_TEMPLATE;'} ." ". $workflow_instance_dir ."/chado2bsml-instance_template.xml";
&execute(\$execution_string);

#
# Set global permissions for the:
# a) chado2bsml-instance.ini
# b) master key-id substitution file
# c) chado2bsml-instance_template.xml
#
$execution_string = "chmod 666 ". $workflow_instance_dir ."/chado2bsml-instance.ini";
&execute(\$execution_string);

$execution_string = "chmod 666 ". $workflow_instance_dir ."/master.subs";
&execute(\$execution_string);

$execution_string = "chmod 666 ". $workflow_instance_dir ."/chado2bsml-instance_template.xml";
&execute(\$execution_string);

#
# Create and execute the chado2bsml workflow
#
$execution_string = $run_wf ." -d ". $workflow_instance_dir ." -c ". $workflow_instance_dir ."/chado2bsml-instance.ini -t ". $workflow_instance_dir ."/chado2bsml-instance_template.xml -i ". $workflow_instance_dir ."/chado2bsml.xml -l ". $workflow_instance_dir ."/log.txt";
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
    if (!defined($$dir)){
	$logger->logdie("dir was not defined");
    }
    if (!-d $$dir){
	$logger->logdie("dir:$dir is not a directory");
    }

#    open (BSMLDIR, "$$dir") or $logger->logdie("Could not open directory:$$dir");
#    my @bsml_files =  readdir BSMLDIR;

#   my @bsml_files = system"ls -1 $$dir/*.bsml";
   my @bsml_files = qx"ls -1 $$dir/*.bsml";
   chomp @bsml_files;
#    foreach my $file (@bsml_files){
#	print STDERR ("$file\n");
#    }
#    print Dumper @bsml_files;
#    die;
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
# print_configuration()
#
#-----------------------------------------------------------------
sub print_configuration {

    my $conf = shift;
    $logger->logdie("conf was not defined") if (!defined($conf));

    foreach my $key (sort keys %$conf){
	$logger->logdie("key was not defined") if (!defined($key));
	print ("$key=$conf->{$key}\n");
    }

    exit;
}#end sub print_configuration()


#---------------------------------------------------------------------------
# retrieve_asmbl_id_hash()
#
#---------------------------------------------------------------------------
sub retrieve_asmbl_id_hash {

    my ($file, $list) = @_;


    my $asmbl_hash;
    #
    # If bsml_file was specified, retrieve all of the bsml documents listed in the bsml file
    # and add unique bsml document names to the bsml_document_hash
    #
    if (defined($file)){

	my $contents = &get_file_contents($file);
	$logger->logdie("contents was not defined") if (!defined($contents));
	
	foreach my $asmbl (@$contents){
	    $logger->logdie("asmbl was not defined") if (!defined($asmbl));

	    $asmbl_hash->{$asmbl}->{'seq'} = $asmbl if (!exists $asmbl_hash->{$asmbl}->{'seq'});

	}
    }
    
    #
    # If a bsml_directory was specified, retrieve all BSML documents contained in the
    # specified directory and add unique bsml document names tot he bsml_document_hash
    #
    if ($list){

	my @array = split(/,/,$list);
	foreach my $asmbl (@array){
	    $logger->logdie("asmbl was not defined") if (!defined($asmbl));
	    $asmbl_hash->{$asmbl}->{'seq'} = $asmbl if (!exists $asmbl_hash->{$asmbl}->{'seq'});
	}
    }

    return $asmbl_hash;


}#end sub retrieve_asmbl_id_hash()


#-----------------------------------------------------------------
# sample_asmbl_file()
#
#-----------------------------------------------------------------
sub sample_asmbl_file {

    print STDERR "\n";
    print STDERR "Each line of the asmbl_file should contain a single chado assembly identifier.\n";
    print STDERR "This identifier should exist in the chado database\n";
    print STDERR "This is an example of how the asmbl_file should be formatted:\n";
    print STDERR "-----------------------------------------------------------------------------------------------\n";
    print STDERR "gbs_799_assembly\n";
    print STDERR "bsp_3839_assembly\n";
    print STDERR "-----------------------------------------------------------------------------------------------\n";
    print STDERR "\n";
    print STDERR "\n";
    exit;
}		      

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


#-----------------------------------------------------------------
# print_usage()
#
#-----------------------------------------------------------------
sub print_usage {

    print STDERR "SAMPLE USAGE:  $0  -U username -P password -D database [-S server] [-f asmbl_file|-l asmbl_list] [-c config_file] [-h] [-m] [-p] [-q] [-v]\n";
    print STDERR "  -U|--username      = database login username\n";
    print STDERR "  -P|--password      = database login password\n";
    print STDERR "  -D|--database      = name of database to extract features from\n";
    print STDERR "  -S|--server        = source server name\n";
    print STDERR "  -f|--asmbl_file    = file containing list of chado asmbl_ids\n";
    print STDERR "  -l|--asmbl_list    = command-line comma separate list of chado asmbl_ids\n";
    print STDERR "  -c|--config_file   = configuration file\n";
    print STDERR "  -h|--help          = print usage\n";
    print STDERR "  -m|--man           = display the pod2usage man page for this utility\n";
    print STDERR "  -p|--printconf     = print out the default key=value configuration pairs\n";
    print STDERR "  -q|--asmbl_format  = print out the formatting rules for the asmbl_file\n";
    print STDERR "  -v|--verbose       = Increases log4perl reporting level to screen from WARN to INFO\n";
    exit 1;

}

