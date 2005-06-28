#!/usr/local/bin/perl

=head1  NAME

run_pipeline.pl

=head1 SYNOPSIS

Run a workflow pipeline from a pipeline configuration file

USAGE: run_pipeline.pl -conf config file [--printconf] [--listconf] [--skiprun] [--wfid id] [--debug level] [--log file]

=head1 OPTIONS

B<--config, -c> Configuration file in INI format.

B<--printconf, -p> Print configuration file.

B<--listconf> List available configuration file.  This option will
list files named *conf.ini in the WORKFLOW_DOCS_DIR installation
directory, which is specified at installation time

B<--wfid id> Override setting for wfid with id.

B<--skiprun> Skip execution of workflow

B<--log, -l> [OPTIONAL] log file

B<--help, -h> [OPTIONAL]  program help

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Workflow::Logger;
use Workflow::Builder;
use Workflow::IteratorBuilder;
use Config::IniFiles;

umask(0000);

# Installation directories are pulled from environment for now
my $WorkflowDocsDir = $ENV{'WORKFLOW_DOCS_DIR'};
my $BinDir = $ENV{'WORKFLOW_WRAPPERS_DIR'};
my $SchemaDocsDir = $ENV{'SCHEMA_DOCS_DIR'};

my %options = ();

my $results = GetOptions (\%options, 
                          'conf|c=s', 
                          'printconf|p', 
			  'listconf',
			  'log|l=s',
			  'ini|i=s',
			  'sectionid|s=s',
			  'wfid|w=s',
			  'pipelineid=s',
                          'debug=s', 
			  'skiprun',
			  'help|h' );

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();

my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

if(!(-e $WorkflowDocsDir)){
    $logger->get_logger()->logdie("Can't read directory $WorkflowDocsDir.  Set environment variable WORKFLOW_DOCS_DIR");
}

if( $options{'listconf'} ){
    &listconf($WorkflowDocsDir);
    exit;
}

if( $options{'conf'} eq ''){
    $logger->get_logger()->logdie("Need to input a workflow configuration file with -conf option.");
}

if(!(-e $options{'conf'})){
    $logger->get_logger()->logdie("Can't read configuration file $options{'conf'}");
}

my $origcfg = new Config::IniFiles( -file => $options{'conf'});

if( $options{'printconf'} ){ 
    my $allcfg = &import_includes($origcfg);
    if($WorkflowDocsDir ne ''){
	$logger->get_logger()->debug('Setting $;WORKFLOWDOCS_DIR$; to '.$WorkflowDocsDir) if($logger->get_logger->is_debug());
	$allcfg->setval('init','$;WORKFLOWDOCS_DIR$;',$WorkflowDocsDir);
    }
    if($BinDir ne ''){
	$logger->get_logger()->debug('Setting $;BIN_DIR$; to '.$BinDir) if($logger->get_logger->is_debug());
	$allcfg->setval('init','$;BIN_DIR$;',$BinDir);
    }
    if($SchemaDocsDir ne ''){
	$logger->get_logger()->debug('Setting $;WORKFLOWDOCS_DIR$; to '.$SchemaDocsDir) if($logger->get_logger->is_debug());
	$allcfg->setval('init','$;SCHEMA_DIR$;',$SchemaDocsDir);
    }
    &printconf($allcfg);
    exit;
}

#import init section if not found
if(!( defined $origcfg->val('init', '$;WFID$;'))){
    $origcfg = &import_includes($origcfg);
}
 
if($origcfg->val('init', '$;WFID$;') =~ /\S/){
    $logger->get_logger()->warn('$;WFID$; set in configuration file to '.$origcfg->val('init', '$;WFID$;'));
}
elsif($options{'wfid'}){
    $logger->get_logger()->debug('Setting $;WFID$; to '.$options{'wfid'}) if($logger->get_logger->is_debug());
    $origcfg->setval('init','$;WFID$;',$options{'wfid'});
}
else{
    $logger->get_logger()->debug('Setting $;WFID$; to '.$$) if($logger->get_logger->is_debug());
    $origcfg->setval('init','$;WFID$;',$$);
}

if($origcfg->val('init', '$;PIPELINEID$;') =~ /\S/){
    $logger->get_logger()->warn('$;PIPELINEID$; set in configuration file to '.$origcfg->val('init', '$;PIPELINEID$;'));
}
elsif($options{'pipelineid'}){
    $logger->get_logger()->debug('Setting $;PIPELINEID$; to '.$options{'pipelineid'}) if($logger->get_logger->is_debug());
    $origcfg->setval('init','$;PIPELINEID$;',$options{'pipelineid'});
}
else{
    $logger->get_logger()->debug('Setting $;PIPELINEID$; to 0') if($logger->get_logger->is_debug());
    $origcfg->setval('init','$;PIPELINEID$;',0);
}



my $cfg = &replace_keys($origcfg);
&check_parameters($cfg);


my @workflows = $cfg->GroupMembers("workflowdocs");


foreach my $workflow (@workflows){
    my $workflowname = $cfg->val($workflow,'$;NAME$;');
    $logger->get_logger()->debug("Create working directory ".$cfg->val( $workflow, '$;WORKFLOW_REPOSITORY$;')) if($logger->get_logger()->is_debug());
    my $wfrepositorydir = $cfg->val( $workflow, '$;WORKFLOW_REPOSITORY$;');
    &create_directory($wfrepositorydir);

    my $instanceconfigfile = "$wfrepositorydir/pipeline.config";
    my $instanceinifile = "$wfrepositorydir/".$workflowname."_pipeline.ini";    
    my $instancexmlfile = "$wfrepositorydir/"."pipeline.xml";
        
    if($options{'ini'} ne "" ){
	#write ini file
	if($options{'sectionid'} ne ""){
	    $logger->get_logger()->debug("Writing ini entry for id=$options{'sectionid'} ini file=$options{'ini'}\n");
	    my $inicfg;
	    if(-e $options{'ini'}){
		$inicfg = new Config::IniFiles(-file=>$options{'ini'});
	    }
	    else{
		$inicfg = new Config::IniFiles();
	    }
	    $inicfg->newval($options{'sectionid'},"param.command","$BinDir/run_wf");
	    $inicfg->newval($options{'sectionid'},'param.--instance',$instancexmlfile);
	    $inicfg->WriteConfig($options{'ini'});
	}
	else{
	    $logger->get_logger()->warn("INI file $options{'sectionid'} specified for writing with no section id");
	}
    }
    my $finalcfg = &replace_keys($cfg);

    #write out config file to repository
    if($cfg->val($workflow, '$;COMPONENT_CONFIG$;') =~ /\S/){
	$logger->get_logger()->warn('$;COMPONENT_CONFIG$; set in configuration file to '.$cfg->val($workflow, '$;COMPONENT_CONFIG$;'));
	$instanceconfigfile = $cfg->val($workflow, '$;COMPONENT_CONFIG$;');
    }
    else{
	$logger->get_logger()->debug('Setting $;COMPONENT_CONFIG$; to '.$instanceconfigfile) if($logger->get_logger->is_debug());
	$instanceconfigfile = "$wfrepositorydir/pipeline.config";
    }

    $cfg->setval($workflow,'$;COMPONENT_CONFIG$;',$instanceconfigfile);
    $finalcfg->WriteConfig($instanceconfigfile);

    
    $logger->get_logger()->debug("Wrote configuration file $instanceconfigfile");

    #create builder
    my $wfmasterobj = new Workflow::Builder('NAME'=>$workflowname); 
    $wfmasterobj->set_config($finalcfg);

    #generate instance ini file
    my $inisuccess = $wfmasterobj->generate_instance_ini($finalcfg->val( $workflow, '$;MASTER_TEMPLATE_INI$;'),$instanceinifile);
    if($inisuccess){
	$logger->get_logger()->debug("Created master workflow ini $instanceinifile")  if($logger->get_logger()->is_debug());
    }
    else{
	$logger->get_logger()->error("Unable to create master workflow ini $instanceinifile. Status $inisuccess")  if($logger->get_logger()->is_debug());
    }

    #generate instance
    my $xmlsuccess = $wfmasterobj->generate_instance_xml($instanceinifile,$finalcfg->val( $workflow, '$;MASTER_TEMPLATE_XML$;'), $instancexmlfile);
    if($xmlsuccess){
	$logger->get_logger()->debug("Created master workflow instance $instancexmlfile")  if($logger->get_logger()->is_debug());
    }
    else{
	$logger->get_logger()->error("Unable to create master workflow instance $instancexmlfile")  if($logger->get_logger()->is_debug());
    }

    #invoke master instance
    if(!($options{'skiprun'})){
	my $wfexec = new Workflow::Run();
	my $runlogfile = "$instancexmlfile.run.log";
	my $runoutfile = "$instancexmlfile.run.out";
	$logger->get_logger()->debug("Running instance file $instancexmlfile") if($logger->get_logger()->is_debug());
	print "Running instance file $instancexmlfile\n";
	print "To launch viewer execute... $BinDir/run_wfmonitor.sh -i $instancexmlfile\n";
	my $runstatus = $wfexec->RunWorkflow($instancexmlfile,$runlogfile,$runoutfile);
	if($runstatus == 0){
	    $logger->get_logger()->debug("Execution complete of $instancexmlfile with status $runstatus") if($logger->get_logger()->is_debug());
	    print "Execution complete\n";
	    exit(0);
	}
	else{
	    $logger->get_logger()->logdie("Unable to run instance file $instancexmlfile. Return code $runstatus");
	    exit ($runstatus);
	}
    }
    else{
	$logger->get_logger()->info("Skipping execution of workflow $instancexmlfile");
	print "Skipping execution of workflow $instancexmlfile\n";
    }

}


sub printconf{
    my $cfg = shift;

    $cfg->WriteConfig("/tmp/$$.config");
    open FILE, ("/tmp/$$.config");
    while (my $line=<FILE>){
	print $line;
    }
    close FILE;
}


sub check_parameters{
    my($cfg) = @_;
    foreach my $param ($cfg->Parameters('init')){
	if($cfg->val( 'init', $param) eq ""){
	    $logger->get_logger()->logdie("Required parameter $param is missing (".$cfg->val( 'init', $param).").  Check init section of configuration file.\n\n");
	}
	else{
	    $logger->get_logger()->debug("Required parameter $param=".$cfg->val( 'init', $param)) if($logger->get_logger()->is_debug());
	}
    }
}

sub import_includes{
    my($cfg) = @_;
    my @includes = $cfg->GroupMembers("include");
    my $currcfg = $cfg;
    foreach my $member (@includes){
	my @parameters = $cfg->Parameters ($member);
	$logger->get_logger()->debug("Scanning parameters in $member") if($logger->get_logger()->is_debug());
	foreach my $param (@parameters){
	    my $includefile = $cfg->val($member,$param);
	    $logger->get_logger()->debug("Found includefile $includefile in section [ $member ] with key $param") if($logger->get_logger()->is_debug());
	    $includefile =~ s/\$;WORKFLOWDOCS_DIR\$;/$WorkflowDocsDir/g;
	    $logger->get_logger()->debug("Set includefile=$includefile using \$;WORKFLOWDOCS_DIR\$;=$WorkflowDocsDir") if($logger->get_logger()->is_debug());
	    my $newcfg = new Config::IniFiles( -file => $includefile, 
					       -import => $currcfg);
	    $currcfg = $newcfg;
	}
    }
    return $currcfg;
}

sub replace_keys{
    my($cfg) = @_;
    my $allkeys = {};
    my $checkvalues = {};
    my @sections = $cfg->Sections();
    my $currcfg = $cfg;
    foreach my $section (@sections){
	my @parameters = $cfg->Parameters ($section);
	foreach my $param (@parameters){
	    my $value = $cfg->val($section,$param);
	    
        ## take spaces off front and back of value, leaving internal ones
	    $value =~ s/^\s*(.+?)\s*$/$1/;
	    
        $cfg->setval($section,$param,$value);
	    $allkeys->{$param}->{'value'} = $value;
	    $allkeys->{$param}->{'section'} = $section;
	    $logger->get_logger()->debug("Scanning $value for key $param in section [ $section ] as candidate for replacement") if($logger->get_logger()->is_debug());
	    if($value =~ /\$;[\w_]+\$;/){
		$logger->get_logger()->debug("Found $value for key $param in section [ $section ] as candidate for replacement") if($logger->get_logger()->is_debug());
		$checkvalues->{$param} = $value;
	    }
	}
    }

    foreach my $key (keys %$checkvalues){
	my $value = $checkvalues->{$key};
	$logger->get_logger()->debug("Replacing $key with $value") if($logger->get_logger()->is_debug());
	while(&checkvalue($value,$allkeys)){
	    $value =~ s/(\$;[\w_]+\$;)/&replaceval($1,$allkeys)/ge;
	    $logger->get_logger()->debug("Value redefined as $value") if($logger->get_logger()->is_debug());
	}
	my $setval = $cfg->setval($allkeys->{$key}->{'section'},$key,$value);
	$logger->get_logger()->debug("Replaced $key in section [ $allkeys->{$key}->{'section'} ] with $value. Return val $setval") if($logger->get_logger()->is_debug());
	if(!$setval){
	    $logger->get_logger()->logdie("Key $key in section [ $allkeys->{$key}->{'section'} ] is not valid");
	}
    }
    return $cfg;
}


sub replaceval{
    my ($val,$keylookup) = @_;
    if(!(exists $keylookup->{$val})){
	$logger->get_logger()->logdie("Bad key $1 in configuration file") if($logger->get_logger()->is_debug());
    }
    elsif($keylookup->{$val}->{'value'} eq ''){
	return $val;
    }
    else{
	return $keylookup->{$val}->{'value'};
    }
}

sub checkvalue{
    my ($val,$keylookup) = @_;
    if($val =~ /\$;[\w_]+\$;/){
	my($lookupval) = ($val =~ /(\$;[\w_]+\$;)/); 
	if($keylookup->{$lookupval}->{'value'} eq ""){
	    return 0;
	}
	else{
	    return 1;
	}
    }
}

sub listconf{
    my ($dir) = @_;
    $logger->get_logger()->debug("Reading directory $dir for conf files") if($logger->get_logger()->is_debug());
    opendir (DIR,$dir) or $logger->get_logger()->logdir("Can't open directory $dir : $!");
    my @conffiles = grep { /conf\.ini$/ }  readdir(DIR);
    foreach my $conf (@conffiles){
	if($conf ne 'sharedconf.ini'){
	    print "$dir/$conf\n";
	}
    }
}

sub create_directory{
    my($dir) = @_;
    my $ret = system("mkdir -m 777 -p $dir");
    $ret >>= 8;
    return $ret;
}
