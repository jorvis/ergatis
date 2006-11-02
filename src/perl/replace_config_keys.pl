#!/usr/local/bin/perl

#Example command line 
#perl -I ../../lib/ replace_config_keys.pl --template_conf=/usr/local/devel/ANNOTATION/angiuoli/code/ergatis/workflow/wait/wait.default.user.config --output_conf=/tmp/test.config --keys=PIPELINEID=23 --debug=5 --log=my.log

use lib (@INC,$ENV{"PERL_MOD_DIR"});
use strict;
use Config::IniFiles;
use Ergatis::Logger;
use Getopt::Long;

my %options;
# Installation directories are pulled from environment for now
my $WorkflowDocsDir = $ENV{'WORKFLOW_DOCS_DIR'};
my $BinDir = $ENV{'WORKFLOW_WRAPPERS_DIR'};
my $SchemaDocsDir = $ENV{'SCHEMA_DOCS_DIR'};

my $results = GetOptions (\%options, 
                          'template_conf|t=s', 
                          'output_conf|o=s' ,
                          'keys|k=s' ,
			  'log=s',
                          'debug=s', 
                          'help|h' );


my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				 'LOG_LEVEL'=>$options{'debug'});

$logger = $logger->get_logger();

my $origcfg = new Config::IniFiles( -file => $options{'template_conf'});
#Import Included files
$origcfg = &import_includes($origcfg);
#Add add'l keys specified via --keys

&add_keys($origcfg,'init',split(/,/,$options{'keys'}));
#Perform initial key replacement
my $cfg = &replace_keys($origcfg);

#Init section is composed of required parameters
&check_parameters($cfg,'init');
#Write the output location of this file as key $;COMPONENT_CONFIG$;
$cfg->setval('workflowdocs','$;COMPONENT_CONFIG$;',$options{'output_conf'});
&check_parameters($cfg,'workflowdocs');
#Perform second key replacement for nested keys
my $finalcfg = &replace_keys($cfg);
$finalcfg->WriteConfig($options{'output_conf'});


$logger->debug("Wrote configuration file $options{'output_conf'}");

exit;

sub check_parameters{
    my($cfg,$section) = @_;
    foreach my $param ($cfg->Parameters($section)){
	if($cfg->val( $section, $param) eq ""){
	    $logger->logdie("Required parameter $param is missing (".$cfg->val( $section, $param).").  Check init section of configuration file.\n\n");
	}
	else{
	    $logger->debug("Required parameter $param=".$cfg->val( $section, $param)) if($logger->is_debug());
	}
    }
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
	    $logger->debug("Scanning $value for key $param in section [ $section ] as candidate for replacement") if($logger->is_debug());
	    if($value =~ /\$;[\w_]+\$;/){
		$logger->debug("Found $value for key $param in section [ $section ] as candidate for replacement") if($logger->is_debug());
		$checkvalues->{$param} = $value;
	    }
	}
    }

    foreach my $key (keys %$checkvalues){
	my $value = $checkvalues->{$key};
	$logger->debug("Replacing $key with $value") if($logger->is_debug());
	while(&checkvalue($value,$allkeys)){
	    $value =~ s/(\$;[\w_]+\$;)/&replaceval($1,$allkeys)/ge;
	    $logger->debug("Value redefined as $value") if($logger->is_debug());
	}
	my $setval = $cfg->setval($allkeys->{$key}->{'section'},$key,$value);
	$logger->debug("Replaced $key in section [ $allkeys->{$key}->{'section'} ] with $value. Return val $setval") if($logger->is_debug());
	if(!$setval){
	    $logger->logdie("Key $key in section [ $allkeys->{$key}->{'section'} ] is not valid");
	}
    }
    return $cfg;
}


sub replaceval{
    my ($val,$keylookup) = @_;
    if(!(exists $keylookup->{$val})){
	$logger->logdie("Bad key $1 in configuration file") if($logger->is_debug());
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

sub import_includes{
    my($cfg) = @_;
    my @includes = $cfg->GroupMembers("include");
    push @includes,"include";
    $logger->debug("Looking for [include] sections...found:".@includes) if($logger->is_debug());
    my $currcfg = $cfg;
    foreach my $member (@includes){
	$logger->debug("Scanning parameters in $member") if($logger->is_debug());
	my @parameters = $cfg->Parameters ($member);
	foreach my $param (@parameters){
	    my $includefile = $cfg->val($member,$param);
	    $logger->debug("Found includefile $includefile in section [ $member ] with key $param") if($logger->is_debug());
	    $includefile =~ s/\$;WORKFLOWDOCS_DIR\$;/$WorkflowDocsDir/g;
	    $logger->debug("Set includefile=$includefile using \$;WORKFLOWDOCS_DIR\$;=$WorkflowDocsDir") if($logger->is_debug());
	    if(-e $includefile){
		my $newcfg = new Config::IniFiles( -file => $includefile, 
						   -import => $currcfg);
		$currcfg = $newcfg;
	    }
	    else{
		$logger->logdie("Can't find included config file $includefile");
	    }
	}
    }
    return $currcfg;
}

sub add_keys{
    my($cfg,$section,@keys) = @_;
	$logger->debug("Adding user defined keys: ".@keys) if($logger->is_debug());
    foreach my $kv (@keys){
	my($key,$value) = split(/=/,$kv);
	$logger->debug("Adding user defined key $key=$value in section [$section]") if($logger->is_debug());
	$cfg->setval($section,'$;'.$key.'$;',$value);
    }
}
