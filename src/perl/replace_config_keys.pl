#!/usr/bin/perl

#Example command line 
#perl -I ../../lib/ replace_config_keys.pl --template_conf=/usr/local/devel/ANNOTATION/angiuoli/code/ergatis/workflow/wait/wait.default.user.config --output_conf=/tmp/test.config --keys=PIPELINEID=23 --debug=5 --log=my.log

use lib (@INC,$ENV{"PERL_MOD_DIR"});
use strict;
use Config::IniFiles;
use Ergatis::Logger;
use Getopt::Long;

my %options;

my $delimeter = '$;';
my $delimeterregex = '\$\;[\w_]+\$\;';
#Keys must match $delimeter[\w_]+$delimeter

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
#Import included files
my $included = {};
$origcfg = &import_includes($origcfg,$included);

#Add add'l keys specified via --keys
&add_keys($origcfg,'component',split(/,/,$options{'keys'}));

#Perform initial key replacement
my $cfg = &replace_keys($origcfg);

#Check that the parameters that need values have them
&check_parameters($cfg,'project');

#Write the output location of this file as key $;COMPONENT_CONFIG$;
my $ret = $cfg->setval('component',$delimeter.'COMPONENT_CONFIG'.$delimeter,$options{'output_conf'});
if(!$ret){
    $logger->logdie("Couldn't add key ${delimeter}COMPONENT_CONFIG$delimeter=$options{'output_conf'}to section [component]");
}
&check_parameters($cfg,'component');

#Perform second key replacement for nested keys
my $finalcfg = &replace_keys($cfg);
$finalcfg->WriteConfig($options{'output_conf'});

$logger->debug("Wrote configuration file $options{'output_conf'}");

exit;

sub check_parameters{
    my ($cfg,$section) = @_;

    my %optional = ( $delimeter . 'PROJECT_CODE' . $delimeter => 1,
                     $delimeter . 'SKIP_WF_COMMAND' . $delimeter => 1);

    for my $param ( $cfg->Parameters($section) ) {
        ## skip checking if this one's optional
        next if defined $optional{$param};
    
        if( $cfg->val( $section, $param) eq "" ){
            $logger->logdie("Required parameter $param is missing (".$cfg->val( $section, $param).").  Check [$section] section of configuration file.\n\n");
        } else{
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
            
            my $ret = $cfg->setval($section,$param,$value);
            if(!$ret){
                $logger->logdie("Couldn't add key $param=$value to section [$section]");
            }
            $allkeys->{$param}->{'value'} = $value;
            $allkeys->{$param}->{'section'} = $section;
            $logger->debug("Scanning $value for key $param in section [ $section ] as candidate for replacement") if($logger->is_debug());
            if($value =~ /$delimeterregex/){
                $logger->debug("Found $value for key $param in section [ $section ] as candidate for replacement") if($logger->is_debug());
                $checkvalues->{$param} = $value;
            }
        }
    }

    use Data::Dumper;
    print "These are the keys with values that need replacing:\n";
    print Dumper( $checkvalues );    

    #Now that we've obtained all possible values, expand...
    foreach my $key (keys %$checkvalues){
        my $value = $checkvalues->{$key};
        $logger->debug("Replacing $key with $value") if($logger->is_debug());
        while(&checkvalue($value,$allkeys)){
            print "\nValue before: $value\n";
            $value =~ s/($delimeterregex)/&replaceval($1,$allkeys)/ge;
            print "Value after: $value\n";
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
    return $keylookup->{$val}->{'value'};
}

sub checkvalue{
    my ($val,$keylookup) = @_;
    my $retval = 0;
    if($val =~ /$delimeterregex/){
        my($lookupval) = ($val =~ /($delimeterregex)/); 
        $retval = 1 if( exists( $keylookup->{$lookupval} ) );
    }
    $retval;
}

#import_includes() - Support for inline import of ini formatted
#configuration files

#eg.
#[include]
#SOMEKEY=file.ini

#The entire contents of file.ini will be imported into the current
#configuration file (SOMKEY is ignored).  Keys with the same name will
#take the value provided in the included file, not the original config
#file.  If there are multiple included ini files, the values in the
#last file have highest precedence.  The file named "software.config"
#is handled special: only values from section [common] or
#[$componentname] are parsed.

sub import_includes{
    my($cfg,$included) = @_;
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
	    if(-e $includefile){
		if(!$included->{$includefile}){
		    if($includefile =~ /software.config$/){
			my $softcfg = new Config::IniFiles( -file => $includefile);
			my $newcfg = new Config::IniFiles(  -import => $currcfg);
			my $componentname = $cfg->val("component",$delimeter.'COMPONENT_NAME'.$delimeter);
			$componentname =~ s/\s//g;
			foreach my $section ($softcfg->Sections()){
			    if($section =~ /^common/){
				$newcfg->AddSection($section);
				foreach my $p ($softcfg->Parameters($section)){
				    $newcfg->newval($section,$p,$softcfg->val($section,$p));
				}
			    }
			    elsif($section =~ /component $componentname$/){
				$newcfg->AddSection($section);
				foreach my $p ($softcfg->Parameters($section)){
				    $newcfg->newval($section,$p,$softcfg->val($section,$p));
				}
			    }
			}
			$included->{$includefile} = 1;
			#Supports nesting of includes
			$currcfg = &import_includes($newcfg,$included);
		    }
		    else{
			my $newcfg = new Config::IniFiles( -file => $includefile, 
							   -import => $currcfg);
			#If keys are not properly delimeted, add delimeter here
			foreach my $section ($newcfg->Sections()){
			    foreach my $p ($newcfg->Parameters($section)){
				if($p =~ /$delimeterregex/){
				}
				else{
				    $logger->logdie("Parameter already contains delimeter") if($p  =~ /^\$\;/);
				    $logger->debug("Updating parameter $p to $delimeter.$p.$delimeter") if($logger->is_debug());
				    $newcfg->newval($section,$delimeter.$p.$delimeter,$newcfg->val($section,$p));
				    $newcfg->delval($section,$p);
				}
			    }
			}
			$included->{$includefile} = 1;
			#Supports nesting of includes
			$currcfg = &import_includes($newcfg,$included);
		    }
		}
	    }
	    else{
		$logger->logdie("Can't find included config file in $member $param $includefile");
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
	my $ret = $cfg->setval($section,$delimeter.$key.$delimeter,$value);
	if(!$ret){
	    $logger->logdie("Couldn't add key $delimeter $key $delimeter=$value to section [$section]");
	}
    }
}
