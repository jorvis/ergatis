#!/usr/local/bin/perl

use strict;

use XML::Twig;

use CGI qw(:standard);

use Tree::DAG_Node;
use Config::IniFiles;
use Data::Dumper;

umask(0000);

my $conffile = param('conffile');
my $dowrite = param('dowrite');
my $limitsect = param('limitsect');
my $ignoresect = param('ignoresect');
my $sharedconf = param('sharedconf');
my $pipeline_id = param('pipeline_id');
my $component_name = param('component_name');
my $repository_root = param('repository_root');

if($dowrite){
    my $newcfg = new Config::IniFiles(-file => $conffile);
    print header();
    my $filewritten = &write_config_from_params($newcfg);
    
    if ($filewritten) {
        print "<html><body>Output file $filewritten written.<br><pre>\n";
        &print_config($newcfg);
        print "</pre><br><a href='javascript:window.close();'>[close]</a></body></html>";
    } else {
        print "<html><body>There was an error creating the output file</body></html>\n";
    }


    exit;
}

my $cfg = new Config::IniFiles(-file => $conffile);

## we have to hide some configurable sections from conf files that have already
##  been written.  These include all ending in .bld.ini
my $preexisting = 0;
if ($conffile =~ /\.bld\.ini$/ ) {
    $preexisting = 1;
}

#print all conf options

my $ignoresectlookup = {};
foreach my $elt (split (/,/,$ignoresect)){
    $ignoresectlookup->{$elt} = 1;
}
my $limitsectlookup = {};
foreach my $elt (split (/,/,$limitsect)){
    $limitsectlookup->{$elt} = 1;
}

my $currcfg = &import_includes($cfg,$sharedconf);

print header();

print <<PageHeadeR;
<html>

<head>

<style type="text/css">
    body {
	    font-family: verdana,sans-serif;
	    font-size: 10pt;
	    line-height: 1.5em;
    }
    
     td, input {
	    font-family: verdana,sans-serif;
	    font-size: 10pt;
     }
     
     input {
        border: 1px solid rgb(150,150,150);
        padding-left: 5px;
     }
     
     h3 {
        color: rgb(100, 100, 100);
        text-decoration: underline;
        margin-bottom: 5px;
        font-size: 11pt;
     }
     
     table {
        margin-left: 5px;
     }
     
     tr.param_comment td {
        font-size: 8pt;
        color: rgb(100, 150, 100);
        padding-top: 3px;

     }
</style>

</head>

<body>

<form name='myform'><input type=hidden name='dowrite' value=1><input type=hidden name='conffile' value='$conffile'>
PageHeadeR

#print "we have access to this: <pre> " . Dumper(\$currcfg) . "</pre>";

$repository_root = $currcfg->val('init', '$;REPOSITORY_ROOT$;') || die "couln't get a REPOSITORY ROOT";

if($pipeline_id ne ""){
    $currcfg->setval('init','$;PIPELINE$;',"${component_name}_$pipeline_id");
}

$currcfg->setval('init','$;PIPELINEID$;', $pipeline_id);

my @sections = $currcfg->Sections();

foreach my $section (@sections){
    if(($limitsect && $limitsectlookup->{$section}) || !$limitsect){
	if($ignoresectlookup->{$section}){
	    print "<h3>$section</h3>\n<table>\n";
	}
	else{
	    print "<h3>$section</h3>\n<table>\n";
	}
	my @parameters = $currcfg->Parameters ($section);
	foreach my $param (@parameters){
	    my $value = $currcfg->val($section,$param);
        my @param_comments = $currcfg->GetParameterComment( $section, $param );
        
        ## print out the comments here
        for my $comment ( @param_comments ) {
            ## only do the comment if it starts with ;;  (double semi-colon)
            if ($comment =~ /^\s*\;\;\s*(\S.*)/ ) {
                print "\t<tr class='param_comment'><td>&nbsp;</td><td>$1</td>\n";
            }
        }
        
        ## trim whitespace off the ends, but leave internal (keeps opts like "-f 10")
        if ($value =~ /^\s*(.+?)\s*$/) {
	        $value = $1;
        }
        
	    if($ignoresectlookup->{$section}){
		print "\t<tr><td>$param</td><td>$value</td></tr>\n";
	    } else { 
            ## don't display editable OUTPUT_TOKEN if this is a pre-existing component build
            if ($preexisting && $param eq '$;OUTPUT_TOKEN$;') {
            print "\t<tr><td>$param</td><td>$value</td></tr>\n";
            } else {
		    print "\t<tr><td>$param</td><td><input type=text name='edit_",$section,"::","$param' value='$value' size=100></td></tr>\n";
            }
	    }
	}
    }
    print "</table>\n";
}

print "<input type=hidden name=repository_root value='$repository_root'>\n";
print "<input type=hidden name=component_name value='$component_name'>\n";
print "<input type=hidden name=pipeline_id value='$pipeline_id'>\n";
print "<input type='submit' value='Write config'>\n</form>\n</body>\n</html>";

sub import_includes{
    my($cfg,$sharedconf) = @_;
    my @includes = $cfg->GroupMembers("include");
    my $currcfg = $cfg;
    foreach my $member (@includes){
	my @parameters = $cfg->Parameters ($member);
	foreach my $param (@parameters){
	    if($param eq '$;SHARED_CONFIG$;'){
		my $includefile =  $sharedconf;#$cfg->val($member,$param);
		if(-e $includefile){
		    my $newcfg = new Config::IniFiles( -file => $includefile, 
						       -import => $currcfg);
		    $newcfg->setval($member,$param,$includefile);
		    $currcfg = $newcfg;
		}
	    }
	}
    }
    return $currcfg;
}

sub write_config_from_params{
    my($newcfg) = @_;
    my @allparams = param();
    foreach my $p (@allparams){
	if($p =~ /^edit_/){
	    #print "Parsing $p<br>\n";
	    my $pval = $p;
	    $pval =~ s/^edit_//;
	    my($section,$secp) = split(/\:\:/,$pval);
	    #print "Writing [$section] $secp=",param($p),"<br><\n>";
	    my $s = $newcfg->setval($section,$secp,param($p));
	    if(!$s){
		    die "Can't write  [$section] $secp".param($p);
	    }
	}
    }
    
    ## if we haven't got a pipeline_id we need to derive it from the conffile,
    ##  this should happen when editing an existing config
    if (! $pipeline_id) {
        if ( $conffile =~ m|Workflow/.+?/(\d+)_| ) {
            $pipeline_id = $1;
        }
    }
    
    my $output_token = $newcfg->val("output $component_name", '$;OUTPUT_TOKEN$;') || die "no OUTPUT_TOKEN found";
    my $outputdir = "$repository_root/Workflow/$component_name/${pipeline_id}_$output_token";
    
    ## make sure the output directory exists
    &create_directory($outputdir);
    
    my $outputfile = $outputdir . "/component.conf.bld.ini";
    
    $newcfg->WriteConfig($outputfile);
    if (-e $outputfile) {
        `chmod 666 $outputfile`;
        return $outputfile;
    } else {
        return 0;
    }
}

sub print_config{
    my($cfg) = @_;
    my @sections = $cfg->Sections();
    foreach my $section (@sections){
	print "[$section]\n";
	my @parameters = $cfg->Parameters ($section);
	foreach my $param (@parameters){
	    my $value = $cfg->val($section,$param);
	    print "$param=$value\n";
	}
    }
}


sub create_directory{
    my($dir) = @_;
    my $ret = system("mkdir -m 777 -p $dir");
    $ret >>= 8;

    return $ret;
}    
