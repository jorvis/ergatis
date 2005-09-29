#!/usr/local/bin/perl

use strict;

use XML::Twig;

use CGI qw(:standard);
#use CGI::Carp qw(fatalsToBrowser set_message);

use Tree::DAG_Node;
use Config::IniFiles;
use File::Copy;
use XMLManip;
use XML::Twig;

umask 000;

my $xmltemplate = param('xmltemplate');
my $node = param('node');
my $conf = param('conf');
my $location = param('location'); #first_child,before,after
my $name = param('name');
my $pipeline_id = param('pipeline_id');

#print "Content-Type: text/html\n\n";

my $component_conf = &parse_conf($conf);
my $component_name = param('component_name');
my $output_token = $component_conf->val("output $component_name", '$;OUTPUT_TOKEN$;');
my $bin_dir = $component_conf->val('init', '$;BIN_DIR$;');
my $repository_root = $component_conf->val('init', '$;REPOSITORY_ROOT$;');
my $name_token = "$component_name.$output_token";

my $fakexml = <<FAKExml;
<commandSet type='serial'>
    <configMapId>component_$name_token</configMapId>
    <command>
        <type>RunUnixCommand</type>
        <configMapId>generate_component_$name_token</configMapId>
        <name>Generate component $name_token</name>
        <param>
            <key>command</key>
            <value>$bin_dir/run_pipeline</value>
        </param>
        <arg>-c $repository_root/Workflow/$component_name/${pipeline_id}_$output_token/component.conf.bld.ini --skiprun --pipelineid=$pipeline_id</arg>
    </command>
    <commandSet type="serial" version="2.2">
        <name>$name_token</name>
        <maxParallelCmds>0</maxParallelCmds>
        <configMapId>run_component_$name_token</configMapId>
        <fileName>$repository_root/Workflow/$component_name/${pipeline_id}_$output_token/pipeline.xml</fileName>
        <state>incomplete</state>
    </commandSet>
</commandSet>
FAKExml

my $t = new XML::Twig;
   $t->parse($fakexml);
my $component2add = $t->root;

XMLManip::addxml($xmltemplate,$node,$component2add,$location);

## add these configMapIds to the component's pipeline ini file
my $pipeline_ini = Config::IniFiles->new( -file => "$xmltemplate.ini" );
$pipeline_ini->AddSection("component_$name_token");
$pipeline_ini->AddSection("generate_component_$name_token");
$pipeline_ini->AddSection("run_component_$name_token");
$pipeline_ini->RewriteConfig();

## make sure the component conf is written into the proper directory under this pipeline id
## first, make sure the directory exists
if (! -d "$repository_root/Workflow/$component_name/${pipeline_id}_$output_token") {
    if (create_directory("$repository_root/Workflow/$component_name/${pipeline_id}_$output_token")) {
        die "failed to create directory $repository_root/Workflow/$component_name/${pipeline_id}_$output_token : $?";
    }
}

## now copy the conf file there (unless it is self)
if ($conf ne "$repository_root/Workflow/$component_name/${pipeline_id}_$output_token/component.conf.bld.ini") {
    copy($conf, "$repository_root/Workflow/$component_name/${pipeline_id}_$output_token/component.conf.bld.ini") or die "couldn't copy conf file to directory: $!";
}

print redirect(-uri=>"./show_pipeline.cgi?xmltemplate=$xmltemplate&edit=1");

exit;


sub get_pipeline_id {
    my($file) = @_;
    my($id) = ($file =~ m|pipeline/(\d+)/pipeline.xml|);
    if($id eq ""){
        $id = 0;
    }
    return $id;
}

sub parse_conf {
    my($file) = @_;

    my $cfg = new Config::IniFiles(-file => $file);

    my @includes = $cfg->GroupMembers("include");
    my $currcfg = $cfg;
    foreach my $member (@includes){
	    my @parameters = $cfg->Parameters ($member);
	    foreach my $param (@parameters){
	        my $includefile = $cfg->val($member,$param);
	        my $newcfg = new Config::IniFiles( -file => $includefile, 
					                           -import => $currcfg);
	        $currcfg = $newcfg;
	    }
    }

    return $currcfg;
#    return $currcfg->val("init",'$;BIN_DIR$;');
}

sub create_directory{
    my($dir) = @_;
    my $ret = system("mkdir -m 777 -p $dir");
    $ret >>= 8;
    return $ret;
}  

__DATA__

example usage:

http://jorvis-lx:8080/tigr-scripts/ergatis/add_component.cgi?
    xmltemplate=/usr/local/scratch/annotation/TGA1/Workflow/pipeline/8840/pipeline.xml&
    node=start&location=last_child&
    conf=/usr/local/scratch/annotation/TGA1/Workflow/aat_aa/8840_default/component.conf.bld.ini&
    name=component.conf.bld.ini&
    component_name=aat_aa&
    pipeline_id=8840

