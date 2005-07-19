#!/usr/local/bin/perl

use strict;

use XML::Twig;

use CGI qw(:standard);
#use CGI::Carp qw(fatalsToBrowser set_message);
use File::Basename;
use File::Copy;
use Tree::DAG_Node;
use Config::IniFiles;
use XMLManip;

## get the structure of the previously run pipeline
##  

my $xmltemplate = param('xmltemplate');
my $node = param('node');
my $pipelinexml = param('pipelinexml');
my $location = param('location'); #first_child,before,after
my $id = param('id');
my $oldpipelineid = &get_pipeline_id($pipelinexml);
my $newpipelineid = &get_pipeline_id($xmltemplate);

my $pipelineini = Config::IniFiles->new( -file => "$xmltemplate.ini" );

#my $newpipelinexml = XMLManip::get_pipeline_xml($pipelineid,$xmltemplate,$pipelinexml);
#################################################
my $newpipelinexml = '';
## this is pretty straightforward - I don't think we need to Twig it
open(my $ifh, "<$pipelinexml") || die "can't read pipeline $pipelinexml: $!";

while (my $line = <$ifh>) {
    ## skip the commandSetRoot
    next if ( $line =~ /<[\/]*commandSetRoot/ );
    
    ## replace the start configMapId with the old pipeline ID
    if ( $line =~ m|<configMapId>(.*?)</configMapId>| ) {
        my $configmapid = $1;
        
        if ($configmapid eq 'start') {
            $configmapid = "pipeline_$oldpipelineid";
            $line =~ s/start/$configmapid/;
        }
    
        ## make sure this configMapId is in the pipeline.xml.ini
        if (! $pipelineini->SectionExists( $configmapid ) ) {
            $pipelineini->AddSection( $configmapid );
        }
    
    ## some arg elements define the component conf
    } elsif ( $line =~ /<arg>-c (.+\/component.conf.bld.ini).+--pipelineid=$oldpipelineid/ ) {
        my $old_conf_loc = $1;
        my $new_conf_loc = $old_conf_loc;
        $new_conf_loc =~ s/$oldpipelineid/$newpipelineid/g;
        
        ## make sure the directory exists where we want to put this file
        create_directory( (fileparse($new_conf_loc, qr{\.ini}))[1] );
        
        ## copy the conf file to its place in this pipeline
        copy($old_conf_loc, $new_conf_loc);
        
        ## point to the new file now
        $line =~ s/$oldpipelineid/$newpipelineid/g;
        
    ## the new pipelineid needs to be reflected in fileName lines    
    } elsif ( $line =~ /<fileName.+pipeline.xml/ ) {
        $line =~ s/$oldpipelineid/$newpipelineid/g;
    }
    
    $newpipelinexml .= $line;
}

## write back out the pipeline.xml.ini
$pipelineini->RewriteConfig();

#################################################
my $t = new XML::Twig;
   $t->parse($newpipelinexml);
my $pipeline2add = $t->root;

XMLManip::addxml($xmltemplate,$node,$pipeline2add,$location);
print redirect(-uri=>"./show_pipeline.cgi?xmltemplate=$xmltemplate&edit=1");

#print "Content-Type: text/plain\n\n";
#print $newpipelinexml;

sub get_pipeline_id{
    my($file) = @_;
    my($id) = ($file =~ m|Workflow/.+?/(\d+)|);
    if($id eq ""){
	    $id = 0;
    }
    return $id;
}

sub create_directory{
    my($dir) = @_;
    my $ret = system("mkdir -m 777 -p $dir");
    $ret >>= 8;
    return $ret;
}  
