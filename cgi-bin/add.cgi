#!/usr/local/bin/perl

use strict;

use XML::Twig;

use CGI qw(:standard);

use Tree::DAG_Node;
use File::Find;
use File::Basename;
use File::stat;
use Config::IniFiles;

my $node = param('node');
my $xmltemplate = param('xmltemplate');
my $location = param('location');


print header();

print "<html><body>";
print "<h2>Edit $xmltemplate node $node</h2>";
print "<a href='show_pipeline.cgi?xmltemplate=$xmltemplate&edit=1'>[view workflow]</a>";
print "<h3>Add commandset</h3>";
print "<a href='add_commandset.cgi?xmltemplate=$xmltemplate&node=$node&location=$location&type=parallel'>[parallel]</a><br>";
print "<a href='add_commandset.cgi?xmltemplate=$xmltemplate&node=$node&location=$location&type=serial'>[serial]</a><br>";
print "<h3>Configured components</h3>";
my ($outputdir) = ($xmltemplate =~ /(.*)\/Workflow/);
my $pipelinedir = "$outputdir/Workflow/pipeline";
$outputdir .= "/workflow_config_files";
my $sharedconf = "$outputdir/sharedconf.ini";
my $WorkflowDocsDir = "/usr/local/devel/ANNOTATION/cas/docs"; #only used to set up shared conf file
my $componentbldconf = &get_component_blds($outputdir);


if( -e $sharedconf){
    print "Shared config <a  target='config' href='view_component.cgi?conffile=$sharedconf'>[view]<a target='config' href='config_component.cgi?conffile=$sharedconf&outputfile=$sharedconf&limitsect=init'>[edit]</a><br>";
}
print "<table>";
foreach my $componentbld (sort {$componentbldconf->{$b}->{'date'} <=> $componentbldconf->{$a}->{'date'}} keys %$componentbldconf){
    print "<tr><td>$componentbldconf->{$componentbld}->{'type'}</td><td>$componentbldconf->{$componentbld}->{'name'}</td><td>$componentbldconf->{$componentbld}->{'user'}</td><td>".localtime($componentbldconf->{$componentbld}->{'date'})."</td><td><a href='add_component.cgi?xmltemplate=$xmltemplate&node=$node&location=$location&conf=$componentbldconf->{$componentbld}->{'file'}&name=$componentbldconf->{$componentbld}->{'name'}'>[add]</a><a  target='config' href='view_component.cgi?conffile=$componentbldconf->{$componentbld}->{'file'}'>[view]</a><a target='config' href='config_component.cgi?conffile=$componentbldconf->{$componentbld}->{'file'}&outputfile=$componentbldconf->{$componentbld}->{'file'}&ignoresect=init&sharedconf=$sharedconf'>[edit]</a><a href='remove_component.cgi?conffile=$componentbldconf->{$componentbld}->{'file'}'>[remove]</a></td></tr>";
}
print "</table>";
print "<h3>Pipelines</h3>";
print "<table>";
my $pipelineconf = &get_pipeline_blds($pipelinedir);
foreach my $pipeline (sort {$pipelineconf->{$b}->{'date'} cmp $pipelineconf->{$a}->{'date'}} keys %$pipelineconf){
    my $subcomponents = $pipelineconf->{$pipeline}->{'components'};
    if(scalar(@$subcomponents)>1){
	print "<tr><td><a href='show_pipeline.cgi?&xmltemplate=$pipeline'>$pipelineconf->{$pipeline}->{'name'}</a></td><td>$pipelineconf->{$pipeline}->{'user'}</td><td>".localtime($pipelineconf->{$pipeline}->{'date'})."</td><td><a href='add_pipeline.cgi?xmltemplate=$xmltemplate&node=$node&location=$location&pipelinexml=$pipeline'>[add]</a></tr>";
	foreach my $component (@$subcomponents){
	    print "<tr><td>&nbsp;</td><td>$component->{'type'}</td><td><a  target='config' href='view_component.cgi?conffile=$component->{'file'}'>$component->{'name'}</a></td><td>$component->{'user'}</td><td>".localtime($component->{'date'})."</td></tr>";
	}
    }
}
print "</table>";
print "<h3>Configure component from template</h3>";


if(! (-e $sharedconf)){
    print "<a target='config' href='config_component.cgi?conffile=$WorkflowDocsDir/sharedconf.ini&outputfile=$sharedconf&limitsect=init>Shared configuration</a><i>This configuration must be set before adding any components</i><br>";
}
else{
    my $workflowdocsdir = &get_workflow_docs($sharedconf);
    print STDERR "Searching $workflowdocsdir parsed from $sharedconf\n";
    my $componentconf = &get_component_conf($workflowdocsdir);
    foreach my $component (sort keys %$componentconf){
	my $outputfile = "$outputdir/$component"."_$$"."conf.bld.ini";
	print "<a target='config' href='config_component.cgi?conffile=$componentconf->{$component}&outputfile=$outputfile&sharedconf=$sharedconf&ignoresect=init&id=$$'>$component</a><br>";
    }
}

print "</body></html>";

sub get_workflow_docs{
    my($file) = @_;
    my $cfg = new Config::IniFiles(-file => $file);
    return $cfg->val("init",'$;WORKFLOWDOCS_DIR$;');
}
sub get_workflow_bin{
    my($file) = @_;
    my $cfg = new Config::IniFiles(-file => $file);
    return $cfg->val("init",'$;BIN_DIR$;');
}

sub get_pipeline_blds{
    my($dir) = @_;
    my $glob = 'pipeline\d+\.xml$';
    my $pipelinefiles = {};
    find(sub {my $file = $File::Find::name;
	      if($file =~ /$glob/){
		  my($date,$user,$name) = &get_pipeline_bld_info($file);
		  my (@components);
		  my ($stale);
		  my $t1 = new XML::Twig( TwigHandlers => { 'command/arg' => 
							    sub {
								my($t,$elt) = @_;
								my $text = $elt->text();
								my($config) = ($text =~ /-c\s+(\S+\.bld\.ini)/);
								print STDERR "Config $config\n";
								if(-e $config){
								    my($type,$date,$user,$name) = &get_component_bld_info($config);
								    push @components,{'type'=>$type,
										      'date'=>$date,
										      'user'=>$user,
										      'name'=>$name,
										      'file'=>$config};
								}
								else{
								    $stale=1;
								}
							    }
							});
		  if($stale!=1){
		      print STDERR "File $file\n";
		      $t1->parsefile($file);			
		      $pipelinefiles->{$file}->{'date'} = $date;
		      $pipelinefiles->{$file}->{'user'} = $user;
		      $pipelinefiles->{$file}->{'name'} = $name;
		      $pipelinefiles->{$file}->{'components'} = \@components;
		  }
	      }
	  },$dir);
    return $pipelinefiles;
}
		  


sub get_component_blds{
    my($dir) = @_;
    my $glob = 'conf\.bld\.ini';
    my $conffiles = {};
    find(sub {my $file = $File::Find::name;
	      if($file =~ /$glob/){
		  my($type,$date,$user,$name) = &get_component_bld_info($file);
		  $conffiles->{$file}->{'date'} = $date;
		  $conffiles->{$file}->{'user'} = $user;
		  $conffiles->{$file}->{'type'} = $type;
		  $conffiles->{$file}->{'file'} = $file;
		  $conffiles->{$file}->{'name'} = $name;

	      }
	  },$dir);
    return $conffiles;
}

sub get_component_conf{
    my($dir) = @_;
    my $glob = 'conf\.ini';
    my $conffiles = {};
    find(sub {my $file = $File::Find::name;
	      if($file =~ /$glob/){
		  my $cfg = new Config::IniFiles(-file => $file);
		  print STDERR "$file\n";
		  if($cfg){
		      my @workflows = $cfg->GroupMembers("workflowdocs");
		      foreach my $workflow (@workflows){
			  my $name = $cfg->val($workflow,'$;NAME$;');
			  $name =~ s/\s//g;
			  if(exists $conffiles->{$name}){
			      print STDERR "Workflow with name $name already stored\n";
			  }
			  else{
			      $conffiles->{$name} = $file;
			  }
		      }
		  }
		  else{
		      print STDERR "Can't parse config file $file\n";
		  }
	      }},
	 $dir);
    return $conffiles;
}


sub get_component_bld_info{
    my($file) = @_;
    my $cfg = new Config::IniFiles(-file => $file);
    my $fullname;
    if($cfg){
	my @workflows = $cfg->GroupMembers("workflowdocs");
	foreach my $workflow (@workflows){
	    my $wfname = $cfg->val($workflow,'$;NAME$;');
	    $wfname =~ s/\s//g;
	    $fullname = $wfname;
	}
    }
    my $st = stat($file);
    my ($name) = fileparse($file); 
    my ($user) = getpwuid($st->uid);
    return ($fullname,$st->mtime,$user,$name);
}

sub get_pipeline_bld_info{
    my($file) = @_;
    my $st = stat($file);
    my ($name) = fileparse($file); 
    my ($user) = getpwuid($st->uid);
    return ($st->mtime,$user,$name);
}

