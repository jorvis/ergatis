#!/usr/local/bin/perl

use strict;

use XML::Twig;

use CGI qw(:standard);
use CGI::Carp qw(fatalsToBrowser set_message);


use Tree::DAG_Node;
use File::Find;
use File::Basename;
use File::stat;
use Date::Manip;
use Data::Dumper;
use view_dag_tree;
use Digest::MD5  qw(md5 md5_hex md5_base64);

my $xmltemplate = param('xmltemplate');
my $conffile = param('conffile');
my $max_level = param('max_level');
my $open_node = param('open');
my $glob = param('glob');
my $editmode = param('edit');
my $summary = param('summary');
my $forcereload = param('forcereload');

my $trackdate=0;
my $dopurge=0;
my $printstderr=0;
my $quotawarn = "85";
my $writesubtrees = 0;
my $subtreecounts = 100;

my $quota = &get_quota($xmltemplate);

if(-d $xmltemplate){
    if($glob eq ''){
	$glob='\.xml$';
    }
    else{
	$glob=$glob.'[^\/]*\.xml';
    }
    &print_directory_page($xmltemplate,$glob,$quota);
}



my $dag_filename = &get_filename($xmltemplate);

my $currdag;
my $workflowtype;
my $olddag = undef;
if(-e $dag_filename){
    print STDERR "Found dag $dag_filename\n" if($printstderr);
    $olddag = view_dag_tree::LockRetrieve($dag_filename);
}
else{
    print STDERR "Can't find dag $dag_filename\n" if($printstderr);
}
my $newdag = Tree::DAG_Node->new();
($workflowtype) = &build_dag($xmltemplate,$newdag);

my $timings;
if($trackdate){
    $timings = &get_timings($newdag);
}
my $st = stat($xmltemplate);
my $age = time - $st->mtime;
my $user = getpwuid($st->uid);
my $last = localtime($st->mtime);
my $refreshrate = &get_refresh_rate($age);
$refreshrate = "60" if($editmode);

print header();

print "<HTML><META HTTP-EQUIV=Refresh CONTENT='$refreshrate; URL=show_pipeline.cgi?xmltemplate=$xmltemplate&summary=".param('summary')."'><title>$xmltemplate $user $last</title>\n";
my $topimage = view_dag_tree::get_commandset_image($newdag->root->attributes->{'state'});
print "<body><table><tr><td colspan=3>$topimage<b>$xmltemplate</b></td></tr><tr><td>User: $user</td><td>Last modified: $last</td><td>Age: $age s</td></tr></table>";
if($quota < $quotawarn){
    print "Project space usage <b>$quota%</b><br>";
}
else{
    print "<blink><font color=red>Project space usage <b>$quota%</b></font></blink><br>";
}
&print_dag_page($newdag,$dag_filename,$workflowtype);
print "</body></HTML>\n";

sub print_dag_page{
    my($dag,$dag_file,$type) = @_;
    $max_level = 9 if($max_level eq "" && $open_node eq "");
    my $dag_tree_obj = new view_dag_tree('link_root' => $0, 
					 'dag_obj' => $dag,
					 'dag_file' => $dag_file, 
					 'max_level' => $max_level, 
					 'xmltemplate'=>$xmltemplate,
					 'editmode'=>$editmode);
    $dag_tree_obj->open_root_nodes();  #only need to use this if have not already filled $node->attributes->{'open'} tag to initial open
    $dag_tree_obj->open_new_node($open_node) if($open_node ne "");
    $dag_tree_obj->open_nodes($dag_tree_obj->{dag_obj}) if($max_level);
    
    my $dirname = dirname($xmltemplate);
    my ($outputdir) = ($xmltemplate =~ /(.*)\/Workflow/);
    $outputdir .= "/Workflow/pipeline";

    print "<a href='./new_pipeline.cgi?&root=$outputdir'>[New]</a>&nbsp";
    print "<a href='$dag_tree_obj->{link_url_root}?&xmltemplate=$xmltemplate&summary=1'>[Summary]</a>&nbsp;";
    print "<a href='$dag_tree_obj->{link_url_root}?&xmltemplate=$xmltemplate'>[Show all nodes]</a>&nbsp;";
    print "<a href='$dag_tree_obj->{link_url_root}?&xmltemplate=$xmltemplate&forcereload=1'>[Force reparse]</a>&nbsp;";
    print "<a href='$dag_tree_obj->{link_url_root}?&xmltemplate=$dirname&glob=pipeline'>[Directory list]</a>&nbsp;";
    print "<a href='http://tools/condor-status/index.cgi'>[Condor status]</a>&nbsp;";
    print "<a href='$dag_tree_obj->{link_url_root}?&xmltemplate=$xmltemplate&edit=1'>[Edit mode]</a>&nbsp;" if($type ne "instance");
    print "Depth: <a href='",$dag_tree_obj->{link_url},"&max_level=",$max_level-1,"&xmltemplate=$xmltemplate'>[--] <a href='",$dag_tree_obj->{link_url},"&max_level=2&xmltemplate=$xmltemplate'>[2]</a> <a href='",$dag_tree_obj->{link_url},"&max_level=3&xmltemplate=$xmltemplate'>[3]</a> <a href='",$dag_tree_obj->{link_url},"&max_level=4&xmltemplate=$xmltemplate'>[4]</a> <a href='",$dag_tree_obj->{link_url},"&max_level=5&xmltemplate=$xmltemplate'>[5]</a> <a href='",$dag_tree_obj->{link_url},"&max_level=6&xmltemplate=$xmltemplate'>[6]</a> <a href='",$dag_tree_obj->{link_url},"&max_level=7&xmltemplate=$xmltemplate'>[7]</a> <a href='",$dag_tree_obj->{link_url},"&max_level=8&xmltemplate=$xmltemplate'>[8]</a> <a href='",$dag_tree_obj->{link_url},"&max_level=9&xmltemplate=$xmltemplate'>[9]</a> <a href='",$dag_tree_obj->{link_url},"&max_level=10&xmltemplate=$xmltemplate'>[10] <a href='",$dag_tree_obj->{link_url},"&max_level=",$max_level+1,"&xmltemplate=$xmltemplate'>[++]</a>";
    if($type eq "instance"){
	print "<br><a href='run_wf.cgi?instancexml=$xmltemplate'>[Run]</a>&nbsp;";
	print "<a href='kill_wf.cgi?instancexml=$xmltemplate'>[Kill]</a>&nbsp;";
    }
    else{
	my $instancexml = "$xmltemplate".".instance";
	if(-e "$instancexml"){
	    print "<br><a href='show_pipeline.cgi?xmltemplate=$instancexml'>[View instance]</a>&nbsp;";
	}
	else{
	    print "<br><a href='run_wf.cgi?xmltemplate=$xmltemplate&inifile=$xmltemplate".".ini"."&instancexml=$instancexml'>[Run]</a>&nbsp;";
	}
    }

    my $states = &get_states($dag);

    if(scalar(keys %{$states->{'commandSet'}})){
	print "<br>Commandset states:&nbsp;";
	my $totalstates = {};
	my $commandsetprogress;
	foreach my $state (keys %{$states->{'commandSet'}}){
	    $totalstates->{'commandSet'}+=$states->{'commandSet'}->{$state}->{'count'};
	}
	foreach my $state (sort {$a cmp $b} (keys %{$states->{'commandSet'}})){
	    if($states->{'commandSet'}->{$state}->{'count'} > 0){
		my $pct = sprintf("%.2f",($states->{'commandSet'}->{$state}->{'count'}/$totalstates->{'commandSet'})*100);
		print "<b>$state:</b>&nbsp;$states->{'commandSet'}->{$state}->{'count'} ($pct%) &nbsp;";
		my $image = view_dag_tree::get_commandset_image($state);
		my $imagepct = $pct*8;
		$image =~ s/height=\d+/height=10/;
		$image =~ s/width=\d+/width=$imagepct/;
		$image =~ s/border=0/border=1/;
		$commandsetprogress .= $image;
	    }
	}
	print "<br>$commandsetprogress";
    }
    if(scalar(keys %{$states->{'command'}})){
	print "<br><hr>Command states:&nbsp;";
	my $totalcommandstates=0;
	my $commandprogress;
	foreach my $state (keys %{$states->{'command'}}){
	    $totalcommandstates+=$states->{'command'}->{$state}->{'count'};
	}
	foreach my $state (sort {$a cmp $b} (keys %{$states->{'command'}})){
	    if($states->{'command'}->{$state}->{'count'} > 0){
		my $pct = sprintf("%.2f",($states->{'command'}->{$state}->{'count'}/$totalcommandstates)*100);
		print "<b>$state:</b>&nbsp;$states->{'command'}->{$state}->{'count'} ($pct%) &nbsp;";
		my $image = view_dag_tree::get_command_image($state);
		my $imagepct = $pct*8;
		$image =~ s/height=\d+/height=10/;
		$image =~ s/width=\d+/width=$imagepct/;
		$image =~ s/border=0/border=1/;
		$commandprogress .= $image;
	    }
	}
	print "<br>$commandprogress";
    }
    print "<br>";

    if($trackdate){
#	my $exectimeperjob = sprintf("%d",($totalexectime/$states->{'command'}->{'complete'}->{'count'}));
#	my $queuetimeperjob = sprintf("%d",($totalqueuetime/$states->{'command'}->{'complete'}->{'count'}));
#	my $doverheadtimeperjob = sprintf("%d",($totaldoverheadtime/$states->{'command'}->{'complete'}->{'count'}));
#	my $currtime = ($topendtime ne "") ? $topendtime : &UnixDate("today","%s");
#	my $elapsedtime = $currtime - $topstarttime;
#	print "<b>Elapsedtime</b> ",&Delta_Format(&ParseDateDelta($elapsedtime),1,"%ht"),"hrs\n";
#	print "&nbsp;<B>Exectime</b> ",&Delta_Format(&ParseDateDelta("$totalexectime"),1,"%ht"),"hrs (",&Delta_Format(&ParseDateDelta("$exectimeperjob"),1,"%mt"),"min/job)\n";
#	print "&nbsp;<b>Queuetime</b> ",&Delta_Format(&ParseDateDelta("$totalqueuetime"),1,"%ht"),"hrs (",&Delta_Format(&ParseDateDelta("$queuetimeperjob"),1,"%mt"),"min/job)\n";
#	print "&nbsp;<b>Distrib overhead</b> N/A\n";
#	print "<br>Number of hosts used ",scalar(keys %$executionhosts)+1,"\n";
#	print "&nbsp;Speedup (exectime/elapsedtime) = ";
#	printf("%.1f", $totalexectime/$elapsedtime);
#	print "x\n";
#	print "<br>";
    }

    print "<table>";
    my $runningimage = view_dag_tree::get_command_image('running');
    my $distribimage = view_dag_tree::get_command_image('running_distributed');
    my $pendingimage = view_dag_tree::get_command_image('pending');
    my $errorimage = view_dag_tree::get_command_image('error');
    my $failedimage = view_dag_tree::get_command_image('failed');
    foreach my $runningjobs (@{$states->{'command'}->{'running'}->{'nodes'}}){
	print "<tr><td>".$runningjobs->attributes->{'label'}."<a href='show_command.cgi?xmltemplate=".$runningjobs->attributes->{'xmlfile'}."&node=".$runningjobs->name()."'>[info]</a></td><td>$runningimage</td></tr>";
    } 
    foreach my $pendingjobs (@{$states->{'command'}->{'pending'}->{'nodes'}}){
	print "<tr><td>".$pendingjobs->attributes->{'label'}."<a href='show_command.cgi?xmltemplate=".$pendingjobs->attributes->{'xmlfile'}."&node=".$pendingjobs->name()."'>[info]</a><a href='show_file.cgi?&file=".$pendingjobs->attributes->{'log'}."' target='_condorlog'>[log]</a><a href='http://htcworker1:8080/antware/htcservice/html/request_display.jsp?RequestID=".$pendingjobs->attributes->{'jobid'}."' target='_condorlog'>[htcinfo]</a></td><td>$pendingimage</td></tr>";
    }
    foreach my $errorjobs (@{$states->{'command'}->{'error'}->{'nodes'}}){
	print "<tr><td>".$errorjobs->attributes->{'label'}."<a href='show_command.cgi?xmltemplate=".$errorjobs->attributes->{'xmlfile'}."&node=".$errorjobs->name()."'>[info]</a></td><td>$errorimage</td></tr>";
    }
    foreach my $failedjobs (@{$states->{'command'}->{'failed'}->{'nodes'}}){
	print "<tr><td>".$failedjobs->attributes->{'label'}."<a href='show_command.cgi?xmltemplate=".$failedjobs->attributes->{'xmlfile'}."&node=".$failedjobs->name()."'>[info]</a></td><td>$failedimage</td></tr>";
    }
    foreach my $distribjobs (@{$states->{'command'}->{'running_distributed'}->{'nodes'}}){
	print "<tr><td>".$distribjobs->attributes->{'label'}."<a href='show_command.cgi?xmltemplate=".$distribjobs->attributes->{'xmlfile'}."&node=".$distribjobs->name()."'>[info]</a></td><td>$distribimage</td></tr>";
    }
    print "</table>";
    print "<hr>";

    if(!$summary){
	$dag_tree_obj->define_dag_html(); 
    }
}
    

sub build_dag{
    my ($xmlfile,$root) = @_;
    
    my $type = "template";

    my $t1 = new XML::Twig( TwigHandlers => { 'commandSetRoot' =>
                                                  sub {
                                                      my ($t, $elt) = @_;
						      &handle_workflow($elt,$root,$xmlfile);
						  },
					      'commandSetRoot/commandSet/id' =>
						  sub {
						      $type = "instance";
						  },
					      'commandSetRoot/commandSet/startTime' =>
						  sub {
						      my ($t, $elt) = @_;
						      $root->root->attributes->{'starttime'} = &UnixDate($elt->text(),"%s");
						  },
					      'commandSetRoot/commandSet/endTime' =>
						  sub {
						      my ($t, $elt) = @_;
						      $root->root->attributes->{'endtime'} = &UnixDate($elt->text(),"%s");
						  },
					      'commandSetRoot/commandSet/state'=>
						  sub {
						      my($t, $elt) = @_;
						      $root->root->attributes->{'state'} = $elt->text();
						  }
					  },
			    ignore_elts => { 
				retryCount => 1,
				retryAttempts => 1,
				status=> 1,
				timeOut=>1,
				parentFileName=>1,
			    });
    print STDERR "Parsing $xmlfile\n" if($printstderr);
    $t1->parsefile($xmlfile);  
    $t1->purge() if($dopurge);

    return ($type);
}


my $xmllookup;

sub handle_workflow{
    my($cs,$node,$xmlfile) = @_;
    my @children = $cs->children();
    foreach my $cs (@children){

	if($cs->gi eq "commandSet"){
	    my $skipparse=0;
	    if($olddag){
		my $statenode = $cs->first_child_text('state');
		my($nodeid,$name) = &get_id($cs);
		$skipparse = &handle_subnode($olddag,$node,$nodeid,$statenode);
	    }
	    else{
	    }
	    if($skipparse == 0){
		&handle_commandset($cs,$node,$xmlfile);
	    }
	}elsif($cs->gi eq "command"){
	    &handle_command($cs,$node,$xmlfile);
	}
    }
}

sub handle_commandset{
    my($cs,$node,$xmlfile) = @_;
    my ($nodeid,$name) = &get_id($cs);

    my $state = $cs->first_child_text('state');

    my $type = $cs->att('type');


    my $new_daughter = $node->new_daughter();
    
    $new_daughter->name("$nodeid");
    #store lookup by _key to node
    $node->root()->attributes()->{'nodelookup'}->{$new_daughter->name()} = $new_daughter;

    $new_daughter->attributes()->{'open'} = 1;
    $new_daughter->attributes()->{'xmlfile'} = "$xmlfile";
    $new_daughter->attributes()->{'state'} = "$state";
    $new_daughter->attributes()->{'elt'} = "commandSet";
    $new_daughter->attributes()->{'type'} = "$type";
    $new_daughter->attributes()->{'label'} = "$name";

    my $file = $cs->first_child('fileName');
    
    if(defined $file){
	my $incxmlfile = $file->text();
	if($incxmlfile ne $xmlfile){
	    &handle_filename($incxmlfile,$new_daughter);
	}
	else{
	    &handle_workflow($cs,$new_daughter,$xmlfile);
	}
    }
    else{
	&handle_workflow($cs,$new_daughter,$xmlfile);
    }
}

sub handle_command{
    my($cs,$node,$xmlfile) = @_;

    my ($nodeid,$name) = &get_id($cs);

    my $state = $cs->first_child_text('state');
    my $type = $cs->first_child_text('type');
    my $configmapid = $cs->first_child_text('configMapId');
    
    my $new_daughter = $node->new_daughter;

    if($trackdate){
	my $starttime= 0;
	my $endtime=0;
	if($cs->first_child('startTime')){
	    $starttime = $cs->first_child('startTime')->text();
	}
	if($cs->first_child('endTime')){
	    $endtime = $cs->first_child('endTime')->text();
	}
	$starttime = &UnixDate($starttime,"%s");
	$endtime = &UnixDate($endtime,"%s");
	$new_daughter->attributes->{'starttime'} = $starttime;
	$new_daughter->attributes->{'endtime'} = $endtime;
    }
    
    my $log;
    my $jobid;
    if($type =~ /Distributed/){ 
	my $dce = $cs->first_child('dceSpec');
	my $host;
	if($dce){
	    $log = $dce->first_child_text('log');
	    $jobid = $dce->first_child_text('jobID');
	    $host = $dce->first_child_text('executionHost');
	    $new_daughter->attributes()->{'log'} = $log;
	    $new_daughter->attributes()->{'host'} = $host;
	    $new_daughter->attributes()->{'jobid'} = $jobid;

	    if(-e $log){
		if($trackdate){
		    my ($submittime,$executionstart,$executionend) = _parse_condor_logs($log);
		    $new_daughter->attributes()->{'submittime'} = "$submittime";
		    $new_daughter->attributes()->{'executionstart'} = "$executionstart";
		    $new_daughter->attributes()->{'executionend'} = "$executionend";
		}		
	    }
	    if($dce && $state eq 'running'){
		$state = 'running_distributed';
	    }
	}
    }

    $new_daughter->name("$nodeid");
    #store lookup by _key to node
    $node->root->attributes->{'nodelookup'}->{$new_daughter->name()} = $new_daughter;
    $new_daughter->attributes()->{'open'} = 1;
    $new_daughter->attributes()->{'xmlfile'} = "$xmlfile"; 
    $new_daughter->attributes()->{'state'} = "$state";
    $new_daughter->attributes()->{'elt'} = "command";
    $new_daughter->attributes()->{'type'} = "$type";
    $new_daughter->attributes()->{'label'} = "$name";

    if($configmapid =~ /^generate_component/){
	my $arg = $cs->first_child_text('arg');
	my ($conf) = ($arg =~ /\-c\s+(\S+)/);
	$new_daughter->attributes()->{'conf'} = $conf;
    }


    if($name =~ /Run/){ 
	my @paramchilds = $cs->children('param');
	foreach my $p (@paramchilds){
	    if($p->first_child('key')->text() eq '--instance'){
		my $dir = dirname($p->first_child('value')->text());
		$new_daughter->attributes()->{'conf'} = "$dir/pipeline.config";
		if(&handle_filename($p->first_child('value')->text(),$new_daughter)){
		    #sub workflow is created
		}
		else{
		}
	    }
	}
    }
}

sub handle_filename{
    my($xmlfile,$node) = @_;
    if(-e $xmlfile){
	$xmllookup->{$xmlfile} = 1;
	my $statenode = $node->attributes->{'state'};
	my $skipparse = &handle_subnode($olddag,$node,$node->name(),$statenode);
	if($skipparse==0){
	    my $dag_subfilename = &get_filename($xmlfile);
	    if(-e $dag_subfilename){
		my $oldsubdag =  view_dag_tree::LockRetrieve($dag_subfilename);
		$skipparse = &handle_subnode($oldsubdag,$node,$node->name(),$statenode);
	    }
	}
	if($skipparse ==0){
	    my $t2 = new XML::Twig( TwigHandlers => { 'commandSetRoot' =>
							  sub {
							      my ($t, $elt) = @_;
							      &handle_workflow($elt,$node,$xmlfile);
							  },
						  },
				    ignore_elts => { 
					retryCount => 1,
					retryAttempts => 1,
					status=> 1,
					timeOut => 1,
					parentFileName=>1,
				    }
				    );
	    print STDERR "Parsing $xmlfile $statenode\n" if($printstderr);
	    $t2->parsefile($xmlfile);  
	    $t2->purge() if($dopurge);
	    
	    if(scalar($node->descendants) > $subtreecounts){
		$node->attributes->{'write'} = 1;
	    }

	    return 1;
	}
	else{
	    return 1;
	}
    }
    else{
	return 0;
    }
}


sub trimfile{
    my($file) = @_;
    my ($trim) = ($file =~ /.*Workflow(.*)/);
    return $trim;
}

sub _parse_condor_logs{
    my($log) = @_;
    my $submittime = 0;
    my $executionstart = 0;
    my $executionend = 0;

    open LOG, "$log" or warn("Can't open log file $log for parsing $?");
    while (my $line = <LOG>){
	if($line =~ /Job submitted/){
	    ($submittime) = ($line =~ /([\d\/]+\s+[\d\:]+)\s+Job submitted/);
	}
	elsif($line =~ /Job executing/){
	    ($executionstart) = ($line =~ /([\d\/]+\s+[\d\:]+)\s+Job executing/);
	}
	elsif($line =~ /Job terminated/){
	    ($executionend) = ($line =~ /([\d\/]+\s+[\d\:]+)\s+Job terminated/);
	}
    }
    close LOG;
    return (&UnixDate($submittime,"%s"),&UnixDate($executionstart,"%s"),&UnixDate($executionend,"%s"));
}
    
sub get_id{
    my($twig) = @_;
    my $cidchild = $twig->first_child('configMapId');
    my $id = $twig->first_child('id');
    my $namenode = $twig->first_child('name');

    my $nodeid;
    my $name;

    if(defined $id){
	$nodeid = $id->text();
    }
    elsif(defined $cidchild){
	$nodeid = $cidchild->text();
    }

    if(defined $namenode){
	$name = $namenode->text();
    }

    my $n = ($name eq "") ? $nodeid : $name;    

    return ($nodeid,$n);
}

sub get_workflows_from_directory{
    my($dir,$glob) = @_;
    my $files = {};
    find(sub {my $file = $File::Find::name;
	      my $currdir = dirname($file);
	      my $currglob = $glob;
	      if($currdir =~ /pipeline$/){
		  $currglob = "$glob"."\.instance";
	      }
	      if($file =~ /$currglob$/){
		  if(-e $file){
		      $files->{$file}->{'filename'} = "$file";
		      my $dir = dirname($file);
		      if($dir =~ /\d+$/){
			  ($dir) = ($dir =~ /(.*)\/\d+/);
		      }
		      $files->{$file}->{'dirname'} = $dir;
		      my $st = stat($files->{$file}->{'filename'});
		      $files->{$file}->{'size'} = $st->size;
		      $files->{$file}->{'date'} = $st->mtime;
		      $files->{$file}->{'user'} = getpwuid($st->uid);
		      my $root = Tree::DAG_Node->new();
		      my $state = 'unknown';
		      my $t1 = new XML::Twig( TwigHandlers => { 'commandSetRoot/commandSet/state' =>
								    sub {
									my ($t, $elt) = @_;
									$state = $elt->text();
								    },
							    });
		      $t1->parsefile($files->{$file}->{'filename'}); 
		      $t1->purge() if($dopurge);
		      
		      $files->{$file}->{'state'} = $state;
		  }
	      }
	  },$xmltemplate);
    return $files;
}

sub print_directory_page{
    my($dir,$glob) = @_;

    print header();
    my $dirname = dirname($dir);
    my ($outputdir) = ($dir =~ /(.*)\/Workflow/);
    $outputdir .= "/Workflow/pipeline";
    print "<html><META HTTP-EQUIV=Refresh CONTENT='120; URL=show_pipeline.cgi?xmltemplate=$dir&glob=".param('glob')."'><h2>Workflows in $dir</h2>";
    if($quota < $quotawarn){
	print "Project space usage <b>$quota%</b><br>";
    }
    else{
	print "<blink><font color=red>Project space usage <b>$quota%</b></font></blink><br>";
    }
    print "<br>Parent directory $dirname&nbsp;<a href='./show_pipeline.cgi?&xmltemplate=$dirname&glob=".param('glob')."'>[view]</a><br>";
    print "<a href='./new_pipeline.cgi?&root=$outputdir'>[New]</a>&nbsp;<hr>";

    my $files = &get_workflows_from_directory($dir,$glob);
    
    print "<table>";
    my $prevdir;
    foreach my $file (sort {
	if($files->{$a}->{'dirname'} eq $files->{$b}->{'dirname'}){
	    $files->{$b}->{'date'} cmp $files->{$a}->{'date'}
	}
	else{
	    $files->{$a}->{'dirname'} cmp $files->{$b}->{'dirname'};
	}} (keys %$files)){ 
	if($files->{$file}->{'dirname'} ne $prevdir){
	    print "<tr><td colspan=6><b>$files->{$file}->{'dirname'}</b>&nbsp;<a href='./show_pipeline.cgi?&xmltemplate=$files->{$file}->{'dirname'}&glob=".param('glob')."'>[view]</a></td></tr>";
	}
	$prevdir = $files->{$file}->{'dirname'};
	my $image = view_dag_tree::get_commandset_image($files->{$file}->{'state'});
	print "<tr><td><a href='show_pipeline.cgi?&xmltemplate=$files->{$file}->{'filename'}'>$files->{$file}->{'filename'}</a></td><td>$image</td><td>$files->{$file}->{'state'}</td><td>$files->{$file}->{'user'}</td><td>".localtime($files->{$file}->{'date'})."</td><td><a href='show_pipeline.cgi?&xmltemplate=$files->{$file}->{'filename'}&summary=1'>[summary]</a>&nbsp;<a href='remove_workflow.cgi?file=$files->{$file}->{'filename'}&removedir=1' target='_new'>[delete]</a>&nbsp;<a href='kill_wf.cgi?instancexml=xmltemplate=$files->{$file}->{'filename'}'>[kill]</a></td></tr>\n";
    }
    print "</table>";
    print "</html>";
    exit;
}

sub get_timings{
    my($dag) = @_;
#queueing time
    my $totalqueuetime = 0;
#est. exectime
    my $totalexectime = 0;
#est. distributed overhead
    my $totaldoverheadtime = 0;
    my $executionhosts = {};
    my $topstarttime;
    my $topendtime;
    my $topstate;
}

sub get_states{
    my ($dag) = @_;
    my $states = {};
    print STDERR "Rewriting states and nodelookup ".scalar($dag->descendants)."\n" if($printstderr);
    $dag->walk_down({ callbackback => sub
		       {
			   my ($node) = @_;
			   if($node && $node->name()){
#			       die "Node ".$node->attributes->{'label'}." ".$node->attributes->{'xmlfile'}." node ".$node->attributes()->{'_key'}." ".$node->attributes->{'state'} if ($node->attributes()->{'_key'} eq "2907766");
			       $node->root->attributes->{'nodelookup'}->{$node->name()} = $node;
			       
			       if($node->attributes->{'elt'} eq "commandSet"){
				   $states->{'commandSet'}->{$node->attributes->{'state'}}->{'count'}++;
				   
			       }
			       elsif($node->attributes->{'elt'} eq "command"){
				   $states->{'command'}->{$node->attributes->{'state'}}->{'count'}++;
				   if(!(exists $states->{'command'}->{$node->attributes->{'state'}}->{'nodes'})){
				       $states->{'command'}->{$node->attributes->{'state'}}->{'nodes'} = [];
				   }
				   push @{$states->{'command'}->{$node->attributes->{'state'}}->{'nodes'}},$node;
				   if($node->attributes->{'host'}){
				       $states->{'executionhosts'}->{$node->attributes->{'host'}}++;
				   }
			       }

			       if($node->attributes->{'write'}){
				   if($writesubtrees){
				       my $dag_subfilename = &get_filename($node->attributes->{'xmlfile'});
				       view_dag_tree::LockStore($node,$dag_subfilename);
				   }
			       }
			   }
		       }
		   });
    return $states;
}

sub get_refresh_rate{
    my($age) = @_;
    if($age < 60){
	$refreshrate = 10;
    }
    elsif($age < 300){
	$refreshrate = 60;
    }
    elsif($age < 600){
	$refreshrate = 120;
    }
    else{
	$refreshrate = 300;
    }
    return $refreshrate;
}

sub get_filename{
    my ($file) = @_;
    my $md5 = md5_hex($file);
    return "/tmp/$md5.dag";
}


sub get_quota{
    my($xmlfile) = @_;
    my($projectdir) = ($xmlfile =~ /(\/usr\/local\/annotation\/[^\/]+)/);
    my $out = `getquota $projectdir`;
    my ($limit) = ($out =~ /Current\s+Quota\s+Limit\s+:\s+(\d+)/gi);
    my ($current) = ($out =~ /Current\s+Space\s+Usage\s+:\s+(\d+)/gi);
    my $pct;
    if($current && $limit){
        $pct = sprintf("%.1f",($current/$limit)*100);
    }
    return $pct;
}

sub handle_subnode{
    my($olddag,$node,$nodeid,$statenode) = @_;
    my $skipparse = 0;
    my $skipstates = {'complete'=>1,
		      'interrupted'=>1,
		      'error'=>1,
		      'failed'=>1,
		      'pending'=>1,
		      'incomplete'=>1};
    if($olddag){
	if($forcereload != 1){
	    if($skipstates->{$statenode}){
		my $replacenode = $olddag->root->attributes->{'nodelookup'}->{$nodeid};
		if($replacenode){
		    print STDERR "Found replace node $nodeid with state ".$replacenode->attributes->{'state'} ." currstate:$statenode\n" if($printstderr);
		    if($statenode eq $replacenode->attributes->{'state'}){
			print STDERR "Skipping subnode $nodeid and replacing with previous parse nodeid:".$replacenode->name()." state:".$replacenode->attributes->{'state'}.". calling node:".$node->name(),"\n" if($printstderr);
			if(!$node->mother){
			    print STDERR "Replacing root node with ".$replacenode->name()." ".scalar($replacenode->descendants())." descendents\n" if($printstderr);
			    if($node->name() eq $nodeid){
				$node->replace_with($replacenode);
			    }
			    else{
				$node->add_daughters($replacenode);
			    }
			}
			else{
			    if($node->name() eq $nodeid){
				$node->replace_with($replacenode);
			    }
			    else{
				$node->add_daughters($replacenode);
			    }
			}
			$skipparse = 1;
		    }
		}
		else{
		    print STDERR "Can't find node ".$nodeid." in oldday lookup\n" if($printstderr);
		}
	    }
	    else{
		print STDERR "Node $nodeid marked $statenode and needs reparse. calling node ".$node->name()."\n" if($printstderr);
	    }
	}
    }
    return $skipparse;
}
