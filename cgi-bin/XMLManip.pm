package XMLManip;
use strict;

use XML::Twig;

sub findxml{
    my($xmlfile,$id) = @_;
    my $myelt;
    my $t1 = new XML::Twig( TwigHandlers => { 'commandSet' =>
					      sub {
						  my ($t, $elt) = @_;
						  if(&_check_node($elt,$id)){
						      $myelt = $elt;
						  }
					      },
					      'command' =>
						  sub {
						      my ($t, $elt) = @_;
						      if(&_check_node($elt,$id)){
							  $myelt = $elt;
						      }
						  }
					  },
			    pretty_print => 'indented'
			    );

    $t1->parsefile($xmlfile);

    return $myelt;
}

sub movexml{
    my($xmlfile,$node,$id,$location) = @_;

    my $t1 = new XML::Twig( TwigHandlers => { 'commandSet' =>
					      sub {
						  my ($t, $elt) = @_;
						  if(&_check_node($elt,$id)){
						      print STDERR "Moving $id $location\n";
						      my $xml = &findxml($xmlfile,$node);
						      $xml->move($location,$elt);
						  }
						  if(&_check_node($elt,$node)){
						      $elt->cut();
						  }
					      },
					      'command' =>
						  sub {
						      my ($t, $elt) = @_;
						      if(&_check_node($elt,$id)){
							  print STDERR "Moving $id $location\n";
							  #$xml->move($location,$elt);
						      }
						  }
					  },
			    pretty_print => 'indented'
			    );

    $t1->parsefile($xmlfile);

    open FILE, "+>$xmlfile" or die("Can't open output xml file $xmlfile : $!");
    $t1->print(\*FILE);
    close FILE;
}

sub removexml{
    my($xmlfile,$id) = @_;
    my $t1 = new XML::Twig( TwigHandlers => { 'commandSet' =>
					      sub {
						  my ($t, $elt) = @_;
						  if(&_check_node($elt,$id)){
						      $elt->delete();
						  }
					      },
					      'command' =>
						  sub {
						      my ($t, $elt) = @_;
						      if(&_check_node($elt,$id)){
							  $elt->delete();
						      }
						  }
					  },
			    pretty_print => 'indented'
			    );

    $t1->parsefile($xmlfile);			
    open FILE, "+>$xmlfile" or die("Can't open output xml file $xmlfile : $!");
    $t1->print(\*FILE);
    close FILE;
}

    

sub addxml{
    my($xmlfile,$id,$newxml,$location) = @_;
    print STDERR "Attempting to add xml $newxml at $id node\n";
    my $t1 = new XML::Twig( TwigHandlers => { 'commandSet' =>
						  sub {
						      my ($t, $elt) = @_;
						      if($location eq ""){
							  $location eq "last_child";
						      }
						      if(&_check_node($elt,$id)){
							  print STDERR "Pasting xml $location $id\n";
							  $newxml->paste($location,$elt);
						      }
						  },
					      'command' =>
						  sub {
						      my ($t, $elt) = @_;
						      if($location eq ""){
							  $location eq "after";
						      }
						      if(&_check_node($elt,$id)){
							  print STDERR "Pasting xml $location $id\n";
							  $newxml->paste($location,$elt);
						      }
						  }
					  },
			    pretty_print => 'indented'
			    );

    $t1->parsefile($xmlfile);			

    open FILE, "+>$xmlfile" or die("Can't open output xml file $xmlfile : $!");
    $t1->print(\*FILE);
    close FILE;
}

sub get_pipeline_xml{
    my($pipelineid,$xmlfile,$pipelinexml) = @_;
    my $topelt;
    my $ok=0;
    my $id=0;
    my $t1 = new XML::Twig( TwigHandlers => { 'commandSet' =>
						  sub {
						      my ($t, $elt) = @_;
						      my $configmapelt = $elt->first_child("configMapId");
						      my $configmapid = $configmapelt->text();
						      if($configmapid eq "start"){
							  $topelt=$elt;
							  $configmapelt->set_text("start_"."$$");
							  my $inifile = "$xmlfile.ini";
							  my $cfg = new Config::IniFiles(-file => $inifile);
							  $cfg->newval("start_"."$$","");
							  $cfg->WriteConfig($inifile);
							  `chmod 666 $inifile`;
						      }
						      elsif($configmapid =~ /^component_/){
							  $id++;
							  my $conf;
							  my $name;
							  my $gen;
							  foreach my $commands ($elt->children('command')){
							      my $cmdconfigmapid = $commands->first_child_text("configMapId");
							      if($cmdconfigmapid =~ /^generate_component/){
								  ($conf) = ($commands->first_child_text("arg") =~ /-c\s+(\S+\.bld\.ini)/);
								  ($name) = ($configmapid =~ /component_(\w+)_\d+/);
								  $gen = $cmdconfigmapid;
							      }
							  }
							  if($conf ne "" && $name ne ""){
							      $ok=1;
							      my $componentid = "$name"."_$$".$id;
							      print STDERR "Parsing config $gen $name $conf\n";
							      my $newxml = XMLManip::get_component_xml($componentid,$pipelineid,$conf,$xmlfile);
							      $newxml->replace($elt);
							  }
							  else{
							      print STDERR "Can't parse generate component\n";
							  }
						      }
						      else{
							  my $inifile = "$xmlfile.ini";
							  my $cfg = new Config::IniFiles(-file => $inifile);
							  $cfg->newval($configmapid,"");
							  $cfg->WriteConfig($inifile);
							  `chmod 666 $inifile`;
						      }
					      }
					  },
			    pretty_print => 'indented'
			);

    $t1->parsefile($pipelinexml);			

    if($ok){
	$topelt->cut();
	return $topelt;
    }
    else{
	print STDERR "Problem finding components in $xmlfile\n";
	return undef;
    }
}


sub get_commandset_xml{
    my($type,$id,$fileprefix) = @_;

    my $inifile = "$fileprefix.ini";
    my $cfg = new Config::IniFiles(-file => $inifile);
    $cfg->newval("$id","");
    $cfg->WriteConfig($inifile);
    `chmod 666 $inifile`;

    return  parse XML::Twig::Elt( "
    <commandSet type='$type'>
    <configMapId>$id</configMapId>
    </commandSet>");
}

sub get_component_xml{
    my($componentid,$pipelineid,$conf,$fileprefix) = @_;
        
    if(&_check_conf($conf)){

	&_write_component_ini($fileprefix,$componentid);
	&_write_component_subflow_xml($fileprefix,$componentid);

	my $WorkflowBinDir = &_get_bin_dir($conf);


	return parse XML::Twig::Elt( "
       <commandSet type='serial'>
         <configMapId>component_$componentid</configMapId>
       <command>
         <type>RunUnixCommand</type>
           <configMapId>generate_component_$componentid</configMapId>
           <name>Generate component $componentid</name>
           <param>
                <key>command</key>
                <value>$WorkflowBinDir/run_pipeline</value>
           </param>
           <arg>-c $conf -i $fileprefix.subflow.$componentid.ini -sectionid run_component_$componentid --skiprun --pipelineid=$pipelineid</arg>
      </command>
      <command>
        <type>RunUnixCommand</type>
        <configMapId>create_component_$componentid</configMapId>
           <name>Create component $componentid</name>
           <param>
                <key>command</key>
                <value>$WorkflowBinDir/CreateWorkflow.sh</value>
           </param>
           <arg> -t $fileprefix.subflow_template.$componentid.xml -i $fileprefix.subflow.$componentid.xml -c $fileprefix.subflow.$componentid.ini -l $fileprefix.subflow.$componentid.log -o $fileprefix.subflow.$componentid.out</arg>
      </command>
      <commandSet type='serial'>
         <configMapId>subflow_$componentid</configMapId>
      </commandSet>
      </commandSet>");
    }
    else{
	die "Bad configuration file $conf for component";
    }
}

sub get_root_commandset{
    my($fileprefix) = @_;

    my $inifile = "$fileprefix.ini";

    my $cfg = new Config::IniFiles(-file => $inifile);
    $cfg->newval("start","");
    $cfg->WriteConfig($inifile);
    `chmod 666 $inifile`;
    
    return  parse XML::Twig::Elt( "<commandSetRoot xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance' xsi:schemaLocation='commandSet.xsd'>
            <commandSet type='serial'>
                <configMapId>start</configMapId>
            </commandSet>
            </commandSetRoot>");
}


sub get_command_xml{
    my($command,$arg,$name,$fileprefix) = @_;
    my $inifile = "$fileprefix.ini";
    my $cfg = new Config::IniFiles(-file => $inifile);
    $cfg->newval("$name","");
    $cfg->WriteConfig($inifile);
    `chmod 666 $inifile`;
    
    return  parse XML::Twig::Elt( "
    <command>
      <type>RunUnixCommand</type>
      <configMapId>$name</configMapId>
       <param>
                <key>command</key>
                <value>$command</value>
           </param>
           <arg>$arg</arg>
    </command>");
}


sub _check_node{
    my($elt,$id) = @_;

    my $configmapid = $elt->first_child("configMapId");
    my $nodeid = $elt->first_child("id");
    
    return(((defined $configmapid) && ($configmapid->text() eq $id)) || ((defined $nodeid) && ($nodeid->text() eq $id)));
}

sub _write_component_ini{
    my($fileprefix,$id) = @_;
    my $inifile = "$fileprefix.ini";
    my $cfg = new Config::IniFiles(-file => $inifile);
    $cfg->newval("generate_component_$id","");
    $cfg->newval("create_component_$id","");
    $cfg->newval("component_$id","");
    $cfg->newval("subflow_$id","fileName","$fileprefix.subflow.$id.xml");
    $cfg->WriteConfig($inifile);
    `chmod 666 $inifile`;

    my $subflowinifile = "$fileprefix.subflow.$id.ini";
    my $inicfg = new Config::IniFiles();
    $inicfg->newval("run_subflow_$id","");
    $inicfg->WriteConfig($subflowinifile);
    `chmod 666 $subflowinifile`;
}

sub _write_component_subflow_xml{
    my($fileprefix,$id) = @_;

    my $subflowxml = "<commandSetRoot xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance' xsi:schemaLocation='commandSet.xsd'>
      <commandSet type='serial'>
        <configMapId>run_subflow_$id</configMapId>
        <command> 
         <type>RunUnixCommand</type>
           <configMapId>run_component_$id</configMapId>
           <name>Run component $id</name>
         </command>
      </commandSet>
</commandSetRoot> 
";
    
    open FILE,"+>$fileprefix.subflow_template.$id.xml" or die "Can't open xml file $fileprefix.subflow_template.xml";
    print FILE "$subflowxml";
    close FILE;
}

sub _check_conf{
    my($conf) = @_;
    if(-e $conf){
	return 1;
    }
    return 0;
}

sub _get_bin_dir{
    my($file) = @_;
    my $cfg = new Config::IniFiles(-file => $file);
    return $cfg->val("init",'$;BIN_DIR$;');
}

sub _get_workflow_docs{
    my($file) = @_;
    my $cfg = new Config::IniFiles(-file => $file);
    return $cfg->val("init",'$;WORKFLOWDOCS_DIR$;');
}

1;
