#!/usr/bin/env perl

eval 'exec /usr/bin/env perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

#Create a complete workflow XML file (--output_xml) from a template
#XML file (--template_xml) that contains a set of placeholders/keys.
#The config file --component_conf provides a list of key=value pairs.
#Each key found in the template file is replaced by its corresponding
#value from the config file

use lib (@INC,$ENV{"PERL_MOD_DIR"});
use strict;
use Config::IniFiles;
use Ergatis::Logger;
use Getopt::Long;
use XML::LibXML;

#--template_xml   Input template XML file for a component 
#--iterator_list  Iterator list.  If iterator list specified then must specify iterator_output_dir
#--component_conf Final component config file
#--output_xml     Output XML

my %options;
my $delimeter = '$;';
my $SKIPTAG='$;SKIP_WF_COMMAND$;';
## TODO: Figure out a better way to do this mapping
my $dce_param_map = { 'OS' => "OS", 'EXECUTION_HOST' => "executionHost", 
                      'WORKINGDIR' => "workingDir", 'reqStartTime' => "reqStartTime" };
my $default_memory='100M';                      
                      
my $results = GetOptions (\%options, 
                          'keys=s',
                          #opts for replacing keys in a single template file
                          'template_xml|t=s', 
                          'template_xml_conf_key|k=s',
                          'output_xml|o=s', 
                          'component_conf|c=s',
                          #opts for replacing keys in an iterator list
                          'output_dir|O=s',
                          'iterator_list|l=s', 
                          'distribopts|d=s',
                          'log=s',
                          'debug=s', 
                          'help|h' );

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

my $cfg = new Config::IniFiles( -file => $options{'component_conf'});
#Add add'l keys specified via --keys

&add_keys($cfg,'component',split(/,/,$options{'keys'}));

my $replacevals = &parseconf($cfg);

if (exists $options{'iterator_list'}) {
    my ($iteratorconf,$iterator_list_key) = &parseiteratorconf($options{'iterator_list'});
    my $outputfiles;
    my @iterator_output_list;
    
    if (!exists $iteratorconf->{$iterator_list_key}) {
        $logger->logdie("Can't find key $iterator_list_key in $options{'iterator_list'}");
    } else {
        $outputfiles = $iteratorconf->{$iterator_list_key};
    }
    
    if (! -d $options{'output_dir'}) {
        $logger->logdie("$options{'output_dir'} is not a valid output directory");
    }
    
    $options{'output_dir'} .= "/" if ($options{'output_dir'} !~ /\/$/);
    
    for (my $i=0; $i < @$outputfiles; $i++) {
        my %ireplacevals = %$replacevals;
        
        foreach my $key (keys %$iteratorconf){
            $ireplacevals{$key} = $iteratorconf->{$key}->[$i];
        }
        
        # If iterator_list is provided, create all group XML using XML template, and replace values,
        my $outputfile = $options{'output_dir'}.$iteratorconf->{$iterator_list_key}->[$i].".xml.gz";
        push @iterator_output_list,$outputfile;
        
        if(! -e $options{'template_xml'}){
	    	$logger->logdie("Can't find template xml file ($options{'template_xml'})");     
		}
        &replacekeys(\%ireplacevals,$options{'template_xml'},$outputfile);
    }
    
    my $distrib;
    my $maxpar;
    my $cstype;
    if ($options{'distribopts'}) {
        if ($options{'distribopts'} =~ /nodistrib=(\d+)/i && $1>0) {
            $distrib = "parallel";
            $maxpar = "<maxParallelCmds>$1</maxParallelCmds>";
            $cstype = "serial";
	} elsif ($options{'distribopts'} =~ /serial/) {
	    $distrib = "serial";
	    $cstype = "serial";
	} else {
            $distrib = "parallel";
            $cstype = "remote-serial";
        }
    }
    
    if ( $options{output_xml} =~ /\.gz$/ ) {
        open(FILE, ">:gzip", $options{output_xml} ) or $logger->logdie("Can't open file $options{'output_xml'}");
    } else {
        open FILE, ">$options{'output_xml'}" or $logger->logdie("Can't open file $options{'output_xml'}");
    }
    
    # Modify iterator XML to include each group XML and other information
    print FILE "<?xml version=\"1.0\" encoding=\"UTF-8\"?>
                <commandSetRoot type=\"instance\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation='commandSet.xsd'>
                <commandSet type=\"$distrib\">
                 <state>incomplete</state>
                 $maxpar\n";
                 
    for(my $i=0;$i<@iterator_output_list;$i++){
        my $file = $iterator_output_list[$i];
        my $name = $iteratorconf->{$iteratorconf->{$iterator_list_key}->[$i]}->[$i] || '';
        print FILE "<commandSet type=\"$cstype\">
                     <name>$name</name>
                     <state>incomplete</state>
                     <id>$i</id>
                     <fileName>$file</fileName>\n";

        print FILE get_dce_spec_parameters($cfg);
        print FILE " </commandSet>\n";
    }
    print FILE "</commandSet>
               </commandSetRoot>\n";
    close FILE; 
}
else{
	# Used to grab component XML from the pipeline
    if($options{'template_xml_conf_key'}){
	my $templatexml = $cfg->val('component',$delimeter.$options{'template_xml_conf_key'}.$delimeter);
	
	if(!$templatexml){
	    $logger->logdie("Can't find key in component config $options{'component_conf'} [component] $options{'template_xml_conf_key'}");     
	}	    
	if(! -e $templatexml){
	    $logger->logdie("Can't find referenced xml file ($templatexml) in component config $options{'component_conf'} [component] $options{'template_xml_conf_key'}");     
	
	}
	&replacekeys($replacevals,$templatexml,$options{'output_xml'});
    }
    else{
    # Used for grabbing nested XML or other non-iterative XML
	if(! -e $options{'template_xml'}){
	    $logger->debug("Can't find template xml file ($options{'template_xml'})") if($logger->is_debug());     
	}
	&replacekeys($replacevals,$options{'template_xml'},$options{'output_xml'});
    }
}

#--------------------------------------------------------------------
# Replaces any DCE spec parameters that have been set by the user 
# in the interface
#--------------------------------------------------------------------
sub get_dce_spec_parameters {
    my $cfg = shift;
    my $dce_block = "<dceSpec type=\"sge\">\n";

    if ($cfg->SectionExists('dce')) {
        for my $parameter ($cfg->Parameters('dce')) {
            my $val = $cfg->val('dce', $parameter);

            if ($val) {
                $parameter =~ s/\$;//g;
                $parameter = exists($dce_param_map->{$parameter}) ? $dce_param_map->{$parameter} : lc($parameter);
                $dce_block .= "<$parameter>$val</$parameter>\n";
            }
        }
    }
        
    ## Kind of kludgy but we want to check if we have the OS parameter 
    ## in our dceSpec block as this is always required. If not we add it
    $dce_block .= "<OS>linux</OS>\n" if ($dce_block !~ /<OS>\w+<\/OS>/);
    
    ## At this point we may also be missing a group (could be critical 
    ## depending on your SGE install) and will want to add it.
    if ($dce_block !~ /<group>.+<\/group>/ && 
        $cfg->val('project', '$;PROJECT_CODE$;')) {
        $dce_block .= "<group>" . $cfg->val('project', '$;PROJECT_CODE$;') . 
                      "</group>\n";
    }

	### SAdkins - 1/11/17 - Commenting this out.  We are worried this would allow users to be lazy and not specify a memory requirement, and grossly underestimate their requirements, and causing chaos on the grid.
	# Required to specify mem_free requirements when submitting to the grid.  This is to ensure jobs without that req have a default value.
	#$dce_block .= "<passthrough>-l mem_free=$default_memory</passthrough>\n" if ($dce_block !~ /<passthrough>[\w=\s-_]+<\/passthrough>/);

    $dce_block .= "</dceSpec>\n";
    return $dce_block;
}

#
#Create substitution lookup for placeholder keys
#Reads an INI formatted config file with key=value
#Lookup is $subs->{$key}=$value
#Key with name "$;SKIP_WF_COMMAND$;" are handled differently
#The key becomes $subs->{$;SKIP_WF_COMMAND$;$value}=1
sub parseconf{
    my($conf) = @_;
    my $subs = {};
    my @sections = $conf->Sections();
    foreach my $section (@sections){
	$logger->debug("Checking section [$section]") if($logger->is_debug());
	my @parameters = $conf->Parameters($section);
	foreach my $parameter (@parameters){
	    my $replacevalue = $conf->val($section,$parameter);
	    $logger->debug("Storing $parameter as $replacevalue") if($logger->is_debug());
	    if($parameter eq $SKIPTAG){
		foreach my $name (split(/,/,$replacevalue)){
		    $subs->{$SKIPTAG.$name}=1;
		}
	    }
	    else{
		$subs->{$parameter} = $replacevalue;
	    }
	}
    }
    return $subs;
}

#
#Process XML $doc with libxml and remove <command> elements
#that are marked in the substitution lookup $subs.
sub skipcommands{
    my($subs,$doc) = @_;

    my $query  = "//command/name";    
    my(@nodes)  = $doc->findnodes($query);
    
    foreach my $node (@nodes){
	my $name = $node->textContent;
	if(exists $subs->{$SKIPTAG.$name}){
	    print STDERR '<!--Skipping command with name ',$name,">\n";
	    my $commandnode = $node->parentNode;
	    my $cmdparentnode = $commandnode->parentNode;
	    $cmdparentnode->removeChild($commandnode);
	}
	else{
	}
    }
}

#
#Parse file $inputfile and make substitutions as defined
#in lookup $subs
#Write output to file $outputfile
sub replacekeys{
    my($subs,$inputfile,$outputfile) = @_;

    $logger->debug("Reading template $inputfile and writing output $outputfile");
    open( INPUTFILE, "$inputfile" ) or $logger->logdie("Could not open input template file $inputfile");
    
    if ( $outputfile =~ /\.gz$/ ) {
        open( OUTPUTFILE, ">:gzip", $outputfile) or $logger->logdie("Could not open output template file $outputfile: $!");
    } else {
        open( OUTPUTFILE, "+>$outputfile") or $logger->logdie("Could not open output template file $outputfile: $!");
    }
    my $run_dist_cmd_flag = 0;
    my $remote_serial_flag = 0;

    my @inputxml;
    #keeping old code that parses XML in a simplistic way... but works
    #TODO refactor to use LibXML
    while( my $line = <INPUTFILE> ){
    print $line;
        ## don't replace vals on comment lines
	if( $line !~ /^\;/ && $line !~ /^\#/ ) {        
	    if($line =~ /\<INCLUDE/){
		$line =~ s/(\$;[\w_]+\$;)([^=])/&replaceval($1,$subs).$2/ge;
	    }
	    else{
		$line =~ s/(\$;[\w_]+\$;)/&replaceval($1,$subs)/ge;	    
	    }
        }
	if($line =~ /\<type\>RunDistributedCommand\<\/type\>/) {
	    my $group = &replaceval('$;PROJECT_CODE$;',$subs);
	    $run_dist_cmd_flag = 1 if( $group || $group == 0 );
	}
	if($run_dist_cmd_flag && $line =~ /<dceSpec/) {
	    my $group = &replaceval('$;PROJECT_CODE$;',$subs);
	    $line .= "            <group>$group</group>\n";
	    $run_dist_cmd_flag = 0;
	}
	if($run_dist_cmd_flag && $line =~ /\<\/command\>/ ) {
	    my $group = &replaceval('$;PROJECT_CODE$;',$subs);
        my $new_line = get_dce_spec_parameters($cfg) . $line;
	    $line = $new_line;
	    $run_dist_cmd_flag = 0;
	}
	
	# To prevent excess dceSpec blocks from being created, just run on remote-serial commandSet types
	if ($line =~/\<commandSet type=\"remote-serial\"\>/ ) {
		$remote_serial_flag = 1;
	}
	
	# If we encounter a commandSeq tag
	if ($remote_serial_flag && $line =~ /\<\/commandSet\>/) {
		
		# "Peek" next line, and if commandSet is just before end of file commandSetRoot, just ignore
		# commandSetRoot is always the end of file, and a commandSet is always nested in it
		#my $pos = tell(INPUTFILE);
		#my $tmp_line = <INPUTFILE>;
		#chomp $tmp_line;
		#seek(INPUTFILE, $pos, 0);
		#if ($tmp_line eq "</commandSetRoot>") {
		#	push @inputxml, $line;
		#	next();
		#}
		
        my $new_line = get_dce_spec_parameters($cfg) . $line;
	    $line = $new_line;		
	    
	    # If dceSpec block preceded commandSet closing tag, then remove the block
	    if ($inputxml[-1] =~ /\<\/dceSpec\>/ ) {
	    	my $pop_value;
	    	do {
	    		$pop_value = pop @inputxml;
	    	} until ($pop_value =~ /\<dceSpec/ );
	    }
	    $remote_serial_flag = 0;
	}
	
	if($line =~ /\<INCLUDE/){
	    if($line !~ />/){
		$logger->logdie("<INCLUDE> directive must be contained on a single line");
	    }
	    my($file) = ($line =~ /file="(.+?)"/);
	    my($keys) = ($line =~ /keys="(.*)"/);
	    if($keys !~ /ITERATOR_RANDOM/){
		$keys .= ',$;ITERATOR_RANDOM$;=0'
	    }
	    if($keys !~ /ITERATOR_TIMESTAMP/){
		$keys .= ',$;ITERATOR_TIMESTAMP$;=0'
	    }
	    $logger->debug("Include file $file using keys $keys: $line") if($logger->is_debug());
	    &import_xml($file,$keys,$subs,\@inputxml);
	    $line = '';
	}
	push @inputxml,$line;
	#print OUTPUTFILE $line;
    }
    
    close INPUTFILE;
    my $parser = XML::LibXML->new();
    my $xmldoc    = $parser->parse_string(join('',@inputxml));
    &skipcommands($subs,$xmldoc);
    print OUTPUTFILE $xmldoc->toString;
    close OUTPUTFILE;
}

sub replaceval{
    my($val,$keylookup) = @_;
    if(!(exists $keylookup->{$val})){
	$logger->logdie("Bad key $val in template file");
    }
    else{
	$logger->debug("Replacing $val with $keylookup->{$val}");
	return $keylookup->{$val};
    }
}

sub parseiteratorconf{
    my ($listfile) = @_;
    my $iteratorconf = {};

    open FILE, $listfile or $logger->logdie("Can't open file $listfile");
    my @keys;
    my $linenum=0;
    my $uniquekeys = {};
    while(my $line=<FILE>){
	chomp $line;
	if($linenum==0){
	    @keys = split(/\t/,$line);
	    foreach my $key (@keys){
		$iteratorconf->{$key} = [];
	    }
	}
	else{
	    my @elts = split(/\t/,$line);
	    for(my $i=0;$i<(@elts);$i++){
		#make sure first element is globally unique
		#or will cause problems
		if($i==0){
		    if(exists $uniquekeys->{$elts[$i]}){
			$logger->logdie("Found duplicate key $elts[$i] in first column. First column must be a unique identifier");
		    }
		    $uniquekeys->{$elts[$i]} = 1;
		}
		push( @{$iteratorconf->{$keys[$i]}}, $elts[$i] );
	    }
	}
	$linenum++;
    }
    close FILE;
    #make first key the name of the iterator
    return ($iteratorconf,$keys[0]);
}

sub add_keys{
    my($cfg,$section,@keys) = @_;
    $logger->debug("Adding user defined keys: ".@keys) if($logger->is_debug());
    foreach my $kv (@keys){
	my($key,$value) = split(/=/,$kv);
	$logger->debug("Adding user defined key $key=$value in section [$section]") if($logger->is_debug());
	if($cfg->setval($section,$delimeter.$key.$delimeter,$value)){
	}
	else{
	    $cfg->newval($section,$delimeter.$key.$delimeter,$value);
	}
    }
}

sub import_xml{
    my($file,$keys,$subs,$outputh) = @_;
    my @keys = split(/,/,$keys);
    foreach my $k (@keys){
	my($tok,$val) = split(/=/,$k);
	if(exists $subs->{$delimeter.$val.$delimeter}){
	    $subs->{$tok} = $subs->{$delimeter.$val.$delimeter};
	}
	else{
	    $subs->{$tok} = $val;
	}
    }
    my $fh;
    open $fh, "$file" or $logger->logdie("Can't open file $file");
    my $run_dist_cmd_flag = 0;
    my $remote_serial_flag = 0;
    my @inputxml;
    while(my $line=<$fh>){
	#keeping old code that parses XML in a simplistic way... but works
	#TODO refactor to use LibXML
	if( $line !~ /^\;/ && $line !~ /^\#/ ) {
	    if($line =~ /\<INCLUDE/){
		$line =~ s/(\$;[\w_]+\$;)([^=])/&replaceval($1,$subs).$2/ge;
	    }
	    else{
		$line =~ s/(\$;[\w_]+\$;)/&replaceval($1,$subs)/ge;	    
	    }
        }
	
	if($line =~ /\<type\>RunDistributedCommand\<\/type\>/) {
	    my $group = &replaceval('$;PROJECT_CODE$;',$subs);
	    $run_dist_cmd_flag = 1 if( $group || $group == 0 );
	}
	if($run_dist_cmd_flag && $line =~ /\<dceSpec/) {
	    my $group = &replaceval('$;PROJECT_CODE$;',$subs);
	    $line .= "            <group>$group</group>\n";
	    $run_dist_cmd_flag = 0;
	}
	if($run_dist_cmd_flag && $line =~ /\<\/command\>/ ) {
        my $new_line = get_dce_spec_parameters($cfg) . $line;
	    $line = $new_line;
	    $run_dist_cmd_flag = 0;
	}
	
	if ($line =~/\<commandSet type=\"remote-serial\"\>/ ) {
		$remote_serial_flag = 1;
	}
	
	# If we encounter a commandSeq tag
	if ($remote_serial_flag && $line =~ /\<\/commandSet\>/) {
		
		# "Peek" next line, and if commandSet is just before end of file commandSetRoot, just ignore
		#my $pos = tell($fh);
		#my $tmp_line = $fh;
		#chomp $tmp_line;		
		#seek($fh, $pos, 0);
		#if ($tmp_line eq "</commandSetRoot>") {
		#	push @inputxml, $line;
		#	next();
		#}
	
        my $new_line = get_dce_spec_parameters($cfg) . $line;
	    $line = $new_line;		
	    
	    # If dceSpec block preceded commandSet closing tag, then remove the block
	    if ($inputxml[-1] =~ /\<\/dceSpec\>/ ) {
	    	my $pop_value;
	    	do {
	    		$pop_value = pop @inputxml;
	    	} until ($pop_value =~ /\<dceSpec/ );
	    }
	    $remote_serial_flag = 0;
	}
		
	if($line =~ /\<INCLUDE/){
	    if($line !~ />/){
		$logger->logdie("<INCLUDE> directive must be contained on a single line");
	    }
	    my($file) = ($line =~ /file="(.+?)"/);
	    my($keys) = ($line =~ /keys="(.*)"/);
	    if($keys !~ /ITERATOR_RANDOM/){
		$keys .= ',$;ITERATOR_RANDOM$;=0'
	    }
	    if($keys !~ /ITERATOR_TIMESTAMP/){
		$keys .= ',$;ITERATOR_TIMESTAMP$;=0'
	    }
	    $logger->debug("Include file $file using keys $keys: $line") if($logger->is_debug());
	    &import_xml($file,$keys,$subs,$outputh);
	    $line = '';
	}
	#print $outfh $line;
	push @$outputh,$line;
    }
    close $fh;
}

1;

