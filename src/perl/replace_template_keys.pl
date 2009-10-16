#!/usr/bin/perl

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
my $SKIPTAG='$;SKIP_WF_COMMAND$;';
my $results = GetOptions (\%options, 
			  'keys=s',
			  #opts for replacing keys in a single template file
                          'template_xml|t=s', 
			  'template_xml_conf_key|k=s',
                          'output_xml|o=s', 
			  'component_conf|c=s',
			  #opts for replacing keys in an iterator list
                          'output_dir|o=s',
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
        
        my $outputfile = $options{'output_dir'}.$iteratorconf->{$iterator_list_key}->[$i].".xml.gz";
        push @iterator_output_list,$outputfile;
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
                     <fileName>$file</fileName>
                     <dceSpec type=\"sge\">
                        <OS>linux</OS>\n";
        
        if ( $cfg->val( 'project', '$;PROJECT_CODE$;' ) ) {
            print FILE "                        <group>", $cfg->val( 'project', '$;PROJECT_CODE$;' ), "</group>\n"; 
        }
        
        print FILE " </dceSpec>
                    </commandSet>\n";
    }
    print FILE "</commandSet>
               </commandSetRoot>\n";
    close FILE; 
}
else{
    if($options{'template_xml_conf_key'}){
	my $templatexml = $cfg->val('component','$;'.$options{'template_xml_conf_key'}.'$;');
	
	if(!$templatexml){
	    $logger->logdie("Can't find key in component config $options{'component_conf'} [component] $options{'template_xml_conf_key'}");     
	}	    
	if(! -e $templatexml){
	    $logger->logdie("Can't find referenced xml file ($templatexml) in component config $options{'component_conf'} [component] $options{'template_xml_conf_key'}");     
	
	}
	&replacekeys($replacevals,$templatexml,$options{'output_xml'});
    }
    else{
	if(! -e $options{'template_xml'}){
	    $logger->debug("Can't find template xml file ($options{'template_xml'})") if($logger->is_debug());     
	}
	&replacekeys($replacevals,$options{'template_xml'},$options{'output_xml'});
    }
}

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
#Process LibXML $doc and remove <command> elements
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

    my @inputxml;
    #keeping old code that parses XML in a simplistic way... but works
    #TODO refactor to use LibXML
    while( my $line = <INPUTFILE> ){
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
        my $new_line = "        <dceSpec type=\"sge\">\n            <group>$group</group>\n        </dceSpec>\n$line";
        $line = $new_line;
        $run_dist_cmd_flag = 0;
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
	if($cfg->setval($section,'$;'.$key.'$;',$value)){
	}
	else{
	    $cfg->newval($section,'$;'.$key.'$;',$value);
	}
    }

}

sub import_xml{
    my($file,$keys,$subs,$outputh) = @_;
    my @keys = split(/,/,$keys);
    foreach my $k (@keys){
	my($tok,$val) = split(/=/,$k);
	if(exists $subs->{'$;'.$val.'$;'}){
	    $subs->{$tok} = $subs->{'$;'.$val.'$;'};
	}
	else{
	    $subs->{$tok} = $val;
	}
    }
    open FILE, "$file" or $logger->logdie("Can't open file $file");
    my $run_dist_cmd_flag = 0;
    my @inputxml;
    while(my $line=<FILE>){
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
	if($run_dist_cmd_flag && $line =~ /<dceSpec/) {
	    my $group = &replaceval('$;PROJECT_CODE$;',$subs);
	    $line .= "            <group>$group</group>\n";
	    $run_dist_cmd_flag = 0;
	}
	if($run_dist_cmd_flag && $line =~ /\<\/command\>/ ) {
	    my $group = &replaceval('$;PROJECT_CODE$;',$subs);
	    my $new_line = "        <dceSpec type=\"sge\">\n            <group>$group</group>\n        </dceSpec>\n$line";
	    $line = $new_line;
	    $run_dist_cmd_flag = 0;
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
    close FILE;
}

1;

