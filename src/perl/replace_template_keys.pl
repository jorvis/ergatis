#!/usr/local/bin/perl

use lib (@INC,$ENV{"PERL_MOD_DIR"});
use strict;
use Config::IniFiles;
use Ergatis::Logger;
use Getopt::Long;

#--template_xml   Input template XML file for a component 
#--iterator_list  Iterator list.  If iterator list specified then must specify iterator_output_dir
#--component_conf Final component config file
#--output_xml     Output XML

my %options;
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
        
        my $outputfile = $options{'output_dir'}.$iteratorconf->{$iterator_list_key}->[$i].".xml";
        push @iterator_output_list,$outputfile;
        &replacekeys(\%ireplacevals,$options{'template_xml'},$outputfile);
    }
    
    my $distrib;
    my $maxpar;
    my $cstype;
    if ($options{'distribopts'}) {
        if ($options{'distribopts'} =~ /nodistrib=1/i) {
            $distrib = "parallel";
            $maxpar = "<maxParallelCmds>1</maxParallelCmds>";
            $cstype = "serial";
        } else {
            $distrib = "parallel";
            $cstype = "remote-serial";
        }
    }
    
    open FILE, ">$options{'output_xml'}" or $logger->logdie("Can't open file $options{'output_xml'}");
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
                        <OS>linux</OS>
                     </dceSpec>
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
	    $subs->{$parameter} = $replacevalue;
	}
    }
    return $subs;
}

sub replacekeys{
    my($subs,$inputfile,$outputfile) = @_;

    $logger->debug("Reading template $inputfile and writing output $outputfile");
    open( INPUTFILE, "$inputfile" ) or $logger->logdie("Could not open input template file $inputfile");
    open( OUTPUTFILE, "+>$outputfile") or $logger->logdie("Could not open output template file $outputfile");
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
	if($line =~ /\<INCLUDE/){
	    if($line !~ />/){
		$logger->logdie("<INCLUDE> directive must be contained on a single line");
	    }
	    my($file) = ($line =~ /file="([\/\w\.\_]+)"/);
	    my($keys) = ($line =~ /keys="(.*)"/);
	    $logger->debug("Include file $file using keys $keys: $line") if($logger->is_debug());
	    &import_xml($file,$keys,$subs,*OUTPUTFILE);
	    $line = '';
	}
	print OUTPUTFILE $line;
    }
    close INPUTFILE;
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
    while(my $line=<FILE>){
	chomp $line;
	if($linenum==0){
	    @keys = split(/,/,$line);
	    foreach my $key (@keys){
		$iteratorconf->{$key} = [];
	    }
	}
	else{
	    my @elts = split(/,/,$line);
	    for(my $i=0;$i<(@elts);$i++){
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
    my($file,$keys,$subs,$outfh) = @_;
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
    while(my $line=<FILE>){
	if( $line !~ /^\;/ && $line !~ /^\#/ ) {
	    if($line =~ /\<INCLUDE/){
		$line =~ s/(\$;[\w_]+\$;)([^=])/&replaceval($1,$subs).$2/ge;
	    }
	    else{
		$line =~ s/(\$;[\w_]+\$;)/&replaceval($1,$subs)/ge;	    
	    }
        }
	if($line =~ /\<INCLUDE/){
	    if($line !~ />/){
		$logger->logdie("<INCLUDE> directive must be contained on a single line");
	    }
	    my($file) = ($line =~ /file="([\/\w\.\_]+)"/);
	    my($keys) = ($line =~ /keys="(.*)"/);
	    $logger->debug("Include file $file using keys $keys: $line") if($logger->is_debug());
	    &import_xml($file,$keys,$subs,*OUTPUTFILE);
	    $line = '';
	}
	print $outfh $line;
    }
    close FILE;
}

1;

