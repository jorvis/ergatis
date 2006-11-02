#!/usr/local/bin/perl

use strict;
use lib (@INC,$ENV{"PERL_MOD_DIR"});
use Ergatis::Logger;
use Getopt::Long;

my %options;
my $results = GetOptions (\%options, 
			  'input_xml|i=s',
                          'output_xml|o=s', 
			  'nodistrib|n=s',
                          'runcmd|r=s', 
                          'log=s',
                          'debug=s', 
                          'help|h' );

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();


if($options{'nodistrib'}){
    exit;
}

open FILE, "<$options{'input_xml'}" or $logger->logdie("Can't open input file $options{'input_xml'}");

my $cmdstr;
my @files;
while(my $line=<FILE>){
    if($line =~ /fileName/){
	my($file) = ($line =~ /\<fileName\>(.*)\<\/fileName\>/);
	push @files,$file;
    }
}

close FILE;

foreach my $file (@files){
	$cmdstr .= "
      <command>
      <type>RunDistributedCommand</type>
      <name>run group</name>
      <id>\$;GROUP_NUMBER\$;</id>
      <state>incomplete</state>
      <executable>$options{'runcmd'}</executable>
      <param>  
	    <key>--template</key>
	    <value>$file</value>
      </param>
      </command>
      ";
    }

close FILE;

open OUT, ">$options{'output_xml'}" or $logger->logdie("Can't open output file $options{'output_xml'}");
print OUT <<_XML
<commandSetRoot xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
				xsi:schemaLocation='commandSet.xsd'>
    <commandSet type="parallel">
    <state>incomplete</state>
    <name>Groups</name>
    $cmdstr
    </commandSet>
</commandSetRoot>
_XML
    ;
close OUT;
