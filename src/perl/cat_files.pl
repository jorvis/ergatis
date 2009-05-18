#!/usr/bin/perl

#!/usr/bin/perl

use lib (@INC,$ENV{"PERL_MOD_DIR"});
use strict;
use Ergatis::Logger;
use Getopt::Long;

my %options;
my $results = GetOptions (\%options, 
			  'list_file|s=s',
			  #opts for replacing keys in a single template file
			  'extension|e=s',
                          'log=s',
                          'debug=s', 
                          'help|h' );

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

open FILE,"$options{'list_file'}" or $logger->logdie("Can't open file\n");
while(my $line=<FILE>){
    if($line =~ /[\/\.]$options{'extension'}$/){
	chomp $line;
	print "cat $line\n";
    }
}
