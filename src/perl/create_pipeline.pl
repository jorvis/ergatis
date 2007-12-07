#!/usr/bin/perl

use lib (@INC,$ENV{"PERL_MOD_DIR"});
use strict;
use Getopt::Long;
use Ergatis::SavedPipeline;
use Ergatis::Pipeline;
use Ergatis::Logger;

my %options;
my $results = GetOptions (\%options, 
			  'template|s=s',
			  'repository_root|r=s',
                          'log=s',
                          'debug=s', 
                          'help|h' );

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

if(!-e $options{'template'}){
    $logger->logdie("Can't read template XML file: $options{'template'}");
}
if(!-d $options{'repository_root'}){
    $logger->logdie("Can't read repository root: $options{'repository_root'}");
}

my $pipe = Ergatis::SavedPipeline->new( 
					template => $options{'template'},
					pipeline_token => "none"
					);
my $pipeline = $pipe->write_pipeline( repository_root => $options{'repository_root'},
				      id_repository => '/usr/local/devel/ANNOTATION/jorvis/global_id_repository/');
print "Created pipeline ",$pipeline->path(),"\n";
print "Invoke with RunWorkflow -i ",$pipeline->path(),"\n";
