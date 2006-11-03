#!/usr/local/bin/perl

use lib (@INC,$ENV{"PERL_MOD_DIR"});
use Ergatis::SavedPipeline;
use Ergatis::Pipeline;
use Ergatis::Logger;

my $logger = new Ergatis::Logger('LOG_FILE'=>"/tmp/create_pipeline.log",
				 'LOG_LEVEL'=>5);
$logger = $logger->get_logger();
$logger->debug("Test");
my $pipe = Ergatis::SavedPipeline->new( 
					template => $ARGV[0],
					pipeline_token => "create_pipeline_test"
					);
my $pipeline = $pipe->write_pipeline( repository_root => $ARGV[1],
				      id_repository => '/usr/local/devel/ANNOTATION/jorvis/global_id_repository/');
print "Created pipeline ",$pipeline->path(),"\n";
