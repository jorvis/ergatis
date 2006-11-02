#!/usr/local/bin/perl

use lib (@INC,$ENV{"PERL_MOD_DIR"});
use Ergatis::SavedPipeline;
use Ergatis::Logger;

my $logger = new Ergatis::Logger('LOG_FILE'=>"/tmp/run_pipeline.log",
				 'LOG_LEVEL'=>5);
$logger = $logger->get_logger();
$logger->debug("Test");
my $pipe = Ergatis::SavedPipeline->new( 
					template => $ARGV[0]
					);
my $pipeline = $pipe->write_pipeline( repository_root => $ARGV[1]);
print "Created pipeline $pipeline\n";
