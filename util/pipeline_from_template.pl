#!/usr/bin/perl

use strict;
#use lib '/usr/local/projects/ergatis/lib/perl5/site_perl/5.8.8';
use Ergatis::SavedPipeline;

my $source_template_dir = shift || print_usage();
my $repository_root = shift || print_usage();
my $id_repository = shift || print_usage();

my $pipeline = Ergatis::SavedPipeline->new( 
    template => "$source_template_dir/pipeline.layout");

$pipeline->write_pipeline( repository_root => $repository_root, id_repository => $id_repository );


sub print_usage {
    print STDERR "\n\nUSAGE: $0 template_directory repository_root id_repository\n\n";
}
