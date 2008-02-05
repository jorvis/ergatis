#!/usr/bin/perl -w

use strict;
#use lib '/usr/local/devel/ANNOTATION/ard/current/lib/5.8.8';
use Ergatis::SavedPipeline;

my $source_pipeline = shift || print_usage();
my $repository_root = shift || print_usage();
my $global_id_repos = shift || print_usage();

my $template = Ergatis::SavedPipeline->new( source => $source_pipeline );

my $temp_area = "/tmp/saved_temp_$$";
#my $temp_area = '/tmp/saved_temp_29820';

## make sure temporary space exists
mkdir("$temp_area");

#print STDERR "writing template into $temp_area\n";
$template->write_template( template => $temp_area );

my $new_pipeline = Ergatis::SavedPipeline->new( template => "$temp_area/pipeline.xml" );
$new_pipeline->write_pipeline( repository_root => $repository_root,
                               id_repository => $global_id_repos );



sub print_usage {
    print STDERR "\n\nUSAGE: $0 source_pipeline.xml repository_root global_id_repository\n\n";
    exit(1);
}
