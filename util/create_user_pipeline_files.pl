#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Config::IniFiles;
use Pod::Usage;
use File::Basename;
use File::stat;
use Storable;

my %opts = &parse_options();
my $ergatis_cfg = new Config::IniFiles(-file => $opts{'ergatis_ini'}) 
    or die ("Could not open ergatis configuration $opts{'ergatis_ini'}");
my $session_dir = $ergatis_cfg->val('authentication', 'session_db_dir');

foreach my $project_name ($ergatis_cfg->Parameters('projects')) {
    my $project_dir = $ergatis_cfg->val('projects', $project_name);
    my $pipeline_dir = File::Spec->catdir($project_dir, 
                                          "workflow/runtime/pipeline");

    my @xml_files = glob("$pipeline_dir/*/pipeline.xml");

    create_pipelines_files($session_dir, @xml_files);
}

##########################################
#            SUBROUTINES                 #
##########################################

sub create_pipelines_files {
    my ($session_dir, @xml_files) = @_;
    my %user_pipelines_map;

    foreach my $xml_file (@xml_files) {
        my $filestat = stat($xml_file);
        my $pipeline_user = getpwuid($filestat->uid);
        my $pipeline_file = File::Spec->catfile($session_dir,
                                                "$pipeline_user.pipelines");

        next if ($pipeline_user eq "ergatis");

        my $user_pipelines = {};            
        if (exists($user_pipelines_map{$pipeline_user})) {
            $user_pipelines = $user_pipelines_map{$pipeline_user}{'pipelines'};
        } elsif (-e $pipeline_file) {
            $user_pipelines = retrieve($pipeline_file);
        }

        my $pipeline_id = basename(dirname($xml_file));
        $user_pipelines->{$pipeline_id} = 1;
        $user_pipelines_map{$pipeline_user}{'pipelines'} = $user_pipelines;
        $user_pipelines_map{$pipeline_user}{'file'} = $pipeline_file;

    }

    write_storable_files(%user_pipelines_map);
}

sub write_storable_files {
    my %pipelines_map = @_;

    foreach my $user (keys %pipelines_map) {
        my $pipelines = $pipelines_map{$user}{'pipelines'};
        my $pipelines_file = $pipelines_map{$user}{'file'};

        store($pipelines, $pipelines_file);            
    }
}

sub parse_options {
    my %options = ();
    my $results = GetOptions (\%options,
                              'ergatis_ini|c=s',
                              'help|h') || pod2usage();

    ## display documentation
    if( $options{'help'} ){
        pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
    }

    unless (exists($options{'ergatis_ini'})) {
        die ("Must provide path to the ergatis ini file");
    }

    return %options;
}
