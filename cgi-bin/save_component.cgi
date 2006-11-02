#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::ConfigFile;

my $q = new CGI;
print $q->header( -type => 'text/plain' );

my $repository_root = $q->param('repository_root') || die "need a repository root";
my $component_name = $q->param('component_name') || die "need a component name";
my $component_id = $q->param('component_id') || die "need a component id";
my $output_token = $q->param($component_id . '_OUTPUT_TOKEN');
my $build_directory = $q->param('build_directory') || die "need a build directory";

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

my $shared_cfg = new Ergatis::ConfigFile( -file => "$repository_root/workflow_config_files/sharedconf.ini" );
my $docs_dir = $shared_cfg->val('init', '$;WORKFLOWDOCS_DIR$;') || die "failed to determine workflowdocs_dir";

## if the build directory doesn't exist yet, create it
if (! -d $build_directory ) {
    mkdir( $build_directory ) || die "can't create build directory $build_directory: $!";
}

## load all of these into a hash so they can be replaced.
my %component_vars;
for ( $q->param ) {
    if ( /${component_id}_(.+)/ ) {
        $component_vars{'$;' . $1 . '$;'} = $q->param($_);
    }
}

## loop through each element of the template and replace with the user's settings.
## open the template for this component
open(my $template_fh, "<$docs_dir/$component_name.config") || die "failed to find template file for $component_name in $docs_dir: $!";

## open the output component file
open(my $user_fh, ">$build_directory/$component_name.$output_token.config") || die "can't create temporary component file; $!";

while (my $line = <$template_fh>) {
    chomp $line;
    
    ## if this isn't a comment and has an X = Y pattern
    if ( $line !~ /^\;/ && $line =~ /^\s*(.+?)\s*\=\s*(.+?)\s*$/ ) {
        my ($key, $val) = ($1, $2);
        
        if ( exists $component_vars{$key} ) {
            print $user_fh "$key = $component_vars{$key}\n";
            
        } else {
            print $user_fh "$line\n";
        }
        
    } else {
        print $user_fh "$line\n";
    }
}


## if we get this far, make sure we exit correctly so the right response code is sent.
exit(0);
