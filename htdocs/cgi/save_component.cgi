#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::ConfigFile;

umask(0000);

my $q = new CGI;
print $q->header( -type => 'text/plain' );

my $repository_root = $q->param('repository_root') || die "need a repository root";
my $component_name = $q->param('component_name') || die "need a component name";
my $component_id = $q->param('component_id') || die "need a component id";
my $build_directory = $q->param('build_directory') || die "need a build directory";

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

my $shared_cfg = new Ergatis::ConfigFile( -file => "$repository_root/workflow/project.config" );
my $docs_dir = $shared_cfg->val('project', '$;DOCS_DIR$;') || die "failed to determine docs_dir in ergatis.ini";

## if the build directory doesn't exist yet, create it
if (! -d $build_directory ) {
    mkdir( $build_directory ) || check_mkdir_error("can't create build directory $build_directory:", $!);
}

## make sure the output token doesn't have any leading or trailing whitespace
$q->param($component_id . '_OUTPUT_TOKEN') =~ /^\s*(.+?)\s*$/;
my $output_token = $1;

## load all of these into a hash so they can be replaced.
my %component_vars;
for ( $q->param ) {
    if ( /${component_id}_(.+)/ ) {
        $component_vars{'$;' . $1 . '$;'} = $q->param($_);
    }
}

use Data::Dumper;
print STDERR Dumper(%component_vars) . "\n";

## disable grid submissions if they are disabled in ergatis.ini
if ( $ergatis_cfg->val('grid', 'grid_enabled') == 0 ) {
    ## reset NODISTRIB if it exists
    if ( exists $component_vars{'$;NODISTRIB$;'} ) {
        $component_vars{'$;NODISTRIB$;'} = 1;
    }
}

## loop through each element of the template and replace with the user's settings.
## open the template for this component
open(my $template_fh, "<$docs_dir/$component_name.config") || die "failed to find template file for $component_name in $docs_dir: $!";
my @config_template_lines = <$template_fh>;

if (-e "$docs_dir/shared_parameters.config") {
    open(my $common_template_fh, "<$docs_dir/shared_parameters.config") || die "failed to find shared_parameters.config in $docs_dir: $!";
    my @common_lines = <$common_template_fh>;
    @config_template_lines = (@config_template_lines, @common_lines);
}

## open the output component file
open(my $user_fh, ">$build_directory/$component_name.$output_token.config") || die "can't create temporary component file; $!";

my $seen_params = {};
for my $line (@config_template_lines) {
    chomp $line;
    
    ## if this isn't a comment and has an X = Y pattern
    if ( $line !~ /^\;/ && $line =~ /^\s*(.+?)\s*\=\s*(.*?)\s*$/ ) {
        my ($key, $val) = ($1, $2);
        
        if (! exists($seen_params->{$key}) ) {
            if ( exists $component_vars{$key} ) {
                print $user_fh "$key=$component_vars{$key}\n";
            } else {
                print $user_fh "$line\n";
            }

            $seen_params->{$key} = 1;
        }
    } else {
        print $user_fh "$line\n";
    }
}

## if we get this far, make sure we exit correctly so the right response code is sent.
exit(0);

## since save_component.cgi is called multiple times simultaneously on the interface in the
#   sample pipeline via AJAX calls (such as when cloning) there can be a race condition
#   when creating directories.  If this fails, it checks to see if the failure is due to the
#   directory already existing.  If so, it is ignored, else it dies.
#
#   This check happens via the error message, which is admittedly tied to the OS and even
#   version.  It can be expanded for other platform messages.  There be dragons in the details
#   of trying to do it other ways.
sub check_mkdir_error {
    my ($usr_msg, $os_msg) = @_;
    
    if ( $os_msg =~ /File exists/i ) { 
        print STDERR "WARN: mkdir race condition caught and handled for build directory $build_directory";
    } else {
        die "$usr_msg : $os_msg";
    } 
}
