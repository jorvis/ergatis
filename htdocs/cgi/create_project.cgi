#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::ConfigFile;
use File::Copy;
use File::Temp qw/ tempfile /;
use File::Spec;
use HTML::Template;
use JSON;

umask(0000);

my $q = new CGI;
print $q->header( -type => 'text/html' );

my $tmpl = HTML::Template->new( filename => 'templates/create_project.tmpl',
                                die_on_bad_params => 0,
                              );

my $repository_root = $q->param('repository_root') || '';
my $user_attempted = $q->param('user_attempted') || 0;
my $group_id = $q->param('group_id') || '';

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

my $username = user_logged_in($ergatis_cfg);
my $auth_method = $ergatis_cfg->val('authentication', 'authentication_method');
unless ($auth_method eq 'open' || defined($username)) {
    print_error_page( ergatis_cfg => $ergatis_cfg,
                      message => "You must be logged in to create projects",
                      links => []
                    );
    exit(0);
}

## message to prompt the user to enter a project root
my $enter_msg = 'enter project root directory';

## don't just add a step here unless the conditionals below are updated accordingly.
my $steps = [
    { description => 'enter a project root', complete => 0, failed => 0, 
      failed_msg => 'please enter a project root above',
    },
    { description => 'check project root existence', complete => 0, failed => 0, 
      failed_msg => 'the project root you entered does not exist',
    },
    { description => 'check project root permissions', complete => 0, failed => 0, 
      failed_msg => 'project root does not have write permissions',
    },
    { description => 'create workflow directory', complete => 0, failed => 0, 
      failed_msg => '',
    },
    { description => 'create runtime directory', complete => 0, failed => 0, 
      failed_msg => '',
    },
    { description => 'create pipeline directory', complete => 0, failed => 0, 
      failed_msg => '',
    },
    { description => 'create lock files directory', complete => 0, failed => 0, 
      failed_msg => '',
    },
    { description => 'create project id repository', complete => 0, failed => 0, 
      failed_msg => '',
    },
    { description => 'validate project id repository', complete => 0, failed => 0, 
      failed_msg => '',
    },
    { description => 'create output directory', complete => 0, failed => 0, 
      failed_msg => '',
    },
    { description => 'create default project configuration', complete => 0, failed => 0, 
      failed_msg => '',
    },
    { description => 'create project templates directory', complete => 0, failed => 0, 
      failed_msg => '',
    },
];

my $project_config = "$repository_root/workflow/project.config";
my $default_config = $ergatis_cfg->val('paths', 'default_project_conf');

## if the project root is one of the defaults prefixes, the user hasn't entered it yet.
if (! check_enter_project_root($repository_root, $ergatis_cfg, $enter_msg) ) {
    if ( $user_attempted ) {
        $steps->[0]->{failed} = 1;
    }
    
    print_template();
}

$steps->[0]->{complete} = 1;

my $run_dir = $ergatis_cfg->val('paths', 'workflow_run_dir') || croak "workflow_run_dir not found in ergatis.ini";
my $sudo_scripts_dir = "$run_dir/scripts";

my $current_user = user_logged_in($ergatis_cfg) || undef;
my $run_string = "";

## If our directory to execute sudo scripts from doesn't exist it should be created.
## This should be created by whatever user executes ergatis CGI scripts
if ( -d $run_dir ) {
    if ( ! -d $sudo_scripts_dir ) {
        ( mkdir -p $sudo_scripts_dir ) || croak "Failed to create sudo scripts directory: $sudo_scripts_dir : $!";
    }
} else {
        croak "Invalid workflow_run_dir (doesn't exist) in ergatis.ini: $run_dir";
}

## We need the directory we are currently executing in (CGI scripts directory)
## to create commands in the shell script we will execute sudo'd as the user
my $create_form_script = File::Spec->rel2abs( __FILE__ );
my $cgi_dir = dirname($create_form_script);

$run_string = "sudo -u $current_user " if ($current_user);
$run_string .= "$cgi_dir/bin/create_project.pl --repository_root=$repository_root --default_project_conf=$default_config";
$run_string .= " --group_id=$group_id" if ($group_id);

my ($fh, $create_proj_script) = tempfile("create_project_XXXX", DIR => $sudo_scripts_dir, SUFFIX => '.sh', UNLINK => 1);
print $fh '#!/bin/bash', "\n\n";
print $fh "$run_string";
$fh->flush();
$fh->close();

chmod 0775, $create_proj_script;

my $output = `$create_proj_script`;
my $rc = ($? >> 8);

if ($rc > 0) {
    croak("Unable to execute create_project.pl script");
}

$steps = parse_create_json($steps, $output);
print_template(1);

#############################################################
#                     SUBROUTINES                           #
#############################################################

###
# Checks if our repository root passed in is one of the defaults or user 
# supplied. If the repository root is a default the page is loaded fresh.
###
sub check_enter_project_root {
    my ($repository_root, $ergatis, $enter_msg) = @_;
    return $repository_root && 
     $repository_root ne $ergatis_cfg->val('paths', 'default_project_root') &&
     $repository_root ne $enter_msg;
}

###
# Prints the template for the create_project_form page
### 
sub print_template {
    my $all_good = shift || 0;

    
    $tmpl->param( STEPS               => $steps );
    $tmpl->param( ENTER_MSG           => $enter_msg );
    $tmpl->param( USER_ATTEMPTED      => $user_attempted );
    $tmpl->param( ALL_GOOD            => $all_good );
    $tmpl->param( PROJECT_CONFIG      => $project_config );
    $tmpl->param( REPOSITORY_ROOT     => $q->param('repository_root') || 
                                         $ergatis_cfg->val('paths', 'default_project_root') || 
                                         $enter_msg );
    $tmpl->param( QUICK_LINKS         => &get_quick_links($ergatis_cfg) );
    $tmpl->param( SUBMENU_LINKS       => [
                                            #{ label => 'create project', is_last => 1, url => './create_project.cgi' },
                                         ] );

    print $tmpl->output;
    exit(0);
}

###
# Parses the returned JSON from the create_project.pl script
###
sub parse_create_json {
    my ($steps, $output) = @_;
    my $decoded_json = decode_json($output);
    
    ## The decoded JSON should return a success or failure, and upon a failure
    ## should give us what step number failed so that we can mark it down as 
    ## a failure as well as why it failed.
    if ($decoded_json->{'success'}) {
        map { $_->{complete} = 1 } @$steps;
    } else {
        my $failure_step = int($decoded_json->{'failure_step'});

        ## Need to mark all steps up to the failure point as complete
        for (my $i = 0; $i < $failure_step; $i++) {
            $steps->[$i]->{complete} = 1;
        }

        $steps->[$failure_step]->{failed} = 1;
        $steps->[$failure_step]->{failed_msg} = ($decoded_json->{'message'} || '');
    }
    
    return $steps;
}
