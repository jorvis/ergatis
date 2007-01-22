#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::ConfigFile;
use File::Copy;
use HTML::Template;

umask(0000);

my $q = new CGI;
print $q->header( -type => 'text/html' );

my $tmpl = HTML::Template->new( filename => 'templates/create_project.tmpl',
                                die_on_bad_params => 0,
                              );

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

## message to prompt the user to enter a project root
my $enter_msg = 'enter project root directory';

my $repository_root = $q->param('repository_root') || '';
my $user_attempted = $q->param('user_attempted') || 0;

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
];

## if the project root is one of the defaults prefixes, the user hasn't entered it yet.
if (! check_enter_project_root() ) {
    if ( $user_attempted ) {
        $$steps[0]{failed} = 1;
    }
    
    print_template();
}
$$steps[0]{complete} = 1;
     
## check that the project root entered exists
if (! -e $repository_root ) {
    $$steps[1]{failed} = 1;
    print_template();
}
$$steps[1]{complete} = 1;

## check that the project root is writeable
if (! -w $repository_root ) {
    $$steps[2]{failed} = 1;
    print_template();
}
$$steps[2]{complete} = 1;

create_dir_and_record("$repository_root/workflow", 3);
create_dir_and_record("$repository_root/workflow/runtime", 4);
create_dir_and_record("$repository_root/workflow/runtime/pipeline", 5);
create_dir_and_record("$repository_root/workflow/lock_files", 6);
create_dir_and_record("$repository_root/workflow/project_id_repository", 7);

## the project id repository needs a file created within it to verify that it's valid.
##  see IdGenerator module for details.
if ( open(my $vfh, ">$repository_root/workflow/project_id_repository/valid_id_repository") ) {
     $$steps[8]{complete} = 1;
     close $vfh;
} else {
    $$steps[8]{failed_msg} = $!;
    $$steps[8]{failed} = 1;
    print_template();
}

create_dir_and_record("$repository_root/output_repository", 9);

my $default_config = $ergatis_cfg->val('paths', 'default_project_conf');
my $project_config = "$repository_root/workflow/project.config";

if (-e $default_config) {
    if (! copy($default_config, $project_config) ) {
        $$steps[10]{failed_msg} = "$!";
        $$steps[10]{failed} = 1;
        print_template();       
    }

    $$steps[10]{complete} = 1;
    
} else {
    $$steps[10]{failed_msg} = "could not find default project config file specified in ergatis.ini : $default_config";
    $$steps[10]{failed} = 1;
    print_template();
}


## if we get this far all steps were successful.  
print_template(1);

exit(0);

sub check_enter_project_root {
    return $repository_root && 
     $repository_root ne $ergatis_cfg->val('paths', 'default_project_root') &&
     $repository_root ne $enter_msg;
}

sub create_dir_and_record {
    my ($dir, $step_num) = @_;
    
    ## if the directory already exists, just mark success and go on.  else try to make it
    if (! -d $dir ) {
        if (! mkdir("$dir") ) {
            $$steps[$step_num]{failed_msg} = $!;
            $$steps[$step_num]{failed} = 1;
            print_template();
        }
    }
    $$steps[$step_num]{complete} = 1;
}

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
