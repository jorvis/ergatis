#!/usr/bin/perl -w

=head1  NAME 

create_project.pl - Creates the project directory structure and places all 
                    required files in the new directories.

=head1 SYNOPSIS

USAGE: create_project.pl
        --repository_root=/usr/local/annotation/AA1
        --default_project_conf=/path/to/default/conf
       [--restrict_permissions=<true or false>
        --log=/path/to/some.log
       ]

=head1 OPTIONS

B<--repository_root,-r> 
    The project directory, just under which we should find the Workflow directory.

B<--default_project_conf, -d>
    The default project configuration defined in the ergatis.ini file

B<--restrict_permissions, -p>
    0 or 1; by default all folders are created with fully open permissions
    but setting this flag to true will restrict permissions to the folder to only
    allow the creator to write and access the project directories.

B<--log,-l> 
    optional.  will create a log file with summaries of all actions performed.

B<--help,-h> 
    This help message/documentation.

=head1   DESCRIPTION

This script is used to archive a pipeline and, optionally, all its associated output.  
This is based on pipeline_id and repository root.  The contents of the following
folders will be compressed in place (in order):

    $repository_root/workflow/runtime/*/$pipelineid_*
    $repository_root/output_repository/*/$pipelineid_*
    $repository_root/workflow/runtime/pipeline/$pipelineid


=head1 INPUT

The input is defined with the --repository_root parameter defining where 
the new project directories will be created.

=head1 OUTPUT

N/A

=head1 CONTACT

    Cesar Arze
    carze@som.umaryland.edu        

=cut

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use File::Copy;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use HTML::Template;
use JSON;
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
              'repository_root|r=s',
              'default_project_conf|d=s',
              'restrict_permissions|p=i',
              'log|l=s',
              'help|h') || pod2usage();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

$options{'restrict_permissions'} ? umask(027) : umask(000);

my $repository_root = $options{'repository_root'};

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

## check that the project root entered exists
if (! -e $repository_root ) {
    print_json_and_exit(1, "Project root does not exist");
}

## check that the project root is writeable
if (! -w $repository_root ) {
    print_json_and_exit(2, "Project root is not writeable");
}

create_dir_and_record("$repository_root/workflow", 3);
create_dir_and_record("$repository_root/workflow/runtime", 4);
create_dir_and_record("$repository_root/workflow/runtime/pipeline", 5);
create_dir_and_record("$repository_root/workflow/lock_files", 6);
create_dir_and_record("$repository_root/workflow/project_id_repository", 7);

## the project id repository needs a file created within it to verify that it's valid.
##  see IdGenerator module for details.
my $validation_file = "$repository_root/workflow/project_id_repository/valid_id_repository";
if ( -f $validation_file ) {
     $$steps[8]{complete} = 1;
} elsif ( open(my $vfh, ">$validation_file") ) {
     $$steps[8]{complete} = 1;
     close $vfh;
} else {
    $$steps[8]{failed_msg} = $!;
    $$steps[8]{failed} = 1;
    print_template();
}

create_dir_and_record("$repository_root/output_repository", 9);

my $default_config = $options{'default_project_conf'};
my $project_config = "$repository_root/workflow/project.config";

if (-e $default_config) {
    if (! copy($default_config, $project_config) ) {
        print_json_and_exit(10, "$|");
    }
} else {
    print_json_and_exit(10, "Could not find default project config file specified in ergatis.ini: $default_config");
}

create_dir_and_record("$repository_root/workflow/project_saved_templates", 11);

## if we get this far all steps were successful.  
print_json_and_exit();

#############################################################                                                                                                                                                                             
#                     SUBROUTINES                           #                                                                                                                                                                             
############################################################# 

###
# Creates a directory passed in and record whether or not the creation failed.
###
sub create_dir_and_record {
    my ($dir, $step_num) = @_;
    
    ## if the directory already exists, just mark success and go on.  else try to make it
    if (! -d $dir ) {
        if (! mkdir("$dir") ) {
            print_json_and_exit($step_num, "$|");
        }
    }
}

sub print_json_and_exit {
    my ($step_num, $msg) = @_;
    $msg = "All steps completed successfully" if (! $msg);
    my $success = $step_num ? JSON::false : JSON::true;

    my $json_dict = {};
    $json_dict->{'success'} = $success;
    $json_dict->{'message'} = $msg;
    $json_dict->{'failure_step'} = $step_num if ($step_num);

    print encode_json($json_dict);
    exit(0);
}

