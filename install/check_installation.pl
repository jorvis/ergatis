#!/usr/local/bin/perl

=head1 NAME

check_installation.pl - This script is meant to perform an exhaustive check of any particular
Ergatis installation including web and project areas, workflow linking and pipeline execution.

=head1 SYNOPSIS

USAGE: check_installation.pl 
            --ergatis_ini=/path/to/ergatis.ini
          [ --log=/path/to/some.log ]

=head1 OPTIONS

B<--ergatis_ini>
    The full path to the ergatis.ini file, usually within the CGI area.  You should have
    customized this file for your installation BEFORE running this validation script.

B<--log,-l> 
    Log file.  If not specified messages will be printed to STDOUT

B<--help,-h>
    This help message

=head1  DESCRIPTION

Ergatis operates based on the values within very customizable configuration files.  It's very
possible (easy) to mess some of these settings up, so this script was written in an attempt
to catch as many of these errors as possible.  It reports both successes and failures, reporting
output like this:

    checking workflow directory path ... PASS
    checking temp space ... FAIL - directory doesn't exist

Counts of each type will be summarized at the end, but case in the warnings makes the output
files easily searched via grep.  There are also SKIPPED and WARNING types.  Warnings are things
that are probably OK but should be verified.  FATAL messages should be rare but are reserved for
those that prevent the script from performing all checks.  Explanations are provided on each type
for guidance.

You should run this script as a user that runs pipelines, else some of the permission checks
may not report correctly.

This script assumes you've run the installer so it doesn't perform all the same checks (such
as the existence of required perl modules)

=head1  INPUT

The main input for this script is the ergatis.ini file, currently found under the Ergatis
directory within your cgi-bin.

=head1  OUTPUT

If the --log option is used this script will write all output to that file, else it will
be written to STDOUT.  Any temporary files created by the script for testing will 
reside in /tmp/ergatis_check

=head1  CONTACT

    Joshua Orvis
    jorvis@users.sf.net

=cut

use strict;
use Config::IniFiles;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
                          'ergatis_ini=s',
                          'log|l=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">$options{log}") || die "can't create log file: $!";
}

## holds the test counts ( pass|fail )
my %results = ( pass => 0, fail => 0, warn => 0 );

my $cfg;

_log("checking ability to parse ergatis.ini ... ");
eval { $cfg = new Config::IniFiles( -file => $options{ergatis_ini} ) };
if ( $@ ) {
    $results{fail}++;
    _logdie("FAIL\n");
} else {
    $results{pass}++;
    _log("PASS\n");
}

## make sure the ini file has all the sections / parameters expected
&check_ergatis_ini_sections();

## each of the following, single-spaced checks map directly to a variable within
##  the ini file, usually checking basic properties.  The subroutines are all named
##  like $section__$parameter.  For example. the 'temp_space' parameter under the 
##  paths section is checked with the subroutine paths__temp_space();
&paths__temp_space();
&paths__pipeline_archive_root();
&paths__workflow_run_dir();
&paths__workflow_log4j();
&paths__workflow_root();
&paths__global_id_repository();
&paths__global_saved_templates();
&paths__pipeline_build_area();
&paths__default_ergatis_dir();
&paths__default_project_root();
&paths__default_project_conf();

&workflow_settings__submit_pipelines_as_jobs();
&workflow_settings__marshal_interval();
&workflow_settings__init_heap();
&workflow_settings__max_heap();
## if we got this far, print results
_log("\n--------------------------------------------------------\n" );
_log( "\npassed $results{pass}/" . ($results{pass} + $results{fail}) . " tests\n" );

if ( $results{warn} == 1 ) {
    _log( "$results{warn} warning\n\n" );
} else {
    _log( "$results{warn} warnings\n\n" );
}

exit(0);

sub _explain {
    my $msg = shift;
    _log("\tEXPLANATION: $msg\n\n");
}


sub _log {
    my $msg = shift;

    if ( $logfh ) {
        print $logfh "$msg";
    } else {
        print "$msg";
    }
}

sub _logdie {
    my $msg = shift;
    
    _log( $msg );
    
    ## print counts
    _log( "\npassed $results{pass}/" . ($results{pass} + $results{fail}) . " tests\n\n" );
    
    exit(1);
}

sub check_ergatis_ini_sections {
    ## here's everything we expect to find
    my $expected = {
        paths => 
            [
                'temp_space', 'pipeline_archive_root', 'workflow_run_dir', 'workflow_log4j', 
                'workflow_root', 'global_id_repository', 'global_saved_templates',
                'pipeline_build_area', 'default_ergatis_dir', 'default_project_root',
                'default_project_conf',
            ],
        workflow_settings => 
            [
                'submit_pipelines_as_jobs', 'marshal_interval', 'init_heap', 'max_heap',
            ],
        display_settings =>
            [
                'pipeline_list_cache_time', 'active_pipeline_age', 'enable_quota_lookup',
                'display_codebase', 'builder_animations',
            ],
        grid =>
            [
                'grid_enabled', 'sge_root', 'sge_cell', 'sge_qmaster_port', 'sge_execd_port', 
                'sge_arch',
            ],
        authentication =>
            [
                'authentication_method', 'kerberos_realm', 'sudo_pipeline_execution',
            ],
        quick_links => [],
        disabled_components => [],
        projects => [],
    };
    
    ## any of these missing variables are fatal.  They need to exist for the other checks to
    #   be performed.
    my $fatal_count = 0;
    
    ## we're only checking that these exist here, no values are examined
    for my $section ( keys %$expected ) {
        _log("checking for section $section in ergatis.ini ... ");
        
        if ( $cfg->SectionExists($section) ) {
            $results{pass}++;
            _log("PASS\n");
            
            
            ## now check each parameter in this section
            for my $param ( @{$$expected{$section}} ) {
                _log("\tchecking for parameter $param ... ");
                
                if ( defined $cfg->val( $section, $param ) ) {
                    $results{pass}++;
                    _log("PASS\n");
                } else {
                    $results{fail}++;
                    $fatal_count++;
                    _log("FAIL - parameter '$param' in section '$section' is missing\n");
                }
            }
            
        } else {
            $results{fail}++;
            $fatal_count++;
            _log("FAIL - section missing\n");
        }
    }
    
    ## if any of these were missing report FATAL
    if ( $fatal_count ) {
        _logdie("FATAL - there were missing sections or parameters in the ergatis.ini file (listed " . 
                "above.  This may be because you are using an old ini file that doesn't have all the " .
                "current parameters.  The script cannot continue unless these exist.\n"
               );
    }
}

sub check_executable {
    my $thing = shift;

    ## check that it executable
    _log("checking that $thing is executable ... ");
    if ( -x $thing ) {
        $results{pass}++;
        _log("PASS\n");
        return 1;
    } else {
        $results{fail}++;
        _log("FAIL\n");
        return 0;
    }
}

sub check_existence {
    my $thing = shift;

    ## check that it exists
    _log("checking that $thing exists ... ");
    if ( -e $thing ) {
        $results{pass}++;
        _log("PASS\n");
        return 1;
    } else {
        $results{fail}++;
        _log("FAIL\n");
        return 0;
    }
}

sub check_is_directory {
    my $thing = shift;

    ## check that it is a directory
    _log("checking that $thing is a directory ... ");
    if ( -d $thing ) {
        $results{pass}++;
        _log("PASS\n");
        return 1;
    } else {
        $results{fail}++;
        _log("FAIL\n");
        return 0;
    }
}

sub check_is_file {
    my $thing = shift;

    ## check that it is a file
    _log("checking that $thing is a file ... ");
    if ( -f $thing ) {
        $results{pass}++;
        _log("PASS\n");
        return 1;
    } else {
        $results{fail}++;
        _log("FAIL\n");
        return 0;
    }
}

sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( ergatis_ini );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            print STDERR "\n--$option is a required option\n\n";
            exit(1);
        }
    }
    
    ## quick check that the ergatis.ini exists and is readable, else no point continuing
    unless ( -e $$options{ergatis_ini} && -r $$options{ergatis_ini} ) {
        print STDERR "\nFAIL: ergatis.ini passed ($$options{ergatis_ini}) doesn't exist or is not readable.\n\n";
        exit(1);
    }
    
    ## handle some defaults
    #$options{optional_argument2}   = 'foo'  unless ($options{optional_argument2});
}

sub check_required_value {
    my ($section, $param) = @_;
    
    _log("checking that a value is defined for $section/$param ... ");
    if ( defined $cfg->val($section, $param) ) {
        $results{pass}++; 
        _log("PASS\n");
    } else {
        $results{fail}++;
        _log("FAIL\n");
    }
}

sub check_writable {
    my $thing = shift;

    ## check that it is writable
    _log("checking that $thing is writable ... ");
    if ( -w $thing ) {
        $results{pass}++;
        _log("PASS\n");
        return 1;
    } else {
        $results{fail}++;
        _log("FAIL\n");
        return 0;
    }
}

sub get_parent_directory {
    my $full = shift;
    
    ## remove any trailing /
    if ( $full =~ m|(.+)/\s*$| ) {
        $full = $1;
    }

    ## get the parent directory
    $full =~ m|(.*)/|;
    
    return $1;
}

sub paths__default_ergatis_dir {
    _log("\nchecking default ergatis directory\n");
    
    my $ergatis_dir = $cfg->val('paths', 'default_ergatis_dir');
    
    ## this directory should already exist
    &check_existence( $ergatis_dir );
    
    ## we can check for a few things within it to make sure
    &check_existence( "$ergatis_dir/bin" );
    &check_existence( "$ergatis_dir/lib" );
}

sub paths__default_project_conf {
    _log("\nchecking default project conf\n");
    
    my $proj_conf = $cfg->val('paths', 'default_project_conf');
    
    &check_existence( $proj_conf );
}

sub paths__default_project_root {
    _log("\nchecking default project root\n");
    
    my $proj_root = $cfg->val('paths', 'default_project_root');
    
    &check_existence( $proj_root );
    &check_writable( $proj_root );
}

sub paths__global_id_repository {
    my $repos = $cfg->val('paths', 'global_id_repository');
    
    _log("\nchecking ID repository\n");
    
    ## the directory needs to exist and contain the validation file
    my $repos_found = &check_existence( $repos );
    &check_is_directory( $repos );
    
    if ( $repos_found ) {
        my $validation_found = &check_existence( "$repos/valid_id_repository" );
        if ( $validation_found ) {
            my $is_file = &check_is_file("$repos/valid_id_repository");
            if ( ! $is_file ) {
                _log("The valid_id_repository is a file that validates any given directory has an ID " .
                     "repository.  You should create an empty file of that name in $repos");
            }
            
        } else {
            _log("checking that $repos/valid_id_repository is a file ... SKIPPED\n");
        }
        
    } else {
        _log("checking for existence of global id repository validation file ... SKIPPED\n");
    }
}

sub paths__global_saved_templates {
    _log("\nchecking global saved templates area\n");

    my $template_dir = $cfg->val('paths', 'global_saved_templates');

    &check_existence( $template_dir );
    &check_writable( $template_dir );
}

sub paths__pipeline_archive_root {
    my $archive_root = $cfg->val('paths', 'pipeline_archive_root');

    ## if the path provided exists, it should be a directory and writable
    if ( -e $archive_root ) {
        &check_is_directory( $archive_root );
        &check_writable( $archive_root );
    
    ## the path provided doesn't have to exist, but the parent does and must be writable
    } else {

        my $parent = get_parent_directory( $archive_root );
        
        _log("WARNING: paths - pipeline_archive_root defined doesn't exist.  This isn't " .
             "necessarily a problem because the interface will create it automatically " .
             "when it needs to.  Instead, we need to check that the parent exists and is " .
             "writable\n");
        $results{warn}++;
        &check_existence( $parent );
        &check_writable( $parent );
    }
}

sub paths__pipeline_build_area {
    _log("\nchecking pipeline build area\n");

    my $build_area = $cfg->val('paths', 'pipeline_build_area');

    ## if the path provided exists, it should be a directory and writable
    if ( -e $build_area ) {
        &check_is_directory( $build_area );
        &check_writable( $build_area );
    
    ## the path provided doesn't have to exist, but the parent does and must be writable
    } else {

        my $parent = get_parent_directory( $build_area );
        
        _log("WARNING: paths - pipeline_build_area defined doesn't exist.  This isn't " .
             "necessarily a problem because the interface will create it automatically " .
             "when it needs to.  Instead, we need to check that the parent exists and is " .
             "writable\n");
        $results{warn}++;
        &check_existence( $parent );
        &check_writable( $parent );
    }
}

sub paths__temp_space {
    &check_required_value( 'paths', 'temp_space' );

    my $temp_space = $cfg->val( 'paths', 'temp_space' );

    ## check that it exists and is a directory
    &check_existence( $temp_space );
    &check_is_directory( $temp_space );
    
    ## check that it's writable
    &check_writable( $temp_space );
}

sub paths__workflow_log4j {
    my $log4j_file = $cfg->val('paths', 'workflow_log4j');

    ## we should later do verification of the contents of this file, but for now just make
    #   sure it exists
    &check_existence( $log4j_file );
}

sub paths__workflow_root {
    my $wf_root = $cfg->val('paths', 'workflow_root');
    
    ## this needs to exist, be a directory and contain a few executables
    &check_existence( $wf_root );
    &check_is_directory( $wf_root );
    
    _log("\nwill now check a few Workflow binaries\n");
    &check_existence( "$wf_root/RunWorkflow" );
    &check_executable( "$wf_root/RunWorkflow" );
    
    &check_existence( "$wf_root/KillWorkflow" );
    &check_executable( "$wf_root/KillWorkflow" );
}

sub paths__workflow_run_dir {
    my $run_dir = $cfg->val('paths', 'workflow_run_dir');

    ## this must exist and be writable
    &check_is_directory( $run_dir );
    &check_writable( $run_dir );
}

sub workflow_settings__init_heap {
    _log("checking workflow_settings:init_heap ... ");
    
    my $init_heap = $cfg->val('workflow_settings', 'init_heap');
    
    ## value should be like 100m
    if ( $init_heap =~ /^\d+m/i ) {
        $results{pass}++;
        _log("PASS\n");
    } else {
        $results{fail}++;
        _log("FAIL\n");
        _explain("workflow_settings:init_heap expected value like 100m");
    }
}

sub workflow_settings__max_heap {
    _log("checking workflow_settings:max_heap ... ");
    
    my $max_heap = $cfg->val('workflow_settings', 'max_heap');
    
    ## value should be like 1024m
    if ( $max_heap =~ /^\d+m/i ) {
        $results{pass}++;
        _log("PASS\n");
    } else {
        $results{fail}++;
        _log("FAIL\n");
        _explain("workflow_settings:max_heap expected value like 1024m");
    }
}

sub workflow_settings__marshal_interval {
    _log("checking workflow_settings:marshal_interval ... ");
    
    my $marshal_interval = $cfg->val('workflow_settings', 'marshal_interval');
    
    ## must be a value greater than 0
    if ( $marshal_interval =~ /^\d+$/ && $marshal_interval > 0 ) {
        $results{pass}++;
        _log("PASS\n");
    } else {
        $results{fail}++;
        _log("FAIL\n");
        _explain("workflow_settings:marshal_interval must be an integer, 0 or greater");
    }
}

sub workflow_settings__submit_pipelines_as_jobs {
    _log("\nchecking workflow_settings:submit_pipelines_as_jobs value ... ");
    my $submit_as_jobs = $cfg->val('workflow_settings', 'submit_pipelines_as_jobs');
    
    ## the only acceptable values are 0 and 1
    if ( $submit_as_jobs !~ /^\d+$/ || ($submit_as_jobs != 0 && $submit_as_jobs != 1) ) {
        $results{fail}++;
        _log("FAIL\n");
        _explain("submit_pipelines_as_jobs value must be either 0 or 1");
        
    } else {
        $results{pass}++;
        _log("PASS\n");
    }
    
    ## if this is set to 1, grid_enabled needs to also be 1
    _log("checking workflow_settings:submit_pipelines_as_jobs grid:grid_enabled agreement ... ");
    if ( $submit_as_jobs == 1 ) {
        if ( $cfg->val('grid', 'grid_enabled') == 1 ) {
            $results{pass}++;
            _log("PASS\n");
        } else {
            $results{fail}++;
            _log("FAIL\n");
            _explain("The workflow_settings:submit_pipelines_as_jobs parameter is set to 1 but " . 
                     "grid:grid_enabled is not 1.  Can't submit jobs without an active grid.");
        }
    } else {
        $results{pass}++;
        _log("PASS\n");
    }
}

















