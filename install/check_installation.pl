#!/usr/bin/perl

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
files easily searched via grep.

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
my %results = ( pass => 0, fail => 0 );

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

## temp space should exist, be a directory, and be writeable
&check_temp_space();


## if we got this far, print counts
_log( "\npassed $results{pass}/" . ($results{pass} + $results{fail}) . " tests\n\n" );

exit(0);


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
                'sge_root', 'sge_cell', 'sge_qmaster_port', 'sge_execd_port', 'sge_arch',
            ],
        authentication =>
            [
                'authentication_method', 'kerberos_realm', 'sudo_pipeline_execution',
            ],
        quick_links => [],
        disabled_components => [],
        projects => [],
    };
    
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
                    _log("FAIL - parameter '$param' in section '$section' is missing\n");
                }
            }
            
        } else {
            $results{fail}++;
            _log("FAIL - section missing\n");
        }
    }
}

sub check_existance {
    my $thing = shift;

    ## check that it exists
    _log("checking that $thing exists ... ");
    if ( -e $thing ) {
        $results{pass}++;
        _log("PASS\n");
    } else {
        $results{fail}++;
        _log("FAIL\n");
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

sub check_temp_space {
    &check_required_value( 'paths', 'temp_space' );

    my $temp_space = $cfg->val( 'paths', 'temp_space' );

    ## check that it exists
    &check_existance( $temp_space );
}











