#!/usr/bin/perl

=head1 NAME

example_run_pipeline_from_template.pl - An example script showing Ergatis pipeline instantiation
and execution via the API.

=head1 SYNOPSIS

USAGE: example_run_pipeline_from_template.pl 
            --template_directory=/path/to/some_dir/
            --repository_root=/path/to/project_dir
            --id_repository=/path/to/foo/id_repository
            --ergatis_ini=/path/to/ergatis.ini


=head1 OPTIONS

B<--template_directory,-t>
    The full path to a directory containing a pipeline.layout file and its associated
    INI configs

B<--repository_root,-r>
    Also known as a project area, under this directory we should find 

B<--id_repository,-i>
    Path to some global-level ID repository.  Will only be used to pull a pipeline ID.

B<--ergatis_ini,-e>
    Path to the ergatis.ini file, usually kept in the cgi directory of an interface instance.

B<--log,-l> 
    Log file

B<--help,-h>
    This help message


=head1  DESCRIPTION

It's not uncommon to want to instantiate and launch a pipeline from the command-line rather
than using the web interface.  This script illustrates how to do that in just a few lines of
code (4 lines, really, the rest is common perl script template code).

=head1  INPUT

Described in the options above, but the main input is the directory containing an Ergatis
pipeline template.  Once the template is loaded, you can use the Ergatis::ConfigFiles module
to customize the component configuration files for your application.  An example is shown in the
code.

=head1  OUTPUT

This script will use the pipeline template to instantiate a full pipeline in the passed project
area (repository root) and then execute it.

=head1  CONTACT

    Joshua Orvis
    jorvis@users.sf.net

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

use FindBin;
use lib "$FindBin::Bin/../lib";

use Ergatis::ConfigFile;
use Ergatis::SavedPipeline;


my %options = ();
my $results = GetOptions (\%options, 
                          'template_directory|t=s',
                          'repository_root|r=s',
                          'id_repository|i=s',
                          'ergatis_ini|e=s',
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

##############################
## important bits here
##############################

my $template = Ergatis::SavedPipeline->new( 
               template => "$options{template_directory}/pipeline.layout");

my $pipeline = $template->write_pipeline( repository_root => $options{repository_root}, 
                                          id_repository => $options{id_repository} );

## here you can use Ergatis::ConfigFiles to edit some of the newly-written
##  component configurations before pipeline execution.  One example is shown.
##  naming and path conventions allow you to know where the component file is
#my $xdformat_config = new Ergatis::ConfigFile( 
#                        -file => "$options{repository_root}/workflow/runtime/xdformat/" . $pipeline->id . 
#                                 "_default/xdformat.default.user.config"
#                      );
#   $xdformat_config->setval( 'parameters', '$;PRECISION$;', 5 );
#   $xdformat_config->RewriteConfig();

my $ergatis_cfg = new Ergatis::ConfigFile( -file => $options{ergatis_ini} );

$pipeline->run( ergatis_cfg => $ergatis_cfg );


##############################
##############################


exit(0);


sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
}


sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( template_directory repository_root id_repository ergatis_ini );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
}
