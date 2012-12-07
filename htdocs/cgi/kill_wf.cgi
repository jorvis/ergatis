#!/usr/bin/perl

use strict;
use CGI qw(:standard);
use Ergatis::Common;
use Ergatis::ConfigFile;

my $q = new CGI;
my $instancexml = param('instancexml');

my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

if(! -d $ergatis_cfg->val('paths','workflow_run_dir')){
    die "Invalid workflow_run_dir in ergatis.ini : " . $ergatis_cfg->val('paths','workflow_run_dir');
}

my $username = user_logged_in($ergatis_cfg);
my $auth_method = $ergatis_cfg->val('authentication', 'authentication_method');
unless ($auth_method eq 'open' || defined($username)) {
    print $q->header( -type => 'text/html' );
    print_error_page( ergatis_cfg => $ergatis_cfg,
                      message => "You must be logged in to kill pipelines.",
                      links => [
                                 { label => "return to pipeline", is_last => 1, url => "./view_pipeline.cgi?instance=$instancexml" },
                               ]
                    );
    exit(0);
}

my $workflow_root = $ergatis_cfg->val('paths', 'workflow_root');

## these WF_ definitions are usually kept in the $workflow_root/exec.tcsh file,
#   which we're not executing.
$ENV{WF_ROOT} = $workflow_root;
$ENV{WF_ROOT_INSTALL} = $workflow_root;
$ENV{WF_TEMPLATE} = "$workflow_root/templates";
$ENV{PATH} = "$ENV{WF_ROOT}:$ENV{WF_ROOT}/bin:$ENV{WF_ROOT}/add-ons/bin:$ENV{PATH}";

## weak attempt to make this a little safer
$instancexml =~ s|[^a-z0-9\.\_\-\/]||ig;

`$ENV{WF_ROOT}/KillWorkflow -i $instancexml &`;

print redirect(-uri=>"./view_pipeline.cgi?instance=$instancexml");
exit(0);

