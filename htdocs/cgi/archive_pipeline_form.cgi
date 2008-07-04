#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::ConfigFile;
use HTML::Template;

my $q = new CGI;
print $q->header( -type => 'text/html' );

umask(0000);

my $repository_root = $q->param('repository_root') || die "didn't get a repository_root";
my $pipeline_id = $q->param('pipeline_id') || die "didn't get a pipeline_id";

my $tmpl = HTML::Template->new( filename => 'templates/archive_pipeline_form.tmpl',
                                die_on_bad_params => 1,
                              );

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

my $archive_root = $ergatis_cfg->val( 'paths', 'pipeline_archive_root' ) || die "pipeline_archive_root not defined in ergatis.ini file";


$tmpl->param( ARCHIVE_ROOT => $archive_root );
$tmpl->param( PIPELINE_ID => $pipeline_id );
$tmpl->param( REPOSITORY_ROOT => $repository_root );

$tmpl->param( QUICK_LINKS         => &get_quick_links($ergatis_cfg) );
$tmpl->param( SUBMENU_LINKS       => [
                                        { label => 'pipeline list', is_last => 1, url => "./pipeline_list.cgi?repository_root=$repository_root" },
                                     ] );

print $tmpl->output;


exit(0);
