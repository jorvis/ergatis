#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::ConfigFile;
use HTML::Template;

my $q = new CGI;
print $q->header( -type => 'text/html' );

my $tmpl = HTML::Template->new( filename => 'templates/build_pipeline.tmpl',
                                die_on_bad_params => 1,
                              );

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );
my $build_area = $ergatis_cfg->val( 'paths', 'pipeline_build_area' ) || die "failed to determine pipeline_build_area";

my $repository_root = $q->param('repository_root') || die "need a repository root";
my $shared_cfg = new Ergatis::ConfigFile( -file => "$repository_root/workflow_config_files/sharedconf.ini" );
my $workflowdocs_dir = $shared_cfg->val( 'init', '$;WORKFLOWDOCS_DIR$;' );

## make sure the build area exists
if (! -d $build_area) {
    mkdir($build_area) || die "failed to make build directory $build_area: $!";
}

## read the available components
my @components;
opendir(my $idh, $workflowdocs_dir) || die "can't read component directory ($workflowdocs_dir): $!";
while ( my $thing = readdir $idh ) {
    if ( $thing =~ /(.+).config/ ) {
        push @components, { name => $1 };
    }
}

## we want the components to be sorted by name
@components = sort {$a->{name} cmp $b->{name}} @components;

$tmpl->param( QUICK_LINKS         => &get_quick_links($ergatis_cfg) );
$tmpl->param( SUBMENU_LINKS       => [
                                        { label => 'run pipeline', is_last => 1, url => 'javascript:checkAndRunPipeline()' },
                                     ] );
$tmpl->param( REPOSITORY_ROOT => $repository_root );
$tmpl->param( COMPONENTS => \@components );
$tmpl->param( BUILD_DIRECTORY => "$build_area/" .temp_pipeline_id() );

print $tmpl->output;


# usage: $string = prettydate( [$time_t] );
# omit parameter for current time/date
sub pretty_date_time {
   @_ = localtime(shift || time);
   return(sprintf("%04d%02d%02d%02d%02d%02d", $_[5]+1900, $_[4]+1, $_[3], @_[2,1], @_[0]));
} 

sub temp_pipeline_id {
    ## currently the date/time with a number between 10 and 100
    return pretty_date_time() . (int(rand(90)) + 10);
}
