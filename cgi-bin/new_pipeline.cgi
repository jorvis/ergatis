#!/usr/local/bin/perl

use strict;

use XML::Twig;
use CGI qw(:standard);
use Tree::DAG_Node;
use Config::IniFiles;
use Ergatis::Common;

umask 000;

## called like: ergatis/new_pipeline.cgi?&root=/usr/local/annotation/AA1/Workflow/pipeline
## no trailing slash
my $repository_root = param('repository_root') || die "I need a repository root to create a pipeline\n";;
my $root = "$repository_root/Workflow/pipeline";

## we need access to the IdGenerator module, but which to use depends on which
##  ergatis installation the project uses.  here we look it up in that project's
##  conf file and then require it.
BEGIN {
    my $repository_root = param('repository_root');
    if ($repository_root) {
        my $lib_dir = get_project_conf_param( $repository_root, 'init', '$;LIB_DIR$;' );

        my $mod_loc = get_module_path( $lib_dir, 'Workflow::IdGenerator' );
        require $mod_loc;
    }
}

my $idgen = new Workflow::IdGenerator;
my $id = $idgen->next_id();

## create the directory for the pipeline
if (create_directory("$root/$id")) {
    die "couldn't create directory $root/$id for pipeline\n";
}

my $xmltemplate = "$root/$id/pipeline.xml";

my $newxml = &get_skeleton_commandset();

open FILE,"+>$xmltemplate" or die("Can't open file $xmltemplate");
print FILE $newxml;
close FILE;

open FILE,"+>$xmltemplate.ini" or die("Can't open file $xmltemplate");
print FILE "[start]\n\n";
close FILE;

system("chmod 666 $xmltemplate");

print redirect(-uri=>"./show_pipeline.cgi?xmltemplate=$xmltemplate");

exit;

sub get_skeleton_commandset{
    return '<commandSetRoot xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="commandSet.xsd">
            <commandSet type="serial">
                <configMapId>start</configMapId>
            </commandSet>
            </commandSetRoot>';
}


sub create_directory{
    my($dir) = @_;
    my $ret = system("mkdir -m 777 -p $dir");
    $ret >>= 8;
    return $ret;
}


