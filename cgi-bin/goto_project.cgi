#!/usr/local/bin/perl

use strict;

use CGI qw(:standard);

my $proj = param('project');

if($proj){
    print redirect(-uri=>'show_pipeline.cgi?&xmltemplate=/usr/local/annotation/'.uc($proj).'/Workflow/pipeline/&glob=pipeline');
}
else{
    print header();
    print "<html><body><form>Project (eg. TRYP) <input type=text name=project></form></body></html>";
}
