#!/usr/local/bin/perl

use strict;

use CGI qw(:standard);

my $proj = param('project');

if($proj =~ /\//){
    print redirect(-uri=>"show_pipeline.cgi?&xmltemplate=$proj&glob=pipeline");
}
elsif($proj){
    print redirect(-uri=>'show_pipeline.cgi?&xmltemplate=/usr/local/annotation/'.uc($proj).'/Workflow/pipeline/&glob=pipeline');
}
else{
    print header();
    print "<html><body><form>Project or directory(eg. TRYP or /usr/local/annotation/TRYP/Workflow) <input type=text name=project size=50></form></body></html>";
}
