#!/usr/local/bin/perl

use strict;
use CGI qw(:standard);


#
# editor:   sundaram@tigr.org
# date:     2005-08-04
# bgzcase:  2030
# comment:  Display dropdown list to choose projects from
#
my $projects = "/usr/local/devel/ANNOTATION/sundaram/datamanagement/project.dat";
my $prj = {};
my $flag=0;

if ((-e $projects) && (-r $projects) && (!-z $projects)){

    open (INFILE, "<$projects") or die "Could not open file '$projects': $!";

    while (my $line = <INFILE>){
	chomp $line;
	next if ($line =~ /^\s+$/ );
	
	$prj->{$line}++;
	$flag++;
    }
}


my $proj = param('project');

if($proj =~ /\//){

    #
    # Add the project path to the control file if not already present
    #
    if (! exists $prj->{$proj}){
	open (PROJECT, ">>$projects") or die "Could not open project control file '$projects' in append-write mode:$!";
	print PROJECT $proj . "\n";
	close PROJECT;
    }

    print redirect(-uri=>"show_pipeline.cgi?&xmltemplate=$proj&glob=pipeline");
}
elsif($proj){
    print redirect(-uri=>'show_pipeline.cgi?&xmltemplate=/usr/local/annotation/'.uc($proj).'/Workflow/pipeline/&glob=pipeline');
}
else{
    print header();
    print "<html>\n".
    "<head>\n".
    "<title>$0</title>\n".
    "<script language=\"JavaScript\">\n".
    "function submitSelection()\n".
    "{\n".
    " newurl = \"show_pipeline.cgi?&xmltemplate=\" + document.form1.selectp.value + \"/Workflow/pipeline/&glob=pipeline\" \n".
    "   window.location = newurl\n".
    "}\n".
    "</script>\n".
    "</head>\n".
    "<body><form name=\"form1\">Project or directory(eg. TRYP or /usr/local/annotation/TRYP/Workflow) <input type=\"text\" name=\"project\" size=\"50\">\n";


    if ($flag > 0){
	print "<br>\n".
	"Or select your project here\n".
	"<select name=\"selectp\" onchange=\"submitSelection()\">\n";
	
	foreach my $option (sort keys %{$prj} ){
	    print "<option>$option\n";
	}

    }


    print "</form>\n".
    "</body></html>\n";
    
}
