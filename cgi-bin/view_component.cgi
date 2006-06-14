#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Date::Manip;
use File::stat;
use Monitor;
use POSIX;
use XML::Twig;


my $q = new CGI;

print $q->header( -type => 'text/html' );

my $pipeline_xml = $q->param("pipeline_xml") || die "pass pipeline_xml";

## it may have been compressed
if (! -e $pipeline_xml && -e "$pipeline_xml.gz") {
    $pipeline_xml .= '.gz';
}

my $pipeline_xml_fh;
if ($pipeline_xml =~ /\.gz/) {
    open($pipeline_xml_fh, "<:gzip", "$pipeline_xml") || die "can't read $pipeline_xml: $!"; 
} else {
    open($pipeline_xml_fh, "<$pipeline_xml") || die "can't read $pipeline_xml: $!";       
}

my $twig = XML::Twig->new( );
$twig->parse($pipeline_xml_fh);

my $parent_commandset = $twig->root->first_child('commandSet');
my $component_state = $parent_commandset->first_child('state')->text || 'unknown';

my $component_name = 'component';
if ( $parent_commandset->first_child('name')->text =~ /(.+) workflow/) {
    $component_name = $1;
}

my $parent_pipeline = '';
if ( $parent_commandset->first_child('parentFileName') ) {
    $parent_pipeline = $parent_commandset->first_child('parentFileName')->text;
}

my $component_project = '?';
my $repository_root = '';
if ( $parent_pipeline =~ /(.+\/(.+?))\/Workflow/ ) {
    $repository_root = $1;
    $component_project = $2;
}

## quota information (only works if in /usr/local/annotation/SOMETHING)
my $quotastring = &quota_string($repository_root);

my ($starttime, $endtime, $lastmodtime, $state, $runtime) = ('n/a', 'n/a', '', 'unknown', 'n/a');

($starttime, $endtime, $runtime) = &time_info($parent_commandset);

my $filestat = stat($pipeline_xml);
my $user = getpwuid($filestat->uid);
$lastmodtime = time - $filestat->mtime;
$lastmodtime = strftime( "%H hr %M min %S sec", reverse split(/:/, DateCalc("today", ParseDate(DateCalc("now", "- ${lastmodtime} seconds")) ) ));


print_header($parent_pipeline);

## look at each of the children of the root
foreach my $child ( $parent_commandset->children() ) {
    
    ## if this is a command, process its details
    if ($child->gi eq 'command') {
        
        &process_command($twig, $child);
    
    ## if it is a commandSet, it should be a file-based subflow
    } elsif ($child->gi eq 'commandSet') {
    
        ## make sure it has a fileName element
        if ( $child->has_child('fileName') ) {
            print "    <div class='subflow'>\n";
            
            ## customize the label for this
            if ( $child->first_child('name')->text() eq 'Iterated subflow' ) {
                print "        <h1>input analysis groups</h1>\n";
            } else {
                print "        <h1>" . $child->first_child('name')->text() . "</h1>\n";
            }
            
            &parse_groups_xml( $child->first_child('fileName')->text() );
            
            print "    </div>\n";
        }
    }
}

print_footer();

exit(0);



sub parse_groups_xml {
    my $filename = shift;
    
    ## it may have been compressed
    if (! -e $filename) {
        if (-e "$filename.gz") {
            $filename .= '.gz';
        } else {
            print "            <div>not yet created</div>\n";    
        }
    }
    
    ## make sure this is a groups.xml file
    if ( $filename !~ /groups.xml/ ) {
        print "            <div>unable to handle $filename</div>\n";
        return;
    }
    
    my ($component_start_time, $component_end_time, $component_state);
    
    my $filename_fh;
    if ($filename =~ /\.gz/) {
        open($filename_fh, "<:gzip", "$filename") || die "can't read $filename: $!"; 
    } else {
        open($filename_fh, "<$filename") || die "can't read $filename: $!";       
    }
    
    ## create the twig
    my $twig = XML::Twig->new( twig_roots => {
                                    'command' => \&process_subflowgroup,
                               }
                             );
    $twig->parse($filename_fh);
    
}

sub process_subflowgroup {
    my ($twig, $command) = @_;
    my $state = 'unknown';
    my $execution_host = '';
    my $name = '';
    my $file = '';
    my $group_num = '?';
    my $workflow_id = '?';
    my $grid_id = '?';
    
    ## need to parse through the params to get the one that references the
    ##  instance file descriptor
    for my $param ( $command->children('param') ) {
        ## it will have a key element --instance
        if ( $param->first_child('key')->text() eq '--instance' ) {
            $file = $param->first_child('value')->text();
        }
    }
    
    ## pull the name out of the subflow_file:
    if ( $file =~ /subflow(\d+)groups(\d+).xml/ ) {
        $name = "subflow$1groups$2";
        $group_num = $2;
    }
    
    ## get the state, if it has one
    if ( $command->first_child('state') ) {
        $state = $command->first_child('state')->text();
    }
    
    ## grab data from the dceSpec if it has one
    if ( $command->first_child('dceSpec') ) {
        if ( $command->first_child('dceSpec')->first_child('executionHost') ) {
            $execution_host = $command->first_child('dceSpec')->first_child('executionHost')->text();
        }
        
        if ( $command->first_child('dceSpec')->first_child('jobID') ) {
            $grid_id = $command->first_child('dceSpec')->first_child('jobID')->text();
        }
    }
    
    if ( $command->first_child('id') ) {
        $workflow_id = $command->first_child('id')->text();
    }
    
    my ($start_time, $end_time, $runtime) = &time_info($command);

    my $ret_value = 'unknown';
    my $message = '';
    
    ## if there is a status and a message, grab it
    if ( $command->first_child('status') ) {
        if ( $command->first_child('status')->first_child('retValue') ) {
            $ret_value = $command->first_child('status')->first_child('retValue')->text;
        }
        
        if ( $command->first_child('status')->first_child('message') ) {
            $message = $command->first_child('status')->first_child('message')->text;
        }
        
        ## don't include 'command finished' messages
        $message =~ s/command finished//;
        
        if ( $message && $ret_value != 0) {
            $message = <<messageBLOCK;
        <div class='messageblock'>
            return value: $ret_value<br>
            message: $message
        </div>
messageBLOCK
        } else {
            $message = '';
        }
    }
    my $hostsrvstr = join(',',split(/\./,$execution_host));
    print <<SubflowGroupBar;
        <div id='${name}_bar' class='subflowgroupbar'>
            <div class='leftside'>
                <img id='${name}_arrow' class='expander' src='/ergatis/arrow_right.gif' onclick='toggle_subflowgroup_display("$name", "$file");' alt='expand' title='expand'>
                <img class='status' src='/ergatis/status_$state.png' title='$state' alt='$state'>
                <span class='minor'>group $group_num</span>
            </div>
            <div class='rightside'>
                <span class='infolabel' id='${name}_infolabel' onclick='toggle_group_info("$name")'>show group info</span>
            </div>
        </div>
$message
        <div id='${name}_info' class='subflowinfo' style='display: none;'>
            <table>
                <tr><th>workflow id:</th><td>$workflow_id</td></tr>
                <tr><th>state:</th><td>$state</td></tr>
                <tr><th>start time:</th><td>$start_time</td></tr>
                <tr><th>end time:</th><td>$end_time</td></tr>
                <tr><th>duration:</th><td>$runtime</td></tr>
                <tr><th>grid id:</th><td>$grid_id</td></tr>
                <tr><th>execution host:</th><td>$execution_host<a href='http://intranet.tigr.org/cgi-bin/sysadmin/hobbit/bb-hostsvc.sh?HOSTSVC=$hostsrvstr.cpu&IP=0.0.0.0&DISPLAYNAME=$execution_host' target='_blank'>[BB]</a><a href='http://enterprise.tigr.org/ganglia/?m=load_one&r=hour&s=descending&c=Main+Cluster&h=$execution_host&sh=1&hc=4' target='_blank'>[Ganglia]</a></td></tr>
                <tr><th colspan='2'>xml:</th></tr>
                <tr><td colspan='2'><a href='./view_formatted_xml_source.cgi?file=$file'>$file</a></td></tr>
            </table>
        </div>
        <div id='${name}_data' class='subflowgroupdata' style='display: none;'></div>
SubflowGroupBar
}



sub print_header {
    my $parent_pipeline = shift;

    print <<HeAdER;
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
"http://www.w3.org/TR/html4/strict.dtd">

<html>

<head>
    <meta http-equiv="Content-Language" content="en-us">
    <meta http-equiv="Content-Type" content="text/html; charset=windows-1252">
    <title>$component_project|$component_name|$component_state</title>
    <style type="text/css">
        body {
            font-family: verdana, sans-serif;
            font-size: 8pt;
            margin: 0;
            padding: 0;
        }
        img.status {
            margin-right: 5px;
        }
        img.expander, img.reloader {
            cursor: pointer;
        }
        img.reloader {
            margin-left: 5px;
        }
        #workflowcontainer {
            margin: 10px;
        }
        #bannertop {
            height: 36px;
            margin: 0px;
            padding: 0px;
            background-color: rgb(76,109,143);
        }
        #bannerbottom {
            background-color: rgb(229,229,220);
            border-bottom: 1px solid rgb(131,150,145);
        }
        a {
            text-decoration: none;
            color: rgb(0,100,0);
            cursor: pointer;
            cursor: hand;
        }
        #pipeline {
            font-weight: bold;
            margin-left: 5px;
        }
        #pipelinesummary {
            padding: 5px 0px 5px 0px;
        }
        #pipelinecommands {
            margin-left: 15px;
        }

        .pipelinestat {
            margin-left: 15px;
            display: inline;
            color: rgb(75,75,75);
        }
        div.subflow {
            margin-bottom: 10px;
        }
        div h1 {
            font-weight: bold; 
            font-size: 100%; 
        }
        div.subflowgroupbar, div.subflowbar {
            margin-bottom: 5px;
            width: 500px;
            padding-bottom: 2px;
            margin-top: 3px;
        }
        div.subflowgroupbar {
            border-bottom: 1px solid rgb(150,150,150);
        }
        div.command {
            margin-bottom: 2px;
            width: 500px;
            padding-bottom: 2px;
        }
        div.subflowgroupdata, div.subflowdata {
            padding: 5px 0px 0px 30px;
            margin-bottom: 15px;
            display: none;
        }
        div.start, div.end {
            border: 1px solid rgb(150,150,150);
            width: 300px;
            background-color: rgb(220,220,220);
            padding: 10px;
            margin: 10px 0px 10px 0px;
            font-weight: bold;
        }
        div.messageblock {
            color: rgb(100,100,100);
            padding-left: 30px;
        }
        span.minor {
            color: rgb(150,150,150);
            font-size: 80%;
            margin-left: 20px;
        }
        span.infolabel {
            color: rgb(0,0,200);
            font-size: 80%;
            cursor: pointer;
            margin-left: 20px;
        }
        div.leftside {
            float: left;
        }
        div.rightside {
            text-align: right;
        }
        div.subflowinfo, div.cmdinfo {
            padding-left: 25px;
        }
        div.subflowinfo table, div.cmdinfo table {
            background-color: rgb(235,235,235);
            padding: 5px;
            width: 500px;
        }
        div.navigation {
            margin-bottom: 10px;
        }
        td, th {
            font-size: 90%;
            text-align: left;
            padding: 0px;
        }
        th {
            color: rgb(100,100,100);
        }
        .hidden {
            display: none;
        }
        img.navbutton {
            border: none;
            margin: 5px 3px 0px 0px;
            padding: 0px;
        }
        div {
            overflow: hidden;
        }
    </style>
    <script type="text/javascript">
        // the following two functions could probably be combined later
        
        function toggle_subflowgroup_display(subflowname, subflowfile) {
            subflownamedata = get_object(subflowname + "_data");
            
            // is it hidden?
            if ( subflownamedata.style.display == 'none' ) {
                subflownamedata.style.display = 'block';
                get_object(subflowname + '_arrow').src = '/ergatis/arrow_down.gif';
                get_object(subflowname + '_data').innerHTML = 'parsing subflow data';
                sendElementUpdateRequest('./subflowgroup_summary.cgi?xml_input=' + encodeURIComponent(subflowfile) + '&nocache=' + no_cache_string(), updateSubflowGroup, subflowname);
                
            } else {
                subflownamedata.style.display = 'none';
                get_object(subflowname + '_arrow').src = '/ergatis/arrow_right.gif';
                get_object(subflowname + '_data').innerHTML = '';
            }
        }
        
        function toggle_subflow_display(subflowname, subflowfile) {
            subflownamedata = get_object(subflowname + "_data");
            
            // is it hidden?
            if ( subflownamedata.style.display == 'none' ) {
                subflownamedata.style.display = 'block';
                get_object(subflowname + '_arrow').src = '/ergatis/arrow_down.gif';
                get_object(subflowname + '_data').innerHTML = 'parsing subflow data';
                sendElementUpdateRequest('./subflow_summary.cgi?xml_input=' + encodeURIComponent(subflowfile) + '&nocache=' + no_cache_string(), updateSubflow, subflowname);
                
            } else {
                subflownamedata.style.display = 'none';
                get_object(subflowname + '_arrow').src = '/ergatis/arrow_right.gif';
                get_object(subflowname + '_data').innerHTML = '';
            }
        }
        
        // threadsafe asynchronous XMLHTTPRequest code
        function sendElementUpdateRequest(url, callback, elementname) {
            // inner functions below.  using these, reassigning the onreadystatechange
            // function won't stomp over earlier requests
            
            function ajaxBindCallback() {
                // progressive transitions are from 0 .. 4
                if (ajaxRequest.readyState == 4) {
                    // 200 is the successful response code
                    if (ajaxRequest.status == 200) {
                        ajaxCallback (ajaxRequest.responseText, elementname);
                    } else {
                        // error handling here
                        
                    }
                }
            }
            
            var ajaxRequest = null;
            var ajaxCallback = callback;
            
            // bind the call back, then do the request
            if (window.XMLHttpRequest) {
                // mozilla, firefox, etc will get here
                ajaxRequest = new XMLHttpRequest();
                ajaxRequest.onreadystatechange = ajaxBindCallback;
                ajaxRequest.open("GET", url, true);
                ajaxRequest.send(null);
            } else if (window.ActiveXObject) {
                // IE, of course, has its own way
                ajaxRequest = new ActiveXObject("Microsoft.XMLHTTP");
                
                if (ajaxRequest) {
                    ajaxRequest.onreadystatechange = ajaxBindCallback;
                    ajaxRequest.open("GET", url, true);
                    ajaxRequest.send();
                }
            }
            
        }
        
        function updateSubflowGroup (sometext, subflowname) {
            get_object(subflowname + '_data').innerHTML = sometext;
        }
        
        function updateSubflow (sometext, subflowname) {
            get_object(subflowname + '_data').innerHTML = sometext;
            
            // get the state of this newly updated subflow
            subflowstate = get_object(subflowname + '_state').innerHTML;
            
            // update the subflow image
            subflowstateimg = get_object(subflowname + '_img');
            subflowstateimg.src   = '/ergatis/status_' + subflowstate + '.png';
            subflowstateimg.title = subflowstate;
            subflowstateimg.alt   = subflowstate;
            
            // set the background color to white
            get_object(subflowname + '_data').style.backgroundColor = 'rgb(255,255,255)';
        }

        function toggle_group_info(subflowname) {
            subflownameinfo = get_object(subflowname + "_info");
        
            // if the group info block is hidden, display it
            if ( subflownameinfo.style.display == 'none' ) {
                subflownameinfo.style.display = 'block';
                get_object(subflowname + '_infolabel').innerHTML = 'hide group info';
            
            // else, hide it
            } else {
                subflownameinfo.style.display = 'none';
                get_object(subflowname + '_infolabel').innerHTML = 'show group info';
            }
        }
        
        function toggle_cmd_info(cmd_id) {
            cmdinfoblock = get_object(cmd_id);
            
            // if the command info block is hidden, display it
            if ( cmdinfoblock.style.display == 'none' ) {
                cmdinfoblock.style.display = 'block';
                get_object(cmd_id + '_infolabel').innerHTML = 'hide info';
            
            // else, hide it
            } else {
                cmdinfoblock.style.display = 'none';
                get_object(cmd_id + '_infolabel').innerHTML = 'show info';
            }
        }
        
        function get_object(name) {
            var ns4 = (document.layers) ? true : false;
            var w3c = (document.getElementById) ? true : false;
            var ie4 = (document.all) ? true : false;
            
            if (ns4) return eval('document.' + name);
            if (w3c) return document.getElementById(name);
            if (ie4) return eval('document.all.' + name);
            return false;
        }
        
        // this provides a random string to append
        //  to the end of URL fetches to prevent browser
        //  caching. 
        function no_cache_string () {
            return (   Math.round(  ( Math.random() * 999999 ) + 1  )   );
        }
        
        function reload_subflow(subflowname, subflowfile) {
            // set the background partially grey
            get_object(subflowname + '_data').style.backgroundColor = 'rgb(225,225,225)';
            
            sendElementUpdateRequest('./subflow_summary.cgi?xml_input=' + subflowfile + '&nocache=' + no_cache_string(), updateSubflow, subflowname);
        }
    </script>
</head>

<body>
<div id='bannertop'>
     <a href='/cgi-bin/ergatis/index.cgi'><img src='/ergatis/banner_main.png' border=0/></a>
</div>
<div id='bannerbottom'>
    <div id='pipelinesummary'>
        <div id='pipeline'><a href='/cgi-bin/ergatis/view_formatted_xml_source.cgi?file=$pipeline_xml' target='_blank'>$pipeline_xml</a></div>
        <div class='pipelinestat' id='pipelinestart'><strong>start:</strong> $starttime</div>
        <div class='pipelinestat' id='pipelineend'><strong>end:</strong> $endtime</div>
        <div class='pipelinestat' id='pipelinelastmod'><strong>last mod:</strong> $lastmodtime</div><br>
        <div class='pipelinestat' id='pipelinestate'><strong>state:</strong> $state</div>
        <div class='pipelinestat' id='pipelineuser'><strong>user:</strong> $user</div>
        <div class='pipelinestat' id='pipelineruntime'><strong>runtime:</strong> $runtime</div><br>
        <div class='pipelinestat'><strong>project:</strong> <span id='projectid'></span></div>
        <div class='pipelinestat' id='projectquota'><strong>quota:</strong> $quotastring</div>
        <div class='timer' id='pipeline_timer_label'></div>
        <div id='pipelinecommands'>
            <a href='./view_workflow_pipeline.cgi?&instance=$parent_pipeline'><img class='navbutton' src='/ergatis/button_blue_pipeline_view.png' alt='pipeline view' title='pipeline view'></a>
            <a href='http://htc.tigr.org/antware/cgi-bin/sgestatus.cgi'><img class='navbutton' src='/ergatis/button_blue_grid_info.png' alt='grid info' title='grid info'></a>
            <a href='/cgi-bin/ergatis/view_formatted_xml_source.cgi?file=$pipeline_xml' target='_blank'><img class='navbutton' src='/ergatis/button_blue_xml.png' alt='View XML' title='View XML'></a>
        </div>

    </div>
</div>
<div id='workflowcontainer'>
    <div class='start'>
        $component_name start
    </div>
HeAdER
}

sub print_footer {
    print <<FooTER;
    <div class='end'>
        end
    </div>
</div> <!-- close the workflow container -->
</body>

</html>
FooTER
}
