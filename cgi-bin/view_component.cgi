#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Date::Manip;
use POSIX;
use XML::Twig;


my $q = new CGI;

print $q->header( -type => 'text/html' );

my $pipeline_xml = $q->param("pipeline_xml") || die "pass pipeline_xml";

my $twig = XML::Twig->new( );
$twig->parsefile($pipeline_xml);

my $parent_commandset = $twig->root->first_child('commandSet');
my $component_state = $parent_commandset->first_child('state')->text || 'unknown';

my $component_name = 'component';
if ( $parent_commandset->first_child('name')->text =~ /(.+) workflow/) {
    $component_name = $1;
}

my $parent_pipeline = $parent_commandset->first_child('parentFileName')->text || '';

my $component_project = '?';
if ( $parent_pipeline =~ /.+\/(.+?)\/Workflow/ ) {
    $component_project = $1;
}

print_header();

## look at each of the children of the root
foreach my $child ( $parent_commandset->children() ) {
    
    ## if this is a command, print its name
    if ($child->gi eq 'command') {
        ## get its state
        my $command_state = 'unknown';
        if ( $child->first_child('state') ) {
            $command_state = $child->first_child('state')->text();
        }
        
        print "    <div class='command'>" .
              "<img class='status' src='/cram/status_${command_state}.png' title='$command_state' alt='$command_state'>" .
              $child->first_child('name')->text() . "</div>\n";
    
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
    
    ## make sure this is a groups.xml file
    if ( $filename !~ /groups.xml/ ) {
        print "            <div>unable to handle $filename</div>\n";
        return;
    }
    
    ## make sure it exists
    if (! -e $filename) {
        print "            <div>not yet created</div>\n";
    }
    
    my ($component_start_time, $component_end_time, $component_state);
    
    ## create the twig
    my $twig = XML::Twig->new( twig_roots => {
#                                    'commandSet/startTime' => 
#                                        sub {
#                                              my ($t, $elt) = @_;
#                                              $component_start_time = $elt->text();
#                                              #print "            <div>start: $component_start_time</div>\n";
#                                        },
#                                    'commandSet/endTime'   => 
#                                        sub {
#                                              my ($t, $elt) = @_;
#                                              $component_end_time = $elt->text();
#                                              #print "            <div>end: $component_end_time</div>\n";
#                                        },
#                                    'commandSet/state'     => 
#                                        sub {
#                                              my ($t, $elt) = @_;
#                                              $component_state = $elt->text();
#                                              #print "            <div>state: $component_state</div>\n";
#                                        },
                                    'command'              => \&process_command,
                               }
                             );
    $twig->parsefile($filename);
    
}

sub process_command {
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
        $execution_host = $command->first_child('dceSpec')->first_child('executionHost')->text();
        $grid_id        = $command->first_child('dceSpec')->first_child('jobID')->text();
    }
    
    if ( $command->first_child('id') ) {
        $workflow_id = $command->first_child('id')->text();
    }
    
    my $start_time_obj = ParseDate( $command->first_child('startTime')->text );
    my $start_time     = UnixDate( $start_time_obj, "%c" );
    my $end_time_obj   = ParseDate( $command->first_child('endTime')->text );
    my $end_time       = UnixDate( $end_time_obj, "%c" );

    ## we can calculate runtime only if start and end time are known, or if start is known and state is running
    my $runtime = '?';
    if ($start_time_obj) {
        if ($end_time_obj) {
            $runtime = strftime( "%H hr %M min %S sec", reverse split(/:/, DateCalc($start_time_obj, $end_time_obj)) );
        } elsif ($state eq 'running') {
            $runtime = strftime( "%H hr %M min %S sec", reverse split(/:/, DateCalc("now", $start_time_obj)) ) . ' ...';
        }
    }
    
    print <<SubflowGroupBar;
        <div id='${name}_bar' class='subflowbar'>
            <div class='leftside'>
                <img id='${name}_arrow' class='expander' src='/cram/arrow_right.gif' onclick='toggle_subflowgroup_display("$name", "$file");' alt='expand' title='expand'>
                <img class='status' src='/cram/status_$state.png' title='$state' alt='$state'>
                <span class='minor'>group $group_num</span>
            </div>
            <div class='rightside'>
                <span class='infolabel' id='${name}_infolabel' onclick='toggle_group_info("$name")'>show group info</span>
            </div>
        </div>
        <div id='${name}_info' class='subflowinfo' style='display: none;'>
            <table>
                <tr><th>workflow id:</th><td>$workflow_id</td></tr>
                <tr><th>state:</th><td>$state</td></tr>
                <tr><th>start time:</th><td>$start_time</td></tr>
                <tr><th>end time:</th><td>$end_time</td></tr>
                <tr><th>duration:</th><td>$runtime</td></tr>
                <tr><th>grid id:</th><td>$grid_id</td></tr>
                <tr><th>execution host:</th><td>$execution_host</td></tr>
                <tr><th colspan='2'>xml:</th></tr>
                <tr><td colspan='2'><a href='./view_formatted_xml_source.cgi?file=$file'>$file</a></td></tr>
            </table>
        </div>
        <div id='${name}_data' class='subflowdata' style='display: none;'></div>
SubflowGroupBar
}



sub print_header {
    print <<HeAdER;
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
"http://www.w3.org/TR/html4/strict.dtd">

<html>

<head>
    <meta http-equiv="Content-Language" content="en-us">
    <meta http-equiv="Content-Type" content="text/html; charset=windows-1252">
    <title>$component_project|$component_name|$component_state</title>
    <!-- <link rel="stylesheet" type="text/css" href="/concept/builder.css"> -->
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
        img.expander {
            cursor: pointer;
        }
        #workflowcontainer {
            margin: 10px;
        }
        div.subflow {
            margin-bottom: 10px;
        }
        div h1 {
            font-weight: bold; 
            font-size: 100%; 
        }
        div.subflowbar, div.itembar {
            margin-bottom: 5px;
            width: 500px;
            padding-bottom: 2px;
            margin-top: 3px;
        }
        div.subflowbar {
            border-bottom: 1px solid rgb(150,150,150);
        }
        div.command {
            margin-bottom: 2px;
            width: 500px;
            padding-bottom: 2px;
        }
        div.subflowdata, div.itemdata {
            padding: 5px 0px 0px 30px;
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
        }
        span.infolabel {
            color: rgb(0,0,200);
            font-size: 80%;
            cursor: pointer;
        }
        div.leftside {
            float: left;
        }
        div.rightside {
            text-align: right;
        }
        div.subflowinfo {
            padding-left: 30px;
        }
        div.subflowinfo table {
            background-color: rgb(235,235,235);
            padding: 5px;
        }
        td, th {
            font-size: 90%;
            text-align: left;
            padding: 0px;
        }
        th {
            color: rgb(100,100,100);
        }
    </style>
    <script type="text/javascript">
        // the following two functions could probably be combined later
        
        function toggle_subflowgroup_display(subflowname, subflowfile) {
            subflownamedata = get_object(subflowname + "_data");
            
            // is it hidden?
            if ( subflownamedata.style.display == 'none' ) {
                subflownamedata.style.display = 'block';
                get_object(subflowname + '_arrow').src = '/cram/arrow_down.gif';
                get_object(subflowname + '_data').innerHTML = 'parsing subflow data';
                sendElementUpdateRequest('./subflowgroup_summary.cgi?xml_input=' + subflowfile + '&nocache=' + no_cache_string(), updateSubflowGroup, subflowname);
                
            } else {
                subflownamedata.style.display = 'none';
                get_object(subflowname + '_arrow').src = '/cram/arrow_right.gif';
                get_object(subflowname + '_data').innerHTML = '';
            }
        }
        
        function toggle_subflow_display(subflowname, subflowfile) {
            subflownamedata = get_object(subflowname + "_data");
            
            // is it hidden?
            if ( subflownamedata.style.display == 'none' ) {
                subflownamedata.style.display = 'block';
                get_object(subflowname + '_arrow').src = '/cram/arrow_down.gif';
                get_object(subflowname + '_data').innerHTML = 'parsing subflow data';
                sendElementUpdateRequest('./subflow_summary.cgi?xml_input=' + subflowfile + '&nocache=' + no_cache_string(), updateSubflow, subflowname);
                
            } else {
                subflownamedata.style.display = 'none';
                get_object(subflowname + '_arrow').src = '/cram/arrow_right.gif';
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
        }

        function toggle_group_info(subflowname) {
            subflownameinfo = get_object(subflowname + "_info");
        
            // if the group info block is hidden, display it
            if ( subflownameinfo.style.display == 'none' ) {
                subflownameinfo.style.display = 'block';
                get_object(subflowname + '_infolabel').innerHTML = 'hide group info';
            } else {
                subflownameinfo.style.display = 'none';
                get_object(subflowname + '_infolabel').innerHTML = 'show group info';
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
    </script>
</head>

<body>

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
