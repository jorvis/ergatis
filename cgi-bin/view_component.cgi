#!/local/perl/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use XML::Twig;


my $q = new CGI;

print $q->header( -type => 'text/html' );

my $pipeline_xml = $q->param("pipeline_xml") || die "pass pipeline_xml";

print_header();

my $twig = XML::Twig->new( );
$twig->parsefile($pipeline_xml);

my $parent_commandset = $twig->root->first_child('commandSet');

## look at each of the children of the root
foreach my $child ( $parent_commandset->children() ) {
    
    ## if this is a command, print its name
    if ($child->gi eq 'command') {
        ## get its state
        my $command_state = 'unknown';
        if ( $child->first_child('state') ) {
            $command_state = $child->first_child('state')->text();
        }
        
        print "    <li class='command'><span class='$command_state'>$command_state</span> " . 
              $child->first_child('name')->text() . "</li>\n";
    
    ## if it is a commandSet, it should be a file-based subflow
    } elsif ($child->gi eq 'commandSet') {
    
        ## make sure it has a fileName element
        if ( $child->has_child('fileName') ) {
            print "    <li>\n";
            print "        <ul class='subflow'>\n";
            print "            <li><h1>" . $child->first_child('name')->text() . "</h1></li>\n";
            
            &parse_groups_xml( $child->first_child('fileName')->text() );
            
            print "        </ul>\n";
            print "    </li>\n";
        }
    }
}

print_footer();

exit(0);



sub parse_groups_xml {
    my $filename = shift;
    
    ## make sure this is a groups.xml file
    if ( $filename !~ /groups.xml/ ) {
        print "            <li>unable to handle $filename</li>\n";
        return;
    }
    
    ## make sure it exists
    if (! -e $filename) {
        print "            <li>not yet created</li>\n";
    }
    
    my ($component_start_time, $component_end_time, $component_state);
    
    ## create the twig
    my $twig = XML::Twig->new( twig_roots => {
                                    'commandSet/startTime' => 
                                        sub {
                                              my ($t, $elt) = @_;
                                              $component_start_time = $elt->text();
                                              print "            <li>start: $component_start_time</li>\n";
                                        },
                                    'commandSet/endTime'   => 
                                        sub {
                                              my ($t, $elt) = @_;
                                              $component_end_time = $elt->text();
                                              print "            <li>end: $component_end_time</li>\n";
                                        },
                                    'commandSet/state'     => 
                                        sub {
                                              my ($t, $elt) = @_;
                                              $component_state = $elt->text();
                                              print "            <li>state: $component_state</li>\n";
                                        },
                                    'command'              => \&process_command,
                               }
                             );
    $twig->parsefile($filename);
    
}

sub process_command {
    my ($twig, $command) = @_;
    my $state = 'unknown';
    my $jobId = '';
    my $execution_host = '';
    my $name = '';
    my $subflow_file = '';
    
    ## need to parse through the params to get the one that references the
    ##  instance file descriptor
    for my $param ( $command->children('param') ) {
        ## it will have a key element --instance
        if ( $param->first_child('key')->text() eq '--instance' ) {
            $subflow_file = $param->first_child('value')->text();
        }
    }
    
    ## pull the name out of the subflow_file:
    if ( $subflow_file =~ /subflow(\d+)groups(\d+).xml/ ) {
        $name = "subflow$1groups$2";
    }
    
    ## get the state, if it has one
    if ( $command->first_child('state') ) {
        $state = $command->first_child('state')->text();
    }
    
    print <<SubflowGroupsBar;
            <li id='${name}_bar' class='subflowbar'>
                <img id='${name}_arrow' src='/cram/arrow_right.gif' onclick='toggle_subflow_display("$name");'>
                group state: <span class='$state'>$state</span>
            </li>
            <li id='${name}_data' class='subflowdata' style='display: none;'>
                subflow data here
            </li>
SubflowGroupsBar
}



sub print_header {
    print <<HeAdER;
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
"http://www.w3.org/TR/html4/strict.dtd">

<html>

<head>
    <meta http-equiv="Content-Language" content="en-us">
    <meta http-equiv="Content-Type" content="text/html; charset=windows-1252">
    <title>component view</title>
    <!-- <link rel="stylesheet" type="text/css" href="/concept/builder.css"> -->
    <style type="text/css">
        body {
            font-family: verdana, helvetica, arial, sans-serif;
            font-size: 10px;
            margin: 0;
            padding: 0;
        }
        
        #workflowcontainer {
            margin-top: 10px;
        }

        ul {
            margin-left: 20px;
            padding-left: 0px;
        }
        
        ul.subflow {
            margin-bottom: 10px;
        }

        li h1 {
            font-weight: bold; 
            font-size: 100%; 
        }
        
        li {
            margin-left: 0px;
            padding-left: 0px;
            list-style-type: none;
        }
        
        li.subflowbar {
            border-bottom: 1px solid grey;
            width: 500px;
            margin-top: 3px;
        }
        
        li.subflowdata {
            padding: 5px 0px 0px 20px;
            display: none;
        }
        
        ul.start, ul.end {
            border: 1px solid grey;
            width: 300px;
            background-color: rgb(220, 220, 220);
            padding-left: 0px;
            margin: 10px 0px 10px 0px;
        }
        
        ul.start li, ul.end li {
            padding-left: 5px;
        }
        

        span.complete {
            color: rgb(0,200,0);
        }
        span.incomplete {
            color: rgb(75,75,75);
        }
        span.failed {
            color: rgb(200,0,0);
        }
        span.pending {
            color: rgb(200,200,0);
        }
        span.errors {
            color: rgb(200,0,0);
        }
        span.error {
            color: rgb(200,0,0);
        }
        span.running {
            color: rgb(0,0,200);
        }
        span.waiting {
            color: rgb(200,200,0);
        }
        span.interrupted {
            color: rgb(200,0,200);
        }
        span.total {
            color: rgb(0,0,0);
        }
        span.unknown {
            color: rgb(0,0,0);
        }
    </style>
    <script type="text/javascript">
        function toggle_subflow_display(subflowname) {
            subflownamedata = get_object(subflowname + "_data");
            
            // is it hidden?
            if ( subflownamedata.style.display == 'none' ) {
                subflownamedata.style.display = 'block';
                get_object(subflowname + '_arrow').src = '/cram/arrow_down.gif';
            } else {
                subflownamedata.style.display = 'none';
                get_object(subflowname + '_arrow').src = '/cram/arrow_right.gif';
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
    </script>
</head>

<body>

<ul id='workflowcontainer'>
    <li>
        <ul class='start'>
            <li><h1>start</h1></li>
        </ul>
    </li>
HeAdER
}

sub print_footer {
    print <<FooTER;
    <li>
        <ul class='end'>
            <li><h1>end</h1></li>
        </ul>
    </li>
</ul> <!-- close the workflow container -->
</body>

</html>
FooTER
}
