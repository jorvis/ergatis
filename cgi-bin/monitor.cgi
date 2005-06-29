#!/local/perl/bin/perl -w

## this shouldn't be used at all yet.  i'm just playing
##  with it right now.

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Date::Manip;
use HTML::Template;

my $q = new CGI;

print $q->header( -type => 'text/html' );

my $xml_input = $q->param("xml_input") || die "pass xml_input";

print_header();

print "<br><br><table>\n";

if (-d $xml_input) {
    my @files = glob( "$xml_input/*.xml.instance" );
    
    for my $file (@files) {
        print "   <tr>\n";
 
        my ($start, $state) = get_start_and_state(\$file);
        my $uid = -o $file;

        print "      <td><a href=\"./xml2tree.cgi?xml_input=$file\">$file</a></td>\n";
        print "      <td class=\"state $state\">$state</td>\n";
        print "      <td>" . getpwuid( (stat $file)[4] ) . "</td>\n";
        print "      <td>" . UnixDate($start, "%b %e, %T") . "</td>\n";
       
        print "   </tr>\n";
    }

} elsif (-f $xml_input) {
    die "can't handle plain files yet";
} else {
    die "don't know what to do with $xml_input";
}


print_footer();

exit(0);

sub get_start_and_state {
    ## we only need a small amount of information here, and this is faster
    ##  than building a twig
    my $file = shift;
    my ($startTime, $state);
    
    open(my $ifh, "<$$file") || die "can't read file $$file : $!";
    
    while (<$ifh>) {
        ## is this a startTime line?
        if (m|<startTime>(.+)</startTime>|) {
            $startTime = $1;
        ## is this a state line?
        } elsif (m|<state>(.+)</state>|) {
            $state = $1;
        }
        
        ## quit once we've found both
        return ($startTime, $state) if ($startTime && $state);
    }
    
}


sub print_header {
    print <<HeAdER;

<html>

<head>
    <meta http-equiv="Content-Language" content="en-us">
    <meta http-equiv="Content-Type" content="text/html; charset=windows-1252">
    <title>pre-pre-alpha workflow interface prototype</title>
    <!-- <link rel="stylesheet" type="text/css" href="/concept/builder.css"> -->
    <style type="text/css">
        body {
            font-family: verdana, helvetica, arial, sans-serif;
            font-size: 10px;
            margin: 0;
            padding: 0;
        }
        td {
            font-family: verdana, helvetica, arial, sans-serif;
            font-size: 10px;
        }
        td.state {
            text-align: center;
            border: 1px solid black;
        }
        td.complete {
            background-color: rgb(0,255,0);
        }
        td.failed {
            background-color: rgb(255,0,0);
        }
        td.interrupted {
            background-color: rgb(255,0,255);
        }
        td.running {
            background-color: rgb(0,0,255);
        }
        a {
            color: rgb(0,0,255);
            text-decoration: none;
        }
        a:hover {
            text-decoration: underline;
        }
    </style>
</head>

<body>

HeAdER
}

sub print_footer {
    print <<FooTER;
</body>

</html>
FooTER
}
