#!/usr/bin/perl -w

use strict;
use CGI qw(:standard);
use CGI::Carp qw(fatalsToBrowser);
use CGI::Cookie;

my $q = new CGI;

## set a new cookie of the same name, but make it already expired
my $c = new CGI::Cookie(-name    =>  'ergatis_user',
                        -value   =>  '',
                        -expires =>  '-1M',);

print "Set-Cookie: $c\n";
print redirect(-uri=>"./index.cgi");
