#!/usr/bin/perl -w

use strict;
use CGI qw(:standard);
use CGI::Carp qw(fatalsToBrowser);
use CGI::Cookie;
use Ergatis::Common;

my $q = new CGI;
my $redirect_url = $q->param('redirect_url') || "./index.cgi";
$redirect_url = parse_referer_url($redirect_url);
my %cookies = fetch CGI::Cookie;

my $session_id = $cookies{'ergatis_user'}->value if ($cookies{'ergatis_user'});

## Kill off our ergatis_users cookie to log us off
my $c = new CGI::Cookie(-name    =>  'ergatis_user',
                        -value   =>  '',
                        -expires =>  '-1M',);

## Create a new cookie storing our CGI session ID so we can make use of it 
## on subsequent logins
if (! $cookies{'CGISESSID'} && $session_id) {
    print STDERR "Setting our CGISESSID cookie!\n";
    my $session_cookie = new CGI::Cookie(-name    => 'CGISESSID',
                                         -value   => $session_id,
                                         -expires => '+3M',);
    print "Set-Cookie: $session_cookie\n";
}

print "Set-Cookie: $c\n";

print redirect(-uri=> $redirect_url);
