#!/usr/bin/perl -w

use strict;
use Authen::Simple::Kerberos;
#use Log::Cabin;
use CGI qw(:standard);
use CGI::Carp qw(fatalsToBrowser);
use CGI::Cookie;
use Ergatis::ConfigFile;

my $q = new CGI;

## read the config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );
my $user_attempted = $q->param('login_user');
my $pass_attempted = $q->param('login_pass');

## for now, this will be a user name
my $valid_user = '';

if ( $ergatis_cfg->val('authentication', 'authentication_method') eq 'kerberos' ) {
    ## validate the user here
    my $realm = $ergatis_cfg->val('authentication', 'kerberos_realm');

#    my $logsimple = new Log::Cabin;
#       $logsimple->level( 8 );
#       $logsimple->set_output(*STDERR);

#    my $logger = $logsimple->get_logger('kerberos');


    my $kerberos = Authen::Simple::Kerberos->new(
        realm => $realm,
#        log => $logger,
    );

    if ( $kerberos->authenticate( $user_attempted, $pass_attempted) ) {
        $valid_user = $user_attempted;
    }
}


if ( $valid_user ) {

    ## set cookie
    my $c = new CGI::Cookie(-name    =>  'ergatis_user',
                            -value   =>  $valid_user,
                            -expires =>  '+1M',);

    print "Set-Cookie: $c\n";
    
    ## redirect to the index page
    print redirect(-uri => "./index.cgi");

## if don't pass
} else {
    ## redirect back to login form
    print redirect(-uri => "./login_form.cgi?failed=1");
}
