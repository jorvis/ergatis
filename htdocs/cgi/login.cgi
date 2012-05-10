#!/usr/bin/perl -w

use strict;

use CGI qw(:standard);
use CGI::Carp qw(fatalsToBrowser);
use CGI::Cookie;
use CGI::Session;
use Ergatis::ConfigFile;
use Ergatis::Common;
use Storable;

my $q = new CGI;
my %cookies = fetch CGI::Cookie;

## read the config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );
my $user_attempted = $q->param('login_user');
my $pass_attempted = $q->param('login_pass');
my $referer = $q->param('redirect_url');

## for now, this will be a user name
my $valid_user = '';

my $auth_module = '';

if ( $ergatis_cfg->val('authentication', 'authentication_method') eq 'kerberos' ) {
    
    $auth_module = "Authen::Simple::Kerberos";
    eval "use $auth_module";
        die "Couldn't load module : $!n" if ($@);

    ## validate the user here
    my $realm = $ergatis_cfg->val('authentication', 'kerberos_realm');

    my $kerberos = Authen::Simple::Kerberos->new(
        realm => $realm,
    );

    if ( $kerberos->authenticate( $user_attempted, $pass_attempted) ) {
        $valid_user = $user_attempted;
    }

} elsif ( $ergatis_cfg->val('authentication', 'authentication_method') eq 'ad' ) {
    
    $auth_module = "Authen::Simple::ActiveDirectory";
    eval "use $auth_module";
        die "Couldn't load module : $!n" if ($@);

    ## validate the user here
    my $host      = $ergatis_cfg->val('authentication', 'ldap_host');
    my $port      = $ergatis_cfg->val('authentication', 'ldap_port');
    my $principal = $ergatis_cfg->val('authentication', 'ad_principal');

    my $ad = Authen::Simple::ActiveDirectory->new(
        host      => $host,
        port      => $port,
        principal => $principal,
    );

    if ( $ad->authenticate( $user_attempted, $pass_attempted) ) {
        $valid_user = $user_attempted;
    }

}  elsif ( $ergatis_cfg->val('authentication', 'authentication_method') eq 'ldap' ) {

    $auth_module = "Authen::Simple::LDAP";
    eval "use $auth_module";
        die "Couldn't load module : $!n" if ($@);

    ## validate the user here
    my $host   = $ergatis_cfg->val('authentication', 'ldap_host');
    my $port   = $ergatis_cfg->val('authentication', 'ldap_port');
    my $basedn = $ergatis_cfg->val('authentication', 'ldap_basedn');

    my $ldap = Authen::Simple::LDAP->new(
        host   => $host,
        port   => $port,
        basedn => $basedn,
    );

    if ( $ldap->authenticate( $user_attempted, $pass_attempted) ) {
        $valid_user = $user_attempted;
    }

} elsif ( $ergatis_cfg->val('authentication', 'authentication_method') eq 'ldap-tls' ) {

    $auth_module = "Net::LDAP";
    eval "use $auth_module";
        die "Couldn't load module : $!n" if ($@);

    ## validate the user here
    my $host   = $ergatis_cfg->val('authentication', 'ldap_host');
    my $basedn = $ergatis_cfg->val('authentication', 'ldap_basedn');

    my $ldap;
    
    eval {
        $ldap = Net::LDAP->new($host);
    };

    if ($@) {
        croak("Unable to establish Net::LDAP connection to $host.");
    }
    
    my $msg = $ldap->start_tls(
        verify => 'require',
        cafile => $ergatis_cfg->val('authentication', 'ldap_tls_certificate'),
    );

    if ( $msg->code() ) {
        croak("Unable to start TLS connection to LDAP server.");
    }

    my $dn_string = "uid=$user_attempted," . $ergatis_cfg->val('authentication', 'ldap_basedn');

    my $login_msg = $ldap->bind("$dn_string", password => $pass_attempted );
    
    if ( ! $login_msg->code() ) {
        $valid_user = $user_attempted;
    }
} 


if ( $valid_user ) {
 
    unless ( $ergatis_cfg->val('authentication', 'session_db_dir') ) {
        croak("A session database directory is required for authentication to function in Ergatis")
    }
    
    ## We will be using sessions to handle recording state that is tied to the user logged in.
    ## We want to check if a session ID already exists for the user attemping to login and if so
    ## retrive the already created session for use.
    my $sid = undef;
    if ($cookies{'CGISESSID'}) {
        $sid = $cookies{'CGISESSID'}->value;
    }

    my $session = get_session($ergatis_cfg, $sid);
    if ( ! $session->param('username') ) {
        $session->param('username', $valid_user);  
        $session->flush();
    }

    my $pipelines_file = $ergatis_cfg->val('authentication', 'session_db_dir') . "/$valid_user.pipelines";
    if ( $ergatis_cfg->val('authentication', 'per_account_pipeline_security') ) {
        my $last_mod = 0;
        if (-e $pipelines_file) {
            $session->param( 'pipelines', retrieve($pipelines_file) );
            $last_mod = (stat ($pipelines_file))[9];
        }

        $session->param('pipelines_last_mod', $last_mod);
        $session->flush();
    }

    ## Use a plain cookie to handle any user state information
    my $c = new CGI::Cookie(-name    =>  'ergatis_user',
                            -value   =>  $session->id,
                            -expires =>  '+3M',);

    print "Set-Cookie: $c\n";
    
    ## Redirect to the page we just came from. Cookies might be a terrible idea here 
    ## as they could be spoofed.
    my $redirect_url = $referer || "./index.cgi";
    $redirect_url = "./index.cgi" if $redirect_url =~ "login_form";
   
    print redirect(-uri => $redirect_url);

## if don't pass
} else {
    ## redirect back to login form
    print redirect(-uri => "./login_form.cgi?failed=1");
}

