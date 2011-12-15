#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use CGI::Cookie;
use Ergatis::Common;
use Ergatis::ConfigFile;
use HTML::Template;

my $q = new CGI;

print $q->header( -type => 'text/html');
my $redirect_url = $q->param('redirect_url') || $ENV{HTTP_REFERER};
$redirect_url = parse_referer_url($redirect_url);

my $tmpl = HTML::Template->new( filename => 'templates/login_form.tmpl',
                                die_on_bad_params => 1,
                              );

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

## check for any reasons why we shouldn't display the form
my $form_ready = 0;
my $form_not_ready_msg = '';

if ( $ergatis_cfg->val('authentication', 'authentication_method') eq 'open' ) {
    $form_not_ready_msg = "Sorry, authentication is set to 'open' in the ergatis.ini file, logins aren't supported in this mode."

} elsif ( $ergatis_cfg->val('authentication', 'authentication_method') eq 'ldap' ) {
    $form_ready = 1;

} elsif ( $ergatis_cfg->val('authentication', 'authentication_method') eq 'ldap-tls' ) {
    $form_ready = 1;

} elsif ( $ergatis_cfg->val('authentication', 'authentication_method') eq 'kerberos' ) {
    $form_ready = 1;
    
} elsif ( $ergatis_cfg->val('authentication', 'authentication_method') eq 'ad' ) {
    $form_ready = 1;
    
} else {
    $form_not_ready_msg = "Error: unrecognized authentication_method value in ergatis.ini.  Please check the documentation in that file.";
}

## check the method

$tmpl->param( FAILED              => $q->param('failed') || 0 );
$tmpl->param( FORM_READY          => $form_ready );
$tmpl->param( FORM_NOT_READY_MSG  => $form_not_ready_msg );
$tmpl->param( QUICK_LINKS         => &get_quick_links($ergatis_cfg) );
$tmpl->param( REDIRECT_URL         => $redirect_url );
$tmpl->param( SUBMENU_LINKS       => [
                                        #{ label => 'create project', is_last => 1, url => './create_project_form.cgi' },
                                     ] );

print $tmpl->output;
