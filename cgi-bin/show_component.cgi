#!/usr/local/bin/perl

use strict;
use CGI qw(:standard);
use CGI::Carp qw(fatalsToBrowser set_message);

use Sybil;

my $logger;

BEGIN {
    $logger = new Coati::Logger('LOG_FILE'=>"/tmp/show_component.log",
				'LOG_LEVEL'=>&param('DEBUG'));  
    set_message(\&Manatee::GetManateeTemplate::get_error_handler);
    &Manatee::GetManateeParameters::import_stored_parameters();
}

my $xmltemplate=param('xmltemplate');

my $template = new Manatee::GetManateeTemplate('FILE'=>"show_component.tt");

