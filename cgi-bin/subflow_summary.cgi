#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Date::Manip;
use File::Basename;
use Monitor;
use POSIX;
use XML::Twig;

my $q = new CGI;

print $q->header( -type => 'text/html' );

my $xml_input = $q->param("xml_input") || die "pass xml_input";
my $subflow_name  = basename($xml_input, '.xml');
my $subflow_state = 'unknown';
my $subflow_start = '';
my $subflow_end   = '';

print "    <h1>analysis steps</h1>\n";

## make sure the file exists, else print message
if (! -e $xml_input) {
    print "<div class='command'>steps not yet available</div>\n";
    exit;
}

my $twig = XML::Twig->new( twig_roots => {
                                'commandSet/state'     => sub {
                                                            my ($t, $elt) = @_;
                                                            $subflow_state = $elt->text;
                                                            print "<span class='hidden' id='${subflow_name}_state'>$subflow_state</span>\n";
                                                          },
                                'commandSet/startTime' => sub {
                                                            my ($t, $elt) = @_;
                                                            $subflow_start = $elt->text;
                                                          },
                                'commandSet/endTime'   => sub {
                                                            my ($t, $elt) = @_;
                                                            $subflow_end = $elt->text;
                                                          },
                                'command'              => \&process_command,
                           },
                      );
$twig->parsefile($xml_input);

