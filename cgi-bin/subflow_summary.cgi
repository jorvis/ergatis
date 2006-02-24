#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use File::Basename;
use Monitor;
use XML::Twig;

my $q = new CGI;

print $q->header( -type => 'text/html' );

my $xml_input = $q->param("xml_input") || die "pass xml_input";
my $subflow_name  = basename($xml_input, ('.xml', '.gz'));
my $subflow_state = 'unknown';
my $subflow_start = '';
my $subflow_end   = '';

print "    <h1>analysis steps</h1>\n";

my $xml_input_fh;
if ($xml_input =~ /\.gz/) {
    open($xml_input_fh, "<:gzip", "$xml_input") || die "can't read $xml_input: $!"; 
} elsif ( ! -e $xml_input ) {
    if ( -e "$xml_input.gz" ) {
        open($xml_input_fh, "<:gzip", "$xml_input.gz") || die "can't read $xml_input: $!";
    } else {
        print "<div class='command'>steps not yet available</div>\n";
        exit;    
    }

} else {
    open($xml_input_fh, "<$xml_input") || die "can't read $xml_input: $!";       
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
$twig->parse($xml_input_fh);

