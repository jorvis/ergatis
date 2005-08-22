#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use XML::Twig;

my $q = new CGI;

print $q->header( -type => 'text/html' );

my $xml_input = $q->param("xml_input") || die "pass xml_input";
my $subflow_state = 'unknown';
my $subflow_start = '';
my $subflow_end   = '';

print "    <h1>analysis steps</h1>\n";

my $twig = XML::Twig->new( twig_roots => {
                                'commandSet/state'     => sub {
                                                            my ($t, $elt) = @_;
                                                            $subflow_state = $elt->text;
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



sub process_command {
    my ($twig, $command) = @_;
    
    my $name  = $command->first_child('name')->text;
    my $state = $command->first_child('state')->text;
    
    print <<CommAnD;
    <div class='command'>
        <div class='leftside'>
            <img class='status' src='/cram/status_$state.png' title='$state' alt='$state'>
            $name
        </div>
        <div class='rightside'>
            <span class='minor'>n/a</span>
        </div>
    </div>
CommAnD

    my $ret_value = 'unknown';
    my $message = '';
    
    ## if there is a status and a message, grab it
    if ( $command->first_child('status') ) {
        if ( $command->first_child('status')->first_child('retValue') ) {
            $ret_value = $command->first_child('status')->first_child('retValue')->text;
        }
        
        if ( $command->first_child('status')->first_child('message') ) {
            $message = $command->first_child('status')->first_child('message')->text;
        }
        
        if ( $message ) {
            print <<messageBLOCK;
        <div class='messageblock'>
            return value: $ret_value<br>
            message: $message
        </div>
messageBLOCK
        }
    }

    $twig->purge;
}











