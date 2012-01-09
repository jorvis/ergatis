#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Monitor;
use File::Basename;
use HTML::Template;
use XML::Twig;

my $q = new CGI;

print $q->header( -type => 'text/html' );

my $tmpl = HTML::Template->new( filename => 'templates/subflow_summary.tmpl',
                                die_on_bad_params => 1,
                              );

my $xml_input = $q->param("xml_input") || die "pass xml_input";

my $pipeline_id = '';
if ( $xml_input =~ m|(.+/(.+?))/workflow/runtime/(.+?)/([A-Z0-9]+)_(.+?)/| ) {
    $pipeline_id = $4;
}

my $subflow_name  = basename($xml_input, ('.xml', '.gz'));
my $subflow_state = 'unknown';
my $subflow_start = '';
my $subflow_end   = '';
my $subflow_found = 0;

my $xml_input_fh;
if ($xml_input =~ /\.gz/) {
    open($xml_input_fh, "<:gzip", "$xml_input") || die "can't read $xml_input: $!"; 
} elsif ( ! -e $xml_input ) {
    if ( -e "$xml_input.gz" ) {
        open($xml_input_fh, "<:gzip", "$xml_input.gz") || die "can't read $xml_input: $!";
    }
} else {
    open($xml_input_fh, "<$xml_input") || die "can't read $xml_input: $!";       
}

$subflow_found = 1 if ($xml_input_fh);

my $elements = [];

if ($subflow_found) {
    my $twig = XML::Twig->new( twig_roots => {
                                    'commandSet/state'     => sub {
                                                                my ($t, $elt) = @_;
                                                                $subflow_state = $elt->text;
                                                              },
                                    'command'              => sub {
                                                                my ($t, $elt) = @_;
                                                                my %parts = &process_command($t, $elt);
                                                                $parts{'pipeline_id'} = $pipeline_id;
                                                                push @$elements, \%parts;
                                                              }
                               },
                          );
    $twig->parse($xml_input_fh);
}
#add links to command
for(my $i=0;$i<@$elements;$i++){
    my $cs_string=$elements->[$i]->{'command_string'};
    my $cs_formatted_string = $cs_string;
    my @cs_string_elts = split(/\s+/,$cs_string);
    print STDERR $cs_string,"\n";
    foreach my $e (@cs_string_elts){
	if( -f "$e"){
	    if($e =~ /\.xml/){
		$cs_formatted_string =~ s/$e/<a href='view_formatted_xml_source.cgi?file=$e&pipeline_id=$pipeline_id'>$e\<\/a\>/;
	    }
	    else{
		$cs_formatted_string =~ s/$e/<a href='view_formatted_log_source.cgi?file=$e&pipeline_id=$pipeline_id'>$e\<\/a\>/;
	    }
	}
	elsif($e =~ /\=/){
	    my($key,$value) = split(/=/,$e);
	    if( -f "$value"){
		if($value =~ /\.xml/){
		    $cs_formatted_string =~ s/$value/<a href='view_formatted_xml_source.cgi?file=$value&pipeline_id=$pipeline_id'>$value\<\/a\>/;
		}
		else{
		    $cs_formatted_string =~ s/$value/<a href='view_formatted_log_source.cgi?file=$value&pipeline_id=$pipeline_id'>$value\<\/a\>/;
		}
	    }
	}
    }
    $elements->[$i]->{'command_string'} = $cs_formatted_string;
}

$tmpl->param( SUBFLOW_NAME  => $subflow_name );
$tmpl->param( SUBFLOW_STATE => $subflow_state );
$tmpl->param( SUBFLOW_FOUND => $subflow_found );
$tmpl->param( ELEMENTS      => $elements );

print $tmpl->output;
