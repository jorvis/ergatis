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

my $tmpl = HTML::Template->new( filename => 'templates/subflowgroup_summary.tmpl',
                                die_on_bad_params => 0,
                              );

## file will be like: $proj_dir/workflow/runtime/wu-blastp/5407_test/i1/g1/g1.xml
my $xml_input = $q->param("xml_input") || die "pass xml_input";

## grab the pipeline ID from our XML input path
my $pipeline_id = '';
if ( $xml_input =~ m|(.+/(.+?))/workflow/runtime/.+?/([A-Z0-9]+)(_.+)?/| ) { 
    $pipeline_id = $3
}

my $xml_input_fh;
if ($xml_input =~ /\.gz/) {
    open($xml_input_fh, "<:gzip", "$xml_input") || die "can't read $xml_input: $!"; 
} elsif ( ! -e $xml_input && -e "$xml_input.gz" ) {
    open($xml_input_fh, "<:gzip", "$xml_input.gz") || die "can't read $xml_input: $!"; 
} else {
    open($xml_input_fh, "<$xml_input") || die "can't read $xml_input: $!";       
}

## slow DOM parse
my $sg_twig = XML::Twig->new();
   $sg_twig->parse($xml_input_fh);

my $parent_commandset = $sg_twig->root->first_child('commandSet');

my $elements = [];

## the gN.xml file will have a series of commands for things like directory creation and variable
#   replacement, then a file-based subflow pointing to the gN.iter.xml

for my $child ( $parent_commandset->children() ) {
    ## if it's a command, just parse the attributes and pass to the template
    if ( $child->gi eq 'command' ) {
        push @$elements, { &process_command( $sg_twig, $child ), is_command => 1, pipeline_id => $pipeline_id };
        
    } elsif ( $child->gi eq 'commandSet' ) {
        process_gN_iter_xml( $child );
    }
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

$tmpl->param( ELEMENTS => $elements );
print $tmpl->output;

exit(0);

sub process_gN_iter_xml {
    my $commandSet = shift;
    my $gN_file = $commandSet->first_child('fileName')->text();
    
    if (! -e $gN_file && ! -e "$gN_file.gz" ) {
        ## we're not currently checking this return at all, but it prevents
        #   errors in attempting to parse a file that wasn't created.
        return 0;
    }
    
    my $gN_fh;
    if ($gN_file =~ /\.gz/) {
        open($gN_fh, "<:gzip", "$gN_file") || die "can't read $gN_file: $!"; 
    } elsif ( ! -e $gN_file && -e "$gN_file.gz" ) {
        open($gN_fh, "<:gzip", "$gN_file.gz") || die "can't read $gN_file: $!"; 
    } else {
        open($gN_fh, "<$gN_file") || die "can't read $gN_file: $!";
    }
    
    my $twig = XML::Twig->new( twig_roots => {
                                        'commandSet/commandSet' => \&process_gN_commandSet,
                                   },
                   );
    $twig->parse($gN_fh);
    
    return 1;
}


sub process_gN_commandSet {
    my ($twig, $commandset) = @_;
    my $file = $commandset->first_child('fileName')->text();
    
    my %props = (
                    file    => $file,
                    name    => basename($file, ('.xml.gz', '.xml', '.gz')),
                    id      => $commandset->first_child('id')->text(),
                    state   => $commandset->first_child('state')->text(),
                    message => '',
                    is_command => 0,
                );
    
    
    ## check for any messages
    my $msg_html = '';
    if ( $commandset->has_child('status') ) {
        if ( $commandset->first_child('status')->has_child('message') ) {
            my $msg = $commandset->first_child('status')->has_child('message')->text();
            
            ## don't display it if it is just a 'finished' message
            unless ( $msg =~ /Command set with name.+ finished/i ) {
                $props{message} = "<div class='messageblock'>$msg</div>\n";
            }
        }
    }

    push @$elements, \%props;   

    $twig->purge;
}
