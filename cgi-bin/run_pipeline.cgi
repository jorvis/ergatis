#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::ConfigFile;
use Ergatis::SavedPipeline;
use XML::Writer;

my $q = new CGI;

print STDERR "very beginning of run_pipeline.cgi\n";

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => 'ergatis.ini' );
my $global_id_repository = $ergatis_cfg->val('paths', 'global_id_repository') || die "global_id_repository not found in ergatis.ini";

## fetch a hashref of all the parameters, since we'll potentially be
##  querying a lot of them
my $qvars = $q->Vars;

my $build_pipeline_path = "$$qvars{build_directory}/pipeline.xml";

## create a skeleton XML file
open(my $ofh, ">$build_pipeline_path") || die "failed to create pipeline at $build_pipeline_path: $!";

print STDERR "step1\n";

my $writer = new XML::Writer( OUTPUT => $ofh,
                              DATA_MODE => 1,
                              DATA_INDENT => 4,
                            );

$writer->startTag('commandSetRoot', 
                        'xmlns:xsi' => 'http://www.w3.org/2001/XMLSchema-instance',
                        'xsi:schemaLocation' => 'commandSet.xsd' );

$writer->startTag('commandSet', type => 'serial');
$writer->startTag('state');
$writer->characters( 'incomplete' );
$writer->endTag('state');
$writer->startTag('name');
$writer->characters('start');
$writer->endTag('name');

print STDERR "step2\n";

my $next_sibling = $$qvars{'pipeline_root_down'};

while ( $next_sibling ne 'pipeline_root_panel' ) {
    ## is this a component?
    if ( $next_sibling =~ /^c\d+/ ) {
        my $name_token = $$qvars{$next_sibling . '_name.token'};
        
        ## value here will be name.token
        if ( $name_token && $name_token =~ /.+?\..+/ ) {
            $writer->startTag('commandSet', type => 'serial');
            $writer->startTag('state');
            $writer->characters( 'incomplete' );
            $writer->endTag('state');
            $writer->startTag('name');
            $writer->characters( "component_$name_token" );
            $writer->endTag('name');
            $writer->endTag('commandSet');
            
        } else {
            die "failed to parse name and token for component $next_sibling";
        }
   
    ## is it a panel?
    } elsif ( $next_sibling =~ /_panel/ ) {
        ## easy, just close a commandSet
        $writer->endTag('commandSet');
    
    ## is it a set?
    } elsif ( $next_sibling =~ /s\d+/ ) {
        ## open a commandSet of the correct type (serial|parallel)
        $writer->startTag('commandSet', type => $$qvars{$next_sibling . '_type'});
        $writer->startTag('state');
        $writer->characters( "incomplete" );
        $writer->endTag('state');
    
    } else {
        die "panic!  I don't know what to do with $next_sibling";
    }
    
    $next_sibling = $$qvars{ $next_sibling . '_down' };
}

print STDERR "step3\n";

$writer->endTag('commandSet'); ## root 'start' command set
$writer->endTag('commandSetRoot');
$writer->end;

## instantiate the pipeline from the template
print STDERR "attempting to read template $build_pipeline_path to repository root ($$qvars{repository_root})\n";
my $pipeline_template = new Ergatis::SavedPipeline( template => $build_pipeline_path );
print STDERR "step4\n";
my $pipeline = $pipeline_template->write_pipeline( repository_root => $$qvars{repository_root},
                                                   id_repository => $global_id_repository );
print STDERR "step5\n";
print STDERR "about to call the run method\n";
$pipeline->run( ergatis_cfg => $ergatis_cfg );

## now redirect to a monitor page
#print $q->redirect( -uri => url_dir_path($q) . "view_pipeline.cgi?instance=" . $pipeline->path() );

exit(0);

sub process_set {
    my $set_id = shift;
    
    my $next_sibling = $$qvars{$set_id . '_down'};
    
    while ( $next_sibling ) {
        
        if ( $next_sibling =~ /_panel/ ) {
            $next_sibling = 0;
        }
    }
}

## returns the path of the current script, including 
sub url_dir_path {
    my $cgi = shift;
    
    my $full = $cgi->url( -full => 1, -query => 0 );
    $full =~ /(.+?)[^\/]+$/;

    return $1;
}
