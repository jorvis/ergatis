#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::ConfigFile;
use HTML::Template;

my $q = new CGI;
print $q->header( -type => 'text/html' );

umask(0000);

my $repository_root = $q->param('repository_root') || die "didn't get a repository_root";
my $pipeline_id = $q->param('pipeline_id') || die "didn't get a pipeline_id";
my $action = $q->param('action') || die "didn't get an action";

my $tmpl = HTML::Template->new( filename => 'templates/archive_pipeline.tmpl',
                                die_on_bad_params => 1,
                              );

my $message = '';
my $message_header = '';

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

my $archive_root = $ergatis_cfg->val( 'paths', 'pipeline_archive_root' ) || die "pipeline_archive_root not defined in ergatis.ini file";
my $temp_space   = $ergatis_cfg->val( 'paths', 'temp_space' ) || die "temp_space not defined in ergatis.ini file";

## create the archive root if it doesn't exist
unless ( -d $archive_root ) {
    mkdir( $archive_root ) || die "failed to create $archive_root: $!";
}

if ($action eq 'delete') {
    
    my $log_file = "$temp_space/ergatis.pipeline.$pipeline_id.delete.log";
    
    my $child_pid = fork;
    
    if ($child_pid) {
        $message_header = "deleting pipeline $pipeline_id";
        $message = "<p>Deletion of pipeline $pipeline_id has started and will run in the background. " .
                   "It may take several minutes.</p><p>A log file is being written to $log_file</p>";
        &print_page( $message_header, $message );
        exit;
    } else {
        ## we have to close everything or it won't detach properly
        close STDOUT;
        close STDERR;
        close STDIN;
        
        ## Fork again.  This helps separate the background process from
        ## the httpd process.  If we're in the original child, $gpid will
        ## hold the process id of the "grandchild", and if we're in the
        ## grandchild it will be zero.
        my $gpid = fork;
        if (! $gpid) {
            
            ## We're in the grandchild.
            my $result = `./bin/delete_pipeline.pl --pipeline_id=$pipeline_id --repository_root=$repository_root --delete_output=1 --log=$log_file --lock_file=$repository_root/workflow/lock_files/pipeline.$pipeline_id.lock`;
            exit;
        }
        
        ## don't do anything, we only forked twice to separate from apache
        exit;
        
    }
    
} else {
    $message_header = "archive/deletion failure";
    $message = "unhandled action: $action.  It hasn't been coded yet?";
    &print_page( $message_header, $message );
}

sub print_page {
    my ( $header, $msg ) = @_;
    
    $tmpl->param( MESSAGE_HEADER => $header );
    $tmpl->param( MESSAGE => $msg );

    print $tmpl->output;
}


exit(0);
