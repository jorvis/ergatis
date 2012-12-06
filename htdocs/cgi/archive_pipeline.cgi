#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::ConfigFile;
use HTML::Template;
use File::Path;

my $q = new CGI;
print $q->header( -type => 'text/html' );

umask(0000);

my $repository_root = $q->param('repository_root') || die "didn't get a repository_root";
my $pipeline_id     = $q->param('pipeline_id') || die "didn't get a pipeline_id";
my $action          = $q->param('action') || die "didn't get an action";
my $process_output  = $q->param('process_output') || 0;

my $tmpl = HTML::Template->new( filename => 'templates/archive_pipeline.tmpl',
                                die_on_bad_params => 1,
                              );

my $message = '';
my $message_header = '';

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

my $archive_root = $ergatis_cfg->val( 'paths', 'pipeline_archive_root' ) || die "pipeline_archive_root not defined in ergatis.ini file";
my $temp_space   = $ergatis_cfg->val( 'paths', 'temp_space' ) || die "temp_space not defined in ergatis.ini file";
my $log_file;

my $username = user_logged_in($ergatis_cfg);
unless ($username) {
    print_error_page( ergatis_cfg => $ergatis_cfg,
                      message => "You must be logged in to archive or delete pipelines",
                      links => [
                                 { label => "pipeline list", is_last => 1, url => "./pipeline_list.cgi?repository_root=$repository_root" },
                               ]
                    );
    exit(0);
}

## create the archive root if it doesn't exist
unless ( -d $archive_root ) {
    mkpath( $archive_root ) || die "failed to create $archive_root: $!";
}

if ($action eq 'delete') {
    
    $log_file = "$temp_space/ergatis.pipeline.$pipeline_id.delete.log";

    $message_header = "deleting pipeline $pipeline_id";
    $message = "<p>Deletion of pipeline $pipeline_id has started and will run in the background. " .
               "It may take several minutes.</p><p>A log file is being written to $log_file</p>";
    
    &fork_and_action('delete_pipeline.pl', $message_header, $message, $log_file);

} elsif ($action eq 'archive_in_place') {

    $log_file = "$temp_space/ergatis.pipeline.$pipeline_id.archive_in_place.log";

    $message_header = "archiving pipeline $pipeline_id";
    $message = "<p>Archival of pipeline $pipeline_id has started and will run in the background. " .
               "It may take several minutes.</p><p>A log file is being written to $log_file</p>";
    
    &fork_and_action('archive_pipeline_in_place.pl', $message_header, $message, $log_file);

} elsif ( $action eq 'archive_to_location' ) {

    $log_file = "$temp_space/ergatis.pipeline.$pipeline_id.archive_to_location.log";

    $message_header = "archiving pipeline $pipeline_id";
    $message = "<p>Archival of pipeline $pipeline_id has started and will run in the background. " .
               "It may take several minutes.</p><p>A log file is being written to $log_file</p>";
    
    &fork_and_action('archive_pipeline_to_location.pl', $message_header, $message, $log_file);

}

sub fork_and_action {
    my ($script, $header, $message, $log) = @_;
    
    my $child_pid = fork;
    
    if ($child_pid) {
        &print_page( $header, $message );
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
            if ( $action eq 'archive_to_location' ) {
                `./bin/$script --pipeline_id=$pipeline_id --repository_root=$repository_root -o=$process_output --log=$log --lock_file=$repository_root/workflow/lock_files/pipeline.$pipeline_id.lock --archive_root=$archive_root`;
            } else {
                `./bin/$script --pipeline_id=$pipeline_id --repository_root=$repository_root -o=$process_output --log=$log --lock_file=$repository_root/workflow/lock_files/pipeline.$pipeline_id.lock`;
            }
            exit;
        }
        
        ## don't do anything, we only forked twice to separate from apache
        exit;
        
    }
}

sub print_page {
    my ( $header, $msg ) = @_;
    my $pipeline_list_url = "./pipeline_list.cgi?repository_root=$repository_root";

#    $tmpl->param( REPOSITORY_ROOT => $repository_root );
    $tmpl->param( MESSAGE_HEADER => $header );
    $tmpl->param( MESSAGE => $msg );

    $tmpl->param( QUICK_LINKS         => &get_quick_links($ergatis_cfg) );
    $tmpl->param( PIPELINE_LIST_URL   => $pipeline_list_url );
    $tmpl->param( SUBMENU_LINKS       => [
                                            { label => 'pipeline list', is_last => 1, url => $pipeline_list_url },
                                         ] );

    print $tmpl->output;
}


exit(0);
