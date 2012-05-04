#!/usr/bin/perl -w

=head1 NAME

pipeline_list_json.cgi - A CGI script that returns a listing of all pipelines 
                         (with metadata) for a given project in JSON format

=head1 SYNOPSIS

    pipeline_list_json.cgi 
        repository_root=/path/to/some/project/repository
        page_num=<Desired page of pipelines to view>

=head1 DESCRIPTION

This CGI script mimics the procedures used by the pipeline_list.cgi script to 
grab pipelines for the provided project but returns the information in JSON
format to be consumed and manipulated by the user.

=head1 AUTHOR

    Cesar Arze
    carze@som.umaryland.edu
    
    Joshu Orvis
    jorvis@som.umaryland.edu

=cut 

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Date::Manip;
use Ergatis::Common;
use Ergatis::ConfigFile;
use Ergatis::Monitor;
use POSIX;
use XML::LibXML;
use JSON;

my $q = new CGI;
my $repository_root = $q->param("repository_root") || die "pass a repository root";
my $view_page_num = $q->param("page_num") || 1;
my $max_pages_to_show = 10;

# Want to grab the URL of this script and manipulate it a bit to give us the location
# of the cgi-bin so we can form links to certain pipeline resources
my $cgi_bin_url = $q->url();
$cgi_bin_url =~ s/\/pipeline_list_json.cgi//;

my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );
my $pipeline_root = "$repository_root/workflow/runtime/pipeline";
my $project_conf_path = "$repository_root/workflow/project.config";
my $project_comment_file = "$repository_root/workflow/project.comment";
my $pipelines_per_page = $ergatis_cfg->val( 'display_settings', 'pipelines_per_page');
my $errors_found  = 0;
my $error_msgs = [];

# Conduct some checks on directories and files that will be used when grabbing
# metadata from workflow XML;
if (! -e $project_conf_path ) {
    generate_err_response($q, 400, "$project_conf_path not found");
}

my $project_conf = new Ergatis::ConfigFile( -file => $project_conf_path );
my $ergatis_dir = $project_conf->val('project', '$;ERGATIS_DIR$;');

my $project_code = 0;
if ( $project_conf->val('project', '$;PROJECT_CODE$;') ) {
    $project_code = $project_conf->val('project', '$;PROJECT_CODE$;');
}

# Check for some required variables in the project.config
my @required_vars = ( '$;PROJECT$;', '$;REPOSITORY_ROOT$;', '$;TMP_DIR$;', '$;PROJECT_ID_REPOSITORY$;',
                      '$;ERGATIS_DIR$;', '$;LIB_DIR$;', '$;BIN_DIR$;', '$;DOCS_DIR$;' );
               
for my $var ( @required_vars ) {
    if (! defined $project_conf->val( 'project', $var ) ) {
        generate_error_response($q, 400, "$var not defined in project's configuration file");
    }
}

# Ensure that our required directories exist and are writable
my @required_dirs = ( $pipeline_root, "$repository_root/workflow", "$repository_root/workflow/runtime",
                      "$repository_root/workflow/lock_files", $project_conf->val( 'project', '$;PROJECT_ID_REPOSITORY$;' ) );

for my $dir ( @required_dirs ) {
    if (! -d $dir ) {
        generate_error_response($q, 400, "directory $dir not found");
    }
}

my $rdh;
if (! opendir $rdh, "$pipeline_root" ) {
    generate_error_response($q, 400, "can't open pipeline root directory: $!");
}

# Do some groundwork to setup everything we need to paginate our pipelines.
my %pipeline_quickstats = get_pipeline_quickstats($rdh, $pipeline_root);
my $pipeline_count = scalar keys %pipeline_quickstats;

my @pagination_vars = prepare_pipeline_pagination($view_page_num, $repository_root, $max_pages_to_show, $pipelines_per_page, $pipeline_count);
my $page_count = $pagination_vars[0];
my $min_pipeline_pos = $pagination_vars[1];
my $max_pipeline_pos = $pagination_vars[2];
my $show_pre_continuation = $pagination_vars[3];
my $show_post_continuation = $pagination_vars[4];

my @pipelines_to_parse = get_pipeline_ids_range('pipeline', $min_pipeline_pos, $max_pipeline_pos, $pipelines_per_page, \%pipeline_quickstats);
my $pipelines = parse_pipeline_xmls($pipeline_root, $repository_root, \%pipeline_quickstats, $cgi_bin_url, \@pipelines_to_parse);

## Add a couple project specific pieces of metadata
$pipelines->{'project_code'} = $project_code;
$pipelines->{'current_page'} = int($view_page_num);
$pipelines->{'total_pages'} = $page_count;
$pipelines->{'next_page'} = "$cgi_bin_url/pipeline_list_json.cgi?repository_root=$repository_root&page_num=" . ($view_page_num + 1) if ($view_page_num  < $page_count);
$pipelines->{'previous_page'} = "$cgi_bin_url/pipeline_list_json.cgi?repository_root=$repository_root&page_num=" . ($view_page_num - 1) if ($view_page_num > 1);

print $q->header( -type => 'application/json' );
print encode_json $pipelines;

#############################################################
#                     SUBROUTINES                           #
#############################################################

###
# Generates an HTTP error response to the request
###
sub generate_error_response {
    my ($q, $error_code, $error_msg) = shift;
    print $q->header('status' => $error_code);
    exit 0;
}

###
# Parses all pipeline ID's in the provided array and grab the corresponding 
# pipeline XML file that will be parsed and have its metadata stored in a 
# hash.
###
sub parse_pipeline_xmls {
    my ($pipeline_root, $repository_root, $quickstats, $cgi_url, $pipeline_ids_ref) = @_;
    my @pipeline_ids = @$pipeline_ids_ref;
    my $pipelines_metadata = {};


    my $parser = XML::LibXML->new();

    foreach my $pipeline_id (@pipeline_ids) {
        # Some variables to house metadata we want to track
        my $state = 'unknown';
        my $start_time = 'unknown';
        my $end_time = 'unknown';
        my $run_time = 'unknown';
        my $component_count = 0;
        my $component_aref = [];
        my $error_message = '';
        my $pipeline_comment = '';

        my $pipeline_file = "$pipeline_root/$pipeline_id/pipeline.xml"; 
        my $archive_link = "$cgi_url/archive_pipeline_form.cgi?repository_root=$repository_root&amp;pipeline_id=$pipeline_id";

        # Check if we have a lock file in place for this pipeline. 
        # A lockfile indicates that the pipeline could be in the archive/delete process
        my $lock_file = "$repository_root/workflow/lock_files/pipeline." .
                        "$pipeline_id.lock";
        if (-e $lock_file) {
            open (my $lockfh, "<$lock_file") || 
                die("Can't read lock file");
            $state = readline $lockfh;
            close ($lockfh);

            chomp $state;
        
            $run_time = ''; 
        } else {
            # If our pipeline file doesn't exists it might exist in gunzip'd
            # format.
            if (! -e $pipeline_file) {
                if (-e "$pipeline_file.gz") {
                    $pipeline_file = "$pipeline_file.gz";
                } else {
                    warn("Could not find pipeline XML for pipeline" .
                                  "$pipeline_id");
                    next;
                }
            }

            if (! -s $pipeline_file) {
                warn("Empty pipeline XML for pipeline $pipeline_id");
                next;
            }

            ## is there a comment?
            if ( -e "$pipeline_file.comment" ) {
                open(my $ifh, "<$pipeline_file.comment") || die "can't read comment file: $!";
                while (<$ifh>) {
                    $pipeline_comment .= $_;
                }
            }

            # Make sure that we actually have components to parse 
            my $doc = $parser->parse_file($pipeline_file);
            my $root = $doc->getDocumentElement;

            next if (! $doc->find('/commandSetRoot/commandSet/state') );
        
            my $commandSet = ($root->findnodes('commandSet'))[0];
            $state = $commandSet->findvalue('state') || "unknown";

            ## Grab timing information about our pipeline.
            ## The run time will need to be de-English'd for our JSON output
            ($start_time, $end_time, $run_time) = time_info_libxml($commandSet);
            $run_time = parse_runtime_to_seconds($run_time);

            ($component_count, $component_aref) = get_component_list($commandSet);

            ## depending on the state, grab the top-level error message if there is one
            ##  there's no use doing this for running pipelines, since workflow won't
            ##  propagate the error until it's finished.
            if ( $state eq 'error' || $state eq 'failed' ) {
                $error_message = $commandSet->findvalue('status/message') || "";
            }
            
            push (@{ $pipelines_metadata->{'pipelines'} }, {
                    'state'             => $state,
                    'run_time'          => $run_time,
                    'components'        => $component_aref,
                    'component_count'   => $component_count,
                    'view_link'         => "$cgi_url/view_pipeline.cgi?instance=$pipeline_file",
                    'clone_link'        => "$cgi_url/clone_pipeline.cgi?instance=$pipeline_file&repository_root=$repository_root",
                    'archive_link'      => $archive_link,
                    'error_message'     => $error_message,
                    'pipeline_comment'  => $pipeline_comment,
                    'pipeline_id'       => int($pipeline_id),
                    'last_mod'          => int($quickstats->{$pipeline_id}->{'last_mod'}),
                    'pipeline_user'     => $quickstats->{$pipeline_id}->{'pipeline_user'}
            });
        }
    }

    return $pipelines_metadata;
}

###
# Takes the english-friendly pipeline run time returned from the time_info
# subroutine of the Monitor.pm module and converts it to a seconds represention
###
sub parse_runtime_to_seconds {
    my $raw_run_time = shift;
    my $run_time = 0;

    # Now for some ugly regex-if stuff...
    if ($raw_run_time =~ /(\d+) day(s)?/) {
        $run_time += $1 * 86400;
    }

    if ($raw_run_time =~ /(\d+) hr/) {
        $run_time += $1 * 3600;
    }

    if ($raw_run_time =~ /(\d+) min/) {
        $run_time += $1 * 60;
    }

    if ($raw_run_time =~ /(\d+) sec/) {
        $run_time += $1;
    }

    # If non of the regex's matched we need to set run time to 1 second
    return ($run_time > 0 ? $run_time : 1);
}
