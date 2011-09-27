#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Date::Manip;
use Ergatis::Common;
use Ergatis::ConfigFile;
use Ergatis::Monitor;
use File::Basename;
use File::stat;
use HTML::Template;
use POSIX;
use XML::LibXML;

my $q = new CGI;

## this template will only be used if we can't find anything.
my $tmpl = HTML::Template->new( filename => 'templates/quick_search.tmpl',
                                die_on_bad_params => 0,
                              );

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

## Setup the information we need to paginate our results
my $view_page_num = $q->param("page_num") || 1;
my $max_pages_to_show = 10;
my @page_links;
my $pipelines_per_page = $ergatis_cfg->val( 'display_settings', 'pipelines_per_page');

## Grab our query term(s). For basis of searching the comments of pipelines we 
## want a unique list of the query terms to avoid over-weighting the results.
my $crit = $q->param('quick_search_crit') || '';
my %terms = map { $_, 1 } split('\s', $crit);
my @query_terms = keys %terms;

my @pipeline_hits = ();
my @search_results = ();
my $total_pipe_hits = 0;

my $parser = XML::LibXML->new();

## search projects:*, paths:default_project_root, pipeline comments, and pipeline IDs for each if it's numeric
for my $project ( $ergatis_cfg->Parameters('projects') ) {
   ## search pipelines
    my $pipeline_dir = $ergatis_cfg->val('projects', $project) . '/workflow/runtime/pipeline';
    if (-d $pipeline_dir) {
        opendir(my $pdh, $pipeline_dir) || die "failed to open pipeline directory $pipeline_dir";
        while (my $dir = readdir $pdh) {
            my $hit_counts = 0;
            my $pipeline_comment_file = "$pipeline_dir/$dir/pipeline.xml.comment";
            my $pipeline_xml_file = "$pipeline_dir/$dir/pipeline.xml";
            my $pipeline_user = "";
            my $last_mod = "";

            if ($dir eq $crit && -e $pipeline_xml_file) {
                print $q->redirect( -uri => url_dir_path($q) . "view_pipeline.cgi?instance=$pipeline_xml_file" );
                exit(0);
            } elsif (-e $pipeline_comment_file && -r $pipeline_comment_file) {
                next if (! -e $pipeline_xml_file );

                ## If we aren't dealing with a pipeline ID go ahead and search through all comments for our pipelines looking for our query terms
                open(PCOMMENT, $pipeline_comment_file) || die "Failed to open pipeline comment $pipeline_comment_file";
                my @comment_lines = <PCOMMENT>;
                my $comment_str = join(" ", @comment_lines);
                close (PCOMMENT);
                
                foreach my $term (@query_terms) {
                    $hit_counts += () = $comment_str =~ /$term/gi;

                    # Highlight the search term in the comment...
                    $comment_str =~ s/($term)/\<b\>$1\<\/b\>/gi;
                }

                ## If our $hit_counts variable is greater than 0 we had a hit in the comment        
                ## and need to pull out the 
                if ($hit_counts > 0) {
                    my $filestat = stat($pipeline_xml_file);
                    $pipeline_user = getpwuid($filestat->uid);
                    $last_mod = $filestat->mtime;

                    push(@pipeline_hits, 
                            [$project, 
                             $dir, 
                             $pipeline_dir, 
                             $pipeline_xml_file, 
                             $comment_str, 
                             $hit_counts,
                             $pipeline_user,
                             $last_mod]);
                }
            }
        }
    }

    ## If we didn't get any hits for a pipeline ID or in the pipeline comments check to see if this is a search for a specific project
    $total_pipe_hits = scalar @pipeline_hits;
    if ( $total_pipe_hits == 0 && lc($project) eq lc($crit) ) {
        ## redirect to project pipeline list and stop
        print $q->redirect( -uri => url_dir_path($q) . "pipeline_list.cgi?repository_root=" . $ergatis_cfg->val('projects', $project) );
        exit(0);
    }
}

## Generate the page links that allow the user to click through different pages of search results.
my ($show_pre_cont, $show_post_cont, $pages_ref) = generate_page_links($view_page_num, $max_pages_to_show, $pipelines_per_page, $total_pipe_hits);
@page_links = @$pages_ref;

## Sort results to move pipelines with the most hits to the top of the group
@pipeline_hits = sort { $b->[5] <=> $a->[5] || $b->[7] <=> $a->[7] } @pipeline_hits;
my @pipeline_range = generate_paginated_pipeline_range($view_page_num, $pipelines_per_page, $total_pipe_hits);

for my $index (@pipeline_range) {
    my @hit_metadata = @{ $pipeline_hits[$index-1] };
    my $metadata = get_pipeline_metadata_libxml($hit_metadata[2], $hit_metadata[1], $hit_metadata[3], $parser);

    push(@search_results, {
            "project"          => $hit_metadata[0],
            "pipeline_id"      => $hit_metadata[1], 
            "pipeline_user"    => $hit_metadata[6],
            "state"            => $metadata->{'state'},
            "is_running"       => $metadata->{'state'} eq "running" ? 1 : 0,
            "run_time"         => $metadata->{'runtime'},
            "components"       => $metadata->{'components'},
            "component_count"  => $metadata->{'component_count'},
            "view_link"        => "./view_pipeline.cgi?instance=$hit_metadata[3]",
            "archive_link"     => "./archive_pipeline_form.cgi?repository_root=$hit_metadata[2]&amp;pipeline_id=$hit_metadata[1]",
            "clone_link"       => "clone_pipeline.cgi?instance=$hit_metadata[3]&repository_root=$hit_metadata[2]",
            "error_message"    => $metadata->{'error_msg'},
            "has_comment"      => length($hit_metadata[4]),
            "pipeline_comment" => $hit_metadata[4],
            "last_mod"         => scalar localtime($hit_metadata[7]),
            "last_mod_raw"     => $hit_metadata[7],
            "links_enabled"    => $metadata->{'links_enabled'},
            "query_hits"       => $hit_metadata[5]
        }
    );
}

## If we didn't have a direct pipeline ID provided or a project name provided 
## we want to attempt to print out the contents of pipeline_hits
print $q->header( -type => 'text/html' );
$tmpl->param( QUICK_LINKS  => get_quick_links($ergatis_cfg) );
$tmpl->param( PIPELINES    => \@search_results);
$tmpl->param( 'NEEDS_PAGINATION' => $total_pipe_hits > $pipelines_per_page ? 1 : 0 );
$tmpl->param( 'PIPELINE_COUNT' => $total_pipe_hits);
$tmpl->param( PAGE_LINKS       => \@page_links ) if (scalar @page_links > 0 );
$tmpl->param( IS_FIRST_PAGE    => $view_page_num == 1 ? 1 : 0 );
$tmpl->param( IS_LAST_PAGE     => ($view_page_num * $pipelines_per_page) > $total_pipe_hits ? 1 : 0 );

## We don't want to pass these to the template if we get 0 results
if (scalar @pipeline_range > 0) {
    $tmpl->param( MIN_PIPELINE_POS => $pipeline_range[0] );
    $tmpl->param( MAX_PIPELINE_POS => $pipeline_range[-1] );
}

$tmpl->param( PREVIOUS_PAGE_URL => "./quick_search.cgi?quick_search_crit=$crit&page_num=" . ($view_page_num - 1) );
$tmpl->param( NEXT_PAGE_URL    =>  "./quick_search.cgi?quick_search_crit=$crit&page_num=" . ($view_page_num + 1) );
$tmpl->param( SHOW_PRE_CONTINUATION => $show_pre_cont);
$tmpl->param( SHOW_POST_CONTINUATION => $show_post_cont);
$tmpl->param('QUERY'       => $crit);
$tmpl->param( RESULT_COUNT => scalar(@search_results) > 0 ? 0 : 1);

print $tmpl->output();
exit(0);

#############################################################
#                     SUBROUTINES                           #
#############################################################

###
# Generates the paginated 'page links' on the search results page
###
sub generate_page_links {
    my ($page_num, $max_pages, $pipes_per_page, $pipe_count) = @_;
    my ($next_lower_interval, $next_higher_interval);
    my @page_links;
    my $page_count = 0;

    if ( $pipe_count > $pipes_per_page ) {
        if ( $pipe_count % $pipes_per_page ) {
            $page_count = int($pipe_count / $pipes_per_page) + 1;
        } else {
            $page_count = int($pipe_count / $pipes_per_page);
        }
    } else {
        $page_count = 1;
    }

    if ($view_page_num == $max_pages_to_show) {
        $next_lower_interval  = ($view_page_num / $max_pages_to_show);
    } elsif ($view_page_num % $max_pages_to_show) {
        $next_lower_interval  = floor($view_page_num/$max_pages_to_show) * $max_pages_to_show + 1;
    } else {
        $next_lower_interval  = $view_page_num - $max_pages_to_show + 1;
    }

    $next_higher_interval = ceil($view_page_num/$max_pages_to_show) * $max_pages_to_show + 1;

    my $show_pre_continuation  = $next_lower_interval > 1 ? 1 : 0;
    my $show_post_continuation = $next_higher_interval < $page_count ? 1 : 0;

    for ( my $i=1; $i<=$page_count; $i++ ) {
        next if $i >= $next_higher_interval || $i < $next_lower_interval;

        push (@page_links, {
            page_num => $i,
            is_active => $i == $view_page_num ? 1 : 0,
            url => "./quick_search.cgi?quick_search_crit=$crit&page_num=$i"
        });
    }

    return ($show_pre_continuation, $show_post_continuation, \@page_links);
}

###
# Generates our range of pipelines that will display in a given page of 
# search results
###
sub generate_paginated_pipeline_range {
    my ($page_num, $pipes_per_page, $hit_counts) = @_;
    my ($min_pipeline_pos, $max_pipeline_pos);

    if ($page_num == 1) {
        $min_pipeline_pos = 1;
    } else {
        $min_pipeline_pos = ($pipes_per_page * ($page_num -1)) + 1;
    }

    if ($page_num * $pipes_per_page >= $hit_counts) {
        $max_pipeline_pos = $hit_counts;
    } else {
        $max_pipeline_pos = $page_num * $pipelines_per_page;
    }

    return ($min_pipeline_pos..$max_pipeline_pos);
}

###
#  Extracts relevant content from a pipeline XML file
###
sub get_pipeline_metadata_libxml {
    my ($project_dir, $pipeline_id, $xml_file, $parser) = @_;
    my $metadata = {};
    my $state = "unknown";
    my $start_time = "unknown";
    my $end_time = "unknown";
    my $run_time = "unknown";
    my $error_message = "";
    my $component_aref = [];
    my $component_count = 0;

    my $links_enabled = 1;

    # This could be a sketchy ploy, but the pipeline ID should be 
    # the name of the folder that house our XML file.            
    my $dirname = dirname($xml_file);       

    my $lock_file = "$project_dir/workflow/lock_files/pipeline." .
                    "$pipeline_id.lock";
    if (-e $lock_file) {
        open (my $lockfh, "<$lock_file") || 
            die("Can't read lock file");
        $state = readline $lockfh;
        close ($lockfh);

        chomp $state;
    
        $run_time = ''; 
        $links_enabled = 0;
    } else {
        # If our pipeline file doesn't exists it might exist in gunzip'd
        # format.
        if (! -e $xml_file) {
            if (-e "$xml_file.gz") {
                $xml_file = "$xml_file.gz";
            } else {
                warn("Could not find pipeline XML for pipeline" .
                              "$pipeline_id");
                next;
            }
        }

        if (! -s $xml_file) {
            warn("Empty pipeline XML for pipeline $pipeline_id");
            next;
        }

        # Make sure that we actually have components to parse 
        my $doc = $parser->parse_file($xml_file);
        my $root = $doc->getDocumentElement;

        next if (! $doc->find('/commandSetRoot/commandSet/state') );
    
        my $commandSet = ($root->findnodes('commandSet'))[0];
        my $state = $commandSet->findvalue('state') || 'unknown';

        ($start_time, $end_time, $run_time) = time_info_libxml($commandSet);

        # Build our component list
        ($component_count, $component_aref) = get_component_list($commandSet);

        ## depending on the state, grab the top-level error message if there is one
        ##  there's no use doing this for running pipelines, since workflow won't
        ##  propagate the error until it's finished.
        if ( $state eq 'error' || $state eq 'failed' ) {
            $error_message = $commandSet->findvalue('status/message') || 0;
            $error_message = CGI::escapeHTML($error_message);
        }

        $metadata->{'error_msg'} = $error_message;
        $metadata->{'state'} = $state;
        $metadata->{'start_time'} = $start_time;
        $metadata->{'end_time'} = $end_time;
        $metadata->{'runtime'} = $run_time;
        $metadata->{'components'} = $component_aref;
        $metadata->{'component_count'} = $component_count . " components";
        $metadata->{'links_enabled'} = $links_enabled;
    }

    return $metadata;
}

###
# Generates a list of all components in this pipeline and a count of each 
# component.
###
sub get_component_list {
    my $command = shift;
    my %components;
    my $component_aref = [];
    my $component_count = 0;

    foreach my $node ($command->findnodes('//commandSet')) {
        my $name = $node->findvalue('name');
        if ($name && $name=~ /^(.+?)\./) {
            $components{$1}{'count'}++;
    
            my $state = $node->findvalue('state');
            if ($state) {
                if ($state eq 'error' || $state eq 'failed') {
                    $components{$1}{'error_count'}++;
                }
            }
        }
    }

    # HTML::Template requires us to pass any complex data structures in an 
    # array of hashes so we need to convert
    foreach my $component (sort keys %components) {
        $component_count++;
        push (@{ $component_aref }, {
            'name'          => $component,
            'count'         => $components{$component}{'count'},
            'error_count'   => $components{$component}{'error_count'}
        });
    }

    return ($component_count, $component_aref);
}

###
# Grabs the relevant time information from our pipeline XML
###
sub time_info_libxml {
    my $command = shift;

    if (! $command->find('startTime') ) {
        return ('unavailable', 'unavailable', 'unavailable');
    }
    
    my $start_time_obj = ParseDate($command->findvalue('startTime'));
    my $start_time = UnixDate($start_time_obj, "%c");

    my ($end_time_obj, $end_time);
    if ( $command->find('endTime') ) {
        $end_time_obj = ParseDate($command->findvalue('endTime'));
        $end_time = UnixDate($end_time_obj, "%c");
    }

    my $runtime = "?";

    if ($start_time_obj) {
        my $diffstring;
        $runtime = '';
        
        if ($end_time_obj) {
            $diffstring = DateCalc($start_time_obj, $end_time_obj);
        } else {
            $diffstring = DateCalc($start_time_obj, "now");
        }
        
       ## take out any non \d: characters
        $diffstring =~ s/[^0-9\:]//g;

        my @parts = split(/:/, $diffstring);
        
        ## years + months + weeks + days
        my $val = ($parts[0] * 365) + ($parts[1] * 30) + ($parts[2] * 7) + ($parts[3]);
        if ($val > 1) {
            $runtime .= "$val days ";
        } elsif ($val == 1) {
            $runtime .= "$val day ";
        }
        
        $runtime .= "$parts[4] hr " if $parts[4];
        $runtime .= "$parts[5] min " if $parts[5];
        $runtime .= "$parts[6] sec";
    }

    $runtime = '&lt; 1 sec' if $runtime eq '0 sec';

    return ($start_time, $end_time, $runtime);
}


###
# Extracts relevant content from a pipeline XML file.
###
sub get_pipeline_metadata_twig {
    my ($project_dir, $pipeline_id, $xml_file, $twig) = @_;
    my $metadata = {};
    my $state = "unknown";
    my $start_time = "unknown";
    my $end_time = "unknown";
    my $run_time = "unknown";
    my $pipeline_user = "";
    my $last_mod = "";
    my $error_message = "";
    my $component_aref = [];
    my $component_count = 0;
    my $links_enabled = 1;

    # This could be a sketchy ploy, but the pipeline ID should be 
    # the name of the folder that house our XML file.            
    my $dirname = dirname($xml_file);       
    
    my $filestat = stat($xml_file);
    $pipeline_user = getpwuid($filestat->uid);
    $last_mod = $filestat->mtime;

    my $lock_file = "$project_dir/workflow/lock_files/pipeline." .
                    "$pipeline_id.lock";
    if (-e $lock_file) {
        open (my $lockfh, "<$lock_file") || 
            die("Can't read lock file");
        $state = readline $lockfh;
        close ($lockfh);

        chomp $state;
    
        $run_time = ''; 
        $links_enabled = 0;
    } else {
        # If our pipeline file doesn't exists it might exist in gunzip'd
        # format.
        if (! -e $xml_file) {
            if (-e "$xml_file.gz") {
                $xml_file = "$xml_file.gz";
            } else {
                warn("Could not find pipeline XML for pipeline" .
                              "$pipeline_id");
                next;
            }
        }

        if (! -s $xml_file) {
            warn("Empty pipeline XML for pipeline $pipeline_id");
            next;
        }

        my $ifh;
        if ($xml_file =~ /\.gz$/) {
            open($ifh, "<:gzip", $xml_file) || 
                die("Couldn't read pipeline XML $xml_file");
        } else {
            open($ifh, "<$xml_file") || 
                die("Couldn't read pipeline XML $xml_file");
        }

        $twig->parse($ifh);
        my $command_set_root = $twig->root;
        my $command_set = $command_set_root->first_child('commandSet');

        next if (! $command_set);

        if ($command_set->has_child('state')) {
            $state = $command_set->first_child('state')->text();
        }

        ($start_time, $end_time, $run_time) = time_info($command_set);

        ## If the state is failed or error we want to grab the top-level error message.
        if ( $state eq 'error' || $state eq 'failed' ) { 
            if ( $command_set->has_child('status') &&  
                $command_set->first_child('status')->has_child('message') ) { 
                $error_message = $command_set->first_child('status')->first_child('message')->text();

                ## handle illegal characters
                $error_message = CGI::escapeHTML($error_message);
            }   
        }   

        ## Grab a list of all our components in the pipeline
        my %pipeline_components = component_count_hash( $xml_file );
        foreach my $component (sort keys %pipeline_components) {
            $component_count += 1;
            push (@{ $component_aref }, { 
                    name        => $component,
                    count       => $pipeline_components{$component}{'count'},
                    error_count => $pipeline_components{$component}{'error_count'},
                }
            );        
        }

        $metadata->{'pipeline_user'} = lc($pipeline_user);
        $metadata->{'error_msg'} = $error_message;
        $metadata->{'last_mod'} = $last_mod;
        $metadata->{'state'} = $state;
        $metadata->{'start_time'} = $start_time;
        $metadata->{'end_time'} = $end_time;
        $metadata->{'runtime'} = $run_time;
        $metadata->{'components'} = $component_aref;
        $metadata->{'component_count'} = $component_count . " components";
        $metadata->{'links_enabled'} = $links_enabled;
    }

    return $metadata;
}
