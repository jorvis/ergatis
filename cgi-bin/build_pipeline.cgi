#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::ConfigFile;
use Ergatis::SavedPipeline;
use HTML::Template;

my $q = new CGI;
print $q->header( -type => 'text/html' );

my $tmpl = HTML::Template->new( filename => 'templates/build_pipeline.tmpl',
                                die_on_bad_params => 1,
                              );

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );
my $build_area = $ergatis_cfg->val( 'paths', 'pipeline_build_area' ) || die "failed to determine pipeline_build_area";

my $repository_root = $q->param('repository_root') || die "need a repository root";
my $shared_cfg = new Ergatis::ConfigFile( -file => "$repository_root/workflow/project.config" );
my $workflowdocs_dir = $shared_cfg->val( 'project', '$;DOCS_DIR$;' );

## make sure the build area exists
if (! -d $build_area) {
    mkdir($build_area) || die "failed to make build directory $build_area: $!";
}

## read the available components
my @components;
my %classes;
opendir(my $idh, $workflowdocs_dir) || die "can't read component directory ($workflowdocs_dir): $!";
while ( my $thing = readdir $idh ) {
    if ( $thing =~ /(.+).config$/ ) {
        my $component_name = $1;
        
        if ( $ergatis_cfg->component_status($component_name) ne 'disabled' ) {
            my $component_cfg = new Ergatis::ConfigFile( -file => "$workflowdocs_dir/$thing" );
            if ( $component_cfg && $component_cfg->val( 'interface', 'classification' ) ) {
                my $class_string = $component_cfg->val( 'interface', 'classification' );
                
                ## this can be a comma-separated list.  split it up
                for ( split(',', $class_string) ) {
                    ## strip leading/tailing whitespace
                    s/^\s+//g;
                    s/\s+$//g;
                    push @{$classes{$_}}, $component_name;
                }
            } else {
                push @{$classes{unclassified}}, $component_name;
            }
        
            push @components, { name => $component_name };
        }
    }
}

closedir $idh;

## we want the components to be sorted by name
@components = sort {$a->{name} cmp $b->{name}} @components;

my @component_classes = ();
for my $class ( sort keys %classes ) {
    push @component_classes, { class => $class, components => [] };
    for my $component ( @{$classes{$class}} ) {
        push @{$component_classes[-1]->{components}}, { name => $component };
    }
}

## 3 columns (groups)
## if you make changes to the group count you'll need to reflect them in the documentation.css
my $component_groups = [ { group_num => 0, component_sections => [] },
                         { group_num => 1, component_sections => [] },
                         { group_num => 2, component_sections => [] },
                         { group_num => 3, component_sections => [] },
                       ];
my $group_counts = [0,0,0,0];

for my $class ( @component_classes ) {
    my $min_group_num = min_group($group_counts);
    
    ## remember how many are in this column (add a few to make up for the header overhead)
    $$group_counts[$min_group_num] += scalar @{$$class{components}} + 2;
    
    ## sort the components within this section
    my @sorted_components = sort { $$a{name} cmp $$b{name} } @{$$class{components}};
    
    ## now add it
    push @{$$component_groups[$min_group_num]{component_sections}}, { section => $$class{class}, components => \@sorted_components };
}


my $recent_pipelines = get_templates( $build_area );
my $project_templates = get_templates( "$repository_root/workflow/project_saved_templates" );

my $build_directory = "$build_area/" .temp_pipeline_id();

$tmpl->param( REPOSITORY_ROOT => $repository_root );
$tmpl->param( WORKFLOWDOCS_DIR => $workflowdocs_dir );
$tmpl->param( COMPONENT_GROUPS => $component_groups );
#$tmpl->param( COMPONENT_CLASSES => \@component_classes );
$tmpl->param( RECENT_PIPELINES => $recent_pipelines );
$tmpl->param( PROJECT_TEMPLATES => $project_templates );
$tmpl->param( BUILD_DIRECTORY => $build_directory );
$tmpl->param( PIPELINE_COMMENT_FILE => "$build_directory/pipeline.xml.comment" );
$tmpl->param( PIPELINE_COMMENT => '' );
$tmpl->param( BUILDER_ANIMATIONS => $ergatis_cfg->val( 'display_settings', 'builder_animations' ) || 0 );
$tmpl->param( QUICK_LINKS         => &get_quick_links($ergatis_cfg) );
$tmpl->param( SUBMENU_LINKS       => [
                                        { label => 'run pipeline', is_last => 0, url => 'javascript:checkAndRunPipeline()' },
                                        { label => 'save pipeline', is_last => 1, url => 'javascript:document.pipeline.skip_run.value=1;checkAndRunPipeline()' },
                                     ] );

print $tmpl->output;

exit(0);

sub min_group {
    my $counts = shift;
    my $min_num = 0;
    my $min_count = 10000000000000;
    
    print STDERR "checking current groups counts: @$counts\n";
    
    for ( my $i=0; $i< scalar @$counts; $i++ ) {
        
        if ( $$counts[$i] < $min_count ) {
            $min_count = $$counts[$i];
            $min_num = $i;
        }
    }
    
    return $min_num;
}

sub get_templates {
    my $dir = shift;
    my @templates = ();

    if ( -d $dir ) {
        opendir( my $recent_dh, $dir ) || die "can't read build area directory: $!";
        while ( my $thing = readdir $recent_dh ) {
            ## these will all have date names
            if ( $thing =~ /^\d+$/ && -e "$dir/$thing/pipeline.layout" ) {
                push @templates, { id => $thing, 
                                   path => "$dir/$thing",
                                   has_comment => 0,
                                   comment => '',
                                   component_count => 0, };

                if ( -e "$dir/$thing/pipeline.xml.comment" ) {
                    $templates[-1]->{has_comment} = 1;

                    open( my $ifh, "$dir/$thing/pipeline.xml.comment" ) || die "can't read comment file: $!";
                    while ( <$ifh> ) {
                        $templates[-1]->{comment} .= $_;
                    }
                }

                my $layout = Ergatis::SavedPipeline->new( template => "$dir/$thing/pipeline.layout" );

                $templates[-1]->{component_count} = $layout->component_count();
            }
        }
    }
    
    return \@templates;
}


# usage: $string = prettydate( [$time_t] );
# omit parameter for current time/date
sub pretty_date_time {
   my @parts = localtime(shift || time);
   return(sprintf("%04d%02d%02d%02d%02d%02d", $parts[5]+1900, $parts[4]+1, $parts[3], $parts[2], $parts[1], $parts[0]));
} 

sub temp_pipeline_id {
    ## currently the date/time with a number between 10 and 100
    return pretty_date_time() . (int(rand(90)) + 10);
}
