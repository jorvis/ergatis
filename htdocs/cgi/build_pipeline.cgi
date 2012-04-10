#!/usr/bin/perl -w

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

## Whether or not to enable email notification by default
my $email_on_default = $shared_cfg->val('project', '$;EMAIL_NOTIFICATION$;') || 0;

## Grab the current logged in user and the email domain (if both of these are specificed)
## and craft together our email.
my $email_user = "";
my $email_domain = $ergatis_cfg->val('workflow_settings', 'email_domain');
my $username = user_logged_in($ergatis_cfg);
if ($username && $email_domain && $email_on_default) {
    $email_user = "$username\@$email_domain";
} else {
    ## If the user isn't logged in and default email notifications is set to true 
    ## we want to behave is if it isn't (i.e. make them add in email).
    $email_on_default = 0;
}

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


my $recent_pipelines = get_pipeline_templates( $build_area );
my $project_templates = get_pipeline_templates( "$repository_root/workflow/project_saved_templates" );

my $build_directory = "$build_area/" .temp_pipeline_id();

my $inputs = &get_inputsets($build_area,10);

$tmpl->param( INPUTSETS => $inputs);
$tmpl->param( REPOSITORY_ROOT => $repository_root );
$tmpl->param( WORKFLOWDOCS_DIR => $workflowdocs_dir );
$tmpl->param( COMPONENT_GROUPS => $component_groups );
$tmpl->param( COMPONENT_NAMES => \@components );
$tmpl->param( RECENT_PIPELINES => $recent_pipelines );
$tmpl->param( PROJECT_TEMPLATES => $project_templates );
$tmpl->param( BUILD_DIRECTORY => $build_directory );
$tmpl->param( AUTOLOAD_TEMPLATE => $q->param('autoload_template') || '' );
$tmpl->param( PIPELINE_COMMENT_FILE => "$build_directory/pipeline.xml.comment" );
$tmpl->param( PIPELINE_COMMENT => '' );
$tmpl->param( EMAIL_ON_DEFAULT => $email_on_default );
$tmpl->param( EMAIL_USER => $email_user );
$tmpl->param( BUILDER_ANIMATIONS => $ergatis_cfg->val( 'display_settings', 'builder_animations' ) || 0 );
$tmpl->param( QUICK_LINKS         => &get_quick_links($ergatis_cfg) );
$tmpl->param( SUBMENU_LINKS       => [
                                        { label => 'run pipeline', is_last => 0, url => 'javascript:checkAndRunPipeline()' },
                                        { label => 'save as a project template', is_last => 0, url => 'javascript:showPipelineNameForm(true)' },
                                        { label => "save, don't run", is_last => 1, url => 'javascript:document.pipeline.skip_run.value=1;checkAndRunPipeline()' },
                                     ] );

print $tmpl->output;

exit(0);

sub get_inputsets{
	my $dir = shift;
	my $limit = shift;
	my @inputs;
	my @inputlists = `find $dir -name "*.config"`;
	
	my $uniqueinputs = {};
	foreach my $file (sort {chomp $a;chomp $b;(stat($b))[9]<=>(stat($a))[9]} (@inputlists)){
	    if(scalar(@inputs)<$limit){
		chomp $file;
		my @fstats = stat($file);
		my $component_cfg = new Ergatis::ConfigFile( -file => "$file" );
		
		## skip this config if it is invalid
		next if( !defined( $component_cfg ) );

		my $cname = $component_cfg->val( 'component', '$;COMPONENT_NAME$;' );
		$cname =~ s/\s//g;
		my $token = $component_cfg->val( 'output', '$;OUTPUT_TOKEN$;');
		my $input_file = $component_cfg->val( 'input', '$;INPUT_FILE$;');
		my $input_list = $component_cfg->val( 'input', '$;INPUT_FILE_LIST$;');
		my $input_dir = $component_cfg->val( 'input', '$;INPUT_DIRECTORY$;');
		my $input_type;
		my $input_val;
		if(defined($input_file) && $input_file ne ""){
		    $input_type = "file";
		    $input_val = $input_file;
		}
		elsif(defined($input_list) && $input_list ne ""){
		    $input_type = "list";
		    $input_val = $input_list;
		}
		elsif(defined($input_dir) && $input_dir ne ""){
		    $input_type = "directory";
		    $input_val = $input_dir;
		}
            	next unless( defined( $input_val ) ); # move here from below
		if(!exists $uniqueinputs->{$input_type}->{$input_val}){
            	# next unless( defined( $input_val ) ); # was here, think this was wrong spot
		    push @inputs,{'NAME'=>"$cname.$token (".scalar localtime($fstats[9]).")",
				  'VALUE'=>$input_val,
				  'TYPE'=>$input_type,
				  'SOURCE'=>"$cname.$token"};
		}
		$uniqueinputs->{$input_type}->{$input_val}++;
	    }
	}
	return \@inputs;
}

sub min_group {
    my $counts = shift;
    my $min_num = 0;
    my $min_count = 10000000000000;
    
#    print STDERR "checking current groups counts: @$counts\n";
    
    for ( my $i=0; $i< scalar @$counts; $i++ ) {
        
        if ( $$counts[$i] < $min_count ) {
            $min_count = $$counts[$i];
            $min_num = $i;
        }
    }
    
    return $min_num;
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
