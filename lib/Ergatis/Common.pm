use strict;
use Ergatis::ConfigFile;
use File::Find;
use HTML::Template;

=head1 DESCRIPTION

This is just a collection of commonly used methods.  No classes.

=cut


=head2 get_project_conf_param( repository_root, section, parameter )

=over 4

finds the project.config file within the passed repository root and
extracts the desired parameter, returning it as a string.

=back

=cut
sub get_project_conf_param {
    my ($repository_root, $section, $parameter) = @_;
    
    my $cfg = new Ergatis::ConfigFile( -file => "$repository_root/workflow/project.config" );
    
    return $cfg->val( $section, $parameter );
}

=head2 get_module_path( root, name )

=over 4

recursively searches a path (root) for a given module and returns the full
path to that module.  the name passed can be like 'IdGenerator', 'IdGenerator.pm',
or 'Workflow::IdGenerator'.

=back

=cut
sub get_module_path {
    my ($root, $name) = @_;

    $name =~ s/::/\//g;
    my $path;

    find(
        sub {
            if ( $File::Find::name =~ /$name/ ) {
                $path = $File::Find::name;
            }
        },
        $root
    );
    
    return $path;
}

=head2 get_pipeline_templates( path )

=over 4

searches each directory within the path passed for pipeline templates and returns a data
structure ideal for use with HTML::Template.  The data structure returned looks like:

    [
        { id => 'pipeline_name', 
            path => "$path/pipeline_name",
            has_comment => 0,
            comment => 'some comment here',
            component_count => 0, },
    ]

=back

=cut
sub get_pipeline_templates {
    my $dir = shift;
    my @templates = ();

    if ( -d $dir ) {
        opendir( my $recent_dh, $dir ) || die "can't read template directory: $!";
        while ( my $thing = readdir $recent_dh ) {

            if ( -e "$dir/$thing/pipeline.layout" ) {
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


=head2 get_quick_links( config_ref )

=over 4

reads the ergatis config file and returns a data structure needed by HTML::Template
to display the quick links at the top of most page headers.  the argument passed
must be a Config::IniFiles object from parsing ergatis.ini.  expects this file to
contain a 'quick_links' section with key (label) value (url) pairs.

=back

=cut
sub get_quick_links {
    my $ergatis_cfg = shift;
    my $quick_links = [];
    
    for my $label ( $ergatis_cfg->Parameters('quick_links') ) {
        push @$quick_links, { 
                                label => $label,
                                url => $ergatis_cfg->val('quick_links', $label),
                                is_last => 0
                            };
    }
    
    $$quick_links[-1]{is_last} = 1;
    
    return $quick_links;
}

=head2 print_error_page( %args )

=over 4

throws generic error handling page.  allows you to pass an error message and links for
continuation.  example:

    print_error_page( ergatis_cfg => $ergatis_cfg,
          message => "The pipeline passed couldn't be found ($xml_input).  " .
                     "It may have been deleted or there could be a network (NFS) problem.",
          links => [ 
                        { label => "$project pipeline list", 
                          is_last => 0, 
                          url => "./pipeline_list.cgi?repository_root=$repository_root" },
                        { label => 'try again', 
                          is_last => 1, 
                          url => "./view_pipeline.cgi?instance=$xml_input" },
                   ],
    );

This assumes that you've printed a header already within the calling CGI and that you'll
perform an explicit exit() afterwards (as appropriate).

=back

=cut
sub print_error_page {
    my %args = @_;
    
    my $tmpl = HTML::Template->new( filename => 'templates/error.tmpl',
                                    die_on_bad_params => 0,
                                  );

    $tmpl->param( MESSAGE => $args{message} );    
    $tmpl->param( QUICK_LINKS         => &get_quick_links($args{ergatis_cfg}) );
    $tmpl->param( SUBMENU_LINKS       => $args{links} );

    print $tmpl->output;
}

=head2 url_dir_path( $cgi )

=over 4

returns the path of the current script, including 

=back

=cut
sub url_dir_path {
    my $cgi = shift;
    
    my $full = $cgi->url( -full => 1, -query => 0 );
    $full =~ /(.+?)[^\/]+$/;

    return $1;
}

1==1;
