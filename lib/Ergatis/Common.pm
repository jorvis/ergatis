use strict;
use CGI;
use CGI::Cookie;
use CGI::Session;
use Cwd qw(abs_path);
use Ergatis::ConfigFile;
use File::Find;
use File::Basename;
use HTML::Template;
use Storable;

=head1 DESCRIPTION

This is just a collection of commonly used methods.  No classes.

=cut

=head2 get_conditional_read_fh( file_path )

=over 4

Accepts a file path as an argument and returns an open, read-only file 
handle for it.  If the file path ends in .gz the file handle is returned
via the gzip IO layer.  It also automatically checks for the existence
of a .gz form of the file you pass, so if you pass 'foo.fsa' and only
'foo.fsa.gz' exists, you'll still get the correct filehandle.

Returns undef if the file and .gz version both don't exist.

=back

=cut

sub get_conditional_read_fh {
    my $path = shift;

    ## auto-handle gzipped files    
    ## if neither version exists return undef
    if (! -e $path) {
        if (-e "$path.gz") {
            $path .= '.gz';
        } else {
            return undef;
        }
    }

    my $fh;
    if ( $path =~ /\.gz$/ ) {
        open($fh, "<:gzip", $path) || die "can't read file $path: $!";
    } else {
        open($fh, "<$path") || die "can't read file $path: $!";
    }
    
    return $fh;
}

=head2 get_conditional_write_fh( file_path )

=over 4

Accepts a file path as an argument and returns an open, writeable file 
handle for it.  If the file path ends in .gz the file handle is returned
via the gzip IO layer.

=back

=cut

sub get_conditional_write_fh {
    my $path = shift;
    my $fh;
    
    if ( $path =~ /\.gz$/ ) {
        open( $fh, ">:gzip", "$path" ) || die "failed to create output file: $!";
    } else {
        open( $fh, ">$path" ) || die "failed to create output file: $!";
    }
    
    return $fh;
}


=head2 parser_referer_url( url )

=over 4

Sanitizes a URL that will be redirected too if needed.

=back

=cut
sub parse_referer_url {
    my $referer = shift;
    my $ret_url = $referer;   

    if ($referer =~ /\/cgi\/(\w+\.cgi.*$)/) {
        $ret_url = $1;
    } 

    return $ret_url;
}

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

Searches each directory within the path passed for pipeline templates and returns a data
structure ideal for use with HTML::Template.  The data structure returned looks like:

    [
        { id => 'pipeline_name', 
            path => "$path/pipeline_name",
            has_comment => 0,
            comment => 'some comment here',
            has_build_guide => 1,
            component_count => 0, },
    ]

If there is any problem parsing the template it is skipped and excluded from the
data structure returned.  This could happen if the pipelines are stored in a temp
directory where some of the files are automatically removed.

=back

=cut
sub get_pipeline_templates {
    my $dir = shift;
    my @templates = ();

    if ( -d $dir ) {
        opendir( my $recent_dh, $dir ) || die "can't read template directory: $!";
        while ( my $thing = readdir $recent_dh ) {
            
            if ( -e "$dir/$thing/pipeline.layout" ) {
                my $layout;
                
                eval { $layout = Ergatis::SavedPipeline->new( template => "$dir/$thing/pipeline.layout" ) };
                
                if (! $@ ) {
                    push @templates, { id => $thing, 
                                       path => "$dir/$thing",
                                       has_comment => 0,
                                       comment => '',
                                       has_build_guide => 0,
                                       component_count => $layout->component_count, };

                    if ( -e "$dir/$thing/pipeline.xml.comment" ) {
                        $templates[-1]->{has_comment} = 1;

                        open( my $ifh, "$dir/$thing/pipeline.xml.comment" ) || die "can't read comment file: $!";
                        while ( <$ifh> ) {
                            $templates[-1]->{comment} .= $_;
                        }
                    }
                    
                    if ( -e "$dir/$thing/pipeline.config" ) {
                        $templates[-1]->{has_build_guide} = 1;
                    }
                    
                } else {
                    print STDERR "failed to parse pipeline template so excluding it from the lists: $@";
                }
            }
        }
    }
    
    return \@templates;
}


=head2 user_logged_in( )

=over 4

Checks if a user is currently logged in - Currently cookie-based..  Returns the
user name or undef.

=back

=cut

sub user_logged_in {
    my $ergatis_cfg = shift;
    my %cookies = fetch CGI::Cookie;
    my $username = undef;
    my $sid = undef;

    if ($cookies{'ergatis_user'}) {
        $sid = $cookies{'ergatis_user'}->value;
    } 

    if ($sid) {
        my $session = get_session($ergatis_cfg, $sid);
        if ($session) {
            $username = $session->param('username') || undef;
        }
    }

    return $username;
}

=head2 get_session($ergatis_cfg)

=over 4

Retrieves the session for the given user. If a session does not exist a 
session will be created.

=back

=cut
sub get_session {
    my ($ergatis_cfg, $sid) = @_;
    my $session = undef;

    if ( $ergatis_cfg->val('authentication', 'session_db_dir') ) {
        my $session_dir = $ergatis_cfg->val('authentication', 'session_db_dir');
        $session = new CGI::Session("driver:File",
                                    $sid,
                                    { Directory => $session_dir }
        );
    }

    return $session;
}

=head2 is_admin_user($ergatis_cfg)

=over 4

If an admin user has been defined in the ergatis configuration this subroutine 
will check whether or not the current logged in user is the admin user 

=back

=cut
sub is_admin_user {
    my $ergatis_cfg = shift;
    my %admin_users = ();
    my $is_admin = 0;

    if ( $ergatis_cfg->val('authentication', 'admin_users') ) {
        my $username = user_logged_in($ergatis_cfg);
        my $admin_users_str = $ergatis_cfg->val('authentication', 'admin_users');
        %admin_users = map { $_ => 1 } split(',', $admin_users_str);
    
        $is_admin = 1 if ( exists($admin_users{$username}) );
    }

    return $is_admin;    
}

=head2 validate_user_authorization($ergatis_cfg, $project, $pipeline_id)

=over 4

If per-account pipeline security is enabled via the ergatis.ini config
file this subroutine will authorize a user to view only pipelines (and 
any components of a pipeline) that the user has created himself/herself.

If a user attempts to view a pipeline (or any components of a pipeline)
that they did not create they will be redirected to an error page indicating
they do not have proper authorization.

=back

=cut
sub validate_user_authorization {                                                                                                                                                                                                       
    my ($ergatis_cfg, $pipeline_id, $defer_exit) = @_;                                                                                                                                                        
    my $authorized_user = 1;                                                                                                                                                                                                             

    ## Authorization should only proceed if we have per-account security                                                                                                                                                                
    ## enabled                                                                                                                                                                                                                           
    if ( $ergatis_cfg->val('authentication','per_account_pipeline_security') && ! is_admin_user($ergatis_cfg) ) {
        my $pipelines = get_account_pipelines($ergatis_cfg);                                                                                                                                                               
                                                                                                                                                                                                                                         
        unless( exists($pipelines->{$pipeline_id}) ) {                                                                                                                                                                                   
            $authorized_user = 0;                                                                                                                                                                                                       
            ## At time we may not want to juse print an error page and exit but                                                                                                                                                          
            ## return a boolean as to whether or not the user has permission to                                                                                                                                                           
            ## view a resource.                                                                                                                                                                                                          
            unless ($defer_exit) {                                                                                                                                                                                                       
                my $username = user_logged_in($ergatis_cfg);
                print_error_page( ergatis_cfg => $ergatis_cfg,                                                                                                                                                                          
                                  message => "User $username does not have authorization to view this resource.",
                                  links => [ { label => "previous page", is_last => 1, url => $ENV{'HTTP_REFERER'} } ]
                                );                                                                                                                                                                                                        
                exit(0);                                                                                                                                                                                                                  
            }                                                                                                                                                                                                                             
        }                                                                                                                                                                                                                                 
    }                                                                                                                                                                                                                                
                                                                                                                                                                                                                                          
    return $authorized_user;                                                                                                                                                                                                             
} 

=head2 get_account_pipelines($ergatis_cfg)

=over 4

Retrieves all pipelines tied to a given account. If a CGI::Session session 
cannot be retrieved or the per_account_pipelines_list flag is not toggled
in the ergatis.ini file an empty hash will be returned.

If the last modified time on the account pipeline list file is newer than
the one recorded in the session metadata the sessions pipeline list will be
refreshed from file.

=back

=cut
sub get_account_pipelines {
    my $ergatis_cfg = shift;
    my $account_pipelines = {};
    my %cookies = fetch CGI::Cookie;
    my $session = undef;

    if ( defined($cookies{'ergatis_user'}) ) {
        my $sid = $cookies{'ergatis_user'}->value;
        $session = get_session($ergatis_cfg, $sid);
    }

    if ( $session && $session->param('username') &&
        $ergatis_cfg->val('authentication', 'per_account_pipeline_security') ) {
        ## Check if our pipelines list file has updated since the last time we updated
        ## our sessions pipeline list and reload if needed
        my $pipelines_file = $ergatis_cfg->val('authentication', 'session_db_dir') .
            "/" . $session->param('username') . ".pipelines";
        my $last_mod = (stat ($pipelines_file))[9] || undef;
        
        if ($session->param('pipelines_last_mod') && $last_mod &&
            $session->param('pipelines_last_mod') < $last_mod) {
            
            $account_pipelines = retrieve($pipelines_file);
            $session->param('pipelines', $account_pipelines);
            $session->param('pipelines_last_mod', $last_mod);

            $session->flush();
        } else {
            $account_pipelines = $session->param('pipelines') || {};
        }
    }

    return $account_pipelines;
}

=head2 add_pipeline_to_user_pipeline_list( $ergatis_cfg, $pipeline_id )

=over 4

Adds a pipeline to the users pipeline list if the per-account pipeline
list flag is enabled in the ergatis.ini

=back

=cut
sub add_pipeline_to_user_pipeline_list {
    my ($ergatis_cfg, $pipeline_id) = @_;
    my %cookies = fetch CGI::Cookie;
    my $session = undef;
    my $pipelines = {};

    if ( defined($cookies{'ergatis_user'}) ) {
        my $sid = $cookies{'ergatis_user'}->value;
        $session = get_session($ergatis_cfg, $sid);
    }

    ## TODO: Figure out whether we need to use file-locking here.
    if ( $session && $session->param('username') ) {
        my $session_db_dir = $ergatis_cfg->val('authentication', 'session_db_dir');
        my $pipelines_file = $session_db_dir . "/" . $session->param('username') .
                             ".pipelines";

        $pipelines = retrieve($pipelines_file) if (-e $pipelines_file);

        ## If this pipeline already exists in our list we are attempting to 
        ## re-run it or something has gone wrong. Either way issue a warning
        ## into our apache log
        if ( exists($pipelines->{$pipeline_id}) ) {
            print STDERR "Duplicate pipeline detected in user pipeline list " .
                         " - $pipeline_id";
        } else {
            ## Add the pipeline ID to our pipeline list.
            $pipelines->{$pipeline_id} = 1;
        }

        $session->param('pipelines', $pipelines);
        store $pipelines, $pipelines_file;                                     

        $session->param('pipelines_last_mod', (stat ($pipelines_file))[9]);
        $session->flush();
    }
}

=head2 delete_pipeline_from_user_pipeline_list( $ergatis_cfg, $pipeline_id )

=over 4

Deletes a pipeline from the users pipeline list.

=back

=cut
sub delete_pipeline_from_user_pipeline_list {
    my ($ergatis_cfg, $pipeline_id) = @_;
    my %cookies = fetch CGI::Cookie;
    my $session = undef;
    my $pipelines = {};

    if ( defined($cookies{'ergatis_user'}) ) {
        my $sid = $cookies{'ergatis_user'}->value;
        $session = get_session($ergatis_cfg, $sid);
    }


    if ( $session && $session->param('username') &&
        $ergatis_cfg->val('authentication', 'per_account_pipeline_security') ) {
        my $session_db_dir = $ergatis_cfg->val('authentication', 'session_db_dir');
        my $pipelines_file = $session_db_dir . "/" . $session->param('username') .
                             ".pipelines";

        $pipelines = retrieve($pipelines_file) if (-e $pipelines_file);
        if ( exists($pipelines->{$pipeline_id}) ) {
            print STDERR "Deleting pipeline ID $pipeline_id from user pipeline list\n";
            delete $pipelines->{$pipeline_id};
        }

        $session->param('pipelines', $pipelines);
        store $pipelines, $pipelines_file;                                     

        $session->param('pipelines_last_mod', (stat ($pipelines_file))[9]);
        $session->flush();
    }
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
    
    ## is the user logged in?
    my $current_user = &user_logged_in($ergatis_cfg);
    
    if ( $current_user ) {

        ## display the logout button
        push @$quick_links, {  
                                label => $current_user . " [logout]",
                                url => './logout.cgi',
                                is_last => 0
                            };
    
    } else {
    
        ## display the login button
        if ( $ergatis_cfg->val('authentication', 'authentication_method') ne 'open' ) {
            push @$quick_links, {  
                                    label => 'login',
                                    url => './login_form.cgi',
                                    is_last => 0
                                };
        }
    }
    
    ## look at the [quick_links] of the ergatis.ini file
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
    my $templates_dir = "templates/";   

    if ( exists($args{templates_dir}) ) {
        $templates_dir = $args{templates_dir};
    } 

    my $tmpl = HTML::Template->new( filename => "$templates_dir/error.tmpl",
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

=head2 run_system_cmd($cmd)

=over 4

Runs a system command wrapped in an eval to catch any issues that 
could occur during the process. If an error occurs this script will die
spitting out an error.

=back

=cut
sub run_system_cmd {
    my ($cmd, $redirect_on_error) = @_;
    my @cmd_output;

    eval {
        @cmd_output = qx{$cmd 2>&1};
        if ( ($? << 8) != 0 ) { 
            die "@cmd_output";
        } 
    };  
    if ($@) {	
	if ($redirect_on_error) {		
            print $redirect_on_error->{'cgi'}->header( -type => 'text/html' );

            my $ergatis_cfg = new Ergatis::ConfigFile( -file => dirname( abs_path($0) ) . "/ergatis.ini" );
            print_error_page( ergatis_cfg => $ergatis_cfg,
                              message => "Command $cmd failed to execute:<br /><br /> <i>$@</i>",
                              links => [ 
                                            { label => $redirect_on_error->{'label'}, is_last => 1, url => $redirect_on_error->{'url'} },
                                       ],  
                              templates_dir => dirname( abs_path($0) ) . "/templates/",
                            );  
            exit(0);
	} else {
        	die "Error executing command $cmd: $@";
	}
    }   
}

1==1;
