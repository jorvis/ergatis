use strict;
use Config::IniFiles;
use File::Find;

=head1 DESCRIPTION

This is just a collection of commonly used methods.  No classes.

=cut


=head2 get_project_conf_param( repository_root, section, parameter )

=over 4

finds the sharedconf.ini file within the passed repository root and
extracts the desired parameter, returning it as a string.

=back

=cut
sub get_project_conf_param {
    my ($repository_root, $section, $parameter) = @_;
    
    my $cfg = new Config::IniFiles( -file => "$repository_root/workflow_config_files/sharedconf.ini" );
    
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

=head2 parse_pipeline_run_lock_file ( lock_file_path )

=over 4

reads a pipeline running lock file and returns a hash the values contained within
it, such as pid, hostname, execuser and retry count.

=back

=cut

sub parse_pipeline_run_lock_file {
    my $lock_file = shift;
    
    ## set defaults
    my %parts = ( 
                  pid => '?',
                  hostname => 'unknown',
                  execuser => 'unknown',
                  retries => '?',
                );
    
    if(-e $lock_file){
        open FILE, "$lock_file" or die "Can't open lock file $lock_file";
        my(@elts) = <FILE>;
        close FILE;
        chomp(@elts);
        $parts{pid} = $elts[0];
        $parts{hostname} = $elts[1];
        $parts{execuser} = $elts[2];
        $parts{retries} = $elts[3];
        
        return %parts;
    }
    
    return %parts;
}

1==1;
