package Workflow::SavedPipeline;

=head1 NAME

SavedPipeline.pm - A module for loading and building saved pipeline templates.

=head1 SYNOPSIS

    full script example here showing common methods

=head1 DESCRIPTION

=head1 METHODS

=over 3

=item I<PACKAGE>->new()

Returns a newly created "Saved" pipeline object.  

=item I<$OBJ>->load_template()

Loads an existing pipeline template.

=item I<$OBJ>->load_pipeline()

Loads an existing pipeline structure, usually in order to create a template for
later use.

=item I<$OBJ>->write_template()

Writes a pipeline template into the saved area, optionally with the associated
component config files.

=item I<$OBJ>->write_pipeline()

Creates a ready-to-execute pipeline xml (not instance xml) in the defined project
space.


=back

=head1 AUTHOR

    Joshua Orvis
    jorvis@tigr.org

=cut


use strict;
use Carp;
use Cwd;
use File::Basename;
use Workflow::IdGenerator;
use XML::Twig;


## class data and methods
{
    my %_attributes = (
                        pipeline_id     => undef,
                        template        => undef,
                        _template_dir   => undef,
                      );

    sub new {
        my ($class, %args) = @_;
        
        ## create the object
        my $self = bless { %_attributes }, $class;
        
        ## set any attribute passed
        for (keys %args) {
            if (exists $_attributes{$_}) {
                $self->{$_} = $args{$_};
            } else {
                croak ("$_ is not a recognized attribute");
            }
        }
        
        ## if the template was defined, load it
        if (defined $self->{template}) {
            $self->load_template();
        }
        
        return $self;
    }

    sub load_template {
        my ($self, %args) = @_;
        
        ## if the template wasn't defined in the constructor, it's required now.
        if (! defined $self->{template} ) {
            if (exists $args{template}) {
                $self->{template} = $args{template};
            } else {
                croak ("template must be defined before it can be loaded");
            }
        }
        
        ## make sure the template file exists.
        if (! -e $self->{template}) {
            croak ("template file passed (" . $self->{template} . ") doesn't exist");
        }
        
        ## parse the directory portion of the template
        $self->{_template_dir} = dirname( $self->{template} );
        
        ## check each of the configMapIds to verify that the local ini files
        ##  exist.
        $self->_check_template_configmaps();
        
    }
    
    sub write_pipeline {
        my ($self, %args) = @_;
        
        ## repository root must be passed
        if (! defined $args{repository_root}) {
            croak("repository_root is a required argument to the write_pipeline method");
        }
        
        ## Workflow/pipeline root must exist
        my $pipeline_root = $args{repository_root} . '/Workflow/pipeline';
        if (! -d $pipeline_root) {
            croak( "$pipeline_root does not exist" );
        }
        
        ## if the pipeline_id is not defined, pull one
        if (! defined $self->{pipeline_id}) {
            my $idgen = new Workflow::IdGenerator;
            $self->{pipeline_id} = $idgen->next_id();
        }
        
        ## create the pipeline directory
        my $pipeline_dir = $pipeline_root . '/' . $self->{pipeline_id};
        mkdir($pipeline_dir) || croak("failed to create pipeline directory $pipeline_dir");
        
        TODO write the component INIs
        
        TODO write the pipeline XML and ini
        
    }
    
    sub _check_template_configmap { 
        ## this just checks to see that the config file defined by the
        ##  name in configMapId actually exists.
        my ($t, $elt) = @_;
        
        my $configmapid = $elt->text();
        
        ## if this is a component, parse the name and check that the file exists.
        if ($configmapid =~ /component_(.+?)\.(.+)/) {
            my ($component_name, $token) = ($1, $2);
            
            ## file will be named like:
            ##  component.token.ini
            my $component_file = "$component_name.$token.ini";
            if (! -e $component_file) {
                croak( "component configuration file $component_file could not be found" );
            }
        }
        
    }
    
    sub _check_template_configmaps {
        my $self = shift;
        
        ## we need to change into the template dir, since I'm not sure
        #   how to pass the directory value to the check_template_configmap sub
        my $origdir = cwd();
        chdir($self->{_template_dir});
        
        my $t = XML::Twig->new(
                                twig_roots => {
                                                configMapId => \&_check_template_configmap,
                                              }
                              );
        $t->parsefile( $self->{template} );
        
        ## change back to the original directory
        chdir($origdir);
    }


} ## end class
