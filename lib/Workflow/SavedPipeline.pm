package Workflow::SavedPipeline;

=head1 NAME

SavedPipeline.pm - A module for loading and building saved pipeline templates.

=head1 SYNOPSIS

    To load an existing pipeline template, and then create a new pipeline
    instance within the AA1 project space:
    
    my $pipe = Workflow::SavedPipeline->new( 
        template => '/path/to/saved/pipelines/some_label/pipeline.xml'
    );
    $pipe->write_pipeline( 
        repository_root => '/usr/local/annotation/FUN' 
    );
    
    To save a pipeline template from an existing pipeline.xml:

    my $pipe = Workflow::SavedPipeline->new(
        source => '/path/PROJECT/Workflow/pipeline/pipeline_id/pipeline.xml'
    );
    $pipe->write_template( template => 'save/pipeline/dir/' );

=head1 DESCRIPTION

=head1 METHODS

=over 3

=item I<PACKAGE>->new()

Returns a newly created "saved" pipeline object.  If loading an existing
template, pass a 'template' argument pointing to the template pipeline.xml

=item I<$OBJ>->load_template()

Loads an existing pipeline template.  This only needs to be called if you
failed to pass a 'template' parameter when using the new method.

=item I<$OBJ>->load_pipeline()

Loads an existing pipeline structure and checks that component.conf.bld.ini 
files are present, usually in order to create a template for later use.
Only needs to be called if a 'source' parameter was not passed with new method.

=item I<$OBJ>->write_template()

Writes a pipeline template into the saved area, 
(optionally with the associated component config files <- not true, all 
component.conf.bld.ini files will be copied into save dir)

=item I<$OBJ>->write_pipeline()

Creates a ready-to-execute pipeline xml (not instance xml) in the defined project
space.  You must pass a repository_root parameter to this method.  You may optionally
pass a shared_config parameter, else it will derive it from the repository root.


=back

=head1 NOTES

Using the same SavedPipeline object for saving a template pipeline and building
a new pipeline from a saved template is not recommended and is prevented.

=head1 AUTHOR

    Joshua Orvis
    jorvis@tigr.org

=cut


use strict;
use Carp;
use Config::IniFiles;
use Cwd;
use File::Basename;
use File::Copy;
BEGIN {
    require '/usr/local/devel/ANNOTATION/ard/chado-v1r5b1/lib/site_perl/5.8.5/Workflow/IdGenerator.pm';
    import Workflow::IdGenerator;
}
use XML::Twig;


## class data and methods
{
    my %_attributes = (
                        pipeline_id             => undef,
                        shared_config           => undef,
                        template                => undef,
		        source                  => undef,
		        _source_twig            => undef,
                        _configmapids           => undef,
                        _template_dir           => undef,
                        _repository_root        => undef,
		        _load_type              => undef,
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
        ## if the source was defined, load it
        if (defined $self->{source}) {
            $self->load_pipeline();
        }
        
        return $self;
    }

    sub load_pipeline {
	my ($self, %args) = @_;

	if (! $self->{_load_type}) {
	    $self->{_load_type} = "pipeline";
	}
	elsif ( $self->{_load_type} ne "pipeline") {
	    croak( "Inappropriate object reuse: cannot load pipelines and load templates with the same object");
	}
        
	# either source must exist from new or they just passed a source
        if (! defined $self->{source} ) {
            if (exists $args{source}) {
                $self->{source} = $args{source};
            } else {
                croak ("source pipeline.xml must be provided either in new or load_pipeline");
            }
        }

	#get repository root and pipeline_id of input xml
	#(this might be incorrect, could result in user overwriting pipeline.xml)
	#(if write pipeline called on same repository_root, pipeline_id would already be stored)
	if ( $self->{source} =~ m|(.+/([^/]+))/Workflow/pipeline/(\d+)| ) {
	    $self->{_repository_root} = $1;
#	    $project_id = $2;
	    $self->{pipeline_id} = $3;
	} else {
	    croak("failed to extract repository root, project id, or pipeline id from ".$self->{source});
	}
	
	#load pipeline from xml
        $self->{_source_twig} = XML::Twig->new(pretty_print => 'indented');
	$self->{_source_twig}->parsefile($self->{source}) || croak "unable to parse ".$self->{source};

	#strip it to the barebones
	$self->_strip_pipeline( $self->{_source_twig}->root );

	#remove version= if present
	$self->{_source_twig}->root->strip_att("version");

	#check that all expected component.conf.bld.ini files are present
	#(not saving configMapIds)
	foreach my $child ( $self->{_source_twig}->root->descendants("configMapId") ) {
	    if ( $child->text =~ m|^component_(.+?)\.(.+)|) {
		my $component = $1;
		my $token = $2;
		#check that conf exists
		my $source_file = $self->{_repository_root}."/Workflow/$component/".$self->{pipeline_id}."_$token/component.conf.bld.ini";
		unless (-e $source_file) {
		    croak("component configuration file $source_file not found");
		}
	    }
	}
    }


    sub write_template {
	my ($self, %args) = @_;
	unless (exists $args{template}) {
	    croak ("template output directory must be provided");
	}

        #print out the pipeline
	unless (-e $args{template}) {
	    mkdir($args{template}) || croak( "unable to create tmeplate directory ".$args{template});
	}
	open(my $FOUT, ">".$args{template}."/pipeline.xml") || croak("unable to create pipeline.xml in ".$args{template});
	$self->{_source_twig}->print($FOUT);
	close($FOUT);

        #copy neccessary .ini files
        #convert from -/PROJECTDIR/Workflow/$component/$pipelineid_$token/component.conf.bld.ini
        #to $template_dir/$component.$token.ini

	#repository_root and pipeline_id should already be defined
	unless (exists $self->{_repository_root} && exists $self->{pipeline_id}) {
	    croak("repository_root or pipeline_id not defined");
	}

	#copy .ini files
	foreach my $child ( $self->{_source_twig}->root->descendants("configMapId") ) {
	    if ( $child->text =~ m|^component_(.+?)\.(.+)|) {
		my $component = $1;
		my $token = $2;
		#copy .ini, setting SHARED_CONFIG = ""
		my $source_file = $self->{_repository_root}."/Workflow/$component/".$self->{pipeline_id}."_$token/component.conf.bld.ini";
		my $dest_file = $args{template}."/$component.$token.ini";
		my $cfg = new Config::IniFiles( -file => $source_file );
		$cfg->setval( "include $component", '$;SHARED_CONFIG$;', "" );        
		$cfg->WriteConfig( $dest_file );
	    }
	}
    }

    #strip all unwanted tags from the pipeline, 
    #leaving only the "core" tags for saving
    sub _strip_pipeline {
	my ($self, $n) = @_; #a tag in the pipeline
	foreach my $child ( $n->children() ) {
#	    $child->gi eq "commandSetRoot" 
	    if ($child->gi eq "commandSet")
	    {
		my $cmi = $child->first_child("configMapId")->text();
		#tags we want to keep and search under
		if ($cmi eq "start" ||
		    $cmi eq "empty" ||
		    $cmi =~ /^component_/ ||
		    $cmi =~ /^pipeline_/ ||
		    $cmi =~ /^\d+/) 
		{
		    $self->_strip_pipeline($child);
		}
		else {
		    $child->delete;
		}
	    }
	    #tags to keep and NOT search under
	    elsif ($child->gi ne "#PCDATA" && #text in the tag is treated as an element
		   $child->gi ne "configMapId") { 
		$child->delete;
	    }
	}
    }



    sub load_template {
        my ($self, %args) = @_;

   	if (! $self->{_load_type}) {
	    $self->{_load_type} = "template";
	}
	elsif ( $self->{_load_type} ne "template") {
	    croak( "Inappropriate object reuse: cannot load pipelines and load templates with the same object");
	}
     
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
    
    sub pipeline_id {
        ## just return the pipeline_id
        return $_[0]->{pipeline_id} || undef;
    }
    
    sub write_pipeline {
        my ($self, %args) = @_;
        
        ## repository root must be passed
        if (! defined $args{repository_root}) {
            croak("repository_root is a required argument to the write_pipeline method");
        }
        
        $self->{_repository_root} = $args{repository_root};
        
        ## was a shared_config passed?
        if (defined $args{shared_config}) {
            $self->{shared_config} = $args{shared_config};
        }

        ## if the shared_config hasn't been defined, derive it
        if (! defined $self->{shared_config} ) {
            $self->{shared_config} = $args{repository_root} . "/workflow_config_files/sharedconf.ini";
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
        $self->_create_dir($pipeline_dir);
        
        ## copy the component INIs
        ##  this copies each file like jaccard.default.ini to
        ##  $repository_root/Workflow/jaccard/$pipeline_id_default/component.bld.conf.ini
        for my $configmapid ( $self->_configmapids() ) {
        
            ## only do the components
            if ($configmapid =~ /component_(.+?)\.(.+)/) {
                my ($component_name, $token) = ($1, $2);
                my $outputdir = "$args{repository_root}/Workflow/$component_name";
                
                ## create the component directory
                $self->_create_dir($outputdir);
                
                ## create the component instance directory
                $outputdir .= '/' . $self->{pipeline_id} . "_$token";
                $self->_create_dir($outputdir);
                
                ## place the component conf into the output directory
                ##  also setting the shared conf
                $self->_place_component_conf( $component_name, $token, $outputdir );
            }
        }
        
        ## create the pipeline XML
        $self->_create_pipeline_xml( $pipeline_dir );
        
        ## create the pipeline INI
        $self->_create_pipeline_ini( $pipeline_dir );
        
        
    } ## end write_pipeline sub
    
    
    sub _create_pipeline_ini {
        my ($self, $pipeline_dir) = @_;
        
        ## create a new Config file
        my $cfg = new Config::IniFiles;
           $cfg->SetFileName( "$pipeline_dir/pipeline.xml.ini" );

        ## add a blank Section for each configMapId
        for my $configmapid ( $self->_configmapids() ) {
            $cfg->AddSection($configmapid);
        }
        
        $cfg->RewriteConfig();
    }

    sub _create_pipeline_xml {
        my ($self, $pipeline_dir) = @_;
        
        ## create a twig and parse through the pipeline template xml
        my $twig = XML::Twig->new( pretty_print => 'indented' );
        $twig->parsefile( $self->{template} );

        ## when you come to a component configMapId, append the new XML after it        
        ## look at each of the configMapIds in the document.
        foreach my $configmapid_elt ( $twig->getElementsByTagName('configMapId') ) {
            my $configmapid =  $configmapid_elt->text();
            
            if ($configmapid =~ /component_(.+?)\.(.+)/) {
                my ($component_name, $token) = ($1, $2);
            
                ## create a twig for the text we need to append and replace
                ##  the template
                my $component_twig = $self->_create_component_twig( $component_name, $token );
                $component_twig->replace( $configmapid_elt->parent );
            }
        }
        
        ## write the twig out to the target directory
        open (my $ofh, ">$pipeline_dir/pipeline.xml") || croak ("can't create $pipeline_dir/pipeline.xml");
        $twig->print($ofh);
    }

    sub _create_component_twig {
        my ($self, $name, $token) = @_;
        my $pipeline_id         = $self->pipeline_id();
        my $bin_dir             = $self->_get_bin_dir();
        my $repository_root     = $self->{_repository_root};
        
        my $xmlfragment = <<XMLfraGMENt;
<commandSet type='serial'>
    <configMapId>component_$name.$token</configMapId>
    <command>
        <type>RunUnixCommand</type>
        <configMapId>generate_component_$name.$token</configMapId>
        <name>Generate component $name.$token</name>
        <param>
            <key>command</key>
            <value>$bin_dir/run_pipeline</value>
        </param>
        <arg>-c $repository_root/Workflow/$name/${pipeline_id}_$token/component.conf.bld.ini --skiprun --pipelineid=$pipeline_id</arg>
    </command>
    <commandSet type="serial" version="2.2">
        <name>$name.$token</name>
        <maxParallelCmds>0</maxParallelCmds>
        <configMapId>run_component_$name.$token</configMapId>
        <fileName>$repository_root/Workflow/$name/${pipeline_id}_$token/pipeline.xml</fileName>
        <state>incomplete</state>
    </commandSet>
</commandSet>
XMLfraGMENt

        my $t = new XML::Twig;
           $t->parse($xmlfragment);

        return $t->root;
    }
    
    sub _get_bin_dir {
        my $self = shift;
        
        my $cfg = new Config::IniFiles( -file => $self->{shared_config} );
        return $cfg->val('init', '$;BIN_DIR$;');
    }
    
    sub _add_configmapid {
        my ($self, $configmapid) = @_;
        
        ## increment a counter for this configmap id
        $self->{_configmapids}{$configmapid}++;
    }
    
    sub _check_template_configmaps {
        my $self = shift;
        
        my $twig = XML::Twig->new(  );
        $twig->parsefile( $self->{template} );
        
        ## look at each of the configMapIds in the document.
        foreach my $configmapid_elt ( $twig->getElementsByTagName('configMapId') ) {
            my $configmapid =  $configmapid_elt->text();

            ## remember the configmapid
            $self->_add_configmapid( $configmapid ); #for the commandset
            
            ## if this is a component, parse the name and check that the file exists.
            if ($configmapid =~ /component_(.+?)\.(.+)/) {
		$self->_add_configmapid( "generate_".$configmapid ); #for command "Generate component"
		$self->_add_configmapid( "run_".$configmapid ); #for commandset that runs subflow 
                my ($component_name, $token) = ($1, $2);

                ## file will be named like:
                ##  component.token.ini
                my $component_file = $self->{_template_dir} . "/$component_name.$token.ini";
                if (! -e $component_file) {
                    croak( "component configuration file $component_file could not be found" );
                }
            }
        }
    }

    sub _configmapids {
        return keys %{$_[0]->{_configmapids}};
    }
    
    sub _create_dir {
        my ($self, $dir) = @_;
        
        if (! -e $dir) {
            mkdir($dir) || croak "failed to create dir $dir";
        }
    }

    sub _place_component_conf {
        my ($self, $component_name, $token, $outputdir) = @_;
        
        ## parse the config file
        my $cfg = Config::IniFiles->new( -file => $self->{_template_dir} . 
                                                  "/$component_name.$token.ini"
                                       );
        
        ## set the shared conf
        $cfg->setval( "include $component_name", '$;SHARED_CONFIG$;', $self->{shared_config} );
        
        ## write it into the proper directory
        $cfg->WriteConfig( "$outputdir/component.conf.bld.ini" );
        
    }

} ## end class
