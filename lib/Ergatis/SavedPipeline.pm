package Ergatis::SavedPipeline;

=head1 NAME

SavedPipeline.pm - A module for loading and building saved pipeline templates.

=head1 SYNOPSIS

    To load an existing pipeline template, and then create a new pipeline
    instance within the AA1 project space:
    
    my $pipe = Ergatis::SavedPipeline->new( 
        template => '/path/to/saved/pipelines/some_label/pipeline.xml'
    );
    $pipe->write_pipeline( 
        repository_root => '/usr/local/annotation/FUN' 
    );
    
    To save a pipeline template from an existing pipeline.xml:

    my $pipe = Ergatis::SavedPipeline->new(
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

Loads an existing pipeline structure and checks that $component_name.$token.user.config
files are present, usually in order to create a template for later use.
Only needs to be called if a 'source' parameter was not passed with new method.

=item I<$OBJ>->write_template()

Writes a pipeline template into the saved area, 
(optionally with the associated component config files <- not true, all 
$component_name.$token.user.config files will be copied into save dir)

=item I<$OBJ>->write_pipeline()

Creates a ready-to-execute pipeline xml in the defined project
space.  You must pass a repository_root parameter to this method.  You may optionally
pass a shared_config parameter, else it will derive it from the repository root.
Returns an Ergatis::Pipeline object.


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
use Ergatis::IdGenerator;
use XML::Twig;
use Ergatis::Logger;
use Ergatis::Pipeline;


## class data and methods
{
    my %_attributes = (
                        pipeline_id             => undef,
                        shared_config           => undef,
                        template                => undef,
                        source                  => undef,
                        instance_path           => undef,  # instance xml
                        _source_twig            => undef,
                        _commandnames           => undef,
                        _template_dir           => undef,
                        _repository_root        => undef,
                        _load_type              => undef,
                      );

    sub new {
        my ($class, %args) = @_;
        
        ## create the object
        my $self = bless { %_attributes }, $class;
        $self->{_logger} =  Ergatis::Logger::get_logger("SavedPipeline.pm");
        ## set any attribute passed
        for (keys %args) {
            if (exists $_attributes{$_}) {
                $self->{$_} = $args{$_};
            } else {
                croak ("$_ is not a recognized attribute");
            }
        }
        
        ## if the template was defined, load it
        ## if the source was defined, load it
        if (defined $self->{template}) {
	    $self->{_logger}->debug("Loading template $self->{template}");
            $self->load_template();
        }
        elsif (defined $self->{source}) {
            $self->load_pipeline();
        }
	else{
	    $self->{_logger}->warn("No template or pipeline specified");
	}
        $self->{_pipeline_dir} = "workflow/runtime/pipeline";
        $self->{_component_dir} = "workflow/runtime";
        $self->{_project_conf} = "workflow/project.config";

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
	if ( $self->{source} =~ m|(.+/([^/]+))/$self->{_pipeline_dir}/(\d+)| ) {
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

	#check that all expected $component_name.$token.user.config files are present
	#(not saving names)
	foreach my $child ( $self->{_source_twig}->root->descendants("name") ) {
	    if ( $child->text =~ m|^component_(.+?)\.(.+)|) {
		my $component = $1;
		my $token = $2;
		#check that conf exists
		my $source_file = $self->{_repository_root}."/$self->{_component_dir}/$component/".$self->{pipeline_id}."_$token/$component.$token.user.config";
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
        #convert from -/PROJECTDIR/Workflow/$component/$pipelineid_$token/$component_name.$token.user.config
        #to $template_dir/$component.$token.ini

	#repository_root and pipeline_id should already be defined
	unless (exists $self->{_repository_root} && exists $self->{pipeline_id}) {
	    croak("repository_root or pipeline_id not defined");
	}

	#copy .ini files
	foreach my $child ( $self->{_source_twig}->root->descendants("name") ) {
	    if ( $child->text =~ m|^component_(.+?)\.(.+)|) {
		my $component = $1;
		my $token = $2;
		#copy .ini, setting SHARED_CONFIG = ""
		my $source_file = $self->{_repository_root}."/$self->{_component_dir}/$component/".$self->{pipeline_id}."_$token/$component.$token.user.config";
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
		my $cmi = $child->first_child("name")->text();
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
		   $child->gi ne "name") { 
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
        $self->{_logger}->debug("Parsed directory name $self->{_template_dir}");
        ## check each of the commandnames to verify that the local ini files
        ##  exist.
        $self->_check_template_commandnames();
        
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
            $self->{shared_config} = $args{repository_root} . "/$self->{_project_conf}";
        }
        
        ## Workflow/pipeline root must exist
        my $pipeline_root = $args{repository_root} . "/$self->{_pipeline_dir}";
        if (! -d $pipeline_root) {
            croak( "$pipeline_root does not exist" );
        }
        
        ## if the pipeline_id is not defined, pull one
        if (! defined $self->{pipeline_id}) {
            my $idgen = new Ergatis::IdGenerator( id_repository => $args{id_repository) );
            $self->{pipeline_id} = $idgen->next_id( type => 'pipeline' );
	        $self->{_logger}->debug("Retrieved pipeline id $self->{pipeline_id}");
        }

        
        ## create the pipeline directory
        my $pipeline_dir = $pipeline_root . '/' . $self->{pipeline_id};
	$self->{_logger}->debug("Creating pipeline dir $pipeline_dir");
        $self->_create_dir($pipeline_dir);
        
        ## copy the component INIs
        ##  this copies each file like jaccard.default.ini to
        ##  $repository_root/Workflow/jaccard/$pipeline_id_default/component.bld.conf.ini
        for my $commandname ( $self->_commandnames() ) {
	    $self->{_logger}->debug("Creating config for $commandname");        
            ## only do the components
            if ($commandname =~ /component_(.+?)\.(.+)/) {
                my ($component_name, $token) = ($1, $2);
                my $outputdir = "$args{repository_root}/$self->{_component_dir}/$component_name";
		$self->{_logger}->debug("Writing to parent outputdir $outputdir");        
                ## create the component directory
                $self->_create_dir($outputdir);
                
                ## create the component instance directory
                $outputdir .= '/' . $self->{pipeline_id} . "_$token";
		$self->{_logger}->debug("Writing to outputdir $outputdir");        
                $self->_create_dir($outputdir);
                
                ## place the component conf into the output directory
                ##  also setting the shared conf
                $self->_place_component_conf( $component_name, $token, $outputdir );
            }
        }
        
        ## create the pipeline XML
        $self->{instance_path} = $self->_create_pipeline_xml( $pipeline_dir );
        return Ergatis::Pipeline->new( id => $self->{pipeline_id}, path => $self->{instance_path} );
        
    } ## end write_pipeline sub
    
    
    sub _create_pipeline_xml {
        my ($self, $pipeline_dir) = @_;
        
        ## create a twig and parse through the pipeline template xml
        my $twig = XML::Twig->new( pretty_print => 'indented' );
        $twig->parsefile( $self->{template} );

        ## when you come to a component commandname, append the new XML after it        
        ## look at each of the commandnames in the document.
        foreach my $commandname_elt ( $twig->getElementsByTagName('name') ) {
            my $commandname =  $commandname_elt->text();
            
            if ($commandname =~ /component_(.+?)\.(.+)/) {
                my ($component_name, $token) = ($1, $2);
            
                ## create a twig for the text we need to append and replace
                ##  the template
                my $component_twig = $self->_create_component_twig( $component_name, $token );
                $component_twig->replace( $commandname_elt->parent );
            }
        }
        
        ## write the twig out to the target directory
        open (my $ofh, ">$pipeline_dir/pipeline.xml") || croak ("can't create $pipeline_dir/pipeline.xml");
        $twig->print($ofh);
	close $ofh;
	return "$pipeline_dir/pipeline.xml";
    }

    sub _create_component_twig {
        my ($self, $component_name, $token) = @_;
        my $pipeline_id         = $self->pipeline_id();
        my $bin_dir             = $self->_get_bin_dir();
        my $repository_root     = $self->{_repository_root};
        
        my $xmlfragment = <<XMLfraGMENt;
  <commandSet type="serial">
    <state>incomplete</state>
    <name>$component_name.$token</name>
    <command>
      <type>RunJavaUnixCommand</type>
      <name>replace_config_keys</name>
      <state>incomplete</state>
      <executable>$bin_dir/replace_config_keys</executable>
      <param>  
	<key>--template_conf</key>
	<value>$repository_root/$self->{_component_dir}/$component_name/${pipeline_id}_$token/$component_name.$token.user.config</value>
      </param>
      <param>  
	<key>--output_conf</key>
	<value>$repository_root/$self->{_component_dir}/$component_name/${pipeline_id}_$token/$component_name.$token.final.config</value>
      </param>   
      <param>
	<key>--keys</key>
        <value>PIPELINEID=$pipeline_id</value>
      </param>
    </command>
    <command>
      <type>RunJavaUnixCommand</type>
      <name>replace_template_keys</name>
      <state>incomplete</state>
      <executable>$bin_dir/replace_template_keys</executable>
      <param>  
	<key>--component_conf</key>
	<value>$repository_root/$self->{_component_dir}/$component_name/${pipeline_id}_$token/$component_name.$token.final.config</value>
      </param> 
      <param>  
	<key>--template_xml_conf_key</key>
	<value>TEMPLATE_XML</value>
      </param>   
      <param>  
	<key>--output_xml</key>
	<value>$repository_root/$self->{_component_dir}/$component_name/${pipeline_id}_$token/component.xml</value>
      </param>   
    </command>
    <commandSet type="serial">
      <name>$component_name</name>
      <state>incomplete</state>
      <fileName>$repository_root/$self->{_component_dir}/$component_name/${pipeline_id}_$token/component.xml</fileName>
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
    
    sub _add_commandname {
        my ($self, $commandname) = @_;
        
        ## increment a counter for this commandname id
        $self->{_commandnames}{$commandname}++;
    }
    
    sub _check_template_commandnames {
        my $self = shift;
        
        my $twig = XML::Twig->new(  );
        $twig->parsefile( $self->{template} );
        $self->{_logger}->debug("Parsing xml file $self->{template}");
        ## look at each of the commandnames in the document.
        foreach my $commandname_elt ( $twig->getElementsByTagName('name') ) {
            my $commandname =  $commandname_elt->text();

	    $self->{_logger}->debug("Found command name $commandname");
            ## remember the commandname
            $self->_add_commandname( $commandname ); #for the commandset
            
            ## if this is a component, parse the name and check that the file exists.
            if ($commandname =~ /component_(.+?)\.(.+)/) {
                my ($component_name, $token) = ($1, $2);
		$self->{_logger}->debug("Found component name $component_name, $token");
                ## file will be named like:
                ##  component.token.ini
                my $component_file = $self->{_template_dir} . "/$component_name.$token.config";
		$self->{_logger}->debug("Looking for file $component_file");
                if (! -e $component_file) {
                    croak( "component configuration file $component_file could not be found" );
                }
            }
        }

	if(scalar(keys %{$self->{_commandnames}}) == 0){
	    $self->{_logger}->warn("No command names found in file $self->{template}");
	}
    }

    sub _commandnames {
        return keys %{$_[0]->{_commandnames}};
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
                                                  "/$component_name.$token.config"
                                       );
        
        ## set the shared conf
        $cfg->setval( "include", '$;SHARED_CONFIG$;', $self->{shared_config} );
	$self->{_logger}->debug("Writing to config file $outputdir/$component_name.$token.user.config");        
        ## write it into the proper directory
        $cfg->WriteConfig( "$outputdir/$component_name.$token.user.config" );
        
    }

} ## end class
