package Ergatis::SavedPipeline;

=head1 NAME

    SavedPipeline.pm - A module for loading and building saved pipeline templates.

    =head1 SYNOPSIS

    To load an existing pipeline template, and then create a new pipeline
    instance within the AA1 project space:
    
    my $pipe = Ergatis::SavedPipeline->new( 
    template => '/path/to/saved/pipelines/some_label/pipeline.layout'
    );
$pipe->write_pipeline( 
    repository_root => '/usr/local/annotation/FUN',
    id_repository => '/path/to/global_id_repository'
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
use Ergatis::ConfigFile;


## class data and methods
{
    my %_attributes = (
                       pipeline_id             => undef,
                       pipeline_token             => undef,
                       shared_config           => undef,
                       template                => undef,
                       source                  => undef,
                       instance_path           => undef,  # instance xml
                       _component_count        => 0,
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

        $self->{_pipeline_dir} = "workflow/runtime/pipeline";
        $self->{_component_dir} = "workflow/runtime";
        $self->{_project_conf} = "workflow/project.config";

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

        return $self;
    }

    ## some accessors
    sub instance_path { return $_[0]->{instance_path} }

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
        if ( $self->{source} =~ m|(.+/([^/]+))/$self->{_pipeline_dir}/([A-Z0-9]+)| ) {
            $self->{_repository_root} = $1;
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
            if ( $child->text =~ m|^(.+?)\.(.+)|) {
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
            mkdir($args{template}) || croak( "unable to create template directory ".$args{template});
        }
        open(my $FOUT, ">".$args{template}."/pipeline.xml") || croak("unable to create pipeline.xml in ".$args{template});
        $self->{_source_twig}->print($FOUT);
        close($FOUT);

        #copy neccessary .ini files
        #convert from -/PROJECTDIR/Workflow/$component/$pipelineid_$token/$component_name.$token.user.config
        #to $template_dir/$component.$token.config

        #repository_root and pipeline_id should already be defined
        unless (exists $self->{_repository_root} && exists $self->{pipeline_id}) {
            croak("repository_root or pipeline_id not defined");
        }

        #copy .ini files
        foreach my $child ( $self->{_source_twig}->root->descendants("name") ) {
            if ( $child->text =~ m|^(.+?)\.(.+)|) {
                my $component = $1;
                my $token = $2;
                #copy .ini, setting SHARED_CONFIG = ""
                my $source_file = $self->{_repository_root}."/$self->{_component_dir}/$component/".$self->{pipeline_id}."_$token/$component.$token.user.config";
                my $dest_file = $args{template}."/$component.$token.config";
                my $cfg = new Config::IniFiles( -file => $source_file );
                my $ret = $cfg->setval("include", '$;PROJECT_CONFIG$;', "" );        
                if(!$ret){
                    croak("can't set \$;PROJECT_CONFIG\$; in section [include]");
                }
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
                if ($cmi =~ /^start pipeline/ ||
                    $cmi eq "empty" ||
                    $cmi eq "serial" ||
                    $cmi eq "parallel" ||
                    $cmi =~ /^\S+\.\S+/ ||  ## matches components
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

        ## preparse the template file and resolve any includes
        my $s;
        $s .= ["a".."z"]->[rand(26)] foreach(1..10);
        my $tmpdir = "/tmp/saved_pipeline_${$}_$s";
        my $t = $self->_resolve_template_includes( $self->{template}, $tmpdir);
        if( $t ne $self->{template} ) {
            $self->{template} = $t;
            $self->load_template();
        } else {
            rmdir( $tmpdir );
        }

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
            my $idgen = new Ergatis::IdGenerator( id_repository => $args{id_repository} );
            $self->{pipeline_id} = $idgen->next_id( type => 'pipeline' );
	        $self->{_logger}->debug("Retrieved pipeline id $self->{pipeline_id}");
        }

        
        ## create the pipeline directory
        my $pipeline_dir = $pipeline_root . '/' . $self->{pipeline_id};
        $self->{_logger}->debug("Creating pipeline dir $pipeline_dir");
        $self->_create_dir($pipeline_dir);
        
        ## copy the layout into the pipeline directory
        copy( $self->{template}, "$pipeline_dir/pipeline.layout" );
        
        ## copy the component INIs
        ##  this copies each file like jaccard.default.conf to
        ##  $repository_root/workflow/runtime/jaccard/$pipeline_id_default/component.token.user.config
        for my $commandname ( $self->_commandnames() ) {
            $self->{_logger}->debug("Creating config for $commandname");        
            ## only do the components
            if ($commandname =~ /(.+?)\.(.+)/) {
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
                $self->_place_component_conf( $component_name, $token, "$outputdir/$component_name.$token.user.config" );
                $self->_place_component_conf( $component_name, $token, "$pipeline_dir/$component_name.$token.config" );
                $self->{user_config} = "$outputdir/$component_name.$token.user.config";
            }
        }
        
        ## create the pipeline XML
        $self->{instance_path} = $self->_create_pipeline_xml( $pipeline_dir );
        return Ergatis::Pipeline->new( id => $self->{pipeline_id}, path => $self->{instance_path} );
        
    } ## end write_pipeline sub
    
    sub _resolve_template_includes {
        my ($self, $template, $tmpdir, $parent_file) = @_;
        
        my $in;
        open($in, "< $template") or croak("Could not open $template: $!");
        my @include_lines = grep( /\<INCLUDE/, <$in> );
        close($in);
        return $template unless( @include_lines );

        my $infh;
        open($infh, "< $template") or croak("Could not open $template: $!");
        my $new_template;
        if( !defined( $parent_file ) ) {
            my $basename = basename( $template );
            $new_template = "$tmpdir/$basename";
            mkdir( $tmpdir ) unless( -d $tmpdir );
        } else {
            my $basename = basename( $parent_file );
            $new_template = "$tmpdir/$basename.child";
        }

        # copy all the configs in the current templates directory to
        # tmpdir
        my $templ_dir = dirname( $template );
        copy( $_, $tmpdir ) foreach( glob("$templ_dir/*config") );

        #where is the other template file in question? In case
        #the include path is relative
        my $temp_basedir = dirname( $template );

        my $outfh;
        open($outfh, "> $new_template") or croak("Could not open $new_template for writing: $!"); 

        #Print out the new pipeline.layout file.
        while( my $line = <$infh> ) {
            chomp $line;
            if( $line =~ /\<INCLUDE.*file\s*=\s*[\"\'](\S+)[\"\']/ ) {
                my $i_file = $1;

                #In case of relative directory
                $i_file = $temp_basedir."/".$i_file unless( $i_file =~ m|^/| );

                my $i_base = basename( $i_file );

                croak("Include file $i_file does not exist") unless( -e $i_file );
                my $i_dir = dirname( $i_file );

                copy($_, $tmpdir) foreach( glob("$i_dir/*config") );
                
                ## Skip printing the commandSetRoot and also skip the outer
                ## command set (name: start pipeline:). This should be 
                ## defined in the original layout.
                my $child_template = $self->_resolve_template_includes( $i_file, $tmpdir, $template );
                my $twig = new XML::Twig( twig_handlers => {
                    'commandSetRoot' => sub {
                        my ($t, $e) = @_;
                        my $start_pipeline_commandSet = $e->first_child('commandSet');
                        croak("Could not find start pipeline commandSet in included pipeline")
                            unless( $start_pipeline_commandSet );
                        foreach my $cs ( $start_pipeline_commandSet->children('commandSet') ) {
                            $cs->print( $outfh );
                        }	
                    }
                });
                $twig->parsefile( $child_template );

                unlink($child_template) if( $child_template ne $i_file );
            } else {
                print $outfh "$line\n";
            }
        }
        close($infh);
        close($outfh);
        return $new_template;
    }  # end _resolve_template_includes

    sub _create_pipeline_xml {
        my ($self, $pipeline_dir) = @_;
        
        ## create a twig and parse through the pipeline template xml
        my $twig = XML::Twig->new( pretty_print => 'indented' );
        $twig->parsefile( $self->{template} );

        ## when you come to a component commandname, append the new XML after it        
        ## look at each of the commandnames in the document.
        foreach my $commandname_elt ( $twig->getElementsByTagName('name') ) {
            my $commandname =  $commandname_elt->text();
            
            if ($commandname =~ /(.+?)\.(.+)/) {
                my ($component_name, $token) = ($1, $2);
                
                ## create a twig for the text we need to append and replace
                ##  the template
                my $component_twig = $self->_create_component_twig( $component_name, $token );
                $component_twig->replace( $commandname_elt->parent );
            }
        }
        foreach my $commandname_elt ( $twig->getElementsByTagName('commandSet') ) {
            if(!$commandname_elt->has_children("name")){
                my $name= XML::Twig::Elt->new('name'=>$commandname_elt->att("type"));
                $name->paste( first_child=>$commandname_elt);
            }
            elsif($commandname_elt->has_children("name")->text() eq 'start'){
                my $name= XML::Twig::Elt->new('name'=>"start pipeline:$self->{pipeline_token}");
                $name->replace( $commandname_elt->has_children("name") );
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
        if(! -d $bin_dir){
            $self->{_logger}->logdie("Can't find directory $bin_dir\n");
        }
        my $repository_root     = $self->{_repository_root};

        
        
        my $xmlfragment;
        my $component_twig_file = $self->_get_component_twig();
        open FILE,"$component_twig_file" or die "Can't open file $component_twig_file";
        my @lines;
        while(my $line=<FILE>){
            if($line !~ /^\<\!\-\-/){
                push @lines,$line;
            }
        }
        close FILE;
        $xmlfragment = join('',@lines);
        
        my $tokens  = {'\$;COMPONENT_DIR\$;'=>"$repository_root/$self->{_component_dir}/$component_name/${pipeline_id}_$token",
                       '\$;COMPONENT_INSTANCE\$;'=>"$component_name.$token",
                       '\$;PIPELINE_DIR\$;' =>"$repository_root/$self->{_pipeline_dir}/$self->{pipeline_id}",
                       '\$;PIPELINE_ID\$;' => "$pipeline_id",
                       '\$;COMPONENT_NAME\$;' => "$component_name",
                       '\$;BIN_DIR\$;'=>"$bin_dir"};
        
        foreach my $key (keys %$tokens){
            $xmlfragment =~ s/$key/$tokens->{$key}/eg;
        }

        my $t = new XML::Twig;
        $t->parse($xmlfragment);

        return $t->root;
    }
    
    sub _get_component_twig{
        my $self = shift;
        
        my $cfg = new Config::IniFiles( -file => $self->{user_config} );
        my $scfg = new Config::IniFiles( -file => $self->{shared_config} );
        my $twigxml;
        if($cfg->val('component', '$;COMPONENT_TWIG_XML$;')){
            $twigxml = $scfg->val('project', '$;DOCS_DIR$;').'/'.$cfg->val('component', '$;COMPONENT_TWIG_XML$;')
            }
        else{
            #default component XML
            $twigxml = $scfg->val('project', '$;DOCS_DIR$;').'/component_template.xml';
        }
        return $twigxml;
    }

    sub _get_bin_dir {
        my $self = shift;
        
        my $cfg = new Config::IniFiles( -file => $self->{shared_config} );
        return $cfg->val('project', '$;BIN_DIR$;');
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
            if ($commandname =~ /(.+?)\.(.+)/) {
                $self->{_component_count}++;
                
                my ($component_name, $token) = ($1, $2);
                
                $self->{_logger}->debug("Found component name $component_name, $token");
                
                ## file will be named like:
                ##  component.token.config
                my $component_file = $self->{_template_dir} . "/$component_name.$token.config";
                
                $self->{_logger}->debug("Looking for file $component_file");
                if (! -e $component_file) {
                    croak( "component configuration file $component_file could not be found" );
                }
            }
        }

        if (scalar(keys %{$self->{_commandnames}}) == 0) {
            $self->{_logger}->warn("No command names found in file $self->{template}");
        }
    }

    sub _commandnames {
        return keys %{$_[0]->{_commandnames}};
    }
    
    sub component_count {
        return $_[0]->{_component_count};
    }
    
    sub _create_dir {
        my ($self, $dir) = @_;
        
        if (! -e $dir) {
            mkdir($dir) || croak "failed to create dir $dir";
        }
    }

    sub _place_component_conf {
        my ($self, $component_name, $token, $output) = @_;
        
        ## parse the config file
        my $cfg = Config::IniFiles->new( -file => $self->{_template_dir} . 
                                         "/$component_name.$token.config"
                                         );
        
        if ( ! defined $cfg ) {
            croak("failed to parse expected ini file: $self->{_template_dir}/$component_name.$token.config with errors:\n".join("\n",@Config::IniFiles::errors));
        }
        
        ## set the shared conf
        my $ret = $cfg->setval( "include", '$;PROJECT_CONFIG$;', $self->{shared_config} );

	    if (!$ret) {
	        croak("can't set \$;PROJECT_CONFIG\$; in section [include] $self->{shared_config} $self->{_template_dir}/$component_name.$token.config");
	    }

	    if ( defined $self->{pipeline_token} && $self->{pipeline_token} ne '') {
	        my $ret = $cfg->setval("component",'$;PIPELINE_TOKEN$;');
	        if(!$ret){
		        $ret = $cfg->newval("component",'$;PIPELINE_TOKEN$;');
		        if(!$ret){
		            croak("can't set \$;PIPELINE_TOKEN\$; in section [component]");
		        }
	        }
	    }
        
        $self->{_logger}->debug("Writing to config file $output");
        ## write it into the proper directory
        $cfg->WriteConfig( "$output" );
    }

    sub configure_saved_pipeline {
        my ($self, $config_file, $repository_root, $id_repository) = @_;

        #Create the new pipeline
        $self->write_pipeline
            ( 'repository_root' => $repository_root,
              'id_repository'   => $id_repository );
        my $pipeline_id = $self->pipeline_id();

        #Replace the internal variables of the input config file
        my $cfg = new Ergatis::ConfigFile( -file => $config_file );

        #Replace some of the variables from the global section.
        my %global_params = map( ($_, $cfg->val('global',$_)), $cfg->Parameters( 'global' ));
        foreach my $section ( $cfg->Sections() ) {
            foreach my $param( $cfg->Parameters( $section ) ) {
                my $val = $cfg->val( $section, $param );
                while( $val =~ /(\$;[\w_]+?\$;)/g ) {
                    if( exists $global_params{$1} ) {
                        my $replacement = $global_params{$1};
                        my $match = $1;
                        $match =~ s/\$/\\\$/g;
                        $val =~ s/$match/$replacement/g;
                        $cfg->setval( $section, $param, $val );
                    }
                }                
            }
        }
        

        #Retrieve all the sections of the input config file
        my @sections = $cfg->Sections( );

        #Check each sections for variables
        foreach my $section ( @sections ) {
            #These are sections that are only used for internal purposes
            next if( $section =~ /global/);

            #Format of section headers contaims the component name and the output token
            my ($component_name, $output_token) = split(/\s+/, $section);

            #This is where the newly written config file should be
            my $ini_file = "$repository_root/workflow/runtime/$component_name/${pipeline_id}_$output_token/$component_name.$output_token.user.config";
            croak("could not find ini_file $ini_file") unless( -e $ini_file );
            
            #The component ini file object
            my $comp_ini = new Ergatis::ConfigFile( -file => $ini_file );

            #Get the parameters fro this section
            my @parameters = $cfg->Parameters( $section );
            foreach my $param ( @parameters ) {

                #Search through each section in the component ini file and replace where necesary.
                #Dies if it can't find the variable ( ie. it doesn't exist, possible typo ).
                my $retval;
                foreach my $comp_section ( $comp_ini->Sections() ) {
                    $retval = $comp_ini->setval( $comp_section, $param, $cfg->val($section, $param) );
                    last if( $retval );
                }
                croak("could not find parameter $param in $ini_file") unless( $retval );
                
            }

            $comp_ini->RewriteConfig();
        }
    } ## end configure_saved_pipeline

} ## end class
