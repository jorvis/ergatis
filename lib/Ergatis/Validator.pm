package Ergatis::Validator;

=head1 NAME

Ergatis::Validator.pm - A class for validating a workflow pipeline.

=head1 SYNOPSIS

    use Ergatis::Validator;

    my $validator = new Ergatis::Validator;
    
    ## $warnings is an array reference of warning messages
    my $warnings = $validator->check_pipeline_template( pipeline => '/path/to/pipeline.xml' );

=head1 DESCRIPTION

This class can be used to check a pipeline for possible errors before it
is run.  The list of things currently checked is given in the documentation
for the check_pipeline_template() method.  No changes are made to the 
pipeline xml or any of the component configuration files.

Much more could be done here.  Right now variable replacement is performed
on the configuration files to the extent possible, but some variables are
not yet populated in the conf files, such as PIPELINEID.  Any values containing
a non-replaced variable are not currently checked.

=head1 METHODS

=over 3

=item I<PACKAGE>->new()

Returns a newly created "Ergatis::Validator" object.  No arguments are necessary.

=item I<$OBJ>->check_pipeline_template( pipeline => I<'/path/to/some_pipeline.xml'> )

This is the only public method of this class.  It is used to perform a series of
checks on a pipeline template to catch things that may prevent it from running correctly.
Current checks include:

=over 6

=item

Two components of the same type (wu-blastp) cannot share the same output
token.

=item

If any OUTPUT_DIRECTORY given already exists and has > 0 non-hidden files a
warning is given with the file count.

=item

Any INPUT_FILE passed must exist and have non-zero size.

=item

Any INPUT_LIST passed must exist and have non-zero size.

=item

Any INPUT_DIRECTORY passed must exist and contain at least 1 file with
the same extension given in INPUT_EXTENSION.

=back

=back

=head1 AUTHOR

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;
use warnings;
use Carp;
use Ergatis::ConfigFile;
use XML::Twig;

## class data and methods
{

    my %_attributes = (
                        _inputs  => {},     # $_inputs->{input} = [ components.token ]
                        _input_dirs => {},  # $_input_dirs->{dir}->{ext} = [ components.token ]
                        _output_dirs => {}, # $_output_dirs->{dir} = [ components.token ]
                        _tokens => {},      # $_tokens->{component}->{token} = count;
                      );
    
    sub new {
        my ($class, %args) = @_;
        
        ## create the object
        my $self = bless { %_attributes }, $class;
        
        ## set any attributes passed, checking to make sure they
        ##  were all valid
        for (keys %args) {
            if (exists $_attributes{$_}) {
                $self->{$_} = $args{$_};
            } else {
                croak ("$_ is not a recognized attribute");
            }
        }
        
        return $self;
    }


    sub check_pipeline_template {
        my ($self, %opts) = @_;
        
        ## check required arguments
        unless (exists $opts{pipeline}) { croak("pipeline is a required option to check_pipeline_template") }

        ## have to load configs of all components in the pipeline to gather
        ##  acceptable inputs and outputs
       
        my $counter = 0;
        my $twig = XML::Twig->new( twig_roots => {  
                    commandSet => sub {
                            my $component_conf = $self->_get_component_conf_ini(@_);
                            $self->_parse_component_conf_ini( $component_conf ) if $component_conf;
                        }
                } 
            );

        $twig->parsefile( $opts{pipeline} );
        
        ##########################
        ## perform the checks now
        my $warnings = [];
        
        ## check that no component/output_token combo was used twice.
        for my $component (keys %{$self->{_tokens}} ) {
        
            for my $output_token ( keys %{ $self->{_tokens}->{$component} } ) {
                if ( $self->{_tokens}->{$component}->{$output_token} > 1 ) {
                    push @$warnings, "The output token '$output_token' was used more than once " .
                                    "for $component components.  No two components of the same " .
                                    "type can have the same output token within the same pipeline.";
                }
            }
        }
        
        ## check each output directory if it doesn't have a run-time variable within its name.  if
        ##  it already exists, warn if it contains any files.
        for my $dir ( keys %{ $self->{_output_dirs} } ) {
            ## does it have a variable in it?
            if ( $dir !~ /\$\;/ ) {
                ## does it exist?
                if ( -d $dir ) {
                    my $count = $self->_directory_content_count($dir);
                    if ( $count ) {
                        push @$warnings, "The output directory $dir specified in component(s) (" .
                                         "@{$self->{_output_dirs}->{$dir}}) already contains $count objects";
                    }
                }
            }
        }
        
        ## check that each input was either previously defined as an output or output_dir, else
        ##  must exist and have non-zero size.
        for my $file (  keys %{ $self->{_inputs} }  ) {
            ## does it have a variable in it?  if so, skip it.
            if ( $file !~ /\$\;/ ) {
                if ( -e $file ) {
                    ## is it non-zero?
                    if (! -s $file ) {
                        push @$warnings, "The input $file specified in component(s) (" .
                        "@{$self->{_inputs}->{$file}}) exists but has zero length";
                    }
                } else {
                    push @$warnings, "The input $file specified in component(s) (" .
                    "@{$self->{_inputs}->{$file}}) does not exist";
                }
            }
        }
        
        ## check that each input directory exists and has >0 files of the specified extension
        for my $dir ( keys %{ $self->{_input_dirs} } ) {
            ## skip the dir if it has a variable in its name.
            if ( $dir !~ /\$\;/ ) {

                for my $ext ( keys %{ $self->{_input_dirs}{$dir} } ) {
                    ## warn if there are no files with the passed extension.
                    if (! $self->_directory_content_count($dir, $ext) ) {
                        push @$warnings, "The input directory $dir specified in component(s) (" .
                        "@{$self->{_input_dirs}->{$dir}->{$ext}}) does not contain any files with the $ext extension";
                    }
                }
            }
        }
        
        return $warnings;
    }
    
    ## this private subroutine accepts a twig and commandSet element that contains
    ##  a call to create a pipeline component.  The reference to that component conf
    ##  is extracted and returned.
    sub _get_component_conf_ini {
        my ($self, $twig, $commandSet) = @_;
        if ( $commandSet->first_child('configMapId')->text() =~ /^component/ ) {

            my $component_cmd = $commandSet->first_child('command');

            if ($component_cmd->first_child('arg')->text() =~ /\-c (\S+) /) {
                return $1;
            } else {
                croak ("can't have a component without an arg pointing to the component conf.");
            }
        }
        else {
            return 0;
        }
    }

    ## this private subroutine accepts the path to a conf file as an argument, parses
    ##  it, and checks for specific parameters we know we want to validate.  The values
    ##  for these parameters are stored within the attribute data structures of the
    ##  object for later validation.
    sub _parse_component_conf_ini {
        my ($self, $conf_file) = @_;
        
        if (! -e $conf_file ) {
            croak("file ($conf_file) not found");
        }
        
        my $conf = Ergatis::ConfigFile->new( -file => $conf_file );
        
        ## replace variables.
        $conf = $conf->import_includes();
        $conf->replace_variables();
        
        my $component_name;
        my $output_token;
        
        foreach my $section ( $conf->Sections() ) {
            
            ## have we got the component name yet?
            if (! defined $component_name) {
                ## if the section name has a space, the last word should be the component name
                if ($section =~ /\s+(\S+)$/) {
                    $component_name = $1;
                
                    ## all components must have an OUTPUT_TOKEN defined
                    $output_token = $conf->val( "output $component_name", '$;OUTPUT_TOKEN$;' );

                    ## remember this combo so we can check for duplicate name/token pairings
                    $self->{_tokens}->{$component_name}->{$output_token}++;
                }
            }
            
            foreach my $param ( $conf->Parameters($section) ) {

                ## skip those that don't have defined values;
                my $value = $conf->val($section, $param) || undef;
                next unless ( defined $value );
                
                ## here we add any parameters and their values into the necessary
                ##  data structures for later checking.
                if ( $param eq '$;OUTPUT_DIRECTORY$;' ) {
                    push @{  $self->{_output_dirs}{ $value }  }, "$component_name.$output_token";
                
                } elsif ( $param eq '$;INPUT_FILE$;' || $param eq '$;INPUT_LIST$;' ) {
                    push @{  $self->{_inputs}{ $value }  }, "$component_name.$output_token";
                
                } elsif ( $param eq '$;INPUT_DIRECTORY$;' ) {
                    ## these must be accompanied by an $;INPUT_EXTENSION$; definition
                    my $extension = $conf->val($section, '$;INPUT_EXTENSION$;');
                    push @{  $self->{_input_dirs}{ $value }{ $extension }  }, "$component_name.$output_token";
                }
            }
        }
    }
    
    ## this private subroutine accepts a directory and, optionally, a file extension and
    ##  returns the count of non-hidden files or folders within that directory.  If a file
    ##  extension is passed, the count is limited to those that match the extension.
    sub _directory_content_count {
        my ($self, $dir, $ext) = @_;
        
        opendir(my $dh, $dir) || return -1;
        
        my $found = 0;
        
        while (my $thing = readdir $dh) {
            ## hidden things don't count, anything else does
            next if ( $thing =~ /^\./ );
            
            ## are we matching an extension?
            if (defined $ext) {
                $found++ if ( $thing =~ /$ext$/ );
                
            } else {
                $found++
            }
        }
        
        return $found;
    }
}


1==1;
