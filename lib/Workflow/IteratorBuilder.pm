package Workflow::IteratorBuilder;

# $Id$

# Copyright (c) 2002, The Institute for Genomic Research. All rights reserved.

=head1 NAME

Builder.pm - A module for building workflow instances

=head1 VERSION

This document refers to version $Name$ of frontend.cgi, $Revision$. 
Last modified on $Date$

=head1 SYNOPSIS



=head1 DESCRIPTION

=head2 Overview


=over 4

=cut


use strict;
use base qw(Workflow::Builder);
use Workflow::Logger;
use Data::Dumper;
use Config::IniFiles;
use File::Basename;
use Storable;

=item new

B<Description:> The module constructor.

B<Parameters:> %arg, a hash containing attribute-value pairs to
initialize the object with. Initialization actually occurs in the
private _init method.

my $builder = new Workflow::Builder('NAME'=>'blastp', #verbose debugging
				    );

B<Returns:> $self (A Workflow::Logger object).

=cut

sub new {
    my ($class) = shift;
    my $self = $class->SUPER::new(@_);
    $self->{_logger} = Workflow::Logger::get_logger(__PACKAGE__);
    return $self;
}


=item $obj->_init([%arg])

B<Description:> Tests the Perl syntax of script names passed to it. When
testing the syntax of the script, the correct directories are included in
in the search path by the use of Perl "-I" command line flag.

B<Parameters:> %arg, a hash containing attributes to initialize the testing
object with. Keys in %arg will create object attributes with the same name,
but with a prepended underscore.

B<Returns:> None.

=cut

sub _init {
    my $self = shift;
    #Each subflow mush have a unique name.  This name is stored in the configuration hash with the following key.
    $self->{_SUBFLOW_NAME_CONF_KEY} = '$;SUBFLOW_NAME$;'; 
    $self->{_CONFIGITERATOR_KEY} = '$;CONFIG_LIST$;'; 	 
    $self->{_INSTANCEITERATOR_KEY} = '$;INSTANCE_LIST$;'; 	 
    $self->{_TEMPLATE_KEY} = '$;TEMPLATE_FILE$;';
    $self->{_NODISTRIB} = 0;

    $self->{_ITERATOR_INI} = "default-iterator.ini";
    $self->{_ITERATOR_TEMPLATE} = "default-iterator_template.xml";
    $self->{_FILENAME_PREFIX} = "iterator";
    
    my %arg = @_;
    foreach my $key (keys %arg) {
        $self->{"_$key"} = $arg{$key}
    }
    if(!(exists $self->{_NAME})){
	$self->{_logger}->logdie("Name not specified via NAME argument");
    }
    if(!(exists $self->{_OUTPUT_DIR})){
	$self->{_logger}->logdie("ouput directory not specified via OUTPUT_DIR argument");
    }
}

sub set_iterator_template{
    my ($self, $iteratorini, $iteratortemplate) = @_;
    $self->{_ITERATOR_INI} = $iteratorini;
    $self->{_ITERATOR_TEMPLATE} = $iteratortemplate;
    return ($self->{_ITERATOR_INI},$self->{_ITERATOR_TEMPLATE});
}

sub get_iterator_template{
    my ($self) = @_;
    return ($self->{_ITERATOR_INI},$self->{_ITERATOR_TEMPLATE});
}

#The configuration hash ($confighash) must contain array(s) of values to be iterated.
#These values are used to create a subflow
#The configuration hash must also contain a unique name for each subflow.  
#This key is name ;_SUBFLOW_NAME; by default.
#Eg.
#$confighash = {';ASMBL;'=>('asmbl_1','asmbl_2','asmbl_3','asmbl_4'),
#               ';_SUBFLOW_NAME;'=>('iterator_1','iterator_2','iterator_3','iterator_4')
#               }

sub generate_iterator_instance{
    my ($self, $inifile, $xmltemplatefile, $confighash, $instancefile) = @_;
    $self->{_logger}->debug("Generating iterator instance $instancefile from inifile=$inifile xmltemplatefile=$xmltemplatefile")  if($self->{_logger}->is_debug());
    if(scalar(keys %$confighash) < 1){
	$self->{_logger}->logdie("No iterator parameter specified when creating iterator");
    }
    my @subflowconfigs;
    foreach my $elt (keys %$confighash){
	if(scalar(@{$confighash->{$elt}}) < 1){
	    $self->{_logger}->logdie(scalar(@{$confighash->{$elt}})." values for option $elt when trying to create iterator");
	}
	$self->{_logger}->debug("Creating subflow config hash for $elt") if($self->{_logger}->is_debug());
	for(my $i=0;$i<@{$confighash->{$elt}};$i++){
	    $self->{_logger}->debug("Storing $elt value $i as $confighash->{$elt}->[$i]") if($self->{_logger}->is_debug());
	    $subflowconfigs[$i]->{$elt} = $confighash->{$elt}->[$i];
	    if($elt eq $self->{_SUBFLOW_NAME_CONF_KEY}){
		if(exists $confighash->{$self->{_SUBFLOW_NAME_CONF_KEY}}){
		    $subflowconfigs[$i]->{$self->{_SUBFLOW_NAME_CONF_KEY}} = "$confighash->{$self->{_SUBFLOW_NAME_CONF_KEY}}->[$i]";
		}
		else{
		    $self->{_logger}->warn("Required key '$self->{_SUBFLOW_NAME_CONF_KEY}' is missing from iterator configuration hash");
		}
	    }
	}
    }
    my @subflowinis;
    my @subflowinstances;
    my $subflowconfigobj =  new Config::IniFiles(-import =>  $self->{_INICONFIGOBJ});
    foreach my $subflowconfighash (@subflowconfigs){
	$self->{_logger}->debug("Creating subflow with config ".Dumper($subflowconfighash)) if($self->{_logger}->is_debug());

	$self->{_logger}->debug("Subflow name ".$subflowconfighash->{$self->{_SUBFLOW_NAME_CONF_KEY}}) if($self->{_logger}->is_debug());
	my $subflowconffile = "$self->{_OUTPUT_DIR}/".$subflowconfighash->{$self->{_SUBFLOW_NAME_CONF_KEY}}.".conf";
	foreach my $key (keys %$subflowconfighash){
	    $self->{_logger}->debug("Storing $key=$subflowconfighash->{$key} in section [input $self->{_NAME}]") if($self->{_logger}->is_debug());
	    my $success = $subflowconfigobj->setval("input $self->{_NAME}",$key,$subflowconfighash->{$key});
	    if(!$success){
		$success = $subflowconfigobj->newval("input $self->{_NAME}",$key,$subflowconfighash->{$key});
	    }
	    if(!$success){
		$self->{_logger}->fatal("Not able to set value $key=$subflowconfighash->{$key} in section [input $self->{_NAME}]");
	    }
	}
	$self->{_logger}->debug("Writing config file $subflowconffile") if($self->{_logger}->is_debug());
	$subflowconfigobj->WriteConfig($subflowconffile) if($self->{_logger}->is_debug());

	my $subflowinifile = "$self->{_OUTPUT_DIR}/$subflowconfighash->{$self->{_SUBFLOW_NAME_CONF_KEY}}.ini";
	my $errors = $self->_create_ini_file($subflowconfigobj, $inifile, $subflowinifile);
	if($errors != 0){
	    $self->{_logger}->logdie("Errors encountered creating ini file $subflowinifile. Check logs");
	}
	my $subflowinstancefile = "$self->{_OUTPUT_DIR}/$subflowconfighash->{$self->{_SUBFLOW_NAME_CONF_KEY}}.xml";
	push @subflowinis,$subflowinifile;
	#note subflow instances are created by the workflow engine at runtime
	push @subflowinstances,$subflowinstancefile;
    }

    my $instancebase = basename($instancefile);

    my $iteratorini = "$self->{_OUTPUT_DIR}/$self->{_FILENAME_PREFIX}.$instancebase.ini";
    my $iteratorconffile = "$self->{_OUTPUT_DIR}/$self->{_FILENAME_PREFIX}.$instancebase.conf";
    
    $self->{_INICONFIGOBJ}->newval("workflowdocs $self->{_NAME}",$self->{_CONFIGITERATOR_KEY},join(', ',@subflowinis));
    $self->{_INICONFIGOBJ}->newval("workflowdocs $self->{_NAME}",$self->{_INSTANCEITERATOR_KEY},join(', ',@subflowinstances));
    $self->{_INICONFIGOBJ}->newval("workflowdocs $self->{_NAME}",$self->{_TEMPLATE_KEY},$xmltemplatefile);

    $self->{_INICONFIGOBJ}->WriteConfig($iteratorconffile);

    $self->_create_ini_file( $self->{_INICONFIGOBJ}, $self->{_ITERATOR_INI}, $iteratorini);
    my $createstatus = $self->_create_instance_file($instancefile,$iteratorini,$self->{_ITERATOR_TEMPLATE});
    if($createstatus == 0){
	$self->{_logger}->debug("Created instance file $instancefile") if($self->{_logger}->is_debug());
	return $instancefile;
    }
    else{
	$self->{_logger}->logdie("Unable to create instance file $instancefile");
	return undef;
    }
}

sub generate_iterator_arrayfile{
    my ($self,$iteratorconf) = @_;
    $self->_create_repository();
    my $iteratorarrayconfile = $self->_build_iteratorarrayconf_filename();
    store $iteratorconf, $iteratorarrayconfile;
    return $iteratorarrayconfile;
}

sub retrieve_iterator_arrayfile{
    my ($self) = @_;
    $self->_create_repository();
    return Storable::retrieve($self->_build_iteratorarrayconf_filename());
}

sub _build_iteratorarrayconf_filename(){
    my ($self) = @_;
    return $self->_build_instanceconf_filename()."array.storable";
}

1;
