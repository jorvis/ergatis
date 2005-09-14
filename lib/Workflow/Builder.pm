package Workflow::Builder;

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
use Workflow::Logger;
use Workflow::Run;
use Workflow::Repository;
use Data::Dumper;

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
    my $self = bless {}, ref($class) || $class;
    $self->{_logger} = Workflow::Logger::get_logger(__PACKAGE__);
    $self->{_KEY_REGEX} = '$;\w+$;';
    $self->{_INICONFIGOBJ} = undef;
    $self->{_NODISTRIB} = 0;
    $self->_init(@_);
    if(!(exists $self->{_NAME})){
	$self->{_logger}->logdie("Name not specified via NAME argument");
    }
    $self->{_logger}->debug("Name specified as $self->{_NAME}") if($self->{_logger}->is_debug());
    
    return $self;
}


=item $obj->_init([%arg])

B<Description:> Initialize object variables

B<Parameters:> %arg, a hash containing attributes to initialize the testing
object with. Keys in %arg will create object attributes with the same name,
but with a prepended underscore.

B<Returns:> None.

=cut

sub _init {
    my $self = shift;
    my %arg = @_;
    foreach my $key (keys %arg) {
        $self->{"_$key"} = $arg{$key}
    }

}

sub set_config{
    my($self,$cfg) = @_;
    $self->{_logger}->debug("Setting configuration");
    $self->{_INICONFIGOBJ} = $cfg;
}

sub get_config{
    my $self = shift;
    $self->{_logger}->debug("Retrieving configuration");
    return $self->{_INICONFIGOBJ};
}

sub generate_instance_ini{
    my($self, $templateini, $instanceini) = @_;
    $self->{_logger}->debug("Generating instance ini file, $instanceini, from template=$templateini");
    my $errors = $self->_create_ini_file($self->{_INICONFIGOBJ}, $templateini, $instanceini);
    if($errors != 0){
	$self->{_logger}->error("Errors ($errors) encountered creating ini file $instanceini. Check logs");
	return ($errors*-1);
    }
    else{
	return 1;
    }
}

sub generate_instance_xml{
    my($self, $instanceini, $xmltemplatefile, $instancefile) = @_;
    $self->{_logger}->debug("Generating instance xml file $instancefile from xmltemplate=$xmltemplatefile and instanceini=$instanceini");
    my $createstatus = $self->_create_instance_file($instancefile,$instanceini,$xmltemplatefile);
    if($createstatus == 0){
	$self->{_logger}->debug("Created instance file $instancefile") if($self->{_logger}->is_debug());
	return $instancefile;
    }
    else{
	$self->{_logger}->logdie("Unable to create instance file $instancefile");
	return undef;
    }
}

sub _create_ini_file{
    my ($self, $iniconfobj, $templateinifile, $instanceinifile) = @_;
    
    $self->{_logger}->debug("Creating ini file $instanceinifile from ini template=$templateinifile");
    my $subs = {};
    my @sections = $iniconfobj->Sections();
    foreach my $section (@sections){
	my $name = $self->{_NAME};
	$self->{_logger}->debug("Checking section [$section] against $name") if($self->{_logger}->is_debug());
	if(($section =~ /$name/) || ($section =~ /init/)){
	    $self->{_logger}->debug("Entering section [$section]") if($self->{_logger}->is_debug());
	    my @parameters = $iniconfobj->Parameters($section);
	    foreach my $parameter (@parameters){
		my $replacevalue = $iniconfobj->val($section,$parameter);
		$self->{_logger}->debug("Storing $parameter as $replacevalue") if($self->{_logger}->is_debug());
		$subs->{$parameter} = $replacevalue;
	    }
	}
    }
    open( TEMPLATEINIFILE, "$templateinifile" ) or $self->{_logger}->logdie("Could not open template ini file $templateinifile");
    open( INSTANCEINIFILE, "+>$instanceinifile") or $self->{_logger}->logdie("Could not open new ini file $instanceinifile");
    
    while( my $line = <TEMPLATEINIFILE> ){
        ## don't replace vals on comment lines
        if ( $line !~ /^\;/ && $line !~ /^\#/ ) {
	    $line =~ s/(\$;[\w_]+\$;)/&replaceval($self,$1,$subs)/ge;
        }
	print INSTANCEINIFILE $line;
    }
    close TEMPLATEINIFILE;
    close INSTANCEINIFILE;
    return $self->_check_subs($instanceinifile);
}

sub _create_instance_file{
    my ($self, $instancefile, $inifile, $templatefile) = @_;
    my $logfile = "$instancefile.create.log";
    my $outfile = "$instancefile.create.out";
    my $wfexec = new Workflow::Run('nodistrib'=>$self->{_NODISTRIB});
    $self->{_logger}->debug("Creating instance file $instancefile from ini file $inifile and templatexml $templatefile");
    return $wfexec->CreateWorkflow($instancefile, $inifile, $templatefile, $logfile, $outfile);
}

sub _write_config{
    my ($self, $conffile, $confighash) = @_;
       
    open( CONFIG, ">$conffile" ) or $self->{_logger}->logdie("Could not open $conffile");

    foreach my $key (keys( %{$confighash} )){
	print CONFIG "$key=$confighash->{$key}\n";
    }
    
    close( CONFIG );
}

sub _check_subs{
    my($self,$file) = @_;
    $self->{_logger}->debug("Checking substitutions in $file") if($self->{_logger}->is_debug());
    open FILE, $file or $self->{_logger}->logdie("Could not open file $file");
    my $missed=0;
    while (my $line=<FILE>){
	while($line =~ /$self->{_KEY_REGEX}/){
	    my($key) = ($line =~ /($self->{_KEY_REGEX})/);
	    $self->{_logger}->warn("Key token $key not found in $file");
	    $line =~ s/($self->{_KEY_REGEX})/_MISSED_/;
	    $missed++;
	}
    }
    $self->{_logger}->debug("Returning $missed keys missed") if($self->{_logger}->is_debug());
    return $missed;
}

sub replaceval{
    my($self,$val,$keylookup) = @_;
    if(!(exists $keylookup->{$val})){
	$self->{_logger}->logdie("Bad key $val in ini file");
    }
    else{
	if($val =~ /TOGGLE\$\;$/){
	    my $val =  ($keylookup->{$val} =~ /\S/) ? $keylookup->{$val} : "";
	}
	else{
	    return $keylookup->{$val};
	}
    }
}

1;
