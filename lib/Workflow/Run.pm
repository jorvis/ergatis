package Workflow::Run;

# $Id$

# Copyright (c) 2002, The Institute for Genomic Research. All rights reserved.

=head1 NAME

Run.pm - A module for building workflow instances


=head1 SYNOPSIS



=head1 DESCRIPTION

=head2 Overview


=over 4

=cut


use strict;
BEGIN {
use Workflow::Logger;
}
use Data::Dumper;

=item new

B<Description:> The module constructor.

B<Parameters:> %arg, a hash containing attribute-value pairs to
initialize the object with. Initialization actually occurs in the
private _init method.

my $builder = new Workflow::Builder('NAME'=>'blastp', #verbose debugging
				    'REPOSITORY_ROOT'=>'/usr/local/annotation',
				    'DATABASE'=>'tryp'
				    );

B<Returns:> $self (A Workflow::Builder object).

=cut

sub new {
    my ($class) = shift;
    my $self = bless {}, ref($class) || $class;
    $self->{_logger} = Workflow::Logger::get_logger(__PACKAGE__);
    $self->_init(@_);
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
    $self->{_WORKFLOW_EXEC_DIR} = "$ENV{'WORKFLOW_WRAPPERS_DIR'}" || ".";
    $self->{_WORKFLOW_CREATE_EXEC} = "CreateWorkflow.sh";
    $self->{_WORKFLOW_RUN_EXEC} = "RunWorkflow.sh";
    $self->{_nodistrib} = 0;

    my %arg = @_;
    foreach my $key (keys %arg) {
        $self->{"_$key"} = $arg{$key}
    }
}

sub CreateWorkflow{
    my($self,$instance, $ini, $template, $log, $outfile) = @_;
    my $execstr = "$self->{_WORKFLOW_EXEC_DIR}/$self->{_WORKFLOW_CREATE_EXEC} -t $template -c $ini -i $instance -l $log -o $outfile";
    $self->{_logger}->debug("Exec via system: $execstr") if ($self->{_logger}->is_debug());
    my $debugstr = "";
    if($self->{_logger}->is_debug()){
	$debugstr = "-v 6";
    }
    
    ## we want to reset the huge ENV so that it doesn't overload XML-RPC
    %ENV = (
                PATH            => $ENV{PATH},
                HOST            => $ENV{HOST},
                USER            => $ENV{USER},
                GROUP           => $ENV{GROUP},
                WORKFLOW_WRAPPERS_DIR => $ENV{WORKFLOW_WRAPPERS_DIR},
                
                ## next two needed by hmmpfam2htab
                HMM_SCRIPTS     => '/usr/local/devel/ANNOTATION/hmm/bin',
                SYBASE          => '/usr/local/packages/sybase',
                
                ## for augustus
                AUGUSTUS_CONFIG_PATH => '/usr/local/devel/ANNOTATION/jorvis/augustus/config',
                #AUGUSTUS_CONFIG_PATH => '/usr/local/annotation/TTA1/augustus/config',
                #AUGUSTUS_CONFIG_PATH => '/usr/local/devel/ANNOTATION/EGC_utilities/AUGUSTUS/augustus/config',
            
                ## for genewise
                WISECONFIGDIR => '/usr/local/devel/ANNOTATION/EGC_utilities/WISE2/wise2.2.0/wisecfg',

                ## for tRNAscan-SE
                ANNOT_DEVEL => $ENV{ANNOT_DEVEL},

                ## should workflow take care of this?
                SGE_ROOT        => '/local/n1ge',
                SGE_CELL => 'tigr',
                SGE_QMASTER_PORT => '536',
                SGE_EXECD_PORT => '537',
                SGE_ARCH => 'lx26-x86',
           );
    
    my $pid;
    if($pid = fork){
	$SIG{'TERM'} = sub {
	    my $signal = shift;
	    $SIG{'TERM'} = '';
	    kill $signal,$pid;
	};
	$SIG{'INT'} = sub {
	    my $signal = shift;
	    $SIG{'INT'} = '';
	    kill $signal,$pid;
	};
	$SIG{'QUIT'} = sub {
	    my $signal = shift;
	    $SIG{'QUIT'} = '';
	    kill $signal,$pid;
	};
	waitpid($pid,0);
	my $ret = $?;
	$ret >>= 8;

	if($self->{_nodistrib} == 1){
	    $self->_replacedistrib($instance,"$instance.nodistrib");
	    `cp $instance $instance.bak`;
	    `cp $instance.nodistrib $instance`;
	}

	return $ret;
    }
    else{
	die "cannot fork: $!" unless defined $pid;
	exec("$execstr $debugstr");
    }
}

sub RunWorkflow{
    my($self,$instance, $log, $outfile) = @_;
    my $execstr = "$self->{_WORKFLOW_EXEC_DIR}/$self->{_WORKFLOW_RUN_EXEC} -i $instance -l $log -o $outfile";
    $self->{_logger}->debug("Exec via system: $execstr") if ($self->{_logger}->is_debug());
    my $debugstr = "";
    if($self->{_logger}->is_debug()){
        $debugstr = "-v 6";
    }
    
    ## we want to reset the huge ENV so that it doesn't overload XML-RPC
    %ENV = (
                PATH            => $ENV{PATH},
                HOST            => $ENV{HOST},
                USER            => $ENV{USER},
                GROUP           => $ENV{GROUP},
                WORKFLOW_WRAPPERS_DIR => $ENV{WORKFLOW_WRAPPERS_DIR},
                
                ## next two needed by hmmpfam2htab
                HMM_SCRIPTS     => '/usr/local/devel/ANNOTATION/hmm/bin',
                SYBASE          => '/usr/local/packages/sybase',

                ## for augustus
                AUGUSTUS_CONFIG_PATH => '/usr/local/devel/ANNOTATION/jorvis/augustus/config',
                #AUGUSTUS_CONFIG_PATH => '/usr/local/annotation/TTA1/augustus/config',
                #AUGUSTUS_CONFIG_PATH => '/usr/local/devel/ANNOTATION/EGC_utilities/AUGUSTUS/augustus/config',

                # for genewise
                WISECONFIGDIR => '/usr/local/devel/ANNOTATION/EGC_utilities/WISE2/wise2.2.0/wisecfg',
                
                ## for tRNAscan-SE
                ANNOT_DEVEL => $ENV{ANNOT_DEVEL},
                
                ## should workflow take care of this?
                SGE_ROOT        => '/local/n1ge',
                SGE_CELL => 'tigr',
                SGE_QMASTER_PORT => '536',
                SGE_EXECD_PORT => '537',
                SGE_ARCH => 'lx26-x86',
           );

    my $pid;
    if($pid = fork){
	$SIG{'TERM'} = sub {
	    my $signal = shift;
	    $SIG{'TERM'} = '';
	    print STDERR "Caught signal $signal. Killing spawned workflow $pid\n";
	    kill $signal,$pid;
	};
	$SIG{'INT'} = sub {
	    my $signal = shift;
	    $SIG{'INT'} = '';
	    print STDERR "Caught signal $signal. Killing spawned workflow $pid\n;";
	    kill $signal,$pid;
	};
	$SIG{'QUIT'} = sub {
	    my $signal = shift;
	    $SIG{'QUIT'} = '';
	    print STDERR "Caught signal $signal. Killing spawned workflow $pid\n";
	    kill $signal,$pid;
	};
	waitpid($pid,0);
	my $ret = $?;
	$ret >>= 8;
	return $ret;
    }
    else{
	die "cannot fork: $!" unless defined $pid;
	exec("$execstr $debugstr");
    }
}

sub _replacedistrib{
    my($self,$file,$outputfile) = @_;
    open( FILEIN, "$file" ) or $self->{_logger}->logdie("Could not open file $file");
    open( FILEOUT, "+>$outputfile") or $self->{_logger}->logdie("Could not open output file $outputfile");
    
    while( my $line = <FILEIN> ){
	$line =~ s/Distributed/Unix/g;
	$line =~ s/commandSet(.*)type="parallel"/commandSet$1type="serial"/g;
	print FILEOUT $line;
    }
    close FILEIN;
    close FILEOUT;
    return $outputfile;
}

1;
