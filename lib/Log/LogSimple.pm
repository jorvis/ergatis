package Log::LogSimple;

# $Id$

# Copyright (c) 2002, The Institute for Genomic Research. All rights reserved.

=head1 NAME

Logger.pm - A module for providing debugging methods to Workflow and Workflow projects.

=head1 VERSION

This document refers to version 1.00 of Logger.pm, released MMMM, DD, YYYY.

=head1 SYNOPSIS

    use base Workflow::Logger;

    I<Module: Running <method>.>

=head1 DESCRIPTION

=head2 Overview

This module is intended to provide Workflow, and Workflow client projects, with an easy to
use api for debugging. Methods for tracing execution and outputting debugging information
based upon a "debug" level will be provided.

=over 4

=cut


use strict;
use Time::HiRes;
use Sys::Hostname;
use IO::Tee;
use Log::LogSimpleLogger;

my $_DEFAULT_LOG_CATEGORY = "LogSimple";
my $_IS_INIT=0;
my $OFF=0;
my $FATAL=1;
my $ERROR=2;
my $WARN=3;
my $INFO=4;
my $DEBUG=5;
my $ALL=6;

my $LEVELTEXT = {0=>'OFF',
		 1=>'FATAL',
		 2=>'ERROR',
		 3=>'WARN',
		 4=>'INFO',
		 5=>'DEBUG',
		 6=>'ALL'};

my $SINGLETON_INSTANCE=undef;
my $OUTPUT_HANDLE = undef;
my @ALLLOGGERS;

=item new

B<Description:> The module constructor.

B<Parameters:> %arg, a hash containing attribute-value pairs to
initialize the object with. Initialization actually occurs in the
private _init method.

my $logger = new Workflow::Logger('LOG_LEVEL'=>2, #verbose debugging
			       'LOG_FILE'=>my.log
			       );

B<Returns:> $self (A Workflow::Logger object).

=cut

sub new {
    my ($class) = shift;
    if($_IS_INIT == 1){
	return $SINGLETON_INSTANCE;
    }
    else{
	
	my $self = bless {}, ref($class) || $class;
	$self->_init(@_);

	if($self->{_SKIP_INIT}){
	}
	else{
	    $_IS_INIT=1;
	    $SINGLETON_INSTANCE = $self;
	    foreach my $loggerinst (@ALLLOGGERS){
		$loggerinst->set_logger_instance($SINGLETON_INSTANCE);
	    }
	}
	return $self;
    }
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

    $self->{_DEFAULT_LOG_LEVEL} = $WARN;
    $self->{_LOG_LEVEL} = $self->{_DEFAULT_LOG_LEVEL};
    $self->{_LOG_FILE} = undef;
    $self->{_HOSTNAME} = hostname;
    $self->{_PID} = $$;
    $self->{_CLOBBER} = 1;
    my %arg = @_;
    foreach my $key (keys %arg) {
        $self->{"_$key"} = $arg{$key};
    }

    if(!defined $self->{_LOG_FILE}){
	open($self->{_OUTPUT_HANDLE},">&=",STDERR) or die "Couldn't copy stderr";
    }
    else{
	$self->set_file_output($self->{_LOG_FILE});
    }
}

sub initialized{
    return $_IS_INIT;
}

sub get_instance{
    return $SINGLETON_INSTANCE;
}

sub set_file_output{
    my($self,$filename) = @_;
    $self->{_LOG_FILE} = $filename;
    my $filehandle;
    if($self->{_CLOBBER}){
	open($filehandle,"+>",$self->{_LOG_FILE})
	    or die "Can't open log file for writing $self->{_LOG_FILE}";
    }
    else{
	open($filehandle,">",$self->{_LOG_FILE})
	    or die "Can't open log file for writing $self->{_LOG_FILE}";
    }
    my $stderr;
    open($stderr,">&=",*STDERR) or die "Couldn't copy stderr";
    $self->{_OUTPUT_HANDLE} = $filehandle;
    *STDERR =  IO::Tee->new($stderr,$filehandle);
}

sub set_output{
    my($self,$handle) = @_;
    $self->{_OUTPUT_HANDLE} = $handle;
}

sub level{
    my($self,$level) = @_;
    if($level =~ /^\-*\d+$/){
	$self->{_LOG_LEVEL} = $level;
    }
    return $self->{_LOG_LEVEL};
}

sub more_logging{
    my($self,$level) = @_;
    if($level =~ /^\-*\d+$/){
	$self->{_LOG_LEVEL} += $level;
    }
}

sub less_logging{
    my($self,$level) = @_;
    if($level =~ /^\-*\d+$/){
	$self->{_LOG_LEVEL} -= $level;
    }
}

sub _output{
    my($self,$msg,$loggername,$level,$package,$filename,$line,$subroutine) = @_;
    my $datestamp = localtime(time());
    if(defined $self->{_OUTPUT_HANDLE}){
	print {$self->{_OUTPUT_HANDLE}} "$loggername $LEVELTEXT->{$level} $datestamp $self->{_HOSTNAME}:$self->{_PID} $filename:$package:$subroutine:$line || $msg\n";
    }
}

##################################################
sub get_logger {
##################################################
    my($class, @args) = @_;

    my $singleton;
    if(defined $SINGLETON_INSTANCE){
	$singleton = $SINGLETON_INSTANCE;
    }
    elsif(ref $class){
	$singleton = $class;
    }
    else{
	#die "get_logger called when Log::LogSimple not instantiated";
	$singleton = new Log::LogSimple('SKIP_INIT'=>1);
    }
    if(!ref $class){

	@args = ($class,@args);
    }
    my $loggerinst = new Log::LogSimpleLogger($singleton,@args);
    push @ALLLOGGERS,$loggerinst;
    return $loggerinst;
}

1;

__END__

=back

=head1 ENVIRONMENT

This module does not use or set any environment variables. The standard
module, File::Basename is required.

=head1 DIAGNOSTICS

=over 4

=item "Error message that may appear."

Explanation of error message.

=item "Another message that may appear."

Explanation of another error message.

=back

=head1 BUGS

Description of known bugs (and any workarounds). Usually also includes an
invitation to send the author bug reports.

=head1 SEE ALSO

List of any files or other Perl modules needed by the file or class and a
brief description why.

=head1 AUTHOR(S)

 The Institute for Genomic Research
 9712 Medical Center Drive
 Rockville, MD 20850

=head1 COPYRIGHT

Copyright (c) 2002, The Institute for Genomic Research. All Rights Reserved.

