package Ergatis::Logger;

# $Id$

# Copyright (c) 2002, The Institute for Genomic Research. All rights reserved.

=head1 NAME

Logger.pm - A module for providing debugging methods to Ergatis and Ergatis projects.

=head1 VERSION

This document refers to version 1.00 of Logger.pm, released MMMM, DD, YYYY.

=head1 SYNOPSIS

    use base Ergatis::Logger;

    I<Module: Running <method>.>

=head1 DESCRIPTION

=head2 Overview

This module is intended to provide Ergatis, and Ergatis client projects, with an easy to
use api for debugging. Methods for tracing execution and outputting debugging information
based upon a "debug" level will be provided.

=over 4

=cut


use strict;
use File::Basename;
use Data::Dumper;
use Log::Cabin;

my $_DEFAULT_LOG_CATEGORY = "Ergatis";
my $_IS_INIT=0;

=item new

B<Description:> The module constructor.

B<Parameters:> %arg, a hash containing attribute-value pairs to
initialize the object with. Initialization actually occurs in the
private _init method.

my $logger = new Ergatis::Logger('LOG_LEVEL'=>2, #verbose debugging
			       'LOG_FILE'=>my.log
			       );

B<Returns:> $self (A Ergatis::Logger object).

=cut

sub new {
    my ($class) = shift;
    my $self = bless {}, ref($class) || $class;
    $_IS_INIT=1;
    $self->{_caller} = join(',',caller());
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

    $self->{_LOG_LEVEL} = 0; #use this parameter to inc/dec log level from default
    $self->{_MODE} = "clobber"; #set to 'append' to append to file
    $self->{_LOG_FILE} = undef;

    my %arg = @_;
    foreach my $key (keys %arg) {
        $self->{"_$key"} = $arg{$key};
    }
    $self->_init_class_loggers();
}

=item $obj->_init_class_loggers($loglevel)

B<Description:> Set up Log4Perl default logger.  All other loggers will inheret from this logger.

B<Parameters:> $loglevel - integer to increase log level. > 1 will produce verbose debugging

B<Returns:> 

=cut

sub _init_class_loggers {
    my($self) = @_;
    
    # make sure that initialization is only done once
    if (Log::Cabin::initialized()) {
	return;
    }

    #set up defaults for root logger
    my $logsimple = new Log::Cabin();
 
    $logsimple->more_logging($self->{_LOG_LEVEL});

   
    if(defined $self->{_LOG_FILE}){
	$logsimple->set_file_output($self->{_LOG_FILE});
	chmod 0666,$self->{_LOG_FILE}; # set to world rw to prevent multi-user access issues
    }else{
	# Use default ouput to stderr
	$logsimple->set_output(*STDERR);
    }
    
    my $logger = $logsimple->get_logger($_DEFAULT_LOG_CATEGORY); 
    $logger->debug('--------------------------') if($logger->is_debug);
    $logger->debug("INIT_LOGGER $self->{_caller}") if($logger->is_debug);
    $logger->debug('--------------------------') if($logger->is_debug);
}


sub get_logger {
    my $name = shift;
    if(defined $name && !ref $name){
	return Log::Cabin::get_logger($_DEFAULT_LOG_CATEGORY."::".$name);
    }
    else{
	return Log::Cabin::get_logger($_DEFAULT_LOG_CATEGORY);
    }
}

sub get_default_logfilename{
    #return "/tmp/".basename($0).".$$.log";
    return '/dev/null';
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

