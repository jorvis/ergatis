#!/usr/local/bin/perl

# $Id$

# Copyright (c) 2002, The Institute for Genomic Research. All rights reserved.

=head1 NAME

frontend.cgi - One line summary of purpose of class (or file).

=head1 VERSION

This document refers to version $Name$ of frontend.cgi, $Revision$. Last modified on $Date$.

=head1 SYNOPSIS

Short examples of code that illustrate the use of the class (if this file is a class).

=head1 DESCRIPTION

=head2 Overview

An overview of the purpose of the file.

=head2 Constructor and initialization.

if applicable, otherwise delete this and parent head2 line.

=head2 Class and object methods

if applicable, otherwise delete this and parent head2 line.

=cut


use strict;
use CGI qw(:standard);
use CGI::Carp qw(fatalsToBrowser set_message);

use Sybil;

my $logger;
BEGIN {
    $logger = new Coati::Logger('LOG_FILE'=>"/tmp/frontend.cgi.log",
				'LOG_LEVEL'=>&param('DEBUG')); 
}

my $user = param("user");
my $password = param("password");
my $db = param("db");
my $testing = param("testing");

my $sybil = new Sybil( user     => $user,
			       password => $password,
			       db       => $db
			       );

# These are Coati methods.
my $retcoati = $sybil->testCoati();
$logger->get_logger()->info("testCoati called and returned value $retcoati");

# These are Sybil specific methods.
my $retproject = $sybil->testSybil();
$logger->get_logger()->info("testCoati called and returned value $retproject");

# This one should fail and therefore AUTOLOAD will step in.
$sybil->shouldfail;


__END__

=head1 ENVIRONMENT

List of environment variables and other O/S related information
on which this file relies.

=head1 DIAGNOSTICS

CHADODUMMY:user=access&password=access&db=pfa1test

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

