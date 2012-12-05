package Ergatis::IdGenerator;

=head1 NAME

IdGenerator.pm - A module for creating unique feature or pipeline IDs.

=head1 SYNOPSIS

    use Ergatis::IdGenerator;

    my $idgen = Ergatis::IdGenerator->new( id_repository => '/path/to/some/id_repository' );

    ## to get a new pipeline id
    my $pipelineid = $idgen->next_id( type => 'pipeline' );
    
    ## optional. set the size of the ID pool for any feature types.  (more efficient)
    $idgen->set_pool_size( exon => 40, transcript => 20 );
    
    ## get get a single id
    ## we have to specify a project here unless it was specified when calling new
    my $transcript_id = $idgen->next_id( type => 'transcript', project => 'aa1' );
    
=head1 DESCRIPTION

This is a module for creating unique pipeline/feature IDs for Ergatis. 

This is an abstract class that instantiates an object of an IdGenerator subclass, defined in Ergatis::IdGenerator::Config.

Subclasses should implement all of the abstract methods present here.

=head1 METHODS

=over 3

=item I<PACKAGE>->new( /%options )

creates the IdGenerator object.  

=over 5

=item B<options>

=item I<id_repository>

required.  this is the path to a directory made to serve as an ID repository for 
unique ID generation.  see the ID REPOSITORY section below for details.

=item I<logging>

Optional. Boolean switch ( 1 or 0 ) to turn on/off logging (default = 1).

=item I<log_dir>

Optional. When logging is turned on, this designates the root where log files will
be written.  Within it, a directory will be created for each host accessing
this module.  Within that directory individual log files will be created
for each process.  By default, this will be created as 'logs' directory under
the id_repository.

=back

=item I<$OBJ>->next_id( )

This returns the next available ID of the passed type and increments the counter 
for subsequent requests.

=over 5

=item B<options>

=item I<type>

Required.  This corresponds to the CV term (SO, SOFA, other) representing the feature
type you want to create.

=item I<project>

Required for non-'pipeline' types.  This is usually a short abbreviation for the
project/database the feature belongs to such as 'aa1' for Aedes aegypti.

=item I<count>

Optional.  Use this is you want to reserve multiple IDs at once.  IDs will be
returned as an array reference.  (default = 1)

=item I<version>

Optional.  The IDs returned by this module contain version information.  Use this
option to overwrite the default version (1).

=back

=item I<$OBJ>->set_pool_size( some_type => N );

This method allows for more efficient ID retrieval by telling the module to pull
blocks of IDs and then use that pool when calling the next_id() method instead of
hitting the file system for each request.  This helps to decrease instances of
heavily locked files and should always be used if possible.  This method can be used
to set (or overwrite) the pool size for as many feature types as necessary by passing
the types and values as a hash (key = type, value = pool size).

Unused IDs in a pool are NOT returned to other processes, so you should attempt to
set a reasonable value here instead of just pooling an arbitrary, overly-large number.

=back

=head1 NAMING CONVENTION

This module will return either pipeline IDs or feature IDs, depending on the
value used in the --type option.  For pipeline IDs, a simple incremented integer
is returned.  Feature IDs contain more information, such as project, feature
type and versioning information.  The current format is like:

    aa1.transcript.14820.1
    
This feature ID indicates version 1 of transcript #14820 in the aa1 project.

=head1 ID REPOSITORY

It is easy to set up an ID repository, but it must be done prior to using the 
module.  This prevents random directories from being accidentally used as
ID repositories.  To set up a directory, just do this:

    mkdir /path/to/some/id_repository
    touch /path/to/some/id_repository/valid_id_repository

The second command will create an empty file that just serves as a marker
to indicate that the directory is meant for ID generation.  The module will
take care of any other necessary file and directory structure needs.

This directory will be checked the first time next_id() is called.

=head1 AUTHOR

    Brett Whitty
    bwhitty@jcvi.org

=cut

use strict;
use warnings;
use Carp;
use Ergatis::IdGenerator::Config;

$|++;

## class data and methods
{
    sub new {
        my ($class, %args) = @_;
        if(! defined($args{logging})) {
            $args{logging} = 0;
        }
        ## create the object
        my $self = new $Ergatis::IdGenerator::Config::class(%args);
        
        return $self;
    }
   
    ## returns the directory for id-generation related files
    sub id_repository {
        return undef;
    }

    ## returns the directory containing log files
    sub log_dir {
        return undef;
    }    

    ## boolean check for whether logging is enabled
    sub logging {
        return 0;
    }
    
    ## returns the next id available
    sub next_id {
       croak "Method next_id is not implemented"; 
    }

    ## sets a the size of the pool of ids to prefetch, if possible
    sub set_pool_size {
        croak "Method set_pool_size is not implemented";
    }

}

1;
