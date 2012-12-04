package Ergatis::IdGenerator::IGSIdGenerator;
@ISA = qw(Ergatis::IdGenerator);

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

This is a module for creating unique pipeline/feature IDs for Ergatis.  It works
based on an incremented file method designed to handle simultaneous ID requests 
from multiple processes, users and hosts over NFS.  The module currently uses
the File::NFSLock module.

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

=head1 TO DO / IMPROVEMENTS

Some signal handling is done to attempt to ensure an ID file isn't left in a locked
or incomplete state if a process is killed.  The File::NFSLock module does occasionally
leave some temp and .NFSLock files behind when this happens, but they do not affect
other processes (in my testing).  More testing and handling could be done here.

There are other modules that claim to provide locking functionality, namely
DotLock and IPC::Locker.  While performance is definitely important, I'm primarily 
interested in getting something working without any race conditions.  I plan to
play with these other modules to see if they are equally stable and test
if they may give a performance advantage.

=head1 AUTHOR

    Kevin Galens
    kgalens@som.umaryland.edu

=cut

use strict;
use warnings;
use Carp;
use Sys::Hostname;
use IGS::MysqlUIDGenerator;

$|++;

umask(0000);

## class data and methods
{

    my %_attributes = (
                        id_repository           => undef,
                        logging                 => 0,
                        log_dir                 => undef,
                        _id_pending             => undef,
                        _id_repository_checked  => 0,
                        _fh                     => undef,
                        _logfh                  => undef,
                        _pool_sizes             => undef,  # will be hashref
                        _pools                  => undef,  # will be hashref
                      );

    ## class variables
    my $host = hostname();

    ## id service
    my $id_service;

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
        
        ## id_repository is required
        #   don't check it until the first next_id call
        if ($self->logging && ! defined $args{id_repository} ) {
            croak ("id_repository is a required argument for the constructor");
        }
        
        ## set the log dir (unless the user has)
        unless ( defined $self->{log_dir} ) {
            $self->{log_dir} = "$args{id_repository}/logs";
        }
        
        ## if we are logging, create a file
        if ($self->logging) {
            ## make sure the log directory exists (if we're logging)
            if(! -d $self->log_dir() ) {
                mkdir($self->log_dir(), 0777);
            }

            ## make sure the host subdir exists
            if (! -d $self->log_dir() . "/$host") {
                mkdir($self->log_dir() . "/$host", 0777);
            }

            open (my $fh, ">>" . $self->log_dir() . "/$host/$$.log") || croak ("can't create log file: $!");
            print $fh "info: starting IdGenerator log for process $$\n";
            $self->{_logfh} = $fh;
        }
       
        ## use EUIDService
        $id_service = new IGS::MysqlUIDGenerator;
        #unless ($id_service->ping()) {
        #    croak "Can't contact GUID service";
        #}
        
        return $self;
    }
    
    sub id_repository {
        ## just return the directory containing the id and lock files
        $_[0]->{id_repository};
    }

    sub log_dir {
        ## just return the directory containing the log files
        $_[0]->{log_dir};
    }    

    sub logging {
        ## just return the logging boolean
        $_[0]->{logging};
    }
    
    sub next_id {
        my ($self, %args) = @_;
        my $current_num = undef;

        ## check some required arguments
        $args{type} || croak "type is a required argument for the next_id method";
        $args{count} = 1 unless ( defined $args{count} );
        
        unless ( $args{type} eq 'pipeline' ) {
            $args{project} || croak "project is a required argument for the next_id method";
            $args{version} = 1 unless ( defined $args{version} );
        }
        
        ## has the id_repository been checked? (jay requested this not be in the constructor)
        if ( $self->logging && ! $self->{_id_repository_checked} ) {
            $self->_check_id_repository();
        }
        
        ## get ID
        
        ## from pool
        $current_num = $id_service->get_next_id();
        
        ## format and return the id we got
        ## if type is pipeline we just return the numerical portion.
        if ( $args{type} eq 'pipeline' ) {
            $self->_log("info: got pipeline ID $current_num");
            return $current_num;
        } else {
            $self->_log("info: got ID $args{project}.$args{type}.$current_num.$args{version}");
            return "$args{project}.$args{type}.$current_num.$args{version}";
        }
        
    }  ## end next_id method block


    ## get next raw GUID
    sub next_guid {
        my ($self, %args) = @_;
        my $current_num = undef;

        $current_num = $id_service->get_next_id();

        $self->_log("info: got GUID $current_num");
        
        return $current_num;
    
    }
    
    sub set_pool_size {
        my ($self, %args) = @_;
        
        for my $type ( keys %args ) {
            ## the value must be a number
            if ( $args{$type} !~ /^\d+$/ ) {
                croak("fatal: failed attempt to set pool size for type $type.  value '$args{$type}' must be an integer");
            }
            
            ## it cannot be 0
            if ( $args{$type} == 0 ) {
                croak("fatal: failed attempt to set pool size for type $type.  value cannot be 0");
            }
            
            $self->_log("info: setting pool size for type $type to $args{$type}");
            $self->{_pool_sizes}->{$type} = $args{$type};
        }

        ## get total id block size and set new size if bigger than current
        my $new_size = 0;
        foreach my $val(values(%{$self->{_pool_sizes}})) {
            $new_size += $val;
        }
        if ($new_size > $id_service->get_batch_size) {
            $id_service->set_batch_size($new_size);
        }
    }

    #####################
    ## private methods ##
    #####################

    sub _check_id_repository {
        my ($self, %args) = @_;
    
        ## make sure it is a valid ID repository
        unless ( -f "$self->{id_repository}/valid_id_repository" ) {
            croak ("fatal: Can't find file $self->{id_repository}/valid_id_repository.  the id repository $self->{id_repository} doesn't appear to be valid.  see the documentation for this module");
        }
        
        $self->{_id_repository_checked} = 1;
    }
    
    sub _log {
        my ($self, $msg) = @_;
        
        my $logfh = $self->{_logfh};
        
        ## print the message to the log file, if logging
        if ($self->logging) {
            print $logfh "$msg\n";
        }
    }
}

1==1;
