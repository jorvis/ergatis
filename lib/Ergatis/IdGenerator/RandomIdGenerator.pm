package Ergatis::IdGenerator::RandomIdGenerator;
@ISA = qw(Ergatis::IdGenerator);

=head1 NAME

RandomIdGenerator.pm - A module for creating unique feature or pipeline IDs.

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

This is a module for creating unique pipeline/feature IDs for Ergatis. It works 
via two methods depending on the type of ID requested.

Pipeline ID - If a pipeline ID is requested a unique ID is generted using the 
              Data::UUID module. This ID will be totally random as prevent users
              from guessing pipeline ID's that they may not have access too.

Feature ID - Works based on an incremented file method designed to handle simultaneous 
             ID requests from multiple processes, users and hosts over NFS.  The module 
             currently uses the File::NFSLock module.

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

    Joshua Orvis
    jorvis@tigr.org

    Cesar Arze
    carze@som.umaryland.edu

=cut

use strict;
use warnings;
use Carp;
use Sys::Hostname;
use POSIX;
use Data::UUID;

$|++;

umask(0000);

## class data and methods
{

    my %_attributes = (
                        id_repository           => undef,
                        logging                 => 1,
                        log_dir                 => undef,
                        _id_pending             => undef,
                        _id_repository_checked  => 0,
                        _fh                     => undef,
                        _logfh                  => undef,
                        _pool_sizes             => undef,  # will be hashref
                        _pools                  => undef,  # will be hashref
                        _uuid_generator         => undef, 
                      );

    ## class variables
    my $host = hostname();

    ## workaround Perl's F_SETLKW bug
    my $f_setlk = 6;
    my $f_setlkw = 7;

    sub new {
        my ($class, %args) = @_;

        ## make sure that the OS is linux (no support for other OS's)
        croak "IdGenerator currently only supports Linux"
            if $^O !~ /linux/;
        
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
        
        ## instantiate a Data::UUID generator for use in generating pipeline ID's
        $self->{_uuid_generator} = new Data::UUID;

        ## id_repository is required
        #   don't check it until the first next_id call
        if ( ! defined $args{id_repository} ) {
            croak ("id_repository is a required argument for the constructor");
        }
        
        ## set the log dir (unless the user has)
        unless ( defined $self->{log_dir} ) {
            $self->{log_dir} = "$args{id_repository}/logs";
        }
        
        ## if we are logging, create a file
        if ($self->logging) {
            ## make sure the log directory exists (if we're logging)
            my $loop = 1;
            while( $loop ) {
                select(undef, undef, undef, 0.3) if( $loop > 1 );
                if(! -d $self->log_dir() ) {
                    croak ("Can't create log_dir (".$self->log_dir()."): $!") if( $loop >= 5 );
                    $loop++;
                    mkdir($self->log_dir(), 0777) || next;
                }
                $loop = 0;
            }

            ## make sure the host subdir exists
            $loop = 1;
            while( $loop ) {
                select(undef, undef, undef, 0.3) if( $loop > 1 );
                if (! -d $self->log_dir() . "/$host") {
                    if( $loop >= 5 ) {
                        croak("Can't create log subdir (".$self->log_dir()."/$host )");
                    }
                    $loop++;
                    mkdir($self->log_dir() . "/$host", 0777) || next;
                }
                $loop = 0;
            }

            open (my $fh, ">>" . $self->log_dir() . "/$host/$$.log") || croak ("can't create log file: $!");
            print $fh "info: starting IdGenerator log for process $$\n";
            $self->{_logfh} = $fh;

        }
        
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
        my $idfh;  ## filehandle for the ID file

        ## check some required arguments
        $args{type} || croak "type is a required argument for the next_id method";
        $args{count} = 1 unless ( defined $args{count} );
        
        unless ( $args{type} eq 'pipeline' ) {
            $args{project} || croak "project is a required argument for the next_id method";
            $args{version} = 1 unless ( defined $args{version} );
        }
        
        ## has the id_repository been checked? (jay requested this not be in the constructor)
        unless ( $self->{_id_repository_checked} ) {
            $self->_check_id_repository();
        }
        
        ## get ID
        
        ## from pool
        if ( defined $self->{_pools}->{ $args{type} } && scalar @{ $self->{_pools}->{ $args{type} } } ) {
            $current_num = shift @{ $self->{_pools}->{$args{type}} };
        
        ## from file system (refresh pool if necessary)
        } else {
            ## if we are generating an ID for a pipeline we want to generate 
            ## our ID making use of Data::UUID
            if ($args{type} eq 'pipeline') {
                $current_num = $self->_generate_UUID();

                ## refresh pool if needed
                if ( defined $self->{_pool_sizes}->{ $args{type} } ) {
                    my $pool_size = $self->{_pool_sizes}->{ $args{type} };
                    
                    if ( $pool_size > 1 ) {
                        for ( 1 .. ($pool_size - 1) ) {
                            push @{ $self->{_pools}->{$args{type}} },
                                $self->_generate_UUID();
                        }
                    }
                } 
            } else {
                ## we need to trap interrupts here so we don't leave the ID file in a *bad* state
                #   i'm doing it locally so in order to minimize messing with any signal trapping the
                #   the calling script may have.
                # sub mymethod { local $SIG{INT} = sub { int-like-stuff }; ... rest of method ... }
                local $SIG{INT} = sub {
                    ## if a lock is open and we have a truncated ID file, we need to write the
                    #   pending ID.
                    if ( defined $self->{_id_pending} ) {
                            sysseek($self->{_fh}, 0, 0);
                            syswrite($self->{_fh}, $self->{_id_pending});
                            close $self->{_fh};
                    }
                    
                    ## if a lock is open, close it
                    $self->_unlock();
                    
                    croak("caught SIGINT. idgen process interrupted");
                };
            
                my $id_file = "$self->{id_repository}/next.$args{type}.id";
                sysopen($self->{_fh}, $id_file, O_RDWR|O_CREAT) ||
                    croak("failed to initialize file $id_file");
                $self->_lock();

                sysseek($self->{_fh}, 0, 0);
                sysread($self->{_fh}, $current_num, 100);
                $current_num = 1 if !$current_num;
                chomp $current_num;

                my $id_to_set;

                    ## do we have a pool to refresh?
                if ( defined $self->{_pool_sizes}->{ $args{type} } ) {
                    my $pool_size = $self->{_pool_sizes}->{ $args{type} };
                    $id_to_set = $current_num + $pool_size;

                            ## reload the pool
                    if ( $pool_size > 1 ) {
                        for ( 1 .. ($pool_size - 1) ) {
                            push @{ $self->{_pools}->{$args{type}} },
                             $_ + $current_num;
                        }
                    }
                } else {
                    $id_to_set = $current_num + 1;
                }
                        
                ## set next id, close
                $self->{_id_pending} = $id_to_set;
                sysseek($self->{_fh}, 0, 0);
                syswrite($self->{_fh}, $id_to_set);

                $self->{_id_pending} = undef;
                $self->_unlock();
                close $self->{_fh};
                undef $self->{_fh};
            }
        }

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
    
    sub _generate_UUID {
        my $self = shift;

        my $uuid = $self->{_uuid_generator}->create_str();
        $uuid = $self->_format_uuid($uuid);

        ## Check if our ID is unique and continue to generate new ID's 
        ## until we reach a unique ID.
        while ( ! $self->_check_uuid_uniqueness($uuid) ) {
            print STDERR "UUID $uuid is not unique; generating new UUID\n";
            $uuid = $self->{_uuid_generator}->create_str();
            $uuid = $self->_format_uuid($uuid);
        }

        return $uuid;
    }

    sub _format_uuid {
        my ($self, $uuid) = @_;
        $uuid =~ s/-//g;
        $uuid = substr($uuid, 0, 12);
        return $uuid;
    }

    sub _check_uuid_uniqueness {
        my ($self, $uuid) = @_;
        my $is_unique = 1;

        ## we want to look up our global ID repository to see if this pipeline
        ## exists prior to generating an ID.
        my $global_id_repo = $self->{'id_repository'};
        my $id_lookup = "$global_id_repo/global_id_lookup";

        if (! -d $id_lookup) {
            mkdir($id_lookup, 0777) || die "Could not create global ID lookup folder: $!";
        }

        my $id_file_subdirectory = "$id_lookup/" . substr($uuid, 0, 1) . "/" . substr($uuid, 1, 1);
        if (! -d $id_file_subdirectory) {
            #mkpath([$id_file_subdirectory, 1, 0777]) || die "Could not create uuid $uuid subdirectory $id_file_subdirectory: $!";
            $self->_run_system_cmd("mkdir -p $id_file_subdirectory");
        }

        ## Check if our pipeline ID exists in our global lookup, if it does not
        ## create a new pipeline file as we will be using this ID.
        if (-e "$id_file_subdirectory/$uuid.id") {
            $is_unique = 0;
        } else {
            ## TODO: find a cleaner way of touching this file.
            $self->_run_system_cmd("touch $id_file_subdirectory/$uuid.id");
        }

        return $is_unique;
    }

    sub _run_system_cmd {
        my ($self, $cmd) = @_;
        my @cmd_output;

        eval {
            @cmd_output = qx{$cmd 2>&1};
            if ( ($? << 8) != 0 ) { 
                die "@cmd_output";
            }   
        };  
        if ($@) {
            die "Error executing command $cmd: $@";
        }   
    }

    sub _log {
        my ($self, $msg) = @_;
        
        my $logfh = $self->{_logfh};
        
        ## print the message to the log file, if logging
        if ($self->logging) {
            print $logfh "$msg\n";
        }
    }

    sub _lock
    {
        my $self = shift;
        my $fl = pack("s! s! l! l! i", F_WRLCK, SEEK_SET, 0, 0, $$);
        fcntl($self->{_fh}, $f_setlkw, $fl) ||
            die "Error locking file: $!";
    }

    sub _unlock
    {
        my $self = shift;
        my $fl = pack("s! s! l! l! i", F_UNLCK, SEEK_SET, 0, 0, $$);
        fcntl($self->{_fh}, $f_setlk, $fl) ||
            die "Error unlocking file: $!";
    }
}

1==1;
