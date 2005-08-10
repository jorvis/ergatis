package Workflow::IdGenerator;

=head1 NAME

IdGenerator.pm - A module for creating pipeline IDs.

=head1 SYNOPSIS

    use Workflow::IdGenerator;

    my $idgen = new Workflow::IdGenerator;

    my $pipelineid = $idgen->new_id();
    
=head1 DESCRIPTION

This is a module for creating unique pipeline IDs for Workflow.  It works
based on an incremented file method designed to handle simultaneous
ID requests from multiple processes, users and hosts.  If one process
tries to pull an ID and the file is already being used, a series of
tests allows it to wait nicely until the file is available.  We have avoided
the use of file-locking methods such as flock which do not function
on NFS.

This module requires a bit of initial setup before it can be used.  Simply
create a directory someplace, and put an empty file in it named like:

    current.1000000.id
    
The integer portion of this can be whatever you like - it will be the first
ID returned to any calling instances of the class.  This file and directory
should both be writeable by those who can use the class.

=head1 METHODS

=over 3

=item I<PACKAGE>->new( /%options )

Returns a newly created "Workflow::IdGenerator" object.  No arguments are necessary.
No individual script should create more than one instance of this class or there
will be a race condition on the log files.

=over 5

=item B<options>

=item I<id_dir>

This specifies the directory where the source id file is.  This is file is
read and incremented by this module to generate the unique IDs.

=item I<logging>

Boolean switch ( 1 or 0 ) to turn off/on logging.  On by default.

=item I<log_dir>

When logging is turned on, this designates the root where log files will
be written.  Within it, a directory will be created for each host accessing
this module.  Within that directory individual log files will be created
for each process.

=back

=item I<$OBJ>->next_id( )

This returns the next available pipeline ID and increments the counter for
subsequent requests.

=item I<$OBJ>->current_id( )

This returns the id most recently grabbed by this object.  Will return undef if
the next_id method hasn't yet been called on the object.

=back

=head1 TO DO / IMPROVEMENTS

A lot needs to be done in the way of error handling if a process dies
prematurely.  This could result in a lock file existing when the process
responsible for it no longer exists.

=head1 AUTHOR

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;
use warnings;
use Carp;
use File::Find;
use Sys::Hostname;

umask(0000);

## class data and methods
{

    my %_attributes = (
                        current_id  => undef,
                        id_dir      => '/usr/local/scratch/jorvis/id_generation',
                        logging  => 1,
                        log_dir     => '/usr/local/scratch/jorvis/id_generation/logs',
                        _id_file    => undef,
                        _logfh      => undef,
                      );

    ## class variables
    my $host = hostname();

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
            print $fh "starting IdGenerator log for process $$\n";
            $self->{_logfh} = $fh;
        }
        
        return $self;
    }
    
    sub current_id {
        ## just return the current_id
        $_[0]->{current_id};
    }
    
    sub id_dir {
        ## just return the directory containing the id and lock files
        $_[0]->{id_dir};
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
        my $self = shift;
    
        ## we want to keep trying until we're able to get an ID or have
        ##  a handled failure
        while (1) {
            ## is there a single id file?
            if ( $self->_get_id_files != 1) {
                $self->_log("debug: no individual id file.  restarting\n");
                next;
            }

            ## are there any locks?
            ## if so, wait 1
            if ( $self->_locks_exist ) {
                $self->_log("debug: existing lock file found.  restarting.\n");
                next;

            ## if not, try to lock
            } else {
                $self->_log("debug: attempting to create " . $self->{_id_file} . ".$$.lock\n");

                ## here we can have this NFS problem:
                ##      Client sends "rename fileA to fileB"
                ##      Server responds "ok, did that"
                ##      (but this response is lost in the network)
                ##      Client times out
                ##      Client re-sends "rename fileA to fileB"
                ##      Server responds "ah!, file 'fileA' not found"
                ## so if the rename reports as failed, we need to make sure it really did before we try again
                if (! rename($self->{_id_file}, $self->{_id_file} . ".$host.$$.lock") ) {
                    ## report failure to lock
                    $self->_log("error: failed to lock, verifying\n");

                    ## verify that the lock file doesn't exist
                    if (! -e $self->{_id_file} . ".$host.$$.lock") {
                        $self->_log("error: verified.  restarting.\n");
                        next;
                    }

                    ## else let it fall through, since it really did work.  NFS reporting problem
                    $self->_log("debug: not verified.  maybe NFS report was lost.  plowing on.\n");
                }

                $self->_log("debug: rename " . $self->{_id_file} . ", " . $self->{_id_file} . ".$host.$$.lock successful\n");

                ## if successful, 
                ## pull id, 
                if ($self->{_id_file} =~ /current\.(\d+)\.id/) {
                    $self->{current_id} = $1;
                    $self->_log("debug: got id " . $self->{current_id} . "\n");
                } else {
                    print "error: failed to get id from file name, dying.\n";
                    exit(1);
                }

                ## create next current_id file,
                $self->_log("debug: creating new id file " . $self->id_dir() . '/current.' . ($self->{current_id} + 1) . ".id\n");
                open(my $fh, ">" . $self->id_dir() . "/current." . ($self->{current_id} + 1) . ".id") || die "couldn't create next file id\n";
                close $fh;

                ## delete locked current_id file used.
                $self->_log("debug: deleting lock file " . $self->{_id_file} . ".$host.$$.lock\n");
                unlink($self->{_id_file} . ".$host.$$.lock") || die "couldn't unlink " . $self->{_id_file} . ".$host.$$.lock";

                ## finish this iteration
                $self->_log("debug: finished\n");
                last;
            }
        }
        
        ## return the id we got
        return $self->current_id();
        
    }  ## end next_id method block


    #####################
    ## private methods ##
    #####################
    sub _get_id_files {
        my $self = shift;
        my $ids_found = 0;

        find ( {    wanted =>   sub {
                                    if (/\.id$/) {
                                        $ids_found++;
                                        $self->{_id_file} = $File::Find::name;
                                        $self->_log("debug: current id_file seems to be " . $self->{_id_file} . "\n");
                                    }
                                }, 
                    no_chdir => 1
                }, $self->id_dir()
             );

        return $ids_found;
    }
    
    sub _locks_exist {
        ## checks to see if there are any locks found within the 
        ##  lock directory.  returns the number of locks found.
        my $self = shift;
        my $locks_found = 0;

        find ( {    wanted =>   sub {
                                    if (/\.lock$/) {
                                        $locks_found++;
                                    }
                                }, 
                    no_chdir => 1
               }, $self->id_dir()
             );

        return $locks_found;    
    }
    
    sub _log {
        my ($self, $msg) = @_;
        
        my $logfh = $self->{_logfh};
        
        ## print the message to the log file, if logging
        if ($self->logging) {
            print $logfh $msg;
        }
    }
}

1==1;
