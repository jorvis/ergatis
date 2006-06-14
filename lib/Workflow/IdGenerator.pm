package Workflow::IdGenerator;

=head1 NAME

IdGenerator.pm - A module for creating unique feature or pipeline IDs.

=head1 SYNOPSIS

    use Workflow::IdGenerator;

    my $idgen = Workflow::IdGenerator->new( id_repository => '/path/to/some/id_repository' );

    ## to get a new pipeline id
    my $pipelineid = $idgen->next_id( type => 'pipeline' );
    
    ## get get a single id
    ## we have to specify a project here unless it was specified when calling new
    my $transcript_id = $idgen->next_id( type => 'transcript', project => 'aa1' );
    
    ## to reserve a range of ids (returns an array ref)
    my $exon_ids = $idgen->next_id( type => 'exon', project => 'rca1', count => 5 );
    
=head1 DESCRIPTION

This is a module for creating unique pipeline/feature IDs for Workflow.  It works
based on an incremented file method designed to handle simultaneous
ID requests from multiple processes, users and hosts.  If one process
tries to pull an ID and the file is already being used, a series of
tests allows it to wait nicely until the file is available.  We have avoided
the use of file-locking methods such as flock which do not function
on NFS.

This file and directory should both be writeable by those who can use the class.

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
module.  This prevents random directories from accidentally being used as
ID repositories.  To set up a directory, just do this:

    mkdir /path/to/some/id_repository
    touch /path/to/some/id_repository/valid_id_repository

The second command will create an empty file that just serves as a marker
to indicate that the directory is meant for ID generation.  The module will
take care of any other necessary file and directory structure needs.

=head1 TO DO / IMPROVEMENTS

A lot needs to be done in the way of error handling if a process dies
prematurely.  This could result in a lock file existing when the process
responsible for it no longer exists.

Sam found another possible method of doing this which is close to what we are
doing here.  I'll paste it below in case we want to experiment using link for
this later.

O_EXCL When used with O_CREAT, if the file already exists it is an error and 
the open will fail. In this context, a symbolic link exists, regardless of 
where its points to.  O_EXCL is broken on NFS file  systems,  programs which  
rely on it for performing locking tasks will contain a race condition.  The
solution for performing atomic file locking using a lockfile is to create a 
unique file on the same fs (e.g., incorporating hostname and pid), use link(2) 
to make a link to the lockfile. If link() returns 0, the lock is successful.  
Otherwise, use stat(2) on the unique file to check if its link  count  has
increased to 2, in which case the lock is also successful.

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
                        id_repository => undef,
                        logging     => 1,
                        log_dir     => undef,
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
        
        ## id_repository is required
        if ( $args{id_repository} ) {
            ## make sure it is a valid ID repository
            unless ( -f "$args{id_repository}/valid_id_repository" ) {
                croak ("the id repository $args{id_repository} doesn't appear to be valid.  see the documentation for this module");
            }
        } else {
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
            print $fh "starting IdGenerator log for process $$\n";
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

        ## check some required arguments
        $args{type} || croak "type is a required argument for the next_id method";
        $args{count} = 1 unless ( defined $args{count} );
        
        unless ( $args{type} eq 'pipeline' ) {
            $args{project} || croak "project is a required argument for the next_id method";
            $args{version} = 1 unless ( defined $args{version} );
        }
    
        ## we want to keep trying until we're able to get an ID or have
        ##  a handled failure
        while (1) {
        
            ## is there a single id file?
            my ($any, $unlocked) = $self->_get_counts_of_type($args{type});
            if ( $any ) {
                if ( $unlocked != 1 ) {
                    $self->_log("debug: no individual id file.  restarting\n");
                    sleep(2);
                    next;                
                }
            } else {
                ## none were found of this type.  initialize
                $self->_log("debug: no arch for type $args{type} found.  initializing.\n");
                open(my $newfh, ">$self->{id_repository}/$args{type}.1.id") || croak("failed to initialize id arch file for type $args{type}: $!");
                $self->{_id_file} = "$self->{id_repository}/$args{type}.1.id";
            }

            ## are there any locks?
            ## if so, wait 1
            if ( $self->_locks_exist ) {
                $self->_log("debug: existing lock file found.  restarting.\n");
                sleep(2);
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
                        sleep(2);
                        next;
                    }

                    ## else let it fall through, since it really did work.  NFS reporting problem
                    $self->_log("debug: not verified.  maybe NFS report was lost.  plowing on.\n");
                }

                $self->_log("debug: rename " . $self->{_id_file} . ", " . $self->{_id_file} . ".$host.$$.lock successful\n");

                ## if successful, 
                ## pull id, 
                if ($self->{_id_file} =~ /$args{type}\.(\d+)\.id/) {
                    $current_num = $1;
                    $self->_log("debug: got id $current_num\n");
                } else {
                    $self->_log("error: failed to get id from file name, dying.\n");
                    exit(1);
                }

                ## create next current_id file,
                $self->_log("debug: creating new id file " . $self->id_repository() . "/$args{type}." . ($current_num + $args{count}) . ".id\n");
                open(my $fh, ">" . $self->id_repository() . "/$args{type}." . ($current_num + $args{count}) . ".id") || die "couldn't create next file id\n";
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
        ## if type is pipeline we just return the numerical portion.
        my $ids = [];

        if ($args{count} > 1 ) {
            for ( my $i=0; $i < $args{count}; $i++ ) {
                if ( $args{type} eq 'pipeline' ) {
                    push @$ids, $current_num + $i;
                } else {
                    push @$ids, "$args{project}.$args{type}." . ($current_num + $i) . ".$args{version}";
                }
            }
            return $ids;

        } else {
            if ( $args{type} eq 'pipeline' ) {
                return $current_num;
            } else {
                return "$args{project}.$args{type}.$current_num.$args{version}";
            }
        }
        
    }  ## end next_id method block


    #####################
    ## private methods ##
    #####################
    sub _get_counts_of_type {
        my ($self, $type) = @_;
        my $found_unlocked = 0;
        my $found_any = 0;

        find ( {    wanted =>   sub {
                                    if (/$type\.\d+\.id$/) {
                                        $found_unlocked++;
                                        $found_any++;
                                        $self->{_id_file} = $File::Find::name;
                                        $self->_log("debug: current id_file seems to be " . $self->{_id_file} . "\n");
                                    } elsif (/$type\.+/) {
                                        $found_any++;
                                    }
                                }, 
                    no_chdir => 1
                }, $self->id_repository()
             );

        return ($found_any, $found_unlocked);
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
               }, $self->id_repository()
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
