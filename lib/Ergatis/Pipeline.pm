package Ergatis::Pipeline;

=head1 NAME

Ergatis::Pipeline.pm - A class for representing pipelines

=head1 SYNOPSIS

    use Ergatis::Pipeline;

    my $conf = Ergatis::Pipeline->new( id => '1234', path => '/path/to/some/pipeline.xml' );

=head1 DESCRIPTION

=head1 METHODS

=over 3

=item I<PACKAGE>->new( path => I<'/path/to/some.xml'> )

Returns a newly created "Ergatis::Pipeline" object.

=item I<$OBJ>->path( )

Returns the path to a given pipeline.

=item I<$OBJ>->id( )

Returns the ID of a given pipeline.

=back

=head1 AUTHOR

    Joshua Orvis
    jorvis@users.sf.net

=cut

use strict;
use warnings;
use Carp;
use Sys::Hostname;
use IO::File;
use Ergatis::Utils qw(build_twig create_progress_bar update_progress_bar handle_component_status_changes);

umask(0000);

## class data and methods
{

    my %_attributes = (
        path       => undef,
        id         => undef,
        debug      => 1,
        debug_file => '/tmp/ergatis.run.debug',
    );

    ## class variables
    ## currently none

    sub new {
        my ( $class, %args ) = @_;

        ## create the object
        my $self = bless {%_attributes}, $class;

        ## set any attributes passed, checking to make sure they
        ##  were all valid
        for ( keys %args ) {
            if ( exists $_attributes{$_} ) {
                $self->{$_} = $args{$_};
            } else {
                croak("$_ is not a recognized attribute");
            }
        }

        return $self;
    }

    sub run {
        my ( $self, %args ) = @_;

        # If the 'block' argument was passed, run the code without the use of forking
        if ($args{block}) {
            my $success = run_non_web($self, %args);
            return $success;
        }

        ## path must be defined and exist to run
        if ( !defined $_[0]->{path} ) {
            croak("failed to run pipeline: no path set yet.");
        } elsif ( !-f $_[0]->{path} ) {
            croak(
                "failed to run pipeline: pipeline file $_[0]->{path} doesn't exist"
            );
        }

        ## ergatis cfg is required
        if ( !defined $args{ergatis_cfg} ) {
            croak("ergatis_cfg is a required option");
        }
        my $run_dir = $args{ergatis_cfg}->val( 'paths', 'workflow_run_dir' )
          || croak "workflow_run_dir not found in ergatis.ini";
        my $pipeline_scripts_dir = "$run_dir/scripts";

        my $current_user = $args{run_as} || '';

        ## create a directory from which to run this pipeline
        if ( -d $run_dir ) {
            ## make sure the scripts directory exists.  this is where the pipeline execution shell
            #   files are written
            if ( !-d $pipeline_scripts_dir ) {
                ( mkdir $pipeline_scripts_dir )
                  || croak
                  "filed to create pipeline scripts directory: $pipeline_scripts_dir : $!";
            }

            # make a subdirectory for this pipelineid if doing data placement
            if ( !$args{ergatis_cfg}->val( 'grid', 'vappio_data_placement' ) ) {
                $run_dir .= '/' . $self->id;

                if ( !-d $run_dir ) {
                    ( mkdir $run_dir )
                      || croak
                      "failed to create workflow_run_dir: $run_dir : $!";
                }
            }

        } else {
            croak
              "Invalid workflow_run_dir (doesn't exist) in ergatis.ini: $run_dir";
        }

        if ( !-e $args{ergatis_cfg}->val( 'paths', 'workflow_log4j' ) ) {
            croak "Invalid workflow_log4j in ergatis.ini : "
              . $args{ergatis_cfg}->val( 'paths', 'workflow_log4j' );
        }


        # If the 'web-based' code is being run, fork processes
        my $child_pid = fork;
        if ($child_pid) {
            while ( !( -e $self->path ) ) {
                sleep 3;
            }
        } else {
            ## these have to stay here or the separation doesn't work right
            close STDOUT;
            close STDERR;
            close STDIN;

            #  Fork again.  This helps separate the background process from
            #  the httpd process.  If we're in the original child, $gpid will
            #  hold the process id of the "grandchild", and if we're in the
            #  grandchild it will be zero.
            my $gpid = fork;
            if ( !$gpid ) {
                ## open the debugging file if needed
                open( my $debugfh, ">>$self->{debug_file}" ) if $self->{debug};

                ##debug
                print $debugfh "debug init\n" if $self->{debug};

                print $debugfh "attempting to chdir to $run_dir\n"
                  if $self->{debug};

                chdir $run_dir
                  || croak "Can't change to running directory $run_dir\n";
                use POSIX qw(setsid);
                setsid() or die "Can't start a new session: $!";

                print $debugfh "got past the POSIX section\n" if $self->{debug};
                $self->_setup_environment( ergatis_cfg => $args{ergatis_cfg} );
                print $debugfh "got past ENV setup section\n" if $self->{debug};

                my $marshal_interval_opt = '';
                if ( $args{ergatis_cfg}
                    ->val( 'workflow_settings', 'marshal_interval' ) )
                {
                    $marshal_interval_opt = "-m "
                      . $args{ergatis_cfg}
                      ->val( 'workflow_settings', 'marshal_interval' );
                }

                my $init_heap =
                  $args{ergatis_cfg}->val( 'workflow_settings', 'init_heap' )
                  || '100m';
                my $max_heap =
                  $args{ergatis_cfg}->val( 'workflow_settings', 'max_heap' )
                  || '1024m';

                my $sudo_prefix = '';
                my $runprefix   = '';
                my $runstring   = '';

                ## should we sudo to a different user?
                if ( $args{run_as} ) {
                    print $debugfh "INFO: run_as parameter was set\n"
                      if $self->{debug};
                    $sudo_prefix = "sudo -E -u $current_user";
                } else {
                    print $debugfh "INFO: run_as parameter not set\n"
                      if $self->{debug};
                }

                print $debugfh "so far runsprefix: $runprefix\n";

                ## are we submitting the workflow as a job?  (CURRENTLY TIED TO SGE)
                if ( $args{ergatis_cfg}
                    ->val( 'workflow_settings', 'submit_pipelines_as_jobs' ) )
                {
                    $runprefix = $args{ergatis_cfg}->val( 'grid', 'sge_qsub' )
                      . " -V -wd $run_dir -b y";

                    my $job_name = "pipeline_" . $self->id;
                    $runprefix .= " -N $job_name";

                    my $pipe_submission_queue = $args{ergatis_cfg}
                      ->val( 'workflow_settings', 'pipeline_submission_queue' );
                    my $pipe_submission_queue_memory = $args{ergatis_cfg}
                      ->val( 'workflow_settings', 'submit_queue_memory' );

                    if ($pipe_submission_queue) {
                        $runprefix .= " -q $pipe_submission_queue";
                    }

                    # If the memory has been specified add that line
                    if ($pipe_submission_queue_memory) {
                        $runprefix .=
                          " -l mem_free=$pipe_submission_queue_memory";
                    }

                    my $pipe_submission_project =
                      $args{ergatis_cfg}->val( 'workflow_settings',
                        'pipeline_submission_project' );

                    if ($pipe_submission_project) {
                        $runprefix .= " -P $pipe_submission_project";
                    }
                }

                $runstring =
                  "$ENV{'WF_ROOT'}/RunWorkflow -i $self->{path} $marshal_interval_opt "
                  . "--init-heap=$init_heap --max-heap=$max_heap ";

                ## If email notification is toggled
                if ( $args{email_user} ) {
                    $runstring .= " --notify $args{'email_user'} ";
                }

                ## Support observer scripts
                if ( $args{ergatis_cfg}
                    ->val( 'workflow_settings', 'observer_scripts' ) )
                {
                    $runstring .=
                      " --scripts="
                      . $args{ergatis_cfg}
                      ->val( 'workflow_settings', 'observer_scripts' );
                }

                $runstring .=
                    " --logconf="
                  . $args{ergatis_cfg}->val( 'paths', 'workflow_log4j' )
                  . " &> $self->{path}.run.out";

                ## write all this to a file
                my $pipeline_script =
                  "$pipeline_scripts_dir/pipeline." . $self->id . ".run.sh";
                print $debugfh "Pipeline script: $pipeline_script\n";
                open( my $pipeline_fh, ">$pipeline_script" )
                  || die "can't write pipeline shell file $pipeline_script: $!";
                print $debugfh "Wrote to pipeline script\n";

                print $pipeline_fh '#!/bin/bash', "\n\n";

                for my $env ( keys %ENV ) {
                    ## We don't want HTTP_* or REMOTE_* env variables to
                    ## be pushed to SGE.
                    if ( $env =~ /(^HTTP_|^REMOTE_)/ ) {
                        print $debugfh "UNSETTING ENV $env\n" if $self->{debug};
                        print $pipeline_fh "unset $env\n";
                        next;
                    }

                    ## don't do HTTP_ variables
                    next if $env =~ /^HTTP_/;

                    print $pipeline_fh "export $env=\"$ENV{$env}\"\n";
                    print $debugfh "ENV $env=\"$ENV{$env}\"\n"
                      if $self->{debug};
                }
                print $pipeline_fh "\n$sudo_prefix $runprefix $runstring";

                close $pipeline_fh;
                print $debugfh
                  "Wrote runstring to $pipeline_script: $runstring\n";

                ## the script needs to be executable
                chmod 0777, $pipeline_script;

                ## create a marker file showing that the pipeline has been started (or attempted to start)
#   we can't rely completely on XML here, since a pipeline submitted as a job won't have any
#   xml changes until the job starts running.  this allows us to show a 'pending' state
#   on the pipeline itself.
                my $final_run_command = "$pipeline_script";

                print $debugfh "preparing to run $final_run_command\n"
                  if $self->{debug};

                my $rc = 0xffff & system($final_run_command);

                print $debugfh
                  "system() returned %#04x: $rc for command $final_run_command\n"
                  #"system(%s) returned %#04x: $rc for command $final_run_command\n"
                  if $self->{debug};
                if ( $rc == 0 ) {
                    print $debugfh "ran with normal exit\n" if $self->{debug};
                } elsif ( $rc == 0xff00 ) {
                    print $debugfh "command failed: $!\n" if $self->{debug};
                    croak
                      "Unable to run workflow command $final_run_command failed : $!\n";
                } elsif ( ( $rc & 0xff ) == 0 ) {
                    $rc >>= 8;
                    print $debugfh "ran with non-zero exit status $rc\n"
                      if $self->{debug};
                    croak
                      "Unable to run workflow command $final_run_command failed : $!\n";
                } else {
                    print $debugfh "ran with " if $self->{debug};
                    if ( $rc & 0x80 ) {
                        $rc &= ~0x80;
                        print $debugfh "coredump from " if $self->{debug};
                    }
                    print $debugfh "signal $rc\n" if $self->{debug};
                }

                close $debugfh if $self->{debug};
            }
            exit;
        }
    }

    sub run_non_web {
        my ($self, %args) = @_;

        ## path must be defined and exist to run
        if ( !defined $_[0]->{path} ) {
            croak("failed to run pipeline: no path set yet.");
        } elsif ( !-f $_[0]->{path} ) {
            croak(
                "failed to run pipeline: pipeline file $_[0]->{path} doesn't exist"
            );
        }

        ## ergatis cfg is required
        if ( !defined $args{ergatis_cfg} ) {
            croak("ergatis_cfg is a required option");
        }
        my $run_dir = $args{ergatis_cfg}->val( 'paths', 'workflow_run_dir' )
          || croak "workflow_run_dir not found in ergatis.ini";
        my $pipeline_scripts_dir = "$run_dir/scripts";

        my $current_user = $args{run_as} || '';

        ## create a directory from which to run this pipeline
        if ( -d $run_dir ) {
            ## make sure the scripts directory exists.  this is where the pipeline execution shell
            #   files are written
            if ( !-d $pipeline_scripts_dir ) {
                ( mkdir $pipeline_scripts_dir )
                  || croak
                  "filed to create pipeline scripts directory: $pipeline_scripts_dir : $!";
            }

            # make a subdirectory for this pipelineid if doing data placement
            if ( !$args{ergatis_cfg}->val( 'grid', 'vappio_data_placement' ) ) {
                $run_dir .= '/' . $self->id;

                if ( !-d $run_dir ) {
                    ( mkdir $run_dir )
                      || croak
                      "failed to create workflow_run_dir: $run_dir : $!";
                }
            }

        } else {
            croak
              "Invalid workflow_run_dir (doesn't exist) in ergatis.ini: $run_dir";
        }

        if ( !-e $args{ergatis_cfg}->val( 'paths', 'workflow_log4j' ) ) {
            croak "Invalid workflow_log4j in ergatis.ini : "
              . $args{ergatis_cfg}->val( 'paths', 'workflow_log4j' );
        }

        ## open the debugging file if needed
        open( my $debugfh, ">>$self->{debug_file}" ) if $self->{debug};

        ##debug
        print $debugfh "debug init\n" if $self->{debug};

        print $debugfh "attempting to chdir to $run_dir\n"
          if $self->{debug};

        chdir $run_dir
          || croak "Can't change to running directory $run_dir\n";

        $self->_setup_environment( ergatis_cfg => $args{ergatis_cfg} );
        print $debugfh "got past ENV setup section\n" if $self->{debug};

        my $marshal_interval_opt = '';
        if ( $args{ergatis_cfg}
            ->val( 'workflow_settings', 'marshal_interval' ) )
        {
            $marshal_interval_opt = "-m "
              . $args{ergatis_cfg}
              ->val( 'workflow_settings', 'marshal_interval' );
        }

        my $init_heap =
          $args{ergatis_cfg}->val( 'workflow_settings', 'init_heap' )
          || '100m';
        my $max_heap =
          $args{ergatis_cfg}->val( 'workflow_settings', 'max_heap' )
          || '1024m';

        my $sudo_prefix = '';
        my $runprefix   = '';
        my $runstring   = '';

        ## should we sudo to a different user?
        if ( $args{run_as} ) {
            print $debugfh "INFO: run_as parameter was set\n"
              if $self->{debug};
            $sudo_prefix = "sudo -E -u $current_user";
        } else {
            print $debugfh "INFO: run_as parameter not set\n"
              if $self->{debug};
        }

        print $debugfh "so far runsprefix: $runprefix\n";

        ## are we submitting the workflow as a job?  (CURRENTLY TIED TO SGE)
        if ( $args{ergatis_cfg}
            ->val( 'workflow_settings', 'submit_pipelines_as_jobs' ) )
        {
            $runprefix = $args{ergatis_cfg}->val( 'grid', 'sge_qsub' )
              . " -V -wd $run_dir -b y";

            my $job_name = "pipeline_" . $self->id;
            $runprefix .= " -N $job_name";

            my $pipe_submission_queue = $args{ergatis_cfg}
              ->val( 'workflow_settings', 'pipeline_submission_queue' );
            my $pipe_submission_queue_memory = $args{ergatis_cfg}
              ->val( 'workflow_settings', 'submit_queue_memory' );

            if ($pipe_submission_queue) {
                $runprefix .= " -q $pipe_submission_queue";
            }

            # If the memory has been specified add that line
            if ($pipe_submission_queue_memory) {
                $runprefix .=
                  " -l mem_free=$pipe_submission_queue_memory";
            }

            my $pipe_submission_project =
              $args{ergatis_cfg}->val( 'workflow_settings',
                'pipeline_submission_project' );

            if ($pipe_submission_project) {
                $runprefix .= " -P $pipe_submission_project";
            }
        }

        $runstring =
          "$ENV{'WF_ROOT'}/RunWorkflow -i $self->{path} $marshal_interval_opt "
          . "--init-heap=$init_heap --max-heap=$max_heap ";

        ## If email notification is toggled
        if ( $args{email_user} ) {
            $runstring .= " --notify $args{'email_user'} ";
        }

        ## Support observer scripts
        if ( $args{ergatis_cfg}
            ->val( 'workflow_settings', 'observer_scripts' ) )
        {
            $runstring .=
              " --scripts="
              . $args{ergatis_cfg}
              ->val( 'workflow_settings', 'observer_scripts' );
        }

        $runstring .=
            " --logconf="
          . $args{ergatis_cfg}->val( 'paths', 'workflow_log4j' )
          . " &> $self->{path}.run.out";

        ## write all this to a file
        my $pipeline_script =
          "$pipeline_scripts_dir/pipeline." . $self->id . ".run.sh";
        print $debugfh "Pipeline script: $pipeline_script\n";
        open( my $pipeline_fh, ">$pipeline_script" )
          || die "can't write pipeline shell file $pipeline_script: $!";
        print $debugfh "Wrote to pipeline script\n";

        print $pipeline_fh '#!/bin/bash', "\n\n";

        for my $env ( keys %ENV ) {
            ## We don't want HTTP_* or REMOTE_* env variables to
            ## be pushed to SGE.
            if ( $env =~ /(^HTTP_|^REMOTE_)/ ) {
                print $debugfh "UNSETTING ENV $env\n" if $self->{debug};
                print $pipeline_fh "unset $env\n";
                next;
            }

            ## don't do HTTP_ variables
            next if $env =~ /^HTTP_/;

            print $pipeline_fh "export $env=\"$ENV{$env}\"\n";
            print $debugfh "ENV $env=\"$ENV{$env}\"\n"
              if $self->{debug};
        }
        print $pipeline_fh "\n$sudo_prefix $runprefix $runstring";

        close $pipeline_fh;
        print $debugfh
          "Wrote runstring to $pipeline_script: $runstring\n";

        ## the script needs to be executable
        chmod 0777, $pipeline_script;

## create a marker file showing that the pipeline has been started (or attempted to start)
#   we can't rely completely on XML here, since a pipeline submitted as a job won't have any
#   xml changes until the job starts running.  this allows us to show a 'pending' state
#   on the pipeline itself.
        my $final_run_command = "$pipeline_script";

        print $debugfh "preparing to run $final_run_command\n"
          if $self->{debug};

		# Run in the background
		$final_run_command .= " &";
		system($final_run_command);
		
    	if ( $? == -1 ) {
        	croak "failed to execute command ($final_run_command): $!\n";
    	} elsif ( $? & 127 ) {
         	my $out = sprintf "command ($final_run_command): child died with signal %d, %s coredump\n",
                    ($? & 127),  ($? & 128) ? 'with' : 'without';
        	croak ($out);
    	}

		# Program should have went through if $? = 0
        print $debugfh "[$final_run_command] ran with normal exit\n" if $self->{debug};

        close $debugfh if $self->{debug};

        # Create a progress bar to visually track progress
        my $p_bar = create_progress_bar($self->{path}, $self->{id}) if $args{show_progress};

       # If 'block' is set to 1, wait for pipeline to return non-running state
        my $p_state = '';
        my $running_components = ();
        my $component_list = ();
		my %prev_component_states = ();
        do {
            $p_state = $self->pipeline_state;
   			# First iteration will be undefined... populate hash with previous states
			if (defined $component_list) {
                foreach my $component (keys %$component_list) {
                    $prev_component_states{$component} = $component_list->{$component}->{'state'};
                }
			}
            $component_list = build_twig($self->{path});
			if ($args{show_progress}) {
			    update_progress_bar($p_bar, $component_list);
		    } else {
				if ( %prev_component_states ) {
                    handle_component_status_changes($component_list, \%prev_component_states);
                }
		    }
            sleep 60 if ( $p_state =~ /(running|pending|waiting|incomplete)/ );
        } while ( $p_state =~ /(running|pending|waiting|incomplete)/ );

        # If end-state is complete, return 1.  Otherwise return 0
        return 1 if ($p_state eq 'complete');
        return 0;
    }

    ## accessors
    sub id   { return $_[0]->{id} }
    sub path { return $_[0]->{path} }

    ## modifiers
    sub set_id   { $_[0]->{id}   = $_[1] }
    sub set_path { $_[0]->{path} = $_[1] }

    sub _setup_environment {
        my ( $self, %args ) = @_;

        ## remove the apache SERVER variables from the environment
        for my $k ( keys %ENV ) {
            if ( $k =~ /^SERVER_/ ) {
                delete $ENV{$k};
            }

            if ( $k =~ /^BASH_FUNC_module/ ) {
                delete $ENV{$k};
            }
        }

        ## these variable seemed to have been causing SGE problems (bug 4565)
        delete $ENV{MC};
        delete $ENV{BASH_FUNC_module};

        $ENV{SGE_ROOT} = $args{ergatis_cfg}->val( 'grid', 'sge_root' );
        $ENV{SGE_CELL} = $args{ergatis_cfg}->val( 'grid', 'sge_cell' );

  #$ENV{SGE_QMASTER_PORT} = $args{ergatis_cfg}->val('grid', 'sge_qmaster_port');
        $ENV{SGE_EXECD_PORT} =
          $args{ergatis_cfg}->val( 'grid', 'sge_execd_port' );
        $ENV{SGE_ARCH} = $args{ergatis_cfg}->val( 'grid', 'sge_arch' );

        ## these WF_ definitions are usually kept in the $workflow_root/exec.tcsh file,
        #   which we're not executing.
        $ENV{WF_ROOT} = $args{ergatis_cfg}->val( 'paths', 'workflow_root' );
        $ENV{WF_ROOT_INSTALL} = $ENV{WF_ROOT};
        $ENV{WF_TEMPLATE}     = "$ENV{WF_ROOT}/templates";

        #$ENV{SYBASE} = '/usr/local/packages/sybase';
        my $sge_bin = $args{ergatis_cfg}->val( 'grid', 'sge_root' ) . '/bin/';
        $ENV{PATH} =
          "$ENV{WF_ROOT}:$ENV{WF_ROOT}/bin:$ENV{WF_ROOT}/add-ons/bin:$sge_bin/$ENV{SGE_ARCH}:$ENV{PATH}";

        $ENV{LD_LIBRARY_PATH} = '';
        $ENV{TERMCAP} =
          '';    #can contain bad characters which crash wrapper shell script
        ## some application-specific env vars
        ############

        ## for the htab.pl script within the hmmpfam component
        $ENV{HMM_SCRIPTS} = '/usr/local/devel/ANNOTATION/hmm/bin';

        ## for the open-source hmmpfam, this ensures that it only runs on a single processor
        $ENV{HMMER_NCPU} = 1;

        ## for the genewise component
        $ENV{WISECONFIGDIR} =
          '/usr/local/devel/ANNOTATION/EGC_utilities/WISE2/wise2.2.0/wisecfg';

        ## for local data placement
		if (defined $args{ergatis_cfg}->val( 'grid', 'vappio_root' )){
            $ENV{vappio_root} = $args{ergatis_cfg}->val( 'grid', 'vappio_root' );
            $ENV{vappio_data_placement} =
                $args{ergatis_cfg}->val( 'grid', 'vappio_data_placement' ); 
		}
        $ENV{PERL5LIB} =
          "/usr/local/packages/perllib/x86_64-linux-thread-multi:$ENV{PERL5LIB}" if (defined $ENV{PERL5LIB});

        ## for overwriting SGEs default
        $ENV{TMPDIR} = "/tmp";

    }

    # Returns the pipeline state of the passed-in XML
    sub pipeline_state {
        my $self = shift;
        my $state;
        open( IN, "<" . $self->path )
          or die("Cannot open pipeline xml for reading $self->{path} ($!)");
        while (<IN>) {

            # Should be first occurrence (outermost nesting)
            if (/<state>(\w+)<\/state>/) {
                $state = $1;
                last;
            }
        }
        close(IN);
        return $state;
    }

}
1 == 1;
