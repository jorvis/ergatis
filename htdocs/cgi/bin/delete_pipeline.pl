#!/usr/bin/perl -w

=head1  NAME

delete_pipeline.pl - delete a workflow pipeline and, optionally, its output data

=head1 SYNOPSIS

USAGE: delete_pipeline.pl
        --pipeline_id=1234
        --repository_root=/usr/local/annotation/AA1
      [ --lock_file=/path/to/somefile.lock
        --delete_output=1
        --log=/path/to/some.log
      ]

=head1 OPTIONS

B<--pipeline_id,-i>
    The ID of the pipeline to delete

B<--repository_root,-r>
    The project directory, just under which we should find the asmbls directory.

B<--lock_file,-k>
    optional.  This file will be created at the start of the run and will contain
    only the word "deleting".  Any existing file will get stomped.  It will be
    deleted once process is finished.

B<--delete_output,-o>
    optional.  If passed, the output in the standard output directory structure
    for this pipeline will also be removed.

B<--exclusion_file, -e>
    optional. If passed, the components listed within this file will not be deleted
    along with any parent directories containing this component.  Each line should
	be in the format of component_name.token, like ncbi-blastx.pre_overlap_analysis

B<--dry_run, -d>
	optional.  When enabled, will show what contents will be removed but will not actually remove them.

B<--log,-l>
    optional.  will create a log file with summaries of all actions performed.

B<--help,-h>
    This help message/documentation.

=head1   DESCRIPTION

This script is used to delete a pipeline and all its associated output.  This
is based on pipeline_id and repository root.  Deleted folders include (in order):

    $repository_root/workflow/runtime/*/$pipelineid_*
    $repository_root/output_repository/*/$pipelineid_*
    $repository_root/workflow/runtime/pipeline/$pipelineid


=head1 INPUT

The input is defined with the --repository_root and --pipeline_id options.

=head1 OUTPUT

This is a deletion script.  There is no output unless you use the --log option.

=head1 CONTACT

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;
use File::Find;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use IO::Handle;
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options,
			  'pipeline_id|i=s',
			  'repository_root|r=s',
			  'exclusion_file|e=s',
			  'lock_file|k=s',
			  'delete_output|o=i',
			  'dry_run|d',
			  'log|l=s',
			  'help|h') || pod2usage();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

my $kept_dir_paths;
my $out_flag = 0;
my $wf_flag = 0;

## play nicely
umask(0000);

## make sure all passed options are peachy
&check_parameters(\%options);

## small sanity check
if ( $options{'pipeline_id'} =~ /[^A-Z0-9]/i ) {
    die "ERROR: encountered a pipeline ID with non-alphanumeric characters.  Cowardly refusing to proceed with requested delete.";
}

## create the lock file if requested
my $lockfh;
my $lock;
if ( $options{lock_file} ) {
	$lock = $options{lock_file};
	if (-e $lock) {
	    _log("Another job is deleting pipeline $options{'pipeline_id'}");
	    mail_to();
	    exit(3); 	#exit if another process is already running
	}
    open($lockfh, ">$lock") || die "failed to create lock file $lock: $!";
    print $lockfh "deleting";
	chmod 0777, $lock;	# in case 3rd-party user needs to remove
	close $lockfh;
}

$SIG{INT} = \&capture_ctrl_c;

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">$options{log}") || die "can't create log file: $!";

    ## i'm setting autoflush here so that logs are written immediately
    ##  i was getting terrible lag otherwise.
    autoflush $logfh 1;
}

&_log("INFO: Dry run mode enabled") if ($options{'dry_run'});

if ($options{'exclusion_file'}) {
    $kept_dir_paths = parse_exclusions($options{'exclusion_file'}, $options{'pipeline_id'});
}

my $wf_dir = "$options{repository_root}/workflow/runtime";
opendir( my $wfdh, $wf_dir );

## get the things in the workflow directory (except the pipeline)
_log("INFO: scanning $wf_dir");
for my $component ( readdir $wfdh ) {
    _log("INFO: scanning $wf_dir/$component");

    for my $rundir ( glob "$wf_dir/$component/$options{pipeline_id}_*" ) {
		$wf_flag = 0;
		# See if this directory matches an excluded component
        foreach my $dir (@$kept_dir_paths){
			$wf_flag++ if ($rundir eq "$wf_dir/$dir");
			last;
		}
        _log("INFO: recursively removing $rundir") if (!$wf_flag);
		next if ($options{'dry_run'});
        &find_remove($rundir) if (!$wf_flag);

    }
}

## the user can optionally delete the output repository too
if ( $options{delete_output} ) {
    my $o_dir  = "$options{repository_root}/output_repository";
    opendir( my $odh, $o_dir );

    for my $component ( readdir $odh ) {
        for my $rundir ( glob "$o_dir/$component/$options{pipeline_id}_*" ) {
            $out_flag = 0;
			# See if this directory matches an excluded component
        	foreach my $dir (@$kept_dir_paths){
				$out_flag++ if ($rundir eq "$o_dir/$dir");
				last;
			}
			_log("INFO: recursively removing $rundir") if (!$out_flag);
			next if ($options{'dry_run'});
            &find_remove($rundir) if (!$out_flag);
        }
    }
}

## delete the pipeline xml explicitly
_log("INFO: recursively removing $wf_dir/pipeline/$options{pipeline_id}");
&find_remove("$wf_dir/pipeline/$options{pipeline_id}") if (! $options{'dry_run'});

## remove the lock file
unlink $lock if ($lock);

_log("INFO: delete operations complete");

exit(0);

#  Parse exclusions file to get the directories that will be spared from deletion
sub parse_exclusions {
    my ($exclusions, $pipe_id) = @_;
    my @kept_dirs;

	open EFH, $exclusions or die "Can't open exclusion file $exclusions for reading: $!\n";
	while (<EFH>) {
		my $line = $_;
		chomp $line;
		my ($component, $token) = split(/\./, $line);
		print $line , "\n";
		my $dir = $component . "/" . $pipe_id . "_" . $token;
		&_log("INFO:  keeping directory $dir from deletion");
		push @kept_dirs, $dir;
	}
	close EFH;
    return \@kept_dirs;
}

# Handle cases where the deletion process was interrupted
sub capture_ctrl_c {
    _log("Ctrl-C was hit.  Archiving will stop");
    _log("removing lock $lock...");
    unlink($lock) if ($lock);
    mail_to();
    exit(2);
}

sub mail_to {
    my (undef, $mon, $day, undef, $year) = split(/\s+/, localtime);
    my $subject = "Pipeline " . $options{'pipeline_id'} . "--" . $mon.'_'.$day.'_'.$year . " Deletion Error Report";
    my $email = 'sadkins@som.umaryland.edu';

    if (defined $logfh) {
   	 my $mail_cmd = "mail -s \"$subject\" $email < $options{'log'}";
   	 system($mail_cmd);
    }

}

sub _log {
    my $msg = shift;
    print "$msg\n";
    print $logfh "$msg\n" if $logfh;
}

sub check_parameters {

    ## pipeline_id and repository_root are required
    unless ( defined $options{pipeline_id} && $options{repository_root} ) {
        print STDERR "pipeline_id and repository_root options are required\n\n";
        pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
    }

    ## make sure the repository root has a Workflow directory
    unless ( -d "$options{repository_root}/workflow/runtime" ) {
        print STDERR "\n\nthe repository root passed doesn't contain a workflow/runtime directory.\n\n";
        exit(1);
    }

    ## make sure the repository root has an output_repository directory
    unless ( -d "$options{repository_root}/output_repository" ) {
        print STDERR "\n\nthe repository root passed doesn't contain an output_repository directory.\n\n";
        exit(1);
    }
}

sub find_remove {
    my $dir = shift;

    finddepth(
        {
            follow => 0,
            no_chdir => 1,
            wanted => \&remove_recursively
        }, $dir
    );
}

sub remove_recursively {
    my $thing = $File::Find::name;

    ## is it a file?
    if (-f $thing) {
        unless ( unlink($thing) ) {
            _log("failed to unlink $thing");
        }
    ## is it a directory?
    } elsif (-d $thing) {
        unless ( rmdir($thing) ) {
            _log("failed to rmdir $thing");
        }
    } else {
        _log("ERROR: failed to delete $thing because I don't know what it is.");
		unlink($lock) if ($lock)
    }

}
