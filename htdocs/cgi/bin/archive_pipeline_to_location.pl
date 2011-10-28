#!/usr/bin/perl -w

=head1  NAME 

archive_pipeline_to_location.pl - archive a workflow pipeline and, optionally, its output data

=head1 SYNOPSIS

USAGE: archive_pipeline_to_location.pl 
        --pipeline_id=1234
        --repository_root=/usr/local/annotation/AA1
        --archive_root=/usr/local/scratch/ergatis/archival/AA1
      [ --lock_file=/path/to/somefile.lock
        --archive_output=1
        --log=/path/to/some.log
      ]

=head1 OPTIONS

B<--pipeline_id,-i> 
    the ID of the pipeline to archive

B<--repository_root,-r> 
    the project directory, just under which we should find the Workflow directory.

B<--archive_root,-a> 
    root of the directory where archivals will be stored.  it will be created if it
    doesn't exist, along with the necessary subdirectory structure.

B<--lock_file,-k> 
    optional.  This file will be created at the start of the run and will contain
    only the word "archiving".  Any existing file will get stomped.  It will be
    deleted once processing is finished.

B<--archive_output,-o> 
    optional.  If passed, the output in the standard output directory structure
    for this pipeline will also be compressed.

B<--log,-l> 
    optional.  will create a log file with summaries of all actions performed.

B<--help,-h> 
    This help message/documentation.

=head1   DESCRIPTION

This script is used to archive a pipeline and, optionally, all its associated output.  
This is based on pipeline_id and repository root.  The contents of each the following
folders will be compressed and migrated to the archive root as individual tarballs:

    $repository_root/workflow/runtime/*/$pipelineid_*
    $repository_root/output_repository/*/$pipelineid_*
    $repository_root/workflow/runtime/pipeline/$pipelineid


=head1 INPUT

The input is defined with the --repository_root, --pipeline_id and --archive_root options.

=head1 OUTPUT

This is a archival script.  There is no output unless you use the --log option.

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
              'archive_root|a=s',
              'lock_file|k=s',
              'archive_output|o=i',
              'log|l=s',
			  'help|h') || pod2usage();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

## play nicely
umask(0000);

## make sure all passed options are peachy
&check_parameters(\%options);

## create the lock file if requested
my $lockfh;
if ( $options{lock_file} ) {
    open($lockfh, ">$options{lock_file}") || die "failed to create lock file $options{lock_file}: $!";
    print $lockfh "archiving";
}

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">$options{log}") || die "can't create log file: $!";
    
    ## i'm setting autoflush here so that logs are written immediately
    ##  i was getting terrible lag otherwise.
    autoflush $logfh 1;
}

my $wf_dir = "$options{repository_root}/workflow/runtime";
opendir( my $wfdh, $wf_dir );

## get the things in the workflow directory (except the pipeline)
_log("INFO: scanning $wf_dir");
for my $component ( readdir $wfdh ) {
    _log("INFO: scanning $wf_dir/$component");
    
    for my $rundir ( glob "$wf_dir/$component/$options{pipeline_id}_*" ) {
        &archive_to_location($rundir);
    }
}

## the user can optionally archive the output repository too
if ( $options{archive_output} ) {
    my $o_dir  = "$options{repository_root}/output_repository";
    opendir( my $odh, $o_dir );

    for my $component ( readdir $odh ) {
        for my $rundir ( glob "$o_dir/$component/$options{pipeline_id}_*" ) {
            &archive_to_location($rundir);
        }
    }
}

## archive the pipeline xml explicitly
&archive_to_location("$wf_dir/pipeline/$options{pipeline_id}");

## remove the lock file
if ( $options{lock_file} ) {
    close $lockfh;
    unlink $options{lock_file};
}

_log("INFO: archive operations complete");

exit(0);


sub _log {
    my $msg = shift;
    
    print $logfh "$msg\n" if $logfh;
}

sub check_parameters {
    
    ## pipeline_id and repository_root are required
    unless ( defined $options{pipeline_id} && $options{repository_root} && $options{archive_root} ) {
        print STDERR "pipeline_id, repository_root and archive_root options are required\n\n";
        pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
    }
    
    ## make sure the repository root has a workflow directory
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

sub archive_to_location {
    my $source = shift;
    
    ## take off any trailing / from the source
    if ( $source =~ m|(.+)/*$| ) {
        $source = $1;
    }
    
    my $arch_dir = "$options{archive_root}/$source";

    $arch_dir =~ m|^(.+)/.+|;
    my $arch_pre_dir = $1;

    &create_directory($arch_pre_dir);

    _log("INFO: recursively archiving $source to $arch_dir.tar.gz");
    
#    _log("INFO: tar -czf $arch_dir.tar.gz $source");
    
    my $ret = system("tar -czf $arch_dir.tar.gz $source");
#    $ret >>= 8;
#    if ( $ret ) {
#        _log( "FATAL: not removing $source because tar returned $ret" );
#    } else {
        _log( "INFO: recursively removing $source" );
        system("rm -rf $source");
#    }
}

sub create_directory{
    my($dir) = @_;
    
    _log("INFO: creating dir $dir");
    print "INFO: creating dir $dir\n";
    
    my $ret = system("mkdir -m 777 -p $dir");
    $ret >>= 8;
    return $ret;
}

