#!/usr/bin/perl

=head1 NAME

create_release_externals.pl - Creates a 'release' directory with svn:externals definitions

=head1 SYNOPSIS

USAGE: template.pl 
            --svn_repository=https://ergatis.svn.sourceforge.net/svnroot/ergatis
            --new_release_tag=ergatis-v2r10b1
          [ --source_release_tag=trunk 
            --log=/path/to/some.log
            --debug=4
          ]

=head1 OPTIONS

B<--svn_repository>
    The root of the svn repository for your project under which we should find
    the 'tags', 'release' and 'trunk' directories.

B<--new_release_tag>
    This is the name of the tag to be created, such as ergatis-v2r10b1

B<--source_release_tag>
    Optional.  This is a release based on which tag?  If not 'trunk' this should
    be a branch label.  (default = trunk)

B<--log,-l> 
    Log file

B<--debug> 
    Debug level.

B<--help,-h>
    This help message

=head1  DESCRIPTION

steps to doing a full release from trunk with example version 'N':

- svn copy trunk to tags/ergatis-N
- mkdir release/ergatis-N
- svn propset svn:externals release/ergatis-N  [ definition file ]

steps to doing a full release from branch with example version N
- svn copy branches/$tag to tags/ergatis-N
- mkdir release/ergatis-N
- svn propset svn:externals release/ergatis-N [ definition file ]

steps to create release externals from an existing tag

- mkdir release/ergatis-N
- svn propset svn:externals release/ergatis-N  [ definition file ]

=head1  INPUT

If the values for --new_release_tag and --source_release_tag are the same, the script won't do
any code copies and will only create the svn:externals definitions under the release directory.
This is probably only needed when the code had been tagged previously by using direct SVN
commands.

=head1  OUTPUT

Just a log - the rest are SVN changes.  Does create a temporary file with the svn:externals
definition within the /tmp directory.


=head1  CONTACT

    Joshua Orvis
    jorvis@users.sf.net

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Ergatis::Logger;
use Pod::Usage;

my %options = ();
GetOptions(\%options, 
           'svn_repository=s',
           'new_release_tag=s',
           'source_release_tag=s',
           'log|l=s',
           'debug=s',
           'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

## Setup the logger.  See perldoc for more info on usage
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE' =>$logfile,
                                 'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

## set this variable to 1 if you just want to log the commands that would have been run
my $testing = 0;

## this is reused
my $cmd;

## if the values for the new and source release are different we need to do an svn copy
if ( $options{new_release_tag} ne $options{source_release_tag} ) {
    my ($src_url, $dest_url);
    
    if ( $options{source_release_tag} =~ /^head$/i || 
         $options{source_release_tag} =~ /^trunk$/i   ) {

        $src_url = "$options{svn_repository}/trunk";

    } else {
        $src_url = "$options{svn_repository}/branches/$options{source_release_tag}";
    }

    $dest_url = "$options{svn_repository}/tags/$options{new_release_tag}";

    $cmd = "svn copy -m \"$options{new_release_tag} release based off $options{source_release_tag}\" $src_url $dest_url";
    run_cmd($cmd, "There was a problem tagging version $options{new_release_tag} from version $options{source_release_tag}");
}

## now create the release directory
$cmd = "svn mkdir -m \"creating release directory for $options{new_release_tag}\" $options{svn_repository}/release/$options{new_release_tag}";
run_cmd($cmd, "There was a problem creating the release/$options{new_release_tag} directory");

## we need to write the properties file
$logger->info("writing svn:externals file to /tmp/$$.externals");
open(my $externals_fh, ">/tmp/$$.externals") || die "failed to create externals file: $!";

print $externals_fh <<SVN_exteRnAls;

$options{new_release_tag} $options{svn_repository}/tags/$options{new_release_tag}/install
$options{new_release_tag}/bin $options{svn_repository}/tags/$options{new_release_tag}/src/perl
$options{new_release_tag}/shell $options{svn_repository}/tags/$options{new_release_tag}/src/shell
$options{new_release_tag}/lib $options{svn_repository}/tags/$options{new_release_tag}/lib
$options{new_release_tag}/components $options{svn_repository}/tags/$options{new_release_tag}/components
$options{new_release_tag}/htdocs $options{svn_repository}/tags/$options{new_release_tag}/htdocs
$options{new_release_tag}/src $options{svn_repository}/tags/$options{new_release_tag}/src/c

SVN_exteRnAls

## can't currently set svn:externals on a remote directory.  check it out first.
$cmd = "svn co $options{svn_repository}/release/$options{new_release_tag} /tmp/$options{new_release_tag}";
run_cmd($cmd, "There was a problem checking out the newly-created release directory.  See the logs for the last command attempted");

## now set the properties
$cmd = "svn propset svn:externals /tmp/$options{new_release_tag}/ -F /tmp/$$.externals";
run_cmd($cmd, "Failed to set svn:externals.  Check log for last executed statement");

## commit the properties
$cmd = "svn commit -m \"set svn:externals on release $options{new_release_tag}\" /tmp/$options{new_release_tag}";
run_cmd($cmd, "Failed to commit svn:externals.  CHeck log for last executed statement");


exit(0);

sub run_cmd {
    my ($command, $err_msg) = @_;

    $logger->info("running: $command");

    if (! $testing ) {
        my $cmd_stdout = `$command`;

        if ( $? ) {
            print STDERR "$err_msg\n";
            exit(1);
        }
    }
}


sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( new_release_tag svn_repository );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
    
    ## the new tag directly shouldn't exist already in /tmp
    if ( -e "/tmp/$options{new_release_tag}" ) {
        die "/tmp/$options{new_release_tag} is in the way.  please remove it.";
    }
        
    
    ## handle some defaults
    $options{source_release_tag} = 'trunk' unless ($options{source_release_tag});
}
