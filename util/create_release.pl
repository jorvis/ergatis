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
    Optional.  This is a release based on which tag?  (default = trunk)

B<--log,-l> 
    Log file

B<--debug> 
    Debug level.

B<--help,-h>
    This help message

=head1  DESCRIPTION

steps to doing a full release with example version 'N':

- svn copy trunk to tags/ergatis-N
- mkdir release/ergatis-N
- svn propset svn:externals release/ergatis-N  [ definition file ]

steps to create release externals from an existing tag

- mkdir release/ergatis-N
- svn propset svn:externals release/ergatis-N  [ definition file ]

=head1  INPUT

If the values for --new_release_tag and --source_release_tag are the same, the script won't do
any code copies and will only create the svn:externals definitions under the release directory.
This is probably only needed when the code had been tagged previously by using direct SVN
commands.

=head1  OUTPUT

Just a log - the rest are SVN changes

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

## if the values for the new and source release 


##
## CODE HERE
##

exit(0);


sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( new_release_tag svn_repository );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
    
    ##
    ## you can do other things here, such as checking that files exist, etc.
    ##
    
    ## handle some defaults
    $options{source_release_tag} = 'trunk' unless ($options{source_release_tag});
}
