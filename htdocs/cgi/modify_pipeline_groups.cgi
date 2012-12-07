#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::ConfigFile;
use HTML::Template;

my $q = new CGI;

my $ids_passed = $q->param('pipeline_ids');
my $action = $q->param('action');
my $repository_root = $q->param('repository_root');

my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

my $username = user_logged_in($ergatis_cfg);
my $auth_method = $ergatis_cfg->val('authentication', 'authentication_method');
unless ($auth_method eq 'open' || defined($username)) {
    print $q->header( -type => 'text/html' );
    print_error_page( ergatis_cfg => $ergatis_cfg,
                      message => "You must be logged in to modify pipeline groups",
                      links => [
                                 { label => "pipeline list", is_last => 1, url => "./pipeline_list.cgi?repository_root=$repository_root&view=group" },
                               ]
                    );
    exit(0);
}

unless ( $ids_passed && $action && $repository_root ) {
    print $q->header( -type => 'text/html' );
    
    print_error_page( ergatis_cfg => $ergatis_cfg,
          message => "There was a problem putting these pipelines into the groups specified. " .
                     "This should not happen, please file a bug report at http://sf.net/projects/ergatis",
          links => [ 
                        { label => "bug report page", 
                          is_last => 1, 
                          url => "http://sf.net/projects/ergatis" },
                   ],
    );
    
    exit();
}

## if we get to here all needed options were passed
my $pipeline_root = "$repository_root/workflow/runtime/pipeline";

my @ids_passed = split(' ', $ids_passed);

for my $id ( @ids_passed ) {
    next unless ($id =~ /^[A-Z0-9]+$/);
    
    if ( -d "$pipeline_root/$id" ) {

        ## start by reading any existing groups    
        my %existing_groups = ();
        my $groupsfh = get_conditional_read_fh("$pipeline_root/$id/pipeline.xml.groups");
        
        if ( $groupsfh ) {
            while (<$groupsfh>) {
                chomp;
                $existing_groups{$_}++;
            }

            close $groupsfh;
        }
    
        
        my @groups_entered = split(',', $q->param('group_labels'));

        for my $group ( @groups_entered ) {

            ## no whitespace on either side of the group
            $group =~ /^\s*(.+?)\s*$/;
            $group = $1;
            
            ## now, are we adding to any existing groups or taking away?
            if ( $action eq 'add' ) {
                $existing_groups{$group}++;
            } elsif ( $action eq 'remove' ) {
                delete $existing_groups{$group};
            } else {
                croak("unrecognized action: $action");
            }
            
        }
        
        ## write the groups out
        open( my $ofh, ">$pipeline_root/$id/pipeline.xml.groups" ) || croak("failed to write groups file: $!");
        for (keys %existing_groups) {
            print $ofh "$_\n";
        }

        
    } else {
        print STDERR "failed to find pipeline directory $pipeline_root/$id";
    }
}

## redirect back to the page the user was editing
print $q->redirect( -uri => url_dir_path($q) . "pipeline_list.cgi?view=group&repository_root=" . $repository_root );
