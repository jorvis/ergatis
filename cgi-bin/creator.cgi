#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Config::IniFiles;
use HTML::Template;

my $q = new CGI;

print $q->header( -type => 'text/html' );

## get the configuration
my $sys_config = Config::IniFiles->new( -file => 'ergatis.ini' );

my %conf = load_shared_conf( $sys_config->val('path', 'workflow_conf') );

## present the login form
my $tmpl = HTML::Template->new( filename => 'templates/creator.tmpl',
                                             die_on_bad_params => 1
                               );

my $components = [];
for ( glob "$conf{SCHEMA_DIR}/*conf.ini" ) {
    m|.+/(.+)conf.ini|;
    my $name = $1;
    my $label = $name;
    $label =~ s/_/ /g;
    
    ## handle any custom label changes here
    $label =~ s/coordinates/coords/g;
    
    push @$components, { 
#                            name  => $name,
                            label => $label,
                       };
}

$tmpl->param( COMPONENTS => $components );

print $tmpl->output;

exit(0);


sub load_shared_conf {
    my $path = shift;
    my %cf;
    
    open(my $ifh, "<$path") || die "couldn't read shared conf file: $!";
    
    while (<$ifh>) {
        if (/^\$\;([A-Z_\-]+)\$\;\s*\=\s*(.+)$/) {
            $cf{$1} = $2;
        }
    }
    
    return %cf;
}
