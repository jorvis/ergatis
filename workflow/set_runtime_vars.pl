#!/usr/local/bin/perl

# set_runtime_vars is a utility for replacing placeholder text strings (keys) with values
# as specified in a simple configuration file of the form...
#
# ;key;=value
#
#

sub help
{
    print " set_runtime_vars is a utility for replacing placeholder text strings (keys) with values\n";
    print " as specified in a simple configuration file of the form...\n\n";
    print " ;key;=value\n\n";

    print "Usage: ./set_runtime_vars.pl -c config_file < INPUT > OUTPUT\n\n";
    exit();
}

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %options = ();
my $results = GetOptions (\%options, 'config|c=s', 'help|h' );

if( $options{'help'} ){ help(); }

open( CONFIGFILE, "$options{'config'}" ) or die "Could not open config file ($options{'config'})\n";

my $subs = {};

while( my $line = <CONFIGFILE> )
{
    chomp $line;
    my( $key, $value ) = split( '=', $line );

    $subs->{$key}=$value;
}

while( my $line = <STDIN> )
{
    foreach my $key (keys(%{$subs}))
    {
	$line =~ s/$key/$subs->{$key}/g;
    }

    print $line;
}
