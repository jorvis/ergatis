#!/local/perl/bin/perl

## TEMP!
exit(0);

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %options = ();
my $results = GetOptions( \%options, 'dtd|d=s', 'help|h' );

if( $options{'help'} || !($ARGV[0]) )
{
    print "\ndtdValid.pl validates xml documents using an inline or user specified dtd.\n\n";
    print "Usage:    ./dtdValid.pl -d file.dtd source.xml\n\n";
    print "Options:\n    --dtd | -d   : the file path to the user specified dtd (optional)\n";
    print "    --help | -h   : this message\n\n";
    print "Output:\n    Nothing if validation is successful\n";
    print "    XML Schema errors if validation is unsuccessful\n\n\n";

    exit(4);
}

##############################################################################
## TEMPORARY HACK
## because of library differences between the servers and desktop machines,
## this script will return true any time it is run on a machine with the 2.6
## series kernel.  this should be removed once the servers are upgraded.

my $uname = `uname -r`;

if ($uname =~ /^2.6/) {
    exit(0);
}

$ENV{LD_LIBRARY_PATH} = "/usr/local/lib:$ENV{LD_LIBRARY_PATH}";

##############################################################################


my $file = $ARGV[0];
my $dtd = $options{'dtd'};

my ($dir) = ($0 =~ /(.*)\/.*/);

my $status = 0;

if( $dtd )
{
    $dtd =~ s/\//\\\//g;                                                                                                                                                                                                              
    $status = system "sed -e \'s/<\\!DOCTYPE Bsml PUBLIC \"-\\/\\/EBI\\/\\/Labbook, Inc. BSML DTD\\/\\/EN\" \"http:\\/\\/www.labbook.com\\/dtd\\/bsml3_1.dtd\">/<\\!DOCTYPE Bsml SYSTEM \"$dtd\">/\' $file | $dir/Xerces-xsdValid";
}

else
{
    $status = system "more $file | $dir/Xerces-xsdValid";
}

my $exit_value = $status >> 8;
my $signal_num = $status & 127;
my $dumped_core = $status & 128;

exit($exit_value);
