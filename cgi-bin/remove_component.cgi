#!/usr/local/bin/perl

use strict;

use XML::Twig;

use CGI qw(:standard);

use Tree::DAG_Node;
use Config::IniFiles;

my $conffile = param('conffile');

unlink($conffile);

print header();
print "<html><body>Removed $conffile</body></html>";

exit;
