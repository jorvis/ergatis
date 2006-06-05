#!/usr/local/bin/perl

=head1	NAME

create_list_file.pl - script to create a list of files given a certain regex

=head1	SYNOPSIS

USAGE: create_list_file.pl
       [  --directory|d=/path/to/some/directory/
          --output_list|o=/path/output.list
          --regex|r=string regular expression
          --include_gz|i ]

=head1	OPTIONS

B<--direcotry,-d>
	Location of the directory to search.  No wild card characters should be
    used in the directory.
    (Default: current directory)

B<--output_list,-o>
    File that will be the list of found files
    (Default: standard out)

B<--regex,-r>
    Regular expression for entire file name (including full path)
    (Default: finds all files)

B<--include_gz|g>
    This flag tells the script to include .gz extensions in the list

B<--help,-h>
	Prints this help message

=head1	DESCRIPTION

    This script will search a directory for file names containing the given
    regular expression.  By default, all .gz or .gzip extensions are stripped
    off.  Programs using these list files should know to search for .gz|.gzip
    versions.

    **Note: No wild card characters should be used in this option.

=head1	INPUT

    There are no required options to the script and there isn't really any
    input.  See above for defaults.  

=head1	OUTPUT

    The script will output a list of files that match the regular expression
    given by option --regex|-r. If no output file is given, this script will 
    print out to standard out.  

=head1	CONTACT

    Kevin Galens
    kgalens@tigr.org

=cut


use strict;
use warnings;
use File::Find;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Cwd;

#GLOBAL VARIABLES##########################################################
my $directory = getcwd;     #The directory to search
my $gz = 0;                 #Flag indicating the inclusion of .gz extension
my $regex = qr/".*"/;           #The regular expression to search for
my $out;                    #Output file handle
###########################################################################

my %options = ();
my $results = GetOptions (\%options, 
                          'directory|d=s',
                          'include_gz|g',
                          'output_list|o=s',
                          'regex|r=s'
                          );

&check_options;
find(\&process, $directory);


#SUB-ROUTINES##############################################################

sub check_options {

    #If directory option is specified and is also a valid directory,
    #then store it.
    if($options{'directory'}) {
        die("$options{directory} is not a valid directory") 
            unless(-d $options{'directory'});
        $directory = $options{'directory'};
    }

    #If the gz option was specified, set this flag to one.
    $gz = 1 if($options{'include_gz'});

    #regular expression setting.  There is no quality check here.
    $regex = qr/$options{'regex'}/ if($options{'regex'});

    #Set the output option to what the user specified or stdout.
    if($options{'output_list'}) {
        open($out, ">$options{'output_list'}") 
            || die ("Unable to open $options{output_list}");
    } else {
        open($out, "> -") || die ("unable to open STDOUT ($!)");
    }
}

#process
#gets called when File::Find finds a file.  This will print out the file
#name to the $out file handle without the gz (or with .gz if the flag is on).
sub process {
    my $filename = $File::Find::name;
    print "$_\n" if(-f);
    (-f && $filename =~ /^$regex(\.gz)?$/) or return;
    $filename =~ s/\.(gz|gzip)$// unless($gz);
    print $out $filename."\n";
}
#EOF###################################################
