#! /usr/local/bin/perl

=head1 NAME

move_niaa_blastp.pl - one-sentence description.

=head1 SYNOPSIS

    USAGE: move_niaa_blastp.pl -D $database -L $list_dir [-b $btab_list -r $raw_list] [--help]

=head1 OPTIONS

=over 4

=item -D, -d, --Database

Annotation database

=item -L, -l  --list_path *

path to the list of results (assuming standard naming for the lists)

=item -b, --btab_list *

path & name of the list of btab files

=item -r, --raw_list *

path & name of the list of raw output files

=item -h, -m, -?, --help, --man

Prints this page and quits

=back

* You need to specify either -L, -l, -list_path or BOTH -b, --bsml_list & -r --raw_list


=head1  DESCRIPTION

This is a wrapper script for Move_Workflow_Transcr_Level_Searches.pl, dedicated to moving the results of the blast searches against Panda (AllGroup.niaa) having simply to specify the name of the database.
In the case the names of the lists of outputs follow Workflow naming conventions, it is necessary only to specify (with the option -L, -l --list_path) the name of the directory where those lists are located.  
Otherwise the names of both the btab and the raw output lists must be specified using respectively the options -b, --btab_list and -r, --raw_list.
The script is able to accept even relative path, transforming them in absolute paths.  If a name of a list is given instead of the path to the lists is given as argument of -L, -l --list_path, the program is able to discard the file name, retaining the path information.
This script assumes that the wrapped script 'Move_Workflow_Transcr_Level_Searches.pl' is located in the very same directory where this very script is located.

In the case one specifies both btab and raw output lists names (the least likely case), the script assumes that you do it correctly and it doesn't currently check for any mixing up between the names.

=head1  INPUT

This script takes as input the name of the database where the files need to be moved to and either the path to the lists of results or the names of both the lists of results (btab & raw output).

=head1  OUTPUT

The only output provided by this script is constituted by a series of messages at STDERR, describing the action that is going to be launched next.

=head1  CONTACT

    Paolo Amedeo
    pamedeo@tigr.org

=begin comment

## legal values for status are active, inactive, hidden, unstable
    status: active
    Comment
=end comment

=cut


use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Spec;

MAIN:{
	my $prN = basename($0);

    my ($db, $list_path, $btab_list, $raw_list, $help);
    
    my $defa_btab_list = 'wu-blastp.btab.list';
    my $defa_raw_list  = 'wu-blastp.raw.list';
    my $target_dir     = 'blastp';
    my $btab_end       = 'nr.btab';
    my $raw_end        = 'nr';
    my $moving_script  = (File::Spec->splitpath(File::Spec->rel2abs($0), 0))[1] . '/Move_Workflow_Transcr_Level_Searches.pl';

    Getopt::Long::Configure("no_ignore_case");
    
    GetOptions('database|D|d=s'     => \$db,
               'list_path|L|l=s'    => \$list_path,
               'btab_list|b=s'      => \$btab_list,
               'raw_list|r=s'       => \$raw_list,
               'help|man|h|m|?'     => \$help) || die "\n\nProblems processing the options\n\n";
    
    pod2usage(-exitval => 0,
              -verbose => 1,
              -output  => \*STDERR) if $help;
    
    pod2usage(-msg     => "\n\nOption -D, -d, --Database (Annotation database) is required\n\n",
              -exitval => 1,
              -verbose => 1,
              -output  => \*STDERR) unless defined $db;
              
    $db = lc($db);

    if (defined $list_path){
        if (defined $btab_list || defined $raw_list){
            pod2usage(-msg     => "\n\nOption -L, -l --list_path is incompatible with -b, --btab_list and -r, --raw_list\n\n",
                      -exitval => 1,
                      -verbose => 1,
                      -output  => \*STDERR);
        }
        my $no_file = -d $list_path;
        $list_path = (File::Spec->splitpath(File::Spec->rel2abs($list_path),  $no_file))[1];
        $btab_list = "$list_path/$defa_btab_list";
        $raw_list  = "$list_path/$defa_raw_list";
    }
    elsif (! defined $btab_list || ! defined $raw_list){
        pod2usage(-msg     => "\n\nOptions -b, --btab_list and -r, --raw_list must be both specified\n\n",
                  -exitval => 1,
                  -verbose => 1,
                  -output  => \*STDERR);   
    }
    die "\n\nImpossible to find the btab list \"$btab_list\"\n\n" unless -f $btab_list;
    die "\n\nImpossible to find the raw list \"$raw_list\"\n\n"   unless -f $raw_list;
    
    open(my $btab, $btab_list) || die "\n\nImpossible to open the file \"$btab_list\" for reading\n\n";
    open(my $raw,  $raw_list)  || die "\n\nImpossible to open the file \"$raw_list\" for reading\n\n";
    
    # Idiot-proofing, checking that the files for matching database name...
    
    foreach my $list ($btab, $raw){
        while (<$list>){
            next unless /\w/;
            chomp();
            
            die "\n\nUnrecognized file name: \"$_\" (wrong database?)\n\n" unless /$db\.model\.\d+_\d+(?:\.\d+)?\.wu-blastp\.[\w.]+/o;
        }
        close($list);
    }
    
    print STDERR "\nMoving the btab files...\n";
    system("$moving_script -D $db -T $target_dir -L $btab_list -t $btab_end")  && die "\nProblems moving btab files: \"$!\"\n\n";
    
    print STDERR "\nMoving the raw output files...\n";
    system("$moving_script -D $db -T $target_dir -L $raw_list -t $raw_end -z")    && die "\nProblems moving raw output files: \"$!\"\n\n";
    
    print STDERR "\nDone\n\n";
}

