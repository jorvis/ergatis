#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

=head1 NAME

merge_lists.pl - Iterate through a directory of lists and merge them into a single list.

=head1 SYNOPSIS

    USAGE: merge_lists.pl 
                --input_dir=/some/dir/
                --output_list=/some/path/something.list
              [ --glob='*.list' ]

=head1 OPTIONS

B<--input_dir,-i>
    The input directory containing the list files.

B<--output_list,-o>
    The file to which the output list will be written.

B<--glob,-g>
    Searches for file names with this pattern within --input_dir

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

Because distributed processes cannot all write to the same list file, each write their
own and we use this script to merge them into a single file, after the iteration steps
are finished.  Of course, this can also be used outside of workflow if you simply have
files that you want to merge.

Why not use just cat with a redirector?  We could, but this allows us to log each file
as it is merged into the final version as well as do any needed sanity checking.

=head1  INPUT

The required input is simply a directory path, specified by --input_dir, containing the files 
you want to merge.  This will merge every file in the directory into the output
file unless you also pass a filter with the --glob option.  Glob allows one to specify
many of the same patterns for file identification as using grep on the command line, such
as * ? [ ] .

=head1  OUTPUT

Each of the files will be concatenated, with no content modification, onto the end of the
output file, specified with the --output_list option.  If the output file already exists
it will be overwritten.

=head1  CONTACT

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
BEGIN {
use Ergatis::Logger;
}

my %options = ();
my $results = GetOptions (\%options, 
                          'input_dir|i=s',
                          'glob|g=s',
                          'output_list|o=s',
                          'log|l=s',
                          'debug|d=s',
                          'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

## open the output file
my $outputfn = $options{output_list};
open(my $ofh, ">$outputfn") || $logger->logdie("can't create output file $outputfn");

## merge each file.
use File::Find;
find sub {
   my $file = $File::Find::name;
   if($file =~ /$options{glob}/){
	$logger->info("merging $file into $outputfn") if ($logger->is_info);
    
        ## open this file, and write it to the output file
        open(my $ifh, "<$file") || $logger->logdie("can't read input file $file");
        while (<$ifh>) { print $ofh $_ }
   }
}, "$options{input_dir}";

exit;

sub check_parameters {
    my $options = shift;
    
    ## input_dir and output_list are required
    unless (defined $options{input_dir} && defined $options{output_list}) {
        pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR}); 
    } 
    
    ## make sure input_dir exists
    if (! -e "$options{input_dir}") {
        $logger->logdie("the input dir passed ($options{input_dir}) cannot be read or does not exist");
    }
    
    ## handle default glob
    if (! defined $options{glob}) { $options{glob} = '*' }
}
