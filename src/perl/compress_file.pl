#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

=head1  NAME 

compress_file.pl - wrapper to gzip an output file

=head1 SYNOPSIS

USAGE: compress_file.pl 
        --file=/path/to/some.file 
      [ --output=/path/to/output
        --compress=1
        --list_file=/path/to/some.list
        --log=/path/to/some.log
        --debug=4
      ]

=head1 OPTIONS

B<--file,-f> 
    The path to any file you wish to compress.

B<--compress,-c>
    Optional.  If 0, the compression will not be performed.

B<--remove_source,-r>
    Optional.  Controls whether source file is removed when compressing (default = 1)

B<--list_file,-s>
    Optional.  Writes a list file of the gzipped file created.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used as a gzip wrapper to compress a file.

=head1 INPUT

The input is any file you wish to compress, usually a text file such
as alignment data.  This file is passed with the --file option.

The --compress option exists so that this script can be incorporated
into workflow.  Since Workflow is not a logic-based system, it cannot
call gzip on output *sometimes*.  So we call this script all the time,
and just pass --compress 1 or 0 depending on whether we actually want
the compression done.

The --list_file options allows us to create a list file containing the
full path to the list file created.  This is useful because a workflow
may create a list file for the raw output, but if we zip it, we want to
overwrite that single-entry list file with the gzip'ed path before it
is merged.

Some checks are performed to handle common problems found handling compression
in production - each handling cases where previous compression attempts
have left partial files.

If neither input or output exists the script will die.  If the input file
exists and output doesn't, normal compression will occur.  If the output file
exists and the input file doesn't, a warning is issued and the script completes
successfully.  

=head1 OUTPUT

The output file will be named the same but with the .gz extension

=head1 CONTACT

    Joshua Orvis
    jorvis@users.sf.net

=cut


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Ergatis::Logger;
use File::Basename;

my %options = ();
my $results = GetOptions (\%options,
			  'file|f=s',
			  'compress|c=s',
                          'list_file|s=s',
			  'output|o=s',
                          'remove_source|r=s',
			  'log|l=s',
			  'debug=s',
			  'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Ergatis::Logger::get_logger();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDERR} );
}

&check_parameters(\%options);

my $return_status = 0;

## we're only going to compress if the user requested it
if ($options{compress}) {
    $logger->debug("preparing to compress file $options{file}") if ($logger->is_debug);
    
    ## do both the input and output already exist?
    my $output = $options{output} || "$options{file}.gz";
    
    ## if the output file already exists but the input file doesn't the operation is 
    #   already completed.  log a warning and continue
    if ( -e $output && ! -e $options{file} ) {
        $logger->warn( "output file already exists, operation appears to be completed" ) if $logger->is_warn;
        
    } else {
    
	if ($options{file} eq $output){
	    $output = $options{output} . '.gz';
	}

        ## if a previous run of this failed in the middle both the input and output files
        #   will already exist.  delete the output file and try again.
        if ( -e $output && -e $options{file} ) {
            unlink($output) || $logger->logdie("failed to clear pre-existing output file: $!");
        }
        
        ## if neither exists Bad Things have happened
        if (! -e $output && ! -e $options{file} ) {
            $logger->logdie("input file doesn't exist and no previous attempts detected.");
        }
        
        my $dirname=dirname($options{file});
	print `mkdir -p $dirname`;
        system("gzip -c -f $options{file} > $output");

        if ( $options{remove_source} ) {
            unlink($options{file}) || die "failed to remove source file after compression: $!";
        }

        ## was it successful?
        $return_status = $? >> 8;
        $logger->debug("return status was $return_status") if ($logger->is_debug);
    }
    
} else {
    if($options{output}){
        my $dirname=dirname($options{output});
        print `mkdir -p $dirname`;
	system("cp $options{file} $options{output}");
    }
    $logger->debug("skipping compression of file $options{file} because --compress=$options{compress}") if ($logger->is_debug);
}

if ($options{list_file}) {
    $logger->debug("creating list file $options{list_file}") if ($logger->is_debug);
    my $dirname=dirname($options{list_file});
    print `mkdir -p $dirname`;
    open(my $lfh, ">$options{list_file}") || $logger->logdie("can't create list file $options{list_file}");
    
    print $lfh "$options{file}";
    print $lfh ".gz" if $options{compress};
    print $lfh "\n";
} else {
    $logger->debug("skipping list file creation because --list_file=$options{list_file}") if ($logger->is_debug);
}

exit($return_status);

sub check_parameters{
    my ($options) = @_;
    
    if ($options{'file'} eq ""){
	    pod2usage({-exitval => 2,  -message => "--file option missing", -verbose => 1, -output => \*STDERR});    
    }
    
    ## handle some defaults
    $options{compress}  = 1 if (! defined $options{compress});
    $options{list_file} = 0 if (! defined $options{list_file});
    $options{remove_source} = 1 if (! defined $options{remove_source});
}
