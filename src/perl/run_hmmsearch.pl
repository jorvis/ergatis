#!/usr/local/packages/perl-5.8.5/bin/perl

eval 'exec /usr/local/packages/perl-5.8.5/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1	NAME

    run_hmmsearch.pl - runs hmmsearch using multiple hmms against one seq_file

=head1	SYNOPSIS

USAGE: bsml2legacydb.pl
        --seq_file|i=/path/to/some/file.fasta
        --hmm_list|l=/path/to/hmm.list
        --hmm_file|f=/path/to/file.hmm
        --output_file|o=/path/output.raw
	[
		--other_opts|r='string of hmmsearch options'
        --help|-h
        
	]

=head1	OPTIONS

B<--seq_file,-i>
    Location of the input sequence file

B<--hmm_list,-l>
    List of hmm files.  Specific formats taken listed below. 

B<--hmm_file,-f>
    The hmm file (if only one is being used).

B<--hmm_dir,-d>
    A directory of hmm files.

B<--hmm_dir_ext,-e>
    The extension of hmm files in the directory to be used to search.  Will only be used
    when --hmm_dir option is provided.  If left blank, the program will try to use all files
    in that directory.

B<--output_file,-o>
    The file to output the results of all the hmm's in list (or single hmm if hmm_file 
    option is given) against the given seq_file.

B<--other_opts,-r>
    Other options which will be passed directly into hmmsearch.

B<--help,-h>
    Prints this message.

=head1	DESCRIPTION

    This script is a wrapper to hmmsearch, allowing multiple hmms to be passed and run against
    one sequence file.  Currently, the list of hmms does not need to include the full path to
    the hmms (although it is much safer if it is included).  If the full path is not included,
    they will be searched for in the same directory that the list resides.  Also, if both of 
    these options do not exist, then the extensions _fwd.HMM and _rev.HMM are added to the end
    of the hmm names in the list and looked for in the directory where the list resides.  

    For example:

    %cat /path/list.file
         
        RF00001
        RF00002
        ...

    %run_hmmsearch.pl -i some_file.fasta -l /path/list.file -o output.raw 
    
    will run these hmmsearch commands (assuming /path/RF00001_fwd.HMM exists)
      
    'hmmsearch /path/RF00001_fwd.HMM some_file.fasta'
        

=head1	INPUT

    The input to this script includes a sequence fasta file (can contain one or multiple fasta
    sequences) and an hmm file or list file.  For requirements of list file, see Description section
    above.

    The hmm files can be input to this program in three ways.  As a list of files (hmm_list), a single file
    (hmm_file) or a directory of files (--hmm_dir) with a specified extensions (--hmm_dir_ext).

=head1	OUTPUT

    The output of this program involves an output file containing all the hmmsearch results from each
    hmm against each sequence in the fasta file.  Each time a new hmm is used, the string '//' to separate
    different outputs.

=head1	CONTACT

    Kevin Galens
    kgalens@tigr.org

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case pass_through);
use Pod::Usage;
use Ergatis;:Logger;
use IPC::Open3;


########### GLOBALS and CONSTANTS####################
use constant HMMSEARCH => '/usr/local/bin/hmmsearch';    ##Holds path to hmmsearch executable
use constant HMMEXT => ['_fwd.HMM','_rev.HMM'];          ##See the find_hmm_file subroutine
use constant ENDHEAD => 'HMMER 2.3';                     ##Holds line to search for end of hmmsearch
                                                            ##output header (error checking)
my $seq_file;                                            ##File name of the seq_file
my @hmms;                                                ##Holds the hmm names
my $output_file_handle;                                  ##File handle for output file
my $other_opts;                                          ##Any other hmmsearch options you may want
######################################################


### Retrieve command line options
my %options = ();
my $results = GetOptions (\%options, 
                          'seq_file|i=s',
                          'hmm_list|l=s',
                          'hmm_file|f=s',
                          'hmm_dir|d=s',
                          'hmm_dir_ext|e=s',
                          'output_file|o=s',
                          'other_opts|r=s',
                          'help|h'
                          ) || pod2usage();


### Setup some logging
my $logfile = $options{'log'} || Ergatis;:Logger::get_default_logfilename();
my $logger = new Ergatis;:Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();


### Make sure all passed options are valid and set global variables
&parse_parameters(\%options);


### We need to run hmmsearch on each of the hmm files against the seq file and 
### print them all to one file.
foreach my $hmm_file(@hmms) {
    chomp $hmm_file;

    ## If the file exists just run hmmsearch
    unless(-e $hmm_file) {
        ## Will return files matching a certain pattern to use as HMMs
        my $tmp_arr_ref = &find_hmm_files($hmm_file);
        $logger->logdie("Could not find a match for hmm $hmm_file") 
            if($tmp_arr_ref == 0);
        push(@hmms,@{$tmp_arr_ref});
        next;
    }

    ## Running hmmsearch
    run_hmmsearch($hmm_file);
}

exit(0);

##################### SUB-ROUTINES ####################################################

#Name:         parse_parameters
#Description:  Checks the command line options for required and invalid options.  Also sets
#              some global variables with values from these options.
sub parse_parameters {
    # Display documentation if the user needs help.
    if( $options{'help'} ){
        pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
    }

    # Make sure at least an hmm_list or an hmm_file name was passed and set @hmms accordingly.
    if($options{hmm_list}) {
        open(IN, "< $options{hmm_list}") || 
            $logger->logdie("Unable to open $options{hmm_list} ($!)");
        @hmms = <IN>;
        close(IN);
    } elsif($options{hmm_file}) {
        $logger->logdie("Value passed in for --hmm_file|-f option ".
                        "($options{hmm_file}) does not exist")
            unless(-e $options{hmm_file});
        push(@hmms, $options{hmm_file});
    } elsif($options{hmm_dir} ) {
        $logger->logdie("Value passed into hmm_dir option ($options{hmm_dir}) must be a valid directory")
            unless(-d $options{hmm_dir});
        my $tmp_ext = "";
        $tmp_ext = $options{hmm_dir_ext} if($options{hmm_dir_ext});
        opendir(DIR, $options{hmm_dir} ) || $logger->logdie("can't opendir $options{hmm_dir}: $!");
        @hmms = grep { -f "$options{hmm_dir}/$_" && /.*$tmp_ext$/ } readdir(DIR);
        closedir DIR;
        $logger->logdie("No files of extensions $tmp_ext were found in the directory $options{hmm_dir}")
            unless(@hmms > 0);
    }else {
        $logger->logdie("One of hmm_list, hmm_file, or hmm_dir options must be passed");
    }

    # Option seq_file must be provided.
    unless($options{seq_file} && -e $options{seq_file}) {
        $logger->logdie("Option seq_file was not passed or does not exist");
    } else {
        $seq_file = $options{seq_file};
    }

    # Option output_file must be provided.
    if($options{output_file}) {
        open($output_file_handle, "> $options{output_file}") ||
            $logger->logdie("Unable to open $options{output_file} for writing ($!)");
    } else {
        $logger->logdie("Option output_file must be passed.");
    }

}


#Name:        run_hmmsearch
#Arguments:   name of an hmm_file
#Description: Runs the hmmsearch command with given hmm_file against seq_file (global variable).
sub run_hmmsearch {
    my $hmm_file = shift;
    my $system_call = HMMSEARCH." $hmm_file $seq_file $other_opts";

    open3(undef,\*OUT,\*ERR,$system_call) ||  #Perl function of the day 2006.21.06
        $logger->logdie("Unable to run $system_call"); 

    my @stdoutFromHmmsearch = <OUT>;
    close(OUT);
    
    my @stderr = <ERR>;
    close(ERR);
    
    $logger->logdie("\n\nhmmsearch call ($system_call) failed.  Message: \"@stderr\"") 
        if(@stderr > 0);

    print $output_file_handle @stdoutFromHmmsearch;
    print $output_file_handle "//\n";

}

#Name:         find_hmm_files
#Agruments:    hmm_file name
#Description:  Searches for an hmm file if the file does not exist
#              (i.e. the full path is not given) it will search for it
#              in the same directory as the list file (if that option is given) and also will search
#              for the same file name with the extensions supplied in the HMMEXT constant.  This
#              was implemented to work with the pre-existing hmm list files.
sub find_hmm_files {
    my $hmm_file = shift;
    my $retval = 0;
    my $retval_array_ref = [];
    my $dir = "";

    if($options{hmm_list} && $options{hmm_list} =~ m|^(.*)/[^/]+|) {
        $dir = $1;
        push(@{$retval}, "$1/$hmm_file") if(-e "$1/$hmm_file");
    }
    my $hmm_extension_ref = HMMEXT;

    foreach my $hmm_ext(@{$hmm_extension_ref}) {
        push(@{$retval_array_ref}, "$hmm_file$hmm_ext") if(-e "$hmm_file$hmm_ext");

        if($dir eq "") {
            $logger->logdie("Could not parse the directory from the list file");
        }
        push(@{$retval_array_ref}, "$dir/$hmm_file$hmm_ext") if(-e "$dir/$hmm_file$hmm_ext" && $dir ne "");
    }
    
    push(@{$retval_array_ref}, [$hmm_file]) if(-e $hmm_file);

    $retval = $retval_array_ref if(@{$retval_array_ref} > 0);
    return $retval;
    
}

## EOF ##############################################################
