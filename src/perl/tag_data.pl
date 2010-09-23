#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

=head1 NAME

transform_WWARN_input.pl - Transforms input data to WWARN to a common format.

=head1 SYNOPSIS

./tag_data.pl 
        --map_file=/path/to/map/file
        --repository_root=/path/to/ergatis/repository/root
        --pipeline_id=<ergatis pipeline ID>
       [--log=/path/to/log/file
        --debug=<debug level>
        --help]
        
=head1 PARAMETERS

B<--map_file, -i>
    A tab delimited file containing the files to be tagged as well as the desired tag-name
    for these set of files.

B<--repository_root, -r>
    The repository where output files to be tagged reside. This value is replaced in the mapping
    file when present
    
B<--pipeline_id, -p>
    The pipeline ID whose files should be downloaded. This value is replaced in the mapping file
    when present    
    
B<--log, -l>
    Optional. Log file.
    
B<--debug, -d>
    Optional. Debug level.
    
B<--help>
    Prints this documentaiton.
    
=head1 DESCRIPTION            

A wrapper script that facilitates the execution of the CloVR tagData.py script through ergatis.
This wrapper script is meant to allow for custom tagging of output data from a pipeline.

=head1 INPUT

Input is a tab delimited file containing all files to be tagged and the associated tag-name for
the files. Files can be a single file, directories, a list containing multiple files or a list of comma-separated
files.

#TAG_NAME\tFILE
test_tag\tfile1.txt,file2.txt,files.list,/usr/local/test_data/

=head OUTPUT

A file containing all tag names will be written for use by the accompanying download_tag component

=head1 CONTACT

    Cesar Arze
    carze@som.umaryland.edu

=cut

use strict;
use warnings;
use warnings;
use Pod::Usage;        
use Ergatis::Logger;
use File::OpenFile qw(open_file);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

## Need to setup the PYTHONPATH env variable here...
$ENV{'PYTHONPATH'} = "/opt/vappio-py:/opt/vappio-py:/opt/opt-packages/bioinf-v1r4b1//Denoiser/" .
		     ":/opt/opt-packages/bioinf-v1r4b1//PyNAST/lib/:/opt/opt-packages/bioinf-v1r4b1//qiime/lib/";

#----------------------------------------------------------
# GLOBALS/COMMAND-LINE OPTIONS
#----------------------------------------------------------
my $logger;
my $TAG_DATA_EXEC = "/opt/vappio-py/vappio/cli/tagData.py"; 

my %options = parse_options();
my $input = $options{'map_file'};
my $repo_root = $options{'repository_root'};
my $pipeline_id = $options{'pipeline_id'};

my $files_to_tag = parse_mapping_file($input);

## Now that we have our files to tag we can go ahead with the tagging
foreach my $tag_name (keys %$files_to_tag) {
    my @files = @{ $files_to_tag->{$tag_name} };
    my $cmd = $TAG_DATA_EXEC . " --tag-name " . $tag_name . " --overwrite --recursive " .
                  join(" ", @files);
    run_system_cmd($cmd);                  
}

###############################################################################
#####                          SUBROUTINES                                #####
###############################################################################

#----------------------------------------------------------
# parse input mapping file and return an array 
# containing hash containing all files and the associated
# tag-name.
#----------------------------------------------------------
sub parse_mapping_file {
    my $map_file = shift;
    my $tag_files = ();
    
    my $map_fh = open_file($map_file, "in");
    while (my $line = <$map_fh>) {
        next if ($line =~ /^#/);
        chomp($line);
        
        my @files = ();
        my ($tag_name, $files_list) = split(/\t/, $line);
        
        ## Currently tagData will break if any spaces are in the tag name
        ## so we need to replace spaces with underscores
        $tag_name =~ s/\s+/_/;
        
        my @file_tokens = split(/,/, $files_list);
        foreach my $file_token (@file_tokens) {
            _verify_file($file_token, \@files);            
        }
        
        push ( @{ $tag_files->{$tag_name} }, @files);
    }    
    
    return $tag_files;
}

#----------------------------------------------------------
# verify whether a file exists, is readable and if it is a 
# directory; if the file is a list then the list of files 
# will be iterated over and verified as well
#----------------------------------------------------------
sub _verify_file {
    my ($file, $files_ref) = @_;
    my $ret_file;
    
    ## If this is being run in conjunction with an ergatis pipeline 
    ## we need to account for files having $;REPOSITORY_ROOT$; and
    ## $;PIPELINEID$; in them.
    $file =~ s/\$\;REPO_ROOT\$\;/$repo_root/;
    $file =~ s/\$\;PIPELINE_ID\$\;/$pipeline_id/;

    ## First check if our file is a directory
    if ( -e $file) {
        
        ## We are dealing with either a single file or a list if 
        ## the "file" is not a directory
        unless (-d $file) {
            my $fh = open_file($file, "in");
            chomp( my $line_peek = <$fh> );
            
            ## Now check if our first line of this file exists, if it does most 
            ## likely we have a list of files here
            if (-e $line_peek) {
                seek ($fh, 0, 0);
                chomp( my @files = <$fh> );
                
                foreach my $list_file (@files) {
		    _verify_file($list_file, $files_ref);
                }
            } else {
                push(@$files_ref, $file);
            }
            
            close ($fh); 
        } else {
            push(@$files_ref, $file);
        }
    } else {
        $logger->logdie("File $file does not exist.");
    }
}

#----------------------------------------------------------
# run unix system command and determine whether or not
# the command executed successfully
#----------------------------------------------------------
sub run_system_cmd {
    my $cmd = shift;
   
    system($cmd);
    my $success = $? >> 8;

    unless ($success == 0) {
        $logger->logdie("Command \"$cmd\" failed.");
    }   
}

#----------------------------------------------------------
# parse command-line options
#----------------------------------------------------------
sub parse_options {
    my %opts = ();
    GetOptions(\%opts,
                'map_file|i=s',
                'repository_root|r:s',
                'pipeline_id|p:s',
                'log|l:s',
                'debug|d:s',
                'help') || pod2usage();
                 
    &pod2usage( {-exitval => 1, -verbose => 1, -output => \*STDOUT}) if ($opts{'help'} );        
    
    ## Configure logger
    my $logfile = $opts{'log'} || Ergatis::Logger::get_default_logfilename();
    $logger = new Ergatis::Logger( 'LOG_FILE'   =>  $opts{'log'},
                                   'LOG_LEVEL'  =>  $opts{'debug'} );
    $logger = Ergatis::Logger::get_logger();
    
    ## Make sure our parameter are declared correctly
    defined ($opts{'map_file'}) || $logger->logdie('Please specify a valid input map file');
    defined ($opts{'repository_root'}) || $logger->logdie('Path to the repository root must be specified');
    defined ($opts{'pipeline_id'}) || $logger->logdie('A valid ergatis pipeline ID must be specified');
    
    return %opts;
}
