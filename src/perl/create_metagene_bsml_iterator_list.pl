#!/usr/bin/perl

use warnings;
use strict;
use File::Basename;
use Pod::Usage;		
use Ergatis::Logger;	
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
my $logger;

my %options = &parse_options();
my $bsml_file_list = $options{'bsml_list'};
my $output_iter = $options{'output'};

## Parse through BSML file list
my @bsml_files = &parse_bsml_file_list($bsml_file_list);

open (OUTFILE, "> $output_iter") or $logger->logdie("Could not open output iterator $output_iter: $!");
print OUTFILE '$;BSML_FILE_BASE$;' . "\t" .
              '$;BSML_FILE_NAME$;' . "\t" .
              '$;BSML_FILE_PATH$;' . "\n";
              
foreach my $file (@bsml_files) {
    my $filename = basename($file);
    my $filebase = fileparse($file, '\.[^\.]*');
    print OUTFILE "$filebase\t$filename\t$file\n";
}              

###############################################################################
#####                          SUBROUTINES                                #####
###############################################################################

#---------------------------------------------------------------------------
# parse through bsml file list, verifying each file exists and is readable
#---------------------------------------------------------------------------
sub parse_bsml_file_list {
    my $file_list = shift;
    my @bsml_files = ();
    
    open(BSML, $file_list) or $logger->logdie("Could not open BSML file list $file_list: $!");
    while (my $file = <BSML>) {
        chomp ($file);
        push (@bsml_files, $file) if ( &verify_file($file) );
    }
    
    return @bsml_files;
}

#---------------------------------------------------------------------
# verifies that a file exists, is readable and is not zero-content.
#---------------------------------------------------------------------
sub verify_file {
    my @files = @_; 
    
    foreach my $file (@files) {
        next if ( (-e $file) && (-r $file) && (-s $file) );
    
        if      (!-e $file) { $logger->logdie("File $file does not exist")   }   
        elsif   (!-r $file) { $logger->logdie("File $file is not readable")  }
        elsif   (!-s $file) { $logger->logdie("File $file has zero content") }
    }
        
    return 1;
}

#--------------------------------------------
# parse command line options
#--------------------------------------------
sub parse_options {
    my %opts = ();
    GetOptions(\%opts, 
                'bsml_list|i=s',
                'output|o=s',
                'log|l=s',
                'debug|d=s',
                'help=s' ) || pod2usage();
                
    ## Intialize logger
    my $logfile = $opts{'log'} || Ergatis::Logger::get_default_logfilename();
    my $debug = $opts{'debug'} ||= 4;
    $logger = new Ergatis::Logger( 'LOG_FILE'      => $logfile,
                                   'LOG_LEVEL'     => $debug );
    $logger = Ergatis::Logger::get_logger();                   
    
    ## Make sure that a list of BSML files was passed in as well as a destination output file..
    defined ($opts{'bsml_list'}) || $logger->logdie("Please specify an input list of BSML files.");
    defined ($opts{'output'}) || $logger->logdie("Please specify an output iterator file.");
    
    return %opts;                                 
}
