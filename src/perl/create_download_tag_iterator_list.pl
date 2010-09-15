#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Basename;
use Ergatis::Logger;
umask(0000);

#-------------------------------------------------
# GLOBALS/DEFAULTS
#-------------------------------------------------
my $logger;

my %options = &parse_options();
my $tag_list = $options{'tag_list'};
my $output_iter = $options{'output'};

my $tags_to_download = &parse_tag_list($tag_list);
print_download_tag_iterator($tags_to_download, $output_iter);

#############################################################
####                  SUBROUTINES                        ####
#############################################################

#--------------------------------------------------
# print output iterator list that will be used 
# by ergatis
#--------------------------------------------------
sub print_download_tag_iterator {
	my ($tags, $output_iter) = @_;
	
	open (ITEROUT, "> $output_iter") or $logger->logdie("Could not open output iterator $output_iter for writing: $!");
	print ITEROUT "\$;I_FILE_BASE\$;\n";    ## Print the header identification line
	foreach my $tag (@$tags) {
		print ITEROUT "$tag\n";
	}
	
	close ITEROUT;
}

#--------------------------------------------------
# parse tag list input file and return all tag
# names to be downloaded.
#--------------------------------------------------
sub parse_tag_list {
	my $tag_list = shift;
	my $ret_tags = ();
	
	open (TAGS, $tag_list) or $logger->logdie("Could not open tag-list input file: $!");
	while (my $line = <TAGS>) {
		next if ($line =~ /^#/);
		chomp ($line);
		
		my ($tag_name, $files) = split(/\t/, $line);
		push (@$ret_tags, $tag_name);
	}
	
	close TAGS;
	return $ret_tags;
}


#--------------------------------------------------
# parse command line arguments
#--------------------------------------------------
sub parse_options {
	my %opts = ();
	GetOptions(\%opts,
	           'tag_list|i=s',
	           'output|o=s',
	           'log|l:s',
	           'debug|d:s',
	           'help') || pod2usage();
	           
    if ( $opts{'help'} ) { pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ); }
    ## Set logger
    my $logfile = $opts{'log'} || Ergatis::Logger::get_default_logfilename();
    $logger = new Ergatis::Logger( 'LOG_FILE'   =>  $opts{'log'},
                                   'LOG_LEVEL'  =>  $opts{'debug'} );
    $logger = Ergatis::Logger::get_logger();
    
    ## We need to make sure our tag_list parameter and output parameter are defined
    defined ($opts{'tag_list'}) || $logger->logdie("Please provide the path to a valid tag list file.");
    defined ($opts{'output'}) || $logger->logdie("Please provide a proper output file iterator path");   	                   
    
    return %opts;
}
