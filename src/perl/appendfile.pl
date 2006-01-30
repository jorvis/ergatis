#!/usr/local/bin/perl

=head1  NAME 

appendfile.pl - 

=head1 SYNOPSIS

USAGE:  appendfile

=head1 OPTIONS

=item *

B<--extension,-e> 

B<--directory,-D> 
 
=item *

B<--log4perl,-l> Log file

=item *

B<--debug,-d> Debug level

=item *

B<--help,-h> This help message

=head1   DESCRIPTION

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

use File::Basename;
BEGIN {
use Workflow::Logger;
}

umask(0000);

my %options = ();

my ($directory, $extension, $log4perl, $help, $debug);

my $results = GetOptions (\%options, 
                          'directory|D=s' => \$directory, 
                          'extension|e=s' => \$extension,
                          'log4perl|l=s'  => \$log4perl, 
			  'debug|d=s'     => \$debug,
                          'help|h'        => \$help ) || pod2usage();



my $log4perl = $options{'log4perl'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$log4perl,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Workflow::Logger::get_logger();


#
#
# 1. Find all .append files in directory
# 2. For each .append file in directory, check for existence of corresponding file
# 3. Append .append file to file
#


&check_directory($directory);


my ($allbcpfiles, $allappendfiles) = &get_all_bcp_files($directory);

&append_files($allbcpfiles, $allappendfiles, $extension);


sub append_files {

    my ($allbcpfiles, $allappendfiles, $extension) = @_;

    foreach my $appendfile (sort keys %{$allappendfiles}){

	if ($appendfile =~ /(\S+)\.$extension/){
	    
	    my $bcpfile = $1;
	    if (( exists $allbcpfiles->{$bcpfile}) && (defined($allbcpfiles->{$bcpfile}))){

		&append_file_contents($bcpfile, $appendfile);
	    }
	    else {
		$logger->logdie("bcpfile did not exist for appendfile '$appendfile'");
	    }

	}
	else{
	    $logger->logdie("Could not extract BCP file from appendfile '$appendfile'");
	}
    }

}





sub append_file_contents {

    my ($bcpfile, $appendfile) = @_;

    print `cat $appendfile >> $bcpfile`;
    print `mv $appendfile $appendfile.$$.bak`;

}






sub check_directory {

    my $directory = shift;



    if (!defined($directory)){
	$logger->logdie("directory was not defined");
    }


    if (!-e $directory){
	$logger->logdie("directory '$directory' does not exist");
    }

    if (!-r $directory){
	$logger->logdie("directory '$directory' does not have read permissions");
    }
    if (!-d $directory){
	$logger->logdie("directory '$directory' is not a directory");
    }
}






sub get_all_bcp_files {

    my ($directory) = shift;
    
    opendir(THISDIR, "$directory") or $logger->logdie("Could not open directory '$directory'");

    my @allfiles = readdir THISDIR;


    my $allbcpfiles = {};
    my $allappendfiles = {};


    foreach my $file (sort @allfiles){

	if ($file =~ /(\S+)\.out$/){
	    
	    my $fullpath = $directory . "/" . $file;

	    $allbcpfiles->{$fullpath} = $fullpath;
	}
	elsif ($file =~ /(\S+)\.out.append$/){
	    
	    my $fullpath = $directory . "/" . $file;

	    $allappendfiles->{$fullpath} = $fullpath;
	}
    }

    return ($allbcpfiles, $allappendfiles);

}



