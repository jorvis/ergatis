#!/usr/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

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
use Ergatis::Logger;
}

umask(0000);

my %options = ();

my ($directory, $extension, $log4perl, $help, $debug, $noAssertFeature);

my $results = GetOptions (\%options, 
                          'directory|D=s' => \$directory, 
                          'extension|e=s' => \$extension,
                          'log4perl|l=s'  => \$log4perl, 
			  'debug|d=s'     => \$debug,
			  'no_assert_feature=s' => \$noAssertFeature,
                          'help|h'        => \$help ) || pod2usage();


my $logger = &getLogger();

#
# 1. Find all .append files in directory
# 2. For each .append file in directory, check for existence of corresponding file
# 3. Append .append file to file
#


&check_directory($directory);

my ($allbcpfiles, $allappendfiles) = &get_all_bcp_files($directory);

&append_files($allbcpfiles, $allappendfiles, $extension);

print "$0 execution completed\n";
print "The log file is '$log4perl'\n";
exit(0);

##--------------------------------------------------------------------
##
##              END OF MAIN -- SUBROUTINES FOLLOW
##
##--------------------------------------------------------------------
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

    if ($bcpfile =~ /feature\.out$/){
	## Check whether the feature.out.append file has no content.
	if (-z $appendfile){
	    if (!$noAssertFeature){
		$logger->logdie("appendfile '$appendfile' has no content");
	    } else {
		$logger->warn("appendfile '$appendfile' has no content");
	    }
	}
    }


    eval {
	print `cat $appendfile >> $bcpfile`;
    };
    if ($@){
	die "Caught exception while attempting to cat $appendfile >> $bcpfile: $!";
    }

    rename($appendfile, "$appendfile.$$.bak") || die "Could not mv $appendfile $appendfile.$$.bak: $!";

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

sub getLogger {


    if (!defined($options{'log4perl'})){
	$log4perl = $directory . '/' . File::Basename::basename($0) . '.' . $$ . '.log';
	print STDERR "--log4perl was not specified and therefore was set to '$log4perl'\n";
    }

    my $logger = new Ergatis::Logger('LOG_FILE'=>$log4perl,
				     'LOG_LEVEL'=>$options{'debug'});

    $logger = Ergatis::Logger::get_logger();

    return $logger;
}
