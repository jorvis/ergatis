#!/usr/local/bin/perl

## free love!
umask(0000);

=head1  NAME

ergatis_test_init - Automated ergatis testing using stored pipeline templates.

=head1 SYNOPSIS

USAGE: ergatis_test_init [--install] [--no-run] --repository-root=/repository/directory/to/run/tests stored_pipeline1 [stored_pipeline2 ...]

=head1 OPTIONS

B<--install,-i>
    Do a complete install of ergatis before running the test pipelines.
    (If omitted will use existing ergatis install.)

B<--no-run,-n>
    Instantiate the test pipeline(s) but do not execute them.
    (For verification purposes to prevent running malformed pipelines.)
   
B<--ignore-lock>
	Disable the lock file system that prevents running multiple instances of this script from the same directory.
	ONLY use this flag if you can be sure that your test pipelines are configured with an ergatis install path, 
    repository root, and database where there is no risk of collision with other testing pipelines.
 
B<--repository-root,-r>
    Repository root in which to create and run the test pipeline instance(s).

B<--workflow-path,-w>
	Path to workflow executables (optional)	

B<--man,-m>
    Display the pod2usage page for this utility.

B<--help,-h>
    Print this help.

=head1   DESCRIPTION

This script is used to do something really magical.

=head1 INPUT

Takes a file of some kind.

=head1 OUTPUT

Outputs another file.

=head1 CONTACT

Aaron Gussman
agussman@tigr.org
Brett Whitty
bwhitty@tigr.org

=head2 Known Bugs

1) All non-switch arguments are interpreted as pipeline names to be run.
Invalid pipelines are run through all steps, including the db query verification steps.
This should be changed so that it chokes on the invalid pipeline name when it finds that
the XML doesn't exist.

=cut
use strict;
use warnings;
use Pod::Usage;
use Cwd qw(getcwd realpath);
use Sys::Hostname;
use File::Basename;
use File::Path;
use File::Copy;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Config::IniFiles;
use Digest::MD5 qw(md5 md5_hex md5_base64);
#use lib '/usr/local/packages/perl-5.8.5/lib/site_perl/5.8.5/'; ## PERL5LIB is set in crontab

my $cwd = getcwd(); ## remember cwd
my $host = hostname();

my %options=();
GetOptions(\%options, 
#	"ergatis-target|t=s", #read from sharedconf.ini
	"install|i",
	"repository-root|r=s", 
	"workflow-path|w=s",
	"ignore-lock",
#	"branch-ergatis=s",  ##branches are now specified in the installer testing_automated.ini file
#	"branch-chado=s", 
#	"branch-prism=s",
#	"branch-bsml=s",
	"man|m",
	"help|h",
) || pod2usage();

if( $options{'man'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 1, -output => \*STDOUT} );
}

#my $repository_root;
unless (defined($options{'repository-root'})) {
	print "USAGE ERROR: Must provide a repository root!\n";
   	pod2usage( {-exitval=>1, -verbose => 1, -output => \*STDOUT} );
}
unless (-d $options{'repository-root'}) {	
   	print "Repository root '".$options{'repository-root'}."' invalid:\n$!\n";
 	pod2usage();
}
$options{'repository-root'} =~ s/\/$//;

if (!defined($options{'workflow-path'})) {
	$options{'workflow-path'} = '/usr/local/devel/ANNOTATION/workflow';
}

#my $install_root = $options{'ergatis-target'};

## read the install_root directory from the sharedconf.ini file 
## in the repository root directory we have been provided
my $cfg = Config::IniFiles->new( -file => $options{'repository-root'}."/workflow_config_files/sharedconf.ini")
	|| print_usage("Repository root (".$options{'repository-root'}.") did not contain a valid workflow_config_files/sharedconf.ini");
	  
my $install_root = $cfg->val('init', '$;BIN_DIR$;') 
    || print_usage($options{'repository-root'}."/workflow_config_files/sharedconf.ini does not contain a valid BIN_DIR");
$install_root =~ s/\/$//; ##strip off trailing slash if it exists
$install_root =~ s/\/bin$//; 

my $install_hash = md5_base64($install_root);

my $tmp_path = "/tmp/".$install_hash."_".allNow()."_".$$;

## make the temp checkout path
mkpath($tmp_path, 1, 0777);

##get absolute path of directory containing this script
my $exec_path = dirname(realpath($0));

##directory for logs
my $log_dir = $exec_path."/logs/".allNow();

## checkout installer script if it hasn't been done already
unless (-e $exec_path."/bin/ergatis_installer.pl") {
	checkout_installer();
}

print "Executing workflow test\n";
print "  logs: $log_dir\n";
print "  repository root: ".$options{'repository-root'}."\n";
print ( $options{'install'} ? "  wipe: yes\n  install: ".$install_root ."\n" : "  wipe: no\n  install: existing ($install_root)\n" );
print "  pipelines: \n";
foreach (@ARGV) {
    print "    - $_\n";
}
if (!$options{'ignore-lock'}) {
if (-e "$exec_path/.lock_$install_hash") {
	##a lock file exists in the testing directory
	print "Lock file exists!!!\n";
	print "  Testing pipeline already running, or previous execution died poorly\n";
	print "  Manual intervention required\n";
	exit(1);
} else {
	##write a lock file
	open (LOCK, ">$exec_path/.lock_$install_hash") || die "Couldn't write lock file!!!";
	print LOCK allNow();
	close LOCK;
}
}

mkpathORdie($log_dir); #should make have no echo??

if ($options{'install'}) {

	open(STDOUT, ">$log_dir/ergatis_install.stdout")|| die "can't create $log_dir/ergatis_install.stdout: $!";
    open(STDERR, ">$log_dir/ergatis_install.stderr")|| die "can't create $log_dir/ergatis_install.stderr: $!";      

	my $dir_name = basename($install_root);
	unless (-e $install_root."/../$dir_name.ini") {
		die $dir_name.".ini must exist in ".realpath($install_root."/../");
	}
	
	doORdie("$exec_path/bin/ergatis_installer.pl -i $install_root -w $tmp_path -U sgc --init --controlfile $install_root/../$dir_name.ini");
	#print("$exec_path/bin/ergatis_installer.pl -i $install_root -w $tmp_path -U sgc --init --controlfile $install_root/../$dir_name.ini\n");
	
	print STDOUT "Copying testing-specific docs to ergatis install...\n";	
	foreach my $file (glob("$exec_path/docs/*")) {
		copy($file, "$install_root/docs");
	}	
    close(STDOUT);
    close(STDERR);

	#### TEMP FIX FOR UMASK PROBLEM
	## change to install root
	chdir($install_root);
	## do a recursive chmod
	rchmod();
	## change back to cwd
	chdir($cwd);
	####
	
}

#run each saved pipeline template passed as an argument
foreach my $pipe (@ARGV) {
	my @badini = ();
	my $exitflag = 0;
    if (!(-d $pipe) && (-d "$exec_path/pipeline_templates/$pipe")) {
		$pipe = "$exec_path/pipeline_templates/$pipe";
    }
    $pipe =~ s/\/$//;
    my $pipe_name = basename($pipe);

    #
    # Run pipeline
    #
	open(STDOUT, ">$log_dir/$pipe_name.pipe.stdout")|| die "can't create $log_dir/$pipe_name.pipe.stdout: $!";
    open(STDERR, ">$log_dir/$pipe_name.pipe.stderr")|| die "can't create $log_dir/$pipe_name.pipe.stdout: $!";

	## $pipe must contain a pipeline.xml file
	unless (-e $pipe."/pipeline.xml") {
		print STDERR "No pipeline.xml in $pipe --- skipping";
		next;
	}

	foreach my $inifile(glob("$pipe/*.ini")) {
		my $exitval = system("$exec_path/bin/ini_variable_check.pl -i $inifile -e $install_root");
		$exitval = $exitval / 256;
		if ($exitval != 0) {
			#print STDERR "Invalid ini file structure in '".basename($inifile)."'.\n";
			print STDERR "Error found validating ini file structure of '$inifile'.\n";
			push(@badini, basename($inifile));
			if ($exitval > 1) {
				print STDERR "This error is a critical error.\n";
				$exitflag = 1;
			} else {
				print STDERR "Error should not interfere with execution of pipeline.\n"
			}
		}
	}

	## if any of the ini files had invalid structure, skip the run
	if ($exitflag) {
		print STDERR "Skipping execution of pipeline '$pipe_name' due to invalid structure in:\n";
		print STDERR join(",", @badini)."\n";
		next;
	}
	
	#clean out stuff in BSML_ and FASTA_ directories 
    #(but not the previous pipelines, because we want to look at them in the future)
	doORdie("find  ".$options{'repository-root'}."/BSML_repository/* -type f -exec rm {} \\; ");
	doORdie("find  ".$options{'repository-root'}."/FASTA_repository/* -type f -exec rm {} \\; ");

    print "\nInstantiating $pipe_name from $pipe\n";
    if (defined($options{'no-run'})) {
    	doORdie("perl $exec_path/bin/instantiate_pipeline.pl --ergatis-install=$install_root --workflow-path=".$options{'workflow-path'}." --template-dir $pipe --repository-root ".$options{'repository-root'});
    } else {
		doORdie("perl $exec_path/bin/instantiate_pipeline.pl --ergatis-install=$install_root --workflow-path=".$options{'workflow-path'}." --execute --template-dir $pipe --repository-root ".$options{'repository-root'});
    }
    close(STDOUT);
    close(STDERR);
		
	## DATABASE QUERY STEP WILL NOW BE ADDED AS LAST STEP TO A SAVED PIPELINE

}

## remove the temp checkout dir
rmtree($tmp_path, 1);

if (!$options{'ignore-lock'}) {
	##remove lock file
	unlink("$exec_path/.lock_$install_hash");
}

exit(0);

#
#  subroutines follow
#
sub doORdie {
    my $cmd = shift;  
    print "$cmd\n";
    system($cmd);
    if ($? == -1) {
       print STDERR "failed to execute: $!\n";
       exit(1);
    } elsif ($? & 127) {
       printf STDERR "child died with signal %d, %s coredump\n",
           ($? & 127),  ($? & 128) ? 'with' : 'without';
       exit(1);
    }
}

#returns the datetime as a string
sub allNow {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    my $time_string = sprintf("%02d%02d%02d%02d%02d%02d\n", substr($year+1900,-2),$mon+1,$mday,$hour,$min,$sec);
    chomp($time_string);
    return $time_string;
}

sub print_usage {
    my $byebye = shift;
    print $byebye;
    pod2usage( {-exitval=>0, -verbose => 1, -output => \*STDOUT} );
}

sub mkpathORdie {
	my $dir = shift;
	eval {mkpath($dir, 1, 0777)};
	if ($@) {
		print STDERR "Couldn't create $dir: $@";
		exit(1);
	}
}

sub rmtreeORdie {
	my $dir = shift;
	eval {rmtree($dir, 1)};
	if ($@) {
		print STDERR "Couldn't remove $dir: $@";
		exit(1);
	}
}


## recursive chmod
sub rchmod {
	dochmod();
	foreach my $item(glob("*")) {
		if (-d $item) {
			chdir($item);
			rchmod();
			chdir(".."); 
		}
	}
}

## does chmod on everything in the current directory
sub dochmod {
	foreach my $item(glob("*")) {
	if (-x $item || -d $item) {
    	chmod(0777, $item);
	} else {
		chmod(0666, $item);
		}
	}
}

sub checkout_installer {
	chdir($tmp_path);
	my $err = system("cvs co -d installer ergatis/ergatis_installer.pl");
	if ($err) {
		die "Checkout of installer script failed.\n";
	}
	copy($tmp_path."/installer/ergatis_installer.pl", $exec_path."/bin/ergatis_installer.pl") || die "couldn't copy installer to bin dir.\n";
	chmod (0777, $exec_path."/bin/ergatis_installer.pl");
	chdir($exec_path);
}
