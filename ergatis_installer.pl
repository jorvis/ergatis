#!/usr/local/bin/perl -w

=head1 NAME

ergatis_installer.pl - Checks out code from CVS repository using specified tags then installs software in specified install directory

=head1 SYNOPSIS

USAGE:  ergatis_installer.pl --installdir=installdir --workingdir=workingdir --username=username [-S server] --cvscode="cvscode" --controlfile=controlfile [-d debug_level] [--init] [-h] [-l log4perl] [-m]

=head1 OPTIONS

=over 8

=item B<--username,-U>
    
    Name of user installing the software

=item B<--installdir,-i>
    
    Directory where software will be installed

=item B<--workingdir,-w>
    
    Optional - Directory where cvs checkouts and other work will take place (Default is current working directory)

=item B<--server,-S>
    
    Optional - For Prism based code - name of Sybase server to store in configuration files (Default is SYBTIGR)

=item B<--cvscode>
    
    Optional - Use this command-line arguement to specify code base and tag e.g. --cvscode="chado_prism=prism-v1r9b1;ergatis=ergatis-v1r2b1"

=item B<--controlfile>
    
    Optional - Configuration file for specifying code base and corresponding tags

=item B<--debug_level,-d>

    Optional: log4perl logging level.  Default threshold is WARN

=item B<--man,-m>

    Display the pod2usage page for this utility


=item B<--init>

    Optional - Specifies that the installdir should be wiped  (Default is to not delete the installdir)

=item B<--chmod>

    Optional - Specifies that the installdir should be changed to world read/write

=item B<--help,-h>

    Print this help

=back

=head1 DESCRIPTION

    eragatis_installer.pl - Checks out code from CVS repository using specified tags then installs software in specified install directory

    Assumptions:
    1. User has appropriate permissions (to execute script, access chado database, write to output directory).
    2. All software has been properly installed, all required libraries are accessible.

    Sample usage:
    ./ergatis_installer.pl --username=sundaram --installdir=/usr/local/devel/ANNOTATION/ergatis_install --workingdir=/tmp --server=SYBIL --log4perl=/tmp/install.log --controlfile=/usr/local/devel/ANNOTATION/ergatis_install/control.ini


=cut

use strict;
use File::Copy;
use File::Basename;
use File::Path;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Cwd qw(realpath);
use Pod::Usage;
use Data::Dumper;
use Log::Log4perl qw(get_logger);


$|=1;

my $tmp_dir;
my %opts;
my $exec_path;


my $datamanagementfile = "/usr/local/devel/ANNOTATION/ard/install_log.txt";
my $htmlfile = "/usr/local/devel/ANNOTATION/ard/install_log.html";

umask(0000);


my ($installdir, $workingdir, $username, $help, $log4perl, $man, $cvscode, $controlfile, $server, $debug_level, $init, $chmod);

my $results = GetOptions (
			  'log4perl|l=s'       => \$log4perl,
			  'debug_level|d=s'    => \$debug_level, 
			  'help|h'             => \$help,
			  'man|m'              => \$man,
			  'installdir|i=s'     => \$installdir,
			  'workingdir|w=s'     => \$workingdir,
			  'username|U=s'       => \$username,
			  'cvscode=s'          => \$cvscode,
			  'controlfile=s'      => \$controlfile,
			  'server|S=s'         => \$server,
			  'init'               => \$init,
			  'chmod|c'            => \$chmod
			  );

&pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if ($man);
&pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT}) if ($help);

print STDERR ("username was not defined\n")     if (!$username); 
print STDERR ("installdir was not defined\n")   if (!$installdir); 


&print_usage if(!$username or !$installdir);

#
# initialize the logger
#
$log4perl = "/tmp/cas_installer.pl.log" if (!defined($log4perl));

my $verbose = 1;

my $screen_threshold = 'WARN';
if ($verbose){
    $screen_threshold = 'INFO';
}

Log::Log4perl->init(
		    \ qq{
#			log4perl.logger                       = WARN, A1, Screen
			log4perl.logger                       = WARN, A1
			log4perl.appender.A1                  = Log::Dispatch::File
			log4perl.appender.A1.filename         = $log4perl
			log4perl.appender.A1.mode             = write
			log4perl.appender.A1.Threshold        = WARN
			log4perl.appender.A1.layout           = Log::Log4perl::Layout::PatternLayout
#			log4perl.appender.A1.layout.ConversionPattern = %d %p> %F{1}:%L %M - %m%n 
			log4perl.appender.A1.layout.ConversionPattern = %p> %L %M - %m%n 
			log4perl.appender.Screen              = Log::Dispatch::Screen
			log4perl.appender.Screen.layout       = Log::Log4perl::Layout::SimpleLayout
                        #log4perl.appender.Screen.layout.ConversionPattern =%d %p> %F{1}:%L %M - %m%n 
			log4perl.appender.Screen.Threshold    = $screen_threshold
			Log::Log4perl::SimpleLayout
		    }
		    );


my $logger = get_logger("INSTALLER");


if (!defined($server)){
	$server = 'SYBTIGR';
	$logger->info("server was set to '$server'");
}


&verify_create_directory($installdir);

&clear_install_dir($installdir, $init);

&create_install_docs($installdir);

if (!defined($workingdir)){
    $workingdir = "/tmp";
    if ($logger->is_debug()){
	$logger->debug("workingdir was set to '$workingdir'");
    }
}

&verify_create_directory($workingdir);

&clear_working_dir($workingdir);


my $cvshash = &determine_cvs_revisions($cvscode, $controlfile);


my $date = `date`;

chomp $date;

&execute_installation($cvshash, $installdir, $workingdir, $server);

&record_installation($installdir, $workingdir, $username, $date, $server, $cvshash, $datamanagementfile, $htmlfile);



#--------------------------------------------------------------------------------------------------------------------------
#
#                                     END OF MAIN    --    SUBROUTINES FOLLOW
#
#--------------------------------------------------------------------------------------------------------------------------


#--------------------------------------------------
# install_ontologies()
#
#--------------------------------------------------
sub install_ontologies {

    my ($workingdir, $installdir, $tag) = @_;

    chdir($workingdir);

    my $execstring = "cvs -Q export -kkv -r $tag ontologies";

    &do_or_die($execstring);

    $execstring = "cp $workingdir/ontologies/*.obo $installdir/docs/obo/.";

    &do_or_die($execstring);
}




#--------------------------------------------------
# install_schema()
#
#--------------------------------------------------
sub install_schema {

    my ($workingdir, $installdir, $tag) = @_;

    chdir($workingdir) || die "Could not chdir $workingdir";

    my $execstring = "cvs -Q export -kkv -d chado_schema_sybase -r $tag ANNOTATION/chado/sybase/ddls";

    &do_or_die($execstring);

    $execstring = "cp $workingdir/chado_schema_sybase/*.ddl $installdir/docs/ddls/sybase/.";

    &do_or_die($execstring);

    $execstring = "cvs -Q export -kkv -d chado_schema_postgresql -r $tag ANNOTATION/chado/postgresql/ddls";

    &do_or_die($execstring);

    $execstring = "cp $workingdir/chado_schema_postgresql/*.ddl $installdir/docs/ddls/postgresql/.";

    &do_or_die($execstring);
}



#--------------------------------------------------
# install_peffect()
#
#--------------------------------------------------
sub install_peffect {

    my ($workingdir, $installdir, $tag) = @_;

    chdir($workingdir);

    my $execstring = "cvs -Q export -kkv -r HEAD -d peffect peffect";

    &do_or_die($execstring);

    chdir("peffect");

    $execstring = "find . -name Makefile -exec perl -i.bak -pe 's|/usr/local/devel/ANNOTATION/cas|$installdir|' \{\} \\;";
    
    &do_or_die($execstring);
    
    $execstring = "make";

    &do_or_die($execstring);

    $execstring = "make install";

    &do_or_die($execstring);


}

#--------------------------------------------------
# install_ergatis_src_c()
#
#--------------------------------------------------
sub install_ergatis_src_c {

    my ($workingdir, $installdir, $tag) = @_;

    chdir($workingdir);

    my $execstring = "cvs -Q export -kkv -r $tag -d ergatis_c ergatis/src/c";

    &do_or_die($execstring);

    chdir("ergatis_c/src/c");

    chdir("SnpClusterer");
    $execstring = "make";
    &do_or_die($execstring);
    $execstring = "cp SnpClusterer $installdir/bin";
    &do_or_die($execstring);

    chdir($workingdir);
    chdir("ergatis_c/src/c");
    chdir("chado_record_uniq");
    $execstring = "make";
    &do_or_die($execstring);
    $execstring = "cp chado_record_uniq $installdir/bin";
    &do_or_die($execstring);
}

#--------------------------------------------------
# install_logcabin()
#
#--------------------------------------------------
sub install_logcabin {

    my ($workingdir, $installdir, $tag) = @_;

    chdir($workingdir);

    my $execstring = "cvs -Q export -kkv -r HEAD -d logcabin logcabin";

    &do_or_die($execstring);

    chdir("logcabin");
    
    &do_or_die($execstring);
    
    $execstring = "perl Makefile.PL PREFIX=$installdir";

    &do_or_die($execstring);

    $execstring = "make";

    &do_or_die($execstring);
    
    $execstring = "make install";

    &do_or_die($execstring);


}



#----------------------------------------------------
# clear_working_dir()
#
#----------------------------------------------------
sub clear_working_dir {
	
	my ($dir) = @_;
	
	
	my $execstring = "rm -rf $workingdir/*_install";
	&do_or_die($execstring);
	$execstring = "rm -rf $workingdir/ergatis_c";
	&do_or_die($execstring);
	$execstring = "rm -rf $workingdir/chado_schema_sybase";
	&do_or_die($execstring);
	$execstring = "rm -rf $workingdir/chado_schema_postgresql";
	&do_or_die($execstring);
	$execstring = "rm -rf $workingdir/peffect";
	&do_or_die($execstring);
	$execstring = "rm -rf $workingdir/ontologies";
	&do_or_die($execstring);
	$execstring = "rm -rf $workingdir/logcabin";
	&do_or_die($execstring);
}
	


#----------------------------------------------------
# clear_install_dir()
#
#----------------------------------------------------
sub clear_install_dir {
	
    my ($dir, $init) = @_;
    
    
    if ((defined($init)) && ($init == 1)){
	
	my $execstring = "rm -rf $dir";
	
	&do_or_die($execstring);
	
	mkdir($dir);
	
    }


    my $docsdir = "$dir/docs";

    if (!-e $docsdir) {
	#
	# Create the docs directory
	#
	mkdir($docsdir);
    }
}



#-----------------------------------------------------
# execute_installation()
#
#-----------------------------------------------------
sub execute_installation {

	my ($cvshash, $installdir, $workingdir, $server) = @_;


	#
	# server in configuration file overrides the server specified as command-line argument
	#
	if ((exists $cvshash->{'server'}) && (defined($cvshash->{'server'})) ) {
	    $server = $cvshash->{'server'};
        }
	
	&install_logcabin($workingdir, $installdir);
	
	foreach my $name (sort keys %{$cvshash} ) {

		my $tag = $cvshash->{$name};

		$logger->logdie("tag was not defined for name '$name'") if (!defined($tag));

		my $installname = $name ."_install";

		if ($name eq 'ontologies'){
		    
		    &install_ontologies($workingdir, $installdir, $tag);
		}
		elsif ($name eq 'chado_schema'){
		    
		    &install_schema($workingdir, $installdir, $tag);
		}
		elsif ($name eq 'peffect'){

		    &install_peffect($workingdir, $installdir, "HEAD");
		}
		elsif( $name ne 'server'){
		    
		    chdir($workingdir);
		    
		    
		    my $execstring = "cvs -Q export -kkv -r $tag $installname";
		    
		    &do_or_die($execstring);
		    
		    chdir($installname);
		    
		    if (($name =~ /prism/) && ($server ne 'SYBTIGR')){
			#
			# Alternative Sybase servers are SYBIL and SYBEST
			#
			my $file = "./conf/Prism.conf";
			
			&edit_prism_conf($file,$server);
		    }
		    
		    
		    my $perlstring = "perl Makefile.PL PREFIX=$installdir WORKFLOW_DOCS_DIR=$installdir/docs SCHEMA_DOCS_DIR=$installdir/docs >& autoinstall.log";

		    &do_or_die($perlstring);
		    
		    my $makestring = "make >> autoinstall.log";
		    
		    &do_or_die($makestring);

		    my $makeinstallstring = "make install >> autoinstall.log";
		    
		    &do_or_die($makeinstallstring);

		    if($name eq 'ergatis'){
			&install_ergatis_src_c($workingdir, $installdir, "HEAD");
		    }

		}else{}

	    }

	if($chmod){
	    my $chmodstrfile = "find $installdir ".'-type f -exec chmod a+rw {} \;';
	    my $chmodstrdir = "find $installdir ".'-type d -exec chmod a+rwx {} \;';
	    my $chmodstrbin = "find $installdir/bin ".'-type f -exec chmod a+rwx {} \;';
	    &do_or_die($chmodstrfile); 
	    &do_or_die($chmodstrdir); 
	    &do_or_die($chmodstrbin); 
	}
}


#---------------------------------------------------------------
# record_installation()
#
#---------------------------------------------------------------
sub record_installation {


    my ($installdir, $workingdir, $username, $date, $server, $cvshash, $datafile, $htmlfile) = @_;
    
    open (DATAFILE, ">>$datafile") or $logger->logdie("Could not open datafile '$datafile' for output: $!");
    
    foreach my $name (sort keys %{$cvshash} ){ 
	
	my $tag = $cvshash->{$name};
	
	print DATAFILE "$installdir||$workingdir||$username||$date||$server||$name||$tag\n";
    }
    

	
    &create_html($datafile, $htmlfile);
    
}


#-------------------------------------------------------------
# create_html()
#
#-------------------------------------------------------------
sub create_html {

	my ($datafile, $htmlfile) = @_;

	if (!defined($datafile)){
	    $logger->logdie("datafile was not defined");
	}
	if (!defined($htmlfile)){
	    $logger->logdie("htmlfile was not defined");
	}
	
	
	open (INFILE, "<$datafile") or $logger->logdie("Could not open file '$datafile': $!");
	open (HTMLFILE, ">$htmlfile") or $logger->logdie("Could not open file '$htmlfile': $!");


	my @contents = <INFILE>;

	chomp @contents;


	my $title = "Installation log";

	&print_header($title);

	&print_body($title, \@contents);

}


#-----------------------------------------------------
# print_body()
#
#-----------------------------------------------------
sub print_body {

    my ($title, $content) = @_;
    
    print HTMLFILE "<body><div align=\"center\"><table border=\"1\" width=\"100%\">\n";
    
    print HTMLFILE "<tr>\n".
    "<th>Install Dir</th>\n".
    "<th>Working Dir</th>\n".
    "<th>Username</th>\n".
    "<th>Date</th>\n".
    "<th>Server</th>\n".
    "<th>CVS codebase</th>\n".
    "<th>CVS tag</th>\n".
    "</tr>\n";
    
    
    foreach my $line (@{$content}) {
	
	my ($installdir, $workingdir, $username, $date, $server, $cvsname, $cvstag) = split(/\|\|/, $line);
	print HTMLFILE "<tr>\n".
	"<td>$installdir</td>\n".
	"<td>$workingdir</td>\n".
	"<td>$username</td>\n".
	"<td>$date</td>\n".
	"<td>$server</td>\n".
	"<td>$cvsname</td>\n".
	"<td>$cvstag</td>\n".
	"</tr>\n";
	
    }

    print HTMLFILE "</table>\n";
    
    my $date = `date`;
    chomp $date;
    
    print HTMLFILE "<br>This page was generated on $date<br>";
    
    print HTMLFILE "</div></body></html>\n";
}



#-----------------------------------------------------
# print_header()
#
#-----------------------------------------------------
sub print_header {

    my ($title) = @_;

    print HTMLFILE <<HeAdER;
	<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
	"http://www.w3.org/TR/html4/strict.dtd">
	
	<html>
	
	<head>
	<meta http-equiv="Content-Language" content="en-us">
	<meta http-equiv="Content-Type" content="text/html; charset=windows-1252">
	<title>$title</title>
	
    <style type="text/css">
	a {
		text-decoration: none;
	  color: rgb(50,50,150);
	  cursor: pointer;
	  cursor: hand;
	}
	td {
		text-align: center;
	}
	</style>
	
HeAdER

}

#-----------------------------------------------------
# do_or_die()
#
#-----------------------------------------------------
sub do_or_die {

    my $cmd = shift;  

    print "$cmd\n";

    my $res =system($cmd);

	if ($res != 0) {
	    $logger->logdie("failed to execute cmd '$cmd': $!");
	    
	} 
}

#-----------------------------------------------------
# edit_prism_conf()
#
#-----------------------------------------------------
sub edit_prism_conf {

    my ($file,$server) = @_;

	my $execstring = "perl -i.bak -pe 's/Chado:BulkSybase:SYBTIGR/Chado:BulkSybase:$server/' $file";

	&do_or_die($execstring);

}



#------------------------------------------------
# verify_create_directory()
#
#------------------------------------------------
sub verify_create_directory {

	my $dir = shift;

	if (-e $dir){

		if (-d $dir){
		   
			if (-w $dir){

				$logger->debug("directory '$dir' exists and has write permissions") if ($logger->is_debug());
			}
			else {
				$logger->logdie("directory '$dir' does not have write permissions");
			}
		}
		else {
			$logger->logdie("'$dir' is not a directory");
		}
	}
	else {
		eval { mkdir($dir); };
		if ($@) {
			$logger->logdie("Could not create directory '$dir': $!");
		}
		else {
			$logger->info("Created directory '$dir'");
		}
	}


}



#-------------------------------------------------
# determine_cvs_revisions()
#
#-------------------------------------------------
sub determine_cvs_revisions {


   
	my ($cvscode, $controlfile) = @_;

	my $cvshash = {};

	if (defined($cvscode)){

		my @groups = split(/;/, $cvscode);

		foreach my $grp (@groups){

			&load_hash($cvshash, $grp);

		}
	}


	if (defined($controlfile)) {
	    print "------------------\nInstalling to $installdir\nwith the following code versions\n\n";
		open (CONTROLFILE, "<$controlfile") or $logger->logdie("Could not open controlfile '$controlfile': $!");

		my @contents = <CONTROLFILE>;

		chomp @contents;


		foreach my $line (@contents){
		    print $line,"\n";
			if ($line =~ /^.+=.+$/) {
				&load_hash($cvshash, $line);
			}
		}
	    print "------------------\n\n\n";
	}


	$logger->debug(Dumper $cvshash) if ($logger->is_debug());
	return $cvshash;
}

#----------------------------------------------------
# load_hash()
#
#----------------------------------------------------
sub load_hash {
	
	my ($cvshash, $grp) = @_;

	my ($name, $tag) = split(/=/, $grp);
	
	if ((defined($name)) && (defined($tag))) {
		
		if (( exists $cvshash->{$name}) && (defined($cvshash->{$name}))) {
			$logger->logdie("name '$name' encountered multiple times!");
		}
		else {
			$cvshash->{$name} = $tag;
		}
	}
	else {
		$logger->logdie("Something was not defined for line '$grp' name '$name' tag '$tag'");
	}
}



#------------------------------------------------------
# print_usage()
#
#------------------------------------------------------
sub print_usage {

    print STDERR "SAMPLE USAGE:  $0 --installdir=installdir --workingdir=workingdir --username=username --server=server --cvscode=\"cvscode\" --controlfile=controlfile [-d debug_level] [-h] [--init] [-l log4perl] [-m]\n".
    "  -i|--installdir          = Installation directory\n".
    "  -w|--workingdir          = Working directory\n".
    "  -U|--username            = Username\n".
    "  -d|--debug_level         = Optional - Coati::Logger log4perl logging level.  Default is 0\n".
    "  -h|--help                = Optional - Display pod2usage help screen\n".
    "  -l|--log4perl            = Optional - Log4perl log file (default: /tmp/cas_installer.pl.log)\n".
    "  -m|--man                 = Optional - Display pod2usage pages for this utility\n".
    "  -S|--server              = Optional - Sybase server for chado databases\n".
    "  --init                   = Optional - if specified the installdir is wiped\n".
    "  --cvscode                = Optional - CVS code and tags\n".
    "  --controlfile            = Optional - control file listing CVS code and tags\n";
    exit 1;

}


sub create_install_docs {

    my ($installdir) = @_;

    &do_or_die("mkdir -p -m 777 $installdir/docs");
    &do_or_die("mkdir -p -m 777 $installdir/docs/ddls/sybase");
    &do_or_die("mkdir -p -m 777 $installdir/docs/ddls/postgresql");
    &do_or_die("mkdir -p -m 777 $installdir/docs/components");
    &do_or_die("mkdir -p -m 777 $installdir/docs/documentation");
    &do_or_die("mkdir -p -m 777 $installdir/docs/obo");
    &do_or_die("mkdir -p -m 777 $installdir/docs/documentation");

}
