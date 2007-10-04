#!/usr/local/bin/perl -w

=head1 NAME

ergatis_installer.pl - Retrieves code from revision control systems using specified tags then installs software in specified install directory

=head1 SYNOPSIS

USAGE:  ergatis_installer.pl --installdir installDirectory --username userInstallingSoftware [ --chmod --controlfile installerControlFile --cvscode cvsCcodeToInstall --debug_level debugLevel -h --init --logfile logFile -m --server serverName --workingdir workingDirectory ]

=head1 REQUIRED ARGUMENTS

=over 8

=item B<--installdir>
    
Directory where software will be installed e.g.: /usr/local/devel/ANNOTATION/ard/

=item B<--username>
    
Name of user installing the software

=head1 OPTIONAL ARGUMENTS

=item B<--chmod>

Specifies that the installdir should be changed to world read/write

=item B<--controlfile>
    
Configuration file for specifying code base and corresponding tags

=item B<--cvscode>
    
Use this command-line argument to specify code base and tag e.g. --cvscode="chado_prism=prism-v1r9b1;ergatis=ergatis-v1r2b1"

=item B<--debug_level>

Log::Log4perl logging level.  Default threshold is WARN

=item B<--help,-h>

Print this help

=item B<--init>

Specifies that the installdir should be wiped  (Default is to not delete the installdir)

=item B<--logfile>

Output Log::Log4perl log filename  (Default is /tmp/ergatis_installer.pl.log)

=item B<--man,-m>

Display the pod2usage page for this utility

=item B<--server>
    
For Prism based code - name of Sybase server to store in configuration files (Default is SYBTIGR)

=item B<--workingdir>
    
Directory where code will be retrieve from the revision control systems and other work will take place (Default is current working directory)

=back

=head1 DESCRIPTION

    ergatis_installer.pl - Checks out code from the various revision control systems using specified tags then installs software in specified install directory

    Assumptions:
    1. User has appropriate permissions (to execute script, access chado database, write to output directory).
    2. All software has been properly installed, all required libraries are accessible.

    Sample usage:
    ./ergatis_installer.pl --username=sundaram --installdir=/usr/local/devel/ANNOTATION/ard --workingdir=/tmp --server=SYBIL --logfile=/tmp/install.log --controlfile=/usr/local/devel/ANNOTATION/ard/chado--v1r14b1.ini

=head1 CONTACT

Jay Sundaram

sundaram@jcvi.org

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
use LWP::UserAgent;

$|=1;

my $tmp_dir;
my %opts;
my $exec_path;

umask(0000);

my $programName = basename($0);

my ($installdir, $workingdir, $username, $help, $logfile, $man, $cvscode,
    $controlfile, $server, $debug_level, $init, $chmod, $datamanagementfile,
    $htmlfile);

my $results = GetOptions (
			  'logfile=s'      => \$logfile,
			  'debug_level=s'  => \$debug_level, 
			  'help|h'         => \$help,
			  'man|m'          => \$man,
			  'installdir=s'   => \$installdir,
			  'workingdir=s'   => \$workingdir,
			  'username=s'     => \$username,
			  'cvscode=s'      => \$cvscode,
			  'controlfile=s'  => \$controlfile,
			  'server=s'       => \$server,
			  'init'           => \$init,
			  'chmod'          => \$chmod,
			  'install_log=s'  => \$datamanagementfile,
			  'install_html=s' => \$htmlfile
			  );

&pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if ($man);
&pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT}) if ($help);

## Check mission critical command-line arguments

my $fatalCtr=0;

if (!defined($username)){
    print STDERR "username was not defined\n";
    $fatalCtr++;
}

if (!defined($installdir)){
    print STDERR "installdir was not defined\n";
    $fatalCtr++;
}

if ($fatalCtr>0){
    &printUsage();
}

if (!defined($logfile)){
    ## Set default log file name
    $logfile = "/tmp/$programName" . '.log';
    print "--logfile was not specified, setting to default value '$logfile'\n";
}

my $logger = &getLogger($logfile, $programName);

if (!defined($server)){
    $server = 'SYBTIGR';
    $logger->info("--server was not specified, setting to default value '$server'");
}

&verify_create_directory($installdir);

&clear_install_dir($installdir, $init);

&create_install_docs($installdir);

if (!defined($workingdir)){
    $workingdir = '/tmp';
    if ($logger->is_debug()){
	$logger->debug("--workingdir was not specified, setting to default value '$workingdir'");
    }
}

&verify_create_directory($workingdir);

&clear_working_dir($workingdir);

my $cvshash = &determine_cvs_revisions($cvscode, $controlfile);

&execute_installation($cvshash, $installdir, $workingdir, $server);

if ((defined($datamanagementfile)) && (defined($htmlfile))){
    &record_installation($installdir, $workingdir, $username, $server, $cvshash, $datamanagementfile, $htmlfile);
}

if (defined($chmod)){
    &setPermissions($installdir);
}

print "$0 execution completed\n";
print "Log file is '$logfile'\n";
exit(0);

#--------------------------------------------------------------------------------------------------------------------------
#
#                                     END OF MAIN    --    SUBROUTINES FOLLOW
#
#--------------------------------------------------------------------------------------------------------------------------

=over 4

=item clear_working_dir()

B<Description:> Will remove all temporary working sub-directories

B<Parameters:> $workingdir (scalar)

B<Returns:> None

=back

=cut

sub clear_working_dir {
	
    my ($workingdir) = @_;
	
    &do_or_die("rm -rf $workingdir/*_install");

    &do_or_die("rm -rf $workingdir/ergatis_c");

    &do_or_die("rm -rf $workingdir/chado_schema_sybase");

    &do_or_die("rm -rf $workingdir/chado_schema_postgresql");

    &do_or_die("rm -rf $workingdir/peffect");

    &do_or_die("rm -rf $workingdir/ontologies");

    &do_or_die("rm -rf $workingdir/logcabin");
}
	
=over 4

=item clear_install_dir()

B<Description:> Re-establish the install directory

B<Parameters:> $dir (scalar), $init (scalar)

B<Returns:> None

=back

=cut

sub clear_install_dir {
	
    my ($dir, $init) = @_;
    
    if ((defined($init)) && ($init == 1)){
	
	&do_or_die("rm -rf $dir");
	
	mkdir($dir);
    }

    my $docsdir = "$dir/docs";

    if (!-e $docsdir) {
	mkdir($docsdir);
    }
}

=over 4

=item execute_installation()

B<Description:> Execute the retrieval from revision systems and then run make install

B<Parameters:> $cvshash (reference to hash), $installdir (scalar), $workingdir (scalar), $server (scalar)

B<Returns:> None

=back

=cut

sub execute_installation {

    my ($cvshash, $installdir, $workingdir, $server) = @_;
    
    foreach my $name (sort keys %{$cvshash} ) {
	
	my $tag = $cvshash->{$name};
	
	if (!defined($tag)){
	    $logger->logdie("tag was not defined for name '$name'");
	}

	my $installname = $name ."_install";

	chdir($workingdir);
	
	if ($name eq 'ontologies'){
	    &install_ontologies($installdir, $tag);
	}
	elsif ($name eq 'chado_schema'){
	    &install_schema($workingdir, $installdir, $tag);
	}
	elsif ($name eq 'peffect'){
	    &install_peffect($installdir);
	}
	elsif( $name eq 'ergatis'){

	    if ($tag eq 'HEAD'){
		$logger->logdie("This version '$tag' is not available as a tarball release.  Contact jorvis\@gmail.com or sundaram\@jcvi.org.");
	    }

	    &downloadErgatisTarball($workingdir, $tag, $installdir, $installname);

	    #&installErgatis($installdir, $tag, $installname);
	    
	    #&installErgatisSrcC($workingdir, $installdir, $tag, $installname);

	}
	elsif( $name eq 'bsml'){
	    &installBsml($installdir, $tag, $installname);
	}
	elsif ($name eq 'coati'){
	    &installCoati($installdir, $tag, $installname);
	}
	elsif ($name =~ /prism/){
	    ## The server value in the configuration file overrides the 
	    ## server value specified in the command-line argument.
	    if ((exists $cvshash->{'server'}) && (defined($cvshash->{'server'})) ) {
		$server = $cvshash->{'server'};
	    }

	    &installPrism($installdir, $tag, $installname, $server);
	}
	else {
	    if (lc($name) eq 'server'){
		next;
	    }
	    $logger->logdie("Don't know how to process project '$name' with tag '$tag'");
	}
    }
}

sub install_ontologies {

    my ($installdir, $tag) = @_;

    &do_or_die("cvs -z3 -d:ext:$username\@ergatis.cvs.sourceforge.net:/cvsroot/ergatis co -r $tag chado");

    &do_or_die("cp chado/ontologies/*.obo $installdir/docs/obo/.");
}

sub install_schema {

    my ($workingdir, $installdir, $tag) = @_;

    if (!-e './chado'){
	## Might have already been retrieved when the ontologies where 'installed'.
	&do_or_die("cvs -z3 -d:ext:$username\@ergatis.cvs.sourceforge.net:/cvsroot/ergatis co -r $tag chado");
    }

    &do_or_die("cp $workingdir/chado/sybase/ddls/*.ddl $installdir/docs/ddls/sybase/.");

    &do_or_die("cp $workingdir/chado/postgresql/ddls/*.ddl $installdir/docs/ddls/postgresql/.");
}


sub install_peffect {

    my ($installdir) = @_;

    &do_or_die("cvs -Q export -kkv -r HEAD -d peffect peffect");
    
    chdir("peffect");
    
    &do_or_die("find . -name Makefile -exec perl -i.bak -pe 's|/usr/local/devel/ANNOTATION/cas|$installdir|' \{\} \\;");
    
    &do_or_die("make");
    
    &do_or_die("make install");
}

sub installErgatis {

    my ($installdir, $tag, $installname) = @_;

    &do_or_die("svn co --non-interactive https://ergatis.svn.sourceforge.net/svnroot/ergatis/$tag $installname");
	
    &perlMakeInstall($installdir, $tag);
}	

sub installErgatisSrcC {

    my ($workingdir, $installdir, $tag) = @_;

    chdir($workingdir);

    if (!-e './ergatis_installer/src/c/SnpClusterer'){

	&do_or_die("svn co --non-interactive https://ergatis.svn.sourceforge.net/svnroot/ergatis/$tag ergatis_installer");
    }

    chdir('ergatis_installer/src/c/SnpClusterer');
    
    &do_or_die("make");
    
    &do_or_die("cp SnpClusterer $installdir/bin/.");
    
    chdir($workingdir);

    if (!-e './ergatis_installer/src/c/chado_record_uniq'){

	&do_or_die("svn co --non-interactive https://ergatis.svn.sourceforge.net/svnroot/ergatis/$tag ergatis_installer");
    }

    chdir('ergatis_installer/src/c/chado_record_uniq');

    &do_or_die("make");

    &do_or_die("cp chado_record_uniq $installdir/bin/.");

}

sub installBsml {

    my ($installdir, $tag, $installname) = @_;

    &do_or_die("cvs -z3 -d:ext:$username\@bsml.cvs.sourceforge.net:/cvsroot/bsml co -r $tag $installname");
	
    &perlMakeInstall($installdir, $installname);

}	

sub installCoati {

    my ($installdir, $tag, $installname) = @_;
    
    &do_or_die("cvs -z3 -d:ext:$username\@manatee.cvs.sourceforge.net:/cvsroot/manatee co -r $tag $installname");
	
    &perlMakeInstall($installdir, $installname);

}	

sub installPrism {

    my ($installdir, $tag, $installname, $server) = @_;

    &do_or_die("cvs -z3 -d:ext:$username\@ergatis.cvs.sourceforge.net:/cvsroot/ergatis co -r $tag $installname");

    if ($server ne 'SYBTIGR'){
	## Alternative Sybase servers are SYBIL and SYBEST

	&do_or_die("perl -i.bak -pe 's/Chado:BulkSybase:SYBTIGR/Chado:BulkSybase:$server/' ./conf/Prism.conf");
    }
    
    &perlMakeInstall($installdir, $installname);
}	

sub perlMakeInstall {

    my ($installdir, $installname) = @_;

    chdir($installname);
	
    &do_or_die("perl Makefile.PL PREFIX=$installdir WORKFLOW_DOCS_DIR=$installdir/docs SCHEMA_DOCS_DIR=$installdir/docs >& autoinstall.log");
    
    &do_or_die("make >> autoinstall.log");

    &do_or_die("make install >> autoinstall.log");
}


=over 4

=item record_installation()

B<Description:> Creates a record of the installation instance

B<Parameters:> $installdir (scalar), $workingdir (scalar), $username (scalar), 
$server (scalar), $cvshash (reference to hash), $datafile (scalar),
$htmlfile (scalar)

B<Returns:> None

=back

=cut

sub record_installation {

    my ($installdir, $workingdir, $username, $server, $cvshash, $datafile, $htmlfile) = @_;
    
    my $date = localtime();

    open (DATAFILE, ">>$datafile") or $logger->logdie("Could not open datafile '$datafile' for output: $!");
    
    foreach my $name (sort keys %{$cvshash} ){ 
	
	my $tag = $cvshash->{$name};
	
	print DATAFILE "$installdir||$workingdir||$username||$date||$server||$name||$tag\n";
    }
    
    &create_html($datafile, $htmlfile);
}


=over 4

=item create_html()

B<Description:> Create HTML log file

B<Parameters:> $datafile (scalar), $htmlfile (scalar)

B<Returns:> None

=back

=cut

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


=over 4

=item print_body()

B<Description:> Print HTML body information

B<Parameters:> $title (scalar), $content (scalar)

B<Returns:> None

=back

=cut

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
    
    my $date = localtime();
    
    print HTMLFILE "<br>This page was generated on $date<br>";
    
    print HTMLFILE "</div></body></html>\n";
}

=over 4

=item print_header()

B<Description:> Print HTML header information

B<Parameters:> $title (scalar)

B<Returns:> None

=back

=cut

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

=over 4

=item do_or_die()

B<Description:> Executes Perl system command

B<Parameters:> $cmd (scalar)

B<Returns:> None

=back

=cut

sub do_or_die {

    my $cmd = shift;  
    
    print "$cmd\n";
    
    my $res =system($cmd);
    
    if ($res != 0) {
	$logger->logdie("failed to execute cmd '$cmd': $!");
	
    } 
}

=over 4

=item verify_create_directory()

B<Description:> Verify the status of the directory or create it is does not exist

B<Parameters:> $dir (scalar)

B<Returns:> None

=back

=cut

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

=over 4

=item determine_cvs_revisions()

B<Description:> Process the cvscode and controlfile values and load hash with key-value pairs of revision name and revision tags

B<Parameters:> $cvscode (scalar), $controlfile (scalar)

B<Returns:> $cvshash (reference to hash)

=back

=cut

sub determine_cvs_revisions {

    my ($cvscode, $controlfile) = @_;
    
    my $cvshash = {};
    
    if (defined($cvscode)){
	## Some command-line value was provided

	my @groups = split(/;/, $cvscode);
	
	foreach my $grp (@groups){
	    
	    &load_hash($cvshash, $grp);
	    
	}
    }

    if (defined($controlfile)) {
	## Some control file was specified

	print "------------------\nInstalling to $installdir\nwith the following code versions\n\n";
	open (CONTROLFILE, "<$controlfile") or $logger->logdie("Could not open controlfile '$controlfile': $!");
	
	my @contents = <CONTROLFILE>;
	
	chomp @contents;
	
	
	foreach my $line (@contents){

	    if ($line =~ m/^\#/){
		next;
	    }

	    print $line,"\n";
	    if ($line =~ /^.+=.+$/) {
		&load_hash($cvshash, $line);
	    }
	}
	print "------------------\n\n\n";
    }



    if ($logger->is_debug()){

	$logger->debug("Here are the revision names and corresponding tags");

	foreach my $name (sort keys %{$cvshash}){
	    $logger->debug("$name \t $cvshash->{$name}\n");
	}
    }

    return $cvshash;
}


=over 4

=item load_hash()

B<Description:> Verifies that the revision name and tag are defined and that the name has not been previously loaded into referenced hash

B<Parameters:> $cvshash (reference to hash), $grp (scalar)

B<Returns:> None

=back

=cut

sub load_hash {
    
    my ($cvshash, $grp) = @_;
    
    my ($name, $tag) = split(/=/, $grp);
    
    if (!defined($name)){
	$logger->logdie("name was not defined while processing '$grp'");
    }

    if (!defined($tag)){
	$logger->logdie("tag was not defined while processing '$grp'");
    }
    
    if ( exists $cvshash->{$name}){
	$logger->logdie("name '$name' encountered multiple times!");
    }

    $cvshash->{$name} = $tag;

}

=over 4

=item printUsage()

B<Description:> Prints the usage information to STDERR and then exits with exit code 0

B<Parameters:> None

B<Returns:> None

=back

=cut

sub printUsage {

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
    exit(1);

}


=over 4

=item create_install_docs()

B<Description:> Will create the docs/ directory and all appropriate sub-directories

B<Parameters:> None

B<Returns:> None

=back

=cut

sub create_install_docs {

    my ($installdir) = @_;

    if (!defined($installdir)){
	$logger->logdie("installdir was not defined!");
    }

    &do_or_die("mkdir -p -m 777 $installdir/docs");
    &do_or_die("mkdir -p -m 777 $installdir/docs/ddls/sybase");
    &do_or_die("mkdir -p -m 777 $installdir/docs/ddls/postgresql");
    &do_or_die("mkdir -p -m 777 $installdir/docs/components");
    &do_or_die("mkdir -p -m 777 $installdir/docs/documentation");
    &do_or_die("mkdir -p -m 777 $installdir/docs/obo");

}

=over 4

=item getLogger()

B<Description:> Will initialize the Log4perl Logger

B<Parameters:> $logfile (scalar), $verbose (scalar), $programName (scalar)

B<Returns:> $logger (reference to Log4perl Logger)

=back

=cut

sub getLogger {

    my ($logfile, $programName) = @_;
    
    ## Initialize the logger
    my $verbose = 1;
    
    my $screen_threshold = 'WARN';
    if ($verbose){
	$screen_threshold = 'INFO';
    }
    
    Log::Log4perl->init(
			\ qq{
			    log4perl.logger                       = WARN, A1
			    log4perl.appender.A1                  = Log::Dispatch::File
			    log4perl.appender.A1.filename         = $logfile
			    log4perl.appender.A1.mode             = write
			    log4perl.appender.A1.Threshold        = WARN
			    log4perl.appender.A1.layout           = Log::Log4perl::Layout::PatternLayout
			    log4perl.appender.A1.layout.ConversionPattern = %p> %L %M - %m%n 
			    log4perl.appender.Screen              = Log::Dispatch::Screen
			    log4perl.appender.Screen.layout       = Log::Log4perl::Layout::SimpleLayout
			    log4perl.appender.Screen.Threshold    = $screen_threshold
			Log::Log4perl::SimpleLayout
		    }
			);

    
    return get_logger($programName);
}

=over 4

=item setPermissions()

B<Description:> Will set appropriate permissions on files and directories nested in install directory

B<Parameters:> $installdir (scalar)

B<Returns:> None

=back

=cut

sub setPermissions{
    
    my ($installdir) = @_;
    
    my $chmodstrfile = "find $installdir ".'-type f -exec chmod a+rw {} \;';
    &do_or_die($chmodstrfile); 

    my $chmodstrdir = "find $installdir ".'-type d -exec chmod a+rwx {} \;';
    &do_or_die($chmodstrdir); 

    my $chmodstrbin = "find $installdir/bin ".'-type f -exec chmod a+rwx {} \;';
    &do_or_die($chmodstrbin); 

}

=over 4

=item downloadErgatisTarball()

B<Description:> Will retrieve ergatis tarball from sourceforge

B<Parameters:> $workingdir (scalar), $tag (scalar), $installdir (scalar), $installname (scalar)

B<Returns:> None

=back

=cut

sub downloadErgatisTarball {
    
    my ($workingdir, $tag, $installdir, $installname) = @_;

    mkdir("$workingdir/$installname");

    chdir("$workingdir/$installname");

    
    my $url = "http://downloads.sourceforge.net/ergatis/$tag" . ".tar.gz";

    &do_or_die("wget $url");

    &do_or_die("tar zxvf $tag.tar.gz");

    &perlMakeInstall($installdir, $tag);

}
