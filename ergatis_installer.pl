#!/usr/local/bin/perl -w
#---------------------------------------------------------------------------
#
#
#
#
#
#---------------------------------------------------------------------------
=head1 NAME

cas_installer.pl - Checks out code from CVS repository using specified tags then installs software in specified install directory

=head1 SYNOPSIS

USAGE:  cas_installer.pl -i installdir --installdir=D database -P password -U username [-a autogen_feat] -b bsmldoc [-d debug_level] [-h] [-i insert_new] [-l log4perl] [-m] [-o outdir] [-p] [-R readonlycache] [-s autogen_seq] [-u update] [-x xml_schema_type] [-y cache_dir] [-z doctype]

=head1 OPTIONS

=over 8

=item B<--username,-U>
    
    Database username

=item B<--password,-P>
    
    Database password

=item B<--database,-D>
    
    Target chado database 

=item B<--bsmldoc,-b>
    
    Bsml document containing pairwise alignment encodings

=item B<--autogen_feat,-a>
    
    Optional - Default behavior is to auto-generate (-a=1) chado feature.uniquename values for all inbound features.  To turn off behavior specify this command-line option (-a=0).

=item B<--autogen_seq,-s>
    
    Optional - Default behavior is to not (-s=0) auto-generate chado feature.uniquename values for all inbound sequences.  To turn on behavior specify this command-line option (-s=1).

=item B<--insert_new,-i>
    
    Optional - Default behavior is to insert (-i=1) insert newly encountered Sequence objects in the BSML document that are not currently present in the Chado database.  To turn off default insert behavior specify this command-line option (-i=0)

=item B<--debug_level,-d>

    Optional: Coati::Logger log4perl logging level.  Default is 0

=item B<--man,-m>

    Display the pod2usage page for this utility

=item B<--outdir,-o>

    Optional: Output directory for the tab delimited .out files.  Default is current directory

=item B<--pparse,-p>

    Optional - turn off parallel load support via global serial identifier replacement (default is ON)

=item B<--update,-u>

    Optional - Default behavior is to not update the database (-u=0).  To turn on update behavior specify this command-line option (-u=1).


=item B<--cache_dir,-y>

    Optional - Query caching directory to write cache files (default is ENV{DBCACHE_DIR})

=item B<--readonlycache,-R>

    Optional - If data file caching is activated and if this readonlycache is == 1, then the tied MLDBM lookup cache files can only be accessed in read-only mode.  Default (-r=0) means cached lookup can be created and access mode is read-write.

=item B<--doctype,-z>

    Optional - If specified, can direct the parser to construct concise lookup - more efficient. One of the following: nucmer, region, promer, pe, blastp, repeat, scaffold, rna, te, coverage

=item B<--help,-h>

    Print this help

=back

=head1 DESCRIPTION

    bsml2chado.pl - Parse BSML document and produce tab delimited .out BCP files for insertion into Chado database

    Assumptions:
    1. The BSML pairwise alignment encoding should validate against the XML schema:.
    2. User has appropriate permissions (to execute script, access chado database, write to output directory).
    3. Target chado database already contains all reference features (necessary to build feature and organism lookups) Review and execute db2bsml.pl if required.
    4. Target chado database contains the necessary controlled vocabulary terms: "match" etc.
    5. All software has been properly installed, all required libraries are accessible.

    Sample usage:
    ./bsml2chado.pl -U access -P access -D tryp -b /usr/local/annotation/TRYP/BSML_repository/blastp/lma2_86_assembly.blastp.bsml  -l my.log -o /tmp/outdir


=cut

use strict;
use File::Copy;
use File::Basename;
use File::Path; #better choice than system calls, rm -rf = rmtree()
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Cwd qw(realpath);
use Pod::Usage;
use Data::Dumper;
use Log::Log4perl qw(get_logger);


$|=1;

my $tmp_dir;
my %opts;
my $exec_path;


my $datamanagementfile = "/usr/local/devel/ANNOTATION/cas/datamanagement/cas_install.txt";
my $htmlfile = "/usr/local/devel/ANNOTATION/cas/datamanagement/cas_install.html";

umask(0000);


my ($installdir, $workingdir, $username, $help, $log4perl, $man, $cvscode, $controlfile, $server, $debug_level);

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
						  'server|S=s'         => \$server
						  );

&pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if ($man);
&pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT}) if ($help);

print STDERR ("username was not defined\n")     if (!$username); 
print STDERR ("installdir was not defined\n")   if (!$installdir); 
print STDERR ("workingdir was not defined\n")   if (!$workingdir); 

&print_usage if(!$username or !$installdir or !$workingdir);

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


my $logger = get_logger("GFF2BSML");


if (!defined($server)){
	$server = 'SYBTIGR';
	$logger->info("server was set to '$server'");
}


&verify_create_directory($installdir);

&clear_install_dir($installdir);

#
# Create the docs directory
#
my $docsdir = "$installdir/docs";
mkdir($docsdir);

&verify_create_directory($workingdir);

&clear_working_dir($workingdir);

my $cvshash = &determine_cvs_revisions($cvscode, $controlfile);


my $date = `date`;

chomp $date;

&execute_installation($cvshash, $installdir, $workingdir, $server);

&install_ontologies($workingdir, $installdir, "HEAD");
&install_schema($workingdir, $installdir, "HEAD");
&install_bcp_files($workingdir, $installdir, "HEAD");
&install_peffect($workingdir, $installdir, "HEAD");

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

    my $execstring = "cvs -Q co -r $tag ontologies";

    &do_or_die($execstring);

    $execstring = "cp ontologies/*.obo $installdir/docs";

    &do_or_die($execstring);
}




#--------------------------------------------------
# install_schema()
#
#--------------------------------------------------
sub install_schema {

    my ($workingdir, $installdir, $tag) = @_;

    chdir($workingdir);

    my $execstring = "cvs -Q co -d chado_schema -r $tag ANNOTATION/chado/tigr_schemas/chado_template/ALL_Modules";

    &do_or_die($execstring);

    $execstring = "cp chado_schema/*.ddl chado_schema/*.sql $installdir/docs";

    &do_or_die($execstring);
}




#--------------------------------------------------
# install_bcp_files()
#
#--------------------------------------------------
sub install_bcp_files {

    my ($workingdir, $installdir, $tag) = @_;

    chdir($workingdir);

    my $execstring = "cvs -Q co -d cvdata -r $tag ANNOTATION/chado/cvdata";

    &do_or_die($execstring);

    $execstring = "cp cvdata/contact.out cvdata/pub.out $installdir/docs";

    &do_or_die($execstring);
}


#--------------------------------------------------
# install_peffect()
#
#--------------------------------------------------
sub install_peffect {

    my ($workingdir, $installdir, $tag) = @_;

    chdir($workingdir);

    my $execstring = "cvs -Q co -d peffect peffect";

    &do_or_die($execstring);

    chdir("peffect");

    $execstring = "find . -name Makefile -exec \"perl -i.bak -pe 's|/usr/local/devel/ANNOTATION/cas|$installdir|' \{\} \;\"";
    
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

}


#----------------------------------------------------
# clear_install_dir()
#
#----------------------------------------------------
sub clear_install_dir {
	
	my ($dir) = @_;
	
	my $execstring = "rm -rf $dir";
	
	&do_or_die($execstring);

	mkdir($dir);

}


#-----------------------------------------------------
# execute_installation()
#
#-----------------------------------------------------
sub execute_installation {

	my ($cvshash, $installdir, $workingdir, $server) = @_;

	foreach my $name (sort keys %{$cvshash} ) {


		my $tag = $cvshash->{$name};

		$logger->logdie("tag was not defined for name '$name'") if (!defined($tag));

		my $installname = $name ."_install";


		chdir($workingdir);
		my $execstring = "cvs -Q co -r $tag $installname";

		&do_or_die($execstring);

		chdir($installname);
		
		if (($name =~ /prism/) && ($server eq 'SYBIL')){

			my $file = "./conf/Prism.conf";
			
			&edit_prism_conf($file);
		}


		my $perlstring = "perl Makefile.PL PREFIX=$installdir WORKFLOW_DOCS_DIR=$installdir/docs SCHEMA_DOCS_DIR=$installdir/docs >& autoinstall.log";

        #die "would have executed $perlstring\n";

		&do_or_die($perlstring);

		my $makestring = "make >> autoinstall.log";

		&do_or_die($makestring);

		my $makeinstallstring = "make install >> autoinstall.log";

		&do_or_die($makeinstallstring);
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



	if (-e $htmlfile){
		
	    my $htmlbak = $htmlfile . ".$$.bak";

	    rename ($htmlfile, $htmlbak);
	    
	    chmod (0666, $htmlbak);

	}

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


	my $title = "CAS Installation";

	&print_header($title);

	&print_body($title, \@contents);


	if (0) {
	    #
	    # Will phase out the invocation of cas2html.pl
	    #


#	die "datafile '$datafile' htmlfile '$htmlfile'";
	    
	    eval { qx{./cas2html.pl $datafile $htmlfile}; };
	    
	    if ($@){
		$logger->logdie("$!");
	    }
	    
	}








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



    if (1){
	#
	# Control testing
	#


	system($cmd);

	if ($? == -1) {
	    $logger->logdie("failed to execute cmd '$cmd': $!");
	    
	} 
	elsif ($? & 127) {
	    printf STDERR "child died with signal %d, %s coredump\n",
	    ($? & 127),  ($? & 128) ? 'with' : 'without';
	    exit(1);
	}
    }


}

#-----------------------------------------------------
# edit_prism_conf()
#
#-----------------------------------------------------
sub edit_prism_conf {

    my ($file) = @_;

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

		open (CONTROLFILE, "<$controlfile") or $logger->logdie("Could not open controlfile '$controlfile': $!");

		my @contents = <CONTROLFILE>;

		chomp @contents;


		foreach my $line (@contents){
			
			&load_hash($cvshash, $line);
		}

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

    print STDERR "SAMPLE USAGE:  $0 --installdir=installdir --workingdir=workingdir --username=username --server=server --cvscode=\"cvscode\" --controlfile=controlfile [-d debug_level] [-h] [-l log4perl] [-m]\n".
    "  -i|--installdir          = Installation directory\n".
    "  -w|--workingdir          = Working directory\n".
    "  -U|--username            = Username\n".
    "  -d|--debug_level         = Optional - Coati::Logger log4perl logging level.  Default is 0\n".
    "  -h|--help                = Optional - Display pod2usage help screen\n".
    "  -l|--log4perl            = Optional - Log4perl log file (default: /tmp/cas_installer.pl.log)\n".
    "  -m|--man                 = Optional - Display pod2usage pages for this utility\n".
    "  -S|--server              = Optional - Sybase server for chado databases\n".
    "  --cvscode                = Optional - CVS code and tags\n".
    "  --controlfile            = Optional - control file listing CVS code and tags\n";
    exit 1;

}




