#!/usr/local/bin/perl -w
use strict;
use File::Copy;
use File::Basename;
use File::Path; #better choice than system calls, rm -rf = rmtree()
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Cwd qw(realpath);

my $exec_name = basename(realpath($0)); 
my $exec_path = dirname(realpath($0));

umask(0000);

my %opts=();
GetOptions(\%opts, "help|h", "target|t=s", "init|i", "branch-ergatis=s", "branch-chado=s", "branch-prism=s","branch-bsml=s");
#-t: target directory of install
#-i: if the current install in that directory should be (re-)initialized
#  if any branchs are to be used
#-branch-ergatis
#-branch-chado
#-branch-prism
#-branch-bsml"

#check that a directory was supplied and that it is writeable
unless ($opts{'target'}) {
    die "No install directory provided\nusage: $exec_name -t <install_directory> [-i]\n";
}
unless (-d $opts{'target'}) {
    die $opts{'target'}." is not a valid directory\n";
}
unless (-w $opts{'target'}) {
    die $opts{'target'}." is not a writable directory\n";
}

my $install_root = $opts{'target'};
my $tmp_dir = "$install_root/cvstmp";
my $workflow_docs_dir = "$install_root/docs";

#check if branch versions were specified
my $branch_ergatis = (exists $opts{'branch-ergatis'}) ? $opts{'branch-ergatis'} : 'HEAD';
my $branch_chado = (exists $opts{'branch-chado'}) ? $opts{'branch-chado'} : 'HEAD';
my $branch_prism = (exists $opts{'branch-prism'}) ? $opts{'branch-prism'} : 'HEAD';
my $branch_bsml = (exists $opts{'branch-bsml'}) ? $opts{'branch-bsml'} : 'HEAD';

print "ergatis cvs install to $install_root\n";

if ($opts{'init'}) {
    print "Initializing install directory\n";
    #0) chmod directory to be removable
    #doORdie("chmod -R 777 $install_root");
    #1) remove old $install_root
	#doORdie("rm -rf $install_root"); ##PLATFORM SPECIFIC
	rmtreeORdie($install_root); #platform independent
    #2) remake new $cvstmp and $workflow_docs_dir
	#doORdie("mkdir -m 777 -p $workflow_docs_dir"); #-p means "parents", make parents if neccessary, no error if existing
	mkpathORdie($workflow_docs_dir);
	#doORdie("mkdir -m 777 -p $tmp_dir");
	mkpathORdie($tmp_dir);
} else {
    print "Not reinitializing previous install\n";
}

#3) #checkout everything we need from CVS
#checkout and install bsml
chdir($tmp_dir);
doORdie("cvs -Q co -r $branch_bsml bsml_install");
chdir('bsml_install');
doORdie("perl Makefile.PL PREFIX=$install_root SCHEMA_DOCS_DIR=$install_root/docs >& autoinstall.log");
doORdie("make >> autoinstall.log");
doORdie("make install >> autoinstall.log");

#checkout and install ergatis
chdir($tmp_dir);
doORdie("cvs -Q co -r $branch_ergatis ergatis_install");
chdir('ergatis_install');
doORdie("perl Makefile.PL PREFIX=$install_root WORKFLOW_DOCS_DIR=$install_root/docs SCHEMA_DOCS_DIR=$install_root/docs >& autoinstall.log");
doORdie("make >> autoinstall.log");
doORdie("make install >> autoinstall.log");

#checkout and install shared_prism
chdir($tmp_dir);
doORdie("cvs -Q co -r $branch_prism shared_prism_install");
chdir('shared_prism_install');
doORdie("chmod 777 conf/Prism.conf");
#chmod(0777, "conf/Prism.conf");
editPrismConf("conf/Prism.conf");
doORdie("perl Makefile.PL PREFIX=$install_root >& autoinstall.log");
doORdie("make >> autoinstall.log");
doORdie("make install >> autoinstall.log");

#checkout and install prok_prism
chdir($tmp_dir);
doORdie("cvs -Q co -r $branch_prism prok_prism_install");
chdir('prok_prism_install');
doORdie("chmod 777 conf/Prism.conf");
#chmod(0777, "conf/Prism.conf");
editPrismConf("conf/Prism.conf");
doORdie("perl Makefile.PL PREFIX=$install_root >& autoinstall.log");
doORdie("make >> autoinstall.log");
doORdie("make install >> autoinstall.log");

#checkout and install euk_prism
chdir($tmp_dir);
doORdie("cvs -Q co -r $branch_prism euk_prism_install");
chdir('euk_prism_install');
doORdie("chmod 777 conf/Prism.conf");
#chmod(0777, "conf/Prism.conf");
editPrismConf("conf/Prism.conf");
doORdie("perl Makefile.PL PREFIX=$install_root >& autoinstall.log");
doORdie("make >> autoinstall.log");
doORdie("make install >> autoinstall.log");

#checkout and install chado_prism
chdir($tmp_dir);
doORdie("cvs -Q co -r $branch_prism chado_prism_install");
chdir('chado_prism_install');
doORdie("chmod 777 conf/Prism.conf");
#chmod(0777, "conf/Prism.conf");
editPrismConf("conf/Prism.conf");
doORdie("perl Makefile.PL PREFIX=$install_root >& autoinstall.log");
doORdie("make >> autoinstall.log");
doORdie("make install >> autoinstall.log");

#checkout and install coati
chdir($tmp_dir);
doORdie("cvs -Q co -r $branch_prism coati_install");
chdir('coati_install');
doORdie("perl Makefile.PL PREFIX=$install_root >& autoinstall.log");
doORdie("make >> autoinstall.log");
doORdie("make install >> autoinstall.log");

#checkout chado so we can get the DDLs
chdir($tmp_dir);
doORdie("cvs -Q co -r $branch_chado chado");
chdir('chado/tigr_schemas/chado_template/ALL_Modules');
for my $file (glob "*.ddl") {
    copy($file, "$workflow_docs_dir/$file");
}

#checkout ontologies so we can get the DDLs
chdir($tmp_dir);
doORdie("cvs -Q co ontologies");
chdir('ontologies');
for my $file (glob "*.obo") {
    copy($file, "$workflow_docs_dir/$file");
}

#checkout, modify install dir and install PEffect
chdir($tmp_dir);
doORdie("cvs -Q co peffect");
chdir('peffect');
my $old_install_dir = "/usr/local/devel/ANNOTATION/cas";
foreach my $old_file (`grep -rl $old_install_dir *`) {
    chomp($old_file);
    doORdie("chmod 777 $old_file");
    #chmod(0777, $old_file);
    editConfFile($old_file, $old_install_dir, $install_root);
}
doORdie("perl Makefile >& autoinstall.log");
doORdie("make >> autoinstall.log");
doORdie("make install >> autoinstall.log");

#need a few files temporarily
chdir($tmp_dir);
doORdie("cvs -Q co ANNOTATION/chado/cvdata");
chdir("ANNOTATION/chado/cvdata");
copy("pub.out", "$workflow_docs_dir/pub.out");
copy("contact.out", "$workflow_docs_dir/contact.out");

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

sub editPrismConf {
    my $file = shift;
    
    open(my $ifh, "<$file") || die "can't read $file : $!";
    open(my $ofh, ">$tmp_dir/prism.conf.temp") || die "can't create prism temp file : $!";
    
    while (<$ifh>) {
        s/Chado:BulkSybase:SYBTIGR/Chado:BulkSybase:SYBIL/;
        print $ofh $_;
    }
    
    ## move the temp file over the original
    doORdie("mv $tmp_dir/prism.conf.temp $file");
    #move($tmp_dir.'/prism.conf.temp', $file) || die "Unable to move $tmp_dir/prism.conf.temp to $file: $!";
    
}

sub editConfFile {
    (my $file, my $old_text, my $new_text) = @_;
    my $tmp_file = "$tmp_dir/editConfFile.temp";
    
    open(my $ifh, "<$file") || die "can't read $file : $!";
    open(my $ofh, ">$tmp_file") || die "can't create temp file : $!";

    while (<$ifh>) {
        s/$old_text/$new_text/;
        print $ofh $_;
    }
    
    ## move the temp file over the original
	#doORdie("mv $tmp_file $file");
    move($tmp_file, $file) || die "Unable to move $tmp_file to $file: $!";
    doORdie("chmod 777 $file"); #as open() gives 666
#chmod($file, 0777);
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
