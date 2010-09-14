use ExtUtils::MakeMaker;
use ExtUtils::MakeMaker::Config;
require ExtUtils::Manifest;
use Term::ANSIColor;
require "perl2wrapper_ergatis.pl";
use strict;

$|++;

## always remake the MANIFEST
ExtUtils::Manifest::mkmanifest();

## let's start by checking for the modules we need
my @modules_to_check = qw( 
    Bio::Seq
    CDB_File
    Class::Struct
    Config::IniFiles
    DB_File
    Data::Dumper
    Date::Manip
    File::Basename
    File::Copy
    File::Find
    File::Mirror
    File::Path
    File::Spec
    Getopt::Long
    HTML::Template
    IO::File
    IO::Tee
    IPC::Open3
    LWP::Simple
    Log::Cabin
    Log::Log4perl
    Mail::Mailer
    Math::Combinatorics
    MLDBM
    PerlIO::gzip
    Storable
    URI::Escape
    XML::LibXML
    XML::Parser
    XML::RSS
    XML::Twig
    XML::Writer
);
my @modules_not_found;

print "checking necessary perl modules:\n";

for ( @modules_to_check ) {
    my $formed = $_;
    $formed =~ s|::|/|g;
    $formed =~ s|$|.pm|g;

    print "\t$_\t:";
    eval {require $formed};

    if ( $@ ) {
        print color 'bold red';
        print " not found\n";
        push @modules_not_found, $_;
    } else {
        print color 'bold green';
        print " found\n";
    }

    print color 'reset';
}

if ( scalar @modules_not_found && !($ARGV[0] eq 'force')) {
    print "\ninstall aborted.  the following modules are missing: \n";
    for ( @modules_not_found ) {
        print "\t$_\n";
    }

    exit(1);
}

## make sure that their ExtUtils::MakeMaker is a high enough version
my $makemaker_found = 0;

eval "require ExtUtils::MakeMaker";

if ( ExtUtils::MakeMaker->VERSION < 6.31 ) {
    print color 'bold red';
    print "ExtUtils::MakeMaker >= 6.31 not found\n";
    print color 'reset';
    exit(1);
}

print "\nall modules found, carrying on\n\n";

my $install_base;
my $build_dir = '.';
my $cbuild_dir = 'src';
my $tbuild_dir = 'global_pipeline_templates';
my $samplebuild_dir = 'samples';
my $sbuild_dir = 'shell';


foreach my $arg (@ARGV){
    if ($arg =~ /INSTALL_BASE/) {
	($install_base) = ($arg =~ /INSTALL_BASE=(\S+)/);
    } elsif ( $arg =~ /PREFIX/ ) {
        print "\nERROR: please use INSTALL_BASE instead of PREFIX.  (This requires a current version of ExtUtils::MakeMaker";
        exit(1);
    } elsif ( $arg =~/LIVE_BUILD/){
	$build_dir = "../";
	$cbuild_dir = "../src/c";
	$tbuild_dir = "../templates";
	$samplebuild_dir = "../samples";
	$sbuild_dir = "../src/shell";
    }
}

my $workflowdocsdir = "$install_base/docs";

if (! -e $workflowdocsdir ) {
    mkdir($workflowdocsdir) || die "failed to create WORKFLOW_DOCS_DIR: $!";
}

my $sampledocsdir = "$install_base/samples";

## workflow and schema docs dirs have always been the same
my $schemadocsdir = $workflowdocsdir;

my($scripts,$wrapperstr) = &get_bins();
print "Writing makefile to install directory $build_dir to $install_base\n";
my $mm = WriteMakefile(
	      'NAME'=>'ergatis',
	      'EXE_FILES'=>$scripts,
	      macro=>{BUILD_DIR => "$build_dir",C_BUILD_DIR=>"$cbuild_dir",T_BUILD_DIR=>"$tbuild_dir",SAMPLE_BUILD_DIR=>"$samplebuild_dir",S_BUILD_DIR=>"$sbuild_dir"},
	      'clean' =>{FILES => "MANIFEST.bak $wrapperstr README.symlinks"}
	      );

package MY;

sub MY::postamble { 
    my $postamble = "WORKFLOW_DOCS_DIR=$workflowdocsdir\nSAMPLE_DOCS_DIR=$sampledocsdir";
    $postamble .= << 'END';

ergatis_install : dirs samples shellwrappers c_apps pipeline_templates 
	@echo "Installing components"
	find $(BUILD_DIR)/components -path '*svn*' -prune -o -print -type f -not -name software.config -not -name "*.tmpl" -exec install  {} ${WORKFLOW_DOCS_DIR} \;
	find $(BUILD_DIR)/components -type f -name "*.tmpl" -exec install  {} ${WORKFLOW_DOCS_DIR}/documentation \;
	@install $(BUILD_DIR)/components/shared/software.config $(PREFIX)

dirs	:
	mkdir -p $(WORKFLOW_DOCS_DIR)
	mkdir -p $(WORKFLOW_DOCS_DIR)/documentation
	mkdir -p $(SAMPLE_DOCS_DIR)

samples	: dirs
	@echo "Installing samples"
	find $(SAMPLES_BUILD_DIR) -type f -exec install {} $(SAMPLE_DOCS_DIR) \;

shellwrappers:
	@echo "Installing shell wrappers"
	find $(S_BUILD_DIR) -type f -exec install {} ${PREFIX}/bin \;

pipeline_templates:
	@echo "Installing global pipeline templates in $(PREFIX)/global_pipeline_templates"
	find $(T_BUILD_DIR) -type d -exec mkdir -p $(PREFIX)/{} \;
	find $(T_BUILD_DIR) -path '*svn*' -prune -o -print -type d -exec sh -c "echo \"{}/* $(PREFIX)/{}\"" \; 

c_apps	:
	@echo "Compiling C applications and installing in $(PREFIX)/bin"
	${MAKE} -C $(C_BUILD_DIR)/cluster
	@install $(C_BUILD_DIR)/cluster/cluster $(PREFIX)/bin
	${MAKE} -C $(C_BUILD_DIR)/chado_record_uniq
	@install $(C_BUILD_DIR)/chado_record_uniq/chado_record_uniq $(PREFIX)/bin
	${MAKE} -C $(C_BUILD_DIR)/pangenome_profile_compare
	@install $(C_BUILD_DIR)/pangenome_profile_compare/pangenome_profile_compare $(PREFIX)/bin
	${MAKE} -C $(C_BUILD_DIR)/SnpClusterer
	@install $(C_BUILD_DIR)/SnpClusterer/SnpClusterer $(PREFIX)/bin

r_scripts:
	@echo "Installing R scripts to $(PREFIX)/docs"
	@tar cf - R |  (cd $(WORKFLOW_DOCS_DIR); tar xf -)

END

    $postamble;
}

sub install { 
    my $self = shift;
    my $installstr = $self->MM::install(@_);
    my $addltarget = 'ergatis_install';
    $installstr =~ s/(pure_install\s+::\s+)(.*)/$1 $2 $addltarget/;
    return $installstr;
}
