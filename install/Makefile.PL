use ExtUtils::MakeMaker;
require ExtUtils::Manifest;
require "perl2wrapper_papyrus.pl";
use strict;

if(-e "MANIFEST"){}
else{
    ExtUtils::Manifest::mkmanifest();
  }


my $workflowdocsdir = '';
my $schemadocsdir = '';
my $isworkflowdocsdir=0;
my $isschemadocsdir=0;

foreach my $arg (@ARGV){
	if($arg =~ /WORKFLOW_DOCS_DIR/){
	    ($workflowdocsdir) = ($arg =~ /WORKFLOW_DOCS_DIR=(.*)/);
	    $isworkflowdocsdir=1;
	}
	if($arg =~ /SCHEMA_DOCS_DIR/){
	    ($schemadocsdir) = ($arg =~ /SCHEMA_DOCS_DIR=(.*)/);
	    $isschemadocsdir=1;
	}
    }

if(!$isworkflowdocsdir){
    die "You must specify a WORKFLOW_DOCS_DIR.\n".
	"Workflow XML and ini files will be stored in this directory\n".
	"Eg. \n\nperl Makefile.PL PREFIX=/usr/local/devel/ANNOTATION/shared WORKFLOW_DOCS_DIR=/usr/local/devel/ANNOTATION/shared/docs\n\n";
}
else{
    print "Ignore the following warning: ";
}
if(!$isschemadocsdir || (!(-e $schemadocsdir))){
    die "You must specify a SCHEMA_DOCS_DIR where BSML XML schemas and the BSML dtd are installed.\n".
	"Eg. \n\nperl Makefile.PL PREFIX=/usr/local/devel/ANNOTATION/shared WORKFLOW_DOCS_DIR=/usr/local/devel/ANNOTATION/shared/docs SCHEMA_DOCS_DIR=/usr/local/devel/ANNOTATION/shared/docs\n\n";
}
else{
    print "Ignore the following warning: ";
}

my($scripts,$wrapperstr) = &get_bins();


WriteMakefile(
	      'NAME'=>'papyrus',
	      'EXE_FILES'=>$scripts,
	      'clean' =>{FILES => "MANIFEST.bak $wrapperstr README.symlinks"}
	      );

package MY;

sub MY::postamble { 
    my $postamble = "WORKFLOW_DOCS_DIR=$workflowdocsdir";
    $postamble .= << 'END';

workflow_install	: dirs shellwrappers promer nucmer blastp ber pe cogs snps legacy2chado chado2bsml bsml2chado jaccard

dirs	:
	mkdir -p $(WORKFLOW_DOCS_DIR)

shellwrappers:
	@cat workflow/run_wf.sh | sed 's?INSTALL_DIR=.*?INSTALL_DIR=$(PREFIX)/bin?' > bin/run_wf.sh
	@install bin/run_wf.sh $(PREFIX)/bin
	@cat workflow/run_wfmonitor.sh | sed 's?INSTALL_DIR=.*?INSTALL_DIR=$(PREFIX)/bin?' > bin/run_wfmonitor.sh
	@install bin/run_wfmonitor.sh $(PREFIX)/bin	
	@install workflow/exec_env.sh $(PREFIX)/bin
	@install workflow/kill_wf.sh $(PREFIX)/bin


promer	:
	@echo "Installing promer workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/promer-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/promer-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/promer.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/promer_template.xml $(WORKFLOW_DOCS_DIR)  

nucmer	:
	@echo "Installing nucmer workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/nucmer-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/nucmer-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/nucmer.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/nucmer_template.xml $(WORKFLOW_DOCS_DIR)  

blastp	:
	@echo "Installing blastp workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/blastp-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/blastp-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/blastp.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/blastp_template.xml $(WORKFLOW_DOCS_DIR)  

ber	:
	@echo "Installing ber workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/ber-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/ber-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/ber.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/ber_template.xml $(WORKFLOW_DOCS_DIR)  

pe	:
	@echo "Installing pe workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/pe-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/pe-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/pe.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/pe_template.xml $(WORKFLOW_DOCS_DIR)  

cogs	:
	@echo "Installing cogs workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/cogs-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/cogs-master_template.xml $(WORKFLOW_DOCS_DIR) 

snps	:
	@echo "Installing snps workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/snps-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/snps-master_template.xml $(WORKFLOW_DOCS_DIR) 

legacy2chado  :
	@echo "Installing legacy2chado workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/legacy2chado-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/legacy2chado-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/legacy2chado.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/legacy2chado_template.xml $(WORKFLOW_DOCS_DIR)  

chado2bsml  :
	@echo "Installing chado2bmsl workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/chado2bsml.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/chado2bsml_template.xml $(WORKFLOW_DOCS_DIR)  

bsml2chado  :
	@echo "Installing bmsl2chado workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/bsml2chado.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/bsml2chado_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/bsml2chado-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/bsml2chado-master_template.xml $(WORKFLOW_DOCS_DIR)  

jaccard  :
	@echo "Installing bmsl2chado workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/jaccard_cluster.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/jaccard_cluster_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/jaccard_cluster-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/jaccard_cluster-master_template.xml $(WORKFLOW_DOCS_DIR)  


END

    $postamble;
}

sub install { 
    my $self = shift;
    my $installstr = $self->MM::install(@_);
    my $addltarget = 'workflow_install';
    $installstr =~ s/(pure_install\s+::\s+)(.*)/$1 $2 $addltarget/;
    return $installstr;
}
