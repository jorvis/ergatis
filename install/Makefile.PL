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

workflow_install	: dirs shellwrappers promer nucmer blastp wait pe chado2bsml cogs region jaccard bsml2chado chunker dechunker delete_analysis initdb backupdb snps
	@install workflow/shared/default-iterator.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/shared/default-serialiterator_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/shared/default-paralleliterator_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/shared/batch-iterator.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/shared/groups-iterator_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/shared/batch-paralleliterator_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/shared/batch-serialiterator_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/shared/groups-iterator.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/shared/sharedconf.ini $(WORKFLOW_DOCS_DIR)

dirs	:
	mkdir -p $(WORKFLOW_DOCS_DIR)

shellwrappers:
	@install workflow/CreateWorkflow.sh $(PREFIX)/bin
	@install workflow/RunWorkflow.sh $(PREFIX)/bin
	@install workflow/run_wfmonitor.sh $(PREFIX)/bin

promer	:
	@echo "Installing promer workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/promer/promer-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/promer/promer-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/promer/promer.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/promer/promer_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/promer/promerconf.ini $(WORKFLOW_DOCS_DIR)  

wait	:
	@echo "Installing wait workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/wait/wait-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/wait/wait-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/wait/wait.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/wait/wait_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/wait/waitconf.ini $(WORKFLOW_DOCS_DIR)

nucmer	:
	@echo "Installing nucmer workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/nucmer/nucmer-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/nucmer/nucmer-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/nucmer/nucmer.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/nucmer/nucmer_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/nucmer/nucmerconf.ini $(WORKFLOW_DOCS_DIR)  

blastp	:
	@echo "Installing blastp workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/blastp/blastp-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/blastp/blastp-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/blastp/blastp.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/blastp/blastp_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/blastp/blastpconf.ini $(WORKFLOW_DOCS_DIR)  

ber	:
	@echo "Installing ber workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/ber/ber-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/ber/ber-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/ber/ber.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/ber/ber_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/ber/berconf.ini $(WORKFLOW_DOCS_DIR)  

pe	:
	@echo "Installing pe workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/pe/pe-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/pe/pe-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/pe/pe.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/pe/pe_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/pe/peconf.ini $(WORKFLOW_DOCS_DIR)  

region	:
	@echo "Installing region workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/region/region-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/region/region-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/region/region.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/region/region_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/region/regionconf.ini $(WORKFLOW_DOCS_DIR)  

cogs	:
	@echo "Installing cogs workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/cogs/cogs-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/cogs/cogs-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/cogs/cogs.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/cogs/cogs_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/cogs/cogsconf.ini $(WORKFLOW_DOCS_DIR)  

snps	:
	@echo "Installing snps workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/snp/snp-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/snp/snp-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/snp/snpconf.ini $(WORKFLOW_DOCS_DIR)

legacy2chado  :
	@echo "Installing legacy2chado workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/legacy2chado-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/legacy2chado-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/legacy2chado.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/legacy2chado_template.xml $(WORKFLOW_DOCS_DIR)  

chado2bsml  :
	@echo "Installing chado2bmsl workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/chado2bsml/chado2bsml-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/chado2bsml/chado2bsml-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/chado2bsml/chado2bsml.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/chado2bsml/chado2bsml_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/chado2bsml/chado2bsmlconf.ini $(WORKFLOW_DOCS_DIR)  

bsml2chado  :
	@echo "Installing bmsl2chado workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/bsml2chado/bsml2chado.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/bsml2chado/bsml2chado_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/bsml2chado/bsml2chado-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/bsml2chado/bsml2chado-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/bsml2chado/bsml2chadoconf.ini $(WORKFLOW_DOCS_DIR)  

jaccard  :
	@echo "Installing jaccard workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/jaccard/jaccard.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/jaccard/jaccard_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/jaccard/jaccard-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/jaccard/jaccard-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/jaccard/jaccardconf.ini $(WORKFLOW_DOCS_DIR)  

chunker  :
	@echo "Installing chunker workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/chunker/chunker.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/chunker/chunker_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/chunker/chunker-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/chunker/chunker-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/chunker/chunkerconf.ini $(WORKFLOW_DOCS_DIR)  

dechunker  :
	@echo "Installing dechunker workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/dechunker/dechunker.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/dechunker/dechunker_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/dechunker/dechunker-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/dechunker/dechunker-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/dechunker/dechunkerconf.ini $(WORKFLOW_DOCS_DIR)  

delete_analysis  :
	@echo "Installing delete_analysis workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/delete_analysis/delete_analysis-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/delete_analysis/delete_analysis-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/delete_analysis/delete_analysisconf.ini $(WORKFLOW_DOCS_DIR)  

initdb  :
	@echo "Installing initdb workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/initdb/initdb-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/initdb/initdb-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/initdb/initdbconf.ini $(WORKFLOW_DOCS_DIR)  

backupdb  :
	@echo "Installing backupdb workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/backupdb/backupdb-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/backupdb/backupdb-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/backupdb/backupdbconf.ini $(WORKFLOW_DOCS_DIR)  


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
