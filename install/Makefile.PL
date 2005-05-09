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


workflow_install	: dirs shellwrappers aat_aa aat_na adjust_bsml_coordinates adjust_btab_coordinates adjust_gap2_coordinates adjust_nap_coordinates backupdb ber bsml2chado bsml2fasta chado2bsml cogs coverage delete_analysis genie glimmerHMM initdb jaccard legacy2bsml nucmer pe promer region repeatmasker snps split_fasta split_multifasta tRNAscan-SE tiling trf wait wu_blastp xdformat

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

aat_aa	:
	@echo "Installing aat_aa workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/aat_aa/aat_aa-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/aat_aa/aat_aa-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/aat_aa/aat_aa.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/aat_aa/aat_aa_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/aat_aa/aat_aaconf.ini $(WORKFLOW_DOCS_DIR)

aat_na	:
	@echo "Installing aat_na workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/aat_na/aat_na-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/aat_na/aat_na-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/aat_na/aat_na.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/aat_na/aat_na_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/aat_na/aat_naconf.ini $(WORKFLOW_DOCS_DIR)

adjust_bsml_coordinates	:
	@echo "Installing adjust_bsml_coordinates workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/adjust_bsml_coordinates/adjust_bsml_coordinates-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/adjust_bsml_coordinates/adjust_bsml_coordinates-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/adjust_bsml_coordinates/adjust_bsml_coordinatesconf.ini $(WORKFLOW_DOCS_DIR) 

adjust_btab_coordinates	:
	@echo "Installing adjust_btab_coordinates workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/adjust_btab_coordinates/adjust_btab_coordinates-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/adjust_btab_coordinates/adjust_btab_coordinates-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/adjust_btab_coordinates/adjust_btab_coordinatesconf.ini $(WORKFLOW_DOCS_DIR) 

adjust_gap2_coordinates	:
	@echo "Installing adjust_gap2_coordinates workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/adjust_gap2_coordinates/adjust_gap2_coordinates-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/adjust_gap2_coordinates/adjust_gap2_coordinates-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/adjust_gap2_coordinates/adjust_gap2_coordinatesconf.ini $(WORKFLOW_DOCS_DIR) 

adjust_nap_coordinates	:
	@echo "Installing adjust_nap_coordinates workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/adjust_nap_coordinates/adjust_nap_coordinates-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/adjust_nap_coordinates/adjust_nap_coordinates-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/adjust_nap_coordinates/adjust_nap_coordinatesconf.ini $(WORKFLOW_DOCS_DIR) 

backupdb  :
	@echo "Installing backupdb workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/backupdb/backupdb-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/backupdb/backupdb-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/backupdb/backupdbconf.ini $(WORKFLOW_DOCS_DIR)  

ber	:
	@echo "Installing ber workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/ber/ber-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/ber/ber-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/ber/ber.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/ber/ber_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/ber/berconf.ini $(WORKFLOW_DOCS_DIR)  

bsml2chado  :
	@echo "Installing bmsl2chado workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/bsml2chado/bsml2chado.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/bsml2chado/bsml2chado_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/bsml2chado/bsml2chado-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/bsml2chado/bsml2chado-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/bsml2chado/bsml2chadoconf.ini $(WORKFLOW_DOCS_DIR)  

bsml2fasta	:
	@echo "Installing bsml2fasta workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/bsml2fasta/bsml2fastaconf.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/bsml2fasta/bsml2fasta-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/bsml2fasta/bsml2fasta-master_template.xml $(WORKFLOW_DOCS_DIR)

chado2bsml  :
	@echo "Installing chado2bmsl workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/chado2bsml/chado2bsml-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/chado2bsml/chado2bsml-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/chado2bsml/chado2bsml.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/chado2bsml/chado2bsml_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/chado2bsml/chado2bsmlconf.ini $(WORKFLOW_DOCS_DIR)  

cogs	:
	@echo "Installing cogs workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/cogs/cogs-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/cogs/cogs-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/cogs/cogs.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/cogs/cogs_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/cogs/cogsconf.ini $(WORKFLOW_DOCS_DIR)  

coverage  :
	@echo "Installing coverage workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/coverage/coverage-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/coverage/coverage-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/coverage/coverageconf.ini $(WORKFLOW_DOCS_DIR) 
	@install workflow/coverage/coverage.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/coverage/coverage_template.xml $(WORKFLOW_DOCS_DIR)  

delete_analysis  :
	@echo "Installing delete_analysis workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/delete_analysis/delete_analysis-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/delete_analysis/delete_analysis-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/delete_analysis/delete_analysisconf.ini $(WORKFLOW_DOCS_DIR)  

genie  :
	@echo "Installing genie workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/genie/genie-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/genie/genie-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/genie/genie.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/genie/genie_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/genie/genieconf.ini $(WORKFLOW_DOCS_DIR)  

glimmerHMM	:
	@echo "Installing glimmerHMM workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/glimmerHMM/glimmerHMM-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/glimmerHMM/glimmerHMM-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/glimmerHMM/glimmerHMM.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/glimmerHMM/glimmerHMM_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/glimmerHMM/glimmerHMMconf.ini $(WORKFLOW_DOCS_DIR)

initdb  :
	@echo "Installing initdb workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/initdb/initdb-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/initdb/initdb-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/initdb/initdbconf.ini $(WORKFLOW_DOCS_DIR)  

jaccard  :
	@echo "Installing jaccard workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/jaccard/jaccard.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/jaccard/jaccard_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/jaccard/jaccard-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/jaccard/jaccard-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/jaccard/jaccardconf.ini $(WORKFLOW_DOCS_DIR)  

legacy2bsml  :
	@echo "Installing legacy2bmsl workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/legacy2bsml/legacy2bsml.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/legacy2bsml/legacy2bsml_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/legacy2bsml/legacy2bsml-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/legacy2bsml/legacy2bsml-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/legacy2bsml/legacy2bsmlconf.ini $(WORKFLOW_DOCS_DIR)  

nucmer	:
	@echo "Installing nucmer workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/nucmer/nucmer-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/nucmer/nucmer-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/nucmer/nucmer.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/nucmer/nucmer_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/nucmer/nucmerconf.ini $(WORKFLOW_DOCS_DIR)  

pe	:
	@echo "Installing pe workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/pe/pe-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/pe/pe-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/pe/pe.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/pe/pe_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/pe/peconf.ini $(WORKFLOW_DOCS_DIR)  

promer	:
	@echo "Installing promer workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/promer/promer-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/promer/promer-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/promer/promer.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/promer/promer_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/promer/promerconf.ini $(WORKFLOW_DOCS_DIR)  

region	:
	@echo "Installing region workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/region/region-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/region/region-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/region/region.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/region/region_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/region/regionconf.ini $(WORKFLOW_DOCS_DIR)  

repeatmasker :
	@echo "Installing repeatmasker workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/repeatmasker/repeatmasker-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/repeatmasker/repeatmasker-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/repeatmasker/repeatmasker.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/repeatmasker/repeatmasker_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/repeatmasker/repeatmaskerconf.ini $(WORKFLOW_DOCS_DIR)

snps	:
	@echo "Installing snps workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/snp/snp-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/snp/snp-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/snp/snpconf.ini $(WORKFLOW_DOCS_DIR)
split_fasta	:
	@echo "Installing split_fasta workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/split_fasta/split_fasta-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/split_fasta/split_fasta-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/split_fasta/split_fasta.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/split_fasta/split_fasta_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/split_fasta/split_fastaconf.ini $(WORKFLOW_DOCS_DIR)

split_multifasta	:
	@echo "Installing split_multifasta workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/split_multifasta/split_multifasta-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/split_multifasta/split_multifasta-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/split_multifasta/split_multifasta.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/split_multifasta/split_multifasta_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/split_multifasta/split_multifastaconf.ini $(WORKFLOW_DOCS_DIR)

tRNAscan-SE	:
	@echo "Installing tRNAscan-SE workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/tRNAscan-SE/tRNAscan-SE-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/tRNAscan-SE/tRNAscan-SE-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/tRNAscan-SE/tRNAscan-SE.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/tRNAscan-SE/tRNAscan-SE_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/tRNAscan-SE/tRNAscan-SEconf.ini $(WORKFLOW_DOCS_DIR)

tiling  :
	@echo "Installing tiling workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/tiling/tiling-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/tiling/tiling-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/tiling/tilingconf.ini $(WORKFLOW_DOCS_DIR) 

trf	:
	@echo "Installing trf workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/trf/trf-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/trf/trf-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/trf/trf.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/trf/trf_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/trf/trfconf.ini $(WORKFLOW_DOCS_DIR)

wait	:
	@echo "Installing wait workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/wait/wait-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/wait/wait-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/wait/wait.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/wait/wait_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/wait/waitconf.ini $(WORKFLOW_DOCS_DIR)

wu-blastp	:
	@echo "Installing wu-blastp workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/wu-blastp/wu-blastp-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/wu-blastp/wu-blastp-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install workflow/wu-blastp/wu-blastp.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/wu-blastp/wu-blastp_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/wu-blastp/wu-blastpconf.ini $(WORKFLOW_DOCS_DIR)


xdformat	:
	@echo "Installing xdformat workflows in $(WORKFLOW_DOCS_DIR)"
	@install workflow/xdformat/xdformat-master.ini $(WORKFLOW_DOCS_DIR)
	@install workflow/xdformat/xdformat-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install workflow/xdformat/xdformatconf.ini $(WORKFLOW_DOCS_DIR) 




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
