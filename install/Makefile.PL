 use ExtUtils::MakeMaker;
require ExtUtils::Manifest;
require "perl2wrapper_ergatis.pl";
use strict;

## always remake the MANIFEST
ExtUtils::Manifest::mkmanifest();

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

if(!$isworkflowdocsdir || !(-e $workflowdocsdir)){
    die "ERROR: Directory '$workflowdocsdir' does not exist.  You must specify a WORKFLOW_DOCS_DIR.\n".
	"Workflow XML and ini files will be stored in this directory\n".
	"Eg. \n\nperl Makefile.PL PREFIX=/usr/local/devel/ANNOTATION/shared WORKFLOW_DOCS_DIR=/usr/local/devel/ANNOTATION/shared/docs\n\n";
}
else{
    print "Ignore the following warning: ";
}
if(!$isschemadocsdir || (!(-e $schemadocsdir))){
    die "ERROR: Directory '$schemadocsdir' does not exist.  You must specify a SCHEMA_DOCS_DIR\n".
	"BSML dtd and schemas are stored in this directory\n".
        "Eg. \n\nperl Makefile.PL PREFIX=/usr/local/devel/ANNOTATION/shared WORKFLOW_DOCS_DIR=/usr/local/devel/ANNOTATION/shared/docs SCHEMA_DOCS_DIR=/usr/local/devel/ANNOTATION/shared/docs\n\n";
}
else{
    print "Ignore the following warning: ";
}

my($scripts,$wrapperstr) = &get_bins();


WriteMakefile(
	      'NAME'=>'ergatis',
	      'EXE_FILES'=>$scripts,
	      'clean' =>{FILES => "MANIFEST.bak $wrapperstr README.symlinks"}
	      );

package MY;

sub MY::postamble { 
    my $postamble = "WORKFLOW_DOCS_DIR=$workflowdocsdir";
    $postamble .= << 'END';


ergatis_install	: dirs shellwrappers wait

ergatis_install_last	: dirs shellwrappers aat_aa aat_aa_masked aat_na aat_na_masked adjust_bsml_coordinates adjust_btab_coordinates adjust_gap2_coordinates adjust_nap_coordinates augustus backupdb ber bsml2chado bsml2chado_fast bsml2fasta bsml2featurerelationships build_supercontig chado2bsml clustalw cogs copydb coverage delete_analysis fasta2bsml fgenesh genbank2bsml genezilla genemark genewise_best_loc genie genscan genscan-plus gff2bsml glimmer3 glimmerHMM hmmpfam hmmpir hmmsearch hmmsmart hmmtigr infernal initdb iprscan jaccard legacy2bsml legacyparafam2bsml mask_by_analysis mepipe mreps muscle ncbi-blastn ncbi-blastp ncbi-blastx ncbi-tblastn ncbi-tblastx ngap nucmer nucmer_snps obo2chado pangenome panther parse_for_ncRNAs pe pepstats phat predotar promer promote_gene_prediction ps_scan psipred region repeatmasker signalp snps snap snoscan split_fasta split_multifasta statistics2chado subfeature_mapping tabula_rasa targetp tmhmm translate_sequence transterm tRNAscan-SE tiling trf wait wu-blastn wu-blastp wu-blastx wu-tblastn wu-tblastx xdformat

	@install component/shared/iterator_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/shared/groups_template.xml $(WORKFLOW_DOCS_DIR)
        @install component/shared/project.config $(WORKFLOW_DOCS_DIR)

dirs	:
	mkdir -p $(WORKFLOW_DOCS_DIR)

shellwrappers:
	@install component/create_ber_db.sh $(PREFIX)/bin
	@install component/pepstats.sh $(PREFIX)/bin
	@install component/transeq.sh $(PREFIX)/bin

aat_aa	:
	@echo "Installing aat_aa workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/aat_aa/aat_aa-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/aat_aa/aat_aa-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/aat_aa/aat_aa.ini $(WORKFLOW_DOCS_DIR)
	@install component/aat_aa/aat_aa_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/aat_aa/aat_aaconf.ini $(WORKFLOW_DOCS_DIR)

aat_aa_masked	:
	@echo "Installing aat_aa_masked workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/aat_aa_masked/aat_aa_masked-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/aat_aa_masked/aat_aa_masked-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/aat_aa_masked/aat_aa_masked.ini $(WORKFLOW_DOCS_DIR)
	@install component/aat_aa_masked/aat_aa_masked_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/aat_aa_masked/aat_aa_maskedconf.ini $(WORKFLOW_DOCS_DIR)

aat_na	:
	@echo "Installing aat_na workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/aat_na/aat_na-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/aat_na/aat_na-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/aat_na/aat_na.ini $(WORKFLOW_DOCS_DIR)
	@install component/aat_na/aat_na_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/aat_na/aat_naconf.ini $(WORKFLOW_DOCS_DIR)

aat_na_masked	:
	@echo "Installing aat_na_masked workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/aat_na_masked/aat_na_masked-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/aat_na_masked/aat_na_masked-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/aat_na_masked/aat_na_masked.ini $(WORKFLOW_DOCS_DIR)
	@install component/aat_na_masked/aat_na_masked_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/aat_na_masked/aat_na_maskedconf.ini $(WORKFLOW_DOCS_DIR)

adjust_bsml_coordinates	:
	@echo "Installing adjust_bsml_coordinates workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/adjust_bsml_coordinates/adjust_bsml_coordinates-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/adjust_bsml_coordinates/adjust_bsml_coordinates-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/adjust_bsml_coordinates/adjust_bsml_coordinatesconf.ini $(WORKFLOW_DOCS_DIR) 
	@install component/adjust_bsml_coordinates/adjust_bsml_coordinates_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/adjust_bsml_coordinates/adjust_bsml_coordinates.ini $(WORKFLOW_DOCS_DIR) 

adjust_btab_coordinates	:
	@echo "Installing adjust_btab_coordinates workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/adjust_btab_coordinates/adjust_btab_coordinates-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/adjust_btab_coordinates/adjust_btab_coordinates-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/adjust_btab_coordinates/adjust_btab_coordinatesconf.ini $(WORKFLOW_DOCS_DIR) 

adjust_gap2_coordinates	:
	@echo "Installing adjust_gap2_coordinates workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/adjust_gap2_coordinates/adjust_gap2_coordinates-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/adjust_gap2_coordinates/adjust_gap2_coordinates-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/adjust_gap2_coordinates/adjust_gap2_coordinatesconf.ini $(WORKFLOW_DOCS_DIR) 

adjust_nap_coordinates	:
	@echo "Installing adjust_nap_coordinates workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/adjust_nap_coordinates/adjust_nap_coordinates-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/adjust_nap_coordinates/adjust_nap_coordinates-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/adjust_nap_coordinates/adjust_nap_coordinatesconf.ini $(WORKFLOW_DOCS_DIR) 

augustus	:
	@echo "Installing augustus workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/augustus/augustus-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/augustus/augustus-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/augustus/augustus.ini $(WORKFLOW_DOCS_DIR)
	@install component/augustus/augustus_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/augustus/augustusconf.ini $(WORKFLOW_DOCS_DIR)

backupdb  :
	@echo "Installing backupdb workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/backupdb/backupdb-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/backupdb/backupdb-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/backupdb/backupdbconf.ini $(WORKFLOW_DOCS_DIR)  

ber	:
	@echo "Installing ber workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/ber/ber-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/ber/ber-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/ber/ber.ini $(WORKFLOW_DOCS_DIR)
	@install component/ber/ber_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/ber/berconf.ini $(WORKFLOW_DOCS_DIR)  

bsml2chado  :
	@echo "Installing bsml2chado workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/bsml2chado/bsml2chado.ini $(WORKFLOW_DOCS_DIR)
	@install component/bsml2chado/bsml2chado_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/bsml2chado/bsml2chado-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/bsml2chado/bsml2chado-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/bsml2chado/bsml2chadoconf.ini $(WORKFLOW_DOCS_DIR)  

bsml2chado_fast  :
	@echo "Installing bsml2chado_fast workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/bsml2chado_fast/bsml2chado_fast.ini $(WORKFLOW_DOCS_DIR)
	@install component/bsml2chado_fast/bsml2chado_fast_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/bsml2chado_fast/bsml2chado_fast-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/bsml2chado_fast/bsml2chado_fast-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/bsml2chado_fast/bsml2chado_fastconf.ini $(WORKFLOW_DOCS_DIR)  

bsml2fasta	:
	@echo "Installing bsml2fasta workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/bsml2fasta/bsml2fastaconf.ini $(WORKFLOW_DOCS_DIR)
	@install component/bsml2fasta/bsml2fasta-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/bsml2fasta/bsml2fasta-master_template.xml $(WORKFLOW_DOCS_DIR)

bsml2featurerelationships	:
	@echo "Installing bsml2featurerelationships workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/bsml2featurerelationships/bsml2featurerelationshipsconf.ini $(WORKFLOW_DOCS_DIR)
	@install component/bsml2featurerelationships/bsml2featurerelationships-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/bsml2featurerelationships/bsml2featurerelationships-master_template.xml $(WORKFLOW_DOCS_DIR)

build_supercontig	:
	@echo "Installing build_supercontig workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/build_supercontig/build_supercontigconf.ini $(WORKFLOW_DOCS_DIR)
	@install component/build_supercontig/build_supercontig-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/build_supercontig/build_supercontig-master_template.xml $(WORKFLOW_DOCS_DIR)

chado2bsml  :
	@echo "Installing chado2bmsl workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/chado2bsml/chado2bsml-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/chado2bsml/chado2bsml-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/chado2bsml/chado2bsml.ini $(WORKFLOW_DOCS_DIR)
	@install component/chado2bsml/chado2bsml_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/chado2bsml/chado2bsmlconf.ini $(WORKFLOW_DOCS_DIR)  

clustalw	:
	@echo "Installing clustalw workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/clustalw/clustalw-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/clustalw/clustalw-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/clustalw/clustalw.ini $(WORKFLOW_DOCS_DIR)
	@install component/clustalw/clustalw_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/clustalw/clustalwconf.ini $(WORKFLOW_DOCS_DIR)  

cogs	:
	@echo "Installing cogs workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/cogs/cogs-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/cogs/cogs-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/cogs/cogs.ini $(WORKFLOW_DOCS_DIR)
	@install component/cogs/cogs_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/cogs/cogsconf.ini $(WORKFLOW_DOCS_DIR)  

copydb	:
	@echo "Installing copydb workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/copydb/copydb-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/copydb/copydb-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/copydb/copydbconf.ini $(WORKFLOW_DOCS_DIR)  

coverage  :
	@echo "Installing coverage workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/coverage/coverage-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/coverage/coverage-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/coverage/coverageconf.ini $(WORKFLOW_DOCS_DIR) 
	@install component/coverage/coverage.ini $(WORKFLOW_DOCS_DIR)
	@install component/coverage/coverage_template.xml $(WORKFLOW_DOCS_DIR)  

delete_analysis  :
	@echo "Installing delete_analysis workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/delete_analysis/delete_analysis-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/delete_analysis/delete_analysis-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/delete_analysis/delete_analysisconf.ini $(WORKFLOW_DOCS_DIR)  

fasta2bsml	:
	@echo "Installing fasta2bsml workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/fasta2bsml/fasta2bsmlconf.ini $(WORKFLOW_DOCS_DIR)
	@install component/fasta2bsml/fasta2bsml-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/fasta2bsml/fasta2bsml-master_template.xml $(WORKFLOW_DOCS_DIR)

fgenesh  :
	@echo "Installing fgenesh workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/fgenesh/fgenesh-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/fgenesh/fgenesh-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/fgenesh/fgenesh.ini $(WORKFLOW_DOCS_DIR)
	@install component/fgenesh/fgenesh_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/fgenesh/fgeneshconf.ini $(WORKFLOW_DOCS_DIR)  

genbank2bsml  :
	@echo "Installing genbank2bsml workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/genbank2bsml/genbank2bsml-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/genbank2bsml/genbank2bsml-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/genbank2bsml/genbank2bsml.ini $(WORKFLOW_DOCS_DIR)
	@install component/genbank2bsml/genbank2bsml_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/genbank2bsml/genbank2bsmlconf.ini $(WORKFLOW_DOCS_DIR)  

genemark  :
	@echo "Installing genemark workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/genemark/genemark-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/genemark/genemark-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/genemark/genemark.ini $(WORKFLOW_DOCS_DIR)
	@install component/genemark/genemark_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/genemark/genemarkconf.ini $(WORKFLOW_DOCS_DIR)  

genewise_best_loc	:
	@echo "Installing genewise_best_loc workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/genewise_best_loc/genewise_best_loc-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/genewise_best_loc/genewise_best_loc-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/genewise_best_loc/genewise_best_loc.ini $(WORKFLOW_DOCS_DIR)
	@install component/genewise_best_loc/genewise_best_loc_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/genewise_best_loc/genewise_best_locconf.ini $(WORKFLOW_DOCS_DIR)

genezilla  :
	@echo "Installing genezilla workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/genezilla/genezilla-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/genezilla/genezilla-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/genezilla/genezilla.ini $(WORKFLOW_DOCS_DIR)
	@install component/genezilla/genezilla_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/genezilla/genezillaconf.ini $(WORKFLOW_DOCS_DIR)  

genie  :
	@echo "Installing genie workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/genie/genie-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/genie/genie-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/genie/genie.ini $(WORKFLOW_DOCS_DIR)
	@install component/genie/genie_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/genie/genieconf.ini $(WORKFLOW_DOCS_DIR)  

genscan  :
	@echo "Installing genscan workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/genscan/genscan-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/genscan/genscan-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/genscan/genscan.ini $(WORKFLOW_DOCS_DIR)
	@install component/genscan/genscan_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/genscan/genscanconf.ini $(WORKFLOW_DOCS_DIR)  

genscan-plus  :
	@echo "Installing genscan-plus workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/genscan-plus/genscan-plus-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/genscan-plus/genscan-plus-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/genscan-plus/genscan-plus.ini $(WORKFLOW_DOCS_DIR)
	@install component/genscan-plus/genscan-plus_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/genscan-plus/genscan-plusconf.ini $(WORKFLOW_DOCS_DIR)  

gff2bsml	:
	@echo "Installing gff2bsml workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/gff2bsml/gff2bsml-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/gff2bsml/gff2bsml-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/gff2bsml/gff2bsml.ini $(WORKFLOW_DOCS_DIR)
	@install component/gff2bsml/gff2bsml_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/gff2bsml/gff2bsmlconf.ini $(WORKFLOW_DOCS_DIR)

glimmer3	:
	@echo "Installing glimmer3 workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/glimmer3/glimmer3-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/glimmer3/glimmer3-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/glimmer3/glimmer3.ini $(WORKFLOW_DOCS_DIR)
	@install component/glimmer3/glimmer3_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/glimmer3/glimmer3conf.ini $(WORKFLOW_DOCS_DIR)

glimmerHMM	:
	@echo "Installing glimmerHMM workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/glimmerHMM/glimmerHMM-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/glimmerHMM/glimmerHMM-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/glimmerHMM/glimmerHMM.ini $(WORKFLOW_DOCS_DIR)
	@install component/glimmerHMM/glimmerHMM_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/glimmerHMM/glimmerHMMconf.ini $(WORKFLOW_DOCS_DIR)

hmmpfam	:
	@echo "Installing hmmpfam workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/hmmpfam/hmmpfam-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/hmmpfam/hmmpfam-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/hmmpfam/hmmpfam.ini $(WORKFLOW_DOCS_DIR)
	@install component/hmmpfam/hmmpfam_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/hmmpfam/hmmpfamconf.ini $(WORKFLOW_DOCS_DIR)

hmmpir	:
	@echo "Installing hmmpir workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/hmmpir/hmmpir-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/hmmpir/hmmpir-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/hmmpir/hmmpir.ini $(WORKFLOW_DOCS_DIR)
	@install component/hmmpir/hmmpir_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/hmmpir/hmmpirconf.ini $(WORKFLOW_DOCS_DIR)

hmmsearch	:
	@echo "Installing hmmsearch workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/hmmsearch/hmmsearch-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/hmmsearch/hmmsearch-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/hmmsearch/hmmsearch.ini $(WORKFLOW_DOCS_DIR)
	@install component/hmmsearch/hmmsearch_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/hmmsearch/hmmsearchconf.ini $(WORKFLOW_DOCS_DIR)

hmmsmart	:
	@echo "Installing hmmsmart workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/hmmsmart/hmmsmart-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/hmmsmart/hmmsmart-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/hmmsmart/hmmsmart.ini $(WORKFLOW_DOCS_DIR)
	@install component/hmmsmart/hmmsmart_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/hmmsmart/hmmsmartconf.ini $(WORKFLOW_DOCS_DIR)

hmmtigr	:
	@echo "Installing hmmtigr workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/hmmtigr/hmmtigr-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/hmmtigr/hmmtigr-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/hmmtigr/hmmtigr.ini $(WORKFLOW_DOCS_DIR)
	@install component/hmmtigr/hmmtigr_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/hmmtigr/hmmtigrconf.ini $(WORKFLOW_DOCS_DIR)

infernal	:
	@echo "Installing infernal workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/infernal/infernal-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/infernal/infernal-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/infernal/infernal.ini $(WORKFLOW_DOCS_DIR)
	@install component/infernal/infernal_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/infernal/infernalconf.ini $(WORKFLOW_DOCS_DIR)

initdb  :
	@echo "Installing initdb workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/initdb/initdb-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/initdb/initdb-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/initdb/initdbconf.ini $(WORKFLOW_DOCS_DIR)  

iprscan	:
	@echo "Installing iprscan workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/iprscan/iprscan-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/iprscan/iprscan-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/iprscan/iprscan.ini $(WORKFLOW_DOCS_DIR)
	@install component/iprscan/iprscan_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/iprscan/iprscanconf.ini $(WORKFLOW_DOCS_DIR)

jaccard  :
	@echo "Installing jaccard workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/jaccard/jaccard.ini $(WORKFLOW_DOCS_DIR)
	@install component/jaccard/jaccard_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/jaccard/jaccard-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/jaccard/jaccard-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/jaccard/jaccardconf.ini $(WORKFLOW_DOCS_DIR)  

legacy2bsml  :
	@echo "Installing legacy2bsml workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/legacy2bsml/legacy2bsml.ini $(WORKFLOW_DOCS_DIR)
	@install component/legacy2bsml/legacy2bsml_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/legacy2bsml/legacy2bsml-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/legacy2bsml/legacy2bsml-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/legacy2bsml/legacy2bsmlconf.ini $(WORKFLOW_DOCS_DIR)  

legacyparafam2bsml  :
	@echo "Installing legacyparafam2bsml workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/legacyparafam2bsml/legacyparafam2bsml.ini $(WORKFLOW_DOCS_DIR)
	@install component/legacyparafam2bsml/legacyparafam2bsml_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/legacyparafam2bsml/legacyparafam2bsml-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/legacyparafam2bsml/legacyparafam2bsml-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/legacyparafam2bsml/legacyparafam2bsmlconf.ini $(WORKFLOW_DOCS_DIR)  

mask_by_analysis  :
	@echo "Installing mask_by_analysis workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/mask_by_analysis/mask_by_analysis.ini $(WORKFLOW_DOCS_DIR)
	@install component/mask_by_analysis/mask_by_analysis_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/mask_by_analysis/mask_by_analysis-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/mask_by_analysis/mask_by_analysis-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/mask_by_analysis/mask_by_analysisconf.ini $(WORKFLOW_DOCS_DIR)  

mepipe  :
	@echo "Installing mepipe workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/mepipe/mepipe.ini $(WORKFLOW_DOCS_DIR)
	@install component/mepipe/mepipe_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/mepipe/mepipe-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/mepipe/mepipe-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/mepipe/mepipeconf.ini $(WORKFLOW_DOCS_DIR)  

mreps  :
	@echo "Installing mreps workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/mreps/mreps.ini $(WORKFLOW_DOCS_DIR)
	@install component/mreps/mreps_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/mreps/mreps-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/mreps/mreps-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/mreps/mrepsconf.ini $(WORKFLOW_DOCS_DIR)  


muscle  :
	@echo "Installing muscle workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/muscle/muscle.ini $(WORKFLOW_DOCS_DIR)
	@install component/muscle/muscle_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/muscle/muscle-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/muscle/muscle-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/muscle/muscleconf.ini $(WORKFLOW_DOCS_DIR)  

ncbi-blastn	:
	@echo "Installing ncbi-blastn workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/ncbi-blastn/ncbi-blastn-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/ncbi-blastn/ncbi-blastn-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/ncbi-blastn/ncbi-blastn.ini $(WORKFLOW_DOCS_DIR)
	@install component/ncbi-blastn/ncbi-blastn_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/ncbi-blastn/ncbi-blastnconf.ini $(WORKFLOW_DOCS_DIR)

ncbi-blastp	:
	@echo "Installing ncbi-blastp workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/ncbi-blastp/ncbi-blastp-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/ncbi-blastp/ncbi-blastp-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/ncbi-blastp/ncbi-blastp.ini $(WORKFLOW_DOCS_DIR)
	@install component/ncbi-blastp/ncbi-blastp_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/ncbi-blastp/ncbi-blastpconf.ini $(WORKFLOW_DOCS_DIR)

ncbi-blastx	:
	@echo "Installing ncbi-blastx workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/ncbi-blastx/ncbi-blastx-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/ncbi-blastx/ncbi-blastx-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/ncbi-blastx/ncbi-blastx.ini $(WORKFLOW_DOCS_DIR)
	@install component/ncbi-blastx/ncbi-blastx_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/ncbi-blastx/ncbi-blastxconf.ini $(WORKFLOW_DOCS_DIR)

ncbi-tblastn	:
	@echo "Installing ncbi-tblastn workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/ncbi-tblastn/ncbi-tblastn-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/ncbi-tblastn/ncbi-tblastn-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/ncbi-tblastn/ncbi-tblastn.ini $(WORKFLOW_DOCS_DIR)
	@install component/ncbi-tblastn/ncbi-tblastn_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/ncbi-tblastn/ncbi-tblastnconf.ini $(WORKFLOW_DOCS_DIR)

ncbi-tblastx	:
	@echo "Installing ncbi-tblastx workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/ncbi-tblastx/ncbi-tblastx-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/ncbi-tblastx/ncbi-tblastx-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/ncbi-tblastx/ncbi-tblastx.ini $(WORKFLOW_DOCS_DIR)
	@install component/ncbi-tblastx/ncbi-tblastx_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/ncbi-tblastx/ncbi-tblastxconf.ini $(WORKFLOW_DOCS_DIR)

ngap	:
	@echo "Installing ngap workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/ngap/ngap-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/ngap/ngap-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/ngap/ngap.ini $(WORKFLOW_DOCS_DIR)
	@install component/ngap/ngap_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/ngap/ngapconf.ini $(WORKFLOW_DOCS_DIR)  

nucmer	:
	@echo "Installing nucmer workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/nucmer/nucmer-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/nucmer/nucmer-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/nucmer/nucmer.ini $(WORKFLOW_DOCS_DIR)
	@install component/nucmer/nucmer_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/nucmer/nucmerconf.ini $(WORKFLOW_DOCS_DIR)  

nucmer_snps	:
	@echo "Installing nucmer_snps workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/nucmer_snps/nucmer_snps-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/nucmer_snps/nucmer_snps-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/nucmer_snps/nucmer_snps.ini $(WORKFLOW_DOCS_DIR)
	@install component/nucmer_snps/nucmer_snps_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/nucmer_snps/nucmer_snpsconf.ini $(WORKFLOW_DOCS_DIR)  

obo2chado	:
	@echo "Installing obo2chado workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/obo2chado/obo2chado-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/obo2chado/obo2chado-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/obo2chado/obo2chado.ini $(WORKFLOW_DOCS_DIR)
	@install component/obo2chado/obo2chado_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/obo2chado/obo2chadoconf.ini $(WORKFLOW_DOCS_DIR)  
pangenome	:
	@echo "Installing pangenome workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/pangenome/pangenome-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/pangenome/pangenome-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/pangenome/pangenome.ini $(WORKFLOW_DOCS_DIR)
	@install component/pangenome/pangenome_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/pangenome/pangenomeconf.ini $(WORKFLOW_DOCS_DIR)  
	
panther	:
	@echo "Installing panther workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/panther/panther-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/panther/panther-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/panther/panther.ini $(WORKFLOW_DOCS_DIR)
	@install component/panther/panther_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/panther/pantherconf.ini $(WORKFLOW_DOCS_DIR)  

parse_for_ncRNAs	:
	@echo "Installing parse_for_ncRNAs workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/parse_for_ncRNAs/parse_for_ncRNA-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/parse_for_ncRNAs/parse_for_ncRNA-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/parse_for_ncRNAs/parse_for_ncRNA.ini $(WORKFLOW_DOCS_DIR)
	@install component/parse_for_ncRNAs/parse_for_ncRNA_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/parse_for_ncRNAs/parse_for_ncRNAconf.ini $(WORKFLOW_DOCS_DIR)  

pe	:
	@echo "Installing pe workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/pe/pe-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/pe/pe-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/pe/pe.ini $(WORKFLOW_DOCS_DIR)
	@install component/pe/pe_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/pe/peconf.ini $(WORKFLOW_DOCS_DIR)  

pepstats	:
	@echo "Installing pepstats workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/pepstats/pepstats-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/pepstats/pepstats-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/pepstats/pepstats.ini $(WORKFLOW_DOCS_DIR)
	@install component/pepstats/pepstats_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/pepstats/pepstatsconf.ini $(WORKFLOW_DOCS_DIR)  

phat	:
	@echo "Installing phat workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/phat/phat-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/phat/phat-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/phat/phat.ini $(WORKFLOW_DOCS_DIR)
	@install component/phat/phat_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/phat/phatconf.ini $(WORKFLOW_DOCS_DIR)  

promer	:
	@echo "Installing promer workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/promer/promer-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/promer/promer-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/promer/promer.ini $(WORKFLOW_DOCS_DIR)
	@install component/promer/promer_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/promer/promerconf.ini $(WORKFLOW_DOCS_DIR)  

predotar	:
	@echo "Installing predotar workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/predotar/predotar-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/predotar/predotar-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/predotar/predotar.ini $(WORKFLOW_DOCS_DIR)
	@install component/predotar/predotar_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/predotar/predotarconf.ini $(WORKFLOW_DOCS_DIR)  

promote_gene_prediction	:
	@echo "Installing promote_gene_prediction workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/promote_gene_prediction/promote_gene_prediction-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/promote_gene_prediction/promote_gene_prediction-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/promote_gene_prediction/promote_gene_prediction.ini $(WORKFLOW_DOCS_DIR)
	@install component/promote_gene_prediction/promote_gene_prediction_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/promote_gene_prediction/promote_gene_predictionconf.ini $(WORKFLOW_DOCS_DIR)  

ps_scan	:
	@echo "Installing ps_scan workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/ps_scan/ps_scan-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/ps_scan/ps_scan-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/ps_scan/ps_scan.ini $(WORKFLOW_DOCS_DIR)
	@install component/ps_scan/ps_scan_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/ps_scan/ps_scanconf.ini $(WORKFLOW_DOCS_DIR)  

psipred	:
	@echo "Installing psipred workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/psipred/psipred-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/psipred/psipred-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/psipred/psipred.ini $(WORKFLOW_DOCS_DIR)
	@install component/psipred/psipred_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/psipred/psipredconf.ini $(WORKFLOW_DOCS_DIR)  

region	:
	@echo "Installing region workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/region/region-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/region/region-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/region/region.ini $(WORKFLOW_DOCS_DIR)
	@install component/region/region_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/region/regionconf.ini $(WORKFLOW_DOCS_DIR)  

repeatmasker :
	@echo "Installing repeatmasker workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/repeatmasker/repeatmasker-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/repeatmasker/repeatmasker-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/repeatmasker/repeatmasker.ini $(WORKFLOW_DOCS_DIR)
	@install component/repeatmasker/repeatmasker_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/repeatmasker/repeatmaskerconf.ini $(WORKFLOW_DOCS_DIR)

signalp	:
	@echo "Installing signalp workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/signalp/signalp-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/signalp/signalp-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/signalp/signalp.ini $(WORKFLOW_DOCS_DIR)
	@install component/signalp/signalp_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/signalp/signalpconf.ini $(WORKFLOW_DOCS_DIR)  

snoscan	:
	@echo "Installing snoscan workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/snoscan/snoscan-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/snoscan/snoscan-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/snoscan/snoscan.ini $(WORKFLOW_DOCS_DIR)
	@install component/snoscan/snoscan_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/snoscan/snoscanconf.ini $(WORKFLOW_DOCS_DIR) 

snps	:
	@echo "Installing snps workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/snp/snp-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/snp/snp-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/snp/snpconf.ini $(WORKFLOW_DOCS_DIR)

snap	:
	@echo "Installing snap workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/snap/snap-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/snap/snap-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/snap/snap.ini $(WORKFLOW_DOCS_DIR)
	@install component/snap/snap_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/snap/snapconf.ini $(WORKFLOW_DOCS_DIR)

snoscan	:
	@echo "Installing snoscan workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/snoscan/snoscan-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/snoscan/snoscan-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/snoscan/snoscan.ini $(WORKFLOW_DOCS_DIR)
	@install component/snoscan/snoscan_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/snoscan/snoscanconf.ini $(WORKFLOW_DOCS_DIR)

split_fasta	:
	@echo "Installing split_fasta workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/split_fasta/split_fasta-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/split_fasta/split_fasta-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/split_fasta/split_fasta.ini $(WORKFLOW_DOCS_DIR)
	@install component/split_fasta/split_fasta_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/split_fasta/split_fastaconf.ini $(WORKFLOW_DOCS_DIR)

split_multifasta	:
	@echo "Installing split_multifasta workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/split_multifasta/split_multifasta-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/split_multifasta/split_multifasta-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/split_multifasta/split_multifasta.ini $(WORKFLOW_DOCS_DIR)
	@install component/split_multifasta/split_multifasta_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/split_multifasta/split_multifastaconf.ini $(WORKFLOW_DOCS_DIR)

statistics2chado	:
	@echo "Installing statistics2chado workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/statistics2chado/statistics2chado-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/statistics2chado/statistics2chado-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/statistics2chado/statistics2chadoconf.ini $(WORKFLOW_DOCS_DIR)

subfeature_mapping  :
	@echo "Installing subfeature_mapping workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/subfeature_mapping/subfeature_mapping-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/subfeature_mapping/subfeature_mapping-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/subfeature_mapping/subfeature_mappingconf.ini $(WORKFLOW_DOCS_DIR)  

tabula_rasa	:
	@echo "Installing tabula_rasa workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/tabula_rasa/tabula_rasa-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/tabula_rasa/tabula_rasa-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/tabula_rasa/tabula_rasa.ini $(WORKFLOW_DOCS_DIR)
	@install component/tabula_rasa/tabula_rasa_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/tabula_rasa/tabula_rasaconf.ini $(WORKFLOW_DOCS_DIR)  

targetp	:
	@echo "Installing targetp workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/targetp/targetp-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/targetp/targetp-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/targetp/targetp.ini $(WORKFLOW_DOCS_DIR)
	@install component/targetp/targetp_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/targetp/targetpconf.ini $(WORKFLOW_DOCS_DIR)  

translate_sequence	:
	@echo "Installing translate_sequence workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/translate_sequence/translate_sequence-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/translate_sequence/translate_sequence-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/translate_sequence/translate_sequence.ini $(WORKFLOW_DOCS_DIR)
	@install component/translate_sequence/translate_sequence_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/translate_sequence/translate_sequenceconf.ini $(WORKFLOW_DOCS_DIR) 

transterm	:
	@echo "Installing transterm workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/transterm/transterm-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/transterm/transterm-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/transterm/transterm.ini $(WORKFLOW_DOCS_DIR)
	@install component/transterm/transterm_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/transterm/transtermconf.ini $(WORKFLOW_DOCS_DIR) 
	
tmhmm	:
	@echo "Installing tmhmm workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/tmhmm/tmhmm-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/tmhmm/tmhmm-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/tmhmm/tmhmm.ini $(WORKFLOW_DOCS_DIR)
	@install component/tmhmm/tmhmm_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/tmhmm/tmhmmconf.ini $(WORKFLOW_DOCS_DIR)  
	
tRNAscan-SE	:
	@echo "Installing tRNAscan-SE workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/tRNAscan-SE/tRNAscan-SE-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/tRNAscan-SE/tRNAscan-SE-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/tRNAscan-SE/tRNAscan-SE.ini $(WORKFLOW_DOCS_DIR)
	@install component/tRNAscan-SE/tRNAscan-SE_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/tRNAscan-SE/tRNAscan-SEconf.ini $(WORKFLOW_DOCS_DIR)

tiling  :
	@echo "Installing tiling workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/tiling/tiling-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/tiling/tiling-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/tiling/tilingconf.ini $(WORKFLOW_DOCS_DIR) 

trf	:
	@echo "Installing trf workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/trf/trf-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/trf/trf-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/trf/trf.ini $(WORKFLOW_DOCS_DIR)
	@install component/trf/trf_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/trf/trfconf.ini $(WORKFLOW_DOCS_DIR)

wait	:
	@echo "Installing wait workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/wait/wait.xml $(WORKFLOW_DOCS_DIR)
	@install component/wait/wait.i1.xml $(WORKFLOW_DOCS_DIR)
	@install component/wait/wait.config $(WORKFLOW_DOCS_DIR)

wu-blastn	:
	@echo "Installing wu-blastn workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/wu-blastn/wu-blastn-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/wu-blastn/wu-blastn-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/wu-blastn/wu-blastn.ini $(WORKFLOW_DOCS_DIR)
	@install component/wu-blastn/wu-blastn_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/wu-blastn/wu-blastnconf.ini $(WORKFLOW_DOCS_DIR)

wu-blastp	:
	@echo "Installing wu-blastp workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/wu-blastp/wu-blastp-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/wu-blastp/wu-blastp-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/wu-blastp/wu-blastp.ini $(WORKFLOW_DOCS_DIR)
	@install component/wu-blastp/wu-blastp_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/wu-blastp/wu-blastpconf.ini $(WORKFLOW_DOCS_DIR)

wu-blastx	:
	@echo "Installing wu-blastx workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/wu-blastx/wu-blastx-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/wu-blastx/wu-blastx-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/wu-blastx/wu-blastx.ini $(WORKFLOW_DOCS_DIR)
	@install component/wu-blastx/wu-blastx_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/wu-blastx/wu-blastxconf.ini $(WORKFLOW_DOCS_DIR)

wu-tblastn	:
	@echo "Installing wu-tblastn workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/wu-tblastn/wu-tblastn-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/wu-tblastn/wu-tblastn-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/wu-tblastn/wu-tblastn.ini $(WORKFLOW_DOCS_DIR)
	@install component/wu-tblastn/wu-tblastn_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/wu-tblastn/wu-tblastnconf.ini $(WORKFLOW_DOCS_DIR)

wu-tblastx	:
	@echo "Installing wu-tblastx workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/wu-tblastx/wu-tblastx-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/wu-tblastx/wu-tblastx-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install component/wu-tblastx/wu-tblastx.ini $(WORKFLOW_DOCS_DIR)
	@install component/wu-tblastx/wu-tblastx_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/wu-tblastx/wu-tblastxconf.ini $(WORKFLOW_DOCS_DIR)

xdformat	:
	@echo "Installing xdformat workflows in $(WORKFLOW_DOCS_DIR)"
	@install component/xdformat/xdformat-master.ini $(WORKFLOW_DOCS_DIR)
	@install component/xdformat/xdformat-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install component/xdformat/xdformatconf.ini $(WORKFLOW_DOCS_DIR) 




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
