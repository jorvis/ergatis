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


ergatis_install : dirs shellwrappers aat_aa aat_aa_masked aat_na aat_na_masked bsml2chado bsml2fasta clustalw fasta2bsml fgenesh formatdb genemark genewise glimmer3 glimmerHMM hmmpfam infernal initdb iprscan j_ortholog_clusters jaccard legacy2bsml ncbi-blastp ngap pangenome parse_for_ncRNA phat repeatmasker runsql signalp split_multifasta targetp tmhmm train_for_glimmer3_iteration trf tRNAscan-SE wait wu-blastp wu-blastn wu-tblastn xdformat

	@install components/shared/iterator_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/shared/file_iterator_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/shared/groups_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/shared/project.config.example $(WORKFLOW_DOCS_DIR)
	@install components/shared/software.config $(PREFIX)

ergatis_install_last : dirs shellwrappers aat_aa aat_aa_masked aat_na aat_na_masked adjust_bsml_coordinates adjust_btab_coordinates adjust_gap2_coordinates adjust_nap_coordinates augustus backupdb ber bsml2chado bsml2chado_fast bsml2fasta bsml2featurerelationships build_supercontig chado2bsml clustalw cogs copydb coverage delete_analysis fasta2bsml fgenesh formatdb genbank2bsml genezilla genemark genie genscan genscan-plus gff2bsml glimmer3 glimmerHMM hmmpfam hmmpir hmmsearch hmmsmart hmmtigr infernal initdb iprscan jaccard legacy2bsml legacyparafam2bsml mask_by_analysis mepipe mreps muscle ncbi-blastn ncbi-blastp ncbi-blastx ncbi-tblastn ncbi-tblastx ngap nucmer nucmer_snps obo2chado pangenome panther parse_for_ncRNAs pe pepstats phat predotar promer promote_gene_prediction ps_scan psipred region repeatmasker signalp snps snap snoscan split_fasta split_multifasta statistics2chado subfeature_mapping tabula_rasa targetp tmhmm train_for_glimmer3_iteration translate_sequence transterm tRNAscan-SE tiling trf wait wu-blastn wu-blastp wu-blastx wu-tblastn wu-tblastx xdformat

dirs	:
	mkdir -p $(WORKFLOW_DOCS_DIR)
	mkdir -p $(WORKFLOW_DOCS_DIR)/documentation

shellwrappers:
	@install bin/create_ber_db.sh $(PREFIX)/bin
	@install bin/pepstats.sh $(PREFIX)/bin
	@install bin/transeq.sh $(PREFIX)/bin

aat_aa	:
	@echo "Installing aat_aa workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/aat_aa/aat_aa.xml $(WORKFLOW_DOCS_DIR)
	@install components/aat_aa/aat_aa.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/aat_aa/aat_aa.config $(WORKFLOW_DOCS_DIR)

aat_aa_masked	:
	@echo "Installing aat_aa_masked workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/aat_aa_masked/aat_aa_masked.xml $(WORKFLOW_DOCS_DIR)
	@install components/aat_aa_masked/aat_aa_masked.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/aat_aa_masked/aat_aa_masked.config $(WORKFLOW_DOCS_DIR)

aat_na	:
	@echo "Installing aat_na workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/aat_na/aat_na.xml $(WORKFLOW_DOCS_DIR)
	@install components/aat_na/aat_na.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/aat_na/aat_na.config $(WORKFLOW_DOCS_DIR)

aat_na_masked	:
	@echo "Installing aat_na_masked workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/aat_na_masked/aat_na_masked.xml $(WORKFLOW_DOCS_DIR)
	@install components/aat_na_masked/aat_na_masked.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/aat_na_masked/aat_na_masked.config $(WORKFLOW_DOCS_DIR)

adjust_bsml_coordinates	:
	@echo "Installing adjust_bsml_coordinates workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/adjust_bsml_coordinates/adjust_bsml_coordinates-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/adjust_bsml_coordinates/adjust_bsml_coordinates-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/adjust_bsml_coordinates/adjust_bsml_coordinatesconf.ini $(WORKFLOW_DOCS_DIR) 
	@install components/adjust_bsml_coordinates/adjust_bsml_coordinates_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/adjust_bsml_coordinates/adjust_bsml_coordinates.ini $(WORKFLOW_DOCS_DIR) 

adjust_btab_coordinates	:
	@echo "Installing adjust_btab_coordinates workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/adjust_btab_coordinates/adjust_btab_coordinates-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/adjust_btab_coordinates/adjust_btab_coordinates-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/adjust_btab_coordinates/adjust_btab_coordinatesconf.ini $(WORKFLOW_DOCS_DIR) 

adjust_gap2_coordinates	:
	@echo "Installing adjust_gap2_coordinates workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/adjust_gap2_coordinates/adjust_gap2_coordinates-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/adjust_gap2_coordinates/adjust_gap2_coordinates-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/adjust_gap2_coordinates/adjust_gap2_coordinatesconf.ini $(WORKFLOW_DOCS_DIR) 

adjust_nap_coordinates	:
	@echo "Installing adjust_nap_coordinates workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/adjust_nap_coordinates/adjust_nap_coordinates-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/adjust_nap_coordinates/adjust_nap_coordinates-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/adjust_nap_coordinates/adjust_nap_coordinatesconf.ini $(WORKFLOW_DOCS_DIR) 

augustus	:
	@echo "Installing augustus workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/augustus/augustus-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/augustus/augustus-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/augustus/augustus.ini $(WORKFLOW_DOCS_DIR)
	@install components/augustus/augustus_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/augustus/augustusconf.ini $(WORKFLOW_DOCS_DIR)

backupdb  :
	@echo "Installing backupdb workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/backupdb/backupdb-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/backupdb/backupdb-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/backupdb/backupdbconf.ini $(WORKFLOW_DOCS_DIR)  

ber	:
	@echo "Installing ber workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/ber/ber-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/ber/ber-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/ber/ber.ini $(WORKFLOW_DOCS_DIR)
	@install components/ber/ber_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/ber/berconf.ini $(WORKFLOW_DOCS_DIR)  

bsml2chado  :
	@echo "Installing bsml2chado workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/bsml2chado/bsml2chado.config $(WORKFLOW_DOCS_DIR)
	@install components/bsml2chado/bsml2chado.xml $(WORKFLOW_DOCS_DIR)  
	@install components/bsml2chado/bsml2chado.parse.xml $(WORKFLOW_DOCS_DIR)
	@install components/bsml2chado/bsml2chado.sql.xml $(WORKFLOW_DOCS_DIR)  
	@install components/bsml2chado/bsml2chado.sqlforce.xml $(WORKFLOW_DOCS_DIR)  
	@install components/bsml2chado/bsml2chado.load.xml $(WORKFLOW_DOCS_DIR)  
	@install components/bsml2chado/bsml2chado.collate.xml $(WORKFLOW_DOCS_DIR)  

runsql  :
	@echo "Installing runsql workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/runsql/runsql.config $(WORKFLOW_DOCS_DIR)
	@install components/runsql/runsql.xml $(WORKFLOW_DOCS_DIR)  
	@install components/runsql/runsql.sql.xml $(WORKFLOW_DOCS_DIR)  

bsml2chado_fast  :
	@echo "Installing bsml2chado_fast workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/bsml2chado_fast/bsml2chado_fast.ini $(WORKFLOW_DOCS_DIR)
	@install components/bsml2chado_fast/bsml2chado_fast_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/bsml2chado_fast/bsml2chado_fast-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/bsml2chado_fast/bsml2chado_fast-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/bsml2chado_fast/bsml2chado_fastconf.ini $(WORKFLOW_DOCS_DIR)  

bsml2fasta	:
	@echo "Installing bsml2fasta workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/bsml2fasta/bsml2fasta.xml $(WORKFLOW_DOCS_DIR)
	@install components/bsml2fasta/bsml2fasta.config $(WORKFLOW_DOCS_DIR)

bsml2featurerelationships	:
	@echo "Installing bsml2featurerelationships workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/bsml2featurerelationships/bsml2featurerelationshipsconf.ini $(WORKFLOW_DOCS_DIR)
	@install components/bsml2featurerelationships/bsml2featurerelationships-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/bsml2featurerelationships/bsml2featurerelationships-master_template.xml $(WORKFLOW_DOCS_DIR)

build_supercontig	:
	@echo "Installing build_supercontig workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/build_supercontig/build_supercontigconf.ini $(WORKFLOW_DOCS_DIR)
	@install components/build_supercontig/build_supercontig-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/build_supercontig/build_supercontig-master_template.xml $(WORKFLOW_DOCS_DIR)

chado2bsml  :
	@echo "Installing chado2bmsl workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/chado2bsml/chado2bsml-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/chado2bsml/chado2bsml-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/chado2bsml/chado2bsml.ini $(WORKFLOW_DOCS_DIR)
	@install components/chado2bsml/chado2bsml_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/chado2bsml/chado2bsmlconf.ini $(WORKFLOW_DOCS_DIR)  

clustalw	:
	@echo "Installing clustalw workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/clustalw/clustalw.xml $(WORKFLOW_DOCS_DIR)
	@install components/clustalw/clustalw.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/clustalw/clustalw.config $(WORKFLOW_DOCS_DIR)

cogs	:
	@echo "Installing cogs workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/cogs/cogs-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/cogs/cogs-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/cogs/cogs.ini $(WORKFLOW_DOCS_DIR)
	@install components/cogs/cogs_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/cogs/cogsconf.ini $(WORKFLOW_DOCS_DIR)  

copydb	:
	@echo "Installing copydb workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/copydb/copydb-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/copydb/copydb-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/copydb/copydbconf.ini $(WORKFLOW_DOCS_DIR)  

coverage  :
	@echo "Installing coverage workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/coverage/coverage-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/coverage/coverage-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/coverage/coverageconf.ini $(WORKFLOW_DOCS_DIR) 
	@install components/coverage/coverage.ini $(WORKFLOW_DOCS_DIR)
	@install components/coverage/coverage_template.xml $(WORKFLOW_DOCS_DIR)  

delete_analysis  :
	@echo "Installing delete_analysis workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/delete_analysis/delete_analysis-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/delete_analysis/delete_analysis-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/delete_analysis/delete_analysisconf.ini $(WORKFLOW_DOCS_DIR)  

fasta2bsml	:
	@echo "Installing fasta2bsml workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/fasta2bsml/fasta2bsml.config $(WORKFLOW_DOCS_DIR)
	@install components/fasta2bsml/fasta2bsml.xml $(WORKFLOW_DOCS_DIR)

fgenesh	:
	@echo "Installing fgenesh workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/fgenesh/fgenesh.xml $(WORKFLOW_DOCS_DIR)
	@install components/fgenesh/fgenesh.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/fgenesh/fgenesh.config $(WORKFLOW_DOCS_DIR) 

formatdb	:
	@echo "Installing formatdb workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/formatdb/formatdb.xml $(WORKFLOW_DOCS_DIR)
	@install components/formatdb/formatdb.config $(WORKFLOW_DOCS_DIR)
	@install components/formatdb/formatdb.tmpl $(WORKFLOW_DOCS_DIR)/documentation

genbank2bsml  :
	@echo "Installing genbank2bsml workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/genbank2bsml/genbank2bsml-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/genbank2bsml/genbank2bsml-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/genbank2bsml/genbank2bsml.ini $(WORKFLOW_DOCS_DIR)
	@install components/genbank2bsml/genbank2bsml_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/genbank2bsml/genbank2bsmlconf.ini $(WORKFLOW_DOCS_DIR)  

genemark	:
	@echo "Installing genemark workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/genemark/genemark.xml $(WORKFLOW_DOCS_DIR)
	@install components/genemark/genemark.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/genemark/genemark.config $(WORKFLOW_DOCS_DIR)

genewise	:
	@echo "Installing genewise workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/genewise/genewise.xml $(WORKFLOW_DOCS_DIR)
	@install components/genewise/genewise.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/genewise/genewise.config $(WORKFLOW_DOCS_DIR)

genezilla  :
	@echo "Installing genezilla workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/genezilla/genezilla-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/genezilla/genezilla-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/genezilla/genezilla.ini $(WORKFLOW_DOCS_DIR)
	@install components/genezilla/genezilla_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/genezilla/genezillaconf.ini $(WORKFLOW_DOCS_DIR)  

genie  :
	@echo "Installing genie workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/genie/genie-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/genie/genie-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/genie/genie.ini $(WORKFLOW_DOCS_DIR)
	@install components/genie/genie_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/genie/genieconf.ini $(WORKFLOW_DOCS_DIR)  

genscan  :
	@echo "Installing genscan workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/genscan/genscan-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/genscan/genscan-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/genscan/genscan.ini $(WORKFLOW_DOCS_DIR)
	@install components/genscan/genscan_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/genscan/genscanconf.ini $(WORKFLOW_DOCS_DIR)  

genscan-plus  :
	@echo "Installing genscan-plus workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/genscan-plus/genscan-plus-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/genscan-plus/genscan-plus-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/genscan-plus/genscan-plus.ini $(WORKFLOW_DOCS_DIR)
	@install components/genscan-plus/genscan-plus_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/genscan-plus/genscan-plusconf.ini $(WORKFLOW_DOCS_DIR)  

gff2bsml	:
	@echo "Installing gff2bsml workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/gff2bsml/gff2bsml-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/gff2bsml/gff2bsml-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/gff2bsml/gff2bsml.ini $(WORKFLOW_DOCS_DIR)
	@install components/gff2bsml/gff2bsml_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/gff2bsml/gff2bsmlconf.ini $(WORKFLOW_DOCS_DIR)

glimmer3	:
	@echo "Installing glimmer3 workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/glimmer3/glimmer3.config $(WORKFLOW_DOCS_DIR)
	@install components/glimmer3/glimmer3.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/glimmer3/glimmer3.xml $(WORKFLOW_DOCS_DIR)

glimmerHMM	:
	@echo "Installing glimmerHMM workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/glimmerHMM/glimmerHMM.xml $(WORKFLOW_DOCS_DIR)
	@install components/glimmerHMM/glimmerHMM.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/glimmerHMM/glimmerHMM.config $(WORKFLOW_DOCS_DIR)

hmmpfam	:
	@echo "Installing hmmpfam workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/hmmpfam/hmmpfam.xml $(WORKFLOW_DOCS_DIR)
	@install components/hmmpfam/hmmpfam.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/hmmpfam/hmmpfam.config $(WORKFLOW_DOCS_DIR)

hmmpir	:
	@echo "Installing hmmpir workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/hmmpir/hmmpir-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/hmmpir/hmmpir-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/hmmpir/hmmpir.ini $(WORKFLOW_DOCS_DIR)
	@install components/hmmpir/hmmpir_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/hmmpir/hmmpirconf.ini $(WORKFLOW_DOCS_DIR)

hmmsearch	:
	@echo "Installing hmmsearch workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/hmmsearch/hmmsearch-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/hmmsearch/hmmsearch-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/hmmsearch/hmmsearch.ini $(WORKFLOW_DOCS_DIR)
	@install components/hmmsearch/hmmsearch_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/hmmsearch/hmmsearchconf.ini $(WORKFLOW_DOCS_DIR)

hmmsmart	:
	@echo "Installing hmmsmart workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/hmmsmart/hmmsmart-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/hmmsmart/hmmsmart-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/hmmsmart/hmmsmart.ini $(WORKFLOW_DOCS_DIR)
	@install components/hmmsmart/hmmsmart_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/hmmsmart/hmmsmartconf.ini $(WORKFLOW_DOCS_DIR)

hmmtigr	:
	@echo "Installing hmmtigr workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/hmmtigr/hmmtigr-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/hmmtigr/hmmtigr-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/hmmtigr/hmmtigr.ini $(WORKFLOW_DOCS_DIR)
	@install components/hmmtigr/hmmtigr_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/hmmtigr/hmmtigrconf.ini $(WORKFLOW_DOCS_DIR)

infernal	:
	@echo "Installing infernal workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/infernal/infernal.config $(WORKFLOW_DOCS_DIR)
	@install components/infernal/infernal.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/infernal/infernal.xml $(WORKFLOW_DOCS_DIR)

initdb	:
	@echo "Installing initdb workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/initdb/initdb.xml $(WORKFLOW_DOCS_DIR)
	@install components/initdb/initdb.sql.xml $(WORKFLOW_DOCS_DIR)
	@install components/initdb/initdb.sqlforce.xml $(WORKFLOW_DOCS_DIR)
	@install components/initdb/initdb.obo.xml $(WORKFLOW_DOCS_DIR)
	@install components/initdb/initdb.config $(WORKFLOW_DOCS_DIR)

iprscan	:
	@echo "Installing iprscan workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/iprscan/iprscan.xml $(WORKFLOW_DOCS_DIR)
	@install components/iprscan/iprscan.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/iprscan/iprscan.config $(WORKFLOW_DOCS_DIR)
	@install components/iprscan/iprscan.tmpl $(WORKFLOW_DOCS_DIR)/documentation

j_ortholog_clusters	:
	@echo "Installing j_ortholog_clusters workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/j_ortholog_clusters/j_ortholog_clusters.xml $(WORKFLOW_DOCS_DIR)
	@install components/j_ortholog_clusters/j_ortholog_clusters.config $(WORKFLOW_DOCS_DIR)

jaccard	:
	@echo "Installing jaccard workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/jaccard/jaccard.xml $(WORKFLOW_DOCS_DIR)
	@install components/jaccard/jaccard.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/jaccard/jaccard.config $(WORKFLOW_DOCS_DIR)

legacy2bsml	:
	@echo "Installing legacy2bsml workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/legacy2bsml/legacy2bsml.xml $(WORKFLOW_DOCS_DIR)
	@install components/legacy2bsml/legacy2bsml.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/legacy2bsml/legacy2bsml.config $(WORKFLOW_DOCS_DIR)
    
legacyparafam2bsml  :
	@echo "Installing legacyparafam2bsml workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/legacyparafam2bsml/legacyparafam2bsml.ini $(WORKFLOW_DOCS_DIR)
	@install components/legacyparafam2bsml/legacyparafam2bsml_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/legacyparafam2bsml/legacyparafam2bsml-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/legacyparafam2bsml/legacyparafam2bsml-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/legacyparafam2bsml/legacyparafam2bsmlconf.ini $(WORKFLOW_DOCS_DIR)  

lipoprotein_motif  :
	@echo "Installing lipoprotein_motif workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/lipoprotein_motif/lipoprotein_motif.ini $(WORKFLOW_DOCS_DIR)
	@install components/lipoprotein_motif/lipoprotein_motif_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/lipoprotein_motif/lipoprotein_motif-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/lipoprotein_motif/lipoprotein_motif-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/lipoprotein_motif/lipoprotein_motifconf.ini $(WORKFLOW_DOCS_DIR)  

mask_by_analysis  :
	@echo "Installing mask_by_analysis workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/mask_by_analysis/mask_by_analysis.ini $(WORKFLOW_DOCS_DIR)
	@install components/mask_by_analysis/mask_by_analysis_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/mask_by_analysis/mask_by_analysis-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/mask_by_analysis/mask_by_analysis-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/mask_by_analysis/mask_by_analysisconf.ini $(WORKFLOW_DOCS_DIR)  

mepipe  :
	@echo "Installing mepipe workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/mepipe/mepipe.ini $(WORKFLOW_DOCS_DIR)
	@install components/mepipe/mepipe_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/mepipe/mepipe-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/mepipe/mepipe-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/mepipe/mepipeconf.ini $(WORKFLOW_DOCS_DIR)  

mreps  :
	@echo "Installing mreps workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/mreps/mreps.ini $(WORKFLOW_DOCS_DIR)
	@install components/mreps/mreps_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/mreps/mreps-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/mreps/mreps-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/mreps/mrepsconf.ini $(WORKFLOW_DOCS_DIR)  


muscle  :
	@echo "Installing muscle workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/muscle/muscle.ini $(WORKFLOW_DOCS_DIR)
	@install components/muscle/muscle_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/muscle/muscle-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/muscle/muscle-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/muscle/muscleconf.ini $(WORKFLOW_DOCS_DIR)  

ncbi-blastn	:
	@echo "Installing ncbi-blastn workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/ncbi-blastn/ncbi-blastn-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/ncbi-blastn/ncbi-blastn-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/ncbi-blastn/ncbi-blastn.ini $(WORKFLOW_DOCS_DIR)
	@install components/ncbi-blastn/ncbi-blastn_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/ncbi-blastn/ncbi-blastnconf.ini $(WORKFLOW_DOCS_DIR)

ncbi-blastp	:
	@echo "Installing ncbi-blastp workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/ncbi-blastp/ncbi-blastp.xml $(WORKFLOW_DOCS_DIR)
	@install components/ncbi-blastp/ncbi-blastp.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/ncbi-blastp/ncbi-blastp.config $(WORKFLOW_DOCS_DIR)

ncbi-blastx	:
	@echo "Installing ncbi-blastx workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/ncbi-blastx/ncbi-blastx-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/ncbi-blastx/ncbi-blastx-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/ncbi-blastx/ncbi-blastx.ini $(WORKFLOW_DOCS_DIR)
	@install components/ncbi-blastx/ncbi-blastx_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/ncbi-blastx/ncbi-blastxconf.ini $(WORKFLOW_DOCS_DIR)

ncbi-tblastn	:
	@echo "Installing ncbi-tblastn workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/ncbi-tblastn/ncbi-tblastn-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/ncbi-tblastn/ncbi-tblastn-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/ncbi-tblastn/ncbi-tblastn.ini $(WORKFLOW_DOCS_DIR)
	@install components/ncbi-tblastn/ncbi-tblastn_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/ncbi-tblastn/ncbi-tblastnconf.ini $(WORKFLOW_DOCS_DIR)

ncbi-tblastx	:
	@echo "Installing ncbi-tblastx workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/ncbi-tblastx/ncbi-tblastx-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/ncbi-tblastx/ncbi-tblastx-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/ncbi-tblastx/ncbi-tblastx.ini $(WORKFLOW_DOCS_DIR)
	@install components/ncbi-tblastx/ncbi-tblastx_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/ncbi-tblastx/ncbi-tblastxconf.ini $(WORKFLOW_DOCS_DIR)

ngap	:
	@echo "Installing ngap workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/ngap/ngap.xml $(WORKFLOW_DOCS_DIR)
	@install components/ngap/ngap.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/ngap/ngap.config $(WORKFLOW_DOCS_DIR)  

nucmer	:
	@echo "Installing nucmer workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/nucmer/nucmer-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/nucmer/nucmer-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/nucmer/nucmer.ini $(WORKFLOW_DOCS_DIR)
	@install components/nucmer/nucmer_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/nucmer/nucmerconf.ini $(WORKFLOW_DOCS_DIR)  

nucmer_snps	:
	@echo "Installing nucmer_snps workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/nucmer_snps/nucmer_snps-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/nucmer_snps/nucmer_snps-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/nucmer_snps/nucmer_snps.ini $(WORKFLOW_DOCS_DIR)
	@install components/nucmer_snps/nucmer_snps_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/nucmer_snps/nucmer_snpsconf.ini $(WORKFLOW_DOCS_DIR)  

obo2chado	:
	@echo "Installing obo2chado workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/obo2chado/obo2chado-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/obo2chado/obo2chado-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/obo2chado/obo2chado.ini $(WORKFLOW_DOCS_DIR)
	@install components/obo2chado/obo2chado_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/obo2chado/obo2chadoconf.ini $(WORKFLOW_DOCS_DIR)  

pangenome	:
	@echo "Installing pangenome workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/pangenome/pangenome.xml $(WORKFLOW_DOCS_DIR)
	@install components/pangenome/pangenome.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/pangenome/pangenome.config $(WORKFLOW_DOCS_DIR)
	@install components/pangenome/pangenome.tmpl $(WORKFLOW_DOCS_DIR)/documentation
    
panther	:
	@echo "Installing panther workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/panther/panther.xml $(WORKFLOW_DOCS_DIR)
	@install components/panther/panther.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/panther/panther.config $(WORKFLOW_DOCS_DIR) 

parse_for_ncRNAs	:
	@echo "Installing parse_for_ncRNAs workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/parse_for_ncRNAs/parse_for_ncRNA.config $(WORKFLOW_DOCS_DIR)
	@install components/parse_for_ncRNAs/parse_for_ncRNA.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/parse_for_ncRNAs/parse_for_ncRNA.xml $(WORKFLOW_DOCS_DIR)  

pe	:
	@echo "Installing pe workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/pe/pe-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/pe/pe-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/pe/pe.ini $(WORKFLOW_DOCS_DIR)
	@install components/pe/pe_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/pe/peconf.ini $(WORKFLOW_DOCS_DIR)  

pepstats	:
	@echo "Installing pepstats workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/pepstats/pepstats.xml $(WORKFLOW_DOCS_DIR)
	@install components/pepstats/pepstats.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/pepstats/pepstats.config $(WORKFLOW_DOCS_DIR)

phat	:
	@echo "Installing phat workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/phat/phat.xml $(WORKFLOW_DOCS_DIR)
	@install components/phat/phat.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/phat/phat.config $(WORKFLOW_DOCS_DIR) 

predotar	:
	@echo "Installing predotar workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/predotar/predotar.xml $(WORKFLOW_DOCS_DIR)
	@install components/predotar/predotar.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/predotar/predotar.config $(WORKFLOW_DOCS_DIR)

promer	:
	@echo "Installing promer workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/promer/promer-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/promer/promer-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/promer/promer.ini $(WORKFLOW_DOCS_DIR)
	@install components/promer/promer_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/promer/promerconf.ini $(WORKFLOW_DOCS_DIR)  

promote_gene_prediction	:
	@echo "Installing promote_gene_prediction workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/promote_gene_prediction/promote_gene_prediction-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/promote_gene_prediction/promote_gene_prediction-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/promote_gene_prediction/promote_gene_prediction.ini $(WORKFLOW_DOCS_DIR)
	@install components/promote_gene_prediction/promote_gene_prediction_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/promote_gene_prediction/promote_gene_predictionconf.ini $(WORKFLOW_DOCS_DIR)  

ps_scan	:
	@echo "Installing ps_scan workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/ps_scan/ps_scan-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/ps_scan/ps_scan-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/ps_scan/ps_scan.ini $(WORKFLOW_DOCS_DIR)
	@install components/ps_scan/ps_scan_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/ps_scan/ps_scanconf.ini $(WORKFLOW_DOCS_DIR)  

psipred	:
	@echo "Installing psipred workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/psipred/psipred-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/psipred/psipred-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/psipred/psipred.ini $(WORKFLOW_DOCS_DIR)
	@install components/psipred/psipred_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/psipred/psipredconf.ini $(WORKFLOW_DOCS_DIR)  

region	:
	@echo "Installing region workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/region/region-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/region/region-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/region/region.ini $(WORKFLOW_DOCS_DIR)
	@install components/region/region_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/region/regionconf.ini $(WORKFLOW_DOCS_DIR)  

repeatmasker	:
	@echo "Installing repeatmasker workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/repeatmasker/repeatmasker.xml $(WORKFLOW_DOCS_DIR)
	@install components/repeatmasker/repeatmasker.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/repeatmasker/repeatmasker.config $(WORKFLOW_DOCS_DIR)

signalp	:
	@echo "Installing signalp workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/signalp/signalp.xml $(WORKFLOW_DOCS_DIR)
	@install components/signalp/signalp.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/signalp/signalp.config $(WORKFLOW_DOCS_DIR)
	@install components/signalp/signalp.tmpl $(WORKFLOW_DOCS_DIR)/documentation

snap	:
	@echo "Installing snap workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/snap/snap-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/snap/snap-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/snap/snap.ini $(WORKFLOW_DOCS_DIR)
	@install components/snap/snap_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/snap/snapconf.ini $(WORKFLOW_DOCS_DIR)

snoscan	:
	@echo "Installing snoscan workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/snoscan/snoscan-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/snoscan/snoscan-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/snoscan/snoscan.ini $(WORKFLOW_DOCS_DIR)
	@install components/snoscan/snoscan_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/snoscan/snoscanconf.ini $(WORKFLOW_DOCS_DIR) 

snps	:
	@echo "Installing snps workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/snp/snp-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/snp/snp-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/snp/snpconf.ini $(WORKFLOW_DOCS_DIR)

split_fasta	:
	@echo "Installing split_fasta workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/split_fasta/split_fasta.xml $(WORKFLOW_DOCS_DIR)
	@install components/split_fasta/split_fasta.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/split_fasta/split_fasta.config $(WORKFLOW_DOCS_DIR)

split_multifasta	:
	@echo "Installing split_multifasta workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/split_multifasta/split_multifasta.xml $(WORKFLOW_DOCS_DIR)
	@install components/split_multifasta/split_multifasta.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/split_multifasta/split_multifasta.config $(WORKFLOW_DOCS_DIR)

statistics2chado	:
	@echo "Installing statistics2chado workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/statistics2chado/statistics2chado-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/statistics2chado/statistics2chado-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/statistics2chado/statistics2chadoconf.ini $(WORKFLOW_DOCS_DIR)

subfeature_mapping  :
	@echo "Installing subfeature_mapping workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/subfeature_mapping/subfeature_mapping-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/subfeature_mapping/subfeature_mapping-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/subfeature_mapping/subfeature_mappingconf.ini $(WORKFLOW_DOCS_DIR)  

tabula_rasa	:
	@echo "Installing tabula_rasa workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/tabula_rasa/tabula_rasa-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/tabula_rasa/tabula_rasa-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/tabula_rasa/tabula_rasa.ini $(WORKFLOW_DOCS_DIR)
	@install components/tabula_rasa/tabula_rasa_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/tabula_rasa/tabula_rasaconf.ini $(WORKFLOW_DOCS_DIR)  

targetp	:
	@echo "Installing targetp workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/targetp/targetp.xml $(WORKFLOW_DOCS_DIR)
	@install components/targetp/targetp.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/targetp/targetp.config $(WORKFLOW_DOCS_DIR) 

tmhmm	:
	@echo "Installing tmhmm workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/tmhmm/tmhmm.xml $(WORKFLOW_DOCS_DIR)
	@install components/tmhmm/tmhmm.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/tmhmm/tmhmm.config $(WORKFLOW_DOCS_DIR)

train_for_glimmer3_iteration	:
	@echo "Installing train_for_glimmer3_iteration workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/train_for_glimmer3_iteration/train_for_glimmer3_iteration.xml $(WORKFLOW_DOCS_DIR)
	@install components/train_for_glimmer3_iteration/train_for_glimmer3_iteration.config $(WORKFLOW_DOCS_DIR)

translate_sequence	:
	@echo "Installing translate_sequence workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/translate_sequence/translate_sequence-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/translate_sequence/translate_sequence-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/translate_sequence/translate_sequence.ini $(WORKFLOW_DOCS_DIR)
	@install components/translate_sequence/translate_sequence_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/translate_sequence/translate_sequenceconf.ini $(WORKFLOW_DOCS_DIR) 

transterm	:
	@echo "Installing transterm workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/transterm/transterm-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/transterm/transterm-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/transterm/transterm.ini $(WORKFLOW_DOCS_DIR)
	@install components/transterm/transterm_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/transterm/transtermconf.ini $(WORKFLOW_DOCS_DIR) 

tiling  :
	@echo "Installing tiling workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/tiling/tiling-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/tiling/tiling-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/tiling/tilingconf.ini $(WORKFLOW_DOCS_DIR) 

trf	:
	@echo "Installing trf workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/trf/trf.xml $(WORKFLOW_DOCS_DIR)
	@install components/trf/trf.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/trf/trf.config $(WORKFLOW_DOCS_DIR)

tRNAscan-SE	:
	@echo "Installing tRNAscan-SE workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/tRNAscan-SE/tRNAscan-SE.xml $(WORKFLOW_DOCS_DIR)
	@install components/tRNAscan-SE/tRNAscan-SE.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/tRNAscan-SE/tRNAscan-SE.config $(WORKFLOW_DOCS_DIR)

wait	:
	@echo "Installing wait workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/wait/wait.xml $(WORKFLOW_DOCS_DIR)
	@install components/wait/wait.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/wait/wait.i2.xml $(WORKFLOW_DOCS_DIR)
	@install components/wait/wait.config $(WORKFLOW_DOCS_DIR)

wu-blastn	:
	@echo "Installing wu-blastn workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/wu-blastn/wu-blastn.config $(WORKFLOW_DOCS_DIR)
	@install components/wu-blastn/wu-blastn.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/wu-blastn/wu-blastn.xml $(WORKFLOW_DOCS_DIR)

wu-blastp	:
	@echo "Installing wu-blastp workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/wu-blastp/wu-blastp.xml $(WORKFLOW_DOCS_DIR)
	@install components/wu-blastp/wu-blastp.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/wu-blastp/wu-blastp.config $(WORKFLOW_DOCS_DIR)
	@install components/wu-blastp/wu-blastp.tmpl $(WORKFLOW_DOCS_DIR)/documentation

wu-blastx	:
	@echo "Installing wu-blastx workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/wu-blastx/wu-blastx-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/wu-blastx/wu-blastx-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/wu-blastx/wu-blastx.ini $(WORKFLOW_DOCS_DIR)
	@install components/wu-blastx/wu-blastx_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/wu-blastx/wu-blastxconf.ini $(WORKFLOW_DOCS_DIR)

wu-tblastn	:
	@echo "Installing wu-tblastn workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/wu-tblastn/wu-tblastn.xml $(WORKFLOW_DOCS_DIR)
	@install components/wu-tblastn/wu-tblastn.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/wu-tblastn/wu-tblastn.config $(WORKFLOW_DOCS_DIR)
	@install components/wu-tblastn/wu-tblastn.tmpl $(WORKFLOW_DOCS_DIR)/documentation

wu-tblastx	:
	@echo "Installing wu-tblastx workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/wu-tblastx/wu-tblastx-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/wu-tblastx/wu-tblastx-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/wu-tblastx/wu-tblastx.ini $(WORKFLOW_DOCS_DIR)
	@install components/wu-tblastx/wu-tblastx_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/wu-tblastx/wu-tblastxconf.ini $(WORKFLOW_DOCS_DIR)

xdformat	:
	@echo "Installing xdformat workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/xdformat/xdformat.xml $(WORKFLOW_DOCS_DIR)
	@install components/xdformat/xdformat.config $(WORKFLOW_DOCS_DIR)
	@install components/xdformat/xdformat.tmpl $(WORKFLOW_DOCS_DIR)/documentation




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
