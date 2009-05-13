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

if ( scalar @modules_not_found ) {
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

foreach my $arg (@ARGV){
    if ($arg =~ /INSTALL_BASE/) {
	($install_base) = ($arg =~ /INSTALL_BASE=(\S+)/);
    } elsif ( $arg =~ /PREFIX/ ) {
        print "\nERROR: please use INSTALL_BASE instead of PREFIX.  (This requires a current version of ExtUtils::MakeMaker";
        exit(1);
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

my $mm = WriteMakefile(
	      'NAME'=>'ergatis',
	      'EXE_FILES'=>$scripts,
	      'clean' =>{FILES => "MANIFEST.bak $wrapperstr README.symlinks"}
	      );

package MY;

sub MY::postamble { 
    my $postamble = "WORKFLOW_DOCS_DIR=$workflowdocsdir\nSAMPLE_DOCS_DIR=$sampledocsdir";
    $postamble .= << 'END';

ergatis_install : dirs samples shellwrappers c_apps pipeline_templates 454_mplex_seq_bin_trim aat_aa aat_aa_masked aat_na aat_na_masked adjust_bsml_coordinates aragorn augustus auto_gene_curation backupdb ber blast2go bsml2chado bsml2gff3 bsml2fasta bsml2featurerelationships bsml2interevidence_fasta bsml2tbl clustalw combine_ncRNA_alignments cvtermpath_closure fasta2bsml fgenesh formatdb genbank2bsml geneid genemark genewise genezilla genscan genscan-plus glimmer3 glimmerHMM hmmpfam infernal initdb intergenic_analysis iprscan j_ortholog_clusters jaccard legacy2bsml lipop lipoprotein_motif load_ber_for_manatee mask_by_analysis metagene mummer mumsgraph muscle ncbi-blastp ncbi-blastn ncbi-blastn_compress ncbi-blastx ngap nucmer_snps obo2chado open_reading_frames outer_membrane_domain overlap_analysis p_func pangenome_blast_parse pangenome_do_analysis pangenome_make_graphs panther parse_evidence parse_for_ncRNAs pepstats percent_gc phage_finder phat phylip_dnadist pipeline_summary predict_prokaryotic_gene_function predotar priam_ec_assignment promote_gene_prediction ps_scan psiblast quicktree rdp_classifier refresh_chado_mart repeatmasker RNAmmer rpsblast runsql seg signalp snap split_fasta split_multifasta stap_aligner stap_classifier start_site_curation targetp tbl2asn tmhmm train_for_glimmer3_iteration translate_sequence trf tRNAscan-SE tRNAscan-SE_compress wait wu-blastp wu-blastn wu-blastx wu-tblastn xdformat

	@install components/shared/iterator_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/shared/file_iterator_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/shared/groups_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/shared/project.config.example $(WORKFLOW_DOCS_DIR)
	@install components/shared/software.config $(PREFIX)

ergatis_install_last : dirs shellwrappers aat_aa aat_aa_masked aat_na aat_na_masked adjust_bsml_coordinates adjust_btab_coordinates adjust_gap2_coordinates adjust_nap_coordinates augustus auto_gene_curation backupdb ber bsml2chado bsml2chado_fast bsml2fasta bsml2featurerelationships build_supercontig chado2bsml clustalw copydb coverage delete_analysis fasta2bsml fgenesh formatdb genbank2bsml genezilla genemark genie genscan genscan-plus gff2bsml glimmer3 glimmerHMM hmmpfam hmmpir hmmsearch hmmsmart hmmtigr infernal initdb iprscan jaccard legacy2bsml legacyparafam2bsml mask_by_analysis mepipe mreps muscle ncbi-blastn ncbi-blastp ncbi-blastx ncbi-tblastn ncbi-tblastx ngap nucmer nucmer_snps obo2chado pangenome parse_for_ncRNAs pe pepstats phat pipeline_summary predotar promer promote_gene_prediction ps_scan psipred rdp_classifier repeatmasker signalp snps snap snoscan split_fasta split_multifasta stap_aligner stap_classifier statistics2chado subfeature_mapping tabula_rasa targetp tmhmm train_for_glimmer3_iteration translate_sequence transterm tRNAscan-SE tiling trf wait wu-blastn wu-blastp wu-blastx wu-tblastn wu-tblastx xdformat

dirs	:
	mkdir -p $(WORKFLOW_DOCS_DIR)
	mkdir -p $(WORKFLOW_DOCS_DIR)/documentation
	mkdir -p $(SAMPLE_DOCS_DIR)

samples	: dirs
	@install samples/proteins.faa $(SAMPLE_DOCS_DIR)
	@install samples/NC_002695.fna $(SAMPLE_DOCS_DIR)
	@install samples/NC_010364.fna $(SAMPLE_DOCS_DIR)

shellwrappers:
	@install shell/create_ber_db.sh $(PREFIX)/bin
	@install shell/pepstats.sh $(PREFIX)/bin

pipeline_templates:
	@echo "Installing global pipeline templates in $(PREFIX)/global_pipeline_templates"
	mkdir -p $(PREFIX)/global_pipeline_templates/Integrated_JOCs-pangenome_pipeline
	@install global_pipeline_templates/Integrated_JOCs-pangenome_pipeline/bsml2fasta.fastamulti.config $(PREFIX)/global_pipeline_templates/Integrated_JOCs-pangenome_pipeline
	@install global_pipeline_templates/Integrated_JOCs-pangenome_pipeline/bsml2fasta.genomicseq.config $(PREFIX)/global_pipeline_templates/Integrated_JOCs-pangenome_pipeline
	@install global_pipeline_templates/Integrated_JOCs-pangenome_pipeline/jaccard.default.config $(PREFIX)/global_pipeline_templates/Integrated_JOCs-pangenome_pipeline
	@install global_pipeline_templates/Integrated_JOCs-pangenome_pipeline/j_ortholog_clusters.default.config $(PREFIX)/global_pipeline_templates/Integrated_JOCs-pangenome_pipeline
	@install global_pipeline_templates/Integrated_JOCs-pangenome_pipeline/muscle.jaccard.config $(PREFIX)/global_pipeline_templates/Integrated_JOCs-pangenome_pipeline
	@install global_pipeline_templates/Integrated_JOCs-pangenome_pipeline/muscle.j_ortholog_clusters.config $(PREFIX)/global_pipeline_templates/Integrated_JOCs-pangenome_pipeline
	@install global_pipeline_templates/Integrated_JOCs-pangenome_pipeline/pangenome.default.config $(PREFIX)/global_pipeline_templates/Integrated_JOCs-pangenome_pipeline
	@install global_pipeline_templates/Integrated_JOCs-pangenome_pipeline/pipeline.layout $(PREFIX)/global_pipeline_templates/Integrated_JOCs-pangenome_pipeline
	@install global_pipeline_templates/Integrated_JOCs-pangenome_pipeline/split_multifasta.default.config $(PREFIX)/global_pipeline_templates/Integrated_JOCs-pangenome_pipeline
	@install global_pipeline_templates/Integrated_JOCs-pangenome_pipeline/wu-blastp.default.config $(PREFIX)/global_pipeline_templates/Integrated_JOCs-pangenome_pipeline
	@install global_pipeline_templates/Integrated_JOCs-pangenome_pipeline/wu-tblastn.default.config $(PREFIX)/global_pipeline_templates/Integrated_JOCs-pangenome_pipeline
	@install global_pipeline_templates/Integrated_JOCs-pangenome_pipeline/xdformat.default.config $(PREFIX)/global_pipeline_templates/Integrated_JOCs-pangenome_pipeline
	@install global_pipeline_templates/Integrated_JOCs-pangenome_pipeline/xdformat.nuc.config $(PREFIX)/global_pipeline_templates/Integrated_JOCs-pangenome_pipeline
	mkdir -p $(PREFIX)/global_pipeline_templates/JOCs_comparative_pipeline
	@install global_pipeline_templates/JOCs_comparative_pipeline/bsml2fasta.fastamulti.config $(PREFIX)/global_pipeline_templates/JOCs_comparative_pipeline
	@install global_pipeline_templates/JOCs_comparative_pipeline/bsml2fasta.fastasingle.config $(PREFIX)/global_pipeline_templates/JOCs_comparative_pipeline
	@install global_pipeline_templates/JOCs_comparative_pipeline/jaccard.default.config $(PREFIX)/global_pipeline_templates/JOCs_comparative_pipeline
	@install global_pipeline_templates/JOCs_comparative_pipeline/j_ortholog_clusters.default.config $(PREFIX)/global_pipeline_templates/JOCs_comparative_pipeline
	@install global_pipeline_templates/JOCs_comparative_pipeline/muscle.jaccard.config $(PREFIX)/global_pipeline_templates/JOCs_comparative_pipeline
	@install global_pipeline_templates/JOCs_comparative_pipeline/muscle.j_ortholog_clusters.config $(PREFIX)/global_pipeline_templates/JOCs_comparative_pipeline
	@install global_pipeline_templates/JOCs_comparative_pipeline/pipeline.layout $(PREFIX)/global_pipeline_templates/JOCs_comparative_pipeline
	@install global_pipeline_templates/JOCs_comparative_pipeline/pipeline.xml.comment $(PREFIX)/global_pipeline_templates/JOCs_comparative_pipeline
	@install global_pipeline_templates/JOCs_comparative_pipeline/split_multifasta.default.config $(PREFIX)/global_pipeline_templates/JOCs_comparative_pipeline
	@install global_pipeline_templates/JOCs_comparative_pipeline/wu-blastp.default.config $(PREFIX)/global_pipeline_templates/JOCs_comparative_pipeline
	@install global_pipeline_templates/JOCs_comparative_pipeline/xdformat.default.config $(PREFIX)/global_pipeline_templates/JOCs_comparative_pipeline
	mkdir -p $(PREFIX)/global_pipeline_templates/nonredundant_sequence_pipeline
	@install global_pipeline_templates/nonredundant_sequence_pipeline/mummer.default.config $(PREFIX)/global_pipeline_templates/nonredundant_sequence_pipeline
	@install global_pipeline_templates/nonredundant_sequence_pipeline/mumsgraph.default.config $(PREFIX)/global_pipeline_templates/nonredundant_sequence_pipeline
	@install global_pipeline_templates/nonredundant_sequence_pipeline/pipeline.layout $(PREFIX)/global_pipeline_templates/nonredundant_sequence_pipeline
	@install global_pipeline_templates/nonredundant_sequence_pipeline/pipeline.xml.comment $(PREFIX)/global_pipeline_templates/nonredundant_sequence_pipeline
	@install global_pipeline_templates/nonredundant_sequence_pipeline/split_multifasta.default.config $(PREFIX)/global_pipeline_templates/nonredundant_sequence_pipeline
	mkdir -p $(PREFIX)/global_pipeline_templates/Pangenome
	@install global_pipeline_templates/Pangenome/pangenome.default.config $(PREFIX)/global_pipeline_templates/Pangenome
	@install global_pipeline_templates/Pangenome/pipeline.layout $(PREFIX)/global_pipeline_templates/Pangenome
	@install global_pipeline_templates/Pangenome/pipeline.xml.comment $(PREFIX)/global_pipeline_templates/Pangenome
	@install global_pipeline_templates/Pangenome/split_multifasta.default.config $(PREFIX)/global_pipeline_templates/Pangenome
	@install global_pipeline_templates/Pangenome/wu-blastp.default.config $(PREFIX)/global_pipeline_templates/Pangenome
	@install global_pipeline_templates/Pangenome/wu-tblastn.default.config $(PREFIX)/global_pipeline_templates/Pangenome
	@install global_pipeline_templates/Pangenome/xdformat.default.config $(PREFIX)/global_pipeline_templates/Pangenome
	@install global_pipeline_templates/Pangenome/xdformat.nuc.config $(PREFIX)/global_pipeline_templates/Pangenome
	@install global_pipeline_templates/Pangenome/xdformat.pepdb.config $(PREFIX)/global_pipeline_templates/Pangenome
    
	mkdir -p $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/auto_gene_curation.overlap_analysis.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/ber.post_overlap_analysis.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/ber.pre_overlap_analysis.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/bsml2fasta.final_cds.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/bsml2fasta.post_overlap_analysis.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/bsml2fasta.prediction_CDS.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/bsml2fasta.pre_overlap_analysis.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/bsml2featurerelationships.post_overlap_analysis.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/bsml2featurerelationships.pre_overlap_analysis.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/bsml2interevidence_fasta.default.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/glimmer3.iter1.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/glimmer3.iter2.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/hmmpfam.post_overlap_analysis.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/hmmpfam.pre_overlap_analysis.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/infernal.find_RNA.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/lipop.default.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/ncbi-blastx.blast_interevidence_regions.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/overlap_analysis.default.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/parse_evidence.ber_post.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/parse_evidence.ber_pre.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/parse_evidence.hmmpfam_post.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/parse_evidence.hmmpfam_pre.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/parse_evidence.hypothetical.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/parse_evidence.lipoprotein.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/parse_evidence.priam_ec.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/parse_evidence.tmhmm.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/p_func.default.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/pipeline.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/pipeline.layout $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/pipeline_summary.default.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/priam_ec_assignment.default.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/promote_gene_prediction.promote_prediction.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/ps_scan.default.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/RNAmmer.default.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/rpsblast.priam.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/signalp.default.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/split_multifasta.default.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
    @install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/split_multifasta.split_interevidence_regions.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/start_site_curation.default.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/tmhmm.default.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/train_for_glimmer3_iteration.train_for_glimmer.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/translate_sequence.default.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/translate_sequence.final_polypeptides.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/translate_sequence.translate.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/translate_sequence.translate_new_models.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/translate_sequence.translate_prediction.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/tRNAscan-SE.find_tRNA.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/wu-blastp.COGS.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/wu-blastp.post_overlap_analysis.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/wu-blastp.pre_overlap_analysis.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/xdformat.post_overlap_analysis.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline
	@install global_pipeline_templates/Prokaryotic_Annotation_Pipeline/xdformat.pre_overlap_analysis.config $(PREFIX)/global_pipeline_templates/Prokaryotic_Annotation_Pipeline

454_mplex_seq_bin_trim	:
	@echo "Installing 454_mplex_seq_bin_trim workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/454_mplex_seq_bin_trim/454_mplex_seq_bin_trim.xml $(WORKFLOW_DOCS_DIR)
	@install components/454_mplex_seq_bin_trim/454_mplex_seq_bin_trim.tmpl $(WORKFLOW_DOCS_DIR)
	@install components/454_mplex_seq_bin_trim/454_mplex_seq_bin_trim.config $(WORKFLOW_DOCS_DIR)

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
	@install components/adjust_bsml_coordinates/adjust_bsml_coordinates.xml $(WORKFLOW_DOCS_DIR)
	@install components/adjust_bsml_coordinates/adjust_bsml_coordinates.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/adjust_bsml_coordinates/adjust_bsml_coordinates.config $(WORKFLOW_DOCS_DIR)


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

aragorn	:
	@echo "Installing aragorn workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/aragorn/aragorn.xml $(WORKFLOW_DOCS_DIR)
	@install components/aragorn/aragorn.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/aragorn/aragorn.config $(WORKFLOW_DOCS_DIR)

augustus	:
	@echo "Installing augustus workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/augustus/augustus.xml $(WORKFLOW_DOCS_DIR)
	@install components/augustus/augustus.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/augustus/augustus.config $(WORKFLOW_DOCS_DIR)
	@install components/augustus/augustus.tmpl $(WORKFLOW_DOCS_DIR)/documentaion

auto_gene_curation	:
	@echo "Installing auto_gene_curation workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/auto_gene_curation/auto_gene_curation.xml $(WORKFLOW_DOCS_DIR)
	@install components/auto_gene_curation/auto_gene_curation.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/auto_gene_curation/auto_gene_curation.config $(WORKFLOW_DOCS_DIR)

backupdb  :
	@echo "Installing backupdb workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/backupdb/backupdb.backuplargetables.xml $(WORKFLOW_DOCS_DIR)
	@install components/backupdb/backupdb.backuptables.xml $(WORKFLOW_DOCS_DIR)  
	@install components/backupdb/backupdb.config $(WORKFLOW_DOCS_DIR)  
	@install components/backupdb/backupdb.xml $(WORKFLOW_DOCS_DIR)  

ber	:
	@echo "Installing ber workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/ber/ber.xml $(WORKFLOW_DOCS_DIR)
	@install components/ber/ber.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/ber/ber.config $(WORKFLOW_DOCS_DIR)

blast2go	:
	@echo "Installing blast2go workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/blast2go/blast2go.xml $(WORKFLOW_DOCS_DIR)
	@install components/blast2go/blast2go.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/blast2go/blast2go.config $(WORKFLOW_DOCS_DIR)

bsml2chado  :
	@echo "Installing bsml2chado workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/bsml2chado/bsml2chado.config $(WORKFLOW_DOCS_DIR)
	@install components/bsml2chado/bsml2chado.xml $(WORKFLOW_DOCS_DIR)  
	@install components/bsml2chado/bsml2chado.parse.xml $(WORKFLOW_DOCS_DIR)
	@install components/bsml2chado/bsml2chado.sql.xml $(WORKFLOW_DOCS_DIR)  
	@install components/bsml2chado/bsml2chado.sqlforce.xml $(WORKFLOW_DOCS_DIR)  
	@install components/bsml2chado/bsml2chado.load.xml $(WORKFLOW_DOCS_DIR)  
	@install components/bsml2chado/bsml2chado.collate.xml $(WORKFLOW_DOCS_DIR)  

bsml2gff3  :
	@echo "Installing bsml2gff3 workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/bsml2gff3/bsml2gff3.config $(WORKFLOW_DOCS_DIR)
	@install components/bsml2gff3/bsml2gff3.xml $(WORKFLOW_DOCS_DIR)  
	@install components/bsml2gff3/bsml2gff3.i1.xml $(WORKFLOW_DOCS_DIR)

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
	@install components/bsml2fasta/bsml2fasta.tmpl $(WORKFLOW_DOCS_DIR)/documentation

bsml2featurerelationships	:
	@echo "Installing bsml2featurerelationships workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/bsml2featurerelationships/bsml2featurerelationships.xml $(WORKFLOW_DOCS_DIR)
	@install components/bsml2featurerelationships/bsml2featurerelationships.config $(WORKFLOW_DOCS_DIR)

bsml2interevidence_fasta	:
	@echo "Installing bsml2interevidence_fasta workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/bsml2interevidence_fasta/bsml2interevidence_fasta.xml $(WORKFLOW_DOCS_DIR)
	@install components/bsml2interevidence_fasta/bsml2interevidence_fasta.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/bsml2interevidence_fasta/bsml2interevidence_fasta.config $(WORKFLOW_DOCS_DIR)

bsml2tbl	:
	@echo "Installing bsml2tbl workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/bsml2tbl/bsml2tbl.xml $(WORKFLOW_DOCS_DIR)
	@install components/bsml2tbl/bsml2tbl.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/bsml2tbl/bsml2tbl.config $(WORKFLOW_DOCS_DIR)
	@install components/bsml2tbl/bsml2tbl.tmpl $(WORKFLOW_DOCS_DIR)/documentation

build_supercontig	:
	@echo "Installing build_supercontig workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/build_supercontig/build_supercontigconf.ini $(WORKFLOW_DOCS_DIR)
	@install components/build_supercontig/build_supercontig-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/build_supercontig/build_supercontig-master_template.xml $(WORKFLOW_DOCS_DIR)

c_apps	:
	@echo "Compiling C applications and installing in $(PREFIX)/bin"
	${MAKE} -C src/cluster
	@install src/cluster/cluster $(PREFIX)/bin
	${MAKE} -C src/chado_record_uniq
	@install src/chado_record_uniq/chado_record_uniq $(PREFIX)/bin
	${MAKE} -C src/pangenome_profile_compare
	@install src/pangenome_profile_compare/pangenome_profile_compare $(PREFIX)/bin
	${MAKE} -C src/SnpClusterer
	@install src/SnpClusterer/SnpClusterer $(PREFIX)/bin

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
	@install components/clustalw/clustalw.tmpl $(WORKFLOW_DOCS_DIR)/documentation

combine_ncRNA_alignments	:
	@echo "Installing combine_ncRNA_alignments workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/combine_ncRNA_alignments/combine_ncRNA_alignments.xml $(WORKFLOW_DOCS_DIR)
	@install components/combine_ncRNA_alignments/combine_ncRNA_alignments.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/combine_ncRNA_alignments/combine_ncRNA_alignments.config $(WORKFLOW_DOCS_DIR)

concatenate_files	:
	@echo "Installing concatenate_files workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/concatenate_files/concatenate_files.xml $(WORKFLOW_DOCS_DIR)
	@install components/concatenate_files/concatenate_files.config $(WORKFLOW_DOCS_DIR)

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

cvtermpath_closure  :
	@echo "Installing cvtermpath_closure workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/cvtermpath_closure/cvtermpath_closure.config $(WORKFLOW_DOCS_DIR)
	@install components/cvtermpath_closure/cvtermpath_closure.xml $(WORKFLOW_DOCS_DIR)  
	@install components/cvtermpath_closure/cvtermpath_closure.process.xml $(WORKFLOW_DOCS_DIR) 

delete_analysis  :
	@echo "Installing delete_analysis workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/delete_analysis/delete_analysis-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/delete_analysis/delete_analysis-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/delete_analysis/delete_analysisconf.ini $(WORKFLOW_DOCS_DIR)  

fasta2bsml	:
	@echo "Installing fasta2bsml workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/fasta2bsml/fasta2bsml.config $(WORKFLOW_DOCS_DIR)
	@install components/fasta2bsml/fasta2bsml.xml $(WORKFLOW_DOCS_DIR)
	@install components/fasta2bsml/fasta2bsml.tmpl $(WORKFLOW_DOCS_DIR)/documentation

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

genbank2bsml	:
	@echo "Installing genbank2bsml workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/genbank2bsml/genbank2bsml.xml $(WORKFLOW_DOCS_DIR)
	@install components/genbank2bsml/genbank2bsml.i1.xml $(WORKFLOW_DOCS_DIR)  
	@install components/genbank2bsml/genbank2bsml.config $(WORKFLOW_DOCS_DIR)
	@install components/genbank2bsml/genbank2bsml.tmpl $(WORKFLOW_DOCS_DIR)/documentation

geneid	:
	@echo "Installing geneid workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/geneid/geneid.xml $(WORKFLOW_DOCS_DIR)
	@install components/geneid/geneid.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/geneid/geneid.config $(WORKFLOW_DOCS_DIR)
	@install components/geneid/geneid.tmpl $(WORKFLOW_DOCS_DIR)/documentation

genemark	:
	@echo "Installing genemark workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/genemark/genemark.xml $(WORKFLOW_DOCS_DIR)
	@install components/genemark/genemark.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/genemark/genemark.config $(WORKFLOW_DOCS_DIR)

genewise	:
	@echo "Installing genewise workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/genewise/genewise.xml $(WORKFLOW_DOCS_DIR)
	@install components/genewise/genewise.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/genewise/genewise.i2.xml $(WORKFLOW_DOCS_DIR)
	@install components/genewise/genewise.config $(WORKFLOW_DOCS_DIR)

genezilla  :
	@echo "Installing genezilla workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/genezilla/genezilla.xml $(WORKFLOW_DOCS_DIR)
	@install components/genezilla/genezilla.i1.xml $(WORKFLOW_DOCS_DIR)  
	@install components/genezilla/genezilla.config $(WORKFLOW_DOCS_DIR)
	@install components/genezilla/genezilla.tmpl $(WORKFLOW_DOCS_DIR)/documentation

genie  :
	@echo "Installing genie workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/genie/genie-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/genie/genie-master_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/genie/genie.ini $(WORKFLOW_DOCS_DIR)
	@install components/genie/genie_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/genie/genieconf.ini $(WORKFLOW_DOCS_DIR)  

genscan  :
	@echo "Installing genscan workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/genscan/genscan.config $(WORKFLOW_DOCS_DIR)
	@install components/genscan/genscan.xml $(WORKFLOW_DOCS_DIR)  
	@install components/genscan/genscan.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/genscan/genscan.tmpl $(WORKFLOW_DOCS_DIR)/documentation

genscan-plus  :
	@echo "Installing genscan-plus workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/genscan-plus/genscan-plus.config $(WORKFLOW_DOCS_DIR)
	@install components/genscan-plus/genscan-plus.xml $(WORKFLOW_DOCS_DIR)  
	@install components/genscan-plus/genscan-plus.i1.xml $(WORKFLOW_DOCS_DIR)

gff2bsml	:
	@echo "Installing gff2bsml workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/gff2bsml/gff2bsml-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/gff2bsml/gff2bsml-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/gff2bsml/gff2bsml.ini $(WORKFLOW_DOCS_DIR)
	@install components/gff2bsml/gff2bsml_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/gff2bsml/gff2bsmlconf.ini $(WORKFLOW_DOCS_DIR)

gff3_to_annotab	:
	@echo "Installing gff3_to_annotab workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/gff3_to_annotab/gff3_to_annotab.xml $(WORKFLOW_DOCS_DIR)
	@install components/gff3_to_annotab/gff3_to_annotab.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/gff3_to_annotab/gff3_to_annotab.config $(WORKFLOW_DOCS_DIR)
	@install components/gff3_to_annotab/gff3_to_annotab.tmpl $(WORKFLOW_DOCS_DIR)/documentation

gff3_to_fasta	:
	@echo "Installing gff3_to_fasta workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/gff3_to_fasta/gff3_to_fasta.xml $(WORKFLOW_DOCS_DIR)
	@install components/gff3_to_fasta/gff3_to_fasta.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/gff3_to_fasta/gff3_to_fasta.config $(WORKFLOW_DOCS_DIR)
	@install components/gff3_to_fasta/gff3_to_fasta.tmpl $(WORKFLOW_DOCS_DIR)/documentation

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
	@install components/glimmerHMM/glimmerHMM.tmpl $(WORKFLOW_DOCS_DIR)/documentation

hmmpfam	:
	@echo "Installing hmmpfam workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/hmmpfam/hmmpfam.xml $(WORKFLOW_DOCS_DIR)
	@install components/hmmpfam/hmmpfam.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/hmmpfam/hmmpfam.config $(WORKFLOW_DOCS_DIR)
	@install components/hmmpfam/hmmpfam.tmpl $(WORKFLOW_DOCS_DIR)/documentation

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
	@install components/initdb/initdb.partition.xml $(WORKFLOW_DOCS_DIR)
	@install components/initdb/initdb.config $(WORKFLOW_DOCS_DIR)

intergenic_analysis	:
	@echo "Installing intergenic_analysis workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/intergenic_analysis/intergenic_analysis.xml $(WORKFLOW_DOCS_DIR)
	@install components/intergenic_analysis/intergenic_analysis.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/intergenic_analysis/intergenic_analysis.config $(WORKFLOW_DOCS_DIR)

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

lipop:
	@echo "Installing lipop workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/lipop/lipop.xml $(WORKFLOW_DOCS_DIR)
	@install components/lipop/lipop.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/lipop/lipop.config $(WORKFLOW_DOCS_DIR)

lipoprotein_motif:
	@echo "Installing lipoprotein_motif workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/lipoprotein_motif/lipoprotein_motif.xml $(WORKFLOW_DOCS_DIR)
	@install components/lipoprotein_motif/lipoprotein_motif.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/lipoprotein_motif/lipoprotein_motif.config $(WORKFLOW_DOCS_DIR)

load_ber_for_manatee:
	@echo "Installing load_ber_for_manatee workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/load_ber_for_manatee/load_ber_for_manatee.xml $(WORKFLOW_DOCS_DIR)
	@install components/load_ber_for_manatee/load_ber_for_manatee.config $(WORKFLOW_DOCS_DIR)

mask_by_analysis	:
	@echo "Installing mask_by_analysis workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/mask_by_analysis/mask_by_analysis.xml $(WORKFLOW_DOCS_DIR)
	@install components/mask_by_analysis/mask_by_analysis.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/mask_by_analysis/mask_by_analysis.config $(WORKFLOW_DOCS_DIR)

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

metagene  :
	@echo "Installing metagene workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/metagene/metagene.xml $(WORKFLOW_DOCS_DIR)
	@install components/metagene/metagene.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/metagene/metagene.config $(WORKFLOW_DOCS_DIR)

mummer  :
	@echo "Installing mummer workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/mummer/mummer.xml $(WORKFLOW_DOCS_DIR)
	@install components/mummer/mummer.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/mummer/mummer.config $(WORKFLOW_DOCS_DIR)

mumsgraph  :
	@echo "Installing mumsgraph workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/mumsgraph/mumsgraph.xml $(WORKFLOW_DOCS_DIR)
	@install components/mumsgraph/mumsgraph.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/mumsgraph/mumsgraph.config $(WORKFLOW_DOCS_DIR)


muscle  :
	@echo "Installing muscle workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/muscle/muscle.xml $(WORKFLOW_DOCS_DIR)
	@install components/muscle/muscle.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/muscle/muscle.config $(WORKFLOW_DOCS_DIR)

ncbi-blastn	:
	@echo "Installing ncbi-blastn workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/ncbi-blastn/ncbi-blastn.xml $(WORKFLOW_DOCS_DIR)
	@install components/ncbi-blastn/ncbi-blastn.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/ncbi-blastn/ncbi-blastn.config $(WORKFLOW_DOCS_DIR)

ncbi-blastn_compress	:
	@echo "Installing ncbi-blastn_compress workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/ncbi-blastn_compress/ncbi-blastn_compress.xml $(WORKFLOW_DOCS_DIR)
	@install components/ncbi-blastn_compress/ncbi-blastn_compress.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/ncbi-blastn_compress/ncbi-blastn_compress.config $(WORKFLOW_DOCS_DIR)

ncbi-blastp	:
	@echo "Installing ncbi-blastp workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/ncbi-blastp/ncbi-blastp.xml $(WORKFLOW_DOCS_DIR)
	@install components/ncbi-blastp/ncbi-blastp.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/ncbi-blastp/ncbi-blastp.config $(WORKFLOW_DOCS_DIR)

ncbi-blastx	:
	@echo "Installing ncbi-blastx workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/ncbi-blastx/ncbi-blastx.xml $(WORKFLOW_DOCS_DIR)
	@install components/ncbi-blastx/ncbi-blastx.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/ncbi-blastx/ncbi-blastx.config $(WORKFLOW_DOCS_DIR)
	@install components/ncbi-blastx/ncbi-blastx.tmpl $(WORKFLOW_DOCS_DIR)/documentation

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
	@install components/nucmer_snps/nucmer_snps.xml $(WORKFLOW_DOCS_DIR)
	@install components/nucmer_snps/nucmer_snps.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/nucmer_snps/nucmer_snps.config $(WORKFLOW_DOCS_DIR)

obo2chado	:
	@echo "Installing obo2chado workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/obo2chado/obo2chado.xml $(WORKFLOW_DOCS_DIR)
	@install components/obo2chado/obo2chado.obo.xml $(WORKFLOW_DOCS_DIR)
	@install components/obo2chado/obo2chado.config $(WORKFLOW_DOCS_DIR)

open_reading_frames	:
	@echo "Installing open_reading_frames workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/open_reading_frames/open_reading_frames.xml $(WORKFLOW_DOCS_DIR)
	@install components/open_reading_frames/open_reading_frames.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/open_reading_frames/open_reading_frames.config $(WORKFLOW_DOCS_DIR)
	@install components/open_reading_frames/open_reading_frames.tmpl $(WORKFLOW_DOCS_DIR)/documentation

outer_membrane_domain :
	@echo "Installing outer_membrane_domain workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/outer_membrane_domain/outer_membrane_domain.xml $(WORKFLOW_DOCS_DIR)
	@install components/outer_membrane_domain/outer_membrane_domain.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/outer_membrane_domain/outer_membrane_domain.config $(WORKFLOW_DOCS_DIR)
	@install components/outer_membrane_domain/outer_membrane_domain.tmpl $(WORKFLOW_DOCS_DIR)/documentation

overlap_analysis :
	@echo "Installing overlap_analysis workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/overlap_analysis/overlap_analysis.xml $(WORKFLOW_DOCS_DIR)
	@install components/overlap_analysis/overlap_analysis.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/overlap_analysis/overlap_analysis.config $(WORKFLOW_DOCS_DIR)

p_func	:
	@echo "Installing p_func workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/p_func/p_func.xml $(WORKFLOW_DOCS_DIR)
	@install components/p_func/p_func.config $(WORKFLOW_DOCS_DIR)

pangenome_blast_parse	:
	@echo "Installing pangenome_blast_parse workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/pangenome_blast_parse/pangenome_blast_parse.xml $(WORKFLOW_DOCS_DIR)
	@install components/pangenome_blast_parse/pangenome_blast_parse.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/pangenome_blast_parse/pangenome_blast_parse.config $(WORKFLOW_DOCS_DIR)
	@install components/pangenome_blast_parse/pangenome_blast_parse.tmpl $(WORKFLOW_DOCS_DIR)/documentation

pangenome_do_analysis	:
	@echo "Installing pangenome_do_analysis workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/pangenome_do_analysis/pangenome_do_analysis.xml $(WORKFLOW_DOCS_DIR)
	@install components/pangenome_do_analysis/pangenome_do_analysis.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/pangenome_do_analysis/pangenome_do_analysis.config $(WORKFLOW_DOCS_DIR)
	@install components/pangenome_do_analysis/pangenome_do_analysis.tmpl $(WORKFLOW_DOCS_DIR)/documentation

pangenome_make_graphs	:
	@echo "Installing pangenome_make_graphs workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/pangenome_make_graphs/pangenome_make_graphs.xml $(WORKFLOW_DOCS_DIR)
	@install components/pangenome_make_graphs/pangenome_make_graphs.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/pangenome_make_graphs/pangenome_make_graphs.config $(WORKFLOW_DOCS_DIR)
	@install components/pangenome_make_graphs/pangenome_make_graphs.tmpl $(WORKFLOW_DOCS_DIR)/documentation

parse_evidence	:
	@echo "Installing parse_evidence workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/parse_evidence/parse_evidence.xml $(WORKFLOW_DOCS_DIR)
	@install components/parse_evidence/parse_evidence.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/parse_evidence/parse_evidence.config $(WORKFLOW_DOCS_DIR)
    
panther	:
	@echo "Installing panther workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/panther/panther.xml $(WORKFLOW_DOCS_DIR)
	@install components/panther/panther.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/panther/panther.config $(WORKFLOW_DOCS_DIR) 
	@install components/panther/panther.tmpl $(WORKFLOW_DOCS_DIR)/documentation

parse_for_ncRNAs	:
	@echo "Installing parse_for_ncRNA workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/parse_for_ncRNA/parse_for_ncRNAs.config $(WORKFLOW_DOCS_DIR)
	@install components/parse_for_ncRNA/parse_for_ncRNAs.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/parse_for_ncRNA/parse_for_ncRNAs.xml $(WORKFLOW_DOCS_DIR)  

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

percent_gc	:
	@echo "Installing percent_gc workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/percent_gc/percent_gc.xml $(WORKFLOW_DOCS_DIR)
	@install components/percent_gc/percent_gc.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/percent_gc/percent_gc.config $(WORKFLOW_DOCS_DIR)

phage_finder	:
	@echo "Installing phage_finder workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/phage_finder/phage_finder.xml $(WORKFLOW_DOCS_DIR)
	@install components/phage_finder/phage_finder.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/phage_finder/phage_finder.config $(WORKFLOW_DOCS_DIR)

phat	:
	@echo "Installing phat workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/phat/phat.xml $(WORKFLOW_DOCS_DIR)
	@install components/phat/phat.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/phat/phat.config $(WORKFLOW_DOCS_DIR) 

phylip_dnadist	:
	@echo "Installing phylip_dnadist workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/phylip_dnadist/phylip_dnadist.xml $(WORKFLOW_DOCS_DIR)
	@install components/phylip_dnadist/phylip_dnadist.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/phylip_dnadist/phylip_dnadist.config $(WORKFLOW_DOCS_DIR) 
	@install components/phylip_dnadist/phylip_dnadist.tmpl $(WORKFLOW_DOCS_DIR) 

pipeline_summary	:
	@echo "Installing pipeline_summary workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/pipeline_summary/pipeline_summary.xml $(WORKFLOW_DOCS_DIR)
	@install components/pipeline_summary/pipeline_summary.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/pipeline_summary/pipeline_summary.config $(WORKFLOW_DOCS_DIR) 

predict_prokaryotic_gene_function	:
	@echo "Installing predict_prokaryotic_gene_function workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/predict_prokaryotic_gene_function/predict_prokaryotic_gene_function.xml $(WORKFLOW_DOCS_DIR)
	@install components/predict_prokaryotic_gene_function/predict_prokaryotic_gene_function.config $(WORKFLOW_DOCS_DIR)
	@install components/predict_prokaryotic_gene_function/predict_prokaryotic_gene_function.tmpl $(WORKFLOW_DOCS_DIR)/documentation

predotar	:
	@echo "Installing predotar workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/predotar/predotar.xml $(WORKFLOW_DOCS_DIR)
	@install components/predotar/predotar.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/predotar/predotar.config $(WORKFLOW_DOCS_DIR)

priam_ec_assignment	:
	@echo "Installing priam_ec_assignment workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/priam_ec_assignment/priam_ec_assignment.xml $(WORKFLOW_DOCS_DIR)
	@install components/priam_ec_assignment/priam_ec_assignment.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/priam_ec_assignment/priam_ec_assignment.config $(WORKFLOW_DOCS_DIR)

promer	:
	@echo "Installing promer workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/promer/promer-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/promer/promer-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/promer/promer.ini $(WORKFLOW_DOCS_DIR)
	@install components/promer/promer_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/promer/promerconf.ini $(WORKFLOW_DOCS_DIR)  

promote_gene_prediction	:
	@echo "Installing promote_gene_prediction workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/promote_gene_prediction/promote_gene_prediction.xml $(WORKFLOW_DOCS_DIR)
	@install components/promote_gene_prediction/promote_gene_prediction.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/promote_gene_prediction/promote_gene_prediction.config $(WORKFLOW_DOCS_DIR) 

ps_scan	:
	@echo "Installing ps_scan workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/ps_scan/ps_scan.xml $(WORKFLOW_DOCS_DIR)
	@install components/ps_scan/ps_scan.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/ps_scan/ps_scan.config $(WORKFLOW_DOCS_DIR) 

psiblast :
	@echo "Installing psiblast workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/psiblast/psiblast.xml $(WORKFLOW_DOCS_DIR)
	@install components/psiblast/psiblast.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/psiblast/psiblast.config $(WORKFLOW_DOCS_DIR)

psipred	:
	@echo "Installing psipred workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/psipred/psipred-master.ini $(WORKFLOW_DOCS_DIR)
	@install components/psipred/psipred-master_template.xml $(WORKFLOW_DOCS_DIR)
	@install components/psipred/psipred.ini $(WORKFLOW_DOCS_DIR)
	@install components/psipred/psipred_template.xml $(WORKFLOW_DOCS_DIR)  
	@install components/psipred/psipredconf.ini $(WORKFLOW_DOCS_DIR)  

quicktree	:
	@echo "Installing quicktree workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/quicktree/quicktree.xml $(WORKFLOW_DOCS_DIR)
	@install components/quicktree/quicktree.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/quicktree/quicktree.config $(WORKFLOW_DOCS_DIR)
	@install components/quicktree/quicktree.tmpl $(WORKFLOW_DOCS_DIR)/documentation

rdp_classifier	:
	@echo "Installing rdp_classifier workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/rdp_classifier/rdp_classifier.xml $(WORKFLOW_DOCS_DIR)
	@install components/rdp_classifier/rdp_classifier.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/rdp_classifier/rdp_classifier.config $(WORKFLOW_DOCS_DIR)
	@install components/rdp_classifier/rdp_classifier.tmpl $(WORKFLOW_DOCS_DIR)/documentation

refresh_chado_mart	:
	@echo "Installing refresh_chado_mart workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/refresh_chado_mart/refresh_chado_mart.xml $(WORKFLOW_DOCS_DIR)
	@install components/refresh_chado_mart/refresh_chado_mart.load.xml $(WORKFLOW_DOCS_DIR)
	@install components/refresh_chado_mart/refresh_chado_mart.config $(WORKFLOW_DOCS_DIR)
	@install components/refresh_chado_mart/refresh_chado_mart.sqlforce.xml $(WORKFLOW_DOCS_DIR)
	@install components/refresh_chado_mart/refresh_chado_mart.tmpl $(WORKFLOW_DOCS_DIR)/documentation

repeatmasker	:
	@echo "Installing repeatmasker workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/repeatmasker/repeatmasker.xml $(WORKFLOW_DOCS_DIR)
	@install components/repeatmasker/repeatmasker.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/repeatmasker/repeatmasker.config $(WORKFLOW_DOCS_DIR)

RNAmmer:
	@echo "Installing RNAmmer workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/RNAmmer/RNAmmer.xml $(WORKFLOW_DOCS_DIR)
	@install components/RNAmmer/RNAmmer.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/RNAmmer/RNAmmer.config $(WORKFLOW_DOCS_DIR)
	@install components/RNAmmer/RNAmmer.tmpl $(WORKFLOW_DOCS_DIR)/documentation

rpsblast :
	@echo "Installing rpsblast workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/rpsblast/rpsblast.xml $(WORKFLOW_DOCS_DIR)
	@install components/rpsblast/rpsblast.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/rpsblast/rpsblast.config $(WORKFLOW_DOCS_DIR)

seg:
	@echo "Installing seg workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/seg/seg.xml $(WORKFLOW_DOCS_DIR)
	@install components/seg/seg.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/seg/seg.config $(WORKFLOW_DOCS_DIR)
	@install components/seg/seg.tmpl $(WORKFLOW_DOCS_DIR)/documentation

signalp	:
	@echo "Installing signalp workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/signalp/signalp.xml $(WORKFLOW_DOCS_DIR)
	@install components/signalp/signalp.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/signalp/signalp.config $(WORKFLOW_DOCS_DIR)
	@install components/signalp/signalp.tmpl $(WORKFLOW_DOCS_DIR)/documentation

snap	:
	@echo "Installing snap workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/snap/snap.xml $(WORKFLOW_DOCS_DIR)
	@install components/snap/snap.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/snap/snap.config $(WORKFLOW_DOCS_DIR)

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
	@install components/split_multifasta/split_multifasta.tmpl $(WORKFLOW_DOCS_DIR)/documentation

stap_aligner	:
	@echo "Installing stap_aligner workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/stap_aligner/stap_aligner.xml $(WORKFLOW_DOCS_DIR)
	@install components/stap_aligner/stap_aligner.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/stap_aligner/stap_aligner.config $(WORKFLOW_DOCS_DIR)
	@install components/stap_aligner/stap_aligner.tmpl $(WORKFLOW_DOCS_DIR)/documentation

stap_classifier	:
	@echo "Installing stap_classifier workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/stap_classifier/stap_classifier.xml $(WORKFLOW_DOCS_DIR)
	@install components/stap_classifier/stap_classifier.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/stap_classifier/stap_classifier.config $(WORKFLOW_DOCS_DIR)
	@install components/stap_classifier/stap_classifier.tmpl $(WORKFLOW_DOCS_DIR)/documentation

start_site_curation	:
	@echo "Installing start_site_curation workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/start_site_curation/start_site_curation.xml $(WORKFLOW_DOCS_DIR)
	@install components/start_site_curation/start_site_curation.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/start_site_curation/start_site_curation.config $(WORKFLOW_DOCS_DIR)
	@install components/start_site_curation/start_site_curation.tmpl $(WORKFLOW_DOCS_DIR)/documentation

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
	@install components/tabula_rasa/tabula_rasa.xml $(WORKFLOW_DOCS_DIR)
	@install components/tabula_rasa/tabula_rasa.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/tabula_rasa/tabula_rasa.config $(WORKFLOW_DOCS_DIR)

targetp	:
	@echo "Installing targetp workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/targetp/targetp.xml $(WORKFLOW_DOCS_DIR)
	@install components/targetp/targetp.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/targetp/targetp.config $(WORKFLOW_DOCS_DIR) 

tbl2asn	:
	@echo "Installing tbl2asn workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/tbl2asn/tbl2asn.xml $(WORKFLOW_DOCS_DIR)
	@install components/tbl2asn/tbl2asn.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/tbl2asn/tbl2asn.config $(WORKFLOW_DOCS_DIR)
	@install components/tbl2asn/tbl2asn.tmpl $(WORKFLOW_DOCS_DIR)/documentation

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
	@install components/translate_sequence/translate_sequence.xml $(WORKFLOW_DOCS_DIR)
	@install components/translate_sequence/translate_sequence.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/translate_sequence/translate_sequence.config $(WORKFLOW_DOCS_DIR)

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

tRNAscan-SE_compress	:
	@echo "Installing tRNAscan-SE_compress workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/tRNAscan-SE_compress/tRNAscan-SE_compress.xml $(WORKFLOW_DOCS_DIR)
	@install components/tRNAscan-SE_compress/tRNAscan-SE_compress.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/tRNAscan-SE_compress/tRNAscan-SE_compress.config $(WORKFLOW_DOCS_DIR)

wait	:
	@echo "Installing wait workflows in $(WORKFLOW_DOCS_DIR)"
	@install components/wait/wait.xml $(WORKFLOW_DOCS_DIR)
	@install components/wait/wait.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/wait/wait.i2.xml $(WORKFLOW_DOCS_DIR)
	@install components/wait/wait.config $(WORKFLOW_DOCS_DIR)
	@install components/wait/wait.tmpl $(WORKFLOW_DOCS_DIR)/documentation

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
	@install components/wu-blastx/wu-blastx.xml $(WORKFLOW_DOCS_DIR)
	@install components/wu-blastx/wu-blastx.i1.xml $(WORKFLOW_DOCS_DIR)
	@install components/wu-blastx/wu-blastx.config $(WORKFLOW_DOCS_DIR)
	@install components/wu-blastx/wu-blastx.tmpl $(WORKFLOW_DOCS_DIR)/documentation

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
