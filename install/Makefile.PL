use ExtUtils::MakeMaker;
require ExtUtils::Manifest;

if(stat "MANIFEST"){}
else{
    ExtUtils::Manifest::mkmanifest();
  }

WriteMakefile(
	      'NAME'=>'papyrus',
	      'EXE_FILES'    => [
				 'bin/add_bit_score.pl',
				 'bin/bit_score_analysis.pl',
				 'bin/make_pep_bsml.pl',
                                 'bin/make_each_seq_bsml.pl',
                                 'bin/do_all_vs_all.pl',
                                 'bin/btab2bsml.pl',
                                 'bin/gene_pos_bsml.pl',
                                 'bin/gene_match_bsml.pl',
                                 'bin/run_pe.pl',
                                 'bin/pe2bsml.pl',
                                 'bin/gene_boundaries_bsml.pl',
                                 'bin/region2bsml.pl',
                                 'bin/do_blastp.pl',
                                 'bin/max_match_post_process.pl'
                                 'bin/make_asmbl_fasta.pl',
                                 'bin/mummer2bsml.pl'
                                ],

	      );
