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
				 'bin/bit_score_analysis.pl'
				 ],

	      );
