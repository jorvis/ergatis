use ExtUtils::MakeMaker;
require ExtUtils::Manifest;
require "perl2wrapper_papyrus.pl";

if(-e "MANIFEST"){}
else{
    ExtUtils::Manifest::mkmanifest();
  }

my($scripts,$wrapperstr) = &get_bins();


WriteMakefile(
	      'NAME'=>'papyrus',
	      'EXE_FILES'=>$scripts,
	      'clean' =>{FILES => "MANIFEST.bak $wrapperstr README.symlinks"}
	      );

