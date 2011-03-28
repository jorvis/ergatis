#!/usr/bin/perl

use strict;

#Tophat needs the BOWTIE_HOME env var set, and bowtie and samtools to be available
my $bowtie_path = $ARGV[0];
my $samtools_path = $ARGV[1];

my $tophat_bin = $ARGV[2];

my @args = @ARGV[3 .. $#ARGV];

umask 0000;

$ENV{"BOWTIE_HOME"} = $bowtie_path;
$ENV{"PATH"} = $ENV{"PATH"} . ":$bowtie_path". ":$samtools_path";

print "\nPATH set to : ".$ENV{"PATH"};
print "\nExecuting $tophat_bin ".join(" ", @args)."\n";

exec ( $tophat_bin , @args);

