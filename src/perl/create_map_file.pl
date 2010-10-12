#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $files = "";
my $file = "";
my $outputFile = "";

GetOptions('input_file_list=s' => \$files,
			'input_file=s' => \$file,
			'output_file=s' => \$outputFile);
unless($outputFile) {
	die "Must provide an output file\n";
}
open(OFH,">$outputFile") or die "Error in writing to the output file, $outputFile, $!\n";
if($files) {
	open(TAG,"<$files") or die "Error in opening the file, $files, $!\n";
	while(my $filePath = <TAG>) {
		chomp $filePath;
		if(trim($filePath) =~ /(\S+)/) {
			printMapInfo($1);
		}
	}
	close TAG;	
} elsif($file) {
	printMapInfo(trim($file));
} else {
	die "Usage: create_pangenome_map_file.pl --input_file_list or --input_file and --output_file\n";
}
close OFH;
exit(0);

sub trim {
	my ($param) = @_;
	$param =~ s/^\s+//;
	$param =~ s/\s+$//;
	return $param;
}

sub printMapInfo {
	my ($genbankFile) = @_;
	open (FH, "<$genbankFile") or die "Error in opening the file, $genbankFile, $!\n";
	my $locus;
	while(my $line = <FH>) {
		if(trim($line) =~ /^LOCUS\s+(\S+)\s+/) {
			$locus = $1;
		}
		elsif(trim($line) =~ /^ORGANISM\s+(.+)\s*$/) {
			unless($locus) {
				die "No locus found in the genbank file, $genbankFile, $!\n";
			}
			print OFH $1." ".$locus."\n";
			last;
		}
	}
	close FH;
}

