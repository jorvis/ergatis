#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $files = "";
my $file = "";
my $outputFile = "";

GetOptions('input_file_list=s' => \$files,
	   'input_file=s' => \$file,
	   'output_file=s' => \$outputFile);

die("Must provide an output file") unless($outputFile);

my @files;
if( $files ) {
    foreach my $list ( split(/[,\s]+/, $files) ) {
	open(IN, "< $list") or die("Could not open list file $list: $!");
	chomp( my @tmp = <IN> );
	close(IN);
	push( @files, @tmp );
    }
}

if( $file ) {
    my @tmp = split(/[,\s]+/, $file);
    push(@files, @tmp);
}

unless( @files ) {
    system("touch $outputFile") && die("Could not run command: touch $outputFile: $@");
    exit(0);
}

my ($file_base,$file_dir,$file_ext) = fileparse($outputFile,qr/\.[^.]*/);
my $listFile = $file_dir."/genbank.list";
open(my $OFL,">$listFile") or die "Error in writing to the output list file, $listFile, $!\n";
open(my $OFH,">$outputFile") or die "Error in writing to the output file, $outputFile, $!\n";
foreach my $f ( @files ) {
    if(-e $f) {
        print $OFL "$f\n";
    } else {
        die "$f does not exist, $!\n";
    }
    printMapInfo(trim($f), $OFH);
}
close $OFL;
close $OFH;

sub trim {
    my ($param) = @_;
    $param =~ s/^\s+//;
    $param =~ s/\s+$//;
    return $param;
}

sub printMapInfo {
    my ($genbankFile, $ofh) = @_;
    open (FH, "<$genbankFile") or die "Error in opening the file, $genbankFile, $!\n";
    my $locus;
    while(my $line = <FH>) {
	if(trim($line) =~ /^LOCUS\s+(\S+)\s+/) {
	    $locus = $1;

	    # if the locus is an IGS internal identifier, parse it.
	    if( $locus =~ /([^\.]+)\.\w+\.\w+/ ) {
		$locus = $1;
	    }

	} elsif(trim($line) =~ /^ORGANISM\s+(.+)\s*$/) {
	    unless($locus) {
		die "No locus found in the genbank file, $genbankFile, $!\n";
	    }
	    print $ofh $1." ".$locus."\n";
	    last;
	}
    }
    close FH;
}
