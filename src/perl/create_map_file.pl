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

die("Either --input_file_list or --input_file is required") unless( @files );

open(my $OFH,">$outputFile") or die "Error in writing to the output file, $outputFile, $!\n";
foreach my $f ( @files ) {
    printMapInfo(trim($f), $OFH);
}
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
	}
	elsif(trim($line) =~ /^ORGANISM\s+(.+)\s*$/) {
	    unless($locus) {
		die "No locus found in the genbank file, $genbankFile, $!\n";
	    }
	    print $ofh $1." ".$locus."\n";
	    last;
	}
    }
    close FH;
}

