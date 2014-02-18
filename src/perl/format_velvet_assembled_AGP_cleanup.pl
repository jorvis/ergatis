#!/usr/bin perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Basename;
use FindBin qw($Bin);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options,
                          "tbl_file|t=s",
                          "output_file|o=s",
                          "log|l=s",
                          "help|h") || pod2usage();
my ($tbl,$head, $tag);
my %contig = ();
## Display documentation
if( $options{'help'} ){
        pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
if (!($options{'output_file'})) {
	die "ERROR : Output file is required\n";
}
if (!($options{'tbl_file'})) {
	die "ERROR : tbl file is required\n";
}
open(FR, "< $options{'tbl_file'}") or die "Could not open $options{'tbl_file'} file for reading";
while($tbl = <FR>) {
	chomp($tbl);
	next if ($tbl =~ /^\s*$/);
	if($tbl =~ /^>Feature\s+(\S+)/) {
		$head = $1;
		next;
	} 
	if($tbl =~ /^\s+locus_tag\s+(.+)/) {
		$tag = substr($1, 0, -1);
		$contig{$tag}{$head} = 1;
	}
}
close(FR);
open(FW, "> $options{'output_file'}") or die "Could not open $options{'output_file'} file for writing";
open(FR, "< $options{'tbl_file'}") or die "Could not open $options{'tbl_file'} file for reading";
while($tbl = <FR>) {
	chomp($tbl);
	next if ($tbl =~ /^\s*$/);
	if($tbl =~ /^\s+locus_tag\s+(.+)/) {
		my $subpart = "";
		my $id = $1;
		$subpart = chop($id) if(substr($1, -1, 1) =~ /[A-Z]/);
		my $n = keys %{$contig{$id}};
		print FW "\t\t\tlocus_tag\t";
		if($n == 1 && length($subpart) > 0) {
			print FW "$id\n";
		} else {
			print FW "$id$subpart\n";
		}
	}
# Removing alphabet from the protein_id if the last half of the split gene is not present in the .tbl file	
	elsif($tbl =~ /^\s+protein_id\s+(.+)/) {
		my $subpart = "";
		my $id = $1;
		my @desc = split(/\|/,$id);
		$subpart = chop($desc[$#desc]) if(substr($desc[$#desc], -1, 1) =~ /[A-Z]/);
		my $n = keys %{$contig{$desc[$#desc]}};
		print FW "\t\t\tprotein_id\t";
		if($n == 1 && length($subpart) > 0) {
			print FW substr($id, 0, -1)."\n";
		} else {
			print FW "$id\n";
		}
# Removing extra notes from CDS if the other half of the gene is not present in .tbl file
	} elsif($tbl =~ /^\s+note\s+(.+)/) {
		$1 =~ m/end is \w+ (\S+) on contig (\S+)/;
		print FW "$tbl\n" if(exists($contig{$1}{$2}));
	} else {
		print FW "$tbl\n";
	}
}
close(FR);
close(FW);
