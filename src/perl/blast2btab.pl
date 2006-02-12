#!/usr/local/bin/perl
use strict;
BEGIN {
use BPlite;
}

# multiple (concatenated) BLAST reports

my $multiple_report = new BPlite::Multi(\*STDIN);
while(my $report = $multiple_report->nextReport) {
    if(defined $report){
	my ($query_name) = ($report->query =~ /(\S+)/);
	my $search_database = $report->database;
	while(my $sbjct = $report->nextSbjct) {
	    if(defined $sbjct){
		my ($dbmatch_accession) = ($sbjct->name =~ /\>(\S+)/);
		my $accession_length = $sbjct->length;
		while(my $hsp = $sbjct->nextHSP) {
		    if(defined $hsp){
			my $date = "N/A";
			my $query_length = $report->queryLength();
			my $blast_program = "N/A";
			my $start_query = $hsp->queryBegin;
			my $stop_query = $hsp->queryEnd;
			my $start_hit = $hsp->sbjctBegin;
			my $stop_hit = $hsp->sbjctEnd;
			my $percent_identity = $hsp->percent;
			my $percent_similarity = $hsp->percent;
			my $bit_score = $hsp->bits;
			my $chain_number = 'N/A';
			my $segment_number = 'N/A';
			my $dbmatch_header = $dbmatch_accession;
			my $unknown1 = 'N/A';
			my $unknown2 = 'N/A';
			my $e_value = $hsp->P;
			my $p_value = $hsp->P;
			print "$query_name\t$date\t$query_length\t$blast_program\t$search_database\t$dbmatch_accession\t$start_query\t$stop_query\t$start_hit\t$stop_hit\t$percent_identity\t$percent_similarity\t$bit_score\t$chain_number\t$segment_number\t$dbmatch_header\t$unknown1\t$unknown2\t$e_value\t$p_value\n";
		    }
		}
	    }
	}
    }
}
