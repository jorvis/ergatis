#!/usr/bin/perl

=head1 NAME

identify_human_contaminants.pl - allows you to identify human contaminants and divide sequences in those similar to human, ref bacterial genomes and those that are similar to nothing or are inconclusive.

=head1 SYNOPSIS

USAGE:identify_human_contaminats.pl 
            --human_contams=/path/to/blast_against_humangenome.blsout',
	    --ref_geno_search=/path/to/blast_against_refgenome.blsout',
	    --fasta_file=/path/to/query_fastafile.fa',
	    --mincov=fraction_coverage',
	    --minpercid=fraction_id',
	    --output_prefix=/path/to/outfiles'
            --list
            --algorithm=search_algorithm
            --for_assembly=1

=head1 OPTIONS

B<--human_contams,-m>
    blast output file from a search of metagenomic sequences against the human genome database (dna)

B<--ref_geno_search,-r>
    blast output file from a search of metagenomic sequences against the reference genome database (dna)

B<--algorithm,-a>
    which algorithm was used to produce search results -- fasta, blast, blasttable for blast -m8
    default blast

B<--fasta_file,-f>
    metagenomic fasta file, that was used at the input for the blast searches done against ref_geno_search and human_contams

B<--mincov, -c>
    This is the minimum coverage acceptable to determine a "organism" hit; default 0.9

B<--minpercid, -p>
    This is the minimum coverage acceptable to determine a "organism" hit; default 0.9

B<--output_prefix,-o>
    output files names, output files would be output_prefix.human.fasta, output_prefix.refgeno.fasta, output_prefix.inconclusive.fasta 
B<--for_assembly,-b>
    1 for assembly and 0 for gene finding

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

B<--list,-t>
    The files are list and not individual files

=head1  DESCRIPTION

this scripts aims to identify human contaminants.

=head1  INPUT

blast results of your sequences against the human genome and a set of bacterial reference genomes.  also the script requires the original blast search input fasta file.

=head1  OUTPUT

the output for this script are several fasta files -- human contaminants,  reference genome hits and unknown hits

=head1  CONTACT

    Brandi Cantarel
    bcantarel@som.umaryland.edu

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use File::Basename;
use Ergatis::Logger;
use MG::ParseSeq;
use MG::RunParse;
use MG::SqlFunc;


my %options = ();
my $results = GetOptions (\%options, 
			  'human_contams|m=s',
			  'ref_geno_search|r=s',
			  'fasta_file|f=s',
			  'mincov|c=s',
			  'minpercid|p=s',
			  'output_prefix|o=s',
			  'help|h',
			  'algorithm|a=s',
			  'list|t',
			 );

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				 'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

unless ($options{human_contams} && $options{ref_geno_search} && $options{fasta_file}) {
    pod2usage({-exitval=>0, -verbose => 2, -output => \*STDOUT});
}
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

$options{'mincov'} = 0.80 unless $options{'mincov'};
$options{'minpercid'} = 0.80 unless $options{'minpercid'};
$options{'algorithm'} = 'nucmer' unless $options{'algorithm'};

use vars qw(%human %seqs %ref %length %bacdbfiles %bin);

if ($options{'list'}) {
    my $human_list = $options{human_contams};
    my $bac_list = $options{ref_geno_search};
    my $seq_list = $options{fasta_file};
    open BLS, "<$human_list" or die $!;
    open TEMP, ">$human_list\.concat" or die $!;
    print TEMP $human_list," ",$seq_list,"\nNUCMER\n";
    while (my $line = <BLS>) {
	chomp($line);
	my $content = `cat $line`;
	my ($compare,$info) = split(/\nNUCMER\s*\n/, $content);
	my ($libfile, $queryfile) = split(/\s+/, $compare);
	print TEMP $info;
    }
    close BLS;
    close TEMP;
    $options{human_contams} = "$human_list\.concat";
    open BLS, "<$bac_list" or die $!;
    open TEMP, ">$bac_list\.concat" or die $!;
    print TEMP $bac_list," ",$seq_list,"\nNUCMER\n";
    while (my $line = <BLS>) {
	chomp($line);
	my $content = `cat $line`;
	my ($compare,$info) = split(/\nNUCMER\s*\n/, $content);
	my ($libfile, $queryfile) = split(/\s+/, $compare);
	$bacdbfiles{$libfile} = 1;
	print TEMP $info;
    }
    close BLS;
    close TEMP;
    $options{ref_geno_search}  = "$bac_list\.concat";
    open BLS, "<$seq_list" or die $!;
    open TEMP, ">$seq_list\.concat" or die $!;
    while (my $line = <BLS>) {
	chomp($line);
	my $content = `cat $line`;
	print TEMP $content;
    }
    close BLS;
    close TEMP;
    $options{fasta_file}  = "$seq_list\.concat";
}
%human = runparse_showcoords($options{human_contams},$options{'minpercid'},$options{'mincov'});
%ref = runparse_showcoords($options{ref_geno_search},$options{'minpercid'},$options{'mincov'});

open RTAX, ">$options{'output_prefix'}\.read_tax.txt" or die $!;
open BAC, ">$options{'output_prefix'}\.refgeno.fsa" or die $!;
open ELSE, ">$options{'output_prefix'}\.unclassified.fsa" or die $!;
open FORASS, ">$options{'output_prefix'}\.forassembly.list" or die $!;

$/ = "\n>";

open SEQIN, "<$options{fasta_file}" or die $!;
while (my $record = <SEQIN>) {
    chomp($record);
    my ($header,$seq) = split(/\n/, $record);
    $header =~ s/^>//;
    my $query = (split(/\s+/, $header))[0];
    my $qname = $query;
    $qname = (split(/_/, $qname))[-1];
    my $length = length($seq);
    next F1 unless ($length > 1);
    if ($human{$query}) {
	my @hits = sort {$a->[2] <=> $b->[2] || $a->[1] <=> $b->[1]} @{$human{$query}};
	my ($hname, $alen, $percid, $per_qcov, $per_hcov, $qlen, $hlen, $lbegin, 
	    $lend, $qbegin, $qend) = @{$hits[0]};
	if ($ref{$query}) {
	    print RTAX join("\t", $qname,0,$hname,$percid,$per_qcov),"\n";
	}else {
	    print RTAX join("\t", $qname,9606,$hname,$percid,$per_qcov),"\n";
	}
    }elsif ($ref{$query}) {
	my @hits = sort {$a->[2] <=> $b->[2] || $a->[1] <=> $b->[1]} @{$ref{$query}};
	my @thit = @{$hits[0]};
	my ($gi, $acc, $taxid) = split(/\|/, $thit[0]);
	print RTAX join("\t", $qname,$taxid,$thit[0],$thit[2],$thit[3]),"\n";
	print BAC ">",$qname,"\n",$seq,"\n";
	push @{$bin{$taxid}}, $qname;
    }else {
	print FORASS $qname,"\n";
	print ELSE ">",$qname,"\n",$seq,"\n";
    }
}
close RTAX;
close ELSE;
close BAC;

$/ = "\n";
my %taxid2file = ();
open FOUT, ">$options{'output_prefix'}\.foramoscmp" or die $!;
open SNGL, ">$options{'output_prefix'}\.singleton.list" or die $!;
foreach my $bacdb (keys %bacdbfiles) {
    my $head = `grep '>' $bacdb|head -n 1`;
    chomp($head);
    $head =~ s/^>//;
    $head =~ s/\n$//;
    my ($gi, $acc, $taxid) = split(/\|/, $head);
    $taxid2file{$taxid} = $bacdb;
}
foreach my $taxid (keys %bin) {
    if (scalar(@{$bin{$taxid}}) > 100) {
	print FOUT $taxid,"\t",$taxid2file{$taxid},"\t",join(",",@{$bin{$taxid}}),"\n";
    }else {
	print SNGL join("\n",@{$bin{$taxid}}),"\n";
    }
}
