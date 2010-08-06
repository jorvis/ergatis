#!/usr/bin/perl -w
#get_seq_by_metagene.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use MG::RunParse;
use MG::ManFasta;
use Bio::Seq::RichSeq;
use Bio::SeqIO;
use MG::Math;

my %options = ();
my $results = GetOptions (\%options,
                          'input_file|i=s',
                          'fasta_input|f=s',
                          'output|o=s',
			  'format|m=s',
			  'prefix|p=s',
                          'help|h');

my $out = Bio::SeqIO->new(-file  => ">$options{'output'}.gbk",
			  -format => $options{'format'});

my ($orfs, $orf_count) = parse_metagene($options{'input_file'});
open IN, "<$options{fasta_input}" or die $!;
open TBL, ">$options{'output'}\.txt" or die $!;
my $j = 0;
$/ = "\n>";
while ($record = <IN>) {
    chomp $record;
    next if $record =~ m/^\s*$/;
    $record =~ s/^>//g;
    my ($defline, @seq) = split(/\n/, $record);
    my $sequence = join("", @seq);
    my ($sequence_id, @other) = split(/\s+/, $defline);
    my $sequence_defline = join(" ", @other);
    $newseq = Bio::Seq::RichSeq->new(-seq=>$sequence,
				     -molecule=>'DNA',
				     -accession_number=>$sequence_id,
				     -locus=>$sequence_id,
				     -primary_id=>$sequence_id,
				     -display_id=>$sequence_id
				    );
    my @tax_class = ('human metagenome','organismal metagenomes,metagenomes','unclassified sequences');
    my $species = Bio::Species->new(-ncbi_taxid=>646099,-classification => \@tax_class);
    $newseq->species($species);
    $newseq->add_date(ncbidate());
    $newseq->desc($sequence_id);
    my $feat = Bio::SeqFeature::Generic->new(-start=>1,-end=>length($sequence),
					     -primary=>'source',
					     -tag=>{organism=>'human metagenome',mol_type=>'genomic DNA',
						    isolation_source=>'Homo sapiens',collection_date=>ncbidate()});
    
    $newseq->add_SeqFeature($feat);
    foreach my $orf(@{$orfs->{$sequence_id}}) {
	$j++;
	my $newacc = $options{'prefix'}.sprintf("%06s",$j);
	$orf_seq = substr($sequence, $orf->{'offset'}, $orf->{'endpos'} - $orf->{'offset'});
	if ($orf->{'complement'} == 1) {
	    $orf_seq = revcomp($orf_seq); 
	}
	$orf_seq = substr($orf_seq, $orf->{'frame'});
	my $feat = Bio::SeqFeature::Generic->new(-start=>$orf->{'startpos'},-end=>$orf->{'endpos'},
						 -strand=>$orf->{'strand'},-primary_tag=>'CDS',-source_tag=>'metagene',
						 -display_name=>$orf->{id},-score=>$orf->{score},-seq_id=>$newacc,
						 -frame=>$orf->{frame},-tag=>{translation=>translate_dna($orf_seq),
									      product=>'hypothetical protein',
									      locus_tag=>$orf->{id},
									      mol_type=>'genomic DNA',
									      organism=>'human metagenome',
									      isolation_source=>'Homo sapiens',
									      protein_id=>$newacc,
									      collection_date=>ncbidate(),
									     }
						);
	$newseq->add_SeqFeature($feat);
	my $complete = 1;
	$complete = 0 if ($orf->{'5_partial'} || $orf->{'3_partial'});
	print TBL join("\t",$sequence_id,$newacc,$orf->{'startpos'},$orf->{'endpos'},$orf->{'strand'},$orf->{frame},$complete,$orf->{score},$orf->{model},'metagene'),"\n";
    }
    $out->write_seq($newseq);
}
$/ = "\n";
close TBL;
