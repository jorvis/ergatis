package MG::ManFasta;

use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use Datastore::MD5;

@ISA = qw(Exporter);
@EXPORT = qw(&recomp &parse_fasta_format &translate_dna);

sub revcomp {
    my $seq = shift(@_);
    my %comp = ("A","T","T","A","C","G","G","C",'N','N','X','X');
    my @seq = split(//, $seq);
    my @revcomp = ();
    foreach $_ (@seq) {
	unshift @revcomp, $comp{$_};
    }
    return join("", @revcomp);
}

sub parse_fasta_format {
    my $infile = shift @_;
    my %hash = ();
    open SEQ, "<$infile" or die $!;
    $/ = "\n>";
    my $i = 0;
    while (my $record = <SEQ>) {
        chomp($record);
	my ($header, @seq) = split(/\n/, $record);
	$header =~ s/^>//g;
	$header = (split(/\s+/, $header))[0];
	$hash{$header} = join("", @seq);
	$hash{$header} = uc($hash{$header});
    }
    $/ = "\n";
    close SEQ;
    return %hash;
}

sub translate_dna {
    my $seq = uc(shift(@_));
    $seq =~ s/U/T/g;
    my %dna_protein = ("GCA"=>"A","GCC"=>"A","GCG"=>"A","GCT"=>"A","TGC"=>"C","TGT"=>"C",
		       "GAC"=>"D","GAT"=>"D","GAA"=>"E","GAG"=>"E","TTC"=>"F","TTT"=>"F",
		       "GGA"=>"G","GGC"=>"G","GGG"=>"G","GGT"=>"G","CAC"=>"H","CAT"=>"H",
		       "ATA"=>"I","ATC"=>"I","ATT"=>"I","AAA"=>"K","AAG"=>"K","TTA"=>"L",
		       "TTG"=>"L","CTA"=>"L","CTC"=>"L","CTG"=>"L","CTT"=>"L","ATG"=>"M",
		       "AAC"=>"N","AAT"=>"N","CCA"=>"P","CCC"=>"P","CCG"=>"P","CCT"=>"P",
		       "CAA"=>"Q","CAG"=>"Q","CGA"=>"R","CGC"=>"R","CGG"=>"R","CGT"=>"R",
		       "AGA"=>"R","AGG"=>"R","TCA"=>"S","TCC"=>"S","TCG"=>"S","TCT"=>"S",
		       "AGC"=>"S","AGT"=>"S","ACA"=>"T","ACC"=>"T","ACG"=>"T","ACT"=>"T",
		       "GTA"=>"V","GTC"=>"V","GTG"=>"V","GTT"=>"V","TGG"=>"W","TAC"=>"Y",
		       "TAT"=>"Y","TAA"=>"*","TAG"=>"*","TGA"=>"*");
    my @dna  = split(//, $seq);
    my $num_codons = int(($#dna + 1)/3);
    my @prot = ();
    foreach (1..$num_codons) {
	my $codon = join("", splice(@dna, 0, 3));
	if ($dna_protein{$codon}) {
	    push @prot, $dna_protein{$codon};   
	}else {
	    push @prot, 'X';
	}
    }
    my $protein = join("", @prot);
    return $protein;
}
return 1;
