package MG::ManFasta;

use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use Data::Dumper;
use XML::Simple;

@ISA = qw(Exporter);
@EXPORT = qw(&revcomp &translate_dna &find_frame);

sub revcomp {
    my $seq = uc(shift(@_));
    my %comp = ("A","T","T","A","C","G","G","C",'N','N','X','X');
    my @seq = split(//, $seq);
    my @revcomp = ();
    foreach $_ (@seq) {
	unshift @revcomp, $comp{$_};
    }
    return join("", @revcomp);
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
    $protein =~ s/\*$//g;
    return $protein;
}
sub find_frame {
    my $ori = shift @_;
    my $prot_seq = shift @_;
    $dna_seq = $ori;
    foreach (0..2) {
	my $trans = translate_dna($dna_seq);
	if ($trans eq $prot_seq) {
	    return ($_);
	}
	$dna_seq = substr($dna_seq,1);
    }
    $dna_seq = $ori;
    foreach (0..2) {
	my $trans = translate_dna($dna_seq);
	unless ($trans =~ m/\*/) {
	    return ($_);
	}
	$dna_seq = substr($dna_seq,1);
    }
    warn "No Frame Matches";
    return (0);
}

return 1;
