#!/usr/bin/perl -w
#BSMLChecker.pm

package BSMLChecker;

use lib ('/usr/local/projects/ergatis/package-latest/lib/perl5');
use Ergatis::Logger;

#  Add more subroutines for checking things as time passes


my $gene_count = 0;
my $cds_count = 0;
my $polypeptide_count = 0;
my $transcript_count = 0;
my $exon_count = 0;
my $tRNA_count = 0;
my %gene_hash = ();
my $id = "";


# Reset gene hash and counts for the next bsml file
sub _reset {
    %gene_hash = ();
    $gene_count = 0;
    $cds_count = 0;
    $polypeptide_count = 0;
    $transcript_count = 0;
    $exon_count = 0;
    $tRNA_count = 0;
    $id = "";
}

sub count_genes {
    my $class = shift;  # class is implicitly stored as the first argument when invoking Class::method
    my $file = shift;
    
    open FH, "<$file" or die "Cannot open $file for reading: $!";
    while (<FH>) {
        chomp;
        $gene_count++ if (/Feature class=\"gene\"/);
    }  
    close FH;
    $gene_count == 0 ? return 1 : return 0;
}

sub count_cds {
    my $class = shift;
    my $file = shift;
    
    open FH, "<$file" or die "Cannot open $file for reading: $!";
    while (<FH>) {
        chomp;
        $cds_count++ if (/Feature class=\"CDS\"/);
    }   
    close FH;
    $cds_count == 0 ? return 1 : return 0;
}

sub count_polypeptides {
    my $class = shift;
    my $file = shift;
    
    open FH, "<$file" or die "Cannot open $file for reading: $!";
    while (<FH>) {
        chomp;
        $polypeptide_count++ if (/Feature class=\"polypeptide\"/);
    }   
    close FH;
    $polypeptide_count == 0 ? return 1 : return 0;
}

sub count_transcripts {
    my $class = shift;
    my $file = shift;
    
    open FH, "<$file" or die "Cannot open $file for reading: $!";
    while (<FH>) {
        chomp;
        $transcript_count++ if (/Feature class=\"transcript\"/);
    }  
    close FH;
    $transcript_count == 0 ? return 1 : return 0;
}

sub count_exons {
    my $class = shift;
    my $file = shift;
    
    open FH, "<$file" or die "Cannot open $file for reading: $!";
    while (<FH>) {
        chomp;
        $exon_count++ if (/Feature class=\"exon\"/);
    }
    close FH;   
    $exon_count == 0 ? return 1 : return 0;
}

sub count_tRNA {
    my $class = shift;
    my $file = shift;
    
    open FH, "<$file" or die "Cannot open $file for reading: $!";
    while (<FH>) {
        chomp;
        $tRNA_count++ if (/Feature class=\"tRNA\"/);
    }
    close FH;   
    $tRNA_count == 0;
    return;
}

sub return_counts {
    my @counts = ($gene_count, $transcript_count, $exon_count, $cds_count, $tRNA_count, $polypeptide_count);
    my @cats = qw(GENE TRANSCRIPT EXON CDS tRNA POLYPEPTIDE);
    my $count_str = "";
    for (my $i = 0; $i < scalar(@counts); $i++) {
	$count_str .= $cats[$i];
	$count_str .= "\t";
	$count_str .= $counts[$i];
	$count_str .= "\n";
    }
    return $count_str;
}

sub divide_by_3 {
    my $class = shift;
    my $file = shift;

    my ($start, $end) = 0;
    my $diff = -1;

    open FH, "<$file" or die "Cannot open $file for reading: $!";
    while (<FH>) {
	chomp;
	if (/id=\"([\w\d]+\.\w+\.\d+.\d+)\"/) {
	    $id = $1;
	}
	if (/startpos=\"(\d+)\" endpos=\"(\d+)\"/) {
	    $start = $1; $end = $2;
	    $diff = $end - $start;
	    $main::logger->warn("WARN, ID: $id has coordinates not divisible by 3\n in file: $file") if ($diff%3 !=0 && $id =~ /transcript/);
	    $id = "";
	}
    }
    return;
}

sub bad_gene_symbols {	#not entirely sure if this is correct but I'm writing based on appearances in the bsml file
    my $class = shift;
    my $file = shift;

    open FH, "<$file" or die "Cannot open $file for reading: $!";
    while (<FH>) {
	chomp;
	if (/id=\"([\w\d]+\.\w+\.\d+.\d+)\"/) {
	    $id = $1;
	}
	unless (/name=\"gene_symbol\" content=\"[a-zA-Z]+\d*\"/) {
	    $main::logger->warn("WARN, ID: $id has a bad gene symbol\n in file: $file") if($id =~ /transcript/);
	    $id = "";
	}
    }
    return;
}

sub duplicate_gene_symbols {
    my $class = shift;
    my $file = shift;

    my $gene_sym = "";

    open FH, "<$file" or die "Cannot open $file for reading: $!";
    while (<FH>) {
	chomp;
	if (/id=\"([\w\d]+\.transcript\.\d+.\d+)\"/) {
	    $id = $1;
	}
	if (/name=\"gene_symbol\" content=\"([a-zA-Z]+\d*)\"/) {
	    $gene_sym = $1;
	    $main::logger->warn("WARN, ID: $id has a duplicate gene symbol ($gene_sym) with an earlier ID\n in file: $file") if ($gene_hash{$gene_sym}++);	# Update gene symbol counts as we see them
	}
    }
    return;
}

sub should_not_have_gs {	#certain gene product names should not contain gene symbols
    my $class = shift;
    my $file = shift;

    my $gene_product ='';

    open FH, "<$file" or die "Cannot open $file for reading: $!";
    while (<FH>) {
	chomp;
	if (/id=\"([\w\d]+\.transcript\.\d+.\d+)\"/) {
	    $id = $1;
	    $gene_product = "";
	}
	if (/\"gene_product_name\" content=\"([a-z]+[^\"]*)\"/) {
	    $gene_product = $1;
	}
	if (/name=\"gene_symbol\" content=\"(\w+\d*)\"/) {
	    $main::logger->warn("WARN, ID: $id should not have a gene symbol as gene product is $gene_product\n in file: $file") if ($gene_product =~ /(hypothetical protein)|(conserved domain protein)/);	
	}
    }
    return;
}

1;
