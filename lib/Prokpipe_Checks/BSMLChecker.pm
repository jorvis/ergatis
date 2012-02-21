#!/usr/bin/perl -w

# BSMLChecker.pm

# This module is designed to run consistency checks on the BSML output from the 
# pipeline_summary component of the Prokaryotic Annotation Pipeline.  If a check
# fails, either a warning will be written to the Ergatis prokpipe_consistency_checks
# component log or an error will be written and kill off the pipeline.

package Prokpipe_Checks::BSMLChecker;

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
my $warn_flag = 0;
my $length_cutoff = 1000;


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

# Determines length of assembly.  If less than 1000 bp, then gene counts are ignored
sub get_length {
    my $class = shift;
    my $file = shift;
    my $length = -1;
    my $title;

    open FH, "<$file" or die "Cannot open $file for reading: $!";
    while (<FH>) {
	chomp;
	if (/Sequence length=\"(\d+)\" class=\"\w+\" title=\"(.+)\"\s/) {
	    $length = $1;
	    $title = $2;
	}
    }
    if ($length < $length_cutoff) {
	$main::logger->warn( "WARN: $title is less than 1000 bp (length: $length) in file $file\n");
	return 1;
    }
    return 0;
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
	    $warn_flag = 1;
	}
    }
    close FH;
    $warn_flag ? return 0 : return 1;
}

sub bad_gene_symbols {	#not entirely sure if this is correct but I'm writing based on appearances in the bsml file
    my $class = shift;
    my $file = shift;

    open FH, "<$file" or die "Cannot open $file for reading: $!";
    while (<FH>) {
	chomp;
	if (/id=\"([a-z]{3}\w*)\"/) {
	    $id = $1;
	}
	unless (/name=\"gene_symbol\" content=\"[a-zA-Z]+\d*\"/) {
	    $main::logger->warn("WARN, ID: $id has a bad gene symbol\n in file: $file") if($id =~ /transcript/);
	    $id = "";
	    $warn_flag = 1;
	}
    }
    close FH;
    $warn_flag ? return 0 : return 1;
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
	    $warn_flag = 1;
	}
    }
    close FH;
    $warn_flag ? return 0 : return 1;
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
	if (/\"gene_product_name\" content=\"([A-Za-z]+[^\"]*)\"/) {
	    $gene_product = $1;
	}
	if (/name=\"gene_symbol\" content=\"(\w+\d*)\"/) {
	    $main::logger->warn("WARN, ID: $id should not have a gene symbol as gene product is $gene_product\n in file: $file") if ($gene_product =~ /(hypothetical protein)|(conserved domain protein)|(family protein)|(putative)/);
	    $warn_flag = 1;	
	}
    }
    close FH;
    $warn_flag ? return 0 : return 1;
}

sub ec_check {	# All gene product names except hypothetical proteins and conserved domain proteins should have EC numbers
    my $class = shift;
    my $file = shift;

    my $gene_product = "";
    my $ec = "";

    open FH, "<$file" or die "Cannot open $file for reading: $!";
    while (<FH>) {
	chomp;
	if (/id=\"([\w\d]+\.transcript\.\d+.\d+)\"/) {
	    $id = $1;
	    $gene_product = "";
	    $ec = "";
	}
	if (/\"gene_product_name\" content=\"([A-Za-z]+[^\"]*)\"/) {
	    $gene_product = $1;
	}
	if (/\"EC\" content=\"(([\d-]+\.){3}[\d-]+)\"/) {
	    $ec = $1;
	    if ($gene_product =~ /(hypothetical protein)|(conserved domain protein)|(family protein)|(putative)/){
	    	$main::logger->warn("WARN, ID: $id should not contain an EC number as gene product is $gene_product\n in file: $file") if ($ec ne "");
	    	$warn_flag = 1;
	    }
	    else {
	    	$main::logger->warn("WARN, ID: $id should contain an EC number as gene product is $gene_product\n in file: $file") if ($ec eq "");
	    	$warn_flag = 1;
	    }
	}

    }
    close FH;
    $warn_flag ? return 0 : return 1;
}

sub TIGR_role_check { #every gene should have a TIGR Role (hypothetical proteins will not be looked at)
    my $class = shift;
    my $file = shift;

    my $gene_product = "";
    my $tigr = "";

    open FH, "<$file" or die "Cannot open $file for reading: $!";
    while (<FH>) {
	chomp;
	if (/id=\"([\w\d]+\.transcript\.\d+.\d+)\"/) {
	    $id = $1;
	    $gene_product = "";
	    $tigr = "";
	}
	if (/\"gene_product_name\" content=\"([A-Za-z]+[^\"]*)\"/) {
	    $gene_product = $1;
	}
	if (/\"TIGR_role\" content=\"(\d+)\"/) {
	    $tigr = $1;
	    if ($tigr eq "" && $gene_product !~ /hypothetical protein/) {
		$main::logger->warn("WARN, ID: $id should contain a TIGR role as gene product is $gene_product\n in file: $file");
		$warn_flag = 1;
	    }
	}

    }
    close FH;
    $warn_flag ? return 0 : return 1;
}

sub valid_start_stop_codons {

    my $class = shift;
    my $file = shift;

    my $path = "";
    my ($header, $sequence);

    open FH, "<$file" or die "Cannot open $file for reading: $!";
    while (<FH>) {
	chomp;
	if (/id=\"([\w\d]+\.\CDS\.\d+.\d+)\"/) {
	    $id = $1;
	}
	if (/source=\"([\/\w]+\.fsa)\"/) {
	    $path = $1;

	    open PATH, "<$path" or $main::logger->warn("WARN, Cannot open $path for reading...see BSML2Fasta component"); 
	    while (<PATH>) {
		chomp;
            	if ($line =~ /^>/) {
                    my $next_header = substr($line, 1);
                    $header = $next_header;
                    $sequence = '';
            	} else {
                    next if ($line =~ /^\s*$/);  #ignore whitespace
                    $sequence .= uc($line);
            	}
	    }
	    close PATH;
	    $path = "";
	    if (substr ($sequence,0,3) ne "ATG" ||
		substr ($sequence,0,3) ne "GTG" ||
		substr ($sequence,0,3) ne "TTG") {
		$main::logger->warn("WARN, Sequence for ID $id does not start with a valid start codon\n see file $path");
		$warn_flag = 1;
	    }

	    if (substr ($sequence,-3,0) ne "TGA" ||
		substr ($sequence,-3,0) ne "TAA" ||
		substr ($sequence,-3,0) ne "TAG") {
		$main::logger->warn("WARN, Sequence for ID $id does not end with a valid stop codon\n see file $path");
		$warn_flag = 1;
	    }
	}
    }
    close FH;
    $warn_flag ? return 0 : return 1;
}

1;
