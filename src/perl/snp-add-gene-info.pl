#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Bio::SeqIO;
use IntervalTree;
use Data::Dumper;
use File::Basename;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
my $iTree = new IntervalTree;
my $seq_cache = {};
my %codons;
my ($P1, $REF_BASE, $QUERY_BASE, $P2, $BUFF, $DIST, $LEN_R, $LEN_Q, $FRM1, $FRM2, $REF_CONFIG, $QUERY_CONFIG) =
    (0 .. 11); #SNPs file columns
my ($GENE_ID, $GENE_START, $GENE_STOP, $POS_IN_GENE, $SYN_NONSYN, $PRODUCT, $GENE_DIRECTION, $REF_CODON, $REF_AMINO_ACID, $QUERY_CODON, $QUERY_AMINO_ACID, $NUM_HOMOPOLYMER) =
    (0 .. 11); # Additional Columns
####################################################

my %options;
my $results = GetOptions (\%options,
			  "reference_genbank|g=s",
			  "show_snps_file|s=s",
			  "query_fasta_list|q=s",
			  "output|o=s",
			  "log|l=s",
			  "debug|d=s",
			  "help|h"
                          );

&check_options(\%options);

# Parse the input fasta files for sequence identifiers and
# cache the path information.
# query_fsa_files:
#   key: sequence identifier
#   val: path to fasta file
my %query_fsa_files;
%query_fsa_files = &map_fasta_inputs( $options{'query_fasta_list'} );

# parse the reference genbank file for gene information
&_log($DEBUG, "Parsing genbank file");
my %reference_genes = &parse_genbank( $options{'reference_genbank'} );

# Build the interval tree. Will be used in searching for location of
# SNPs within genes.
&_log($DEBUG, "Building Tree");
$iTree->buildTree();

# read in the input show_snps_file
my $snpfh;
open($snpfh, "< $options{'show_snps_file'}") or &_log($ERROR, "Can't open file: $options{'show_snps_file'}");
my $outfh;
open($outfh, "> $options{'output'}") or &_log($ERROR, "Can't open file: $options{'output'} for writing");

&parse_and_print_snps( $snpfh, $outfh, \&additional_columns );

close($snpfh);
close($outfh);

print "See $options{'output'}\n";

#p1, ref_base, query_base, p2, buff, dist, len_r, len_q, frm1, frm2, ref_contig, query_contig
sub parse_and_print_snps {
    my ($infh, $outfh, $subref) = @_;
    
    print $outfh join("\t", qw(p1 ref_base query_base p2 buff dist len_r len_q frm1 frm2 ref_contig query_contig) )."\t";
    print $outfh join("\t", qw(gene_id gene_start gene_stop position_in_gene syn_nonsyn product gene_direction ref_codon ref_amino_acid) )."\t";
    print $outfh join("\t", qw(query_codon query_amino_acid num_homopolymer) )."\n";
    
    while( my $line = <$infh> ) {
	chomp($line);
	my @col = split(/\t/, $line);
	my $additional = $subref->(@col);
	if (exists($query_fsa_files{$col[11]})) {
	    my $qcontig = $query_fsa_files{$col[11]};
	    my ($file_base,$file_dir,$file_ext) = fileparse($qcontig,qr/\.[^.]*/);
	    $col[11] = $file_base;
	}
	print $outfh join("\t", @col)."\t$additional\n";
    }
}

# gene_id, gene_start, gene_stop, position_in_gene, syn_nonsyn, product, gene_direction, ref_codon, ref_amino_acid, query_codon, query_amino_acid
sub additional_columns {
    my ($rsnploc, $qsnploc, $refbase, $querybase, $qorient, $ref, $query) = ($_[0], $_[3], $_[1], $_[2], $_[9], $_[10], $_[11]);
    my @retval;

    # Is the SNP within a gene?
    my @ols = $iTree->searchInterval( $rsnploc, $rsnploc );

    # If the SNP is not within a gene, all NA
    return join("\t", qw(intergenic NA NA NA NA NA NA NA NA NA NA) ) if( @ols == 0 );

    my ($gene_id, $gene_start, $gene_stop, $pos_in_gene, $syn, $prod, $strand, $rcodon, $raa, $qcodon, $qaa, $numh) = 
	([],[],[],[],[],[],[],[],[],[],[],[]);
    foreach my $ol ( @ols ) {
	my @ret = &process_overlap( $rsnploc, $qsnploc, $refbase, $ref, $querybase, $query, $qorient, $ol );
	push( @{$gene_id}, $ret[$GENE_ID] );
	push( @{$gene_start}, $ret[$GENE_START] );
	push( @{$gene_stop}, $ret[$GENE_STOP] );
	push( @{$pos_in_gene}, $ret[$POS_IN_GENE] );
	push( @{$syn}, $ret[$SYN_NONSYN] );
	push( @{$prod}, $ret[$PRODUCT] );
	push( @{$strand}, $ret[$GENE_DIRECTION] );
	push( @{$rcodon}, $ret[$REF_CODON] );
	push( @{$raa}, $ret[$REF_AMINO_ACID] );
	push( @{$qcodon}, $ret[$QUERY_CODON] );
	push( @{$qaa}, $ret[$QUERY_AMINO_ACID] );
	push( @{$numh}, $ret[$NUM_HOMOPOLYMER] );
    }
    my @cols;
    my $count = 0;
    map { 
	
	if ( @{$_} == 0 || !defined( $_->[0] ) ) {
	    print "BAD THINGS [$count]\n";
	    die Dumper( ($gene_id, $gene_start, $gene_stop, $pos_in_gene, $syn, $prod, $strand, $rcodon, $raa, $qcodon, $qaa, $numh) );
	    exit(1);
	}
	push(@cols, join("/", @{$_}) );
	$count++;
    } ($gene_id, $gene_start, $gene_stop, $pos_in_gene, $syn, $prod, $strand, $rcodon, $raa, $qcodon, $qaa, $numh);
    return join("\t", @cols);
}


# gene_id, gene_start, gene_stop, position_in_gene, syn_nonsyn, product, gene_direction, ref_codon, ref_amino_acid, query_codon, query_amino_acid
sub process_overlap {
    my ($rsnploc, $qsnploc, $refbase, $ref, $querybase, $query, $qorient, $ol) = @_;

    my $gene_id = $ol->[2];
    &_log($ERROR, "Could not get gene id from searchInterval results")
	unless( $gene_id );
    &_log($ERROR, "Could not find gene_id [$gene_id] in reference_genes")
	unless( exists( $reference_genes{$gene_id} ) );

    # Grab info from parsed reference genbank file 
    my $start = $reference_genes{$gene_id}->{'start'};
    my $end = $reference_genes{$gene_id}->{'end'};
    my $strand = $reference_genes{$gene_id}->{'strand'};


    # Calculate position in gene
    my $pos_in_gene;
    if( $strand == 1 ) {
	$pos_in_gene = $rsnploc - ( $start - 1 );
    } elsif( $strand == -1 ) {
	$pos_in_gene = ( $end + 1 ) - $rsnploc;
    } else {
	&_log($ERROR, "Do not understand strand $strand");
    }

    # Sanity check. Just check the SNP bases which aren't indels
    if( $refbase ne '.' ) {

	# The reference should always be on the forward strand.
	my $seqbase = &get_sequence_from_parent( $ref, $rsnploc, $rsnploc, 1 );
	&_log($ERROR, "Reference nucleotide parsed from fasta does not match SNP predicted in snps file. ".
	      "Possible incorrect fasta provided. [SNP: $rsnploc, ($refbase != $seqbase)]")
	    if( $refbase ne $seqbase );
    }

    if( $querybase ne '.' ) {
	my $seqbase = &get_sequence_from_parent( $query, $qsnploc, $qsnploc, $qorient );
	&_log($ERROR, "Query nucleotide parsed from fasta does not match SNP predicted in snps file. ".
	      "Possible incorrect fasta provided. [SNP: $rsnploc, ($querybase != $seqbase)]")
	    if( $querybase ne uc($seqbase));
    }


    my $numh = "NA";
    if( $refbase eq '.' ) {
	my $tmpquerybase = $querybase;
	$tmpquerybase = &reverse_complement( $tmpquerybase ) if( $qorient == -1 );
	$numh = &get_homopolymer_run( $qsnploc, $query, 2, $tmpquerybase );
    } elsif( $querybase eq '.' ) {
	$numh = &get_homopolymer_run( $rsnploc, $ref, 2, $refbase );
    }

    # At this point, make sure we have a CDS feature. Otherwise, grabbing codons/translations
    # doesn't really make sense. Do the same if this is an indel.
    if( !exists( $reference_genes{$gene_id}->{'CDS'} ) || $refbase eq '.' || $querybase eq '.') {
	my $product = $reference_genes{$gene_id}->{'product'} || "NA";
	return ($gene_id, $start, $end, $pos_in_gene, "NA", $product, 
		$reference_genes{$gene_id}->{'strand'}, "NA", "NA", "NA", "NA", $numh);
    }

    # Make sure everything is okay with the parsed information
    unless( exists( $reference_genes{$gene_id}->{'product'} ) &&
	    defined( $reference_genes{$gene_id}->{'product'} ) ) {
	&_log($ERROR, "Could not find product for gene: $gene_id");
    }		

    # Grab the codons
    my $pos_in_codon = ( ( $pos_in_gene - 1 ) % 3 ); # Will be 0, 1 or 2
    my $refcodon = &get_sequence_region( $reference_genes{$gene_id}->{'seq'},
					 $pos_in_gene - $pos_in_codon,
					 $pos_in_gene + ( 2 - $pos_in_codon ),
					 1 );

    # Calculate query coordinates
    my ($qstart, $qstop);
    if( $reference_genes{$gene_id}->{'strand'} == 1 ) {
	$qstart = $qsnploc - $pos_in_codon;
	$qstop = $qsnploc + 2 - $pos_in_codon;
    } elsif( $reference_genes{$gene_id}->{'strand'} == -1 ) {
	$qstart = $qsnploc - 2 + $pos_in_codon;
	$qstop = $qsnploc + $pos_in_codon;
    } else {
	&_log($ERROR, "Could not understand strand: $reference_genes{$gene_id}->{'strand'}. Expected 1 or -1");
    }

    my $querycodon = &get_sequence_from_parent( $query, $qstart, $qstop, $reference_genes{$gene_id}->{'strand'} );


    if( $refcodon =~ /\./ || $querycodon =~ /\./ ) {
	print "Either reference or query codons contain a '.'. This is not expected.\n";
	print "ref: $ref, query: $query\n";
	print "refcodon: $refcodon, querycodon: $querycodon\n";
	&_log($ERROR, "Ruh roh");
    }

    my $refaa = &translate_codon( $refcodon );
    my $qaa = &translate_codon( $querycodon );
    my $syn = ( $refaa eq $qaa ) ? "SYN" : "NSYN";

    return ($gene_id, $start, $end, $pos_in_gene, $syn, $reference_genes{$gene_id}->{'product'}, 
	    $reference_genes{$gene_id}->{'strand'}, $refcodon, $refaa, $querycodon, $qaa, $numh);
    
}

sub get_homopolymer_run {
    my ( $snploc, $seq, $direction, $base ) = @_;
    my $buffer = 15;
    my $run = 0;
    
    # both directions
    if( $direction == 2 ) {

	$run =  &get_homopolymer_run( $snploc, $seq,  1, $base );
	$run += &get_homopolymer_run( $snploc - 1, $seq, -1, $base );
	
	# going forward
    } elsif( $direction == 1 ) {

	my $region = &get_sequence_from_parent( $seq, $snploc, $snploc + $buffer, 1 );
	my $parsed_base = substr( $region, 0, 1 );
	$run = ( $parsed_base ne $base ) ? $run = 0 : $run = length($1) if( $region =~ /^($parsed_base+)/ );
	$run += &get_homopolymer_run( $snploc + 10, $seq, $direction, $base ) if( $run == $buffer );

	# and backwards
    } elsif( $direction == -1 ) {

	my $region = &get_sequence_from_parent( $seq, $snploc - $buffer, $snploc, 1 );
	my $parsed_base = substr( $region, length($region) - 1, 1 );
	$run = ( $parsed_base ne $base ) ? $run = 0 : $run = length($1) if( $region =~ /($parsed_base+)$/ );
	$run += &get_homopolymer_run( $snploc - $buffer, $seq, $direction, $base ) if( $run == $buffer );
    }

    return $run;
}



# Just need to parse the gebank file for
# genes, coords, function and locus_tag (and also if it's pseudo or not)
# Also build an IntervalTree while parsing to search with later.

sub parse_genbank {
    my ($gbk_file) = @_;
    my %retval;
    my $seq = new Bio::SeqIO( -file => $gbk_file )->next_seq();
    my @features = $seq->get_SeqFeatures();


    my $accession = $seq->display_id;
    $seq_cache->{$accession} = $seq->seq;


    foreach my $feat ( @features ) {
	next if( $feat->primary_tag eq 'source' );
	

	# Get the coordinates from the gene feature
	if( $feat->primary_tag eq 'gene' ) {
	    
	    # Grab the locus_tag
	    my $locus_tag = &get_locus_tag( $feat );

	    # Exit if we've already seen a gene with this locus_tag
	    $retval{$locus_tag} = {} unless( exists( $retval{$locus_tag} ) );
	    &_log($ERROR, "Found multiple genes with same locus_tag: $locus_tag") 
		if( exists( $retval{$locus_tag}->{'gene'} ) && $retval{$locus_tag}->{'gene'} );

	    # Grab the coordinates and strand. Script will die if there is a SplitLocation because we don't handle
	    # those (mostly because I don't think it will happen and I'm lazy. If we get this error, I'll have to
	    # handle them.)
	    my $loc = $feat->location();
	    &_log($ERROR, "This program doesn't handle SplitLocations.") if( $loc->isa('Bio::Location::SplitLocationI') );
	    
	    $retval{$locus_tag}->{'start'} = $loc->start;
	    $retval{$locus_tag}->{'end'} = $loc->end;
	    $retval{$locus_tag}->{'strand'} = $loc->strand;
	    $retval{$locus_tag}->{'gene'} = 1;
	    $retval{$locus_tag}->{'seq'} = &get_sequence_from_parent( $accession, $loc->start, $loc->end, $loc->strand );
	    $retval{$locus_tag}->{'pseudo'} = 1 if( $feat->has_tag('pseudo') );

	} else {
	    
	    # Grab the locus tag
	    my $locus_tag = &get_locus_tag( $feat );

	    # Grab the feature type
	    my $type = $feat->primary_tag;

	    # If this is a repeat region, we want to remove the entry we just made 
	    # for the 'gene'. Why are repeat regions annotated as genes?
	    if( $type eq 'repeat_region' ) {
		delete( $retval{$locus_tag} ) if( exists( $retval{$locus_tag} ) );
		next;
	    }

            # There might already be an entry in retval from gene
	    $retval{$locus_tag} = {} unless( exists( $retval{$locus_tag} ) );

	    # Make sure we haven't seen a feature of this type for this locus tag already
	    &_log($ERROR, "Found multiple $type features with locus_tag: $locus_tag")
		if( exists( $retval{$locus_tag}->{$type} ) && $retval{$locus_tag}->{$type} );

	    # Get the product. If there is no product, use the note.
	    my $product = "No product";
	    if( $feat->has_tag('product') ) {
		$product = join(", ", $feat->get_tag_values( 'product' ) );
	    } elsif( $feat->has_tag('note') ) {
		my @ns = join(", ", $feat->get_tag_values( 'note' ) );
	    }
	    
	    $retval{$locus_tag}->{'product'} = $product;
	    $retval{$locus_tag}->{$type} = 1;	    

	    my $loc = $feat->location();
	    &_log($ERROR, "Start was greater than stop. Didn't think this would happen. [$locus_tag]")
		if( $loc->start > $loc->end );
	    $iTree->addInterval( $locus_tag, $loc->start, $loc->end );

	}
    }

    return %retval;
}

sub get_locus_tag {
    my ($feat) = @_;
    
    unless( $feat->has_tag('locus_tag') ) {
	my $type = $feat->primary_tag;
	print Dumper( $feat );
	&_log($ERROR, "$type feature does not have locus tag");
    }

    # Grab the locus_tag as the primary gene id. Make sure there is one and only one.
    my @lts = $feat->get_tag_values( 'locus_tag' );
    &_log($ERROR, "Could not parse locus_tag from feature in reference genbank file") unless( @lts );
    &_log($ERROR, "Found more than one locus_tag for feature [@lts] in reference genbank file") if( @lts > 1 );
    $lts[0];
}

sub get_sequence_from_parent {
    my ($parent, $start, $stop, $strand) = @_;
    &_log($ERROR, "Start should always be less that stop for get_sequence_region [$start $stop]") if( $start > $stop );
    # Have we parsed this sequence yet?
    unless( exists( $seq_cache->{$parent} ) ) {
	&_log($ERROR, "Could not find query sequence fasta file for sequence $parent") unless( exists( $query_fsa_files{$parent} ) );
	open(IN, "< $query_fsa_files{$parent}") or &_log($ERROR, "Could not open file: $query_fsa_files{$parent} $!");
	my $seq;
	my $first = 1;
	my $head;
	while(<IN>) {
	    next if( /^\s*$/ );
	    chomp;
	    if($_ =~ /^>(\S+)/) {
		unless ($first) {
		    $seq_cache->{$head} = $seq;
		    $seq = '';
		}
		$first = 0;
		$head = $1;
	    } else {
		$seq .= $_ ;
	    }
	}
	$seq_cache->{$head} = $seq if(defined( $head ));
	close(IN);
    }
    &get_sequence_region( $seq_cache->{$parent}, $start, $stop, $strand );
    
}

sub get_sequence_region {
    my ($seq, $start, $stop, $strand) = @_;	
    &_log($ERROR, "Start should always be less that stop for get_sequence_region [$start $stop]")
	if( $start > $stop );
    &_log($ERROR, "Expected strand to be 1 or -1. Do not understand: $strand [get_sequence_region]")
	unless( $strand == 1 || $strand == -1 );

    my $length = ( $stop - $start ) + 1;
    my $subseq = substr( $seq, $start - 1, $length );
    $subseq = &reverse_complement( $subseq ) if( $strand == -1 );
    uc($subseq);

}

sub reverse_complement {
    my ($seq) = @_;
    my $rev = reverse $seq;
    $rev =~ tr/ATCGRYWSKMVBHDatcgrywskmvbhd/TAGCYRWSMKBVDHtagcyrwsmkbvdh/;
    $rev;
}

sub translate_codon {
    my ($codon) = @_;
    &init_codons unless( keys %codons > 0 );
    my $retval = "UNK";
    $retval = $codons{$codon} if( exists( $codons{$codon} ) );
    return $retval;
}

sub map_fasta_inputs {
    my ($list) = @_;
    my %retval;
    open(IN, "<$list") or &_log($ERROR, "Can't open file $list: $!");
    map {
	chomp;
	open(FSA, "< $_") or &_log($ERROR, "Can't open fsa file $_: $!");
	while(my $fl = <FSA>) {
	    chomp($fl);
	    if( $fl =~ /^>(\S+)/ ) {
		&_log($ERROR, "Couldn't parse header from defline from fasta file: $_") unless( $1 );
		$retval{$1} = $_;
	    }
	}
	close(FSA);
    } <IN>;
    return %retval;
}

sub check_options {

    my $opts = shift;

    if( $opts->{'help'} ) {
	&_pod;
    }

    if( $opts->{'log'} ) {
	open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
    }

    foreach my $req ( qw(reference_genbank show_snps_file output query_fasta_list) ) {
	&_log($ERROR, "Option $req is required") unless( $opts->{$req} );
    }
    
}

sub _log {
    my ($level, $msg) = @_;
    my $t = localtime();
    print STDOUT "$msg\n" if( $level <= $debug );
    print $logfh "[$t]: $msg\n" if( defined( $logfh ) );
    exit(1) if( $level == $ERROR );
}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

sub init_codons {
    %codons = ( 
		'AAA'=>'K',
		'AAC'=>'N',
		'AAG'=>'K',
		'AAT'=>'N',
		'ACA'=>'T',
		'ACC'=>'T',
		'ACG'=>'T',
		'ACT'=>'T',
		'AGA'=>'R',
		'AGC'=>'S',
		'AGG'=>'R',
		'AGT'=>'S',
		'ATA'=>'I',
		'ATC'=>'I',
		'ATG'=>'M',
		'ATT'=>'I',
		'CAA'=>'Q',
		'CAC'=>'H',
		'CAG'=>'Q',
		'CAT'=>'H',
		'CCA'=>'P',
		'CCC'=>'P',
		'CCG'=>'P',
		'CCT'=>'P',
		'CGA'=>'R',
		'CGC'=>'R',
		'CGG'=>'R',
		'CGT'=>'R',
		'CTA'=>'L',
		'CTC'=>'L',
		'CTG'=>'L',
		'CTT'=>'L',
		'GAA'=>'E',
		'GAC'=>'D',
		'GAG'=>'E',
		'GAT'=>'D',
		'GCA'=>'A',
		'GCC'=>'A',
		'GCG'=>'A',
		'GCT'=>'A',
		'GGA'=>'G',
		'GGC'=>'G',
		'GGG'=>'G',
		'GGT'=>'G',
		'GTA'=>'V',
		'GTC'=>'V',
		'GTG'=>'V',
		'GTT'=>'V',
		'TAA'=>'Stop',
		'TAC'=>'Y',
		'TAG'=>'Stop',
		'TAT'=>'Y',
		'TCA'=>'S',
		'TCC'=>'S',
		'TCG'=>'S',
		'TCT'=>'S',
		'TGA'=>'Stop',
		'TGC'=>'C',
		'TGG'=>'W',
		'TGT'=>'C',
		'TTA'=>'L',
		'TTC'=>'F',
		'TTG'=>'L',
		'TTT'=>'F'
		);

}


=head1 NAME

 snp-add-gene-info.pl - Adds gene information to show-snps [nucmer] output: intergenic vs genic, ref vs query codons/amino acids, etc.

=head1 SYNOPSIS

  USAGE: snp-add-gene-info.pl
    --reference_genbank=/path/to/ref.gbk
    --show_snps_file=/path/to/file.snps
    --query_fasta_list=/path/to/queries.fsa.list
    --output=/path/to/output.snps
    [ --log=/path/to/file.log
    --debug=3
    --help
    ]

=head1 OPTIONS

B<--reference_genbank,-g>
    Input genbank file for reference sequence.

B<--show_snps_file,-s>
    This is the output file from the show-snps program. Expects tabular output [using the -T options to show-snps].

B<--add_homopolymer,-a>
    This will include homopolymer stretches in indel regions in output tab file. If this is zero, do not need to include the --query_fasta.

B<--query_fasta_list,-q>
    List of fasta files for the query genome(s).

B<--output,-o>
    The output is a tab file. 

B<--log,-l>
    Logfile.

B<--debug,-d>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

    It can be very useful to look at the output of show-snps in the context of annotation. Looking at
    synonomous vs non-synonomous SNPs or indels which may have been created by homopolymer runs [possible
    sequencing issue.]
    
=head1  INPUT

  Required: 
    Referenece genome genbank file.
    show-snps output [tabular output, -T]
    Output file path.
    Query fasta list

=head1 OUTPUT

    The first 12 columns are the same as those output from show-snps [see show-snps documentation for 
    description of columns]:
    p1, ref_base, query_base, p2, buff, dist, len_r, len_q, frm1, frm2, ref_contig, query_contig

    This script will add these additional columns:
    gene_id, gene_start, gene_stop, position_in_gene, snps_per_gene, syn_nonsyn, product, gene_direction, ref_codon, ref_amino_acid, query_codon, query_amino_acid, num_homopolymer

=head1  CONTACT

    Rewrite of Elliot Drabek's process-snps script

    Kevin Galens
    kgalens@gmail.com

=cut
