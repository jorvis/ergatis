#!/usr/bin/env perl

=head1 NAME

	snp-verify.pl - Description

	=head1 SYNOPSIS

  USAGE: snp-verify.pl
	--input_file=/path/to/some/input.file
	--output=/path/to/transterm.file
	[ --log=/path/to/file.log
	--debug=3
	--help
	]

	=head1 OPTIONS

	B<--input_file,-i>

	B<--output_file,-o>

	B<--log,-l>
    Logfile.

	B<--debug,-d>
    1,2 or 3. Higher values more verbose.

	B<--help,-h>
    Print this message

	=head1  DESCRIPTION

	DESCRIPTION
	
	=head1  INPUT
    Describe the input

	=head1 OUTPUT
    Describe the output

	=head1  CONTACT

    Kevin Galens
    kgalens@gmail.com

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Data::Dumper;
use File::Basename;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
my %aa;
####################################################

my %options;
my $results = GetOptions (\%options,
						  "input_list|i=s",
						  "output_dir|o=s",
						  "log|l=s",
						  "debug|d=s",
						  "help|h"
                          );

&check_options();

my $coords = &parse_coords_file( $options{'coords_file'} );
my $queries = &parse_queries_from_table( $options{'blast_list'} );
my $snps = &parse_all_blast_results( $options{'blast_list'}, $queries);

## open the merged table and the output file.
open(TBL, "< $options{'merged_table'}") or die("Could not open $options{'merged_table'}: $!");
my ($gb_base,$gb_dir,$gb_ext) = fileparse($options{'coords_file'},qr/\.[^.]*/);
my $res_file = $options{'output_dir'}."/merged_table_".$gb_base.".verified";
my $outfh;
open($outfh, "> $res_file") or die("Could not open $options{'output_file'} for writing: $!");

## Print the header for the output file
my $header = &make_header($queries);
print $outfh $header."\n";

while( my $line = <TBL> ) {
	next if( $line =~ /^\s*$/ || $line =~ /^refpos/ );
	chomp $line;

	my @c = split(/\t/, $line );
	my ($refpos, $gene, $syn, $refbase, $product) = 
		($c[0], $c[1], $c[2], $c[3], $c[4]);
	$refpos =~ s/^\s+|\s+$//;
	$gene =~ s/\s+$//;
	$gene =~ s/^\s+//;
	unless( exists( $snps->{$refpos} ) ) {
		print STDERR "No query blast hits at SNP location $refpos\n";
	}

	print $outfh "$refpos\t$gene\t$syn\t$refbase\t";

	my $num_hits = "";
	if( !exists( $snps->{$refpos} ) ) {
		my $num_hits = "";
		map { 
			print $outfh "No Hit\t";
			$num_hits .= "0\t";
		} @{$queries};
	} else {
		foreach my $q ( @{$queries} ) {
			my $snp = &get_snp_from_hits( $snps->{$refpos}, $q );
			die("SNP was blank [$snp], $refpos, $q") 
				if( !defined( $snp ) || $snp =~ /^\s+$/ );
			print $outfh "$snp\t";
			
			if( exists( $snps->{$refpos}->{$q} ) ) {
				my $num = scalar(@{$snps->{$refpos}->{$q}->{'hits'}});
				$num = ($num == 1 ) ? $num : "$num";
				$num_hits .= "$num\t";
			} else {
				$num_hits .= "0\t";
			}

		}
	}

	print $outfh "$product\t$num_hits";

	if( $gene ne 'None' && exists( $snps->{$refpos} ) ) {
		my $codons = &get_codon_info( $snps->{$refpos}, $refpos, $refbase, $gene );
		if($codons != 0) { 
			print $outfh "$codons->{'gene_length'}\t";
			print $outfh "$codons->{'pos_in_gene'}\t";
			print $outfh "$codons->{'ref_codon'}\t";
			print $outfh "$codons->{'ref_aa'}\t";
			print $outfh "$codons->{'query_codon'}\t";
			print $outfh "$codons->{'query_aa'}";
		} else {
			print $outfh "NA\tNA\tNA\tNA\tNA\tNA";
		}
			
	} else {
		print $outfh "NA\tNA\tNA\tNA\tNA\tNA";
	}

	if( $gene ne 'None' && exists( $coords->{$gene} ) ) {
		print $outfh "\t$coords->{$gene}->{'start'}\t$coords->{$gene}->{'stop'}";
	}

	print $outfh "\n";
	
}

close(TBL);
close($outfh);

##################################################
#    subroutines                                 #
##################################################


sub aa{
    my ($codon) = @_;

    if (scalar(keys %aa) == 0){
	$aa{"ATT"} = "I";
	$aa{"ATC"} = "I";
	$aa{"ATA"} = "I";
	$aa{"CTT"} = "L";
	$aa{"CTC"} = "L";
	$aa{"CTA"} = "L";
	$aa{"CTG"} = "L";
	$aa{"TTA"} = "L";
	$aa{"TTG"} = "L";
	$aa{"GTT"} = "V";
	$aa{"GTC"} = "V";
	$aa{"GTA"} = "V";
	$aa{"GTG"} = "V";
	$aa{"TTT"} = "F";
	$aa{"TTC"} = "F";
	$aa{"ATG"} = "M";
	$aa{"TGT"} = "C";
	$aa{"TGC"} = "C";
	$aa{"GCT"} = "A";
	$aa{"GCC"} = "A";
	$aa{"GCA"} = "A";
	$aa{"GCG"} = "A";
	$aa{"GGT"} = "G";
	$aa{"GGC"} = "G";
	$aa{"GGA"} = "G";
	$aa{"GGG"} = "G";
	$aa{"CCT"} = "P";
	$aa{"CCC"} = "P";
	$aa{"CCA"} = "P";
	$aa{"CCG"} = "P";
	$aa{"ACT"} = "T";
	$aa{"ACC"} = "T";
	$aa{"ACA"} = "T";
	$aa{"ACG"} = "T";
	$aa{"TCT"} = "S";
	$aa{"TCC"} = "S";
	$aa{"TCA"} = "S";
	$aa{"TCG"} = "S";
	$aa{"AGT"} = "S";
	$aa{"AGC"} = "S";
	$aa{"TAT"} = "Y";
	$aa{"TAC"} = "Y";
	$aa{"TGG"} = "W";
	$aa{"CAA"} = "Q";
	$aa{"CAG"} = "Q";
	$aa{"AAT"} = "N";
	$aa{"AAC"} = "N";
	$aa{"CAT"} = "H";
	$aa{"CAC"} = "H";
	$aa{"GAA"} = "E";
	$aa{"GAG"} = "E";
	$aa{"GAT"} = "D";
	$aa{"GAC"} = "D";
	$aa{"AAA"} = "K";
	$aa{"AAG"} = "K";
	$aa{"CGT"} = "R";
	$aa{"CGC"} = "R";
	$aa{"CGA"} = "R";
	$aa{"CGG"} = "R";
	$aa{"AGA"} = "R";
	$aa{"AGG"} = "R";
	$aa{"TAA"} = "*";
	$aa{"TAG"} = "*";
	$aa{"TGA"} = "*";
    }
    if (($codon =~ /N/) || ($codon =~ /X/)){
		return "X";
    }
    if (exists($aa{$codon})){
		return $aa{$codon};
    }
    if (length($codon) == 3){
		return "?";
    }
	die("BAD CODON");
}
sub get_codon_info {
	my ($hits, $refpos, $refbase, $gene) =@_;
	my ($gene_length, $pos_in_gene, $ref_codon, $ref_aa, $query_codon, $query_aa);
	my $retval;

	if( $gene eq 'None' ) {
		$retval = {
			'gene_length' => 'NA',
			'pos_in_gene' => 'NA',
			'ref_codon' => 'NA',
			'ref_aa' => 'NA',
			'query_codon' => 'NA',
			'query_aa' => 'NA'
			};
		return $retval;
	}
	
	## Does the gene exist?
#	die("Could not find gene $gene in coord file") unless( exists( $coords->{$gene} ) );
	unless( exists( $coords->{$gene} ) ) {
		print "Could not find gene $gene in coords file\n";
		return 0;
	} 
	$retval->{'gene_length'} = abs( $coords->{$gene}->{'start'} - $coords->{$gene}->{'stop'} ) + 1;

	## find the position within the gene.
	my $forward = ( $coords->{$gene}->{'start'} < $coords->{$gene}->{'stop'} ) ? 1 : 0;
	if( $forward ) {
		$retval->{'pos_in_gene'} = $refpos - $coords->{$gene}->{'start'} + 1;
	} else {
		$retval->{'pos_in_gene'} =  $coords->{$gene}->{'start'} - $refpos + 1;
	}

	my $pos_in_codon = (($retval->{'pos_in_gene'} + 1) % 3) + 1;;
	if( !$forward ) {
		if( $pos_in_codon == 3 ) {
			$pos_in_codon = 1;
		} elsif( $pos_in_codon == 1 ) {
			$pos_in_codon = 3;
		}
	} 
	my ($rcodon, $qcodon);
  QUERY:
	foreach my $q ( keys %{$hits} ) {
	  HIT:
		foreach my $hit ( @{$hits->{$q}->{'hits'}} ) {
			next HIT if( $hit->{'indel'} );

			my $qline = $hit->{'lines'}->[0];
			my $sline = $hit->{'lines'}->[1];

			my ($qstart, $qseq, $qend) = ($1, $2, $3) if( $qline =~ /Query:\s*(\d+)\s*([\w\-]+)\s*(\d+)/ );
			my ($sstart, $sseq, $send) = ($1, $2, $3) if( $sline =~ /Sbjct:\s*(\d+)\s*([\w\-]+)\s*(\d+)/ );
			next HIT if( $qstart > 21 || $qend < 21 );

			my $offset = 20 - ($qstart - 1);
			$offset = $offset - ($pos_in_codon - 1);

			my $refcodon = substr( $qseq, $offset, 3 );
			my $query_codon = substr( $sseq, $offset, 3 );

			if( !$forward ) {
				$refcodon = &reverse_complement( $refcodon );
				$query_codon = &reverse_complement( $query_codon );
			}

			next HIT if( length($refcodon) != 3 || length($query_codon) != 3 );

			$rcodon = $refcodon unless( defined( $rcodon ) );
			$qcodon = $query_codon if( !defined( $qcodon ) || $query_codon ne $rcodon );

			if( $rcodon ne $refcodon ) {
				print Dumper( $hits );
				print "Rcodon: $rcodon vs $refcodon\n";
				print "REFBASE: $refbase\n";
				die("Found two different reference codons");
			} elsif( $qcodon ne $query_codon && $query_codon ne $rcodon ) {
				print Dumper( $hits );
				print "SNP POS: $refpos\n";
				print "QUERY: $q\n";
				print "Refbase: $refbase, $rcodon\n";
				print "$qcodon vs $query_codon\n";
				die("Found two different query codons");
			}
			
		} # HIT
	} # QUERY

	$retval->{'ref_codon'} = "NA";
	$retval->{'query_codon'} = "NA";
	$retval->{'ref_aa'} = "NA";
	$retval->{'query_aa'} = "NA";
	if( defined( $rcodon ) ) {
		$retval->{'ref_codon'} = uc($rcodon);
		$retval->{'ref_aa'} = &aa(uc($rcodon));
	}
	if( defined( $qcodon ) ) {
		$retval->{'query_codon'} = uc($qcodon);
		$retval->{'query_aa'} = &aa(uc($qcodon));
	}
	return $retval;
}
sub reverse_complement {
	my ($codon) = @_;
	my $rev = reverse($codon);
	$rev =~ tr/acgt/tgca/;
	return $rev;
}

sub get_snp_from_hits {
	my ($data, $q) = @_;
	my $snp_string = "No Hit";
	if( exists( $data->{$q} ) ) {
		my %uniq;
		foreach my $hit ( @{$data->{$q}->{'hits'}} ) { 
			$uniq{$hit->{'snp'}} = 1;
		}
		$snp_string = join("", keys %uniq);
	}
	return $snp_string;
	
}
sub make_header {
	my ($queries) = @_;
	my $header = join("\t", ('refpos','gene','syn?','refbase',@{$queries},'product'));
	map { $header .= "\tnum_hits:$_" } @{$queries};
	$header .= "\t".join("\t", ('gene_length', 'pos_in_gene', 'ref_codon', 'ref_aa', 'query_codon', 'query_aa') );
	return $header;
}
sub parse_all_blast_results {
	my ($blast_list, $queries) = @_;
	print "Queries : @$queries\n";
	## make sure all the files are there
	my %files;
	open(FR, "$blast_list") or die "Could not open file $blast_list for reading.\nReason : $!\n";
	while(<FR>) { 
		chomp($_);
		$_ =~ s/^\s+//;
		$_ =~ s/\s+$//;
		my $bf = $_;
		if( -e $bf ) {
			my ($bf_base,$bf_dir,$bf_ext) = fileparse($bf,qr/\.[^.]*/);
			$bf_base =~ /(.+)_refmol_/;
			my $query_name = $1;
			foreach my $q (@{$queries}) {
				chomp($q);
				$q =~ s/^\s+//;
				$q =~ s/\s+$//;
				if($q eq $query_name) {
					$files{$q} = $bf;
					last;
				} 
			}	
		} else {
			die("Blast file does not exist");
		}
	} 

	my $results = {};
	foreach my $query ( keys %files ) {
		print "Parsing $query\n";
		$results = &parse_blast_file( $files{$query}, $query, $results );
	}

	return $results;
}

sub parse_blast_file {
	my ($blast_file, $query, $data) = @_;

	open(IN, "< $blast_file") or die("Couldn't open $blast_file: $!");

	## DEBUG
	my $flag = 0;

	my $count = 0;
	my $total = 0;

	my $snp_pos;
	my $lines;
	while( my $line = <IN> ) {
		next if( $line =~ /^\s*$/ );
		chomp($line);

		if( $line =~ /Query= SNP_(\d+)/ ) {
			## Store the snp position
			$snp_pos = $1;
		} elsif( $line =~ /Sbjct:/ ) {
			## Push the current line onto the list of hits
			push(@{$lines}, $line);

			## If we haven't stored any hits for this snp yet, 
			## create an array ref to store them
			$data->{$snp_pos}->{$query}->{'hits'} = [] unless( exists( $data->{$snp_pos}->{$query}->{'hits'} ) );
			my $parsed_hit = &parse_blast_hit( $lines, $flag );
			push( @{$data->{$snp_pos}->{$query}->{'hits'}}, $parsed_hit ) if( defined( $parsed_hit ) );

			if( $flag ) {
				print Dumper( $lines );
				print Dumper( $data->{$snp_pos}->{$query} );
				exit(0);
			}

			## Reset the array which holds the hit
			$lines = [];

			## DEBUG
			$count++;
			last if( $total != 0 && $count == $total );

		} elsif( $line =~ /^[\s\|]+$/ ) {
			## Skip the line which contains the vertical bars
			next;
		} else {
			## Keep capturing the hit line. Pretty sure this will only grab the
			## Query: 1 agata ... line.
			push(@{$lines}, $line);
		}
	}

	close(IN);
	return $data;
}

sub parse_blast_hit {
	my ($lines, $flag) = @_;
	my $retval = {};
	
	# find pos to parse from query
	my $query_line = $lines->[0];
	my ($qstart,$seq,$qend) = ($1,$2,$3) if( $query_line =~ /^Query: (\d+)\s+([\w\-]+)\s+(\d+)$/ );
	die("Could not parse query start or query end from query line $query_line [start: $qstart, end: $qend]")
		unless( defined( $qstart ) && defined( $qend ) );

	## Considered a No Hit if qstart doesn't contain the SNP
	## (i.e. the qstart is greater than 21)
	return if( $qstart > 21 || $qend < 21 );

	## calculate the offset to use in the substr routine
	my $offset = 20 - ($qstart - 1);
	print "Offset: $offset\n" if $flag;
	
	## Grab the subject line
	my $sub_line = $lines->[1];
	
	## parse start, stop and the sequence
	my ($sstart,$sseq,$send) = ($1,$2,$3) if( $sub_line =~ /^Sbjct: (\d+)\s+([\w\-]+)\s+(\d+)$/ );
	die("Could not parse subject start or end from subject line $sub_line [$sstart, $send]") 
		unless( defined( $sstart ) && defined( $send ) );

	## grab the nucleotide in the SNP position
	$retval->{'snp'} = uc(substr( $sseq, $offset, 1 ));
	## is it a partial match?
	$retval->{'partial'} = ( $qstart > 1 || $qend < 41 ) ? 1 : 0;

	## is it an indel?
	$retval->{'indel'} = 0;
	$retval->{'indel'} = 1 if( $sseq =~ /-/ || $seq =~ /-/ );

	$retval->{'lines'} = $lines;
	
	return $retval;
}

sub parse_queries_from_table {
	my ($table) = @_;
	my $queries = [];

	open( IN, "< $table") or die("Couldn't open $table: $!");
	while(<IN>) {
		my $line = $_;
		chomp($line);
		next if ($line =~ /^\s*$/);
		next if ($line =~ /^#/);
		if (-e $line) {
			my ($file_base,$file_dir,$file_ext) = fileparse($line,qr/\.[^.]*/);
			$file_base =~ /(.+)_refmol_/;
			push (@{$queries}, $1);
		}		
	}
	close(IN);

#	my @c = split(/\t/, $line);
#	@{$queries} = @c[4..(scalar(@c)-8)];
#	print " @$queries\n";
	return $queries;
}

sub parse_coords_file {
	my ($file) = @_;
	my $retval = {};

	my ($locus, $fivep, $threep, $common_name) = (0..3);

	open(COORDS, "< $file" ) or die("couldn't open file $file for reading: $!");
	while( my $line = <COORDS> ) {
		chomp($line);

		my @c = split(/\t/, $line );
		die("Found repeat in coords file for gene $c[$locus]")
			if( exists( $retval->{$c[$locus]} ) );

		$retval->{$c[$locus]} = {
			'start' => $c[$fivep],
			'stop'  => $c[$threep],
			'name'  => $c[$common_name]
			};

	}
	close(COORDS);

	return $retval;
}

sub check_options {
	if( $options{'help'} ) {
		&_pod;
	}

	if( $options{'log'} ) {
		open( $logfh, "> $options{'log'}") or die("Can't open log file ($!)");
	}

	foreach my $req ( qw(input_list output_dir) ) {
		&_log($ERROR, "Option $req is required") unless( $options{$req} );
	}

	if (-e $options{'input_list'}) {
		my @ifiles;
		open(FLIST, "$options{'input_list'}") or &_log($ERROR, "Could not open file $options{'input_list'} for reading.\nReason : $!");
		my @file_list = <FLIST>;
		foreach my $line (@file_list) {
			chomp($line);
			next if ($line =~ /^\s*$/);
			next if ($line =~ /^#/);
			push @ifiles, $line;
		}
		close(FLIST);
		if(-e $ifiles[0]) {
			$options{'coords_file'} = $ifiles[0];
		} else {
			&_log($ERROR, "Coords file $ifiles[0] does not exist\n");
		}
		if (-e $ifiles[1]) {
			$options{'merged_table'} = $ifiles[1];
		} else {
			&_log($ERROR, "SNP positions file $ifiles[1] does not exist\n");
		}
		if (-e $ifiles[2]) {
			$options{'blast_list'} = $ifiles[2];
		} else {
			&_log($ERROR, "Parsed BLAST list file $ifiles[2] does not exist\n");
		}
	} else {
		&_log($ERROR, "Input list file $options{'input_list'} does not exist\n");
	}
}

sub _log {
	my ($level, $msg) = @_;
	if( $level <= $debug ) {
		print STDOUT "$msg\n";
	}
	print $logfh "$msg\n" if( defined( $logfh ) );
	exit(1) if( $level == $ERROR );
}

sub _pod {
	pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
