#!/usr/bin perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Basename;
use FindBin qw($Bin);
use Pod::Usage;

=head1 NAME

format_velvet_assembled.pl - Program to format velvet assembled, annotated genomes to split the contigs at the sequence of N (gaps) 
			     and modify the original .tbl and .fsa files. 

=head1 SYNOPSIS

	format_velvet_assembled.pl --tbl_file <tbl_file> --fsa_file <fsa_file> --split_param <num_of_N> --output_dir <output_directory> --min_contig_len <minimum_contig_length_allowed>

=head1 OPTIONS
	
	--tbl_file <tbl_file>		= /path/to/tbl_file. Path to the original .tbl file generated
				  	  for submission after annotation.  Optional

	--fsa_file <fsa_file>		= /path/to/fsa_file. Path to the original multi FASTA file
				          containing assembled contig sequences with Ns.

	--split_param <num_of_N>	= Integer. Number of Ns to be used for splitting the contigs.
					  Default is 10.

	--output_dir <output_directory>	= /path/to/output_dir. Output directory to store newly formatted
					  files.

	--min_contig_len <min_length>	= Integer. Minimum allowed length of the contig. All contigs less than
					  this specified length and the annotation associated will be deleted.
					  Default is 200bp.

	--log <log_file>		= /path/to/log_file. Log file. Optional

=head1 DESCRIPTION

This program reads in originally created .tbl and .fsa files and formats them to remove a series of Ns (gaps)
from the contig sequences, splitting them further into contigs and adjusts the annotation coordinates accordingly.
It also creates a .agp containing the mapping where contigs were split and number of Ns removed.

=head1  CONTACT

Sonia Agrawal
sagrawal@som.umaryland.edu

=cut

my ($ctbl, $corrected_tbl, $cfsa, $agp, $num_of_n, $tbl_base, $agp_file, $corrected_fsa, $min_contig_len);
my (@old_tbl, @old_fsa);
my $no_table = 0;
my (%options, %contig_hash);
my $results = GetOptions (\%options,
			  "tbl_file|t=s",
			  "fsa_file|f=s",
			  "split_param|n=i",
			  "output_dir|o=s",
			  "min_contig_len|m=i",
			  "log|l=s",
			  "help|h") || pod2usage();

## Display documentation
if( $options{'help'} ){
	pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

&check_options();

$agp_file = $options{'output_dir'}."/".$tbl_base.".agp";
open($agp, "> $agp_file") or die "Could not open $agp_file for writing\n";

&create_contig_hash();

if ((-e $corrected_fsa) && (-s $corrected_fsa)) {
	system("$Bin/clean_fasta.pl $corrected_fsa");
} else {
	print STDERR "Not able to clean parse newly created fasta file $corrected_fsa\n";
}

unless ($no_table) {
	&edit_tbl_file($options{'tbl_file'});

	if ((-e $corrected_tbl) && (-s $corrected_tbl)) {
		my $ctbl_cleaned = $corrected_tbl.".tmp";
		system("perl $Bin/format_velvet_assembled_AGP_cleanup.pl -t $corrected_tbl -o $ctbl_cleaned");
		system("mv $ctbl_cleaned $corrected_tbl");
	}
}
print STDERR "\n";

## Subroutine to check the options passed at the command line
sub check_options {
	if (!($options{'output_dir'})) {
		die "ERROR : Output directory is required\n";	
	}
	if ($options{'tbl_file'}) {
		my ($tbl_dir, $tbl_ext);
		($tbl_base, $tbl_dir, $tbl_ext) = fileparse($options{'tbl_file'}, qr/\.[^.]*/);
		$corrected_tbl = $options{'output_dir'}."/".$tbl_base."_corrected".$tbl_ext;
		open($ctbl, "> $corrected_tbl") or die "Could not open $corrected_tbl for writing\n";
	} else {
		print STDERR "No table file provided.  Skipping table adjustment.\n";
		$no_table = 1;
	}

	if ($options{'fsa_file'}) {
		my ($fsa_base, $fsa_dir, $fsa_ext);
		($fsa_base, $fsa_dir, $fsa_ext) = fileparse($options{'fsa_file'}, qr/\.[^.]*/);
		$corrected_fsa = $options{'output_dir'}."/".$fsa_base."_corrected".$fsa_ext;
		open($cfsa, "> $corrected_fsa") or die "Could not open $corrected_fsa for writing\n";
	} else {
		die "ERROR : fsa_file is required\n";
	}

	if ($options{'split_param'}) {
		$num_of_n = $options{'split_param'};
	} else {
		$num_of_n = 10;
	}
	if ($options{'min_contig_len'}) {
		$min_contig_len = $options{'min_contig_len'};
	} else {
		$min_contig_len = 200;
	}
} 


## Subroutine to alter .tbl file by adjusting the coordinates of the contigs after they are split
sub edit_tbl_file {
	my ($tbl_file) = @_;
	my ($header, $tbl);
	my ($feat_start, $feat_stop, $feat, $id, $desc);
	my $counter = 0;
	my %contig_genes;
	open (FH, "< $tbl_file") or die "Could not open $tbl_file file for reading. $!\n";
	while(<FH>) {
		$tbl = $_;
		chomp($tbl);
		next if ($tbl =~ /^\s*$/);
		if($tbl =~ /^>Feature\s+(\S+)/) {
			$header = $1;
			print STDERR "\nProcessing Feature $header .....\n";
			$counter = 0;
			%contig_genes = ();
			next;
		} 
		# If the contig was not split just print the features in it
		if ((!(exists $contig_hash{$header}{'split_num'}))) {
			if (!(exists($contig_hash{$header}{'seq_len'}))) {
				print $ctbl ">Feature $header"."_1\n";
				print $ctbl "$tbl\n";
				while (<FH>) {
					$tbl = $_;
					chomp($tbl);
					if ($tbl =~ /^>Feature\s+(\S+)/) {
						$header = $1;
						print STDERR "\nProcessing Feature $header .....\n";
						$counter = 0;
						%contig_genes = ();
						last;
					}
					next if ($tbl =~ /^\s*$/);
					print $ctbl "$tbl\n";
				}
			}
		} else {
			if ($tbl =~ /^(<*\d+)\s+(>*\d+)\s+(\w+)/) {
					$feat_start = $1;
					$feat_stop = $2;
					$feat = $3;
					if ($feat =~ /gene/) {
						$counter++;
#						print STDERR "\t$counter";
						$contig_genes{$counter}{'start'} = $feat_start;
						$contig_genes{$counter}{'stop'} = $feat_stop;
					}
					$contig_genes{$counter}{$feat} = "";
			} elsif ($tbl =~ /^\s+(\w+)\s+(.+)/) {
					$id = $1;
					$desc = $2;
					$contig_genes{$counter}{$feat} .= $id.":::".$desc.";";
			}
			while (<FH>) {
				$tbl = $_; 
				chomp($tbl);
				if ($tbl =~ /^>Feature\s+(\S+)/) {
#						print STDERR "\n";
						# Process each contig to adjust the coordinates according to the split 
						&process_contig($header,\%contig_genes);
						$header = $1;
						print STDERR "\nProcessing Feature $header .....\n";
						$counter = 0;
						%contig_genes = (); 
						last;
				}
				next if ($tbl =~ /^\s*$/);
				if ($tbl =~ /^(<*\d+)\s+(>*\d+)\s+(\w+)/) {
					$feat_start = $1;
					$feat_stop = $2;
					$feat = $3;
					if ($feat =~ /gene/) {
						$counter++;
#						print STDERR "\t$counter";
						$contig_genes{$counter}{'start'} = $feat_start;
						$contig_genes{$counter}{'stop'} = $feat_stop;
					}
					$contig_genes{$counter}{$feat} = "";
				} elsif ($tbl =~ /^\s+(\w+)\s+(.+)/) {
					$id = $1;
					$desc = $2;
					$contig_genes{$counter}{$feat} .= $id.":::".$desc.";";
				}
			}
		}
	}
# Processing last contig's annotation
	if ((exists $contig_hash{$header}{'split_num'})) {
#		print STDERR "\n";
		&process_contig($header,\%contig_genes);
	}
	close(FH);
}

## Subroutine to adjust a contig's annotation coordinates when it is split  
sub process_contig {
	my ($head, $split_hash) = @_;
	my %gene_counter = ();
	my %gene_index = ();
	my (@gene, @tags, @sorted_gene_counter);
	my ($gene_start, $gene_stop, $feat, $id, $desc, $start, $stop, $locus);
	my ($count, $split_head, $diff, $subpart);
	my ($partial_start, $partial_stop);

	foreach $count (sort {$a <=> $b} keys %{$split_hash}) {
		$gene_start = $split_hash->{$count}{'start'};
		if ($gene_start =~ /^<(\d+)/) {
			$gene_start = $1;
		}
		$gene_stop = $split_hash->{$count}{'stop'};
		if ($gene_stop =~ /^>(\d+)/) {
			$gene_stop = $1;
		}
		
		print STDERR "\tGene $count:$gene_start-$gene_stop"; 
		
		foreach $split_head (keys %{$contig_hash{$head}}) {
			next if (!($split_head =~ /^$head/));
			my $temp_start = $gene_start;
			my $temp_stop = $gene_stop;
			if ($gene_start > $gene_stop) {
                                $temp_start = $gene_stop;
                                $temp_stop = $gene_start;
                        }
#			print STDERR "($temp_start <= $contig_hash{$head}{$split_head}->{'end_pos'}) && ($contig_hash{$head}{$split_head}->{'start_pos'} <= $temp_stop)\n";
			if (($temp_start <= $contig_hash{$head}{$split_head}->{'end_pos'}) && ($contig_hash{$head}{$split_head}->{'start_pos'} <= $temp_stop)) {
#				print STDERR "(($contig_hash{$head}{$split_head}->{'start_pos'} <= $gene_start) && ($gene_start <= $contig_hash{$head}{$split_head}->{'end_pos'})) OR $contig_hash{$head}{$split_head}->{'start_pos'} <= $gene_stop AND $gene_stop <= $contig_hash{$head}{$split_head}->{'end_pos'}\n";
				## Modification made by Sonia Agrawal on 2012/05/25 to split a gene when it falls in two split contigs
#				if ((($contig_hash{$head}{$split_head}->{'start_pos'} <= $gene_start) && ($gene_start <= $contig_hash{$head}{$split_head}->{'end_pos'})) || (($contig_hash{$head}{$split_head}->{'start_pos'} <= $gene_stop) && ($gene_stop <= $contig_hash{$head}{$split_head}->{'end_pos'}))) {
					if (!(defined $gene_index{$split_head})) {
						$gene_index{$split_head} = "";
					}
					$gene_index{$split_head} .= "$count,";
					## Modification made by Sonia Agrawal on 2012/06/05 to assign a unique locus_tag to split genes when they fall in two or more split contigs
					# Storing occurence of gene in various contigs
					$split_head =~ m/^$head\_(\d+)/;
					if(exists($gene_counter{$count})) {
						push(@{$gene_counter{$count}}, $1);
					} else {
						$gene_counter{$count}[0] = $1;
					}
					print STDERR " --> $split_head .....";
#					last;
#				}
			}
		}
		print STDERR "\n";
	}
	for (my $i=1;$i<=$contig_hash{$head}{'split_num'};$i++) {
		$split_head = $head."_".$i;
		if ($contig_hash{$head}{$split_head}->{'length'} < $min_contig_len) {
			next;
		}
		print STDERR "\n\tWriting $split_head .....";
		print $ctbl ">Feature $split_head\n";
		next if(!defined($gene_index{$split_head}));
		@gene = split(/,/, $gene_index{$split_head});
		for (my $j=0;$j<@gene;$j++) {
			$count = $gene[$j];
			$partial_start = $partial_stop = "";
			$gene_start = $split_hash->{$count}{'start'};
			if ($gene_start =~ /^<(\d+)/) {
				$gene_start = $1;
				$partial_start = "<";
			}
			$gene_stop = $split_hash->{$count}{'stop'};
			if ($gene_stop =~ /^>(\d+)/) {
				$gene_stop = $1;
				$partial_stop = ">";
			}
			$start = $gene_start - $contig_hash{$head}{$split_head}->{'start_pos'} + 1;
			$stop = $gene_stop - $contig_hash{$head}{$split_head}->{'start_pos'} + 1;
			## Modification made by Sonia Agrawal on 2012/05/25 to split a gene when it falls in two split contigs
			if ($gene_start < $gene_stop) {
				if ($contig_hash{$head}{$split_head}->{'end_pos'} < $gene_stop) {
					$stop = $contig_hash{$head}{$split_head}->{'length'};
					$partial_stop = ">";
				} 
				if ($gene_start < $contig_hash{$head}{$split_head}->{'start_pos'}) {
					$partial_start = "<";
					$start = 1;
				}
			} else {
				if ($gene_stop < $contig_hash{$head}{$split_head}->{'start_pos'}) {
					$stop = 1;
					$partial_stop = ">";
				} 
				if($contig_hash{$head}{$split_head}->{'end_pos'} < $gene_start) {
					$partial_start = "<";
					$start = $contig_hash{$head}{$split_head}->{'length'};
				}
			}
			$subpart = "";
			@sorted_gene_counter = ();
			my $flag = 0;
			my ($t, $alpha, $codon_start, $frame, $frag_len, $gl_flag);
			$codon_start = 1;
			# Convert ASCII code to capital alphabetic character for genes that are split across contigs
			if(@{$gene_counter{$count}} > 1) {
				@sorted_gene_counter = sort {$a <=> $b} @{$gene_counter{$count}} if($start < $stop);
				@sorted_gene_counter = sort {$b <=> $a} @{$gene_counter{$count}} if($stop < $start);
				@sorted_gene_counter = map { $head."_".$_ } @sorted_gene_counter;
				for($t=0 ; $t<@sorted_gene_counter ; $t++) {
					if($sorted_gene_counter[$t] eq $split_head) {
						$alpha = $t;
						last;
					}
				}
				$subpart = chr ($alpha + 65);
			}
#				if($alpha > 0) {			
					$codon_start = $1 if((defined($split_hash->{$count}{'CDS'}) && $split_hash->{$count}{'CDS'} =~ m/codon_start:::(\d)/) || (defined($split_hash->{$count}{'tRNA'}) && $split_hash->{$count}{'tRNA'} =~ m/codon_start:::(\d)/) || (defined($split_hash->{$count}{'rRNA'}) && $split_hash->{$count}{'rRNA'} =~ m/codon_start:::(\d)/));		
#					print STDERR "$split_hash->{$count}{'CDS'}\t$split_head\t$contig_hash{$head}{$split_head}->{'start_pos'}\t$gene_start\t$gene_stop\t$codon_start\t";
					$frag_len = ($contig_hash{$head}{$split_head}->{'start_pos'} - $gene_start + 1) - ($codon_start - 1) if($gene_start < $gene_stop);
					$frag_len = ($gene_start - $contig_hash{$head}{$split_head}->{'end_pos'} + 1) - ($codon_start - 1) if($gene_stop < $gene_start);
					$frame = $frag_len%3;
					if($frame == 0) {
						$codon_start = 2;
					} elsif($frame == 2) {
						$codon_start = 3;
					} else {
						$codon_start = 1;
					}
#					print STDERR "$frag_len\t$codon_start\n";
#				}
#			}
			$gl_flag = check_gene_length($start, $stop, $partial_stop, $codon_start);
			next if($gl_flag == 1);
			print $ctbl $partial_start.$start."\t".$partial_stop.$stop."\t"."gene\n";
			@tags = split(/;/, $split_hash->{$count}{'gene'});
			my $locus_tag;
			my $attrib = "gene";
			for (my $i=0;$i<@tags;$i++) {
				($id, $locus) = split(/:::/, $tags[$i]);
				$locus_tag = $locus if($id eq "locus_tag");
				print $ctbl "\t\t\t$id\t$locus";
				print $ctbl "$subpart" if($id eq "locus_tag");
				print $ctbl "\n";
			}
			foreach $feat (sort keys %{$split_hash->{$count}}) {
				next if (($feat =~ /start/) || ($feat =~ /stop/) || ($feat =~ /gene/));
				print $ctbl $partial_start.$start."\t".$partial_stop.$stop."\t"."$feat\n";
				@tags = split(/;/, $split_hash->{$count}{$feat});
				for (my $i=0;$i<@tags;$i++) {
					($id, $desc) = split(/:::/, $tags[$i]);
					next if($id eq "codon_start" && @{$gene_counter{$count}} > 1 && $alpha > 0);
					$flag = 1 if($id eq "codon_start"); 
					print $ctbl "\t\t\t$id\t$desc";
					print $ctbl "$subpart" if($id eq "protein_id");
					print $ctbl "\n";
				}
				print $ctbl "\t\t\tcodon_start\t$codon_start\n" if((@{$gene_counter{$count}} > 1 && length($partial_start) > 0 && ($alpha > 0)) || (length($partial_start) > 0 && $flag != 1));		
				$attrib = $feat if($feat =~ /RNA$/);	
				print $ctbl "\t\t\tnote\t5' end; 3' end is $attrib $locus_tag on contig ".$sorted_gene_counter[-1]."\n" if(length($partial_stop) > 0 && @{$gene_counter{$count}} > 1 && $sorted_gene_counter[-1] ne $split_head && $sorted_gene_counter[0] eq $split_head );
				print $ctbl "\t\t\tnote\t3' end; 5' end is $attrib $locus_tag on contig $sorted_gene_counter[0]\n" if(length($partial_start) > 0 && @{$gene_counter{$count}} > 1 && $sorted_gene_counter[0] ne $split_head && $sorted_gene_counter[-1] eq $split_head);
			}
		}
	}
#	print STDERR "\n";
}

# Subroutine to check the length of a gene after split. Minimum length of a gene in 3 bases to form a codon if it is partial on 3'. 
# Minimum length is 6 if it is complete on 3' to find a stop codon and another codon before stop
sub check_gene_length {
	my ($gs, $gsp, $ps, $cs) = @_;
	my $flag = 0;
	my $glen = ($gsp - $gs + 1) if($gsp >= $gs);
	$glen = ($gs - $gsp + 1) if($gs > $gsp);
# If partial on 3' end
	if(length($ps) > 0 && $glen < ($cs + 2)) {
		$flag = 1;
	}
# If complete on 3' end
	if(length($ps) == 0 && $glen < ($cs + 5)) {
		$flag = 1;
	}	
	return($flag);
}


## Subroutine to read a multi FASTA file and create the altered multi FASTA file by splitting the contigs by a series of Ns
## It also creates a .agp file to store the mapping of the split contigs and gap fragments for the original contigs
sub create_contig_hash {
	my $flag = 1; 
	my $seq = "";
	my $header;
	@old_fsa = read_file($options{'fsa_file'});

	foreach my $line (@old_fsa) {
                chomp($line);
                if($line =~ /^>(.*)/) {
                        unless($flag) {
                                chomp($seq);
                                my $seq_len = length($seq);
#				print "$header\t$seq_len\n";
				my $annot_flag = &check_contig_len($header, $seq_len);
				if ($annot_flag == 1) {
					my $frag_size = &split_contig($seq, $header);
					if ($frag_size == 1) { 
						print $agp "$header\t1\t$seq_len\t1\tW\t".$header."_1\t1\t$seq_len\t+\n";
#                                		$contig_hash{$header}{'seq_len'} = $seq_len;
						print $cfsa ">$header"."_1\n$seq\n";
					}
				} else {
					$contig_hash{$header}{'seq_len'} = $seq_len;
				}
                                $seq = ""; 
                        }
                        $flag = 0;
                        $header = $1; 
                } else {
                        # skip it if it is just whitespace
                        next if ($line =~ /^\s*$/);
                        $seq .= $line;
                }
        }
# Concatenate the last contig sequence
        chomp($seq);
	my $seq_len = length($seq);
	if (defined($header)) {
		my $annot_flag = &check_contig_len($header, $seq_len);
		if ($annot_flag == 1) {
			my $frag_size = &split_contig($seq, $header);
			if ($frag_size == 1) {
				print $agp "$header\t1\t$seq_len\t1\tW\t".$header."_1\t1\t$seq_len\t+\n";
#			$contig_hash{$header}{'sequence'} = $seq;
				print $cfsa ">$header"."_1\n$seq\n";
			}
		} else {
			$contig_hash{$header}{'seq_len'} = $seq_len;
		}
	}
}

sub check_contig_len {
	my ($head, $seq_length) = @_;
#	my $tbl_head;
	my $annot_yes = 1;
	if ($seq_length < $min_contig_len) {
		$annot_yes = 0;
#		$tbl_head = `grep -o -w $head $options{'tbl_file'}`;
#		chomp($tbl_head);
#		$tbl_head =~ s/^\s+|\s+$//;
#		if ($tbl_head eq $head) {
#			$annot_yes = 1;
#		} else {
		print STDERR "Removing $head from FASTA file as it is less than $min_contig_len bp in length\n";
#		}
	}
	return($annot_yes);
}



## Subroutine to split a sequence by a series of Ns and return the number of fragments generated
sub split_contig {
	my ($sequence, $head) = @_;
	my @print_contigs = ();
	my @frags = split(/(N{$num_of_n,})/i, $sequence);
	my $frags_size = scalar(@frags);
	if ($frags_size > 1) {
		my $i = 1;
		my $start_pos = 1;
		my $end_pos;
		foreach my $fr (@frags) {
			$end_pos = $start_pos + length($fr) - 1;
			if($fr !~ /^N+$/i) {
				my $split_head = $head."_".$i;
				my $split_len = length($fr);
				my $contig_frag = {
					'start_pos' => $start_pos,
					'end_pos' => $end_pos,
					'length' => $split_len,
				};
				$contig_hash{$head}{$split_head} = $contig_frag;
				my $annot_flag = &check_contig_len($split_head, $split_len);
				if ($annot_flag == 1) {		
					print $cfsa ">$split_head\n$fr\n";
					push(@print_contigs,$split_head);
				}
				$i++;
			}
			$start_pos += length($fr);
		}
		$contig_hash{$head}{'split_num'} = ($i - 1);
		$start_pos = 1;
		my $part_num = 1;
		my $j;
		my $sub_offset = 0;
## Modifications made by Sonia Agrawal on May 29, 2012 to remove gaps from the beginning and end of a split scaffold in the .agp file and always start the first contig in the split scaffold with 1
		for($j=0;$j < @print_contigs;$j++) {
			if(($contig_hash{$head}{$print_contigs[$j]}{'start_pos'} > $start_pos) && ($start_pos != 1)) {
				$end_pos = $contig_hash{$head}{$print_contigs[$j]}{'start_pos'} - $sub_offset - 1;
				print $agp "$head\t$start_pos\t$end_pos\t$part_num\tN\t".($end_pos - $start_pos + 1)."\tscaffold\tyes\tpaired-ends\n";
				$part_num++;  
			} else {
				$sub_offset = $contig_hash{$head}{$print_contigs[$j]}{'start_pos'} - 1;
			} 
			print $agp "$head\t".($contig_hash{$head}{$print_contigs[$j]}{'start_pos'} - $sub_offset)."\t".($contig_hash{$head}{$print_contigs[$j]}{'end_pos'} - $sub_offset)."\t$part_num\tW\t$print_contigs[$j]\t1\t$contig_hash{$head}{$print_contigs[$j]}{'length'}\t+\n";
			$part_num++;
			$start_pos = $contig_hash{$head}{$print_contigs[$j]}{'end_pos'} - $sub_offset + 1;
		}
#		if($contig_hash{$head}{$print_contigs[$j-1]}{'end_pos'} < length($sequence)) {
#			print $agp "$head\t$start_pos\t".length($sequence)."\t$part_num\tN\t".(length($sequence) - $start_pos + 1)."\tscaffold\tyes\tpaired-ends\n";
#		}
	}
	return($frags_size);
}

## Subroutine to read files
sub read_file {
        my $filename = shift;
        my @lines;
        open(FH , "< $filename")  || die("Could not open $filename file for reading.$!");
        @lines = <FH>;
        close(FH);
        return(@lines);
}

__END__ 
