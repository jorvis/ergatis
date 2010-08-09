package MG::RunParse;

use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use Data::Dumper;
use XML::Simple;
use MG::ManFasta;
use MG::Common;
use Bio::SearchIO;
use Bio::Seq::RichSeq;
use MG::ParseSeq;
use MG::Math;
@ISA = qw(Exporter);
@EXPORT = qw(&run_search &parse_search &run_orthomcl &parse_orthomcl &runparse_lalign &runparse_ssearch &activesite_ssearch &parse_metagene &remove_overlaps &runparse_showsnp &runparse_showcoords &runparse_fastx &create_db &runparse_showtiling &parse_hmmscan);

my %global = MG::Common::initialize();

sub parse_search {
    my $algo = shift @_;
    my $report = shift @_;
    my $eval_min = shift @_;
    my $dna2prot = shift @_;
    my $rev = shift @_;
    my $one_hsp = shift @_;
    my %bhits = ();
    my %trans = ();
    my $sio = new Bio::SearchIO(-format => $algo,
				-file   => $report);
    
 W1:while (my $result = $sio->next_result) {
	my $query_length = $result->query_length;
	my $query_name = $result->query_name;
    W2:while (my $hit = $result->next_hit) {
	    my $hit_name = $hit->name; #could also be $hit->description
	    my $per_qcov = $hit->frac_aligned_query if ($result->query_length);
	    my $per_hcov = $hit->frac_aligned_hit if ($hit->length);
	    my $hit_length = $hit->length;
	    my $hit_desc = $hit->description;
	W3:while (my $hsp = $hit->next_hsp) {
		my $eval = $hsp->significance();
		$eval =~ s/,//g;
		$eval = "1".$eval if $eval =~ m/^e/;
		my $bits = $hsp->bits;
		my $alen = $hsp->length("hit");
		my $percid = sprintf("%.2f", $hsp->frac_identical);
		my $percsim = sprintf("%.2f", $hsp->frac_conserved);
		my $qbegin = $hsp->start('query');
		my $qend = $hsp->end('query');
		my $lbegin = $hsp->start('hit');
		my $lend = $hsp->end('hit');	
		next W2 unless ($eval < $eval_min);
		if ($dna2prot) {
		    my $strand = '+';
		    $strand = '-' if ($hsp->strand == -1);
		    if ($qbegin > $qend) {
			$strand = '-';
			my $t = $qbegin;
			$qbegin = $qend;
			$qend = $t;
		    }
		    my $query_seq =  $hsp->seq_str;
		    $query_seq =~ s/-|\\|\///g;
		    $query_seq =~ s/\*/X/g;
		    my $gap = $hsp->cigar_string;
		    $trans{$query_name}{$strand}{$hit_name} .=  $query_seq;
		    push @{$bhits{$query_name}{$strand}}, [$hit_name, $eval, $percid, $strand, 
							   $qbegin, $qend, $query_length, $lbegin, 
							   $lend, $hit_length];
		}elsif ($rev) {
		    push @{$bhits{$hit_name}}, [$query_name, $eval, $bits, $alen, $percid, $percsim, 
						$per_qcov, $per_hcov, $query_length, $hit_length, 
						$lbegin,$lend, $hit_desc,$qbegin, $qend,];
		}else {
		    push @{$bhits{$query_name}}, [$hit_name, $eval, $bits, $alen, $percid, $percsim, 
						  $per_qcov, $per_hcov, $query_length, $hit_length, 
						  $lbegin,$lend, $hit_desc,$qbegin, $qend];
		}
		next W1 if $one_hsp;
	    }
	}
    }
    if ($dna2prot) {
	return(\%bhits, \%trans);
    }else {
	return %bhits;
    }
}
sub runparse_showsnp {
    my $deltafile = shift @_;
    my %snps;
    system("show-snps -I -T -H $deltafile > $deltafile\.snps") unless (-e "$deltafile\.snps");
    open IN, "<$deltafile\.snps" or die $!;
    while (my $line = <IN>) {
	chomp($line);
	my($chrpos,$orint,$newnt,$readpos,$buff,$dist,$r,$q,
	   $from1,$from2,$chr,$read) = split(/\s+/, $line);
	next if ($orint eq '.' || $newnt eq '.');
	$snps{$chr}{$read}{$chrpos}{$newnt} ++;
    }
    return %snps;
}
sub runparse_showtiling {
    my $deltafile = shift @_;
    my $minpercid = shift @_;
    $minpercid = $minpercid *100 if ($minpercid < 1);
    my $mincov = shift @_;
    $mincov = $mincov *100 if ($mincov < 1);
    my %tiles = ();
    system("show-tiling -c -i $minpercid -V 0 -R -v $mincov $deltafile > $delta.tiling") unless (-e "$deltafile\.tiling");
    open IN, "<$deltafile\.tiling" or die $!;
    $/ = "\n>";
    while (my $record = <IN>) {
	chomp($record);
	$record =~ s/^>//g;
	my ($ref_info, @reads) = split(/\n/, $record);
	$ref_info =~ m/(.+)\s+(\d+)\s+bases/;
	my ($header, $ref_len) = ($1,$2);
	my $ref_name = (split(/\s+/, $header))[0];
	foreach $line (@reads) {
	    my ($ref_start,$ref_end,$gap,$q_len,$q_cov,
		$q_percid,$strand,$q_name) = split(/\s+/,$line);
	    push @{$tiles{$ref_name}}, [$q_name,$q_percid,$q_cov,$q_len,
					$strand,$ref_len,$ref_start,$ref_end];
	}
    }
    $/ = "\n";
    return %tiles;
}

sub runparse_showcoords {
    my $deltafile = shift @_;
    my $minpercid = shift @_;
    my $mincov = shift @_;
    my $revname = shift @_;
    $minpercid = $minpercid *100 if ($minpercid < 1);
    $mincov = $mincov *100 if ($mincov < 1);
    my %coords;
    system("/usr/local/bin/show-coords -H -l -I $minpercid -c -T $deltafile > $deltafile\.coords") unless (-e "$deltafile\.coords");
    open IN, "<$deltafile\.coords" or die $!;
    while (my $line = <IN>) {
	chomp($line);
	my($ref_begin, $ref_end, $qbegin, $qend, $ref_alen, $qalen, $percid,
	   $rlen,$qlen, $rcov, $qcov, $ref_name, $query_name) = split(/\s+/, $line);
	next if ($percid < $minpercid || $qcov < $mincov);
	$perc_qcov = sprintf("%.2f",$qcov/100);
	$perc_hcov = sprintf("%.2f",$rcov/100);
	$percid = sprintf("%.2f",$percid/100);
	if ($revname) {
	    push @{$coords{$ref_name}}, [$query_name,$qalen, $percid,$perc_qcov,$perc_hcov,$qlen,$rlen, 
					 $ref_begin,$ref_end,$qbegin, $qend];
	}else {
	    push @{$coords{$query_name}}, [$ref_name,$qalen,$percid,$perc_qcov,$perc_hcov,$qlen,$rlen, 
					   $ref_begin,$ref_end,$qbegin, $qend];
	}
    }
    return %coords;
}
sub remove_overlaps {
    my @ali = sort {$a->[4] <=> $b->[4]} @{shift @_};
    my $j = 0;
    my %limits = ();
    my %hash = ();
 F1:foreach $seg (@ali) {
	my ($hit,$eval, $percid, $strand, $qbegin, $qend, $qlen, $lbegin, $lend, $llen) = @{$seg};
    F2:foreach my $i (keys %limits) {
	    my ($bhid, $min_eval, $min_qbegin, $max_qend) = @{$limits{$i}};
	    if ($qbegin < ($max_qend - 20)) {
		if ($eval < $min_eval) {
		    @{$hash{$i}} = ($hit,$eval, $percid, $strand, $qbegin, $qend, $qlen, $lbegin, $lend, $llen);
		    @{$limits{$i}} = ($hit, $eval,$qbegin,$qend);
		}
		next F1;
	    }
	}
	$j ++;
	@{$hash{$j}} = ($hit,$eval, $percid, $strand, $qbegin, $qend, $qlen, $lbegin, $lend, $llen);
	@{$limits{$j}} = ($hit, $eval,$qbegin,$qend);
    }
    return values %hash;
}
sub create_db {
    my ($outfile,$seqfile,@acc) = @_;
    my %h = ();
    foreach (@acc) {
	$h{$_} = 1;
    }
    open IN, "<$seqfile" or die $!;
    open OUT, ">$outfile" or die $!;
    $/ = "\n>";
    while (my $record = <IN>) {
	chomp($record);
	my ($header, @seq) = split(/\n/, $record);
	$header =~ s/^>//g;
	my $acc = (split(/\s+/, $header))[0];
	if ($h{$acc}) {
	    print OUT ">".$header."\n".join("",@seq)."\n";
	}
    }
    close OUT;
    close IN;
    return 1;
}
sub runparse_lalign {
    my ($ref, $test) = @_;
    system("$global{local_bin}\/lalign35 -I $ref $test > $global{tmp}\/lalign.out");
    open IN, "<$global{tmp}\/lalign.out" or die $!;
    my $hit;
    my %hash = ();
    while (my $line = <IN>) {
	if ($line =~ m/^>>(\w\S+)/) {
	    $hit = $1;
	}elsif ($line =~ m/([\d\.]+%\s+identity).+\((\d+)-(\d+):(\d+)-(\d+)\)/) {
	    my @array = ($1,$2,$3,$4,$5);
	    $hash{$hit} = \@array;
	}
    }
    return %hash;
}
sub runparse_ssearch {
    my ($ref, $test) = @_;
    system("$global{local_bin}\/ssearch35 -m9c -d 0 -f 16 -H -g 4 $ref $test > $global{tmp}\/lalign.out");
    my $file = `cat $global{tmp}\/lalign.out`;
    open IN, "<$global{tmp}\/lalign.out" or die $!;
    my $hit;
    $/ = "\n";
    my %hash = ();
    while (my $line = <IN>) {
	if ($line=~ m/^The best scores are:/) {
	INNER:while ($line = <IN>) {
		last INNER if ($line =~ m/^\s*$/o );
		chomp($line);
		my $idnote;
		my $desc = (split(/\t/,$line))[0];
		my $scores = (split(/\t/,$line))[1];
		my $alignment = (split(/\t/,$line))[2];
		my $hit = (split(/\s+/, $desc))[0];
		my $percid = 100 * (split(/\s+/, $scores))[0];
		my $alen = (split(/\s+/, $scores))[3];
		my $qbegin = (split(/\s+/, $scores))[4];
		my $qend = (split(/\s+/, $scores))[5];
		my $lbegin = (split(/\s+/, $scores))[8];
		my $lend = (split(/\s+/, $scores))[9];
		my $qlen = (split(/\s+/, $scores))[7];
		my $llen = (split(/\s+/, $scores))[11];
		my $aa_gaps =(split(/\s+/, $scores))[12]+(split(/\s+/, $scores))[13];
		my $gaps=0;
		my $indels=0;
		$gaps++ if($aa_gaps);
		foreach ($alignment =~ m/\+|-/g) {
		    $indels ++;
		}
		my $match = int((split(/\s+/, $scores))[0] * $alen);
		$fpercid = $match/($alen-$aa_gaps);
		my $fake_percid = sprintf("%.2f",100*$fpercid);
		if ($llen > $qlen) {
		    $idnote = "Warning: $llen aa ";
		}elsif ($fake_percid <= 95) {
		    $idnote = "Warning: $fake_percid ";
		}
		$idnote .= "$fake_percid% identity";
		if ($gaps) {
		    $idnote .= " with $indels gaps"
		}
		my @array = ($idnote,$qbegin,$qend,$lbegin,$lend);
		$hash{$hit} = \@array;
	    }
	}
    }
    return %hash;
}
sub runparse_fastx {
    my ($ref, $test,$outname,$done) = @_;
    my ($report,$sfile)= @{$done};
    unless ($done) {
	system("mv $ref $global{scratch}\/$outname\.query");
	system("mv $test $global{scratch}\/$outname\.lib");
	system("$global{local_bin}\/fastx36 -m 9c -H $global{scratch}\/$outname\.query $global{scratch}\/$outname\.lib > $global{scratch}\/$outname\.fxout");
	
	$report = "$global{scratch}\/$outname\.fxout";
	$sfile  = "$global{scratch}\/$outname\.query";
    }
    my ($nt, $oname) = parsefasta($sfile);
    
    my %seq = %{$nt};
    foreach (keys %seq) {
	$acc = $seq{$_}{db_acc};
	$seq{$acc} = $seq{$_};
	delete $seq{$_};
    }
    my %seen;
    my %seqobj;
    my $j = 0;
    my $in = Bio::SearchIO->new(-file  => $report,
				-format => 'fasta');
 W1:while (my $r = $in->next_result) {
	my $qname = $r->query_name;
	my $i = 0;
	unless ($seq{$r->query_name}{sequence}) {
	    warn "Trouble in Denmark: Seqname in Search and Fasta file do not match!\n";
	    next W1;
	}
	$newseq = Bio::Seq::RichSeq->new(-seq=>$seq{$r->query_name}{sequence},
					 -accession_number=>$qname,-id=>$qname,
					 -primary_id=>$qname,-display_id=>$qname);
	$newseq->desc($qname);
	my @tax_class = ('human metagenome','organismal metagenomes,metagenomes','unclassified sequences');
	my $species = Bio::Species->new(-ncbi_taxid=>646099,-classification => \@tax_class);
	$newseq->species($species);
	$newseq->add_date(ncbidate());
	my $feat = Bio::SeqFeature::Generic->new(-start=>1,-end=>length($seq{$r->query_name}{sequence}),
						 -primary=>'source',
						 -tag=>{organism=>'human metagenome',
							mol_type=>'genomic DNA',
							isolation_source=>'Homo sapiens',
							collection_date=>ncbidate()});
	$newseq->add_SeqFeature($feat) unless (length($seq{$r->query_name}{sequence}) < 20);
    W2:while (my $h = $r->next_hit) {
	W3:while (my $hsp = $h->next_hsp) {
		if ($seen{$qname}{$hsp->strand}) {
		    foreach (@{$seen{$qname}{$hsp->strand}}) {
			next W2 if ($hsp->start('query') >= $_->{begin} && $hsp->start('query') <= $_->{end});
			next W2 if ($hsp->end('query') >= $_->{begin} && $hsp->end('query') <= $_->{end});
			next W2 if ($_->{begin} >= $hsp->start('query') && $_->{begin} <= $hsp->end('query'));
			next W2 if ($_->{end} >= $hsp->start('query') &&  $_->{end} <= $hsp->end('query'));

		    }
		}
		my $prot_seq = $hsp->seq_str;
		my $location = new Bio::Location::Split->new();
		if ($hsp->strand eq -1) {
		    my $start = $hsp->end('query');
		    my $offset = 0;
		    my $len;
		    while ($prot_seq ne '') {
			if ($prot_seq =~ m/^(\w+)/) {
			    $len = length($1);
			    $plen = 3*$len - $offset;
			    $location->add_sub_Location(Bio::Location::Simple->new(-end=>$start,
										   -start=>$start - $plen +1,
										   -strand=>$hsp->strand
										  ));
			    $start = $start - $plen;
			    $offset = 0;
			}elsif ($prot_seq =~ m/^(\*+)/) {
			    $len = length($1);
			    $plen = 3*$len - $offset;
			    $start = $start - $plen;
			    $offset = 0;
			}elsif ($prot_seq =~ m/^(\/+)/) {
			    $len = length($1);
			    $offset = $len;
			}elsif ($prot_seq =~ m/^(\\+)/) {
			    $len = length($1);
			    $start -= $len;
			}elsif ($prot_seq =~ m/(-+)/) {
			    $len = length($1);
			}else {
			    warn "Forgot about this Character: $prot_seq";
			}
			substr($prot_seq,0,$len,'');
		    }
		} else {
		    my $start = $hsp->start('query');
		    my $offset = 0;
		    my $len;
		    while ($prot_seq ne '') {
			if ($prot_seq =~ m/^(\w+)/) {
			    $len = length($1);
			    $plen = 3*length($1) - $offset;
			    $location->add_sub_Location(Bio::Location::Simple->new(-start=>$start,
										   -end=>$start + $plen - 1,
										   -strand=>$hsp->strand
										  ));
			    $start = $start + $plen;
			    $offset = 0;
			}elsif ($prot_seq =~ m/^(\*+)/) {
			    $len = length($1);
			    $plen = 3*length($1) - $offset;
			    $start = $start + $plen;
			    $offset = 0;
			}elsif ($prot_seq =~ m/^(\/+)/) {
			    $len = length($1);
			    $offset = $len;
			}elsif ($prot_seq =~ m/^(\\+)/) {
			    $len = length($1);
			    $start += $len;
			}elsif ($prot_seq =~ m/(-+)/) {
			    $len = length($1);
			}else {
			    warn "Forgot about this Character: $prot_seq";
			}
			substr($prot_seq,0,$len,'');
		    }
		}
		my $fasta_pseq = $hsp->seq_str;
		$fasta_pseq =~ s/\\|\/|-//g;
		$i++;
		my $did = $qname.'|fastx|gene_'.$i;
		my $feat = Bio::SeqFeature::Generic->new(-start=>$hsp->start('query'),-end=>$hsp->end('query'),
							 -strand=>$hsp->strand,-primary=>'CDS',-source_tag=>'fastx',
							 -display_name=>$did,-score=>$hsp->bits,
							 -frame=>0,-tag=>{translation=>$fasta_pseq,
									  organism=>'human metagenome',
									  mol_type=>'genomic DNA',
									  locus_tag=>$did,
									  product=>'hypotetical protein',
									  isolation_source=>'Homo sapiens',
									  best_hit=>$h->name,
									  collection_date=>ncbidate(),
									 },
							);
		$feat->location($location);
		$newseq->add_SeqFeature($feat);
		push @{$seen{$qname}{$hsp->strand}}, {begin=>$hsp->start('query'),end=>$hsp->end('query')};
	    }
	}
	$seqobj{$r->query_name} = $newseq;
    }
    system("rm $global{scratch}\/$outname.fxout $global{scratch}\/$outname\.query $global{scratch}\/$outname\.lib");
    return %seqobj;
}
sub activesite_ssearch {
    my ($ref, $test) = @_;
    system("$global{local_bin}\/ssearch35 -m9c -E 1e-6 -H -d 0 -V '&*@#%' $ref $test > $global{tmp}\/lalign.out");
    warn "ssearch completed\n";
    my $file = `cat $global{tmp}\/lalign.out`;
    open IN, "<$global{tmp}\/lalign.out" or die $!;
    my $hit;
    $/ = "\n";
    my %hash = ();
    my $annotation;
    while (my $line = <IN>) {
	if ($line=~ m/Annotation:\s*(\S+)/) {
	    $annotation = $1;
	}elsif ($line=~ m/^The best scores are:/) {
	INNER:while ($line = <IN>) {
		last INNER if ($line =~ m/^\s*$/o );
		chomp($line);
		
		my $desc = (split(/\t/,$line))[0];
		my $scores = (split(/\t/,$line))[1];
		my $alignment = (split(/\t/,$line))[2];
		my $activesite_aln = (split(/\t/,$line))[3];
		
		my $hit = (split(/\s+\(/, $desc))[0];
		my $percid = 100 * (split(/\s+/, $scores))[0];
		my $alen = (split(/\s+/, $scores))[3];
		my $eval = (split(/\s+/, $desc))[-1];
		my $qbegin = (split(/\s+/, $scores))[4];
		my $qend = (split(/\s+/, $scores))[5];
		my $lbegin = (split(/\s+/, $scores))[8];
		my $lend = (split(/\s+/, $scores))[9];
		my $qlen = (split(/\s+/, $scores))[7];
		my $llen = (split(/\s+/, $scores))[11];
		
		my @array = ($eval,$percid,$qbegin,$qend,$lbegin,$lend,$alignment,$activesite_aln);
		$hash{$hit} = \@array;
	    }
	}
    }
    return ($annotation, \%hash);
}
sub parse_hmmscan {
    my $in = shift @_;
    $/ = "\n//\n";
    my %results;
    open IN, "<$in" or die $!;
 W1:while (my $record = <IN>) {
	chomp($record);
	next W1 if ($record =~ m/No hits detected/);
	my ($query_info, @hits) = split(/\n\s*>>\s*/, $record);
	$query_info =~ m/Query:\s+(\S+)\s+\[L=(\d+)\]/;
	my ($query,$qlen) = ($1,$2);
    F1:foreach $hit (@hits) {
	    my ($hit_descript, @lines) = split(/\n/, $hit);
	    my ($hit_acc, @hit_desc) = split(/\s+/, $hit_descript);
	    my $hit_desc = join(" ", @hit_desc);
	F2:foreach $line (@lines) {
		next F2 unless ($line =~ m/^\s+(\d+)\s+\S+\s+(.+)/);
		my ($num, $scores) = ($1,$2);
		my ($bit, $bias, $eval1, $eval2,$hmmstart, $hmmend,$sym,$qstart,$qend,$sym2) = split(/\s+/, $scores);
		next F2 unless ($eval1 <= 1e-5 and $eval2 <= 1e-5);
		push @{$results{$query}{$hit_acc}}, [$hit_desc,$bit,$eval1,$eval2,$qstart,$qend,$hmmstart,$hmmend]
	    }
	}
    }
    $/ = "\n";
    return %results;
}
sub parse_metagene {
    my $infile = shift @_;
    my $orfs;
    my $fixed = `grep -c \">>>\" $infile`;
    chomp($fixed);
    unless ($fixed > 0) {
	system("perl -pi -e \'s/\# gc/<gc/g\' $infile");
	system("perl -pi -e \'s/\# +self/<self/g\' $infile");
	system("perl -pi -e \'s/\# />>>\n/g\' $infile");
    }
    $/ = ">>>\n";
    open IN, "<$infile" or die $!;
 W1:while (my $record = <IN>) {
	chomp $record;
    W1:next if ($record =~ m/^\s*$/);
	my ($defline, @lines) = split(/\n/, $record);
	$defline =~ s/#\s*//g;
	my $id = (split(/\s+/, $defline))[0];
	
	## read % gc line
	my $temp = shift @lines;
	chomp $temp;
	$temp =~ /^<gc\s*=\s*(\S+)/ || die "Failed parsing gc line";
	my $gc = $1;
    F1:foreach (@lines) {
	    chomp;
	    next F1 if ($_ =~ m/\</ || $_ =~ m/^\s*$/);
	    my ($gid,$startpos,$endpos,$strand,
		$frame,$state,$score,$model,
		$rbsstart,$rbsstop,$rbsscore) = split("\t", $_);
	    ## set flags for partial ORFs
	    my ($five_prime_complete, $three_prime_complete) = split(//, $state);
	    
	    push (@{$orfs->{$id}},  {
				     'gc'          =>  $gc,
				     'model'       =>  $model,
				     'startpos'    =>  $startpos,
				     'offset'      =>  $startpos - 1,
				     'endpos'      =>  $endpos,
				     'complement'  =>  ($strand eq '+') ? 0 : 1,
				     'frame'       =>  $frame,
				     'score'       =>  $score,
				     'strand'      =>  $strand,
				     '5_partial'   =>  ($five_prime_complete == 1) ? 1 : 0,
				     '3_partial'   =>  ($three_prime_complete == 1) ? 1 : 0,
				     'id'          =>   join("|",$id,'metagene',$gid)
});
	    $orf_count++;
	}
    }
    $/ = "\n";
    return ($orfs, $orf_count);
}
sub run_orthomcl {
    my @seqfiles = @_;
    my $fref = join(",", @seqfiles);
    $command = "$global{orthotmp}/ORTHOMCLV1.4/orthomcl.pl --mode 1 --fa_files \"$fref\"";
    warn $command,"\n";
    my $output = `$command`;
    $output =~ m/($global{orthotmp}\/(\S+?)\/)/;
    return $1;
}
sub parse_orthomcl {
    my $seqfiledir = shift @_;
    my $orthomcldir = shift @_;
    my $opt_modo = shift @_;
    my $sfileref = shift @_;
    my %cname2file = %{shift @_};
    my %cname2orginfo = %{shift @_};
    my $orthomclout = $orthomcldir."all_orthomcl.out";
    my $allseqsfile = $orthomcldir."tmp/all.fa";
    my @files = @{$sfileref};
    my (%seq2file, %orthomcl, %geneannot);
    my $fref = join(",", @files);
    $command = "$global{orthotmp}/ORTHOMCLV1.4/orthomcl.pl --mode 1 --fa_files \"$fref\"";
    warn $command,"\n";
    open MCL, "<$orthomclout" or die "ORTHOMCL FAILED";
    while (my $line = <MCL>) {
	chomp($line);
	$i++;
	my $id = sprintf('%04s', $i);
	my $gname = "C$id";
	my ($name, $other)= split(/\t+/, $line);
	my @ary = split(/\s+/, $other);
	$name =~ m/\w+\((\d+)\s+genes,(\d+)\s+taxa\)/;
	$orthomcl{$gname}{num_genes} = $1;
	$orthomcl{$gname}{num_taxa} = $2;
	my %taxarep = ();
	my %cat;
	my %all;
	my $has_sign = 0;
	my $has_gpi = 0;
	foreach $g (@ary) {
	    $g =~ s/_\((\w+)\)/$1/;
	    my ($cname, $fam,$coor, $org, $file, $start, $end);
	    if ($opt_modo) {
		($cname,$fam,$coor,$org,$file) = split(/[\(\|]/, $g);
		($start, $end) = split(/-/, $coor);
		$uniname = join("\|", $cname, $fam, $coor);
	    }else {
		($cname, $org, $file) = split(/[\(\|]/, $g);
		$uniname = $cname;
	    }
	    if ($file) {
		$db = 'cazy_7';
		$cname =~ m/(\d+)/;
		$cnum = $1;
		$file = $cname2file{$cname};
		push @{$orthomcl{$gname}{members}}, $uniname;
		push @{$orthomcl{$gname}{files}{$file}}, $uniname;
		push @{$orthomcl{$gname}{dbs}{$db}}, $cnum;
		my $orgname = $cname2orginfo{$cname}{orgname};
		my $taxid = $cname2orginfo{$cname}{taxid};
		$seq2file{$uniname} = $file;
		$geneannot{$uniname}{db} = $db;
		$geneannot{$uniname}{czog} = $gname;
		$geneannot{$uniname}{org} = $orgname;
		$geneannot{$uniname}{taxid} = $taxid;
		if (not $opt_modo) {
		    my $get_modules = get_modules_by_entryid($cnum);
		    my @modules = ();
		    my @catmodules = ();
		    foreach my $row (sort {$a->[2] <=> $b->[2]} @{$get_modules}) {
			my ($fam, $sfam, $start, $end, $rstart, $rend) = @{$row};
			if ($fam eq 'SIGN') {
			    $has_sign ++;
			}elsif ($fam eq 'GPI') {
			    $has_gpi ++;
			}
			if (not $fam =~ m/SIGN|TM|UNK|LNK|GPI|SORT|PROP/) {
			    if ($sfam && $sfam ne '') {
				push @modules, join("_", $fam, $sfam);
			    }else {
				push @modules, $fam;
			    }
			}if ($fam =~ m/GH|GT|PL|CE/) {
			    $orthomcl{$gname}{fams}{$fam} = 1;
			    if ($sfam && $sfam ne '') {
				push @catmodules, join("_", $fam, $sfam);
			    }else {
				push @catmodules, $fam;
			    }
			}
		    }
		    $cat{$uniname} = join("\|", @catmodules);
		    $all{$uniname} = join("\|", @modules);
		}else {
		    $orthomcl{$gname}{fams}{(split(/_/,$fam))[0]} = 1;
		    $cat{$uniname} = $fam;
		    $all{$uniname} = $fam;
		}
		$orthomcl{$gname}{hassign} = $has_sign;
		$orthomcl{$gname}{hasgpi} = $has_gpi;
	    }
	}
	my %catstr;
	my %allstr;
	foreach $gene (keys %cat) {
	    $geneannot{$gene}{catmod} = $cat{$gene};
	    $geneannot{$gene}{allmod} = $all{$gene};
	    $catstr{$cat{$gene}} ++;
	    $allstr{$all{$gene}} ++;
	}
	$orthomcl{$gname}{numcatstr} = scalar(keys %catstr);
	$orthomcl{$gname}{numallstr} = scalar(keys %allstr);
	my @catstrdata;
	foreach (keys %catstr) {
	    push @catstrdata, "$_\($catstr{$_}\)"
	}
	$orthomcl{$gname}{posscatstr} = join("<br>", @catstrdata);
	my @allstrdata;
	foreach (keys %allstr) {
	    push @allstrdata, "$_\($allstr{$_}\)"
	}
	$orthomcl{$gname}{possallstr} = join("<br>", @allstrdata);
    }
    close MCL;
    my %singrp = ();
    open FA, "<$allseqsfile" or die $!;
    while (my $line = <FA>) {
	chomp($line);
	if ($line =~ m/>(\S+)/) {
	    if ($opt_modo) {
		($cname,$fam,$coor,$org) = split(/\|/, $1);
		($start, $end) = split(/-/, $coor);
		$uniname = join("\|", $cname, $fam, $coor);
	    }else {
		($cname, $org) = split(/\|/, $1);
		$uniname = $cname;
	    }
	    if (not $geneannot{$uniname}) {
		$db = 'cazy_7';
		$cname =~ m/(\d+)/;
		$cnum = $1;
		my $orgname = $cname2orginfo{$cname}{taxid};
		my $taxid = $cname2orginfo{$cname}{orgname};
		$geneannot{$uniname}{org} = $orgname;
		$geneannot{$uniname}{taxid} = $taxid;
		$file = $cname2file{$cname};
		if ($opt_modo) {
		    $geneannot{$uniname}{catmod} = $fam;
		    $geneannot{$uniname}{allmod} = $fam;
		    $gname = "S_$fam";
		    $orthomcl{$gname}{fams}{(split(/_/,$fam))[0]} = 1;
		    push @{$orthomcl{$gname}{files}{$file}}, $uniname;
		    push @{$orthomcl{$gname}{dbs}{$db}}, $cnum;
		    push @{$orthomcl{$gname}{members}}, $uniname;
		    $orthomcl{$gname}{catstr}{$fam} = 1;
		    $orthomcl{$gname}{allstr}{$fam} = 1;
		    $singrp{$gname} = 1;
		}else {
		    my $get_modules = get_modules_by_entryid($cnum);
		    my @modules = ();
		    my @catmodules = ();
		    my $has_sign = 0;
		    my $has_gpi = 0;
		    my %fams = ();
		    foreach my $row (sort {$a->[2] <=> $b->[2]} @{$get_modules}) {
			my ($fam, $sfam, $start, $end, $rstart, $rend) = @{$row};
			my $gname = "S_$fam";
			if ($fam eq 'SIGN') {
			    $has_sign ++;
			}elsif ($fam eq 'GPI') {
			    $has_gpi ++;
			}
			if (not $fam =~ m/SIGN|TM|UNK|LNK|GPI|SORT|PROP/) {
			    if ($sfam && $sfam ne '') {
				push @modules, join("_", $fam, $sfam);
			    }else {
				push @modules, $fam;
			    }
			}if ($fam =~ m/GH|GT|PL|CE/) {
			    $fams{$fam} = 1;
			    if ($sfam && $sfam ne '') {
				push @catmodules, join("_", $fam, $sfam);
			    }else {
				push @catmodules, $fam;
			    }
			}
		    }
		    $cat{$uniname} = join("\|", @catmodules);
		    $all{$uniname} = join("\|", @modules);
		    $geneannot{$uniname}{catmod} = join("\|", @catmodules);
		    $geneannot{$uniname}{allmod} = join("\|", @catmodules);
		    foreach $fam (keys %fams) {
			$gname = "S_$fam";
			$orthomcl{$gname}{fams}{$fam} = 1;
			$orthomcl{$gname}{hassign} += $has_sign;
			$orthomcl{$gname}{hasgpi} += $has_gpi;
			push @{$orthomcl{$gname}{files}{$file}}, $uniname;
			push @{$orthomcl{$gname}{dbs}{$db}}, $cnum;
			push @{$orthomcl{$gname}{members}}, $uniname;
			$orthomcl{$gname}{catstr}{$cat{$cname}} = 1;
			$orthomcl{$gname}{allstr}{$all{$cname}} = 1;
			$singrp{$gname} = 1;
		    }
		}
	    }
	}
    }
    foreach $singrp (keys %singrp) {
	$orthomcl{$singrp}{num_genes} = scalar(@{$orthomcl{$singrp}{members}});
	$orthomcl{$singrp}{num_taxa} = scalar(keys %{$orthomcl{$singrp}{files}});
	$orthomcl{$singrp}{numcatstr} = scalar(keys %{$orthomcl{$singrp}{catstr}});
	$orthomcl{$singrp}{numallstr} = scalar(keys %{$orthomcl{$singrp}{allstr}});
	$orthomcl{$singrp}{posscatstr} = join("<br>", keys %{$orthomcl{$singrp}{catstr}});
	$orthomcl{$singrp}{possallstr} = join("<br>", keys %{$orthomcl{$singrp}{allstr}});
    }
    my @lines = ();
    my %outhash = ();
    foreach $og(keys %orthomcl) {
	my @ostats = ();
	foreach $f (@files) {
	    my $file = (split(/[\.\/]/, $f))[-2];
	    if ($orthomcl{$og}{files}{$file}) {
		push @ostats, scalar(@{$orthomcl{$og}{files}{$file}});
		$outhash{$file}{$og} = scalar(@{$orthomcl{$og}{files}{$file}});
	    }else {
		push @ostats, 0;
		$outhash{$file}{$og} = 0;
	    }
	}
	$ids1 = join(",",@{$orthomcl{$og}{dbs}{cazy_7}});
	foreach $fg (keys %{$orthomcl{$og}{fams}}) {
	    $fg =~ m/(\D+)(\d+)/;
	    $fgp = $1;
	    $fnm = $2;
	    if (not ($fgp && $fnm)) {
		$fg =~ m/(\w+)\*/;
		$fgp = $1;
		$fnm = "9999";
	    }
	    push @lines, [$pgrps{$fgp},$fgp, $fnm, $og, $orthomcl{$og}{num_genes},$orthomcl{$og}{num_taxa}, $orthomcl{$og}{numcatstr}, $orthomcl{$og}{numallstr}, $orthomcl{$og}{hassign}, $orthomcl{$og}{hasgpi},$orthomcl{$og}{posscatstr},$orthomcl{$og}{possallstr}, $ids1,@ostats];
	}
    }
    @lines = sort {$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] || $a->[3] <=> $b->[3]} @lines;
    return(\@lines, \%outhash);
}

return 1;
