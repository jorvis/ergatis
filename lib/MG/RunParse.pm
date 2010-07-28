package MG::RunParse;

use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use Data::Dumper;
use XML::Simple;
use MG::ManFasta;
use MG::Common;
use Bio::SearchIO;
@ISA = qw(Exporter);

@EXPORT = qw(&run_search &parse_search &run_orthomcl &parse_orthomcl &runparse_lalign &runparse_ssearch &activesite_ssearch &runparse_metagene &remove_overlaps);

my %global = MG::Common::initialize();

sub parse_search {
    my $algo = shift @_;
    my $report = shift @_;
    my $eval_min = shift @_;
    my $dna2prot = shift @_;
    my %bhits = ();
    my %trans = ();
    my $sio = new Bio::SearchIO(-format => $algo,
				-file   => $report);
    
 W1:while (my $result = $sio->next_result) {
	my $query_length = $result->query_length;
	my $query = $result->query_name;
    W2:while (my $hit = $result->next_hit) {
	    my $hit_name = $hit->name; #could also be $hit->description
	    my $per_qcov = $hit->frac_aligned_query;
	    my $per_hcov = $hit->frac_aligned_hit;
	    my $hit_length = $hit->length;
	    my $hit_desc = $hit->description;
	W3:while (my $hsp = $hit->next_hsp) {
		my $eval = $hsp->significance();
		$eval =~ s/,//g;
		$eval = "1".$eval if $eval =~ m/^e/;
		my $bits = $hsp->bits;
		my $alen = $hsp->length("hit");
		my $percid = 100 * sprintf("%.2f", $hsp->frac_identical);
		my $percsim = sprintf("%.2f", $hsp->frac_conserved);
		my $qbegin = $hsp->start('query');
		my $qend = $hsp->end('query');
		my $lbegin = $hsp->start('hit');
		my $lend = $hsp->end('hit');	
		next W2 unless ($eval < $eval_min);
		if ($dna2prot) {
		    my $strand = "+";
		    if ($qbegin > $qend) {
			$strand = '-';
			my $t = $qbegin;
			$qbegin = $qend;
			$qend = $t;
		    }
		    my $query_seq =  $hsp->seq_str;
		    $query_seq =~ s/-|\\|\///g;
		    $query_seq =~ s/\*/X/g;
		    my $gap = $hsp->generate_cigar_string;
		    $trans{$query}{$strand}{$hit_name} = $query_seq;
		    push @{$bhits{$query}{$strand}}, [$hit_name, $eval, $percid, $strand, 
						      $qbegin, $qend, $query_length, $lbegin, 
						      $lend, $hit_length];
		}else {
		    push @{$bhits{$query}}, [$hit_name, $eval, $bits, $alen, $percid, $percsim, 
					     $per_qcov, $per_hcov, $query_length, $hit_length, 
					     $lbegin, $lend,$hit_desc,$qbegin, $qend];
		}
	    }
	}
    }
    if ($dna2prot) {
	return(\%bhits, \%trans);
    }else {
	return %bhits;
    }
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
				     '5_partial'   =>  ($five_prime_complete == 1) ? 1 : 0,
				     '3_partial'   =>  ($three_prime_complete == 1) ? 1 : 0,
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
