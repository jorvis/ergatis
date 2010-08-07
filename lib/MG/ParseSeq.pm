package MG::ParseSeq;

use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use FileHandle;
use MG::Math;
@ISA = qw(Exporter);
@EXPORT = qw(&parsegenpt &parsegb &parsefasta &parsesp &parsegb4dna);

use vars qw(%month);
%month = ("JAN"=>'01',"FEB"=>'02',"MAR"=>'03',
	  "APR"=>'04',"MAY"=>'05',"JUN"=>'06',
	  "JUL"=>'07',"AUG"=>'08',"SEP"=>'09',
	  "OCT"=>'10',"NOV"=>'11',"DEC"=>'12');


sub parsefasta {
    $/ = "\n>";
    my $infile = shift @_;
    my $source = shift @_;
    my $geno =   shift @_;
    open IN, "<$infile" or die "$infile  ??? ".$!;
    my %length;
    my %results;
    my $i = 0;
 W1:while (my $record = <IN>) {
	chomp($record);
	my @rec = split(/\n/, $record);
	my $header = shift @rec;
	my $seq = join("",@rec);
	next W1 unless $seq;
	$i ++;
	my ($seqname,@desc) = split(/\s+/, $header);
	$desc = join(" ", @desc);
	$header =~ s/^[\s*|>]//g;
	$seqname =~ s/^[\s*|>]//g;
	$seqname =~ s/\s*$//g;
	$seqname =~ s/:|-/_/g;
	$seqname =~ s/,\s*//g;
	$desc = (split(/\x1/, $desc))[0] if $desc;
	$desc =~ s/,\s*/ /g if $desc;
	$results{$i}{ori_head} = $header;
	$results{$i}{db_desc} = $desc;
	$results{$i}{is_confidential} = 1;
	$results{$i}{is_auto} = 1;
	$results{$i}{db_name} = "other";
	if ($seqname =~ m/^gi\|/) {
	    $results{$i}{gi} = (split(/\|/, $seqname))[1];
	    $results{$i}{db_acc} = (split(/\|/, $seqname))[3]; 
	    my @ary = split(/\[/, $desc);
	    $results{$i}{species} = pop @ary;
	    $results{$i}{species} =~ s/\]//g;
	    $results{$i}{db_desc} = join(" ", @ary);
	    $results{$i}{db_name} = "ncbi";
	    $results{$i}{description} = $results{$i}{db_acc};
	    $results{$i}{is_confidential} = 0;
	    $results{$i}{db_name} = "ncbi";
	}elsif ($seqname =~ m/^jgi\|(.+)/) {
	    ($project,$gpid,$pname) = split(/\|/, $1);
	    $results{$i}{db_acc} = $gpid;
	    $results{$i}{db_desc} = $header;
	    $results{$i}{locus_tag} = $pname;
	    $results{$i}{db_name} = "jgi"
	}elsif ($seqname =~ m/^gnl\|(.+)/) {
	    $results{$i}{db_desc} = $1;
	    $results{$i}{db_acc} = (split(/\s+/, $1))[0];
	    $results{$i}{db_acc} =~ s/\|/_/g;
	}else {
	    $results{$i}{db_desc} = $seqname;
	    $results{$i}{db_acc} = (split(/\s+/, $seqname))[0]; 
	}
	$results{$i}{description}= $results{$i}{db_acc} if not $results{$i}{description};
	$results{$i}{annot_note} = $source;
	$results{$i}{organism} = $results{$i}{species};
	$results{$i}{organism} .= " ".$results{$i}{strain} if ($results{$i}{strain});
	$seq =~ s/%([\da-f]{1,2})/pack(C,hex($1))/eig;
	$seq =~ s/[\+]/ /g;
	$seq =~ s/\r//g;
	$seq =~ s/\*//g;
	$seq =~ s/J|j|U|u/X/g;
	$results{$i}{sequence} = uc($seq);
	$results{$i}{length} = length($seq);
	$length{$seqname} = length($seq);
    }
    $/ = "\n";
    return(\%results,\%length);
}
sub parsegenpt {
    $/ = "//\n";
    my $infile = shift @_;
    my $source = shift @_;
    my $genoid = shift @_;
    my $geno =   shift @_;
    open IN, "<$infile" or die $!;
    my %orgnames;
    my %results;
    my $i = 0;
    while (my $record = <IN>) {
	chomp($record);
	$i ++;
	my @pubmed;
	next if ($record =~ m/^\s$/);
	my ($data, $seq) = split(/ORIGIN\s*\n/, $record);
	$data =~ s/\*/\\\*/g;
	$seq =~ s/[\s0-9\/\n]//g;
	$seq =~ s/%([\da-f]{1,2})/pack(C,hex($1))/eig;
	$seq =~ s/[\+]/ /g;
	$seq =~ s/\r//g;
	$seq =~ s/\*//g;
	$seq =~ s/J|j/X/g;
	$results{$i}{length} = length($seq);
	$results{$i}{sequence} = uc($seq);
	$results{$i}{is_confidential} = 0;
	$results{$i}{db_name} = 'ncbi';
	$results{$i}{is_auto} = 1;
	$results{$i}{genome_id} = $genoid if $genoid;
	while ($data =~ /^[A-Z].*\n(^\s.*\n)*/gm) {
            my $value = $&;
            (my $key = $value) =~ s/^([A-Z]+).*/$1/s;
            $value =~ s/\s*\n\s*/ /g;
            chomp($value);
            $value =~ s/^[A-Z]+//;
            $value =~ s/\n\s+/ /g;
            $value =~ s/\s+\//\n/g;
            $value =~ s/^\s*(.+)\s*$/$1/;
            if ($key eq 'REFERENCE') {
                while ($value =~ m/PUBMED\s+(\d+)/g) {
                    push @pubmed, $1;
                }
            }elsif ($key eq 'VERSION') {
                my($gpid, $gi) = split(/\s+/, $value);
                $gi =~ s/GI://;
                my ($acc, $version) = split(/\./, $gpid);
                %{$results{$i}} = (%{$results{$i}},gi=>$gi,db_acc=>$gpid);
	    }elsif ($key eq 'LOCUS') {
		my $date = (split(/\s+/, $value))[-1];
		my ($day, $month, $year) = split(/-/, $date);
		my $nmonth = $month{$month};
		$results{$i}{db_date} = join("-", $year, $nmonth, $day);
	    } elsif ($key eq 'DBSOURCE') {
		$value =~ m/accession\s+(\S+)/;
		$results{$i}{extra_acc} = $1;
	    } elsif ($key eq 'FEATURES') {
		if ($value =~ m/organism=\"(.+)\"/) {
		    $results{$i}{species} = $1;
		}
		if ($value =~ m/strain=\"(.+)\"/) {
		    $results{$i}{strain} = $1;
		}
		if ($value =~ m/db_xref=\"taxon:(\d+)\"/) {$results{$i}{tax_id} = $1;}
		if ($value =~ m/gene=\"(\S+)\"/) {
		    $results{$i}{gene} = lc($1);
		    $results{$i}{gene} =~ s/^(\w)/uc($1)/eg;
		    $results{$i}{gene} =~ s/(i+)$/uc($1)/eg;
		    $results{$i}{gene} =~ s/(\w)$/uc($1)/eg;
		}
		if ($value =~ m/chromosome=\"(\S+)\"/) {
		    $results{$i}{chromosome} = $1;
		}
		while ($value =~ m/EC_number=\"([\w\.-]+)\"/g) {push @{$results{$i}{ec}}, $1;}
		if ($value =~ m/product=\"(.+)\"/) {$results{$i}{product} = $1;}
		if ($value =~ m/coded_by=\"(.+)\"/) {$results{$i}{coded_by} = $1;} 
		if ($value =~ m/locus_tag=\"(.+)\"/) {$results{$i}{locus_tag} = $1;} 
		if ($value =~ m/\/name=\"(.+)\"/) {
		    $results{$i}{pname} = $1;
		    if ($results{$i}{pname} =~ m/\w\w(\w{4}?)/) {
			$results{$i}{pname} = $1;
		    }elsif ($results{$i}{pname} =~ m/\w*\w*(\w{3}?\d{1,3}?\w)/) {
			$results{$i}{pname} = $1;
		    }
		}
		if ($results{$i}{coded_by} && $results{$i}{coded_by} =~ m/\<.+\>/) {
		    $results{$i}{entry_note} = "Fragment N-Term/C-Term"
		}elsif ($results{$i}{coded_by} && $results{$i}{coded_by} =~ m/\</) {
		    if ($results{$i}{coded_by} =~ m/complement/) {
			$results{$i}{entry_note} .= "Fragment C-Term";
		    }else {
			$results{$i}{entry_note} .= "Fragment N-Term";
		    }
		}elsif ($results{$i}{coded_by} && $results{$i}{coded_by} =~ m/\>/) {
		    if ($results{$i}{coded_by} =~ m/complement/) {
			$results{$i}{entry_note} .= "Fragment N-Term";
		    }else {
			$results{$i}{entry_note} .= "Fragment C-Term";
		    }
		}
	    } else {
		$results{$i}{lc($key)} = $value;
	    }
	}
	if (not $results{$i}{gene}) {
	    $results{$i}{gene} = $results{$i}{pname};
	    $results{$i}{gene} = lc($1);
	    $results{$i}{gene} =~ s/^(\w)/uc($1)/eg;
	    $results{$i}{gene} =~ s/(i+)$/uc($1)/eg;
	    $results{$i}{gene} =~ s/(\w)$/uc($1)/eg;
	}
	$results{$i}{db_desc} = (split(/\s*\[/, $results{$i}{definition}))[0];
	if ($results{$i}{product} && $results{$i}{product} =~ m/neuraminidase/) {
	    $results{$i}{description} = "neuraminidase";
	}elsif ($results{$i}{locus_tag}) {
	    $results{$i}{description} = $results{$i}{locus_tag};
	}else {
	    $results{$i}{description} = $results{$i}{db_desc};
	}if ($results{$i}{gene}) {
	    $results{$i}{description} .= " \($results{$i}{gene}\)";
	}
	$results{$i}{organism} = $results{$i}{species};

	if ($results{$i}{strain} && $results{$i}{species} !~ m/\Q$results{$i}{strain}\E/) {
	    $results{$i}{org_xtra} .= $results{$i}{strain};
	    $results{$i}{organism} .= " ".$results{$i}{strain};
	}

	$results{$i}{pubmed} = \@pubmed;
	if ($results{$i}{db_desc} =~ m/(putative|potential|potencial|homolog|like|hypothetical|similar\s+to)/) {
	    delete $results{$i}{ec};
	    if (not $results{$i}{locus_tag}) {
		$results{$i}{description} = 'ORF';
	    }
	}
	delete $results{$i}{source};
	$orgnames{$results{$i}{organism}} ++;
    }
    return (\%results, \%orgnames);
}
sub parsesp {
    $/ = "//\n";
    my $infile = shift @_;
    my $source = shift @_;
    my $genoid = shift @_;
    my $geno =   shift @_;
    open IN, "<$infile" or die $!;
    my %orgnames;
    my %results;
    my $i = 0;
    while (my $record = <IN>) {
	chomp($record);
	$i ++;
	my @pubmed;
	my ($data, $seq) = split(/\nSQ\s.+\n/, $record);
	$data =~ s/\*/\\\*/g;
	$data =~ s/[\(\)]//g;
	$seq =~ s/[\s0-9\/\n]//g;
	$seq =~ s/%([\da-f]{1,2})/pack(C,hex($1))/eig;
	$seq =~ s/[\+]/ /g;
	$seq =~ s/\r//g;
	$seq =~ s/\*//g;
	$seq =~ s/J|j/X/g;
	$results{$i}{length} = length($seq);
	$results{$i}{sequence} = uc($seq);
	$results{$i}{is_confidential} = 0;
	$results{$i}{db_name} = 'uniprot';
	$results{$i}{is_auto} = 1;
	$results{$i}{genome_id} = $genoid if $genoid;
	my %data = ();
	foreach my $line (split(/\n/, $data)) {
	    $line =~ m/^(\w{2}?)\s+(\S.+)/;
	    $data{$1} .= " ".$2;
	}
	if ($data{'RX'}) {
	    while ($data{'RX'} =~ m/PubMed=(\d+)/g) {
		push @pubmed, $1;
	    }
	}if ($data{ID} =~ m/(\w+)/) {
	    $results{$i}{locus_tag} = $1;
	}if ($data{AC} =~ m/(\w+)/) {
	    $results{$i}{db_acc} = $1;
	}while ($data{DT} =~ m/(\d+-\w+-\d+)/g) {
	    my ($day, $month, $year) = split(/-/, $1);
	    my $nmonth = $month{$month};
	    $results{$i}{db_date} = join("-", $year, $nmonth, $day);
	}if ($data{DE}) {
	    $results{$i}{definition} = $data{DE};
	    $results{$i}{definition} =~ s/^\s*//;
	    $results{$i}{db_desc} = $results{$i}{definition};
	}if ($data{GN} =~ m/Name=(\w+);/) {
	    $results{$i}{gene} = lc($1);
	    $results{$i}{gene} =~ s/^(\w)/uc($1)/eg;
	    $results{$i}{gene} =~ s/(i+)$/uc($1)/eg;
	    $results{$i}{gene} =~ s/(\w)$/uc($1)/eg;
	}if ($data{OS}) {

	    $results{$i}{species} = $data{OS};
	    $results{$i}{species} =~ s/^\s*//;
	}if ($data{DE} =~ m/EC ([\d\.]+\d)/) {
	    push @{$results{$i}{ec}}, $1;
	}if ($data{DR} =~ m/EMBL; (\w+);/) {
	    $results{$i}{extra_acc} = $1;
	}if ($data{OX} =~ m/NCBI_TaxID=(\d+)/) {
	    $results{$i}{tax_id} = $1;
	}
	$results{$i}{pubmed} = \@pubmed;
	$results{$i}{organism} = $results{$i}{species};
	$orgnames{$results{$i}{organism}} ++;
    }
    $/="\n";
    return (\%results, \%orgnames);
}
sub parsegb4dna {
    $/ = "//\n";
    my $infile = shift @_;
    my $source = shift @_;
    my $genoid = shift @_;
    my $geno =   shift @_;
    open IN, "<$infile" or die $!;
    my %gbdata = ();
    my $i = 0;
 W1:while (my $record = <IN>) {
	chomp($record);
	my ($data, $seq) = split(/ORIGIN\s*\n/, $record);
	next W1 unless $seq;
	$data =~ s/\*/\\\*/g;
	$seq =~ s/[\s0-9\/\n]//g;
	$seq =~ s/%([\da-f]{1,2})/pack(C,hex($1))/eig;
	$seq =~ s/[\+]/ /g;
	$seq =~ s/\r//g;
	$seq =~ s/\*//g;
	$seq =~ s/J|j/X/g;
	$i ++;
	$gbdata{$i}{sequence} = $seq;
	my @section = split(/\n([A-Z]+)/, $data);
	unshift @section, "LOCUS";
    F2:foreach (1..(scalar(@section)/2)) {
	    $key = shift @section;
	    $value = shift @section;
	    $value =~ s/\s*\n\s*/\n/g;
	    $value =~ s/^\s*(.+)\s*$/$1/;
	    chomp($value);
	    if ($key eq 'REFERENCE') {
		while ($value =~ m/PUBMED\s+(\d+)/g) {
		    push @pubmed, $1;
		}
	    }elsif ($key eq 'VERSION') {
		my($gpid, $gi) = split(/\s+/, $value);
		$gbdata{$i}{db_acc} = $gpid;
		$gbdata{$i}{gi} = $gi;
		$gbdata{$i}{gi} =~ s/GI://g;
	    }elsif ($key eq 'DBLINK') {
		if ($value =~ m/Project:(\d+)/) {
		    $gbdata{$i}{project_id} = $1;
		}
	    }elsif ($key eq 'SOURCE') {
		my @source = split(/\n/, $value);
		shift @source;
		shift @source;
		$gbdata{$i}{taxtree} = join("", @source);
		$gbdata{$i}{taxtree} =~ s/\s+/ /g;
	    } elsif ($key eq 'LOCUS') {
		my $date = (split(/\s+/, $value))[-1];
		next F2 unless($date);
		my ($day, $month, $year) = split(/-/, $date);
		my $nmonth = $month{$month};
		$gbdata{$i}{db_date} = join("-", $year, $nmonth, $day);
	    } elsif ($key eq 'FEATURES') {
		$value =~ s/gene\s+/\/gene\n\/gene_coor=/g;
		$value =~ s/mRNA\s+/\/mRNA\n\/mrna_coor=/g;
		$value =~ s/CDS\s+/\/CDS\n\/cds_coor=/g;
		$value =~ s/\n\s*\//\n>>><<</g;
		$value =~ s/\s*\n\s*/ /g;
		$value =~ s/\s+>>><<</\n/g;
		$value =~ s/\n\s+/ /g;
		if ($value =~ m/organism=\"(.+)\"/) {
		    $gbdata{$i}{species} = $1;
		}if ($value =~ m/strain=\"(.+)\"/) {
		    $gbdata{$i}{strain} = $1;
		}if ($value =~ m/mol_type=\"(.+)\"/) {
		    $gbdata{$i}{mol_type} = $1;
		}if ($value =~ m/db_xref=\"taxon:(\d+)\"/) {
		    $gbdata{$i}{tax_id} = $1;
		}
		$gbdata{$i}{pubmed} = \@pubmed;
		$gbdata{$i}{organism} = $gbdata{$i}{species};
		$gbdata{$i}{organism} .= " ".$gbdata{$i}{strain} if ($gbdata{$i}{species} !~ m/$gbdata{$i}{strain}/);
	    }
	}
    }
    $/ = "\n";
    return %gbdata;
}
sub parsegb {
    $/ = "//\n";
    my $infile = shift @_;
    open IN, "<$infile" or die $!;
    my %results;
    my $i = 0;
    my $g = 0;
 W1:while (my $record = <IN>) {
	chomp($record);
	my %gbdata = ();
	my @pubmed;
	my ($data, $scafseq) = split(/ORIGIN\s*\n/, $record);
	next W1 unless $data =~ m/\s+CDS\s+/;
	$scafseq =~ s/[\s0-9\/\n]//g;
	$scafseq =~ s/%([\da-f]{1,2})/pack(C,hex($1))/eig;
	$scafseq =~ s/[\+]/ /g;
	$scafseq =~ s/\r//g;
	$scafseq =~ s/\*//g;
	$scafseq =~ s/J|j/X/g;
	$gbdata{scafseq} = uc($scafseq);
	$data =~ s/\*/\\\*/g;
	my @section = split(/\n([A-Z]+)/, $data);
	unshift @section, "LOCUS";
	foreach (1..(scalar(@section)/2)) {
	    $key = shift @section;
	    $value = shift @section;
	    $value =~ s/\s*\n\s*/\n/g;
	    $value =~ s/^\s*(.+)\s*$/$1/;
	    chomp($value);
	    if ($key eq 'REFERENCE') {
		while ($value =~ m/PUBMED\s+(\d+)/g) {
		    push @pubmed, $1;
		}
	    }elsif ($key eq 'VERSION') {
		my($gpid, $gi) = split(/\s+/, $value);
		$gbdata{extra_acc} = $gpid;
	    }elsif ($key eq 'LOCUS') {
		my @ary = split(/\s+/, $value);
		my ($day, $month, $year) = split(/-/, $ary[$#ary]);
		my $nmonth = $month{$month};
		$gbdata{db_date} = join("-", $year, $nmonth, $day) if ($year && $nmonth && $day);
	    } elsif ($key eq 'FEATURES') {
		$value =~ s/gene\s+/\/gene\n\/gene_coor=/g;
		$value =~ s/mRNA\s+/\/mRNA\n\/mrna_coor=/g;
		$value =~ s/CDS\s+/\/CDS\n\/cds_coor=/g;
		$value =~ s/\n\s*\//\n>>><<</g;
		$value =~ s/\s*\n\s*/ /g;
		$value =~ s/\s+>>><<</\n/g;
		$value =~ s/\n\s+/ /g;
		if ($value =~ m/organism=\"(.+)\"/) {
		    $gbdata{species} = $1;
		}if ($value =~ m/strain=\"(.+)\"/) {
		    $gbdata{strain} = $1;
		}if ($value =~ m/mol_type=\"(.+)\"/) {
		    $gbdata{mol_type} = $1;
		}if ($value =~ m/db_xref=\"taxon:(\d+)\"/) {
		    $gbdata{tax_id} = $1;
		}
		$gbdata{pubmed} = \@pubmed;
		$gbdata{organism} = $gbdata{species} if ($gbdata{species});
		$gbdata{organism} .= " ".$gbdata{strain} if ($gbdata{species} && $gbdata{species} !~ m/$gbdata{strain}/);
		my @genes = split(/\nCDS\s*\n/, $value);
	    FCDS:foreach $cds (@genes[1..$#genes]) {
		    next FCDS unless $cds =~ m/translation/;
		    $i ++;
		    foreach (keys %gbdata) {
			$results{$i}{$_} = $gbdata{$_};
		    }
		    print $i,"\n";
		    $results{$i}{genome_id} = $genoid if $genoid;
		    if ($cds =~ m/db_xref=\"GI:(\d+)\"/) {$results{$i}{gi} = $1;}
		    if ($cds =~ m/protein_id=\"(\S+)\"/) {
			$results{$i}{db_acc} = $1;
		    }if ($cds =~ m/translation=\"(.+)\"/) {
			$results{$i}{sequence} = uc($1);
			$results{$i}{sequence} =~ s/ //g;
		    }
		    if ($cds =~ m/gene=\"(\S+)\"/) {$results{$i}{gene} = $1}
		    if ($cds =~ m/chromosome=\"(\S+)\"/) {
			$results{$i}{chromosome} = $1;
		    }
		    while ($cds =~ m/EC_number=\"([\w\.-]+)\"/g) {push @{$results{$i}{ec}}, $1;}
		    if ($cds =~ m/product=\"(.+)\"/) {
			$results{$i}{product} = $1;
		    }
		    if ($cds =~ /note=\"frame=(\d+)/) {
			$results{$i}{frame} = $1;
		    }if ($cds =~ m/score=\"(.+)\"/) {
			$results{$i}{score} = $1;
		    }
		    if ($cds =~ m/cds_coor=(.+)/) {
			$results{$i}{coded_by} = $1;
			$results{$i}{strand} = '+';
			$results{$i}{strand} = '-' if ($1 =~ m/complement/);
			$results{$i}{coded_by} =~ m/(\d+)[\.>]+(\d+)/;
			$results{$i}{genome_start} = $1;
			$results{$i}{genome_end} = $2;
			$b = $results{$i}{genome_start} - 1;
			$l = $results{$i}{genome_end} - $b;
			$results{$i}{dnaseq} = substr($gbdata{scafseq},$b,$l)
		    }
		    if ($cds =~ m/locus_tag=\"(.+)\"/) {$results{$i}{locus_tag} = $1;} 
		    if ($cds =~ m/\/note=\"(.+)\"/ && not $results{$i}{gene}) {
			$results{$i}{xtradesc} = $1;
			if ($results{$i}{xtradesc} =~ m/\w\w(\w{4}?)/) {
			    $results{$i}{gene} = $1 unless $results{$i}{gene};
			}elsif ($results{$i}{xtradesc} =~ m/\w*\w*(\w{3}?\d{1,3}?\w)/) {
			    $results{$i}{gene} = $1 unless $results{$i}{gene};
			}
		    }
		    if ($results{$i}{locus_tag}){
			$results{$i}{description} = $results{$i}{locus_tag};
		    }elsif ($results{$i}{desc} =~ m/(putative|potential|potencial|homolog|like|hypothetical|similar\s+to)/) {
			$results{$i}{description} = 'ORF';
			delete $results{$i}{ec} if $results{$i}{ec};
		    } elsif ($results{$i}{product}) {
			$results{$i}{description} = $results{$i}{product};
		    }else {
			$results{$i}{description} = 'ORF';
		    }
		    if ($results{$i}{gene} && $results{$i}{gene} =~ m/\w+/) {
			$results{$i}{gene} = lc($results{$i}{gene});
			$results{$i}{gene} =~ s/^(\w)/uc($1)/eg;
			$results{$i}{gene} =~ s/(i+)$/uc($1)/eg;
			$results{$i}{gene} =~ s/(\w)$/uc($1)/eg;
			$results{$i}{description} .= " \($results{$i}{gene}\)";
		    }
		    if ($results{$i}{coded_by} =~ m/\<.+\>/) {
			$results{$i}{entry_note} = "Fragment N-Term/C-Term"
		    }elsif ($results{$i}{coded_by} =~ m/\</) {
			if ($results{$i}{coded_by} =~ m/complement/) {
			    $results{$i}{entry_note} .= "Fragment C-Term";
			}else {
			    $results{$i}{entry_note} .= "Fragment N-Term";
			}
		    }elsif ($results{$i}{coded_by} =~ m/\>/) {
			if ($results{$i}{coded_by} =~ m/complement/) {
			    $results{$i}{entry_note} .= "Fragment N-Term";
			}else {
			    $results{$i}{entry_note} .= "Fragment C-Term";
			}
		    }
		    if ($results{$i}{entry_note}) {
			$results{$i}{db_desc} = join(" ",$results{$i}{description},$results{$i}{entry_note});
		    }
		    $results{$i}{length} = length($results{$i}{sequence});
		    $results{$i}{is_auto} = 1;
		    $results{$i}{db_desc} = $results{$i}{description};
		}
	    }
	}
    }
    $/ = "\n";
    return (\%results);
}
return 1;
