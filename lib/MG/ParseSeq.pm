package Annot::ParseSeq;

use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use FileHandle;
use MG::Math;
@ISA = qw(Exporter);
@EXPORT = qw(&parsegenpt &parsegb &parsefasta &parsesp &clean_header &delpat &parsegb4dna);

use vars qw(%month);
%month = ("JAN"=>'01',"FEB"=>'02',"MAR"=>'03',
	  "APR"=>'04',"MAY"=>'05',"JUN"=>'06',
	  "JUL"=>'07',"AUG"=>'08',"SEP"=>'09',
	  "OCT"=>'10',"NOV"=>'11',"DEC"=>'12');


sub parsefasta {
    $/ = "\n>";
    my $infile = shift @_;
    my $source = shift @_;
    my $genoid = shift @_;
    my $geno =   shift @_;
    open IN, "<$infile" or die "$infile  ??? ".$!;
    my %unique_gpid;
    my %orgnames;
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
	    $results{$i}{db_acc} =~ s/\||\./_/g;
	}
	$results{$i}{description}= $results{$i}{db_acc} if not $results{$i}{description};
	if ($geno) {
	    my ($taxid, $species, $strain) = split(/:/, $geno);
	    $results{$i}{tax_id} = $taxid;
	    $results{$i}{species} = $species;
	    $results{$i}{org_xtra} = $strain;
	    $results{$i}{mol_type} = "genomic DNA";
	}
	die "Locus Tag not Unique\n$seqname\n$results{$i}{db_acc}\n" if ($unique_gpid{$results{$i}{db_acc}});
	$unique_gpid{$results{$i}{db_acc}} = 1;
	$results{$i}{annot_note} = $source;
	$results{$i}{genome_id} = $genoid if $genoid;
	$results{$i}{length} = length($seq);
	$results{$i}{organism} = $results{$i}{species};
	$results{$i}{organism} .= " ".$results{$i}{strain} if ($results{$i}{strain});
	$seq =~ s/%([\da-f]{1,2})/pack(C,hex($1))/eig;
	$seq =~ s/[\+]/ /g;
	$seq =~ s/\r//g;
	$seq =~ s/\*//g;
	$seq =~ s/J|j|U|u/X/g;
	$results{$i}{sequence} = uc($seq);
	$orgnames{$results{$i}{organism}} ++ if $results{$i}{organism};
    }
    if ($geno) {
	my ($taxid, $species, $strain) = split(/:/, $geno);
	my $date = currdate();
	$date =~ s/-//g;
	%results = clean_header(\%results, $geno);
	my $newname = "cg--$date--$species--$strain--$taxid--$source.fa";
	$newname =~ s/\s+/_/g;
	system("cp $infile /tera/data/cazydata/complete_genomes/genomes_rename_file/$newname");
	system("mv $infile /tera/data/cazydata/complete_genomes/genomes_original_file/");
    }
    $/ = "\n";
    return(\%results,\%orgnames);
}

sub parsegenpt {
    $/ = "//\n";
    my $infile = shift @_;
    my $source = shift @_;
    my $genoid = shift @_;
    my $geno =   shift @_;
    open IN, "<$infile" or die $!;
    my %unique_gpid;
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
	$data =~ s/alpha-/a-/g;
	$data =~ s/beta-/b-/g;
	$data =~ s/gamma-/g-/g;
	$data =~ s/kappa-/k-/g;
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
		    $results{$i}{chromosome} =~ s/IX/9/g;
		    $results{$i}{chromosome} =~ s/VIII/8/g;
		    $results{$i}{chromosome} =~ s/VII/7/g;
		    $results{$i}{chromosome} =~ s/VI/6/g;
		    $results{$i}{chromosome} =~ s/V/5/g;
		    $results{$i}{chromosome} =~ s/IV/4/g;
		    $results{$i}{chromosome} =~ s/III/3/g;
		    $results{$i}{chromosome} =~ s/II/2/g;
		    $results{$i}{chromosome} =~ s/I/1/g;
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
    my %unique_gpid;
    my %orgnames;
    my %results;
    my $i = 0;
    while (my $record = <IN>) {
	chomp($record);
	$i ++;
	my @pubmed;
	my ($data, $seq) = split(/\nSQ\s.+\n/, $record);
	$data =~ s/\*/\\\*/g;
	$data =~ s/alpha-/a-/g;
	$data =~ s/beta-/b-/g;
	$data =~ s/gamma-/g-/g;
	$data =~ s/kappa-/k-/g;
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
	$data =~ s/alpha-/a-/g;
	$data =~ s/beta-/b-/g;
	$data =~ s/gamma-/g-/g;
	$data =~ s/kappa-/k-/g;
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

sub clean_header {
    my %seqhash = %{shift @_};
    my $geno = shift @_;
    my ($taxid, $species, $strain) = split(/:/, $geno);
    my $date = currdate();
    my $date2 = join("", split(/-/, $date));
    
    #---------------
    # Declarations
    # --------------
    my @headtmp=();
    my @testlet=();
    my $val= "";
    my $next="no";
    my $war;
    
    
 #---------------------------------
 # Traitement du header gi or jgi
 # --------------------------------

    $orf=scalar(keys %seqhash);
    $x = $seqhash{1}{ori_head};
    $length_x=length($x); 
    my @head = ();
    if($x=~m/^jgi/){# 2 patterns  JGI >jgi|blabla|Locus_id >jgi|blabla|Locus_id|Prot_name
	my @jgi=split(/\|/,$x); 
	foreach my $i (keys %seqhash) {
	    if($#jgi == 2) { if ($seqhash{$i}{ori_head}=~/^jgi\|.*\|(\d+)$/) { $seqhash{$i}{db_acc}=$1; } } 
	    elsif($#jgi == 3) { if ($seqhash{$i}{ori_head}=~/^jgi\|.*\|(\d+)\|(.*)$/) { $seqhash{$i}{db_acc}=$1; $seqhash{$i}{locus_tag}=$2; } }     
	} 
    }elsif($x=~m/^gi/){ # 1 pattern GI		>gi|Locus_id|ref|GP_ac|_Model_id
	foreach my $i (keys %seqhash) {
	    my @gi=split(/\|/,$seqhash{$i}{ori_head}); 
	    if ($gi[4]=~m/^_(.*)$/) {$seqhash{$i}{db_desc}=$1;}
	    elsif ($gi[4]=~m/^\s(.*)\s\[.*\]\s?$/)  {$seqhash{$i}{db_desc}=$1;}
	    $seqhash{$i}{db_acc}=$gi[3];
	} 
    } else { 
	#-------- pretraitement :  it removes all that one does not -------------- 
	foreach my $i (sort {$a <=> $b} keys %seqhash) {
	    $seqhash{$i}{ori_head}=~s/\|//g; 
	    $seqhash{$i}{ori_head}=~s/\(\d+\saa\)//;
	    $seqhash{$i}{ori_head}=~s/\d+\saa//;
	    $seqhash{$i}{ori_head}=~s/MW:\d+//;
	    $seqhash{$i}{ori_head}=~s/\(translation\)//;
	    $seqhash{$i}{ori_head}=~s/forward//;
	    $seqhash{$i}{ori_head}=~s/direct//;
	    $seqhash{$i}{ori_head}=~s/reverse//;
	    $seqhash{$i}{ori_head}=~s/org:[^ ]+//i;
	    $seqhash{$i}{ori_head}=~s/strain\:[^ ]+//i; 
	    $seqhash{$i}{ori_head}=~s/source\:[^ ]+//i;
	    $seqhash{$i}{ori_head}=~s/tax:\d+//i; 
	    $seqhash{$i}{ori_head}=~s/\[\s?$species\s?\]//i; 
	    $seqhash{$i}{ori_head}=~s/$species//i;
	    my @species=split(/\s/,$species);
	    for ( my $j=0; $j<=$#species; $j++) { $seqhash{$i}{ori_head}=~s/$species[$j]//i; }
	    $seqhash{$i}{ori_head}=~s/sp\.//i;  
	    if ($strain && $strain !~ m/^\s+$/) { $seqhash{$i}{ori_head}=~s/$strain//; }
	    $seqhash{$i}{ori_head}=~s/[^\s]*:?complement\(.*\)$//;
	    $seqhash{$i}{ori_head}=~s/[^\s]*:?join\(.*\)$//; 
	    $seqhash{$i}{ori_head}=~s/[^\s]*:?\d+\.\.\d+.*$//;
	    $seqhash{$i}{ori_head}=~s/\(//g; $seqhash{$i}{ori_head}=~s/\)//g;
	    $seqhash{$i}{ori_head}=~s/\{//g; $seqhash{$i}{ori_head}=~s/\}//g;
	    $seqhash{$i}{ori_head}=~s/\[//g; $seqhash{$i}{ori_head}=~s/\]//g;
	    $seqhash{$i}{ori_head}=~s/\s\s*?\d+\s\s*-\s\s*?\d+\s*?$//;
	    $seqhash{$i}{ori_head}=~s/\s\s*?\d+:\d+\s*?$//;
	    $seqhash{$i}{ori_head}=~s/,$//;
	    $seqhash{$i}{ori_head}=~s/_+\s?$//;
	    $seqhash{$i}{ori_head}=~s/\s+/ /g;
	    $seqhash{$i}{ori_head}=~s/\s$//;
	    $seqhash{$i}{ori_head}=~s/^\s//;         
	    push @head, $seqhash{$i}{ori_head}; 
	}
	my ($x,@headtmp) = @head;
	$length_x=length($x); 
	#--------fin du  pretraitement -----------------------------------------------
	
	######### IF pas d'espace / without space
	if($x !~ m/\s/) {
	    if($x=~m/contig/i){
		foreach my $i (keys %seqhash) {
		    if($seqhash{$i}{ori_head}=~m/^.*(contig.*)$/i) { $seqhash{$i}{db_acc}=$1;} 
		}
	    }
	    elsif($x=~m/:/) { 
		@twopts=split(/:/,$x);
		@locus=delpat($orf,scalar(@twopts),@twopts,@head);
	    }elsif($x=~m/_/){ 
		@underscore=split(/_/,$x);
		@locus=delpat($orf,scalar(@underscore),@underscore,@head); 
	    }else {
		for( my $j=0; $j<$length_x; $j++ ){ 
		    $ct=1;
		    $testlet[$j]=substr($x,0,$j+1);
		    for( my $i=0; $i<=$#headtmp; $i++ ){
			if($testlet[$j]=~m/.*\d/) {
			    if ($x=~m/.*\.\w+\..*/){ goto NEXT; } # cas vraiment particulier >bma1.CDS.218966.0 voir si ne fait pas buger le reste
			    $val=$testlet[$j]; $next="yes";
			    goto SORT; 
			}
		    NEXT:if ($headtmp[$i]=~m/^$testlet[$j]/) {$ct++;}
		    } 
		    if($orf != $ct) {$val=$testlet[$j]; $next="yes"; goto SORT; }
		}
	    }
	SORT:   
	}
	######### ELSE des espaces / with spaces
	else {
	    if($x=~m/^(\d+)/){
		foreach my $i (sort {$a <=> $b} keys %seqhash) {
		    if($seqhash{$i}{ori_head}=~m/^(\d+),?_?\s/) { $seqhash{$i}{db_acc}=$1; $seqhash{$i}{ori_head}=~s/$seqhash{$i}{db_acc}//; $seqhash{$i}{ori_head}=~s/^,\s?//; $seqhash{$i}{ori_head}=~s/^_//;  }
		    elsif($seqhash{$i}{ori_head}=~m/^([^\s]+)/){ $seqhash{$i}{db_acc}=$1; $seqhash{$i}{ori_head}=~s/$seqhash{$i}{db_acc}//; $seqhash{$i}{ori_head}=~s/^,\s?//; $seqhash{$i}{ori_head}=~s/^_//;  }
		    if($seqhash{$i}{ori_head}=~m/^(\d+),?\s(.*)/) {$seqhash{$i}{locus_tag}=$1; $seqhash{$i}{db_desc}=$2;}
		} 
		goto SORT2;
	    } for( my $j=0; $j<$length_x; $j++ ) {
		$ct=1;
		$testlet[$j]=substr($x,0,$j+1);
		if ($testlet[$j]=~m/.*_(\d)/) {$val=$testlet[$j]; $next="yes"; goto SORT2; }
		for (my $i=0; $i<=$#headtmp; $i++ ){ if ($headtmp[$i]=~m/^$testlet[$j]/) {$ct++;} }        
		if($orf != $ct)  { $val=$testlet[$j]; $next="yes"; goto SORT2; }
	    } 
	SORT2:
	} 
	######### Var NEXT=YES traitement commun des 2 cas : avec et sans espaces  / both with/without spaces
	if ($next eq "yes") {
	    $val=$testlet[-2]; 
	    if($val=~m/(.*)\d\d*?$/) {$val=$1;}
	    foreach my $t (keys %seqhash) {
		if ($seqhash{$t}{ori_head}=~m/^$val([^\s,]+),?\s(.*)\n?$/){
		    $next="no";
		    $seqhash{$t}{db_acc}=$1;
		    $seqhash{$t}{db_desc}=$2;
		} elsif ($seqhash{$t}{ori_head}=~m/^$val(.*)\s?\n?$/){ 
		    $next="no";
		    $seqhash{$t}{db_acc}=$1;
		}
	    }
	}
    }
    #------------------------------------------------
    #          Verification du nb de seq traitee
    # -----------------------------------------------
    
    my $verif=0;
    foreach my $i (keys %seqhash) {
	$seqhash{$i}{db_acc}=~s/^_//; $seqhash{$i}{db_acc}=~s/_$//; $seqhash{$i}{db_acc}=~s/^\s//;
	if ($seqhash{$i}{locus_tag}) {
	    $seqhash{$i}{locus_tag}=~s/^_//; $seqhash{$i}{locus_tag}=~s/_$//; $seqhash{$i}{locus_tag}=~s/^\s//;
	}
	if ($seqhash{$i}{db_desc}) {
	    $seqhash{$i}{db_desc}=~s/^_//; $seqhash{$i}{db_desc}=~s/_$//; $seqhash{$i}{db_desc}=~s/^\s//;
	}
	$seqhash{$i}{description} = $seqhash{$i}{db_acc};
	$verif++;    
    }
    
    if ($orf != $verif) {
	return "Clean Geno Failed  Nb ORF = ".$orf." <--DIFFERENT--> Nb of processed sequences = ".$verif."<br>";
    }elsif ($seqhash{1}{db_acc} eq '' ) {
	return "Clean Geno Failed $x <br>";
    }
    return %seqhash;
}
sub delpat {
    my ($orf) = shift @_;
    my ($lgth) = shift @_;
    my (@head) = @_;
    my @pattern = ();
    my @pat = ();
    my $ctpat=0;
    for (my $i=0; $i<$lgth; $i++) { $pattern[$i] = $head[$i]; }
    for (my $i=0; $i<$lgth; $i++) { shift @head; }

    for (my $j=0; $j<=$#pattern; $j++)
     {  $ctpat=0;
        for( my $i=0; $i<=$#head; $i++ ) { if($head[$i] =~ m/$pattern[$j]/){$ctpat++;} }
        if ($ctpat==$orf)
	    { $head[0]=~s/$pattern[$j]:|:$pattern[$j]//;
	      $head[0]=~s/$pattern[$j]_|_$pattern[$j]//;
	      for( my $l=1; $l<=$#head; $l++ )
		  { $head[$l]=~s/$pattern[$j]:|:$pattern[$j]//;
		    $head[$l]=~s/$pattern[$j]_|_$pattern[$j]//; 
		} 
	  } 
    }         
    return @head;    
}
