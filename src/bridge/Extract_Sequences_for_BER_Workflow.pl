#! /usr/local/bin/perl

=head2
	Extract_Sequences_for_BER_Workflow.pl
		
=text

Paolo Amedeo, June 28th 2005

This script extracts protein sequence and underlining, extended nucleotide sequence
for running BER searches.
It also output a table with the corresponding names of the sequences of the two fasta
files (in this case identical).

=cut

use strict;
use warnings;
use lib '/local/perl/lib/site_perl';
use DBI;
use Getopt::Std;

our ($opt_D, $opt_a, $opt_L, $opt_l, $opt_e, $opt_O, $opt_h);
getopts('D:a:L:l:e:O:h');


MAIN:{
	my ($prN) = ($0 =~ /\/?([^\/]+)$/);
	my $message = "\n\nUsage:  $prN  <options>\n\n";
	my $bad = 0;
	my $ext_ln = 300;
	my @asmbls = ();
	my @models = ();
	
	my ($db, $out_name);

	if (defined $opt_D && $opt_D =~ /\w/){
		$db = $opt_D;
	} else {
		$message .= "Option -D is required\n\n" unless $opt_h;
		++$bad;
	}

	if (defined $opt_a && defined $opt_L || defined $opt_a && defined $opt_l || defined $opt_l && defined $opt_L){
		$message .= "Options -a and -L of -l  are incompatible: chose one out of the three\n\n";
		++$bad;
	}

	if (defined $opt_a && $opt_a =~ /\d/ && $opt_a !~ /\D/){
		@asmbls = ($opt_a);
	} 
	elsif (defined $opt_a){
		$message .= "Bad value for option -a (assembly ID)\n\n";
		++$bad;
	}

	if (defined $opt_L && open(LIST, $opt_L)){
		while (<LIST>){
			push(@asmbls, $1) if /(\d+)/;
		}
		close(LIST);
	}
	elsif (defined $opt_L){
		$message .= "Bad value for option -L (list of assemblies) or impossible to open the file\n\n";
		++$bad;
	}

	if (defined $opt_l && open(LIST, $opt_l)){
		while (<LIST>){
			push(@models, $1) if /(\d+\.m\d+)/;
		}
		close(LIST);
	}
	elsif (defined $opt_l){
		$message .= "Bad value for option -l (list of models) or impossible to open the file\n\n";
		++$bad;
	}

	if (defined $opt_O && $opt_O =~ /\S/ && $opt_O !~ /\s/){
		$out_name = $opt_O;
	}
	elsif (defined $opt_O){
		$message .= "Bad value for option -O (Name for the Output fasta files)\n\n";
		++$bad;
	} else {
		$message .= "Option -O (Name for the Output fasta files) is required\n\n" unless $opt_h;
		++$bad;
	}


	$message .= "


###############################################################################
#
# -D annotation database (mandatory) 
#
# -a assembly ID
#
# -L list of assembly IDs
#
# -l list of models
#  **if none of -a, -l or -L are specified, it will act on all current assemblies
# 
# -e Extension length (default: $ext_ln nt)
#
# -O Name for the Output files
#
# -h print this option menu and quit
#
###############################################################################

";


	die $message if $opt_h || $bad;
	
	$ext_ln = $opt_e if defined $opt_e && $opt_e =~ /\d/ && $opt_e !~ /\D/;
	
	open(DB, "$ENV{EGC_SCRIPTS}/egc_password") || die "\n\nImpossible to open the file with information for database connection\n\n";
	chomp(my ($user, $pass) = <DB>);
	close(DB);
	
	my $dbh1 = DBI->connect("dbi:Sybase:server=SYBTIGR; packetSize=8092",$user, $pass) || die "\n\nProblems with database connection\n\n";
	$dbh1->do("use $db");
	
	open(my $aa, ">$out_name.aa.fsa")   || die "\n\nImpossible to open the file $out_name.aa.fasta for writing\n\n";
	open(my $nt, ">$out_name.nt.fsa")   || die "\n\nImpossible to open the file $out_name.nt.fasta for writing\n\n";
	open(my $table, ">$out_name.table") || die "\n\nImpossible to open the file $out_name.table for writing\n\n";

	my $modl_done;
	
	if (@asmbls){
		my $get_models_qry = "SELECT a.feat_name, a.end5, a.end3, a.protein FROM asm_feature a, phys_ev p WHERE a.feat_name = p.feat_name  AND a.asmbl_id = ? AND a.feat_type = 'model' AND p.ev_type = 'working'";
		my $get_models = $dbh1->prepare($get_models_qry) || die "\n\nProblems preparing the following query:\n\"$get_models_qry\"\n\n";
		$modl_done = &WriteModels($dbh1, $db, $nt, $aa, $table, $get_models, $ext_ln, \@asmbls) || die "\n\nNo models written\n\n";
	}
	elsif (@models){
		my $get_models_qry = "SELECT a.feat_name, a.end5, a.end3, a.protein FROM asm_feature a WHERE a.feat_name =  ?";
		my $get_models = $dbh1->prepare($get_models_qry) || die "\n\nProblems preparing the following query:\n\"$get_models_qry\"\n\n";
		$modl_done = &WriteModels($dbh1, $db, $nt, $aa, $table, $get_models, $ext_ln, \@models) || die "\n\nNo models written\n\n";
	} else {
		my $get_models_qry = "SELECT a.feat_name, a.end5, a.end3, a.protein FROM asm_feature a, phys_ev p, clone_info c WHERE a.feat_name = p.feat_name  AND a.asmbl_id = c.asmbl_id AND c.is_public = 1 AND a.feat_type = 'model' AND p.ev_type = 'working'";
		my $get_models = $dbh1->prepare($get_models_qry) || die "\n\nProblems preparing the following query:\n\"$get_models_qry\"\n\n";
		$modl_done = &WriteModels($dbh1, $db, $nt, $aa, $table, $get_models, $ext_ln) || die "\n\nNo models written\n\n";
	}
		
	print "\n\nWritten $modl_done models\n\n";
	
	close($aa);
	close($nt);
	close($table);
	
	$dbh1->disconnect;
}

sub WriteModels {
	my ($dbh, $annot_db, $nt_file, $aa_file, $table_file, $sth, $ext, $r_list) = @_;
	my $processed = 0;
	
	if (defined $r_list){
		foreach my $elem (@{$r_list}){
			$sth->execute($elem) || die "\n\nProblems executing the query\n\n";
			
			$processed += &Fetch_n_Write($dbh, $sth, $annot_db, $nt_file, $aa_file, $table_file, $ext) || warn "Problems fetching the results or no models for \"$elem\"\n";
		}
	} else {
		$sth->execute() || die "\n\nProblems executing the query\n\n";
		$processed += &Fetch_n_Write($dbh, $sth, $annot_db, $nt_file, $aa_file, $table_file, $ext) || warn "\n\nProblems fetching the results or no models\n\n";
	}
	return($processed);
}

sub Fetch_n_Write {
	my ($dbh, $prepped, $dbase, $ntfile, $aafile, $tblfile, $ext_ln) = @_;
	my $written = 0;
	my ($mod_name, $e5, $e3, $prot);
	my %modl = ();
	
	$prepped->bind_columns(undef, \$mod_name, \$e5, \$e3, \$prot);
	
	while ($prepped->fetch()){
		unless (defined $prot && $prot =~ /\w/){
			warn "Model $mod_name lacks Protein sequence - skept\n";
			next;
		}
		my ($asmbl) = ($mod_name =~ /(\d+)/);
		push(@{$modl{$asmbl}}, [$mod_name, $e5, $e3, $prot]);
	}
	
	my $max_ln_txt;
	my $qry = "SELECT MAX(length) FROM clone_info WHERE is_public = 1";
	my $sth = $dbh->prepare($qry) || die "\n\nProblem preparing the following query:\n\"$qry\"\n\n";
	$sth->execute() || die "\n\nProblem executing the following query:\n\"$qry\"\n\n";
	$sth->bind_columns(undef, \$max_ln_txt);
	$sth->fetch();
	$sth->finish();
	die "\n\nProblems retrieving info from the database\n\n" unless defined $max_ln_txt;
	$max_ln_txt += 10;
	$dbh->do("SET TEXTSIZE $max_ln_txt") || die "\n\nImpossible to set the textsize variable to $max_ln_txt\n\n";
	my $get_asmbl_seq_qry = "SELECT sequence FROM assembly WHERE asmbl_id = ?";
	my $get_asmbl_seq = $dbh->prepare($get_asmbl_seq_qry)  || die "\n\nProblem preparing the following query:\n\"$get_asmbl_seq_qry\"\n\n";
	
	while (my ($asmbl, $models) = each %modl){
		my $asmbl_seq;
		$get_asmbl_seq->execute($asmbl) || die "\n\nProblem preparing the following query:\n\"$get_asmbl_seq_qry\" ('$asmbl')\n\n";
		$get_asmbl_seq->bind_columns(undef, \$asmbl_seq);
		$get_asmbl_seq->fetch();
		$get_asmbl_seq->finish();
		die "Sequence for assembly $asmbl is not in the database\n\n" unless defined $asmbl_seq;
		
		foreach my $model (@{$models}){
			my ($start, $end, $minusstrand) = $model->[1] < $model->[2] ? (@{$model}[1,2], 0) : (@{$model}[2,1], 1);
			(my $model_name = $model->[0]) =~ s/\.m/_/;
			my $header = "$dbase.model.$model_name";
			print {$tblfile} "$header\t$header\n";
			
			$start = $start < $ext_ln ? 0 : $start - $ext_ln - 1;
			my $ln = 1 + $end + $ext_ln - $start;
			my $nucl = substr($asmbl_seq, $start, $ln);
			&ComplementDNA(\$nucl) || die "\n\nProblems complementing the nucleotide sequence of model $model->[0]\n\n" if $minusstrand;
		
			foreach my $seq ($nucl, $prot){
				$seq =~ s/\W+//g;
				$seq =~ s/(.{1,60})/$1\n/g;
			}
			print {$ntfile} ">$header\n$nucl";
			print {$aafile} ">$header\n$prot";
			++$written;
		}
	}
	return($written);
}



sub ComplementDNA {
	my ($r_seq) = @_;
	${$r_seq} =~ tr/AaCcGgTtMmRrWwSsYyKkVvHhDdBb/TtGgCcAaKkYyWwSsRrMmBbDdHhVv/;
	${$r_seq} = reverse(${$r_seq});
	return (length(${$r_seq}));
}

