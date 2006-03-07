#! /usr/local/bin/perl

=head2
	Extract_Models_for_WorkFlow.pl
		
=text

Paolo Amedeo, April 22nd 2005

This script is based on Extract_Model_List.pl with the following few differences:

1) usage of options (Getopt::Std)
2) writes as header the names of the models according with WorkFlow specs:
	db.molecule.feat# (i.e. $db.model.$asmbl+5_digits_model_number

=cut

use strict;
use warnings;
use lib '/local/perl/lib/site_perl';
use DBI;
use Getopt::Std;

our ($opt_D, $opt_a, $opt_L, $opt_l, $opt_O, $opt_h, $opt_s);
getopts('D:a:L:l:O:h:s');


MAIN:{
	my ($prN) = ($0 =~ /\/?([^\/]+)$/);
	my $message = "\n\nUsage:  $prN  <options>\n\n";
	my $bad = 0;
	my @asmbls = ();
	my @models = ();
	
	my ($db, $out_name);

    if (! defined $opt_s ) {
        ## default length cutoff
        $opt_s = 1;
    }

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
		my %seen = ();
        
        while (<LIST>){
			next unless /(\d+)/ &! exists $seen{$1};
            push(@asmbls, $1);
            undef $seen{$1};
		}
		close(LIST);
	}
	elsif (defined $opt_L){
		$message .= "Bad value for option -L (list of assemblies) or impossible to open the file\n\n";
		++$bad;
	}

	if (defined $opt_l && open(LIST, $opt_l)){
		my %seen = ();
        
        while (<LIST>){
			next unless /(\d+\.m\d+)/ &! exists $seen{$1};
            push(@models, $1);
            undef $seen{$1};
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
# -O Name for the Output multi-fasta file
#
# -s size cutoff.  sequences with protein seqs shorter than this will not be included.
#
# -h print this option menu and quit
#
###############################################################################

";


	die $message if $opt_h || $bad;

	open(DB, "$ENV{EGC_SCRIPTS}/egc_password") || die "\n\nImpossible to open the file with information for database connection\n\n";
	chomp(my ($user, $pass) = <DB>);
    #($user, $pass) = ('access', 'access');
	close(DB);
	
	my $dbh1 = DBI->connect("dbi:Sybase:server=SYBTIGR; packetSize=8092",$user, $pass) || die "\n\nProblems with database connection\n\n";
	$dbh1->do("use $db");
	
	open(my $aa, ">$out_name.aa.fsa") || die "\n\nImpossible to open the file $out_name.aa.fasta for writing\n\n";
	open(my $nt, ">$out_name.nt.fsa") || die "\n\nImpossible to open the file $out_name.nt.fasta for writing\n\n";

	my $modl_done;
	
	if (@asmbls){
		my $get_models_qry = "SELECT a.feat_name, a.sequence, a.protein FROM asm_feature a, phys_ev p WHERE a.feat_name = p.feat_name  AND a.asmbl_id = ? AND a.feat_type = 'model' AND p.ev_type = 'working'";
		my $get_models = $dbh1->prepare($get_models_qry) || die "\n\nProblems preparing the following query:\n\"$get_models_qry\"\n\n";
		$modl_done = &WriteModels($db, $nt, $aa, $get_models, \@asmbls) || die "\n\nNo models written\n\n";
	}
	elsif (@models){
		my $get_models_qry = "SELECT a.feat_name, a.sequence, a.protein FROM asm_feature a WHERE a.feat_name =  ?";
		my $get_models = $dbh1->prepare($get_models_qry) || die "\n\nProblems preparing the following query:\n\"$get_models_qry\"\n\n";
		$modl_done = &WriteModels($db, $nt, $aa, $get_models, \@models) || die "\n\nNo models written\n\n";
	} else {
		my $get_models_qry = "SELECT a.feat_name, a.sequence, a.protein FROM asm_feature a, phys_ev p, clone_info c WHERE a.feat_name = p.feat_name  AND a.asmbl_id = c.asmbl_id AND c.is_public = 1 AND a.feat_type = 'model' AND p.ev_type = 'working'";
		my $get_models = $dbh1->prepare($get_models_qry) || die "\n\nProblems preparing the following query:\n\"$get_models_qry\"\n\n";
		$modl_done = &WriteModels($db, $nt, $aa, $get_models) || die "\n\nNo models written\n\n";
	}
		
	print "\n\nWritten $modl_done models\n\n";
	
	close($aa);
	close($nt);
	
	$dbh1->disconnect;
}

sub WriteModels {
	my ($annot_db, $nt_file, $aa_file, $sth, $r_list) = @_;
	my $processed = 0;
	
	if (defined $r_list){
		foreach my $elem (@{$r_list}){
			$sth->execute($elem) || die "\n\nProblems executing the query\n\n";
			
			$processed += &Fetch_n_Write($sth, $annot_db, $nt_file, $aa_file) || warn "\n\nProblems fetching the results or no models for \"$elem\"\n\n";
		}
	} else {
		$sth->execute() || die "\n\nProblems executing the query\n\n";
		$processed += &Fetch_n_Write($sth, $annot_db, $nt_file, $aa_file) || warn "\n\nProblems fetching the results or no models\n\n";
	}
	return($processed);
}

sub Fetch_n_Write {
	my ($prepped, $dbase, $ntfile, $aafile) = @_;
	my $written = 0;
	my ($mod_name, $nucl, $prot);
	
	$prepped->bind_columns(undef, \$mod_name, \$nucl, \$prot);
	
	while ($prepped->fetch()){
		unless (defined $nucl && $nucl =~ /\w/ && defined $prot && $prot =~ /\w/){
			warn "Model $mod_name lacks Nucleotide and / or Protein sequence - skipped\n";
			next;
		}
		$mod_name =~ s/\.m/_/;
		my $header = ">$dbase.model.$mod_name";
		
		foreach my $seq ($nucl, $prot){
			$seq =~ s/\W+//g;
			$seq =~ s/(.{1,60})/$1\n/g;
		}
        
        if ($opt_s && length $prot >= $opt_s ) {
		    print $ntfile "$header\n$nucl";
		    print $aafile "$header\n$prot";
		    ++$written;
        }
	}
	return($written);
}
