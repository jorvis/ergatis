#! /usr/local/bin/perl

=head2
	Create_Assemblies_Fasta.pl
		
=text

Paolo Amedeo, March 31st 2005

This program extracts the sequence of a list of assemblies from the database

=cut

use strict;
use warnings;
use lib ('/local/perl/lib/site_perl');
use DBI;
use Getopt::Std;

our ($opt_D, $opt_a, $opt_L, $opt_o, $opt_h);
getopts('D:a:L:o:h');

MAIN:{
	my ($prN) = ($0 =~ /([^\/]+)$/);
	my $db_file = "$ENV{EGC_SCRIPTS}/egc_password";
	my $message = "\n\nUsage:  $prN  <options>\n\n";
	my $bad = 0;
	my @asmbls = ();
	my ($db, $out_file);

	if (defined $opt_D && $opt_D =~ /\w/){
		$db = $opt_D;
	} else {
		$message .= "Option -D is required\n\n" unless $opt_h;
		++$bad;
	}
	if (defined $opt_a){
		if (defined $opt_L){
			$message .= "Options -a and -L are incompatible: chose one out of the two\n\n";
			++$bad;
		}
		elsif($opt_a =~ /\d/ && $opt_a !~ /\D/){
			@asmbls = ($opt_a);
		} else {
			$message .= "Bad value for option -a (assembly ID)\n\n";
			++$bad;
		}
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
	if (defined $opt_o && $opt_o =~ /\S/){
		$out_file = $opt_o;
	} else {
		$message .= "Options -o (Name for output file) is mandatory\n\n" unless $opt_h;
		++$bad;
	}

	$message .= "


###############################################################################
#
# -D annotation database (required) 
#
# -a assembly ID
#
# -L list of assembly IDs
#  **if neither -a or -L are specified, it will act on all current assemblies
# 
# -o output file (required)
#
# -h print this option menu and quit
#
###############################################################################

";


	die $message if $opt_h || $bad;

	open(OUT, ">$out_file") || die "\n\nImpossible to open the file $out_file for writing\n\n";

	open(DB, $db_file) || die "\n\nImpossible to open the file $db_file freading\n\n";
	chomp(my ($user, $pass) = <DB>);
	close(DB);
	
	my $dbh1 = DBI->connect("dbi:Sybase:server=SYBTIGR; packetSize=8092",$user, $pass) || die "\n\nProblems with database connection\n\n";
	$dbh1->do("use $db");

	unless (@asmbls){
		my $asmbl_id;
		my $qry = "SELECT asmbl_id from clone_info where is_public = 1";
		my $sth = $dbh1->prepare($qry);
		$sth->execute() || die "\n\nProblems executing the following query:\"$qry\"\n\n";
		$sth->bind_columns(undef, \$asmbl_id);
		
		while ($sth->fetch()){
			push(@asmbls, $asmbl_id);
		}
	}
	
	# Getting the length and name of the clone....
	
	my $get_asmbl_info_qry = "SELECT clone_name, length FROM clone_info WHERE asmbl_id = ?";
	my $get_asmbl_info = $dbh1->prepare($get_asmbl_info_qry) || die "\n\nProblems preparing the following query:\n\"$get_asmbl_info_qry\"\n\n";
	my $max_ln = 0;
	
	foreach my $asmbl (@asmbls){
		my ($clone_name, $length);
		$get_asmbl_info->execute($asmbl) || die "\n\nProblems executing the following query:\n\"$get_asmbl_info_qry\"\n\n";
		$get_asmbl_info->bind_columns(undef, \$clone_name, \$length);
		$get_asmbl_info->fetch();
		$get_asmbl_info->finish();
		$max_ln = $length if $length > $max_ln;
		$asmbl = [$asmbl, $length, $clone_name];
	}
	
	# Getting the sequence of the assemblies...
	
	$dbh1->do("set textsize $max_ln");
	
	my $get_asmbl_qry = "SELECT sequence FROM assembly WHERE asmbl_id = ?";
	my $get_asmbl_seq = $dbh1->prepare($get_asmbl_qry) || die "\n\nProblems preparing the following query:\n\"$get_asmbl_qry\"\n\n";
	
	foreach my $asmbl (@asmbls){
		print STDERR "Retrieving the sequence for assembly $asmbl->[0]\n";
		my $seq;
		$get_asmbl_seq->execute($asmbl->[0]) || die "\n\nProblems executing the following query:\n\"$get_asmbl_qry\"\n\n";
		$get_asmbl_seq->bind_columns(undef, \$seq);
		$get_asmbl_seq->fetch();
		$get_asmbl_seq->finish();
		
		$seq =~ s/\W+//g;
		$seq =~ s/(.{1,60})/$1\n/g;
		
		print OUT ">$db.assembly.$asmbl->[0]\n$seq";
	}
	$dbh1->disconnect;
	close(OUT);
}
