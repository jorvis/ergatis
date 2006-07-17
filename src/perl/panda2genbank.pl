#!/usr/local/bin/perl
# modified version of generate_signature_sets.pl

# Sample execution:
# % panda2gb.pl --query_taxid=83334,83333 --background_taxid=543 --output_dir=/tmp

# taxon_file has two tab separated columns of comma-delimited taxon_ids
# the first column is the target set, the second column is the background, e.g.
# target1, target2, target3      background1, background2
# target3       background2
# etc...

# 0) Pull list of children taxon ids from SYBPANDA.prometheus
#    eg. get_gb_children 562;
#    for target and background sets
#    Separate target from background set.

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use DBI qw(:sql_types); #for SQL_NUMERIC datatype
umask(0000);

my %options = &parse_options();
my $storage = $options{output_dir};

if (defined ($options{input_file})) {
    print "Input from: ".$options{input_file}."\n";
    open(my $FIN, $options{input_file}) || print_usage("Unable to open taxon file ".$options{input_file});
    while (my $row = <$FIN>) {
	unless ($row =~ /^#/) { # use # for comments
		chomp $row;
		my @junk = split("\t",$row);
		(@junk == 2) || print_usage("Incorrect format on line $. in taxon file ".$options{input_file});
		generate_accession_set($junk[0], $junk[1]); #divide columns on tab, send to generate
	    }
    }
}
else {
    print "Input from: command line\n";
    generate_accession_set($options{query_taxid}, $options{background_taxid});
}

#
# Subroutines follow
#

# input: target and background taxon sets
# output: a signature.input file in $storage
sub generate_accession_set {
    my @target_source = (split(',',$_[0]));
    my @background_source = (split(',',$_[1]));
    my $setname = "t".join(".",@target_source)."_b".join(".",@background_source);
    my $ofile = "$storage/$setname.signature.input";
    print "Output file: $ofile\nset: $setname\n";
    my @target = get_descendents_from_taxon_id(@target_source);
    my @gestalt = get_descendents_from_taxon_id(@background_source);
    my @background = a_not_b(\@gestalt, \@target);

    #get and print out target taxons and accesions
    open(my $FOUT, ">$ofile") or die "Unable to create signature.input file";
    print {$FOUT} "\#[target: @target_source]";
    get_and_print_accessions($FOUT, @target);
    print {$FOUT} "\n\#[background: @background_source]";
    get_and_print_accessions($FOUT, @background);
}


#Connect to the Panda database to obtain taxonomy info
#returns tax_id node and all its descendents
#for each tax_id node passed to function
sub get_descendents_from_taxon_id {
    unless (@_ > 0) {
	die "Must provide at least 1 taxon_id to get_descendents_from_taxon_id";
    }
    my $db = "prometheus";
    my $dbh = connect_db('loadpanda');
    $dbh->do("use $db");

    my $query = 'exec get_gb_children @tax_id=?';
    my $sth = $dbh->prepare($query);

    #get descendents from all... order doesn't matter...don't want dupes
    my %kids;
    foreach (@_) {
	#skip if something in the input set was a child of something else already looked for
	unless (defined $kids{$_}) {
	    $sth->bind_param(1,$_,SQL_NUMERIC);#stored procedure needs numeric data type
		$sth->execute();
	    while (my $row = $sth->fetchrow_hashref) {
		$kids{$row->{'child_tax_id'}}++;
	    }
	}
    }
    $sth=undef;
    $dbh->disconnect; 

    return (keys %kids);
}


#input: two arrayrefs of arrays of *unique* items
#returns: array of elements from A not in B
sub a_not_b {
    my @A = @{$_[0]};
    my @B = @{$_[1]};

    #stolen from perl cookbook
    my %seen = ();  # lookup table to test membership of B
    my @aonly = (); # answer

    # build lookup table
    foreach my $item (@B) { $seen{$item} = 1 }

    # find only elements in @A and not in @B
    foreach my $item (@A) {
	unless ($seen{$item}) {
	    # it's not in %seen, so add to @aonly
	    push(@aonly, $item);
	}
    }  
    return @aonly;
}


#since we print out accessions for target and background the same way
#just use a function to standardize it
#format is
#taxon_id:acc1, acc2, ...
sub get_and_print_accessions {
    (my $FILE, my @taxons) = @_;
    my %rethash = %{get_best_accessions(@taxons)};
    foreach my $k (keys %rethash) {
	my @acc_files = ();
	foreach my $v (@{$rethash{$k}}) {
	    unless (-e "$storage/$v.gbk") {
		`/usr/local/common/yank_panda --genbank -a=$v > $storage/$v.gbk`;
	    }
	    if (`cat $storage/$v.gbk` =~ /^Found 0 results/) {
		warn "t:Unable to find gbk for $v.  Removing empty file.";
		`rm $storage/$v.gbk`;
	    }
	    else {
		push(@acc_files, $v);
	    }
	}
	if (@acc_files > 0) {
	    print {$FILE} "\n$k:".join(", ",@acc_files);
	}
    }
}

# Get the "best set" of accessions associated with a taxon id
# input: taxon_id(s)
# output: hashref, key is taxonid, values are arrays of accessions

# NOTE: Upon completion of bug 2717 this should be rewritten
# http://jorvis-lx:8080/bugzilla/show_bug.cgi?id=2717
# all of the get_XXX_genomes will be almost identical functions,
# varying only in the WHERE entry_info.xxx = 1 clause

sub get_best_accessions {
    my %tax_acc; #keys are input taxids, values are arrays of accessions
    foreach (@_) {
	$tax_acc{$_} = [];
    }

    foreach my $tid (@_) {
	#1: Check for complete genomes
	if (my @cg_accs = get_panda_cg($tid)) {
	    #print "\n\tcg:";
	    foreach my $acc (@cg_accs) {
		#print " $acc,";
		push(@{$tax_acc{$tid}}, $acc);
	    }
	}

	#CURRENTLY NO DATA
	#2: Check for wgs sequences from genome projects
#	elsif (my @wgs_accs = get_panda_wgs($tid)) {
#	    foreach my $acc (@wgs_accs) {
#		push(@{$tax_acc{$tid}}, $acc);
#	    }
#	}

	#3: Check for genbank records from genome projects
	#   ie: AE, CP, CY,(genbank) AL, BX, CR, CT, (embl) or AP (ddjb) accession prefixes
	elsif ( my @wgp_accs = get_panda_wgp($tid)) {
	    #print "\n\twgp:";
	    foreach my $acc (@wgp_accs) {
		#print " $acc,";
		push(@{$tax_acc{$tid}}, $acc);
	    }
	}

	else {
	    #???
	}
    }

    return \%tax_acc;
}

#1: get full refseq genomes
#   return accessions of full genomes associated w/ a taxon id(s)
sub get_panda_cg {
    my $dbh = connect_db('loadpanda');
    
    my %sth = (
	       #in the future, this will be the call
#	       pull_cg => 'SELECT gi, accession
#                          FROM genbank..entry_info
#                          WHERE cg = 1
#                          AND tax_id=?'
	       #welcome to the present...
	       pull_cg => 'SELECT n.accession_id, n.accession, count(*) as C
                           FROM prometheus..protein p, prometheus..protein_taxon pt, prometheus..cds c, prometheus..nucleotide n 
                           WHERE p.cg = 1
                           AND pt.prot_id = p.prot_id
                           AND c.cds_id = p.cds_id
                           AND n.nucl_id = c.nucl_id
                           AND pt.taxon_id = ?
                           GROUP BY n.accession_id, n.accession
                           ORDER BY C DESC'
	       );
    foreach (keys %sth) {
	$sth{$_} = $dbh->prepare($sth{$_});
    }

    my %genomes; #use hash because only want unique (not caring about versioning)
    #probably just one passed value, but maybe more
    foreach (@_) {
	#print "\npull_cg($_)";
	$sth{pull_cg}->execute($_);
	foreach my $row (@{$sth{pull_cg}->fetchall_arrayref()}) {
	    #row->[0] is gi, row->[1] is accession
	    $genomes{$row->[1]}++;
	}
    }

    #close up shop
    foreach (keys %sth) {$sth{$_}->finish();}
    $dbh->disconnect;

    return keys %genomes;
}

#2: WGS sequence
#   Currently known to not contain any data
# input: taxon_id(s)
# output: accessions of all entries in the 
# wgs_genomeprj table from given taxon_id(s)
sub get_panda_wgs {
    my $dbh = connect_db('loadpanda');
    
    my %sth = (
	       #in the future, this will be the call
	       pull_wgs => 'SELECT gi, accession
                            FROM genbank..entry_info
                            WHERE wgs = 1
                            AND tax_id=?'
	       );
    foreach (keys %sth) {
	$sth{$_} = $dbh->prepare($sth{$_});
    }

    my %genomes; #use hash because only want unique (not caring about versioning)
    #probably just one passed value, but maybe more
    foreach (@_) {
	$sth{pull_wgs}->execute($_);
	foreach my $row (@{$sth{pull_wgs}->fetchall_arrayref()}) {
	    #row->[0] is gi, row->[1] is accession
	    $genomes{$row->[1]}++;
	}
    }

    #close up shop
    foreach (keys %sth) {$sth{$_}->finish();}
    $dbh->disconnect;

    return keys %genomes;
}


#3: genbank records from genome projects
#goes against panda database, which may not be as complete as it should?
#also, is slow
sub get_panda_wgp {
 my $dbh = connect_db('loadpanda');
    
    my %sth = (
	       #in the future, this will be the call
#	       pull_wgp => 'SELECT gi, accession
#                          FROM genbank..entry_info
#                          WHERE cg = 1
#                          AND tax_id=?'
    # Keep accessions beginning with letters AE, CP, CY,(genbank) AL, BX, CR, CT, (embl) AP (ddjb)
    # http://www.ncbi.nlm.nih.gov/Sequin/acc.html
	       pull_wgp => 'SELECT gi, accession 
                           FROM genbank..entry_info
                           WHERE tax_id = ? 
                           AND (accession like "AE%" OR
              	            accession like "AE%" OR 
   	                    accession like "CP%" OR
   	                    accession like "CY%" OR
   	                    accession like "AL%" OR 
   	                    accession like "BX%" OR
   	                    accession like "CR%" OR 
   	                    accession like "CT%" OR
 	           	    accession like "AP%" )'
	       );
    foreach (keys %sth) {
	$sth{$_} = $dbh->prepare($sth{$_});
    }

    my %genomes; #use hash because only want unique (not caring about versioning)
    #probably just one passed value, but maybe more
    foreach (@_) {
	#print "\npull_wgp($_)";
	$sth{pull_wgp}->execute($_);
	foreach my $row (@{$sth{pull_wgp}->fetchall_arrayref()}) {
	    #row->[0] is gi, row->[1] is accession
	    $genomes{$row->[1]}++;
	}
    }

    #close up shop
    foreach (keys %sth) {$sth{$_}->finish();}
    $dbh->disconnect;

    return keys %genomes;
}

sub print_usage
{
        my $byebye = shift;
        my $progname = $0;
        die << "END";
$byebye
usage: $progname --output_dir <directory> --input_file <file>
or provide taxids separated by commas (no space):
$progname --output_dir <directory> --query_taxid=<taxid1,taxid2..> --background_taxid=<taxidA,taxidB,...>

Input is target and background taxid(s).
Output is a signature.input file listing the genbank id(s) for each taxid.

Input taxids can be given on the command line like:
  panda2gb.pl --query_taxid=83334,83333 --background_taxid=543 --output_dir=/tmp
or via a file:
  $progname --output_dir=/tmp --input_file=./set_lists
where set_lists would be formatted to describe one or more runs:
  #target background
  83334,83333    543
  562    543

The output signature.input file will have the format:
  #[target: 132564]
  <taxid>:<acc>,<acc>,<acc>...
  <taxid>:<acc>,<acc>...
  #[background: 193568]
  <taxid>:<acc>,<acc>...
  etc...
END
}

sub parse_options {
    my %options = ();
    GetOptions( \%options,
                'help|h',
                'output_dir|o=s',
                'input_file|i=s',
                'query_taxid|q=s',
		'background_taxid|b=s'
                ) || &print_usage("Unprocessable option");

    # check for required parameters
    (defined $options{help}) && print_usage("How to use the program:");
    (defined $options{output_dir}) || print_usage("output_dir required");
    (-r $options{output_dir}) || print_usage("output_dir ($options{outout_dir}) not readable");
 
    if (defined $options{input_file}) {
	    (-r $options{input_file}) || print_usage("input_file ($options{input_file}) not readable");
    }

    print "Executing $0 with options\n";
    foreach (keys %options) { print "  $_: $options{$_}\n";}
 
    return %options;
}


sub connect_db {
    my $dbname = shift;
    if (lc($dbname) eq 'loadpanda') {
	my $server = "LOADPANDA";
	my $username = "access";
	my $password = "access";
    
	my $dbh = DBI->connect("dbi:Sybase:$server", $username,$password)
            or die "Unable to connect to $dbname";

        return $dbh;
    }
    else {
        die "Unexpected database, cannot connect";
    }
}
