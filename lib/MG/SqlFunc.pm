package MG::SqlFunc;

use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use Data::Dumper;
use MG::Common;
use DBI;

@ISA = qw(Exporter);
@EXPORT = qw(&populate_table &get_general_taxclass &get_general_tax &get_cds_seq &get_data &update_table &get_lca_by_gi &get_lca_by_taxname);

sub populate_table {
    my ($dbh,$table, $data) = @_;
    foreach (@{$data}) {
	my %d = %{$_};
	my $columns = join(",", keys %d);
	my $values = join("','", values %d);
	$insert = $dbh->do(qq{insert into $table ($columns) VALUES (\'$values\')});
    }
    return 1;
}
sub update_table {
    my ($dbh,$table, $data,$where) = @_;
    my %where = %{$where};
    my @w = ();
    foreach my $k (keys %where) {
	push @w, "$k = $where{$k}"
    }
    $wsment = join(" and ", @w);
    my @set = ();
    foreach my $k (keys %{$data}) {
	push @set, "$k = $data->{$k}"
    }
    $set = join(",", @set);
    $update = $dbh->do(qq{update $table set $set where $wsment});
}
sub get_data {
    my ($dbh,$table, $data,$where) = @_;
    my %where = %{$where};
    my @w = ();
    foreach my $k (keys %where) {
	push @w, "$k = $where{$k}"
    }
    $wsment = join(" and ", @w);
    my $select = join(",",@{$data});
    my $info = $dbh->selectall_arrayref(qq{select $select from $table where $wsment});
    return @{$info};
}

sub get_general_taxclass {
    my ($dbh, $tax_id) = @_;
    $get_class = $dbh->selectall_arrayref(qq{select tax_class.tax_id, tax_class.name from tax_class inner join tax_node as t1 using(tax_id), tax_node as t2 where t2.tax_left > t1.tax_left and t2.tax_left < t1.tax_right and t2.tax_id = $tax_id});
    my ($taxid, $name) = @{$get_class->[0]};
    return ($taxid, $name);
}

sub get_general_tax {
    my $dbh = shift @_;
    my @tids = @_;
    my @taxa = ();
    my @where;
    foreach (@tids) {
	push @where, "tax_id = $_";
    }
    $where = join(" or ", @where);
    my ($min_left, $max_right) = $dbh->selectrow_array(qq{select min(tax_left), max(tax_right) from tax_node where $where});
    my ($tax_id, $tax_name) = $dbh->selectrow_array(qq{select tax_node.tax_id, tax_name.name from tax_node inner join tax_name using(tax_id) where tax_name.class = 'scientific name' and tax_left <=$min_left and tax_right >= $max_right order by tax_left DESC limit 1});
    return ($tax_id, $tax_name)
}
sub get_cds_seq {
    my ($dbh, $where) = @_;
    my @cond;
    my %info = ();
    foreach $col (keys %{$where}) {
	if ($where->{$col} =~ m/^like|>|</) {
	    push @cond, "$col $where->{$col}";
	}else {
	    push @cond, "$col = $where->{$col}";
	}
    }
    my $cond = join(" and ", @cond);
    my $info  = $dbh->selectall_arrayref(qq{select extra_acc,db_acc,genome_start,genome_end,genome_strand,seq_dna.sequence,seq_prot.sequence from ref_cds inner join seq_dna using(seq_dna_id) inner join seq_prot using (seq_prot_id) where $cond});
    my $i = 0;
    foreach $row (@{$info}) {
	$i ++;
	my ($dna_acc, $prot_acc, $start, $end, $strand, $dna_seq, $prot_seq) = @{$row};
	$info{$prot_acc} = {prot_acc=>$prot_acc,dna_acc=>$dna_acc,start=>$start,end=>$end,strand=>$strand,dna_seq=>$dna_seq, prot_seq=>$prot_seq};
    }
    return %info;
}
sub get_all_cds {
    my ($dbh, $where) = @_;
    my @cond = ();
    my @info = ();
    foreach $col (keys %{$where}) {
	if ($where->{$col} =~ m/^like|>|</) {
	    push @cond, "$col $where->{$col}";
	}else {
	    push @cond, "$col = $where->{$col}";
	}
    }
    my $cond = join(" and ", @cond);
    my $info  = $dbh->selectall_arrayref(qq{select extra_acc,db_acc,genome_start,genome_end,genome_strand,seq_dna.sequence,seq_prot.sequence from ref_cds using(ref_assembly_id) inner join seq_dna using(seq_dna_id) inner join seq_prot using (seq_prot_id) where $cond});
    my $i = 0;
    foreach $row (@{$info}) {
	$i ++;
	my ($dna_acc, $prot_acc, $start, $end, $strand, $dna_seq, $prot_seq) = @{$row};
	push @info, {prot_acc=>$prot_acc,dna_acc=>$dna_acc,start=>$start,end=>$end,strand=>$strand,dna_seq=>$dna_seq, prot_seq=>$prot_seq};
    }
    return @info;
}

sub get_lca_by_gi {
    my $dbh = shift @_;
    my @gis = @_;
    my @where;
    foreach (@gis) {
	push @where, "gi = $_" if $_;
    }
    $where = join(" or ", @where);
    my @tids = ();
    #warn("select distinct tax_id from gi_tax where $where");
    $result = $dbh->selectall_arrayref(qq{select distinct tax_id from gi_tax where $where});
    foreach (@{$result}) {
	push @tids, $_->[0];
    }
    if ($tids[0]) {
	my ($lca_tid, $lca_name) = get_general_tax($dbh,@tids);
	return($lca_tid, $lca_name);
    }else {
	return (0,'no tids');
    } 
}
sub get_lca_by_taxname {
    my $dbh = shift @_;
    my @names = @_;
    my @where;
    foreach (@names) {
	$_ =~ s/\'/\\\'/g;
	push @where, "name = \'$_\'"if $_;
    }
    $where = join(" or ", @where);
    my @tids = ();
    #warn("select distinct tax_id from tax_name where $where");
    $result = $dbh->selectall_arrayref(qq{select distinct tax_id from tax_name where $where});
    foreach (@{$result}) {
	push @tids, $_->[0];
    }
    if ($tids[0]) {
	my ($lca_tid, $lca_name) = get_general_tax($dbh,@tids);
	return($lca_tid, $lca_name);
    }else {
	return (0,'no tids');
    }
}

return 1;
