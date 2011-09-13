=head1 NAME

GiTaxon

=head1 DESCRIPTION

Utility to look up taxon information in a mongo database.

=head1 AUTHOR

David Riley

driley@som.umaryland.edu

=cut

package GiTaxon;
use strict;
use MongoDB;
use Bio::DB::Taxonomy;



sub new {
    my ($class, $args) = @_;

    my $self = {};
    $self->{'nodes'} = $args->{'nodes'} ? $args->{'nodes'} : '/local/db/by_source/ncbi/taxonomy/latest/nodes.dmp';
    $self->{'names'} = $args->{'names'} ? $args->{'names'} : '/local/db/by_source/ncbi/taxonomy/latest/names.dmp';
    $self->{'gi2tax'} = $args->{'gi2tax'} ? $args->{'gi2tax'} : '/local/db/by_source/ncbi/taxonomy/latest/gi_taxid_nucl.dmp';
    $self->{'chunk_size'} = $args->{'chunk_size'} ? $args->{'chunk_size'} : 10000;
    $self->{'idx_dir'} = $args->{'idx_dir'} ? $args->{'idx_dir'} : '/tmp/';
    $self->{'host'} = $args->{'host'} ? $args->{'host'} : 'tettelin-lx.igs.umaryland.edu';
    $self->{'gi_db'} = $args->{'gi_db'} ? $args->{'gi_db'} : 'lgt';
    $self->{'gi_coll'} = $args->{'gi_coll'} ? $args->{'gi_coll'} : 'gi2taxonnuc'; 
    print STDERR "Reading the NCBI Taxonomy $self->{idx_dir}\n";
    $self->{'db'} = Bio::DB::Taxonomy->new(-source => 'flatfile',
                                 -nodesfile => $self->{'nodes'},
                                 -namesfile => $self->{'names'},
                                 -directory => $self->{'idx_dir'});
    $self->{'cache'} = {};
    bless $self;
    $self->{'gi2taxon'} = $self->getgi2taxon($self->{'gi2tax'});

    return $self;
}

sub getTaxon {
    my($self,$acc) = @_;
    $acc =~ s/^\s+//;
    $acc =~ s/\s/\t/;
    my $gi = $acc;
    if($acc =~ /\|/) {
        my @fields = split(/\|/, $acc);
        $gi = $fields[1]; 
    }
    my $taxonid = '';
    my $retval = {};

    # First check the cache
    if($self->{cache}->{$gi}) {
        $retval = $self->{cache}->{$gi};
    }
    else {
        my $taxon_lookup = $self->{'gi2taxon'}->find_one({'gi' => "$gi"}, {'taxon' => 1});
        
        
        if($taxon_lookup) {
            $taxonid = $taxon_lookup->{'taxon'};
        }
        else {
#        print STDERR "Unable to find taxon for $gi\n";
        }
        
        my $taxon = $self->{'db'}->get_taxon(-taxonid => $taxonid);
        if(!$taxon) {
#        print STDERR "Unable to find taxon for $gi\n";
        }
        if(!$taxon) {
            $retval = {
                'acc'=> $acc,
                'gi'=> $gi,
                'taxon_id'=> $taxonid
            };
        }
        elsif ($taxon->isa('Bio::Taxon')) {
            my $name = $taxon->scientific_name;
            my $c = $taxon;
            my @lineage = ($name);
            while (my $parent = $self->{'db'}->ancestor($c)) {
                unshift @lineage, $parent->scientific_name;
                $c = $parent;
            }
            $retval = {
                'gi' => $gi,
                'acc'=> $acc,
                'taxon_id'=> $taxonid, 
                'name'=> $name, 
                'lineage'=> join(";", @lineage)
            };
        }else {
            print STDERR "Had something other than a Bio::Taxon\n";
#	print join("\t", $taxonid),"\n";
        }
        $self->{cache}->{$gi} = $retval;
    }
    return $retval;
}
sub getgi2taxon {
    my($self, $data_file) = @_;

    my $mongo = $self->get_mongodb_connection($self->{'gi_db'},$self->{'host'});
    my $coll = $mongo->get_collection($self->{'gi_coll'});

    print STDERR "Checking collection $self->{'gi_coll'}\n";
    if(!$coll->find_one()) {
        print STDERR "No records found in $self->{'gi_coll'}\nGetting the line count\n";
        my $lc = `wc -l $data_file`;
        chomp $lc;
        $lc =~ s/\s.*//;
        print "Got the line count\n";
        open IN, "<$data_file" or die "Unable to open $data_file\n";
        my $num_in_chunk = 0;
        my $total = 0;
        my @chunk;
        while(<IN>) {
            chomp;
            my ($gi, $taxon) = split(/\t/, $_);
            $num_in_chunk++;
            push(@chunk, {'gi' => $gi, 'taxon' => $taxon});
            if($num_in_chunk==$self->{'chunk_size'}) {
                $total += $num_in_chunk;
                print join("", ("\r",(sprintf('%.2f', (($total/$lc)*100))),"% complete"));
                $self->insert_chunk($coll, \@chunk);
                @chunk = ();
                $num_in_chunk =0;
            }
        }
        $self->insert_chunk($coll,\@chunk);

        close IN;

        print "Indexing the gis...\n";
        $coll->ensure_index({'gi' => 1},{'safe' => 1, 'background' => 0});
        print "Done indexing the gis\n";
    }
    else {
        print STDERR "Records found in $self->{'gi_coll'}\n";
    }
    return $coll;
}
sub insert_chunk {
    my $self =shift;
    my $coll = shift;
    my $chunk = shift;
    $coll->batch_insert($chunk);
}
sub get_mongodb_connection {
    my($self,$dbname,$host) = @_;
    print STDERR "Connecting to $dbname db on $host\n";
    # First we'll establish our connection to mongodb
    my $conn = MongoDB::Connection->new(host => $host);
    return $conn->get_database($dbname);
}

1;
