
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
use warnings;
use MongoDB;
use Bio::DB::Taxonomy;
use Bio::DB::EUtilities;
use File::Find;

my $version;

my $NODES = '/local/db/repository/ncbi/blast/20120414_001321/taxonomy/taxdump/nodes.dmp';
my $NAMES =  '/local/db/repository/ncbi/blast/20120414_001321/taxonomy/taxdump/names.dmp';
my $GI2TAX = '/local/db/repository/ncbi/blast/20120414_001321/taxonomy/gi_taxid_nucl.dmp';
my $CHUNK_SIZE = 10000;
my $HOST = 'revan.igs.umaryland.edu:10001';
my $DB = 'gi2taxon';
my $COLL = 'gi2taxonnuc';
my $TMP = '/tmp';

sub new {
    my ( $class, $args ) = @_;

    # Get version of MongoDB being used
    # Version is 'v#.#.#' so let's just get the first digit only.
    $version = substr($MongoDB::VERSION, 1, 1);

    my $self = {};
    # If we choose to not provide any taxonomy flatfiles, enable this property
    # This will cause all taxon entries to be added from the Entrez database
    $self->{'no_flatfiles'} =
        $args->{'no_flatfiles'}
        ? $args->{'no_flatfiles'}
        : 0;

    # Nodes/Names/Gi2Tax will only matter if no_flatfiles = 0
    $self->{'nodes'} =
        $args->{'nodes'}
        ? $args->{'nodes'}
	    : $NODES;
    $self->{'names'} =
        $args->{'names'}
        ? $args->{'names'}
	    : $NAMES;

    $self->{'gi2tax'} =
        defined $args->{'gi2tax'}
        ? $args->{'gi2tax'}
        : '';

    $self->{'chunk_size'} =
      $args->{'chunk_size'} ? $args->{'chunk_size'} : $CHUNK_SIZE;
	$self->{'host'} = 'mongodb://';
    $self->{'host'} .=
      $args->{'host'} ? $args->{'host'} : $HOST;
	$self->{'gi_db'} = $args->{'gi_db'} ? $args->{'gi_db'} : $DB;
    $self->{'gi_coll'} =
      $args->{'gi_coll'} ? $args->{'gi_coll'} : $COLL;
    $self->{'taxonomy_dir'} = $args->{'taxonomy_dir'} ? $args->{'taxonomy_dir'} : $TMP;

    $self->{'cache'} = {};

    $self->{'type'} = $args->{'type'} ? $args->{'type'} : 'nucleotide';
    my $gi_tax_file = 'gi_taxid_nucl.dmp';
    if ( $self->{'type'} eq 'protein' ) {
        $gi_tax_file = 'gi_taxid_prot.dmp';
    }

# This option can be used if the user want's to override all the nodes/names params at once
    if ( $args->{'taxon_dir'} ) {
        print STDERR "Here with a taxon directory $args->{'taxon_dir'}\n";

        # Find the nodes, names and nucleotide mapping file
        find(
            sub {
                if ( $File::Find::name =~ /nodes.dmp/ ) {
                    $self->{'nodes'} = $File::Find::name;
                }
                elsif ( $File::Find::name =~ /names.dmp/ ) {
                    $self->{'names'} = $File::Find::name;
                }
                elsif ( $File::Find::name =~ /$gi_tax_file/ ) {
                    $self->{'gi2tax'} = $File::Find::name;
                }
            },
            $args->{'taxon_dir'}
        );
        if ( !$args->{'gi_coll'} ) {
            if ( $args->{'taxon_dir'} =~ /(\d+\_\d+)/ ) {
                my $date = $1;
                if ( $self->{'type'} eq 'nucleotide' ) {
                    $self->{'gi_coll'} = "gi2taxonnuc_$date";
                }
                else {
                    $self->{'gi_coll'} = "gi2taxonprot_$date";
                }
            }
        }
    }

    if (! $self->{'no_flatfiles'}) {
        $self->{'db'} = Bio::DB::Taxonomy->new(
            -source    => 'flatfile',
            -nodesfile => $self->{'nodes'},
            -namesfile => $self->{'names'},
            -directory => $self->{'taxonomy_dir'}
        );
    }

    if ( $args->{verbose} ) {
        if ($self->{'no_flatfiles'}) {
            print STDERR "======== &Gi2Taxon - Not using NCBI taxonomy flatfiles\n";
        } else {
            print STDERR "======== &Gi2Taxon - Using $self->{nodes}\n";
            print STDERR "======== &Gi2Taxon - Using $self->{names}\n";
            print STDERR "======== &Gi2Taxon - Using $self->{gi2tax}\n";
        }
        print STDERR "======== &Gi2Taxon - Using $self->{'taxonomy_dir'}\n";
        print STDERR "======== &Gi2Taxon - Using $self->{'gi_coll'}\n";
        print STDERR "======== &Gi2Taxon - Using $self->{'host'}\n";
    }
    bless $self;
    $self->{'gi2taxon'} = $self->getgi2taxon( $self->{'gi2tax'} );

    return $self;
}

sub getTaxon {
    my ( $self, $acc ) = @_;
    $acc =~ s/^\s+//;
    $acc =~ s/\s/\t/;
    my $gi = $acc;
    if ( $acc =~ /\|/ ) {
        my @fields = split( /\|/, $acc );
        $gi = $fields[1];
    }
    my $taxonid = '';
    my $retval  = {};
    # First check the cache
    if ( $self->{cache}->{$gi} ) {
        return $self->{cache}->{$gi};
    }
    else {
        my $taxon_lookup =
          $self->{'gi2taxon'}->find_one( { 'gi' => "$gi" }, { 'taxon' => 1 } );
        if ($taxon_lookup) {
            $taxonid = $taxon_lookup->{'taxon'};
        }
        else {
            print STDERR
"*** GiTaxon-getTaxon: Unable to find taxon for $gi, Checking NCBI\n";
            my $factory = Bio::DB::EUtilities->new(
                -eutil => 'esummary',
                -email => 'example@foo.bar',
                -db    => $self->{'type'},
                -id    => [$gi]
            );
            while ( my $ds = $factory->next_DocSum ) {
                my @res = $ds->get_contents_by_name('TaxId');
                if (@res) {
                    $taxonid = $res[0];
                }
                if ( !$taxonid ) {
                    print STDERR "Unable to find taxonid at NCBI\n";
                }
                else {
                    my $res = $self->{'gi2taxon'}->update(
                        { 'gi'     => "$gi" },
                        { '$set' => { 'gi'     => "$gi", 'taxon' => $taxonid } },
                        { 'upsert' => 1, 'safe' => 1 }
                    );
                    print STDERR
                      "*** GiTaxon-getTaxon: Added $gi\t$taxonid to the db\n";
                }
            }

        }

        ## NEW VVV 01.08.15 KBS v1.07
        ## I added this so that if the gi isn't in our DB we pull the data from NCBI
        if ( ! $self->{'no_flatfiles'} &&
            ( my $taxon = $self->{'db'}->get_taxon( -taxonid => $taxonid ) ) ) {
            if ( $taxon->isa('Bio::Taxon') ) {
                my $name    = $taxon->scientific_name;
                my $c       = $taxon;
                my @lineage = ($name);
                while ( my $parent = $self->{'db'}->ancestor($c) ) {
                    unshift @lineage, $parent->scientific_name;
                    $c = $parent;
                }
                $retval = {
                    'gi'       => $gi,
                    'acc'      => $acc,
                    'taxon_id' => $taxonid,
                    'name'     => $name,
                    'lineage'  => join( ";", @lineage )
                };
            }
        }
        else {
			# SAdkins 1/13/16 - URL in Entrez module is obsolete and incorrect
			# We really should update BioPerl
            my $db = Bio::DB::Taxonomy->new( -source => 'entrez', -location => 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/');
            my $taxon = $db->get_taxon( -taxonid => $taxonid );
            if (defined($taxon) && $taxon->isa('Bio::Taxon') ) {
                my $name    = $taxon->scientific_name;
                my $c       = $taxon;
                my @lineage = ($name);
                while ( my $parent = $db->ancestor($c) ) {
                    unshift @lineage, $parent->scientific_name;
                    $c = $parent;
                }
                $retval = {
                    'gi'       => $gi,
                    'acc'      => $acc,
                    'taxon_id' => $taxonid,
                    'name'     => $name,
                    'lineage'  => join( ";", @lineage )
                };
            }
            else {
                print STDERR "**GiTaxon unable to find taxon for taxon_id: $taxonid & gi:$gi\n";
            }
        }

        ## NEW ^^^ 01.08.15 KBS v1.07
        $self->{cache}->{$gi} = $retval;
    }
    return $retval;
}

sub getgi2taxon {
    my ( $self, $data_file ) = @_;

    my $mongo =
      $self->get_mongodb_connection( $self->{'gi_db'}, $self->{'host'} );
    my $coll = $mongo->get_collection( $self->{'gi_coll'} );
	# If collection not found in database, update the db using the datafile
    if ($data_file) {
	       if ( !$coll->find_one() ) {
            print
"Found nothing in database $self->{gi_db} collection $self->{gi_coll} on $self->{host}\n";
            print "Getting the line count\n";
            my $lc = `wc -l $data_file`;
            chomp $lc;
            $lc =~ s/\s.*//;
            print "Got the line count\n";
            open IN, "<$data_file" or die "Unable to open $data_file\n";
            my $num_in_chunk = 0;
            my $total        = 0;
            my @chunk;

		      # In the data file, add new data to MongoDB database in chunks
            while (<IN>) {
                chomp;
                my ( $gi, $taxon ) = split( /\t/, $_ );
                $num_in_chunk++;
                push( @chunk, { 'gi' => $gi, 'taxon' => $taxon } );
                if ( $num_in_chunk == $self->{'chunk_size'} ) {
                    $total += $num_in_chunk;
                    print join(
                        "",
                        (
                            "\r", ( sprintf( '%.2f', ( ( $total / $lc ) * 100 ) ) ),
                            "% complete"
                        )
                    );
                    $self->insert_chunk( $coll, \@chunk );
                    @chunk        = ();
                    $num_in_chunk = 0;
                }
            }
            $self->insert_chunk( $coll, \@chunk );

            close IN;
            $coll->ensure_index( { 'gi' => 1 }, { 'safe' => 1 } );
        }
    }
    return $coll;
}

sub insert_chunk {
    my $self  = shift;
    my $coll  = shift;
    my $chunk = shift;

    if ($version < 1) {
        $coll->batch_insert( $chunk, { 'safe' => 1 } );    # Uses older MongoDB
    } else {
        $coll->insert_many( $chunk, { 'safe' => 1 } );      # Uses newer MongoDB
    }
}

sub get_mongodb_connection {
    my ( $self, $dbname, $host ) = @_;

    # First we'll establish our connection to mongodb
    my $conn;
    if ($version < 1) {
        print STDERR "Using version of MongoDB before version 1.0\n";
        require MongoDB::Connection;
        $conn = MongoDB::Connection->new( host => $host );
        $conn->query_timeout(60000000);
    } else {
        print STDERR "Using version of MongoDB after version 1.0\n";
        $conn = MongoDB->connect( $host );
    }
    return $conn->get_database($dbname);
}

1;
