#!/usr/bin/perl

use strict;
use Bio::Annotation::Collection;
use Bio::Annotation::Comment;
use Bio::Seq::RichSeq;
use Bio::SeqFeature::Generic;
use Bio::SeqIO;
use Bio::Species;
#use Bio::Tools::GFF;  ## this is a pain.
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

use Data::Dumper;

my %options = ();
my $results = GetOptions (\%options,
                          'database=s',
                          'user=s',
                          'password=s',
                          'password_file=s',
                          'output_directory=s',
                          'format=s',
                          'locus_db=s',
                          'locus_db_version=s',
                          'translation_table=i',
                          'intergenic_regions=i',
                          'server=s',
                          'database_type=s',
                          'comment=s',
                          'source=s',
                          'organism_id=i',
                          'use_assembly_names',
			              'add_definition',
			              'unique_gene_symbols',
			              'include_pmarks',
                          'cgi_mode=i',
                          'log=s',
                          'help') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## Will either be acquired directly as an option or from the password file
my $password;
my $download_file_name;

## make sure everything passed was peachy
&check_parameters(\%options);

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">$options{log}") || die "can't create log file: $!";
}

my $dsn = '';

if ( $options{database_type} eq 'mysql' ) {
    $dsn = "dbi:mysql:database=$options{database};host=$options{server}";

} elsif ( $options{database_type} eq 'postgresql' ) {
    $dsn = "DBI:Pg:dbname=$options{database};host=$options{server}";
}

_log("attempting to create database connection");
my $dbh = DBI->connect($dsn, $options{user}, $password, {PrintError=>1, RaiseError=>1} );
$dbh->{mysql_auto_reconnect} = 1;

## cvterm lookup
my %cvterm;

## will be reused throughout
my ($qry, @qry_args);

## start by creating the sequences
$cvterm{assembly} = get_cvterm_id('assembly');

## this script was meant for prok genomes, so we should be OK storing thes in memory
##  structure like $h->{$feature_id} = { uniquename => ?, seq_obj =>  )
my $assemblies = {};

$qry = qq{
    SELECT feature_id, uniquename, residues, name
	FROM feature
	WHERE type_id = $cvterm{assembly}
    AND is_obsolete = ?
	AND is_analysis = ?
    };

@qry_args = (0,0);

if ( $options{organism_id} ) {
    $qry .= " AND organism_id = ? ";
    push @qry_args, $options{organism_id};
}

my $assembly_selector = $dbh->prepare($qry);
$assembly_selector->execute( @qry_args );

#get the organism information
my ($genus, $species, $sub_species, $com_name) = &get_organism_information();
my $bio_species = new Bio::Species();
$bio_species->genus( $genus );
$bio_species->species( $species );
$bio_species->sub_species( $sub_species ) if( $sub_species );
$bio_species->common_name( $com_name );
$bio_species->classification( $com_name, $genus." ".$species, $species );

while ( my $row = $assembly_selector->fetchrow_hashref ) {
    _log("INFO: found assembly: $$row{uniquename}");
    $$assemblies{$$row{feature_id}} = { uniquename => $$row{uniquename}, name => $$row{name} };

    my $asmbl_id = ($options{'use_assembly_names'}) ? $$row{'name'} : $$row{'uniquename'};

    my %seq_obj_opts = ( -seq => $$row{residues},
			 -molecule => 'DNA',
			 -display_id => $asmbl_id,
			 -division => 'BCT',
			 -dates => &today()
			 );

    if( $options{'add_definition'} ) {
	my $org_name = $bio_species->binomial('FULL');
	$seq_obj_opts{'description'} = $org_name;
    }

    $$assemblies{ $$row{feature_id} }{seq_obj} = Bio::Seq::RichSeq->new ( %seq_obj_opts );
    $$assemblies{ $$row{feature_id} }{seq_obj}->species( $bio_species );

    if ( $options{comment} ) {
	my $collection = $$assemblies{ $$row{feature_id} }{seq_obj}->annotation();
        my $comment = Bio::Annotation::Comment->new( -text => $options{comment} );
        $collection->add_Annotation( 'comment', $comment );
    }

    my $location = Bio::Location::Simple->new( -start => 1,
					       -end => length( $$row{residues} ),
					       -strand => 1 );

    my $source = Bio::SeqFeature::Generic->new( -primary => 'source',
						-location => $location );

    my $org_name = $bio_species->binomial('FULL');
    $source->add_tag_value('organism', $org_name );
    $source->add_tag_value('molecule', $asmbl_id );
    $source->add_tag_value('strain', $bio_species->sub_species() ) if( $bio_species->sub_species() );
    $source->add_tag_value('mol_type', 'genomic DNA');


    $$assemblies{ $$row{feature_id} }{seq_obj}->add_SeqFeature( $source );

}

$assembly_selector->finish();

## qry for all genes on an assembly
$cvterm{gene} = get_cvterm_id('gene');
$cvterm{transcript} = get_cvterm_id('transcript');
$cvterm{tRNA} = get_cvterm_id('tRNA');
$cvterm{rRNA} = get_cvterm_id('rRNA');
$cvterm{exon} = get_cvterm_id('exon');
$cvterm{frameshift} = get_cvterm_id('frameshift');

## query to fetch all dbxrefs for a given feature id
$qry = qq{
    SELECT db.name, dbx.accession, dbx.version
	FROM db, dbxref dbx, feature f, feature_dbxref fdbx
	WHERE f.feature_id=fdbx.feature_id
	AND dbx.dbxref_id=fdbx.dbxref_id
	AND dbx.db_id=db.db_id
	AND f.feature_id = ?;
};
my $dbxref_selector = $dbh->prepare($qry);

## query to fetch gene product names for a given feature
$cvterm{gene_product_name} = get_cvterm_id('gene_product_name');
$qry = qq{
    SELECT value
	FROM featureprop
	WHERE type_id = $cvterm{gene_product_name}
    AND feature_id = ?
    };
my $product_selector = $dbh->prepare($qry);

## query to fetch gene symbols for a given feature
$cvterm{gene_symbol} = get_cvterm_id('gene_symbol');
$qry = qq{
    SELECT value
	FROM featureprop
	WHERE type_id = $cvterm{gene_symbol}
    AND feature_id = ?
    };
my $gene_symbol_selector = $dbh->prepare($qry);

## query for EC numbers
$qry = qq{
    SELECT t.uniquename, d.accession
	FROM feature t
        JOIN feature_cvterm fc ON t.feature_id = fc.feature_id
        JOIN cvterm c ON fc.cvterm_id = c.cvterm_id
        JOIN cv ON c.cv_id = cv.cv_id
        JOIN cvterm_dbxref cd ON fc.cvterm_id = cd.cvterm_id
        JOIN dbxref d ON cd.dbxref_id = d.dbxref_id
	WHERE t.feature_id = ?
	AND t.type_id =  $cvterm{transcript}
    AND cv.name = 'EC';
};
my $ec_num_selector = $dbh->prepare($qry);

## query to select all primary GO terms for a given feature id
$qry = qq{
    SELECT t.uniquename, db.name as db_name, cvt.name, go.accession
	FROM feature t
	JOIN feature_cvterm fc ON t.feature_id = fc.feature_id
	JOIN cvterm cvt ON fc.cvterm_id = cvt.cvterm_id
	JOIN dbxref go ON go.dbxref_id = cvt.dbxref_id
	JOIN db ON db.db_id = go.db_id
	WHERE db.name IN ('process', 'component', 'function', 'GO' )
	AND t.type_id = $cvterm{transcript}
    AND t.feature_id = ?
    };
my $go_term_selector = $dbh->prepare($qry);

## query to select all (supplemental) GO terms for a given feature id
## NOT USED CURRENTLY
$qry = qq{
    SELECT t.uniquename, db.name as db_name, cv.name, go.accession
	FROM feature t
	JOIN feature_cvterm fc ON t.feature_id = fc.feature_id
	JOIN cvterm cv ON fc.cvterm_id = cv.cvterm_id
	JOIN cvterm_dbxref cd ON cv.cvterm_id = cd.cvterm_id
	JOIN dbxref go ON go.dbxref_id = cd.dbxref_id
	JOIN db ON db.db_id = go.db_id
	WHERE db.name IN ('process', 'component', 'function', 'GO' )
	AND t.type_id = $cvterm{transcript}
    AND t.feature_id = ?
    };
my $supplemental_go_term_selector = $dbh->prepare($qry);

## query to fetch all features of a given type located on a feature
$qry = qq{
    SELECT f1.feature_id, f1.uniquename, f2.feature_id as srcfeature_id, ftl.fmin, ftl.fmax, ftl.strand, ftl.locgroup, ftl.rank
	FROM featureloc ftl, feature f1, feature f2
	WHERE ftl.feature_id=f1.feature_id
	AND ftl.srcfeature_id=f2.feature_id
	AND f1.is_obsolete = ?
	AND f1.is_analysis = ?
	AND f2.feature_id = ?
	AND f1.type_id = ?
    };
my $feature_on_assembly_selector = $dbh->prepare($qry);

## query to fetch all features of a given type located on a feature
$qry = qq{
    SELECT f1.feature_id, f1.uniquename, f2.feature_id as srcfeature_id, ftl.fmin, ftl.fmax, ftl.strand, ftl.locgroup, ftl.rank
	FROM featureloc ftl, feature f1, feature f2
	WHERE ftl.feature_id=f1.feature_id
	AND ftl.srcfeature_id=f2.feature_id
	AND f1.is_obsolete = ?
	AND f2.feature_id = ?
	AND f1.type_id = ?
    };
my $feature_on_assembly_selector_nonfiltered = $dbh->prepare($qry);

## selects CDS related to a passed transcript
$cvterm{derives_from} = get_cvterm_id('derives_from');
$qry = qq{
    SELECT fr.subject_id, f1.uniquename
	FROM feature_relationship fr, feature f1, feature f2
	WHERE fr.subject_id = f1.feature_id
	AND fr.object_id  = f2.feature_id
	AND fr.type_id    = $cvterm{derives_from}
    AND fr.object_id = ?
    };
my $cds_on_transcript_selector = $dbh->prepare($qry);

## selects gene related to a passed transcript
$qry = qq{
    SELECT f1.feature_id, f1.uniquename
	FROM feature_relationship fr, feature f1, feature f2
	WHERE fr.object_id = f1.feature_id
	AND fr.subject_id  = f2.feature_id
	AND fr.type_id    = $cvterm{derives_from}
    AND f2.feature_id = ?
    };
my $gene_on_transcript_selector = $dbh->prepare($qry);

## selects exon(s) related to a passed transcript and their coordinates on an assembly
$cvterm{part_of} = get_cvterm_id('part_of');
$qry = qq{
    SELECT fr.subject_id, exon.uniquename, fl.fmin, fl.fmax, fl.strand
	FROM feature_relationship fr
        JOIN feature transcript ON transcript.feature_id = fr.object_id
        JOIN feature exon ON fr.subject_id=exon.feature_id
        JOIN featureloc fl ON exon.feature_id=fl.feature_id
        JOIN feature assembly ON fl.srcfeature_id=assembly.feature_id
	WHERE exon.type_id = $cvterm{exon}
    AND assembly.feature_id = ?
	AND transcript.uniquename = ?
    };
my $exon_on_transcript_selector = $dbh->prepare($qry);

## selects all pmarks on given assembly.
$qry = qq{
    SELECT fl.fmin, fl.fmax
	FROM featureloc fl, cvterm c, feature assembly, feature pmark
	WHERE assembly.feature_id = ?
	AND fl.srcfeature_id = assembly.feature_id
	AND fl.feature_id = pmark.feature_id
	AND pmark.type_id = c.cvterm_id
	AND c.name = 'pmark_spacer'
};
my $pmarks_on_asmbl = $dbh->prepare($qry);


## track parents for GFF export
## structure like $h{$cds_uniquename} = $transcript_uniquename;
my %cds_parents;
my %transcript_parents;
my $frameshifts = {};
my $point_mutations = {};


for my $feature_id ( keys %$assemblies ) {

    ## This will hold all gene symbols to detect repeats per assembly.
    my %gene_symbols;

    ## until I can find a feature sorting method within the bioperl objects I'm using
    #   this to store the feature references so they can be added to the SeqIO in order
    my @assembly_feats;

    ## get all the validated frameshifts on this molecule
    $feature_on_assembly_selector->execute( 0, 0, $feature_id, $cvterm{frameshift} );
    while ( my $row = $feature_on_assembly_selector->fetchrow_hashref ) {
        $$frameshifts{ $feature_id }{ $$row{uniquename} } = {
            fmin => $$row{fmin},
            fmax => $$row{fmax},
            strand => $$row{strand},
        };
    }

    ## get the transcript features.  these correspond to the CDS entries for proks in the gbk file
    ##  using 'CDS' would include the tRNAs, which are queried separately
    _log("INFO: getting transcripts by executing: feature_on_assembly_selector->execute( 0, 0, $feature_id, $cvterm{transcript} )");
    $feature_on_assembly_selector->execute( 0, 0, $feature_id, $cvterm{transcript} );

    while (my $row = $feature_on_assembly_selector->fetchrow_hashref ) {

        ## make a gene entry for this too
	my $location = Bio::Location::Simple->new( -start => $$row{fmin} + 1,
						   -end => $$row{fmax},
						   -strand => $$row{strand} );

        my $gene = Bio::SeqFeature::Generic->new(
						 -seq_id => $$row{uniquename},
						 -primary => 'gene',
						 -location => $location
						 );
        my $transcript_uniquename = $$row{uniquename};

	## add any CDS related to this gene (fails if more than one is found)
	my $cds;
	$cds_on_transcript_selector->execute( $$row{feature_id} );
	while ( my $cds_row = $cds_on_transcript_selector->fetchrow_hashref ) {
	    if ( defined $cds ) {
		die "found more than one CDS on transcript feature ID $$row{feature_id}\n";
	    }

	    $cds_parents{$$cds_row{uniquename}} = $$row{uniquename};

	    my $location = Bio::Location::Simple->new( -start => $$row{fmin} + 1,
						       -end => $$row{fmax},
						       -strand => $$row{strand} );
	    $cds = Bio::SeqFeature::Generic->new(
						 -seq_id => $$cds_row{uniquename},
						 -primary => 'CDS',
						 -location => $location,
						 -tag => {
						     codon_start => 1,
						     transl_table => $options{translation_table},
						 },
						 );
	}

        ## examine any dbxrefs this has
        my $gene_feature_id;
        $gene_on_transcript_selector->execute( $$row{feature_id} );
        while( my $gene_row = $gene_on_transcript_selector->fetchrow_hashref ) {
            if( defined( $gene_feature_id ) ) {
                die("found more than on gene on transcript feature ID $$row{feature_id}\n");
            }
            $gene_feature_id = $$gene_row{feature_id};
            $transcript_parents{$transcript_uniquename} = $$gene_row{uniquename};
        }
        $dbxref_selector->execute( $gene_feature_id );

        while ( my $dbxref_row = $dbxref_selector->fetchrow_hashref ) {
            ## don't do anything if the accession value is empty
            next unless $$dbxref_row{accession};

            ## if this row source matches the locus_db and version specified by the user, add
            ##  it as a locus, else just a dbxref.
            if ( defined $options{locus_db} &&
                 defined $options{locus_db_version} &&
                 $options{locus_db} eq $$dbxref_row{name} &&
                 $options{locus_db_version} eq $$dbxref_row{version} ) {
                $cds->add_tag_value( 'locus_tag', "$$dbxref_row{accession}" );
                $gene->add_tag_value( 'locus_tag', "$$dbxref_row{accession}" );
            } else {
                $cds->add_tag_value( 'db_xref', $$dbxref_row{name}."_".$$dbxref_row{version}.":".$$dbxref_row{accession} );
            }
        }

        ## look for a gene product name
        $product_selector->execute( $$row{feature_id} );
        while ( my $product_row = $product_selector->fetchrow_hashref ) {
            $cds->add_tag_value( 'product', $$product_row{value} );
            #$gene->add_tag_value( 'product', $$product_row{value} );	# Gene tags should not have products

            if ( $$product_row{value} =~ /point mutation/i ) {
                $$point_mutations{ $gene->seq_id }++;
            }

            last;
        }

        ## gene symbol
        $gene_symbol_selector->execute( $$row{feature_id} );
        while ( my $genesym_row = $gene_symbol_selector->fetchrow_hashref ) {
	    my $gs = $genesym_row->{'value'};
	    if( $options{'unique_gene_symbols'} && exists( $gene_symbols{$gs} ) ) {
		$gs .= "_".$gene_symbols{$gs}++;
	    } else {
		$gene_symbols{$gs} = 1;
	    }
            $cds->add_tag_value( 'gene', $gs );
            $gene->add_tag_value( 'gene', $gs );
            last;
        }

        ## EC number
        $ec_num_selector->execute( $$row{feature_id} );
        while ( my $ec_num_row = $ec_num_selector->fetchrow_hashref ) {
            $cds->add_tag_value( 'EC_number', $$ec_num_row{accession} );
            #last;
        }

        ## Go Terms
        $go_term_selector->execute( $$row{feature_id} );
        my @go_terms;
        while ( my $go_term_row = $go_term_selector->fetchrow_hashref ) {
            #we don't want to include the most generic of go terms
            #(biological_process, molecular_function, cellular_component)
            next if( $$go_term_row{name} eq 'biological_process' ||
                     $$go_term_row{name} eq 'molecular_function' ||
                     $$go_term_row{name} eq 'cellular_component' );

            push( @go_terms, [$$go_term_row{db_name}, $$go_term_row{accession}, $$go_term_row{name}] );

        }
        my @values;
        map {
            my $namespace = 'GO';
            $namespace .= "_".$_->[0] unless( $_->[0] eq 'GO' );
            push( @values, $namespace.": $_->[1] - $_->[2]" );
        } @go_terms;
        $cds->add_tag_value( 'note', join( "; ", @values ) ) unless (scalar @values == 0);

        push @assembly_feats, $gene;
        push @assembly_feats, $cds;

    }

    ## now get tRNAs
    _log("INFO: getting tRNAs by executing: feature_on_assembly_selector->execute( 0, 0, $feature_id, $cvterm{tRNA} )");
    $feature_on_assembly_selector->execute( 0, 0, $feature_id, $cvterm{tRNA} );

    while (my $row = $feature_on_assembly_selector->fetchrow_hashref ) {

	my $location = Bio::Location::Simple->new( -start => $$row{fmin} + 1,
						   -end => $$row{fmax},
						   -strand => $$row{strand} );

        my $cds = Bio::SeqFeature::Generic->new(
						-seq_id => $$row{uniquename} || 'UNKN0WN',
						-primary => 'tRNA',
						-location => $location,
						);

	# I'll makea  new one just so I know they won't share the same
	$location = Bio::Location::Simple->new( -start => $$row{fmin} + 1,
						-end => $$row{fmax},
						-strand => $$row{strand} );

        my $gene = Bio::SeqFeature::Generic->new(
						 -seq_id => $$row{uniquename} || 'UNKN0WN',
						 -primary => 'gene',
						 -location => $location
						 );

        ## look for locus ids (dbxref)
        my $gene_feature_id;
        $gene_on_transcript_selector->execute( $$row{feature_id} );
        while( my $gene_row = $gene_on_transcript_selector->fetchrow_hashref ) {
            if( defined( $gene_feature_id ) ) {
                die("found more than on gene on tRNA feature ID $$row{feature_id}\n");
            }
            $gene_feature_id = $$gene_row{feature_id};
        }
        $dbxref_selector->execute( $gene_feature_id );

        while ( my $dbxref_row = $dbxref_selector->fetchrow_hashref ) {
            ## don't do anything if the accession value is empty
            next unless $$dbxref_row{accession};

            ## if this row source matches the locus_db and version specified by the user, add
            ##  it as a locus, else just a dbxref.
            if ( defined $options{locus_db} &&
                 defined $options{locus_db_version} &&
                 $options{locus_db} eq $$dbxref_row{name} &&
                 $options{locus_db_version} eq $$dbxref_row{version} ) {
                $cds->add_tag_value( 'locus_tag', "$$dbxref_row{accession}" );
                $gene->add_tag_value( 'locus_tag', "$$dbxref_row{accession}" );
            } else {
                $cds->add_tag_value( 'db_xref', $$dbxref_row{name}."_".$$dbxref_row{version}.":".$$dbxref_row{accession} );
            }

        }

        ## look for a gene product name
        $product_selector->execute( $$row{feature_id} );
        while ( my $product_row = $product_selector->fetchrow_hashref ) {
            $cds->add_tag_value( 'product', $$product_row{value} );
            $gene->add_tag_value( 'product', $$product_row{value} );

            if ( $$product_row{value} =~ /point mutation/i ) {
                $$point_mutations{ $gene->seq_id }++;
            }

            last;
        }

        push @assembly_feats, $cds;
        push @assembly_feats, $gene;
    }

    ## now get all rRNAs
    _log("INFO: getting rRNAs by executing: feature_on_assembly_selector_nonfiltered->execute( 0, $feature_id, $cvterm{rRNA} )");
    $feature_on_assembly_selector_nonfiltered->execute( 0, $feature_id, $cvterm{rRNA} );

    while( my $row = $feature_on_assembly_selector_nonfiltered->fetchrow_hashref ) {
        my $cds = Bio::SeqFeature::Generic->new(
						-seq_id => $$row{uniquename} || 'UNKN0WN',
						-primary => 'rRNA',
						-start => $$row{fmin} + 1,
						-end => $$row{fmax},
						-strand => $$row{strand},
						);

        my $gene = Bio::SeqFeature::Generic->new(
						 -seq_id => $$row{uniquename} || 'UNKN0WN',
						 -primary => 'gene',
						 -start => $$row{fmin} + 1,
						 -end => $$row{fmax},
						 -strand => $$row{strand},
						 );

        ## look for locus ids (dbxref)
        my $gene_feature_id;
        $gene_on_transcript_selector->execute( $$row{feature_id} );
        while( my $gene_row = $gene_on_transcript_selector->fetchrow_hashref ) {
            if( defined( $gene_feature_id ) ) {
                die("found more than on gene on rRNA feature ID $$row{feature_id}\n");
            }
            $gene_feature_id = $$gene_row{feature_id};
        }
        $dbxref_selector->execute( $gene_feature_id );

        while ( my $dbxref_row = $dbxref_selector->fetchrow_hashref ) {
            ## don't do anything if the accession value is empty
            next unless $$dbxref_row{accession};

            ## if this row source matches the locus_db specified by the user, add
            ##  it as a locus, else just a dbxref.
            if ( defined $options{locus_db} && $options{locus_db} eq $$dbxref_row{name} ) {
                $cds->add_tag_value( 'locus_tag', "$$dbxref_row{accession}" );
                $gene->add_tag_value( 'locus_tag', "$$dbxref_row{accession}" );
            } else {
                $cds->add_tag_value( 'db_xref', "$$dbxref_row{name}:$$dbxref_row{accession}" );
            }
        }

        ## look for a gene product name
        $product_selector->execute( $$row{feature_id} );
        while ( my $product_row = $product_selector->fetchrow_hashref ) {
            $cds->add_tag_value( 'product', $$product_row{value} );
            last;
        }

        push @assembly_feats, $gene;
        push @assembly_feats, $cds;

    }

    ## if they asked for pmarks
    if( $options{'include_pmarks'} ) {
	$pmarks_on_asmbl->execute( $feature_id );
	while( my $row = $pmarks_on_asmbl->fetchrow_hashref ) {
	    my $misc = Bio::SeqFeature::Generic->new(
						     -primary => 'misc_feature',
						     -start => $$row{fmin} + 1,
						     -end => $$row{fmax},
						     -strand => 1,
						     );
	    $misc->add_tag_value('note', 'pmark_spacer');
	    push( @assembly_feats, $misc );

	}
    }

	add_intergenic_regions(\@assembly_feats, $feature_id) if ($options{'intergenic_regions'});

    ## now add all the features
    for my $feat ( sort { $a->start <=> $b->start } @assembly_feats ) {
	    $$assemblies{$feature_id}{seq_obj}->add_SeqFeature( $feat );
	    my $whole_seq = $$assemblies{$feature_id}{seq_obj};
        my $seq_len = length($whole_seq->seq );
	    if ( $feat->primary_tag eq 'CDS' ) {
	    # Assign start and end coordinates, but adjust if partial gene is present.  Add translation tag to the sequence
            eval {
                my $temp_start = ($feat->start < 1) ? 1 : $feat->start;	#for partial genes extending beyond coord 1
                my $temp_end = ($feat->end > $seq_len) ? $seq_len : $feat->end;	# for partial genes extending beyond the end
                my $subseq = $whole_seq->trunc($temp_start, $temp_end);	#get specific sequence from this region
                $subseq = $subseq->revcom if ($feat->strand eq "-1");	#get compliments of sequences that need them

                #print STDERR "\nINFO: seq " . $feat->seq_id . "  has a length of: " . length($subseq->seq) . "\n";
                #print $temp_start, " - " , $temp_end, "\n";
                if ( length $subseq->seq  ) {# should always have a length
                   # print STDERR "INFO: exporting sequence " . $feat->seq_id . "\n";
                    if ( $options{format} eq 'gbk' ) {
                	$feat->add_tag_value( 'translation', $subseq->translate(-complete => 1)->seq );
                    } else {
                    	$feat->add_tag_value( 'translation', $subseq->translate()->seq );
                    }
                }
            };

            if ($@) {
                print STDERR "WARN: failure to export sequence " . $feat->seq_id . " for feature at coordinates (" . $feat->start . '/' . $feat->end . ")\n";
            }
	    }
    }
}


## don't need the database anymore
$exon_on_transcript_selector->finish();
$cds_on_transcript_selector->finish();
$ec_num_selector->finish();
$gene_symbol_selector->finish();
$product_selector->finish();
$dbxref_selector->finish();
$feature_on_assembly_selector->finish();
$feature_on_assembly_selector_nonfiltered->finish();
$go_term_selector->finish();
$supplemental_go_term_selector->finish();
$pmarks_on_asmbl->finish();

$dbh->disconnect();


## if this is a CGI, we need to print a header:
if ( $options{cgi_mode} == 1 ) {

    ## quicker way than this to get YYYYMMDD?
    my ($year, $mon, $day) = (localtime time)[5,4,3];
    my $datestamp = sprintf("%d%02d%02d", $year + 1900, ++$mon, $day);

    if ( $options{format} eq 'gbk' ) {
        $download_file_name = "$options{database}.annotation.$datestamp.gbk";
    } elsif ( $options{format} eq 'gff' ) {
        $download_file_name = "$options{database}.annotation.$datestamp.gff3";
    } elsif ( $options{format} eq 'tbl' ) {
        $download_file_name = "$options{database}.annotation.$datestamp.tbl";
    }

    print "Content-Type:application/x-download\n";
    print "Content-Disposition:attachment;filename=$download_file_name\n\n";
}

## if writing tbl format, all output goes into a single file
my $tbl_fh;
if ( $options{format} eq 'tbl' ) {

    if ( $options{cgi_mode} == 1 ) {
	    (open $tbl_fh, ">". $options{output_directory}."/".$download_file_name) || die "can't create TBL output file: $!";
    } else {
        _log("INFO: writing $options{output_directory}/$options{database}.tbl");
        open($tbl_fh, ">$options{output_directory}/$options{database}.tbl") || die "can't create TBL output file: $!";
    }
}

my $gbk;
if ( $options{cgi_mode} == 1 ) {
    $gbk = Bio::SeqIO->new( -format => 'genbank',
                            -file => ">" . $options{output_directory}."/".$download_file_name );
}

## write out the annotations in either GBK, GFF or tbl format
for my $feature_id ( keys %$assemblies ) {
    my $assembly_length = length( $$assemblies{$feature_id}{seq_obj}->seq );

    if ( $options{format} eq 'gbk' ) {
        if (! $options{cgi_mode}) {
            _log("INFO: writing $options{output_directory}/$$assemblies{$feature_id}{uniquename}.gbk");
            my $file;
            if( $options{'use_assembly_names'} ) {
                $file = ">$options{output_directory}/$$assemblies{$feature_id}{name}.gbk";
            } else {
                $file = ">$options{output_directory}/$$assemblies{$feature_id}{uniquename}.gbk";
            }
            $gbk = Bio::SeqIO->new( -format => 'genbank', -file => $file );
        }

        foreach my $feat ( $$assemblies{$feature_id}{seq_obj}->get_SeqFeatures ) {
            # is this a partial?
            $feat->location->start("<1") if ($feat->start <= 0);
            $feat->location->end("$assembly_length>") if ($feat->end > $assembly_length);
        }

        $gbk->write_seq( $$assemblies{$feature_id}{seq_obj} );

    } elsif ( $options{format} eq 'tbl' ) {

        if( $options{'use_assembly_names'} ) {
            print $tbl_fh ">Feature $$assemblies{$feature_id}{name}\n";
        } else {
            print $tbl_fh ">Feature $$assemblies{$feature_id}{uniquename}\n";
        }


        foreach my $feat ( $$assemblies{$feature_id}{seq_obj}->get_SeqFeatures ) {

            my ($is_5_partial, $is_3_partial) = (0,0);

            ## strand: 1 = plus, strand 2 = minus
            if ($feat->strand == 1 ) {
                ## is this a partial?
                if ( $feat->start <= 0 ) {
                    print $tbl_fh "<1";
                    $is_5_partial = 1;
                } else {
                    print $tbl_fh $feat->start;
                }

                print $tbl_fh "\t";

                if ( $feat->end > $assembly_length ) {
                    print $tbl_fh ">$assembly_length";
                    $is_3_partial = 1;
                } else {
                    print $tbl_fh $feat->end;
                }

                print $tbl_fh "\t", $feat->primary_tag, "\n";

            } else {

                if ( $feat->end > $assembly_length ) {
                    print $tbl_fh "<$assembly_length";
                    $is_5_partial = 1;
                } else  {
                    print $tbl_fh $feat->end;
                }

                print $tbl_fh "\t";

                if ( $feat->start <= 0 ) {
                    print $tbl_fh ">1";
                    $is_3_partial = 1;
                } else {
                    print $tbl_fh $feat->start;
                }

                print $tbl_fh "\t", $feat->primary_tag, "\n";
            }

            my $locus_tag = '';

            if ( $feat->primary_tag eq 'gene' ) {

                if ( $feat->has_tag('gene') ) {
                    print $tbl_fh "\t\t\tgene\t", ($feat->get_tag_values('gene'))[0], "\n";
                }

                if ( $feat->has_tag('locus_tag') ) {
                    $locus_tag = ($feat->get_tag_values('locus_tag'))[0];
                    print $tbl_fh "\t\t\tlocus_tag\t$locus_tag\n";
                }

                if ( has_frameshift($feature_id, $feat) ||
                     exists $$point_mutations{ $feat->seq_id } ) {

                    print $tbl_fh "\t\t\tpseudo\n";
                }

            } elsif ( $feat->primary_tag eq 'CDS' ) {

                my $min_coord = $feat->start < $feat->end ? $feat->start : $feat->end;

                if ( $feat->has_tag('locus_tag') ) {
                    $locus_tag = ($feat->get_tag_values('locus_tag'))[0];
                }

                if ( $locus_tag eq '' ) {

                }

                ## we need to specify a codon_start if it's a partial gene
                if ( $is_5_partial ) {
                    if ( $feat->strand != 1 ) {

                        my $offset = $assembly_length - $feat->end;
                        #print STDERR "INFO: exporting a different codon_start for reverse locus ($locus_tag) ($offset), offset calculated by: $assembly_length - (" . $feat->end . ")\n";

                        if ( $offset == -1 ) {
                            print $tbl_fh "\t\t\tcodon_start\t3\n";  ## tested on gss1617_contigs.SD1617_0004
                        } elsif ( $offset == -2 ) {
                            print $tbl_fh "\t\t\tcodon_start\t2\n";
                        } elsif ( $offset == -3 ) {
                            #print $tbl_fh "\t\t\tcodon_start\t3\n";
                        }

                    } else {

                        print STDERR "INFO: exporting a different codon_start for forward locus ($locus_tag) ($min_coord)\n";

                        if ( $min_coord == 0 ) {
                            print $tbl_fh "\t\t\tcodon_start\t3\n";
                        } elsif ( $min_coord == -1 ) {
                            print $tbl_fh "\t\t\tcodon_start\t2\n";
                        } elsif ( $min_coord == -2 ) {
                            print $tbl_fh "\t\t\tcodon_start\t1\n";
                        }
                    }
                }

                if ( $feat->has_tag('product') ) {
                    my $product_name = ($feat->get_tag_values('product'))[0];

                    if ( has_frameshift($feature_id, $feat) && $product_name !~ /frameshift/ ) {
                        print $tbl_fh "\t\t\tproduct\t$product_name, frameshift\n";
                    } else {
                        print $tbl_fh "\t\t\tproduct\t$product_name\n";
                    }
                }

                if ( $feat->has_tag('EC_number') ) {
                    for my $ec_num ( $feat->get_tag_values('EC_number') ) {
                        print $tbl_fh "\t\t\tEC_number\t$ec_num\n";
                    }
                }

                if ( $feat->has_tag('locus_tag') ) {
                    print $tbl_fh "\t\t\tprotein_id\tgnl|IGS|", ($feat->get_tag_values('locus_tag'))[0], "\n";
                }

            } elsif ( $feat->primary_tag eq 'rRNA' || $feat->primary_tag eq 'tRNA' ) {

                if ( $feat->has_tag('product') ) {
                    my $product_name = ($feat->get_tag_values('product'))[0];

                    if ( has_frameshift($feature_id, $feat) && $product_name !~ /frameshift/ ) {
                        print $tbl_fh "\t\t\tproduct\t$product_name, frameshift\n";
                    } else {
                        print $tbl_fh "\t\t\tproduct\t$product_name\n";
                    }
                }
            }
        }

    } elsif ( $options{format} eq 'gff' ) {

        my $gff_fh;

        if ( $options{cgi_mode} == 1 ) {
            open($gff_fh, ">".$options{output_directory}."/".$download_file_name) || die "can't create GFF output file: $!";
        } else {
            _log("INFO: writing $options{output_directory}/$$assemblies{$feature_id}{uniquename}.gff3");
            my $file;
            if( $options{'use_assembly_names'} ) {
                $file = "$options{output_directory}/$$assemblies{$feature_id}{name}.gff3";
            } else {
                $file = "$options{output_directory}/$$assemblies{$feature_id}{uniquename}.gff3";
            }
            open($gff_fh, "> $file") || die "can't create GFF output file [$file]: $!";
        }

		my %gene_locus;

		# First let's map transcript to locus tags
        foreach my $feat ( $$assemblies{$feature_id}{seq_obj}->get_SeqFeatures ) {
        	if ( $feat->primary_tag eq 'gene' ) {
                if ( $feat->has_tag('locus_tag') ) {
                    $gene_locus{$feat->seq_id} = ($feat->get_tag_values('locus_tag'))[0];
                } else {
                	print $feat->seq_id . " does not have a locus tag\n" if (defined $options{locus_db});
                }
            }
        }

        foreach my $feat ( $$assemblies{$feature_id}{seq_obj}->get_SeqFeatures ) {

            ## writing my own GFF output.  A Bio::SeqIO layer would be cool eventually.
            my $asmbl_id = ($options{'use_assembly_names'}) ? $$assemblies{$feature_id}{name} : $$assemblies{$feature_id}{uniquename};
            my @columns = ( $asmbl_id,
                            $options{source},
                            $feat->primary_tag,
                            $feat->start,
                            $feat->end,
                            '.'
			    );

            push @columns, ($feat->strand == 1) ? '+' : '-';
            push @columns, ($feat->primary_tag eq 'CDS') ? 0 : '.';

            if ( $feat->primary_tag eq 'gene' ) {

                ## print the first 8 columns
                print $gff_fh join("\t", @columns);
                print $gff_fh "\t";

                if (defined $options{locus_db}){
                    die ("No locus tag found for gene " . $transcript_parents{$feat->seq_id}) if (! defined $gene_locus{$feat->seq_id});
                    print_gff_col9_attribute($gff_fh, 'ID', $gene_locus{$feat->seq_id} );
                } else {
                	print_gff_col9_attribute($gff_fh, 'ID', $transcript_parents{$feat->seq_id} );
		}

                ## now handle any column 9 attributes
                if ( $feat->has_tag('product') ) {
                    print_gff_col9_attribute($gff_fh, 'description', ($feat->get_tag_values('product'))[0] );
                }

                if ( $feat->has_tag('gene') ) {
                    print_gff_col9_attribute($gff_fh, 'gene_symbol', ($feat->get_tag_values('gene'))[0] );
                }

                if ( $feat->has_tag('locus_tag') ) {
                    print_gff_col9_attribute($gff_fh, 'locus_tag', ($feat->get_tag_values('locus_tag'))[0] );
                }

                print $gff_fh "\n";

		## GFF requires exons too.  write an mRNA feature, then the CDS ones
            } elsif ( $feat->primary_tag eq 'CDS' ) {

                my $gene_id = $transcript_parents{$cds_parents{ $feat->seq_id}};
                my $transcript_id = $cds_parents{ $feat->seq_id};

                ## the mRNA feature matches the gene one
                $columns[2] = 'mRNA';
                $columns[7] = '.';
                print $gff_fh join("\t", @columns);
                print $gff_fh "\t";
                if (defined $options{locus_db}){
                    die ("No locus tag found for mRNA" . $transcript_id) if (! defined $gene_locus{$transcript_id});
                    print_gff_col9_attribute($gff_fh, 'ID', $gene_locus{$transcript_id} );
                } else {
                	print_gff_col9_attribute($gff_fh, 'ID', $transcript_id );
                }
                print_gff_col9_attribute($gff_fh, 'Parent', $gene_id );
                print $gff_fh "\n";

                $exon_on_transcript_selector->execute( $feature_id, $transcript_id );

                ## we'll write one CDS per exon, as described in the GFF3 specs
                while ( my $row = $exon_on_transcript_selector->fetchrow_hashref ) {
                    print $gff_fh "$$assemblies{$feature_id}{uniquename}\t" .
			"$options{source}\tCDS\t" .
			($$row{fmin} + 1) . "\t$$row{fmax}\t.\t";

                    print $gff_fh ($$row{strand} == 1) ? '+' : '-';
                    print $gff_fh "\t0\t";

                    # Changed to print locus tag in place of CDS_ID 10/22/14 -- Shaun Adkins  (JIRA AE-621)
                    if (defined $options{locus_db}){
                	die ("No locus tag found for CDS parent " . $transcript_id) if (! defined $gene_locus{$transcript_id});
                    	print_gff_col9_attribute($gff_fh, 'ID', $gene_locus{$transcript_id} );
                    } else {
                    	print_gff_col9_attribute($gff_fh, 'ID', $feat->seq_id);
                    }
                    print_gff_col9_attribute($gff_fh, 'Parent', $transcript_id );
                    print $gff_fh "\n";

                    # Write out exon feature now
                    print $gff_fh "$$assemblies{$feature_id}{uniquename}\t" .
			"$options{source}\texon\t" .
			($$row{fmin} + 1) . "\t$$row{fmax}\t.\t";

                    print $gff_fh ($$row{strand} == 1) ? '+' : '-';
                    print $gff_fh "\t.\t";
                    if (defined $options{locus_db}){
                	die ("No locus tag found for exon parent " . $transcript_id) if (! defined $gene_locus{$transcript_id});
                    	print_gff_col9_attribute($gff_fh, 'ID', $gene_locus{$transcript_id} );
                    } else {
                    	print_gff_col9_attribute($gff_fh, 'ID', $$row{uniquename} );
                    }
                    print_gff_col9_attribute($gff_fh, 'Parent', $transcript_id );
                    print $gff_fh "\n";

                }

            }
        }
    }
}


exit(0);

sub has_frameshift {
    my ($assembly_id, $feat) = @_;

    ## check the coordinates of this feature against known frameshift coordinates
    if ( exists $$frameshifts{$assembly_id} ) {

        for my $frameshift_id ( keys %{ $$frameshifts{$assembly_id} } ) {
            my $frameshift = $$frameshifts{$assembly_id}{$frameshift_id};

            if ( $$frameshift{fmin} >= $feat->start() &&
                 $$frameshift{fmin} <= $feat->end() ) {

                return 1;
            }
        }
    }

    return 0;
}

sub add_intergenic_regions {
	my $a_feats = shift;
	my $feature_id = shift;

	# Create array to store intergenic region features and initialize the source area
	my @igr_feats;
	my $igr_coord_start = 1;
	my $igr_coord_end;
	my $igr_locus_start = "Source0";
	my $igr_locus_end;
    ## Sort and iterate through all 'gene' features
    for my $feat ( sort { $a->start <=> $b->start } @{$a_feats} ) {
        if ($feat->primary_tag eq 'gene'){
            # Creating the intergenic region before the current gene feature
        	$igr_coord_end = $feat->start - 1;
            if ($feat->has_tag('locus_tag')){
        	    $igr_locus_end = ($feat->get_tag_values('locus_tag'))[0];
	        } else {
                $igr_locus_end = "coord-" . $igr_coord_end;
            }

            my $igr = Bio::SeqFeature::Generic->new(
						     -primary => 'intergenic',
						     -start => $igr_coord_start,
						     -end => $igr_coord_end,
						     -strand => 1,
						     );
	    	$igr->add_tag_value( 'locus_tag', "ig-$igr_locus_start-$igr_locus_end" );

            # If we have genes that overlap then do not write an intergenic region
            if ($feat->start <= $igr_coord_start){
                $igr_coord_start = $feat->end + 1;  # Want to update start coord AFTER comparison check
                $igr_locus_start = $igr_locus_end;
                next;
            }

        	$igr_coord_start = $feat->end + 1;
        	$igr_locus_start = $igr_locus_end;

            push @igr_feats, $igr;
        }

            if ($@) {
                print STDERR "WARN: failure to export sequence " . $feat->seq_id . " for feature at coordinates (" . $feat->start . '/' . $feat->end . ")\n";
            }
	}

	# Add final intergenic region at the end
	my $whole_seq = $$assemblies{$feature_id}{seq_obj};
    my $seq_len = length($whole_seq->seq );
	$igr_coord_end = $seq_len;
	$igr_locus_end = "end";
	my $igr = Bio::SeqFeature::Generic->new(
						     -primary => 'intergenic',
						     -start => $igr_coord_start,
						     -end => $igr_coord_end,
						     -strand => 1,
						     );
	$igr->add_tag_value( 'locus_tag', "ig-$igr_locus_start-$igr_locus_end" );

    push @igr_feats, $igr;


	push @{$a_feats}, $_ foreach (@igr_feats);
	return;
}

sub get_organism_information {
    my $query = "SELECT genus, species, common_name FROM organism";
    if( $options{'organism_id'} ) {
        $query .= " WHERE organism_id = ".$options{'organism_id'};
    }

    my $sth = $dbh->prepare($query);
    $sth->execute();

    my ($genus, $species, $common_name);
    while( my $row = $sth->fetchrow_arrayref ) {
        next if( $row->[0] =~ /not known/ );
        if( $genus || $species || $common_name ) {
            die("Found multiple organism in database. Please specify an organism_id");
        }
        ($genus, $species, $common_name) = @{$row};
    }
    $sth->finish();

    my $sub_species;
    my @parts = split(/\s+/, $species);
    if( @parts > 1 ) {
        $species = shift @parts;
        $sub_species = join(" ", @parts);
    }

    return ($genus, $species, $sub_species, $common_name);
}

sub get_cvterm_id {
    my $name = shift;
    my $qry = qq{
        SELECT cvterm_id
	    FROM cvterm
	    WHERE name = ?
	};

    my $cvterm_selector = $dbh->prepare($qry);
    $cvterm_selector->execute($name);

    my $cvterm_id = ( $cvterm_selector->fetchrow_array )[0];
    $cvterm_selector->finish();

    if ( $cvterm_id ) {
        _log("INFO: got cvterm_id $cvterm_id for name $name");
        return $cvterm_id;
    } else {
        _log("ERROR: failed to retrieve cvterm_id for name $name");
        die "ERROR: failed to retrieve cvterm_id for name $name\n";
    }
}


## gets the numerical db id for a passed db.name.  returns undef if not found.
sub get_db_id {
    my $db_name = shift;

    my $db_id;

    ## name is a unique key, so we don't have worry about multiples
    my $qry = qq{
        SELECT db_id
	    FROM db
	    WHERE name = ?
	};
    my $db_selector = $dbh->prepare($qry);
    $db_selector->execute( $db_name );

    while ( my $db = $db_selector->fetchrow_hashref ) {
        $db_id = $db->{db_id};
        last;
    }

    $db_selector->finish();

    return $db_id;
}

sub print_gff_col9_attribute {
    my ($fh, $key, $val) = @_;

    ## we could use general encoding methods here, such as:
    #   $val =~ s/([^A-Za-z0-9])/sprintf("%%%02X", ord($1))/seg;
    #  but the GFF3 spec only calls for encoding of:
    #
    #   , -> %2C
    #   = -> %3D
    #   ; -> %3B

    $val =~ s/\,/\%2C/g;
    $val =~ s/\=/\%3D/g;
    $val =~ s/\;/\%3B/g;

    print $fh "$key=\"$val\"\;";
}

sub today {
    my @time = localtime();
    my ($day, $mon, $year) = @time[3..5];
    $day = sprintf( "%02d", $day );
    my @m = qw(JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC);
    $mon = $m[$mon];
    $year += 1900;
    return "$day-$mon-$year";
}

sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
}

sub check_parameters {
    my $options = shift;

    ## make sure required arguments were passed
    my @required = qw( database user output_directory format );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }

    if (!defined($$options{password}) && !defined ($$options{password_file}) ){
    	print STDERR ("Neither a password nor a path to a password file were defined.  Please provide one or the other\n");
    }

    if (defined $$options{password}) {
    	$password = $$options{password};
    }

    #Assign password to be read from the file if it exists.
    if (defined ($$options{password_file}) && -s $$options{password_file} ) {
    	my $pass_file = $$options{password_file};
	open PASS, $pass_file or die ("Cannot open password file $pass_file : $!\n");
	print STDERR ("Password from file will take priority over --password option\n") if (defined $password);
	$password= <PASS>;
	chomp $password;
	close PASS;
    }

    ## handle some defaults
    $$options{server} = 'localhost' unless defined $$options{server};
    $$options{database_type} = 'mysql' unless defined $$options{database_type};
    $$options{translation_table} = 11 unless defined $$options{translation_table};
    $$options{source} = '.' unless defined $$options{source};
    $$options{cgi_mode} = 0 unless defined $$options{cgi_mode};
    $$options{intergenic_regions} = 0 unless defined $$options{intergenic_regions};

    ## format should be either gff or gbk
    if ( $$options{format} ne 'gbk' &&
         $$options{format} ne 'gff' &&
         $$options{format} ne 'tbl' ) {
        die "value for --format option must be either 'gbk' or 'gff' or 'tbl' ";
    }

    ## reset in case the user actually passes 'all'
    $$options{organism_id} = undef if $$options{organism_id} eq 'all';

    ## database type must be either mysql or postgresql
    if ( $$options{database_type} ne 'mysql' && $$options{database_type} ne 'postgresql' ) {
        die "value for option --database_type must be either 'mysql' or 'postgresql'";
    }
}

=head1 NAME

    chado_aengine_dumper.pl - creates a genbank or GFF file from a prokaryotic annotation stored
    in chado.

=head1 SYNOPSIS

  USAGE: assign_dbxref_ids.pl
    --database=abs2
    --user=someuser
    --password=somepass
    --output_directory=/some/dir
    --format=gbk
    [ --locus_db=IGS_bba1
    --locus_db_version=current
    --translation_table=11
    --server=khan.igs.umaryland.edu
    --database_type=postgreqsql
    --comment="This is a temporary/private release"
    --source=IGS
    --organism_id=all
    --use_assembly_names
    --cgi_mode=1
    --log=/path/to/some.log ]


=head1 OPTIONS

B<--database>
    Database name to connect to.

B<--user>
    User account with select, insert, and update privileges on the specified database.

B<--password>
    Password for user account specified.

B<--password_file>
    File which contains the password you want to provide.  Useful if you do not want to publicly reveal your password

B<--output_directory>
    Output directory.

B<--format>
    Output format desired.  Value should be either 'gbk', 'gff' or 'tbl'.

B<--locus_db>
    Optional.  If passed, corresponds to the db.name entry that should be used as the source
    for locus identifiers.  All other dbxref entries will be encoded as db_xref GBK attributes

B<--locus_db_version>
    Optional. Must specify if using locus_db option. This will specify the version of the
    dbxrefs used for locus_tag.

B<--translation_table>
    Optional.  Numeric value for translations table (default = 11)

B<--intergenic_regions>
	Optional.  Add intergenic region features to the output when enabled.

B<--server>
    Optional.  Server to connect to (default = localhost).

B<--database_type>
    Optional.  Database type (vendor.)  Currently supports 'mysql' or 'postgresql' (default = mysql).

B<--comment>
    Optional.  Populates the COMMENT field within each file written. (default = none).

B<--source>
    Optional.  Only used when exporting GFF - populates column 2 for all gene and CDS features. (default = none).

B<--organism_id>
    Optional.  Only process a specific organism_id (default = all).

B<--use_assembly_names>
    Optional flag. If included, will use the feature.name field for ids for assemblies rather than uniquenames.

B<--cgi_mode>
    Optional.  Transforms output mode of script to be suitable for CGI execution.  All output is written as a
    single stream to STDOUT, with a header that forces the user to download as a file.

B<--log>
    Optional.  Full path to a log file to create.

B<--help>
    This help message

=head1  DESCRIPTION

    description needed

=head1  INPUT

    Schema assumptions: Generally, this expects that an Annotation Engine genome has been run and
    loaded using Ergatis with representation in chado as defined by Coati/Prism loader conventions.

    Specifically, reference sequences are of type 'assembly' and CDS  feature 'derived_from' their
    transcripts via the feature_relationship table.

    Almost all features are only exported if the the values of their is_obsolete and is_analysis fields
    in the feature table are 0.  The only exception currently are rRNA features, which can be either.

    =head1  OUTPUT

    One file is written per assembly/sequence and is named using the feature.uniquename of the
    sequence with the gbk|gff extension added.  The Genbank output is probably a bit more complete -
    the addition is GFF is a recent feature and usage conventions are still being worked out.

    Example output of each feature type (double-spaces between each block added only for ease of reading
    in this documentation.)

    1830	2966	gene
    gene	dnaN
    locus_tag     OBB_0002
    1830	2966	CDS
    product	DNA-directed DNA polymerase III beta chain
    EC_number	2.7.7.7
    protein_id	gnl|ncbi|OBB_0002

    91493	93058	gene
    gene	rrsA
    locus_tag     OBB_0089
    91493	93058	rRNA
    product	16S ribosomal RNA


    96468	96543	gene
    gene	trnV
    locus_tag     OBB_0092
    96468	96543	tRNA
    product	tRNA-Val


    This is taken from:

  http://www.ncbi.nlm.nih.gov/Genbank/genomesubmit-Examples.html#fig2



=head1 TEST INFO

    bioperl howto:

http://www.bioperl.org/wiki/HOWTO:Feature-Annotation#Building_Your_Own_Sequences

    Entrez entry using as reference:

  http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nuccore&id=157844830

    commands used for testing:

    ./chado_aengine_dumper.pl --database=b_theta --database_type=postgresql --user=driley --password=whatever --organism_id=1 --server=driley-lx.igs.umaryland.edu --output_directory=/tmp/gbktest --locus_db='Btheta_' --source=IGS
    ./chado_aengine_dumper.pl --database=hik|cgsp --database_type=mysql --user=jorvis --password=whatever --organism_id=1 --server=localhost --output_directory=/tmp/gbktest --locus_db='TIGR_moore' --source=IGS --format=gff

    db_id 1 feature counts:

    name     | count
    -------------+-------
    gene        |  6028
    exon        |  6028
    CDS         |  6028
    tRNA        |    77
    assembly    |     2
    polypeptide |  5951
    transcript  |  5951

    =head1 TO DO

    - need to respect is_obsolete for all feature queries.  mysql and postgres treat these differently

    - BRC column 9 convention example:

    ID=bhbCDS.3517814;Name=FTF0002;gene_symbol=dnaN;Ontology_term=GO:0003677,GO:0003887,GO:0005737,GO:0006260,GO:0006261,GO:0008408,GO:0009360;Dbxref=NCBI_NP:YP_666216.1,Swiss-Prot:Q14K60,EC:2.7.7.7;Parent=bhb|gene.3517811;description=DNA polymerase III beta chain

    should test with gff3_to_annotab (at aaron's request)

    CVTERMS ( in gcj1099 )

    gene_product_name = 68
    frameshift = 34754
    point_mutation = 35139


=head1  CONTACT

    Joshua Orvis
    jorvis@users.sf.net

=cut
