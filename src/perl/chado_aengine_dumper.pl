#!/usr/bin/perl

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
            --translation_table=11
            --server=khan.igs.umaryland.edu
            --database_type=postgreqsql
            --comment="This is a temporary/private release"
            --source=IGS
            --organism_id=all
            --cgi_mode=1
            --log=/path/to/some.log ]


=head1 OPTIONS

B<--database>
    Database name to connect to.

B<--user>
    User account with select, insert, and update privileges on the specified database.

B<--password>
    Password for user account specified.

B<--output_directory>
    Password for user account specified.

B<--format>
    Output format desired.  Value should be either 'gbk', 'gff' or 'tbl'.

B<--locus_db>
    Optional.  If passed, corresponds to the db.name entry that should be used as the source
    for locus identifiers.  All other dbxref entries will be encoded as db_xref GBK attributes

B<--translation_table>
    Optional.  Numeric value for translations table (default = 11)

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

=head1  OUTPUT

One file is written per assembly/sequence and is named using the feature.uniquename of the
sequence with the gbk|gff extension added.  The Genbank output is probably a bit more complete -
the addition is GFF is a recent feature and usage conventions are still being worked out.

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
ID=bhb|CDS.3517814;Name=FTF0002;gene_symbol=dnaN;Ont
ology_term=GO:0003677,GO:0003887,GO:0005737,GO:0006260,GO:0006261,GO:0008408,GO:0009360;Dbxref=NCBI_NP:YP_666216.1,Swiss-Prot:Q14K60
,EC:2.7.7.7;Parent=bhb|gene.3517811;description=DNA polymerase III beta chain

should test with gff3_to_annotab (at aaron's request)

=head1  CONTACT

    Joshua Orvis
    jorvis@users.sf.net

=cut

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
$|++;

my %options = ();
my $results = GetOptions (\%options, 
                          'database=s',
                          'user=s',
                          'password=s',
                          'output_directory=s',
                          'format=s',
                          'locus_db=s',
                          'translation_table=i',
                          'server=s',
                          'database_type=s',
                          'comment=s',
                          'source=s',
                          'organism_id=i',
                          'cgi_mode=i',
                          'log=s',
                          'help') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

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
my $dbh = DBI->connect($dsn, $options{user}, $options{password}, {PrintError=>1, RaiseError=>1} );

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
    SELECT feature_id, uniquename, residues
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
    $$assemblies{$$row{feature_id}} = { uniquename => $$row{uniquename} };
    
    $$assemblies{ $$row{feature_id} }{seq_obj} = Bio::Seq::RichSeq->new (
                                                        -seq => $$row{residues},
                                                        -display_id => $$row{uniquename},
                                                        -division => 'BCT',
                                                 );

    $$assemblies{ $$row{feature_id} }{seq_obj}->species( $bio_species );

    my $collection = new Bio::Annotation::Collection;
    
    if ( $options{comment} ) {
        my $comment = Bio::Annotation::Comment->new( -text => $options{comment} );
        $collection->add_Annotation( 'comment', $comment );
        
    }

    $$assemblies{ $$row{feature_id} }{seq_obj}->annotation($collection);
}

$assembly_selector->finish();

## qry for all genes on an assembly
$cvterm{gene} = get_cvterm_id('gene');
$cvterm{transcript} = get_cvterm_id('transcript');
$cvterm{tRNA} = get_cvterm_id('tRNA');
$cvterm{rRNA} = get_cvterm_id('rRNA');
$cvterm{exon} = get_cvterm_id('exon');

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
    SELECT t.uniquename, ec.accession
      FROM feature t
      JOIN feature_cvterm fc ON t.feature_id = fc.feature_id
      JOIN cvterm cv ON fc.cvterm_id = cv.cvterm_id
      JOIN cvterm_dbxref cd ON cv.cvterm_id = cd.cvterm_id
      JOIN dbxref ec ON cd.dbxref_id = ec.dbxref_id
      JOIN db ON ec.db_id = db.db_id
     WHERE db.name = 'EC'
       AND t.type_id = $cvterm{transcript}
       AND t.feature_id = ?
};
my $ec_num_selector = $dbh->prepare($qry);

## query to select all GO terms for a given feature id
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
my $go_term_selector = $dbh->prepare($qry);

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


## track parents for GFF export
## structure like $h{$cds_uniquename} = $transcript_uniquename;
my %cds_parents;
my %transcript_parents;

for my $feature_id ( keys %$assemblies ) {

    ## until I can find a feature sorting method within the bioperl objects I'm using
    #   this to store the feature references so they can be added to the SeqIO in order
    my @assembly_feats;

    ## get the transcript features.  these correspond to the CDS entries for proks in the gbk file
    ##  using 'CDS' would include the tRNAs, which are queried separately
    $feature_on_assembly_selector->execute( 0, 0, $feature_id, $cvterm{transcript} );
    
    while (my $row = $feature_on_assembly_selector->fetchrow_hashref ) {
        
        ## make a gene entry for this too
        my $gene = Bio::SeqFeature::Generic->new(
            -seq_id => $$row{uniquename},
            -primary => 'gene',
            -start => $$row{fmin} + 1,
            -end => $$row{fmax},
            -strand => $$row{strand},
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
            
            $cds = Bio::SeqFeature::Generic->new(
                -seq_id => $$cds_row{uniquename},
                -primary => 'CDS',
                -start => $$row{fmin} + 1,
                -end => $$row{fmax},
                -strand => $$row{strand},
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
            $gene->add_tag_value( 'product', $$product_row{value} );
            last;
        }

        ## gene symbol
        $gene_symbol_selector->execute( $$row{feature_id} );
        while ( my $genesym_row = $gene_symbol_selector->fetchrow_hashref ) {
            $cds->add_tag_value( 'gene', $$genesym_row{value} );
            $gene->add_tag_value( 'gene', $$genesym_row{value} );
            last;
        }
        
        ## EC number
        $ec_num_selector->execute( $$row{feature_id} );
        while ( my $ec_num_row = $ec_num_selector->fetchrow_hashref ) {
            $cds->add_tag_value( 'EC_number', $$ec_num_row{accession} );
            last;
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
        $cds->add_tag_value( 'note', join( "; ", @values ) );


        push @assembly_feats, $gene;
        push @assembly_feats, $cds;
    }

    ## now get tRNAs
    $feature_on_assembly_selector->execute( 0, 0, $feature_id, $cvterm{tRNA} );
    
    while (my $row = $feature_on_assembly_selector->fetchrow_hashref ) {
        
        my $cds = Bio::SeqFeature::Generic->new(
            -seq_id => $$row{uniquename} || 'UNKN0WN',
            -primary => 'tRNA',
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
                die("found more than on gene on tRNA feature ID $$row{feature_id}\n");
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
            $gene->add_tag_value( 'product', $$product_row{value} );
            last;
        }

        push @assembly_feats, $cds;
        push @assembly_feats, $gene;
    }    
    
    ## now get all rRNAs
    $feature_on_assembly_selector->execute( 0, 0, $feature_id, $cvterm{rRNA} );
    
    while( my $row = $feature_on_assembly_selector->fetchrow_hashref ) {
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
            $gene->add_tag_value( 'product', $$product_row{value} );
            last;
        }
        
        push @assembly_feats, $cds;
        push @assembly_feats, $gene;
    }
    
    ## now add all the features
    for my $feat ( sort { $a->start <=> $b->start } @assembly_feats ) {
         $$assemblies{$feature_id}{seq_obj}->add_SeqFeature( $feat );
         
         if ( $feat->primary_tag eq 'CDS' ) {
            
            ## bioperl could fail here for any number of reasons, mostly if our coordinates on contig ends
            ##  are messed up.  warn here, but don't die.
            eval {
                if ( length $feat->seq ) {
                    $feat->add_tag_value( 'translation', $feat->seq->translate->seq );
                }
            };
            
            if ($@) {
                print STDERR "WARN: failure to export sequence for feature at coordinates (" . $feat->start . '/' . $feat->end . ")\n";
            }
         }
    }
}

## if this is a CGI, we need to print a header:
if ( $options{cgi_mode} == 1 ) {

    ## quicker way than this to get YYYYMMDD?
    my ($year, $mon, $day) = (localtime time)[5,4,3];
    my $datestamp = sprintf("%d%02d%02d", $year + 1900, ++$mon, $day);

    my $download_file_name;
    
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
            $tbl_fh = \*STDOUT;
    } else {
        _log("INFO: writing $options{output_directory}/$options{database}.tbl");
        open($tbl_fh, ">$options{output_directory}/$options{database}.tbl") || die "can't create TBL output file: $!";
    }
}

## write out the annotations in either GBK, GFF or tbl format
for my $feature_id ( keys %$assemblies ) {
    if ( $options{format} eq 'gbk' ) {

        my $gbk;
        
        if ( $options{cgi_mode} == 1 ) {
            $gbk = Bio::SeqIO->new( -format => 'genbank', 
                                    -fh => \*STDOUT );        
        } else {
            _log("INFO: writing $options{output_directory}/$$assemblies{$feature_id}{uniquename}.gbk");
            $gbk = Bio::SeqIO->new( -format => 'genbank', 
                                    -file => ">$options{output_directory}/$$assemblies{$feature_id}{uniquename}.gbk" );

        }
        
        $gbk->write_seq( $$assemblies{$feature_id}{seq_obj} );
    
    } elsif ( $options{format} eq 'tbl' ) {
    
        print $tbl_fh ">Feature $$assemblies{$feature_id}{uniquename}\n";
        
        foreach my $feat ( $$assemblies{$feature_id}{seq_obj}->get_SeqFeatures ) {
            
            if ( $feat->primary_tag eq 'gene' ) {
                if ($feat->strand == 1 ) {
                    print $tbl_fh $feat->start, "\t", $feat->end, "\tgene\n";
                } else {
                    print $tbl_fh $feat->end, "\t", $feat->start, "\tgene\n";
                }

                if ( $feat->has_tag('gene') ) {
                    print $tbl_fh "\t\t\tgene\t", ($feat->get_tag_values('gene'))[0], "\n";
                }
                
                if ( $feat->has_tag('locus_tag') ) {
                    print $tbl_fh "\t\t\tlocus_tag\t", ($feat->get_tag_values('locus_tag'))[0], "\n";
                }

            } elsif ( $feat->primary_tag eq 'CDS' ) {
                if ($feat->strand == 1 ) {
                    print $tbl_fh $feat->start, "\t", $feat->end, "\tCDS\n";
                } else {
                    print $tbl_fh $feat->end, "\t", $feat->start, "\tCDS\n";
                }
                
                if ( $feat->has_tag('product') ) {
                    print $tbl_fh "\t\t\tproduct\t", ($feat->get_tag_values('product'))[0], "\n";
                }
                
                if ( $feat->has_tag('locus_tag') ) {
                    print $tbl_fh "\t\t\tprotein_id\tgnl|IGS|", ($feat->get_tag_values('locus_tag'))[0], "\n";
                }
            }
        }
    
    } elsif ( $options{format} eq 'gff' ) {
        
        my $gff_fh;
        
        if ( $options{cgi_mode} == 1 ) {
            $gff_fh = \*STDOUT;
        } else {
            _log("INFO: writing $options{output_directory}/$$assemblies{$feature_id}{uniquename}.gff3");
            open($gff_fh, ">$options{output_directory}/$$assemblies{$feature_id}{uniquename}.gff3") || die "can't create GFF output file: $!";
        }

        foreach my $feat ( $$assemblies{$feature_id}{seq_obj}->get_SeqFeatures ) {
            
            ## writing my own GFF output.  A Bio::SeqIO layer would be cool eventually.
            my @columns = ( $$assemblies{$feature_id}{uniquename}, 
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

                print_gff_col9_attribute($gff_fh, 'ID', $transcript_parents{$feat->seq_id} );

                ## now handle any column 9 attributes
                if ( $feat->has_tag('product') ) {
                    print_gff_col9_attribute($gff_fh, 'description', ($feat->get_tag_values('product'))[0] );
                }

                if ( $feat->has_tag('gene') ) {
                    print_gff_col9_attribute($gff_fh, 'gene_symbol', ($feat->get_tag_values('gene'))[0] );
                }

                print $gff_fh "\n";
            
            ## GFF requires exons too.  write an mRNA feature, then the CDS ones
            } elsif ( $feat->primary_tag eq 'CDS' ) {
                
                ## the mRNA feature matches the gene one
                $columns[2] = 'mRNA';
                $columns[7] = '.';
                print $gff_fh join("\t", @columns);
                print $gff_fh "\t";
                print_gff_col9_attribute($gff_fh, 'ID', $cds_parents{$feat->seq_id} );
                print_gff_col9_attribute($gff_fh, 'Parent', $transcript_parents{$cds_parents{ $feat->seq_id }} );
                print $gff_fh "\n";
            
                $exon_on_transcript_selector->execute( $feature_id, $cds_parents{ $feat->seq_id } );
                
                ## we'll write one CDS per exon, as described in the GFF3 specs
                while ( my $row = $exon_on_transcript_selector->fetchrow_hashref ) {
                    print $gff_fh "$$assemblies{$feature_id}{uniquename}\t" . 
                                  "$options{source}\tCDS\t" .
                                  ($$row{fmin} + 1) . "\t$$row{fmax}\t.\t";
                    
                    print $gff_fh ($$row{strand} == 1) ? '+' : '-';
                    print $gff_fh "\t0\t";
                    print_gff_col9_attribute($gff_fh, 'ID', $feat->seq_id );
                    print_gff_col9_attribute($gff_fh, 'Parent', $cds_parents{ $feat->seq_id } );
                    print $gff_fh "\n";
                    
                    print $gff_fh "$$assemblies{$feature_id}{uniquename}\t" . 
                                  "$options{source}\texon\t" .
                                  ($$row{fmin} + 1) . "\t$$row{fmax}\t.\t";
                    
                    print $gff_fh ($$row{strand} == 1) ? '+' : '-';
                    print $gff_fh "\t.\t";
                    print_gff_col9_attribute($gff_fh, 'ID', $$row{uniquename} );
                    print_gff_col9_attribute($gff_fh, 'Parent', $cds_parents{ $feat->seq_id } );
                    print $gff_fh "\n";
                    
                }
            
            }
        }
    }
}

$exon_on_transcript_selector->finish();
$cds_on_transcript_selector->finish();
$ec_num_selector->finish();
$gene_symbol_selector->finish();
$product_selector->finish();
$dbxref_selector->finish();
$feature_on_assembly_selector->finish();

$dbh->disconnect();

exit(0);

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

sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
}

sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( database user password output_directory format );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
    
    ## handle some defaults
    $$options{server} = 'localhost' unless defined $$options{server};
    $$options{database_type} = 'mysql' unless defined $$options{database_type};
    $$options{translation_table} = 11 unless defined $$options{translation_table};
    $$options{source} = '.' unless defined $$options{source};
    $$options{cgi_mode} = 0 unless defined $$options{cgi_mode};
    
    ## format should be either gff or gbk
    if ( $$options{format} ne 'gbk' && 
         $$options{format} ne 'gff' && 
         $$options{format} ne 'tbl' ) {
        die "value for --format option must be either 'gbk' or 'gff'";
    }
    
    ## reset in case the user actually passes 'all'
    $$options{organism_id} = undef if $$options{organism_id} eq 'all';
    
    ## database type must be either mysql or postgresql
    if ( $$options{database_type} ne 'mysql' && $$options{database_type} ne 'postgresql' ) {
        die "value for option --database_type must be either 'mysql' or 'postgresql'";
    }
}











