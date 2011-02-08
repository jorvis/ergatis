#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

chado2flatfile.pl - Description

=head1 SYNOPSIS

 USAGE: chado2flatfile.pl
       --database=pak26
       --host=khan.igs.umaryland.edu
       --username=user
       --password=password
       --pep_file=/path/to/file.pep
       --seq_file=/path/to/file.seq
       --coords_file=/path/to/file.coords
     [ --locus_ids=TIGR_moore
       --log=/path/to/file.log
       --debug=4
       --help
     ]

=head1 OPTIONS

B<--database,-d>
    The name of the database

B<--host,-s>
    Server where the database is located

B<--username,-u>
    Username for database

B<--password,-p>
    Password for database

B<--pep_file,-f>
    The output pep file name.

B<--seq_file,-e>
    The output seq file name.

B<--coords_file,-c>
    The output coords file name.

B<--locus_ids>
    The name of the db of the locus identifiers to be used in the output

B<--help,-h>
    Print this message

=head1  DESCRIPTION

    Reads information from a chado database and creates a pep file, seq file and/or coords file
    for the features described in the database.
 
=head1  INPUT
    
    The only input is the name of database, host, username and password.

=head1 OUTPUT

    All output files are sorted first by molecule located on and then by featureloc.fmax
    
    pep file:
    Multi fasta file of peptides. Headers have identifiers (feature.uniquename, unless --locus_ids option is used)
    followed by the common name (featureprop.value for 'gene_product_name' cvterm) then followed by the
    molecule uniquenae on which the feature is located.

    ex:
    >hik_1522	glyceraldehyde-3-phosphate dehydrogenase, type I hik.assembly.1
    MAIKIGINGFGRIGRIVFRAAQHRDDIEVVGINDLIDVEYMAYMLKYDSTHGRFDGTVEV
    KDGNLVVNGKTIRVTAERDPANLNWGAIGVDIAVEATGLFLTDETARKHITAGAKKVVLT
    GPSKDATPMFVRGVNFNAYAGQDIVSNASCTTNCLAPLARVVHETFGIKDGLMTTVHATT
    ATQKTVDGPSAKDWRGGRGASQNIIPSSTGAAKAVGKVLPALNGKLTGMAFRVPTPNVSV
    VDLTVNLEKPASYDAIKQAIKDAAEGKTFNGELKGVLGYTEDAVVSTDFNGCALTSVFDA
    DAGIALTDSFVKLVSWYDNETGYSNKVLDLVAHIYNYKG
    >hik_1013	hypothetical protein
    LKKLTALRSVFYYLEFNLRSKNEQGITKNNLKQSIFF
    >hik_593	AMP-binding enzyme family protein
    MNLDLHFVHRIQQQAKTRTNMTALRYKEHGLWRDISWKNFQEQLNQLSRALLAHNIDVQD
    KIAIFAHNMERWTIVDIATLQIRAITVPIYATNTAQQAEFILNHADVKILFVGDQEQYDQ
    TLEIAHHCPKLQKIVAMKSTIQLQQDPLSCTWESFIKTGSNAQQDELTQRLNQKQLSDLF
    TIIYTSGTTGEPKGVMLDYANLAHQLETHDLSLNVTDQDISLSFLPFSHIFERAWAAYIL
    HRGAILCYLEDTNQVRSALTEIRPTLMCAVPRFYEKIYAAVLDKVQKAPKLRQIMFHWAI
    SVGQKYFDLRANNKAIPFLLKKQFALADKLVLSKLRQLLGGRIKMMPCGGAKLEPAIGLF
    FHAIGINIKLGYGMTETTATVSCWHDFQFNPNSIGTLMPKAEVKIGENNEILVRGGMVMK
    GYYKKPEETAQAFTEDGFLKTGDAGEFDEQGNLFITDRIKELMKTSNGKYIAPQYIESKI
    GKDKFIEQIAIIADAKKYVSALIVPCFDSLEEYAKQLNIKYHDRLELLKNSDILKMFEHR
    INAVQKELAHFEQVKKFTLLSQAFSIKLGEITPTLKLRRKVILERYRKQIEAMYHSQEA
    ...

    seq file:
    Same as pep file but with nucleotide sequences

    coords file:
    tab separated file with identifiers (uniquename unless --locus_ids option is used)
    followed by molecule name, end5 and end3.

    ex:
    hik_1522	hik.assembly.1   1	1021
    hik_1013	hik.assembly.1   1029	1143
    
    ...
    

=head1  CONTACT

    Kevin Galens
    kgalens@som.umaryland.edu

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use POSIX;
BEGIN {
use Ergatis::Logger;
}
use DBI;
use Data::Dumper;
use File::OpenFile qw(open_file);

my $database;
my $host;
my $username;
my $password;
my ($pep, $seq, $coords);
my %options;
my $peptide_info;
my $locus = 0;

my $results = GetOptions (\%options,
                          'database|d=s',
                          'host|h=s',
                          'username|u=s',
                          'password|p=s',
                          'seq_file|e=s',
                          'pep_file|f=s',
                          'coords_file|c=s',
                          'locus_ids|l=s',
                          'help'
                          );

&check_options( \%options );

my $dbh = &connect_to_db( $database, $host, $username, $password );
#print "querying peptide info\n";
$peptide_info = &get_polypeptide_info();

my $locus_lookup = {};
$locus_lookup = &get_polypeptide_to_locus_lookup( ) if( $locus );

if( $pep ) {
#    print "creating .pep file\n";
    &print_pep_file( $peptide_info, $pep );
}

if( $seq ) {
#    print "creating .seq file\n";
    my $cds_info = &get_cds_info( );
    &print_seq_file( $peptide_info, $cds_info, $seq );
}

if( $coords ) {
#    print "creating .coords file\n";
    &print_coords_file( $peptide_info, $coords );
}

sub get_polypeptide_to_locus_lookup {

    #check to see if the db is there
    my $d_query = "SELECT db_id ".
        "FROM db ".
        "WHERE name = '$locus'";

    my $dth = $dbh->prepare($d_query);
    $dth->execute();

    my $tmp_aref = $dth->fetchall_arrayref();
    $dth->finish();
    die("Could not find db $locus in database. Please check the database name for the locus identifiers desired ".
        "and change the value of the --locus_ids option") 
        if( @{$tmp_aref} == 0 );

    #make gene -> polypeptide lookup
    my $t_query = "SELECT g.uniquename, p.uniquename ".
        "FROM feature g, feature t, feature p, cvterm cg, cvterm ct, cvterm cp, feature_relationship frgt, feature_relationship frtp ".
        "WHERE g.type_id = cg.cvterm_id ".
        "AND cg.name = 'gene' ".
        "AND g.feature_id = frgt.object_id ".
        "AND t.type_id = ct.cvterm_id ".
        "AND ct.name = 'transcript' ".
        "AND t.feature_id = frgt.subject_id ".
        "AND t.feature_id = frtp.object_id ".
        "AND p.type_id = cp.cvterm_id ".
        "AND cp.name = 'polypeptide' ".
        "AND p.feature_id = frtp.subject_id ".
        "AND t.is_obsolete = 0";

    my $tth = $dbh->prepare( $t_query );
    $tth->execute();
    
    my $lookup = {};
    while( my $row = $tth->fetchrow_arrayref() ) {
        $lookup->{$row->[0]} = $row->[1];
    }
    $tth->finish();

    #grab the loci
    my $l_query = "SELECT g.uniquename, x.accession ".
        "FROM feature g, cvterm cg, dbxref x, feature_dbxref fd, db ".
        "WHERE g.type_id = cg.cvterm_id ".
        "AND cg.name = 'gene' ".
        "AND g.feature_id = fd.feature_id ".
        "AND fd.dbxref_id = x.dbxref_id ".
        "AND x.db_id = db.db_id ".
        "AND db.name = '$locus' ".
        "AND g.is_obsolete = 0 ";

    my $lth = $dbh->prepare( $l_query );
    $lth->execute();
    
    my $retval = {};
    while( my $row = $lth->fetchrow_arrayref() ) {
        my $pid = $lookup->{$row->[0]};
        next unless( defined( $pid ) ); #we don't want RNAs
        $retval->{$pid} = $row->[1];
    }

    $lth->finish();

    return $retval;
    
}

sub print_pep_file {
    my ($data, $file) = @_;

    my $pep = open_file( $file, 'out' );

    foreach my $polypeptide ( sort by_molecule_then_fmax values %{$data} ) {
        my $uniquename = $polypeptide->{'uniquename'};
        my $print_id = $uniquename;
        $print_id = $locus_lookup->{$uniquename} if( $locus );

        my $seq = $polypeptide->{'residues'};
        
        #skip if we don't have any sequence information
        unless( defined( $seq ) && length($seq) > 0 ) {
            warn("Could not find sequence information for polypeptide $uniquename");
            next;
        }

        print $pep ">$print_id\t".$polypeptide->{'value'}."\t".$polypeptide->{'molname'}."\n";
        print $pep $1."\n" while( $seq =~ /(\w{1,60})/g );
    }
    
    close($pep);
    
}

sub by_molecule_then_fmax {
    my $retval;
    unless( $retval = ( $a->{'molname'} cmp $b->{'molname'} ) ) {
        $retval = ( $a->{'fmax'} <=> $b->{'fmax'} );
    }
    return $retval;
}

sub print_coords_file {
    my ($data, $file) = @_;
    
    my $crd = open_file( $file, 'out' );

    foreach my $polypeptide ( sort by_molecule_then_fmax values %{$data} ) {
        my $uniquename = $polypeptide->{'uniquename'};
        my $print_id = $uniquename;
        $print_id = $locus_lookup->{$uniquename} if( $locus );
        print $crd "$print_id\t".$polypeptide->{'molname'}."\t";
        my ($start, $end) = ( $data->{$uniquename}->{'strand'} == -1 ) ? 
            ( $data->{$uniquename}->{'fmax'}, $data->{$uniquename}->{'fmin'} ) :
            ( $data->{$uniquename}->{'fmin'}, $data->{$uniquename}->{'fmax'} );
        print $crd "$start\t$end\n";
    }

    close($crd);
    
}

sub print_seq_file {
    my ($pep_data, $cds_data, $file) = @_;
    my $seq = open_file( $file, 'out' );

    foreach my $polypeptide ( sort by_molecule_then_fmax values %{$pep_data} ) {
        my $uniquename = $polypeptide->{'uniquename'};
        my $cds_unq = $cds_data->{$uniquename}->{'CDS'} || die "CDS not defined for uniquename $uniquename";
        my $cds_seq = $cds_data->{$uniquename}->{'seq'};

        unless( defined( $cds_seq ) && length( $cds_seq ) > 0 ) {
            warn( "Could not find nucleotide sequence for $uniquename ($cds_unq). Skipping..." );
            next;
        }

        my $print_id = $cds_unq."\t".$uniquename;
        $print_id = $locus_lookup->{$uniquename} if( $locus );
        die("Could not find locus id for $uniquename") if( !defined( $print_id ) );

        print $seq ">$print_id\t".$pep_data->{$uniquename}->{'value'}."\t".$cds_data->{$uniquename}->{'molecule'}."\n";
        print $seq $1."\n" while( $cds_seq =~ /(\w{1,60})/g );
    }
    
    close($seq);
    
}

sub print_data_files {
    my ($basename, $data) = @_;

    my $pep = $basename.".pep";
    my $coords = $basename.".coords";

    open( PEP, ">$pep") or die("Could not open $pep for writing ($!)");
    open( CRD, ">$coords") or die("Could not open $coords for writing ($!)");

    foreach my $uniquename ( sort { $data->{$a}->{'fmin'} <=> $data->{$b}->{'fmax'} } keys %{$data} ) {
        print PEP ">$uniquename\t".$data->{$uniquename}->{'value'}."\n";
        my $seq = $data->{$uniquename}->{'residues'};
        print PEP $1."\n" while( $seq =~ /(\w{1,60})/g );

        print CRD "$uniquename\t";
        my ($start, $end) = ( $data->{$uniquename}->{'strand'} == -1 ) ? 
            ( $data->{$uniquename}->{'fmax'}, $data->{$uniquename}->{'fmin'} ) :
            ( $data->{$uniquename}->{'fmin'}, $data->{$uniquename}->{'fmax'} );
        print CRD "$start\t$end\n";
    }

    close(PEP);
    close(CRD);

    print "$pep\n";
    print "$coords\n";

}

sub connect_to_db {
    my ($db, $host, $user, $pass) = @_;
    my $retval = DBI->connect("dbi:mysql:host=$host;packetSize=8092", $user, $pass) or
        die("Could not connect to server");
    $retval->do("use $db");
    return $retval;
}

sub get_cds_info {
    
    my $query = "SELECT p.uniquename, c.uniquename, c.residues, m.uniquename ".
        "FROM feature c, feature p, feature t, cvterm ct, cvterm cc, cvterm cp, feature_relationship frpt, feature_relationship frct, feature m, featureloc fl ".
        "WHERE frpt.subject_id = p.feature_id ".
        "AND p.type_id = cp.cvterm_id ".
        "AND cp.name = 'polypeptide' ".
        "AND frpt.object_id = t.feature_id ".
        "AND t.type_id = ct.cvterm_id ".
        "AND ct.name = 'transcript' ".
        "AND frct.object_id = t.feature_id ".
        "AND frct.subject_id = c.feature_id ".
        "AND c.type_id = cc.cvterm_id ".
        "AND cc.name = 'CDS' ".
        "AND c.is_obsolete = 0 ".
        "AND p.feature_id = fl.feature_id ".
        "AND fl.srcfeature_id = m.feature_id";

    my $sth = $dbh->prepare( $query );
    $sth->execute;

    my $retval = {};

    while( my $row = $sth->fetchrow_arrayref ) {
        $retval->{$row->[0]}->{'CDS'} = $row->[1];
        $retval->{$row->[0]}->{'seq'} = $row->[2];
        $retval->{$row->[0]}->{'molecule'} = $row->[3];
    }

    return $retval;
}

sub get_polypeptide_info {

    my $query = " SELECT p.uniquename, p.residues, fp.value, fl.fmin, fl.fmax, fl.strand, m.uniquename as molname ".
        "FROM featureloc fl, feature p, feature t, featureprop fp, feature_relationship fr, ".
        "cvterm cg, cvterm cp, cvterm ct, cvterm cfr, feature m ".
        "WHERE p.type_id = cp.cvterm_id ".
        "AND cp.name = 'polypeptide' ".
        "AND fr.subject_id = p.feature_id ".
        "AND t.feature_id = fr.object_id ".
        "AND fr.type_id = cfr.cvterm_id ".
        "AND cfr.name = 'part_of' ".
        "AND t.type_id = ct.cvterm_id ".
        "AND ct.name = 'transcript' ".
        "AND t.feature_id = fp.feature_id ".
        "AND fp.type_id = cg.cvterm_id ".
        "AND cg.name = 'gene_product_name' ".
        "AND p.feature_id = fl.feature_id ".
        "AND p.is_obsolete = 0 ".
        "AND fl.srcfeature_id = m.feature_id";

    my $sth = $dbh->prepare( $query );
    $sth->execute();

    my $retval = $sth->fetchall_hashref('uniquename');
#    print "Fetched the rows\n";
    return $retval;
    
}

sub check_options {
    my $opts = shift;

    my @required_options = qw( database host username password );
    foreach my $req( @required_options ) {
        die("Option $req is required") unless( exists( $opts->{$req} ) );
    }

    $database = $opts->{'database'};
    $host = $opts->{'host'};
    $username = $opts->{'username'};
    $password = $opts->{'password'};
    $pep = $opts->{'pep_file'} if( $opts->{'pep_file'} );
    $seq = $opts->{'seq_file'} if( $opts->{'seq_file'} );
    $coords = $opts->{'coords_file'} if( $opts->{'coords_file'} );
    $locus = $opts->{'locus_ids'} if( $opts->{'locus_ids'} );
}
