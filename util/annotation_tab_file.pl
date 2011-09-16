#!/usr/local/bin/perl

=head1 NAME

annotation_tab_file.pl - generates tab delimited file of annotation including common name, gene symbol, EC number, GO ids, and TIGR role ids

=head1 USAGE: 

    chado2flatfile.pl
    --database=pak26
    --host=khan.igs.umaryland.edu
    --username=user
    --password=password
    --output_file=/path/to/output_file.txt
    --filter_complete

=head1 OPTIONS

B<--database,-d>
    The name of the database

B<--host,-s>
    Server where the database is located

B<--username,-u>
    Username for database

B<--password,-p>
    Password for database

B<--output_file,-o>
    The output tab delimited file name.

B<--filter_complete,-c>
    Optional, if used, will only provide annotations for features which are 
    marked as complete.

=head1  OUTPUT

One file is written per organism and is named using the database name with the ".txt" extension added and provided as download option.

=head1  DESCRIPTION

* The out file will a tab delimited file and will some look something like this ...

    cgsp_3829       diguanylate cyclase (GGDEF) domain protein                      GO:0009975,GO:0009966   264,710
    cgsp_3830       AFG1-like ATPase family protein                         157
    cgsp_3831       transposase, Mutator family protein                             154
    cgsp_3832       hypothetical protein                            856
    cgsp_3833       alkaline phosphatase                            703
    cgsp_3834       ATP-dependent Clp protease, ATP-binding subunit ClpX            clpX    GO:0008462,GO:0051082,GO:0006510,GO:0009368,GO:0005524  95,138
    cgsp_3835       hypothetical protein                            856
    cgsp_3836       conserved hypothetical protein                          156
    cgsp_3837       conserved hypothetical protein                          156
    cgsp_3838       HTH-type transcriptional regulator prtR (Pyocin repressor protein)                              129
    cgsp_3839       RNA polymerase sigma-H factor (Sigma-30)                                703
    cgsp_3840       pyrimidine-specific ribonucleoside hydrolase rihB (Cytidine/uridine-specific hydrolase)                         703
    cgsp_3841       tryptophan synthase, alpha subunit      4.2.1.20        trpA    GO:0004834,GO:0000162   70
    cgsp_3842       hypothetical protein                            856
    ...
    ...
    ...

* Will print logging information to terminal if no log file is specified
    
    
=head1  CONTACT
    
    Kevin galens
    kgalens@som.umaryland.edu

=cut

use strict;
use warnings;
use lib "../lib";
use Aengine::AnnotationDB;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use DBI;
use Data::Dumper;
use File::OpenFile qw(open_file);

my $database;
my $host;
my $username;
my $password;
my ($output);
my %options;
my $filter_complete = 0;

my $results = GetOptions (\%options,
                          'database|d=s',
                          'host|h=s',
                          'username|u=s',
                          'password|p=s',
                          'output_file|o=s',
                          'filter_complete|c'
                          );

&check_options( \%options );

############### VARS ###############
my $flags = {
    'common_name' => 1,
    'ec_number' => 1,
    'gene_symbol' => 1,
    'go_terms' => 1,
    'tigr_roles' => 1,
    'transcript' => 0,
    'product_name_source' => 1,
    'fmin' => 1,
    'fmax' => 1,
    'strand' => 1
    };
my $indices = {
    'common_name' => 0,
    'gene_symbol' => 1,
    'ec_number' => 2,
    'go_terms' => 3,
    'tigr_roles' => 4,
    'transcript' => 5,
    'product_name_source' => 6,
    'fmin' => 7,
    'fmax' => 8,
    'strand' => 9
    };
####################################

my $db = new Aengine::AnnotationDB( "username" => "$username", "password" => "$password", "database" => "$database", "host" => "$host" );

my $transcripts = $db->get_all_transcripts;

my %molecules;
foreach my $t ( @{$transcripts} ) {
    my $fid = $t->[0];

    next if( $filter_complete && !$db->feature_is_complete( $fid ) );

    my $locus = $db->get_locus_id( $fid, 'TIGR_moore' );
    next unless( defined( $locus ) );
    my $tmp = [$locus, $fid];

    my ($mol) = $db->get_molecule_and_coords( $fid );
    
    $molecules{$mol} = [] unless( exists( $molecules{$mol} ) );
    my $ext_locus = $db->get_external_locus_id( $fid, 'NCBI_locus_tag' );
    push(@{$tmp}, $ext_locus) if( defined( $ext_locus ) );

    push( @{$molecules{$mol}}, $tmp );
}

my $file = open_file( $output, 'out' );

print $file "internal locus\texternal locus\tmolecule id\tcommon_name\tgene_symbol\tec_number\tgo_terms\tproduct_name_source\tfmin\tfmax\tstrand\tTIGR_roles\n";

foreach my $molecule ( keys %molecules ) {
    my @loci = @{$molecules{$molecule}};
    foreach my $locus ( sort by_locus @loci ) {

        print $file "$locus->[0]";
        print $file "\t";
        print $file $locus->[2] if( defined( $locus->[2] ) );
        print $file "\t$molecule";
        my @annotation = $db->get_annotation( $locus->[1] );
        die("Expected 10 columns got ".scalar(@annotation)) unless( @annotation == 10 );
        foreach my $term ( qw(transcript common_name gene_symbol ec_number go_terms product_name_source fmin fmax strand tigr_roles) ) {
            die("$term not in index") unless( defined( $indices->{$term} ) );
            print $file "\t$annotation[$indices->{$term}]" if( $flags->{$term} );
        }
        print $file "\n";
    }
}

sub by_locus {
    my $a_number = $1 if( $a->[0] =~ /(\d+)$/ );

    if( !defined( $b->[0] ) ) {
        
        print Dumper( $a );
        print Dumper( $b );
        exit(0);
    }

    my $b_number = $1 if( $b->[0] =~ /(\d+)$/ );
    return $a_number <=> $b_number;
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
    $output = $opts->{'output_file'} if( $opts->{'output_file'} );
    $filter_complete = 1 if( $opts->{'filter_complete'});

    if( $filter_complete ) {
        print "Filtering only complete features\n";
    }
}
