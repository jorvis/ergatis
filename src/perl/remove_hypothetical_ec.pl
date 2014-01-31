#!/usr/bin/env perl

=head1 NAME

remove_hypothetical_ec.pl - Will remove EC numbers from genes annotated as hypothetical proteins

=head1 SYNOPSIS

 USAGE: remove_hypothetical_ec.pl
	--database_file=/path/to/db_file.txt
	--username=username
	--password=password
	--server=server
	--update=1
	--help

=head1 OPTIONS

B<--username,-u> 
	Database username

B<--password,-p>
	Database password

B<--server,-s>
	Database hostname
	
B<--database_file,-d>
    Either a comma-separated or tab-separated file containing the following columns:
    1) Name of Database
    2) Locus_tag ID prefix
    3) Path to a curate_common_names rules file
    4+) Any other DB related information (to come later)	

B<--update,-a> 
	Default = 0. Set to non-zero to make changes to database. Will not change anything by default

B<--help,-h>
	Will print this message

=head1  CONTACT

    Kevin Galens
    kgalens@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use DBI;
use Pod::Usage;
use Data::Dumper;

##############################
my $dbh;
my $update = 0;
my $prepared_queries = {}; #Used to cache some queries we will be using multiple times
##############################

&check_options();

my $data = remove_ec_from_hypothetical();

$dbh->disconnect();

sub remove_ec_from_hypothetical {
	
	#Get hypo genes
	my @hypo_genes_with_ec = grep { @{$_} > 0 } map { &has_ec_number( $_ ) } &get_hypo_genes();

	if( $update ) {
		&remove_ec_annotations( $_ ) for( @hypo_genes_with_ec );
		$dbh->commit();
		print "Removed ".scalar(@hypo_genes_with_ec)." EC annotations\n";
	} else {
		print "Found ".scalar(@hypo_genes_with_ec)." hypothetical genes with EC annotations\n";
	}
}

sub remove_ec_annotations {
	my ($fc_ids) = @_;
	my $sth = $prepared_queries->{'delete_fc'};
	$sth->execute( $_ ) for( @{$fc_ids} );
}

sub has_ec_number {
	my ($feature_id) = @_;
	my $sth = $prepared_queries->{"has_ec_number_sth"};
	$sth->execute( $feature_id );
	my $res = $sth->fetchall_arrayref();
	
	[map { $_->[0] } @{$res}];
}

sub get_hypo_genes {

	my $gpn_term = &get_cvterm_id( 'gene_product_name' );
	my $transcript = &get_cvterm_id( 'transcript' );

	my @hypo_terms = ('hypothetical protein', 'conserved hypothetical protein');
	
	my $q = 
		"SELECT f.feature_id ".
		"FROM feature f, featureprop fp ".
		"WHERE f.feature_id = fp.feature_id ".
		"AND f.type_id = ? ".
		"AND fp.type_id = ? ".
		"AND fp.value IN (?,?)";
	my $sth = $dbh->prepare( $q );
	$sth->execute( $transcript, $gpn_term, @hypo_terms );

	my $res = $sth->fetchall_arrayref();
	map { $_->[0] } @{$res};
}

sub get_cvterm_id {
	my ($name) = @_;

	my $sth = $prepared_queries->{'get_cvterm_sth'};
	$sth->execute( $name );
	my $res = $sth->fetchall_arrayref();
	
	if( @{$res} > 1 ) {
		die("Found more than one cvterm for $name");
	} elsif( @{$res} == 0 ) {
		die("Could not find term for $name");
	}

	$res->[0]->[0];
}

sub cache_queries {
	my $get_cvterm_query = 
		"SELECT cvterm_id FROM cvterm ".
		"WHERE name = ?";
	$prepared_queries->{'get_cvterm_sth'} = $dbh->prepare( $get_cvterm_query );

	my $has_ec_number = 
		"SELECT fc.feature_cvterm_id ".
		"FROM feature_cvterm fc, cvterm c, cv ".
		"WHERE cv.name = 'EC' ".
		"AND cv.cv_id = c.cv_id ".
		"AND c.cvterm_id = fc.cvterm_id ".
		"AND fc.feature_id = ?";
	$prepared_queries->{'has_ec_number_sth'} = $dbh->prepare( $has_ec_number );

	my $delete_fc = 
		"DELETE FROM feature_cvterm WHERE feature_cvterm_id = ?";
	$prepared_queries->{'delete_fc'} = $dbh->prepare( $delete_fc );
}

sub check_options {
    
    my %options;
    my $results = GetOptions (\%options,
                              'database_file|d=s',
                              'server|s=s',
                              'username|u=s',
                              'password|p=s',
                              'update|a=s',
							  'help|h'
                              );
    
	if( $options{'help'} ) {
		&_pod;
	}

    my @reqs = qw(database_file server username password);
    foreach my $req( @reqs ) {
        die("Option $req is required") unless( $options{$req} );
    }

    open DB, $options{'database_file'} or die ("Cannot open database file " . $options{'database_file'} . "$!\n");
   	my $line = <DB>;
   	chomp $line;
   	my ($database, $rest) = split(/,|\t/, $line, 2);
    close DB;

    $update = 1 if( $options{'update'} );

    $dbh = &_connect( $options{'username'}, $options{'password'},
                      $database, $options{'server'} );

	&cache_queries();

}

sub _connect {
    my ($user, $password, $db, $server) = @_;

    my $dbh;
    eval {
        $dbh = DBI->connect("DBI:mysql:$db:$server", "$user", "$password",
                            { 
                                'RaiseError' => 1,
                                'AutoCommit' => 0,
                          } );
    };
    if( $@ ) {
        carp("Could not connect to database ".DBI->errstr);
    }

    $dbh->do("use $db");

    return $dbh;
}

sub _pod {
	pod2usage( {-exitval => 1, -verbose => 2, -output => \*STDERR });
}
