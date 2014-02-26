#!/usr/bin/env perl

=head1 NAME

curate_common_names.pl - 

=head1 SYNOPSIS

 USAGE: curate_common_names.pl
       --input_map_file=/path/to/map.txt
       --username=kgalens
       --password=pass
       --host=server
       --no_change=1
     [ --password_file=/path/to/password.txt
       --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_map_file,-i>
REQUIRED. Either a comma-separated or tab-separated file containing the following columns:
    1) Name of Database
    2) Locus_tag ID prefix
    3) Path to a curate_common_names rules file
    4+) Any other DB related information (to come later)

B<--username,-u>
REQUIRED. Database username

B<--password,-p>
Database password.  Must provide either --password or --password_file

B<--password_file,-P>
Password stored in a text file.  Useful for keeping password non-visible.

B<--host,-s>
REQUIRED. Database server name

B<--no_change,-n>
Optional. Default = 0. If set to non-zero, will not change database. 

B<--log,-l>
    Logfile.

B<--debug,-D>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

 There were name that were being assigned like this:

    functional name (silly info).

 Turns out, NCBI doesn't like stuff in parens. So they
 should be removed. Also removing trailing periods.
 
=head1  INPUT
    Just the database name, username, password and server name

=head1 OUTPUT
    None. Except for log file.

=head1  CONTACT

    Kevin Galens
    kgalens@gmail.com

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use DBI;
use Text::Balanced qw(extract_bracketed);
use Pod::Usage;
use Data::Dumper;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
my $dbh;
my $no_change = 0;
my %checks = ();
my %rules = ();
my %gene_sym_checks = ();
my $chp = 'conserved hypothetical protein';
my $conserved_hypo = {
    'gene_product_name' => "conserved hypothetical protein",
    'GO'                => ['GO:0008150', 'GO:0003674', 'GO:0005575'],
    'TIGR_Role'         => 156,
};
my ($total_changed, $total_ch_hypo) = (0,0);
my $password;
my $database;
my $rules_file;
####################################################

my %options;
my $results = GetOptions (\%options,
                          'username|u=s',
                          'password|p=s',
                          'password_file|P=s',
                          'host|s=s',
			  'input_map_file|i=s',
                          'no_change|n=s',
                          "log|l=s",
                          "debug|D=s",
                          "help|h"
                          );

######################################################################
##                              Main                                ##
######################################################################
&check_options(\%options);

if(-e $rules_file) {
	open(RF, "< $rules_file") or die "Could not open $rules_file file for reading\n";
	my @contents = <RF>;
	close(RF);
	foreach my $line (@contents) {
		chomp($line);
		next if($line =~ /^#/);
		next if ($line =~ /^\s*$/);
		my ($correct,$incorrect) = split(/\t/,$line,2);
		#print $correct . "-----".$incorrect."\n";
		if(!exists($rules{$incorrect})) {
			$rules{$incorrect} = $correct;
		}
	}
} else {
	die "$rules_file file does not exist";
}

%checks = &register_common_name_checks();
%gene_sym_checks = &register_gene_symbol_checks();

print STDERR "Now curating database $database...\n";
$dbh = &_connect( $options{'username'}, $password, 
                 $database, $options{'host'} );

#Grab all transcripts, with common names
my @transcripts = &get_all_transcripts($dbh, 'gene_product_name');
# Get all transcripts with gene symbol
my @gene_trans = &get_all_transcripts($dbh, 'gene_symbol');

print STDERR "Found ".scalar(@transcripts)."\n";
print STDERR "Found ".scalar(@gene_trans)."\n";
#Find the ones to be changed
my %to_change = &change_com_name('gene_product_name', @transcripts);
my %gene_to_change = &change_com_name('gene_symbol', @gene_trans);

print STDERR "Changing ".scalar(keys %to_change)." common name\n";
print STDERR "Changing ".scalar(keys %gene_to_change)." gene symbol\n";

#And make the changes
unless( $no_change ) {
    	&load_changes( $dbh, \%to_change, 'gene_product_name');
    	&load_changes($dbh, \%gene_to_change, 'gene_symbol');
}else {
    	print Dumper( \%to_change );
    	print Dumper(\%gene_to_change);
}

print STDERR "Changed the annotation of $total_changed transcripts\n";
print STDERR "Changed $total_ch_hypo to conserved hypothetical\n";

$dbh->disconnect();
exit(0);

######################################################################
##                           Subs                                   ##
######################################################################
sub load_changes {
    my ($dbh, $data, $filter) = @_;

    my $query = 
        "UPDATE featureprop SET value = ? ".
        "WHERE featureprop_id = ?";
    my $sth = $dbh->prepare( $query );

    foreach my $tid ( keys %{$data} ) {
        $total_changed++;
        if(($data->{$tid}->{'new'} eq $chp) || ($data->{$tid}->{'new'} eq "hypothetical protein")) {
	    if($filter eq "gene_product_name") {
            	$total_ch_hypo++;
            	print $data->{$tid}->{'old'}."\t=>\tCONSERVED HYPO\n";
            	&load_conserved_hypothetical( $tid, $data->{$tid}->{'featureprop_id'}, $data->{$tid} );
	    }
        } else {
            print $data->{$tid}->{'old'}."\t=>\t".$data->{$tid}->{'new'}."\n";
            $sth->execute( $data->{$tid}->{'new'}, $data->{$tid}->{'featureprop_id'} );
        }
    }
    $sth->finish();
    $dbh->commit();
}
sub load_conserved_hypothetical {
    my ($tid, $fid, $data) = @_;
    die("Transcript id was not defined") unless( defined( $tid ) );
    die("Featureprop id was not defined") unless( defined( $fid ) );
    die("data was not defined") unless( defined( $data ) );
    
    
    #delete go terms
    &delete_go_terms( $tid );
    
    #delete tigr roles
    &delete_tigr_roles( $tid );
    
    #delete gene symbol
    &delete_gene_symbol( $tid );

    #delete ec
    &delete_ec_number( $tid );

    #change common name
    &update_com_name( $fid, $data->{'new'} );

    #add go terms
    foreach my $go ( @{$conserved_hypo->{'GO'}} ) {
        &add_go_term( $tid, $go );
    }

    #add tigr role
    &add_tigr_role( $tid, $conserved_hypo->{'TIGR_Role'} );
    
}
sub add_go_term {
    my ($tid, $go) = @_;
    
    #Find cvterm id for go term
    my $query = 
        "SELECT c.cvterm_id FROM cvterm c, dbxref d, db ".
        "WHERE d.accession = ? and d.dbxref_id = c.dbxref_id ".
        "AND db.db_id = d.db_id ".
        "AND db.name IN ('biological_process', 'molecular_function', 'cellular_component', 'GO')";
    my $sth = $dbh->prepare( $query );
    $sth->execute( $go );
    my $res = $sth->fetchall_arrayref();
    my $cvterm_id;
    if( @{$res} > 0 ) {
        $cvterm_id = $res->[0]->[0];
    } else {
        die("Could not find go term $go");
    }

    #add the proper entry into feature_cvterm
    my $next_id = &next_id( 'feature_cvterm', 'feature_cvterm_id' );
    $query = 
        "INSERT INTO feature_cvterm (feature_cvterm_id, feature_id, cvterm_id, pub_id, is_not) ".
        "VALUES( ?, ?, ?, 1, 0)";
    $sth = $dbh->prepare( $query );
    $sth->execute( $next_id, $tid, $cvterm_id );
}
sub add_tigr_role {
    my ($tid, $role) = @_;
    
    #Find cvterm id for go term
    my $query = 
        "SELECT c.cvterm_id FROM cvterm c, dbxref d, db, cvterm_dbxref cd ".
        "WHERE d.accession = ? AND d.db_id = db.db_id ".
        "AND db.name = 'TIGR_role' AND d.dbxref_id = cd.dbxref_id ".
        "AND cd.cvterm_id = c.cvterm_id";
    my $sth = $dbh->prepare( $query );
    $sth->execute( $role );
    my $res = $sth->fetchall_arrayref();
    my $cvterm_id;
    if( @{$res} > 0 ) {
        $cvterm_id = $res->[0]->[0];
    } else {
        die("Could not find tigr role $role");
    }

    #add the proper entry into feature_cvterm
    my $next_id = &next_id( 'feature_cvterm', 'feature_cvterm_id' );
    $query = 
        "INSERT INTO feature_cvterm (feature_cvterm_id, feature_id, cvterm_id, pub_id, is_not) ".
        "VALUES( ?, ?, ?, 1, 0)";
    $sth = $dbh->prepare( $query );
    $sth->execute( $next_id, $tid, $cvterm_id );
}
sub next_id {
    my ($table, $col) = @_;
    my $sth = $dbh->prepare( "SELECT max( $col ) FROM $table" );
    $sth->execute();
    my $res = $sth->fetchall_arrayref();
    my $max = $res->[0]->[0];
    ++$max;
}
sub delete_gene_symbol {
    my ($tid) = @_;
    my $sth = $dbh->prepare( "SELECT fp.featureprop_id ".
                             "FROM featureprop fp, cvterm c ".
                             "WHERE fp.feature_id = $tid ".
                             "AND fp.type_id = c.cvterm_id ".
                             "AND c.name = 'gene_symbol'" );
    $sth->execute();
    my $delete_gs = $dbh->prepare( "DELETE FROM featureprop ".
                                   "WHERE featureprop_id = ?" );
    while( my $row = $sth->fetchrow_arrayref() ) {
        $delete_gs->execute( $row->[0] );
    }
    $sth->finish();          
}
sub delete_ec_number {
    my ($tid) = @_;
    my $query = 
        "SELECT fc.feature_cvterm_id ".
        "FROM feature_cvterm fc, cvterm c, cv ".
        "WHERE fc.feature_id = $tid ".
        "AND fc.cvterm_id = c.cvterm_id ".
        "AND c.cv_id = cv.cv_id ".
        "AND cv.name = 'EC'";
    my $sth = $dbh->prepare( $query );

    my $delete_fc = $dbh->prepare( "DELETE FROM feature_cvterm WHERE feature_cvterm_id = ?");
    $sth->execute();
    while( my $row = $sth->fetchrow_arrayref() ) {
        $delete_fc->execute( $row->[0] );
    }
    $sth->finish;
}
sub delete_go_terms {
    my ($tid) = @_;
    my $query = 
        "SELECT fc.feature_cvterm_id ".
        "FROM feature_cvterm fc, cvterm c, cv ".
        "WHERE fc.feature_id = $tid ".
        "AND fc.cvterm_id = c.cvterm_id ".
        "AND c.cv_id = cv.cv_id ".
        "AND cv.name IN ('biological_process', 'molecular_function', 'cellular_component')";
    my $sth = $dbh->prepare( $query );

    my $delete_fc = $dbh->prepare( "DELETE FROM feature_cvterm WHERE feature_cvterm_id = ?");
    $sth->execute();
    while( my $row = $sth->fetchrow_arrayref() ) {
        $delete_fc->execute( $row->[0] );
    }
    $sth->finish;
}
sub delete_tigr_roles {
    my ($tid) = @_;
    my $query = 
        "SELECT fc.feature_cvterm_id ".
        "FROM feature_cvterm fc, cvterm c, cv ".
        "WHERE fc.feature_id = $tid ".
        "AND fc.cvterm_id = c.cvterm_id ".
        "AND c.cv_id = cv.cv_id ".
        "AND cv.name ='TIGR_role'";
    my $sth = $dbh->prepare( $query );

    my $delete_fc = $dbh->prepare( "DELETE FROM feature_cvterm WHERE feature_cvterm_id = ?");
    $sth->execute();
    while( my $row = $sth->fetchrow_arrayref() ) {
        $delete_fc->execute( $row->[0] );
    }
    $sth->finish;
}
sub update_com_name {
    my ($fpid, $name) = @_;
    my $query = 
        "UPDATE featureprop SET value = ? ".
        "WHERE featureprop_id = ?";
    my $sth = $dbh->prepare( $query );
    eval {
        $sth->execute( $name, $fpid );
    };
    if( $@ ) {
        die("Couldn't execute query");
    }
}

sub change_com_name {
    my ($filter, @transcripts) = @_;
    my %retval;
    foreach my $transcript (@transcripts) {
        my ($tid, $fid, $com_name) = @{$transcript};
        my $orig = $com_name;

        my $prev = $com_name;
        do {
            $prev = $com_name;
            $com_name = &run_curate_subs( $com_name, $filter );
	    $com_name = &additional_rules($com_name) if (-e $rules_file);
        } while( $prev ne $com_name );

        unless( ($com_name eq $orig) && ($orig ne "hypothetical protein") && ($orig ne $chp) ) {
            $retval{$tid}->{'old'} = $orig;
            $retval{$tid}->{'new'} = $com_name;
            $retval{$tid}->{'featureprop_id'} = $fid;
        }

    }    
    %retval;
}

sub run_curate_subs {
    my ($com_name, $filter) = @_;
    ## cycle through all the register common name handlers
    if ($filter eq "gene_product_name") {
	foreach my $name ( keys %checks ) {
        	$com_name = $checks{$name}->($com_name);
        	die("Error running subroutine $name") unless( defined( $com_name ) );
    	}
    } elsif ($filter eq "gene_symbol") {
	foreach my $name ( keys %gene_sym_checks ) {
		$com_name = $gene_sym_checks{$name}->($com_name);
		die("Error running subroutine $name") unless( defined( $com_name ) );
	}
    }
    $com_name;
}
sub get_all_transcripts {
    my ($dbh, $filter) = @_;
    my $query = 
        "SELECT f.feature_id, fp.featureprop_id, fp.value ".
        "FROM feature f, cvterm c, featureprop fp, cvterm cfp ".
        "WHERE f.type_id = c.cvterm_id ".
        "AND c.name = 'transcript' ".
        "AND f.feature_id = fp.feature_id ".
        "AND fp.type_id = cfp.cvterm_id ".
        "AND cfp.name = '$filter'";
    print "$query\n";
    my $sth = $dbh->prepare( $query );
    $sth->execute;

    my $res = $sth->fetchall_arrayref();
    my @retval;
    map { push(@retval, $_) } @{$res};
    @retval;
}

sub register_common_name_checks {
    my %checks = (
                  'remove_trailing_parens' => \&remove_parens,
                  'remove_trailing_period' => \&remove_trailing_period,
                  'remove_protein_protein' => \&remove_protein_protein,
                  'remove_family_family' => \&remove_family_family,
                  'begins_with_orf' => \&begins_with_orf,
                  'find_putative_transposase' => \&find_putative_transposase,
                  'contains_homolog' => \&contains_homolog,
                  'similar_to' => \&similar_to,
                  'phospholipase_active_site' => \&phospholipase_active_site,
                  'contains_duf_upf' => \&contains_duf_upf,
                  'contains_golgi' => \&contains_golgi,
                  'contains_uncharacterised' => \&contains_uncharacterised,
                  'contains_underscore' => \&contains_underscore,
                  'bad_plurals' => \&contains_bad_plural,
                  'hypothetica' => \&hypothetica,
                  'americanise_some_words' => \&americanise_some_words,
                  'correct_putaive' => \&correct_putaive,
                  'correct_asparate' => \&correct_asparate,
                  'correct_unnamed' => \&correct_unnamed,
                  'only_the_word_protein' => \&only_protein,
                  'begins_with_crappy_word' => \&begins_with_ZB_locus,
                  'begins_with_residues' => \&begins_with_residues,
                  'short_bogus_name' => \&short_bogus_name,
                  'ttg_start' => \&ttg_start,
                  'gene_number_protein' => \&gene_number_protein,
                  'misc_name_changes' => \&misc_name_changes,
		  'remove_trailing_symbols' => \&remove_symbols,
		  'ends_in_binding' => \&ends_in_binding,
		  'remove_terminus' => \&remove_terminus,
		  'remove_colon' => \&remove_colon
                  );
    return %checks;
}

sub register_gene_symbol_checks {
	my %gene_sym_checks = (
			       'starts_with_orf' => \&starts_with_orf,
			       'contains_digits' => \&contains_digits
			      );
	return %gene_sym_checks;
}

######################################################################

sub misc_name_changes {
    my ($com_name) = @_;
    
    if( $com_name eq "transcriptional regulatory protein, C terminal family protein" ) {
        $com_name = 'putative transcriptional regulator';
    } elsif( $com_name eq "bacterial regulatory proteins, luxR family protein" ) {
        $com_name = "transcriptional regulator, LuxR family";
    } elsif( $com_name eq "bacterial regulatory helix-turn-helix proteins, AraC family protein" ) {
        $com_name = "transcriptional regulator, AraC family";
    } elsif( $com_name eq "bacterial regulatory proteins, tetR family protein" ) {
        $com_name = "transcriptional regulator, TetR family";
    } elsif( $com_name eq "bacterial regulatory proteins, gntR family protein" ) {
        $com_name = "transcriptional regulator, GntR family";
    } elsif( $com_name eq "bacterial regulatory proteins, lacI family protein" ) {
        $com_name = "transcriptional regulator, LacI family";
    } elsif( $com_name eq "FGGY family of carbohydrate kinases, C-terminal domain protein" ) {
        $com_name = "carbohydrate kinase, FGGY family";
    } elsif( $com_name eq "tripartite ATP-independent periplasmic transporters, DctQ component family protein" ) {
        $com_name = "tripartite ATP-independent periplasmic transporter, DctQ family";
    } elsif( $com_name eq "thiamin/thiamin pyrophosphate ABC transporter, thiamin/thiamin pyrophospate-binding protein" ) {
        $com_name = "thiamin/thiamine pyrophosphate ABC transporter, thiamin/thiamine pyrophospate-binding protein";
    } elsif( $com_name eq "GSPII_E N-terminal domain protein" ) {
        $com_name = "bacteriophage N4 adsorption protein B";
    } elsif( $com_name eq "transcriptional activator of defense systems" ) {
        $com_name = "multiple antibiotic resistance protein MarA";
    } elsif($com_name eq "type IV secretory pathway VirD2 components") {
        $com_name = "type IV secretory pathway protein";
    } elsif($com_name eq "hydro-lases, Fe-S type, tartrate/fumarate subfamily, beta region domain protein") {
        $com_name = "fumarate hydratase family protein";
    } elsif($com_name eq "glycogen/starch synthases, ADP-glucose type family protein") {
        $com_name = "glycogen/starch synthase";
    } elsif($com_name eq "glutamate synthases, NADH/NADPH, small subunit domain protein") {
        $com_name = "glutamate synthase, NADH/NADPH, small subunit";
    } elsif($com_name eq "K+ potassium transporter family protein") {
        $com_name = "potassium uptake protein";
    } elsif($com_name eq "domain related to MnhB subunit of Na+/H+ antiporter family protein") {
        $com_name = "Na+/H+ antiporter family protein";
    } elsif($com_name eq "arginine-tRNA-transferase, C terminus family protein") {
        $com_name = "putative arginine-tRNA-transferase";
    } elsif($com_name eq "cytochrome b(C-terminal)/b6/petD family protein") {
        $com_name = "cytochrome b family protein";
    } elsif($com_name eq "traG-like , N-terminal region family protein") {
        $com_name = "putative traG protein";
    } elsif($com_name eq "cyclic di-GMP binding protein VCA0042") {
        $com_name = "cyclic di-GMP binding protein";
    } elsif($com_name eq "alr5027 protein") {
        $com_name = "heme-binding protein HutZ";
    } elsif($com_name eq "putative 2-hydroxyacid dehydrogenase HI_1556") {
        $com_name = "putative 2-hydroxyacid dehydrogenase";
    } elsif($com_name eq "SULFATE TRANSPORTER SULFATE TRANSPORTER FAMILY PROTEIN") {
        $com_name = "sulfate permease family protein";
    } elsif($com_name eq "conserved protein with nucleoside triphosphate hydrolase domain") {
        $com_name = "putative ATP-dependent endonuclease";
    } elsif($com_name eq "gene 25-like lysozyme family protein") {
        $com_name = "lysozyme family protein";
    } elsif($com_name eq "bordetella uptake gene (bug) product family protein") {
        $com_name = "bug family protein";
    } elsif($com_name eq "phage/plasmid replication , gene II/X family protein") {
        $com_name = "phage/plasmid replication protein, gene II/X family";
    } elsif($com_name eq "invasion gene expression up-regulator, SirB family protein") {
        $com_name = "invasion gene expression up-regulator";
    } elsif($com_name eq "FAD linked oxidases, C-terminal domain protein") {
	$com_name = "FAD linked oxidase domain protein";
    } elsif(($com_name eq "PIII") || ($com_name eq "zn-dependent hydrolase of the beta-lactamase fold")) {
        $com_name = $chp;
    }
    
    $com_name;
}

sub gene_number_protein {
    my ($com_name) = @_;
    $com_name = $chp if( $com_name =~ /gene \d+ protein/ );
    $com_name;
}

sub ttg_start {
    my ($com_name) = @_;
    $com_name = $chp if( $com_name =~ /^ttg start/i );
    $com_name;
}
sub short_bogus_name {
    my ($com_name) = @_;
    $com_name = $chp if( $com_name =~ /^\w{1,2}\d{1,3}$/ );
    $com_name;
}
sub begins_with_residues {
    my ($com_name) = @_;
    $com_name = $chp if( $com_name =~ /^residues\b/i );
    $com_name;
}
sub begins_with_ZB_locus {
    my ($com_name) = @_;
    $com_name = $chp if( $com_name =~ /^\s*(Z|B)\d+\s+gene\s+product/ );
    $com_name;
}
sub only_protein {
    my ($com_name) = @_;
    $com_name = $chp if( $com_name =~ /^\s*protein\s*$/ );
    $com_name;
}
sub correct_unnamed {
    my ($com_name) = @_;
    if( $com_name =~ /\bunnamed\b/i ) {
        $com_name = $chp;
    }
    return $com_name;
}
sub correct_asparate {
    my ($com_name) = @_;
    if( $com_name =~ /\basparate\b/ ) {
        $com_name =~ s/asparate/aspartate/;
    }
    return $com_name;
}
sub correct_putaive {
    my ($com_name) = @_;
    $com_name =~ s/^(predicted|possible|potential|probable)/putative/i;
    $com_name =~ s/putaive/putative/gi;
    return $com_name;
}
sub americanise_some_words {
    #hah, americanise isn't americanized. ugh...
    my ($com_name) = @_;
    $com_name =~ s/utilisation/utilization/g;
    $com_name =~ s/utilising/utilizing/g;
    $com_name =~ s/dimerisation/dimerization/g;
    $com_name =~ s/disulphide/disulfide/g;
    $com_name =~ s/sulphur/sulfur/g;
    $com_name =~ s/mobilisation/mobilization/g; 
    if( $com_name =~ /\bhaem(\w+)/ ) {
        my $p = $1;
        $com_name =~ s/haem$p/hem$p/;
    }
    return $com_name;
    
}
sub hypothetica {
    my ($com_name) = @_;
    $com_name =~ s/hypotheical|hypothetic\b/hypothetical/g;
    if( $com_name =~ /conserved hypothetica\b/ ) {
        $com_name = $chp;
    } elsif( $com_name =~ /hypothetica\b/ )  {
        $com_name = 'hypothetical protein';
    }
    return $com_name;
}
sub contains_bad_plural {
    my ($com_name) = @_;
    $com_name =~ s/desulfurases/desulfurase/g;
    $com_name =~ s/synthases/synthase/g;
    $com_name =~ s/kinases/kinase/g;
    $com_name =~ s/decarboxylases/decarboxylase/g;
    $com_name =~ s/oxidases/oxidase/g;
    $com_name =~ s/\bgenes\b/gene/g;
    if( $com_name =~ /^(cons|no significant matches)$/i ) {
        $com_name = $chp;
    }
    return $com_name;
}
sub contains_underscore {
    my ($com_name) = @_;
    if( $com_name =~ /_/ ) {
        if( $com_name eq 'menC_gamma/gm+: o-succinylbenzoic acid (OSB) synthetase' ) {
            $com_name = 'o-succinylbenzoic acid (OSB) synthetase';
        }
    }
    return $com_name;
}
sub contains_uncharacterised {
    my ($com_name) = @_;
    if( $com_name =~ /\buncharacteri(z|s)ed\b/i ) {
        $com_name = $chp;
    }
    return $com_name;
}
sub contains_golgi {
    my ($com_name) = @_;
    if( $com_name =~ /golgi/i ) {
        $com_name = $chp;
    }
    return $com_name;
}
sub contains_duf_upf {
    my ($com_name) = @_;
    if( $com_name =~ /\b(DUF|UPF)/ ) {
        $com_name = $chp;
    }
    return $com_name;
}
sub phospholipase_active_site {
    my ($com_name) = @_;
    if( $com_name =~ /phospholipase d active site motif family protein/i ) {
        $com_name = 'phospholipase D family protein';
    }
    return $com_name;
}
sub similar_to {
    my ($com_name) = @_;
    if( $com_name =~ /similar to/ ) {
        $com_name = $chp;
    }
    return $com_name;
}
sub contains_homolog {
    my ($com_name) = @_;

    if( $com_name =~ /homolog/ ) {
        if( $com_name =~ /shiA homolog/ ) {
            $com_name = 'shiA protein';
        } elsif( $com_name =~ /virulence factor mviM homolog/ ) {
            $com_name = 'virulence factor mviM';
        } elsif( $com_name =~ /protein phnA homolog/ ) {
            $com_name = 'phnA protein';
        } elsif( $com_name =~ /protein seqA homolog/ ) {
            $com_name = 'seqA protein';
        } else {
            $com_name = $chp;
        }
    }
    return $com_name;
}
sub find_putative_transposase {
    my ($com_name) = @_;
    
    if( $com_name eq 'transposase and inactivated derivative' ) {
        $com_name = 'putative transposase';
    } elsif( $com_name =~ /^IS/ ) {
        $com_name = 'putative transposase';
    }
    return $com_name;
}
sub begins_with_orf {
    my ($com_name) = @_;
    if( $com_name =~ /^(orf)/i ) {
        $com_name = $chp;
    }
    return $com_name;
}

sub ends_in_binding {
    my ($com_name) = @_;
    if ($com_name =~ /(binding|like)\s*$/) {
	$com_name =~ s/^(protein)\s//;
	$com_name .= ' protein';  
    }
    return $com_name;
}

sub remove_terminus {
   my ($com_name) = @_;
   $com_name =~ s/\d*\s*C\-terminus//;
   $com_name =~ s/\d*\s*N\-terminus//;
   return $com_name;
}

sub remove_parens {
    my ($com_name) = @_;
    my $orig = $com_name;
    
    if( $com_name =~ /^(.*([\)\]]))\s*(domain protein|family protein|protein)?$/ ) {
        my ($name, $delim, $p) = ($1, $2, $3);
        my $reversed_name = reverse($name);
        $reversed_name =~ tr/\[\]\(\)/\]\[\)\(/;
        my ($extracted, $remainder) = extract_bracketed( $reversed_name, $delim );
        $remainder =~ tr/\[\]\(\)/\]\[\)\(/;
        $remainder =~ s/^\s+//;
        $com_name = reverse($remainder);
        $com_name .= " $p" if( $p );
    }
    $com_name = &remove_parens( $com_name ) unless( $orig eq $com_name );
    $com_name = $orig if( $com_name eq '' );
    return $com_name;
}
sub remove_trailing_period {
    my ($com_name) = @_;
    if( $com_name =~ /\.\s*(domain protein|family protein|protein)?$/ ) {
        my $dp = $1 if( $1 );
        $com_name =~ s/\s*\.\s*(domain protein|family protein|protein)?$//;
        $com_name .= " $dp" if( $dp );
    }
    return $com_name;
}
sub remove_protein_protein {
    my ($com_name) = @_;
    $com_name =~ s/protein\s+protein/protein/g;
    return $com_name;
}
sub remove_family_family {
    my ($com_name) = @_;
    $com_name =~ s/family\s+family/family/g;
    return $com_name;
}

sub remove_symbols {
    my ($com_name) = @_;
    if ( $com_name =~ /^(.+?)(\.|\,|\-|\_|\:|\/)$/ ) {
	$com_name = $1;
    }
    return $com_name;
}

# Gene product name contains gene symbol in the beginning. eg. rapA: rapA protein
# Remove the beginning gene symbol
sub remove_colon {
	my ($com_name) = @_;
	if ($com_name =~ /^(.+?):\s*(.*)/) {
		$com_name = $2;
	}
	return $com_name;
}

# If additional rules file is passed
sub additional_rules {
	my ($com_name) = @_;
	my $corrected_name = $com_name;
	if(exists($rules{$com_name})) {
		$corrected_name = $rules{$com_name};
	}
	return $corrected_name;
}

# Gene symbol should not start with "orf"
sub starts_with_orf {
    my ($gene_symbol) = @_;
    if ($gene_symbol =~ /^orf/i) {
	$gene_symbol = "";
    }
    return $gene_symbol;
}

# Gene symbol should not contain 4 consecutive digits
sub contains_digits {
    my ($gene_symbol) = @_;
    if ($gene_symbol =~ /\d{4}/) {
	$gene_symbol = "";
    }
    return $gene_symbol;
}

######################################################################
##                         Utility                                  ##
######################################################################
sub check_options {

   my $opts = shift;

   if( $opts->{'help'} ) {
       &_pod;
   }

   if( $opts->{'log'} ) {
       open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
   }

   foreach my $req ( qw(input_map_file username host) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }

   $no_change = 1 if( $opts->{'no_change'} );
   $debug = $opts->{'debug'} if( $opts->{'debug'} );

    #Assign password to be read from the file if it exists.
    if (defined ($options{'password_file'} && -s $options{'password_file'} ) ) {
	  open PASS, $options{'password_file'} or die ("Cannot open password file ". $options{'password_file'} ." : $!\n");
	  print STDERR ("Password from file will take priority over --password option\n") if (defined $options{'password'});
	  $password= <PASS>;
	  chomp $password;
	  close PASS;
    } elsif (defined ($options{'password'}) ){
	  $password = $options{'password'};
    } else {
	  die("Neither a password or a password file were supplied.  Please supply one or the other");
    }

    open MAP, $options{'input_map_file'} or die ("Cannot open input mapping file for reading: $!\n");
    my $line = <MAP>;
    chomp $line;
    my ($unused, $rest);
    ($database, $unused, $rules_file, $rest) = split(/,|\t/ , $line, 4);
    close MAP;


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
        die("Could not connect to database ".DBI->errstr);
    }
    $dbh->do("use $db");
    return $dbh;
}
sub _log {
   my ($level, $msg) = @_;
   if( $level <= $debug ) {
      print STDOUT "$msg\n";
      print $logfh "$msg\n" if( defined( $logfh ) );
      exit(1) if( $level == $ERROR );
   }
}
sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
######################################################################
##                             End                                  ##
######################################################################
