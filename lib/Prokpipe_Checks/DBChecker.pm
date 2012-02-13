#!/usr/bin/perl -w

#DBChecker.pm

# This module is designed to run consistency checks on the database created after
# Chado loading in the Prokaryotic Annotation Pipeline.  If a check fails,
# either a warning will be written to the Ergatis prokpipe_consistency_checks
# component log or an error will be written and kill off the pipeline.

package Prokpipe_Checks::DBChecker;

use Ergatis::Logger;

my $query = "";
my @queries;
our $dbh;
my $row;

#  Add more subroutines for checking things as time passes

# This test checks for GO role links that do not have GO evidence listed
sub no_GO_evidence {

	$query = "SELECT f.uniquename, d.accession, fcp.value
    		FROM feature f
        	JOIN feature_cvterm fc ON f.feature_id = fc.feature_id
        	JOIN cvterm c ON fc.cvterm_id = c.cvterm_id
        	JOIN cv ON c.cv_id = cv.cv_id
        	JOIN dbxref d ON c.dbxref_id = d.dbxref_id
        	JOIN feature_cvtermprop fcp ON fc.feature_cvterm_id = fcp.feature_cvterm_id
        	JOIN cvterm c2 ON fcp.type_id = c2.cvterm_id
        	JOIN cvtermsynonym cs ON c2.cvterm_id = cs.cvterm_id
		WHERE NOT EXISTS (SELECT c2.name cs.synonym)";

   	$queries->{'no_go_evidence'} = $dbh->prepare( $query );
	$queries{'no_go_evidence'}->execute();
	while ($row = $queries{'no_go_evidence'}->fetchrow_arrayref() ) { 
		return 1;
	}
}

# This test checks for no gene symbol or EC number in conserved hypothetical, plain hypothetical, or conserved domain proteins
sub no_GS_EC {

	$query = "SELECT i.feat_name, i.com_name, isnull(i.gene_sym, \"NULL\"), isnull(i.ec#, \"NULL\"), r.role_id " .
	"FROM ident i, role_link r " .
	    "WHERE i.feat_name = r.feat_name " .
		"AND r.role_id in (156, 270, 704) " .
		    "and (i.ec# is not null or i.gene_sym is not null) " .
			"and i.gene_sym !=\"\" " .
			    "and i.complete is not null";

	$query = 'SELECT f.uniquename 
		FROM feature f';

}

1;
