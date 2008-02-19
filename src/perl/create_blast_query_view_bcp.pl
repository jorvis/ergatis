#!/usr/bin/perl

=head1 NAME

create_blast_query_view_bcp.pl - reads blast results from a chado instance and creates a
bcp file to load a table (cm_blast) more appropriate for fast querying.

=head1 SYNOPSIS

USAGE: create_blast_query_view_bcp.pl 
            --database=strep
            --output_file=/path/to/somefile.bcp
            --analysis_id=8
          [ --first_id=1
            --user=access
            --pass=access
            --server=SYBTIGR
          ]

=head1 OPTIONS

B<--database,-d>
    Sybase database name to connect to.

B<--output_file,-o>
    Output BCP file that will be created.  See OUTPUT section for more details.

B<--analysis_id,-a>
    ID of a blast analysis to gather results for.  This should correspond to some
    analysis.analysis_id

B<--first_id,-f>
    Optional.  The first column of the output is an integer and can serve as the
    row ID.  Here you can specify the first value to use and all others will
    increment from it.  (default = 1)

B<--server,-s>
    Optional.  Sybase server to connect to (default = SYBTIGR).

B<--user,-u>
    Optional.  User account with select privileges on the specified database. (default = 'access')

B<--password,-p>
    Optional.  Password for user account specified. (default = 'access')

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

Chado is meant to be used as a flexible, extensible database schema and was not designed
with query speed in mind.  To query all BLAST results, including properties of all
sequences involved and HSP scoring statistics, requires fairly large queries that
join across multiple tables.  These can be prohibitively slow, especially when driving
the display of a web app where users can only be expected to wait a few seconds for a
page load.

A collection of materialized views, designed for fast querying, can solve many of these
problems.  Because Sybase doesn't support materialized views we have to generate them
on the fly.  This script generates the BCP file needed to load a materialized view
for fast BLAST querying.

=head1  INPUT

BLAST storage convention expected should be entered here.

=head1  OUTPUT

The output is written in BCP format to the file defined by the --output_file option.  It
is made up of these columns:

    1. unique row ID (auto-generated)
    2. query feature_id
    3. query organism_id
    4. hit feature_id
    5. hit organism_id
    6. percent identity
    7. percent similarity
    8. p value
    9. match feature_id

Columns 1-5 are straight-forward.  Columns 5 and 6 are averages of the percent identity and 
percent similarity, respectively, of each of the HSPs of the match.  P value is the lowest
(best) p value of the HSPs of the match.

The columns are delimited by a null-byte character followed by a tab, and each line ends with 
null-byte character and a new line (\n).

=head1 LOADING

Here is an example command to load this bcp file into a sybase instance:

    /usr/local/devel/ANNOTATION/ard/current/bin/flatFileToChado 
        --username=chado_admin 
        --password=chado_admin99 
        --database=entamoeba 
        --server=SYBTIGR 
        --database_type=sybase 
        --bcp_ext=bcp 
        --bcpmode=in 
        --batchsize=3000 
        --infile=cm_blast.bcp 
        --debug_level=0 
        --ignore_empty_bcp=1 
        --logfile=cm_blast.load.log

=head1  CONTACT

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
                          'database|d=s',
                          'analysis_id|a=s',
                          'first_id|f=s',
                          'user|u=s',
                          'password|p=s',
                          'server|s=s',
                          'output_file|o=s',
                          'log|l=s',
                          'help|h') || pod2usage();

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

_log("connecting to the database ($options{server}.$options{database})");
my $dbh = DBI->connect("dbi:Sybase:server=$options{server}; packetSize=8092", $options{user}, $options{password}, {PrintError=>1, RaiseError=>1});
$dbh->do("use $options{database}");

## reused in multiple preps
my $qry;

## open the output file
open (my $ofh, ">$options{output_file}") || die "can't create output file: $!";

## gather some terms
my %terms = (
    computed_by => get_cvtermid('computed_by'),
    input_of => get_cvtermid('input_of'),
    match => get_cvtermid('match'),
    match_part => get_cvtermid('match_part'),
    percent_similarity => get_cvtermid('percent_similarity'),
    p_value => get_cvtermid('p_value'),
);

## get all match features involved in the passed analysis
$qry = qq{
    SELECT af.analysisfeature_id, query.organism_id query_organism_id, hit.organism_id hit_organism_id,
           match.feature_id match_feature_id, query.feature_id query_feature_id, hit.feature_id hit_feature_id, 
           af.pidentity, af.significance
      FROM feature match, feature query, feature hit, featureloc q2mfl, featureloc h2mfl, analysisfeature af
     WHERE match.feature_id=af.feature_id
       AND q2mfl.srcfeature_id=query.feature_id
       AND q2mfl.feature_id=match.feature_id
       AND q2mfl.rank=0
       AND h2mfl.srcfeature_id=hit.feature_id
       AND h2mfl.feature_id=match.feature_id
       AND h2mfl.rank=1
       AND af.analysis_id = ?
       AND match.type_id = ?
};

my $match_selector = $dbh->prepare($qry);
   $match_selector->execute( $options{analysis_id}, $terms{match} );

## the pval.rank here should be unnecessary, but there was a BSML encoding bug that 
##  entered two featureprops for each p_value
$qry = qq{
    SELECT match_part.feature_id, percsim.featureprop_id, percsim.value percent_similarity, 
           pval.featureprop_id, pval.value p_value, af.pidentity percent_identity
      FROM feature match_part, feature match, analysisfeature af, feature_relationship fr,
           featureprop percsim, featureprop pval
     WHERE match_part.feature_id=af.feature_id
       AND fr.subject_id=match_part.feature_id
       AND fr.object_id=match.feature_id
       AND percsim.feature_id=match_part.feature_id
       AND percsim.type_id=?
       AND pval.feature_id=match_part.feature_id
       AND pval.type_id=?
       AND pval.rank=0
       AND af.analysis_id=?
       AND match_part.type_id=?
       AND match.feature_id=?
};
## execute with $match_part_selector->execute( percsim_type_id, pval_type_id, analysis_id, match_part_type_id, specific_match_feature_id );
my $match_part_selector = $dbh->prepare($qry);

my @cols;
my $next_row_id = $options{first_id};

while ( my $match = $match_selector->fetchrow_hashref ) {
    ## reset for this row
    @cols = (
                $$match{query_feature_id}, 
                $$match{query_organism_id}, 
                $$match{hit_feature_id}, 
                $$match{hit_organism_id}, 
                '', 
                '', 
                '', 
                $$match{match_feature_id},
            );
    
    $match_part_selector->execute( $terms{percent_similarity}, $terms{p_value}, $options{analysis_id}, 
                                   $terms{match_part}, $$match{match_feature_id} );
    
    my $min_p_value = 10000;  ## fake
    my $p_ident_sum = 0;
    my $p_ident_count = 0;
    my $p_sim_sum = 0;
    my $p_sim_count = 0;
    
    while ( my $match_part = $match_part_selector->fetchrow_hashref ) {
        $p_ident_sum += $$match_part{percent_identity};
        $p_ident_count++;
        $p_sim_sum = $$match_part{percent_similarity};
        $p_sim_count++;
        
        if ( $$match_part{p_value} < $min_p_value ) {
            $min_p_value = $$match_part{p_value};
        }
    }
    
    $cols[4] = $p_ident_sum / $p_ident_count;
    $cols[5] = $p_sim_sum / $p_sim_count;
    $cols[6] = $min_p_value;
    
    print $ofh "$next_row_id\0\t", join("\0\t", @cols), "\0\n";
    $next_row_id++;
} 

$match_part_selector->finish();
$match_selector->finish();
$dbh->disconnect();


exit(0);

sub get_cvtermid {
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

    _log("got cvterm_id $cvterm_id for name $name");
    return $cvterm_id;
    
}


sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
}


sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( database output_file analysis_id );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
    
    ##
    ## you can do other things here, such as checking that files exist, etc.
    ##
    
    ## handle some defaults
    $options{user} = 'access' unless ($options{user});
    $options{password} = 'access' unless ($options{password});
    $options{server} = 'SYBTIGR' unless ($options{server});
    $options{first_id} = 1 unless ($options{first_id});
}
