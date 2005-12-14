#!/usr/local/bin/perl -w

=head1  NAME 

iprscan2legacydb.pl - load iprscan data into the legacy database schema.

=head1 SYNOPSIS

USAGE: iprscan2bsml.pl 
        --input_list=/path/to/somefile.raw.list
        --database=sma1
      [ --log=/path/to/some/log ]

=head1 OPTIONS

B<--input_list,-i> 
    Raw list file from an iprscan workflow component run.

B<--database,-d> 
    Sybase project database ID.

B<--log,-d> 
    optional.  will create a log file with summaries of all actions performed.

B<--help,-h> 
    This help message/documentation.

=head1   DESCRIPTION

This script is used to load iprscan data into the legacy database schema.  For
each model loaded, the script first checks the database to see if there are
any existing matches marked as "curated".  If the same model -> method -> accession
combination is also found among the new data it will be marked as curated.

As the evidence for each model is loaded all previous iprscan evidence is removed
for that model.

=head1 INPUT

The input is the raw, tab-delimited output from an iprscan component run.  For example:

    bma1.model.14040_00266  B679F3E794F518D5        810     HMMPfam PF00058.7       Low-density lipoprotein receptor repeat 675     715        0.017   T       03-Nov-2005     IPR000033       Low-density lipoprotein receptor, YWTD repeat   Cellular Component: membrane(GO:0016020)
    bma1.model.14040_00266  B679F3E794F518D5        810     HMMSmart        SM00192 no description  1       36      2.5e-07 T       03-Nov-2005  IPR002172       Low density lipoprotein-receptor, class A
    bma1.model.14040_00266  B679F3E794F518D5        810     FPrintScan      PR00261 LDLRECEPTOR     8       29      6.9e-48 T       03-Nov-2005  IPR002172       Low density lipoprotein-receptor, class A
    bma1.model.14040_00266  B679F3E794F518D5        810     ScanRegExp      PS00010 ASX_HYDROXYL    393     404     8e-5    T       03-Nov-2005  IPR000152       Aspartic acid and asparagine hydroxylation site
    bma1.model.14040_00266  B679F3E794F518D5        810     ProfileScan     PS50068 LDLRA_2 40      77      14.250  T       03-Nov-2005        IPR002172       Low density lipoprotein-receptor, class A
    bma1.model.14040_00266  B679F3E794F518D5        810     superfamily     SSF57424        LDL receptor-like module        162     201        4e-13   T       03-Nov-2005     NULL    NULL
    bma1.model.14040_00266  B679F3E794F518D5        810     Seg     seg     seg     17      25      NA      ?       03-Nov-2005     NULL        NULL

=head1 OUTPUT

    there is no output unless you use the --log option.

=head1 CONTACT

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
			  'input_list|i=s',
              'database|d=s',
              'log|l=s',
			  'help|h') || pod2usage();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

## make sure all passed options are peachy
&check_parameters(\%options);

## get the current user
my ($user) = getpwuid($<);

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">$options{log}") || die "can't create log file: $!";
}

## create the database connection
my $dbh = DBI->connect("dbi:Sybase:server=SYBTIGR; packetSize=8092", 'egc', 'egcpwd', {RaiseError => 1, PrintError => 1});
my $result = $dbh->do("use $options{database}");

##
## prepare all the SQL statements to be used later.
##
my $qry = 'INSERT INTO evidence (feat_name, ev_type, accession,  end5, end3,  ' .
                               '  curated, m_lend, m_rend, rel_end5, rel_end3, date, assignby, ' .
                               ' method, change_log, save_history, score, total_score) ' .
          "VALUES (?, ?, ?, ?, ?, ?, -1, -1, 0, 0, GETDATE(), '$user', 'workflow', 0, 0, ?, ?)";
my $insert = $dbh->prepare($qry);

$qry = 'SELECT ev_type, accession ' .
       'FROM evidence ' .
       "WHERE ev_type in ('BlastProDom', 'Coil', 'FPrintScan', 'HMMPfam', 'HMMpir', " .
                         "'HMMSmart', 'HMMTigr', 'ProfileScan', 'ScanPrositeC', 'ScanRegExp', 'Seg', 'superfamily') " .
       '  AND feat_name = ? ' .
       '  AND curated   = 1';
my $curated_check = $dbh->prepare($qry);

$qry = 'DELETE FROM evidence ' .
       "WHERE ev_type in ('BlastProDom', 'Coil', 'FPrintScan', 'HMMPfam', 'HMMpir', " .
                         "'HMMSmart', 'HMMTigr', 'ProfileScan', 'ScanPrositeC', 'ScanRegExp', 'Seg', 'superfamily') " .
       ' and feat_name = ?';
my $delete = $dbh->prepare($qry);

open (my $listfh, $options{input_list}) || die "can't read list file\n";

my %curated;  ## gets overwritten for each list file.

while (<$listfh>) {
    chomp;
    next if (/^\s*$/);

    _log("processing file $_");

    ## extract the feat_name from the file name:
    ##  bma1.model.14990_07908.iprscan.raw -> 14990.m07908
    my $feat_name;
    if ( /(\d+)_(\d+)\.iprscan\.raw/ ) {
        $feat_name = "$1.m$2";
    } else {
        die "improperly named input file.  couldn't extract feat name from $_";
    }

    ## query each program and remember the ev_type->accession relationships that have curated set to 1.
    ## if these are encountered, they will be continue to be marked as curated.
    undef %curated;
    
    $curated_check->execute($feat_name);
    while (my $res = $curated_check->fetchrow_hashref) {
        $curated{ $res->{ev_type} }{ $res->{accession} } = 1;
    }

    ## now delete the existing interpro evidence for this feat_name
    _log("deleting all interpro evidence for $feat_name");
    $delete->execute($feat_name);

    ## now we load the new predictions.

    ## open this file
    open (my $resfh, "<$_") || die "can't read file $_ : $!";

    while (<$resfh>) {
        chomp;
        next if (/^\s*$/);
        my @cols = split("\t", $_);

        if (scalar(@cols) > 12) {
            ## column 8 (eval) can be null
            $cols[8] = 'NULL' if ($cols[8] eq 'NA');
            $cols[3] = 'HMMpir' if ($cols[3] eq 'HMMPIR');
            
            ## two possible naming conventions
            ## reformat the feat_name from pva1.model.2817_00001 to 2817.m00001
            if ($cols[0] =~ /.+\..+\.(\d+)_(\d+)/) {
                $cols[0] = "$1.m$2";
                store_results_for_model( $cols[0], $cols[3],  $cols[4], $cols[6], $cols[7], $cols[8] );

            ## reformat the feat_name from bba1.model.1806356 to 18.m06356
            ## DEPRECATED for legacydb
            } elsif ($cols[0] =~ /.+\..+\.(\d+)(\d{5})/) {
                $cols[0] = "$1.m$2";
                store_results_for_model( $cols[0], $cols[3],  $cols[4], $cols[6], $cols[7], $cols[8] );
            } else {
                print STDERR "skipping unrecognized id $cols[0]\n";
            }
        } else {
            unless (/FAILED/ || /failed and exit/) {
                print STDERR "skipping the following line because of incorrect column count:\n$_\n";
            }
        }
    }

}

$insert->finish();
$curated_check->finish();
$dbh->disconnect();

sub store_results_for_model {
    ## this feat_name should be named properly, like 2817.m00001
    my ($feat_name, $ev_type, $accession, $end5, $end3, $score) = @_;
    
    if ( exists $curated{ $ev_type }{ $accession } ) {
        _log("inserting ( $feat_name, $ev_type, $accession, $end5, $end3, 1, $score, $score)");
        $insert->execute( $feat_name, $ev_type, $accession, $end5, $end3, 1, $score, $score);
    } else {
        _log("inserting ( $feat_name, $ev_type, $accession, $end5, $end3, 0, $score, $score)");
        $insert->execute( $feat_name, $ev_type, $accession, $end5, $end3, 0, $score, $score);    
    }
}

sub _log {
    my $msg = shift;
    
    print $logfh "$msg\n" if $logfh;
}

sub check_parameters {
    
    ## database and input_list are required
    unless ( defined $options{database} && $options{input_list} ) {
        print STDERR "database and input_list options are required\n\n";
        pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
    }
    
    ## make sure input list exists
    unless ( -e $options{input_list} ) {
        print STDERR "\n\ninput_list $options{input_list} does not exist\n\n";
        exit(1);
    }
}
