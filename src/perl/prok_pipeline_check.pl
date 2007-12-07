#!/usr/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

data_consistency_check.pl - Checks data produced in the prokaryotic pipeline

=head1 SYNOPSIS

USAGE: data_consistency_check.pl
            --pipeline_id=9923
            --repository_root=/usr/local/annotation/MOORE
          [ --log=/path/to/file.log
            --debug=4
            --help
          ]

=head1 OPTIONS

B<--pipeline_id,-p>
    REQUIRED. The pipeline id to review data

B<--repository_root,-r>
    REQUIRED. Location of the pipeline ( Ergatis file structure )

B<--log,-l>
    REQUIRED. Logfile.

B<--debug,-d>
    OPTIONAL. Larger number is more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

    This script will review the data produced by a normal prokaryotic annotation pipeline run.
    Specifics for fail conditions will be added here.

=head1  INPUT

    The data ( and where to find it ) requirements will be put here.


=head1 OUTPUT

    Output is in the form of a log file.  

=head1  CONTACT

    Kevin Galens
    kgalens@jcvi.org

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Ergatis::Logger;
use Prokaryotic::Pipeline::Statistics;
use Data::Dumper;
use DBI;

####### GLOBALS AND CONSTANTS ###########
my $prok;
my $pipeline_id;
my $repository_root;
my $transcript_level = 0;
#########################################

my %options = ();
my $results = GetOptions (\%options, 
                          'pipeline_id|p=s',
                          'repository_root|r=s',
                          'include_transcript_level|t',
                          'log|l=s',
                          'debug=s',
                          'help|h') || &_pod;

#Setup the logger
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

&check_parameters( \%options );


print "\nSTART: ".localtime(time)."\n\n";

#Get the genome information
my $organism;
eval {$organism = $prok->get_genome_name;};
unless( $@ ) {
    print "\tOrganism name: $organism\n";
} else {
    &_warn("\tOrganism name: Unknown. Check pipeline_summary.default component configuration ($@)");
}

my $abbreviation;
eval { $abbreviation = $prok->get_genome_abbreviation; };
unless( $@ ) {
    print "\tOrganism abbreviation: $abbreviation\n";
} else {
    &_warn("\tOrganims abbreviation: Unknown.  Check pipeline_summary.default component configuration ($@)");
}

#See if it exists in the common database
my @exists_in_common = &exists_in_common( $organism );
if( @exists_in_common == 0 ) {
    &_warn( "\tDid not find the organism in common..genomes table" );
} elsif( @exists_in_common > 1 ) {
    &_warn("\tFound multiple entries in common..genomes table for this organism");
    foreach my $result ( @exists_in_common ) {
        print "\t\tname: ".$result->[0]."\tstate: ".$result->[1]."\n";
    }
} else {
    print "\tFound the organism in the common..genomes table. Stage: ".$exists_in_common[0]->[1]."\n";
}

my $genome_length;
eval { $genome_length = $prok->get_genome_length; };
if( $@ ) {
    &_warn("\tGenome length: Could not get genome length ($@)");
} else {
    print "\tGenome length: $genome_length residues\n";
}

#Assembly Info
print "\nASSEMBLY INFO:\n";
my $asmbls;
eval{ $asmbls = $prok->get_asmbl_info; };
&_warn( $@ ) if( $@ );

my $gene_count;
eval { $gene_count = $prok->get_gene_count; };
&_warn( $@ ) if( $@ );

my @sorted_keys = sort by_assembly_num keys %{$asmbls};

foreach my $key ( @sorted_keys ) {
    
    print "\t$key:\n\t\t".($asmbls->{$key})." residues, ".
        $prok->{'gene_count_assemblies'}->{$key}." ORFs\n";
}

#ORF Info
print "\nORF INFO:\n";
my $min = (($genome_length/1000) * 0.8) - 20;
my $max = (($genome_length/1000) * 1.2) + 20;
print "\torf count: $gene_count\n";

if( $gene_count < $max && $gene_count > $min ) {
    print "\torf count in acceptable range [ ".int($min)." - ".int($max)." ]\n";
} else {
    die("\torf count is not in acceptable range.  Check gene calls\n\t".
        "Expected range: [".int($min)." - ".int($max)."]" );
    
}

print "\n\tEvidence (BER):\n";
my $ber_results_count;
eval { $ber_results_count = $prok->get_ber_results_count;};
&_warn("\t$@") if( $@ );
print "\tNon-zero BER results: $ber_results_count\n";
die("Could not find any non-zero sized ber results") unless( $ber_results_count );

print "\n\tEvidence (HMM):\n";
my $hmm_results_count;
eval{ $hmm_results_count = $prok->get_hmm_results_count;};
&_warn("\t$@") if( $@ );
print "\tNon-zero HMM results: $hmm_results_count\n";
die("Could not find any non-zero sized hmm results") unless( $hmm_results_count );

print "\ntRNA:\n";
my $tRNA_count;
eval { $tRNA_count = $prok->get_rna_count };
&_warn("\t$@") if( $@ );
print "\ttRNA predictions: $tRNA_count\n";
&_warn("tRNA-scan should be run on this genome") unless( $tRNA_count );


if( $transcript_level ) {
    print "\nTranscript Level Computes\n";

    #The transcript level computes.  Set up a data structure to do processing in a loop
    my @transcript_computes = 
        (
         [ 'lipoprotein_motif', 'default', 'bsml' ],
         [ 'tmhmm', 'default', 'raw' ],
         [ 'targetp', 'default', 'raw' ],
         [ 'signalp', 'default', 'raw' ],
         [ 'ps_scan', 'default', 'bsml' ],
         [ 'wu-blastp', 'COGS', 'bsml' ] 
         );

    foreach my $compute( @transcript_computes ) {
        my $count;
        eval {
            $count = $prok->_count_evidence_files( $compute->[0],
                                                   $compute->[1],
                                                   $compute->[2] );
        };
        &_warn( "\t$@" ) if( $@ );
        print "\t".$compute->[0].": $count\n";
        &_warn( "\tThe number of ".$compute->[0]." output files should equal the".
                " gene count") unless( $count == $gene_count );
    }
}

print "\nDone ".localtime()."\n";
print "\n";
#################### SUB ROUTINES ######################
sub by_assembly_num {

    my $a_num = $1 if( $a =~ /assembly\.(\d+)/ );
    my $b_num = $1 if( $b =~ /assembly\.(\d+)/ );

    return $a_num <=> $b_num;
}

sub check_parameters {
    my ($opts) = @_;

    &_pod if( $opts->{'help'} );

    #Pipeline Id is required
    if( $opts->{'pipeline_id'} ) {
        $pipeline_id = $opts->{'pipeline_id'};
    } else {
        $logger->logdie("option --pipeline_id is required");
    }

    #repository root is required
    if( $opts->{'repository_root'} ) {
        $repository_root = $opts->{'repository_root'};
    } else {
        $logger->logdie("option --repository_root is required");
    }

    if( $opts->{'include_transcript_level'} ) {
        $transcript_level = 1;
    }

    $prok = new Prokaryotic::Pipeline::Statistics( {
        'pipeline_id' => $pipeline_id,
        'repository_root' => $repository_root,
    });

}

sub connect_to_db {
    my ($user, $pass) = @_;
    $user = $pass = "access" unless( $user && $pass );

    my $retval = DBI->connect('dbi:Sybase:server=SYBTIGR; packetSize=8092', $user, $pass)
        or $logger->logdie("Unable to connect to database".${DBI->errstr});

    return $retval;
}

sub exists_in_common {
    my ($organism) = @_;

    my $dbh = &connect_to_db;
    
    my $query = 
        "SELECT name, stage ".
        "FROM common..genomes ".
        "WHERE name = \"$organism\"";


    my $sth = $dbh->prepare( $query );
    $sth->execute;

    my @results;

    while( my $ar = $sth->fetchrow_arrayref() ) {
        push( @results, $ar );
    }

    
    $sth->finish;
    $dbh->disconnect;

    return @results;
    
}

sub _pod {   
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

sub _warn {
    my ($msg) = @_;
    print "$msg\n";
    $logger->warn( $msg );
}
