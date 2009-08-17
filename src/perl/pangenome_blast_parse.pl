#!/usr/bin/perl

=head1  NAME

pangenome_blast_parse.pl - Parses BLAST BSML output files and extracts data needed for pangenome analysis.

=head1 SYNOPSIS

USAGE: pangenome_query_list.pl
        --input=/path/to/somefile.bsml
        --output=/path/to/output.dat
      [ --log=/path/to/some/log ]

=head1 OPTIONS

B<--input_bsml,-i>
    BSML file containing BLAST results to parse. 

B<--output_path,-o>
    Desired prefix of output files.

B<--log,-d>
    optional. Will create a log file with summaries of all actions performed.

B<--debug>
    optional. Will display verbose status messages to STDERR if logging is disabled.

B<--help,-h>
    This help message/documentation.

=head1   DESCRIPTION

    This script parses the BSML output from WU-BLASTP and WU-TBLASTN runs for pangenome analysis
    and serializes the prepared data so that it can be later unserialized for the analysis step.

=head1 INPUT

    The input should be one of the set of BSML files containing BLAST results from a pangenome
    pipeline run.

=head1 OUTPUT

    There is no output unless you use the --log option.

=head1 CONTACT

    Brett Whitty
    bwhitty@tigr.org

=cut

use File::Basename;
use Pod::Usage;
use Storable qw(nstore retrieve);
use XML::Twig;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use strict;

my %options = ();
my $results = GetOptions (  \%options,
                            'input|i=s',
                            'output_path|o=s',
#                           'filter|f=s',
                            'organism_to_db_mapping|odb:s',
                            'db_list|dl:s',
                            'debug|d=s',
                            'command_id=s',       ## passed by workflow
                            'logconf=s',          ## passed by workflow (not used)
                            'log|l=s',
                            'help|h',
                         ) || pod2usage();

my @results = ();
my $dups_temp = {};
my %dups = ();
my $db_to_org = {};
my $db_filter = undef;
my $dbs = {};

my $twig = XML::Twig->new(
                            twig_roots  => { 
                                                'Seq-pair-alignment' => \&processSeqPairAlignment
                                           }
                         );

if($options{'db_list'}) {
    &read_db_list();
}
if($options{'organism_to_db_mapping'}) {
    &read_organism_to_db_mapping();
}

if (!-d $options{'output_path'}) {
    die "must specify output path as a directory";
}

$options{'output_path'} =~ s/\/$//;

my $input_prefix = basename($options{'input'},".bsml");

my $ifh;

if (-e $options{'input'}.".gz" ) {
    $options{'input'} .= ".gz";
} elsif (-e $options{'input'}.".gzip") {
    $options{'input'} .= ".gzip";
}
if ($options{'input'} =~ /\.(gz|gzip)$/) {
    open ($ifh, "<:gzip", $options{'input'}) || die "couldn't open '$options{input}' for reading: $!";
} else {
    open ($ifh, "<".$options{'input'}) || die "couldn't open '$options{input}' for reading: $!";
}

$twig->parse($ifh);

close $ifh;


foreach my $key (keys %$dups_temp) {
    if (scalar(@{$dups_temp->{$key}}) > 1) {
        $dups{$key} = join(" ", @{$dups_temp->{$key}});
    }
}

nstore([\@results,\%dups], $options{'output_path'}."/".$input_prefix.".blast.stored") || die "couldn't serialize results";

exit(0);

sub processSeqPairAlignment { 
    my ($twig, $feat) = @_;
    my ($p_sim, $p_ident, $p_cov_ref, $p_cov_comp, $this_chain);
    
    my $query_id   = $feat->{'att'}->{'refseq'};
    my $subject_id = $feat->{'att'}->{'compseq'};
    print STDERR "$query_id $subject_id\n";
    my $len_total = 0;
    my $pident_sum = 0;
    my $psim_sum = 0;
    my $all_seg_p_ident = 0;
    my $all_seg_p_sim = 0;

    my $parent_chain = 0;
    
    for my $att ($feat->children('Attribute')) {
        if($att->{att}->{'name'} eq 'percent_coverage_compseq') {
            $p_cov_comp = $att->{att}->{'content'};
        }
        elsif($att->{att}->{'name'} eq 'percent_coverage_refseq') {
            $p_cov_ref = $att->{att}->{'content'};
        }
        elsif($att->{att}->{'name'} eq 'percent_similarity') {
            $all_seg_p_sim = $att->{att}->{'content'};
        }
        elsif($att->{att}->{'name'} eq 'percent_identity') {
            $all_seg_p_ident = $att->{att}->{'content'};
        }        
    }
#print STDERR "$p_cov_comp $p_cov_ref $all_seg_p_sim $all_seg_p_ident\n";

#    for my $seqpairrun ( $feat->children('Seq-pair-run') ) {
#        my $run_length = $seqpairrun->{att}->{runlength};
#        $len_total += $run_length;
#        for my $attribute ( $seqpairrun->children('Attribute') ) {
#            if ($attribute->{att}->{'name'} eq 'percent_identity') {
#                $p_ident = $attribute->{att}->{'content'};
#            } elsif ($attribute->{att}->{'name'} eq 'percent_similarity') {
#                $p_sim = $attribute->{att}->{'content'};
#            } elsif ($attribute->{att}->{'name'} eq 'percent_coverage_refseq') {
#                #$p_cov_ref = $attribute->{att}->{'content'};
#            } elsif ($attribute->{att}->{'name'} eq 'percent_coverage_compseq') {
#                #$p_cov_comp = $attribute->{att}->{'content'};
#            }
#        }
#        $pident_sum += (($p_ident/100.0) * $run_length); 
#        $psim_sum += (($p_sim/100.0) * $run_length); 
#    }
#    my $all_seg_p_ident = sprintf("%.1f", $pident_sum / $len_total * 100);
#    my $all_seg_p_sim = sprintf("%.1f", $psim_sum / $len_total * 100);

    # $qprot may be an assembly for TBLASTN, not a protein
    my($qdb, $qprot) = ($2, $1);
    if ($query_id =~ /^(([^\.]+)\..*)$/) {
	($qdb, $qprot) = ($2, $1);
    }
    # handle BSML sequences (for NC GenBank entries) created via genbank2bsml
    elsif ($query_id =~ /^ref_(NC_\d+)_(\S+)$/) {
	($qdb, $qprot) = ($2, $1);
    } else {
	die "couldn't parse out db and sequence id from query\n$query_id\n";
    }

    # $sprot may be an assembly for TBLASTN, not a protein
    my($sdb, $sprot);
    if ($subject_id =~ /^(([^\.]+)\..*)$/) {
        ($sdb, $sprot) = ($2, $1);
    } 
    # handle BSML sequences (for NC GenBank entries) created via genbank2bsml
    elsif ($subject_id =~ /^ref_(NC_\d+)_(\S+)$/) {
	($sdb, $sprot) = ($2, $1);
    } else {
	die "couldn't parse out db and sequence id from subject\n$subject_id\n";
    }
    
    if(!$options{'db_list'} && !$options{'organism_to_db_mapping'}) {
        $db_to_org->{$qdb} = $qdb;
        $db_to_org->{$sdb} = $sdb;
    }
    
    if(!$options{'db_list'} || ($dbs->{$qdb} && $dbs->{$sdb})) {
#    if ($skip_filter || ($filter->{$qdb}->{$qprot} && $filter->{$sdb}->{$sprot})) {

#        if (!$db_filter || ($db_filter->{$db_to_org->{$qdb}}->{$qprot} && $db_filter->{$db_to_org->{$sdb}}->{$sprot})) {

        if (!$db_filter || ($db_filter->{$qdb} && $db_filter->{$sdb})) {

            #   if ($all_seg_p_ident >= 50.0 && ($p_cov_ref >= 50.0 || $p_cov_comp >= 50.0)) {  ## switched from p_ident to p_sim
            #       push (@results, [$qdb,$qprot,$sdb,$sprot,$all_seg_p_ident,$p_cov_ref]);
            if ($all_seg_p_sim >= 50.0 && ($p_cov_ref >= 50.0 || $p_cov_comp >= 50.0)) {
                #print STDERR "$db_to_org->{$qdb},$qprot,$db_to_org->{$sdb},$sprot,$all_seg_p_sim,$p_cov_ref\n";
                push (@results, [$db_to_org->{$qdb},$qprot,$db_to_org->{$sdb},$sprot,$all_seg_p_sim,$p_cov_ref]);
            }
            else {
                print STDERR "filtering 1 $qprot $sprot $p_cov_ref $p_cov_comp $all_seg_p_sim\n";
            }
            #   if ($qdb eq $sdb && $qprot ne $sprot && $all_seg_p_indent == 100 && $p_cov_ref == 100 && $p_cov_comp == 100.0) {
            #   I won't exclude the self hit here because we do want the query id in the group of dups
            if ($db_to_org->{$qdb} eq $db_to_org->{$sdb} && $all_seg_p_ident == 100 && $p_cov_ref == 100 && $p_cov_comp == 100) {
            push(@{$dups_temp->{$db_to_org->{$sdb}}}, $sprot);
            print $sprot."\n";
            print $p_cov_ref." ".$p_cov_comp."\n";
            }
#            else {
#                print STDERR "filtering 2 $qprot $sprot $p_cov_ref $p_cov_comp $all_seg_p_sim\n";
#            }
        }
        elsif($db_filter) {
            print STDERR "Filtered from organism_to_prefix ".$db_filter->{$qdb}." or ".$db_filter->{$sdb}."\n";
        }
    }
    elsif($options{'db_list'}) {
        print STDERR "Filtered $qdb or $sdb\n";
    }
}

sub read_db_list {
    open IN, "<".$options{'db_list'} || die "couldn't open '$options{db_list}' for reading: $!";
    while(<IN>) {
        chomp;
        $dbs->{$_} = 1;
        $db_to_org->{$_} = $_;
    }close IN;
    
}

sub read_organism_to_db_mapping {
   open IN, "<".$options{'organism_to_db_mapping'} or die "couldn't open '$options{organism_to_db_mapping}' for reading: $!";
    while(<IN>) {
        chomp;
        $_ =~ /^(.*)\s(\S+)$/;
        my $org = $1;
        my $db = $2;
        $org =~ s/\s/_/g;
        $db_filter->{$db} = 1;
        $db_to_org->{$db} = $org;
    }close IN;
  
}
