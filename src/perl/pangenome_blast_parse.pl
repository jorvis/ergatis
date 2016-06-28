#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

=head1  NAME

pangenome_blast_parse.pl - Parses BLAST BSML output files and extracts data needed for pangenome analysis.

=head1 SYNOPSIS

USAGE: pangenome_query_list.pl
        --input=/path/to/somefile.bsml
        --output_path=/path/to/output/directory/
        --coverage_cutoff=50
        --similarity_cutoff=50
      [ --log=/path/to/some/log
        --organism_to_db_mapping=/path/to/organism_to_db_mapping.txt
        --db_list=/path/to/database.list
      ]

=head1 OPTIONS

B<--input,-i>
    BSML file containing BLAST results to parse.

B<--output_path,-o>
    Directory to write stored file

B<--coverage_cutoff>
    Coverage cutoff. Hits below this will not be considered hits.

B<--similarity_cutoff>
    Similarity cutoff. Hits below this will not be considered hits.

B<--organism_to_db_mapping>
    Mapping file of the format:
    Genus Species strain uniquename_prefix
    Streptococcus pneumoniae TIGR4 sptigr4

B<--db_list>
    List of uniquename prefixes that should be included in the analysis.
    It is not necessary to include this list if you are filtering with the
    organism_to_db_mapping file.

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

my $TOP_HIT_CUTOFF = 1;	# Max number of good hits from query gene to subject genome to keep

my %options = ();
my $results = GetOptions (  \%options,
                            'input|i=s',
                            'output_path|o=s',
                            'organism_to_db_mapping|odb:s',
                            'db_list|dl:s',
                            'coverage_cutoff:s',
                            'similarity_cutoff:s',
                            'debug|d=s',
                            'command_id=s',       ## passed by workflow
                            'logconf=s',          ## passed by workflow (not used)
                            'log|l=s',
                            'help|h',
                         ) || pod2usage();

my %hit_results = ();
my $dups_temp = {};
my %dups = ();
my $db_to_org = {};
my $db_filter = undef;
my %good_hit_counts = {};
my $coverage_cutoff = $options{'coverage_cutoff'} ne '' ? $options{'coverage_cutoff'} : 50;
my $similarity_cutoff = $options{'similarity_cutoff'} ne '' ? $options{'similarity_cutoff'} : 50;


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

# Foreach gene that has genuine duplicates besides itself, add to dups hash.
foreach my $org (keys %$dups_temp) {
	foreach my $prot ( keys %{$dups_temp->{$org}}) {
    	if (scalar(@{$dups_temp->{$org}->{$prot}}) > 0) {
    		# Add the query gene in the dup string to complete the set
    		push (@{$dups_temp->{$org}->{$prot}}, $prot);
    		# Sort array to make all identical arrays have the same order
    		my @dup_array = sort {$a cmp $b} @{$dups_temp->{$org}->{$prot}};
    		my $dup_set = join(" ", @dup_array);
    	    $dups{$org}{$dup_set} = 1;
    	}
    }
}

#foreach my $q_genome (keys %good_hit_counts){
#	foreach my $q_gene (keys %{$good_hit_counts{$q_genome}}){
#		foreach my $s_genome (keys %{$good_hit_counts{$q_genome}{$q_gene}}){
#			print STDERR "$q_genome\t$q_gene\t$s_genome\t" . $good_hit_counts{$q_genome}{$q_gene}{$s_genome} . "\n";
#		}
#	}
#}

nstore([\%hit_results,\%dups], $options{'output_path'}."/".$input_prefix.".blast.stored") || die "couldn't serialize results";

exit(0);

sub processSeqPairAlignment {
    my ($twig, $feat) = @_;
    my ($p_sim, $p_ident, $p_cov_ref, $p_cov_comp, $this_chain);

    #This is used if we want to split via underscore
    #my $query_id   = $feat->{'att'}->{'refseq'};
    #my $subject_id = $feat->{'att'}->{'compseq'};
    #chop $query_id;	#Query ID has a final "_" a the end.  Want to make uniform with Query ID for regex later

    #Finding it difficult to break up the sequence on underscores since some of the DB names and hit names have them
    my $query_id_line = $feat->{'att'}->{'refxref'};
    my $subject_id_line = $feat->{'att'}->{'compxref'};
    my ($nope, $query_id) = split(':', $query_id_line, 2);
    my ($nada, $subject_id) = split(':', $subject_id_line, 2);
    
    my $method = $feat->{'att'}->{'method'};

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

    # $qprot may be an assembly for TBLASTN, not a protein
    my($qdb, $qprot) = ($2, $1);

    # For general database identifiers (gnl|qdb|qprot)
    if ($query_id =~/^(gnl|gb)\|([^\|]+)\|([^\.\|]+)[\.\|]?.*$/) {
    	($qdb, $qprot) = ($2, $3);
#    	print "QDB:\t$qdb\tQPROT:\t$qprot\n";
    }
    elsif ($query_id =~ /^(([^\.]+)\..*)$/) {
		($qdb, $qprot) = ($2, $1);
    }
    # handle BSML sequences (for NC GenBank entries) created via genbank2bsml
    elsif ($query_id =~ /^ref\|(NC_\d+)\|(\S+)$/) {
    #elsif ($query_id =~ /^ref_(NC_\d+)_(\S+)$/) {
		($qdb, $qprot) = ($2, $1);
    } else {
	die "couldn't parse out db and sequence id from query\n$query_id\n";
    }

    # $sprot may be an assembly for TBLASTN, not a protein
    my($sdb, $sprot);
    # For general database identifiers (gnl|qdb|qprot)
    if ($subject_id =~/^(gnl|gb)\|([^\|]+)\|([^\.\|]+)[\.\|]?.*$/){
    	($sdb, $sprot) = ($2, $3);
#    	print "SDB:\t$sdb\tSPROT:\t$sprot\n";
    }
    elsif ($subject_id =~ /^(([^\.]+)\..*)$/) {
        ($sdb, $sprot) = ($2, $1);
    }
    # handle BSML sequences (for NC GenBank entries) created via genbank2bsml
    elsif ($subject_id =~ /^ref\|(NC_\d+)\|(\S+)$/) {
    #elsif ($subject_id =~ /^ref_(NC_\d+)_(\S+)$/) {
		($sdb, $sprot) = ($2, $1);
    } else {
	die "couldn't parse out db and sequence id from subject\n$subject_id\n";
    }

    if(!$options{'db_list'} && !$options{'organism_to_db_mapping'}) {
        $db_to_org->{$qdb} = $qdb;
        $db_to_org->{$sdb} = $sdb;
    }

    # If we are doing any filtering then we will check that both the database names are good.
    if (!$db_filter || ($db_filter->{$qdb} && $db_filter->{$sdb})) {
		# Filter out the self-hits so they aren't run through cutoff filtering (blastp results)
		next if ($db_to_org->{$qdb} eq $db_to_org->{$sdb} && ($qprot eq $sprot));
		# Filter protein queries that are self-hits with the nucleotide genome/contig (tblastn results)
		next if ($db_to_org->{$qdb} eq $db_to_org->{$sdb} && $method == "TBLASTN" && $all_seg_p_ident == 100 && $p_cov_ref == 100 && $p_cov_comp == 100);
		$good_hit_counts{$qdb}{$qprot}{$sdb} = 0 unless ($good_hit_counts{$qdb}{$qprot}{$sdb});
		# Does the sequence hit meet the cutoff for being good?
		if ($all_seg_p_sim >= $similarity_cutoff && ($p_cov_ref >= $coverage_cutoff || $p_cov_comp >= $coverage_cutoff)) {
        	if ($db_to_org->{$qdb} eq $db_to_org->{$sdb} && $all_seg_p_ident == 100 && $p_cov_ref == 100 && $p_cov_comp == 100) {
			    # If we have a non-self (different gene) hit on the same genome that exceeds the cutoffs, send to duplicates hash
            	push(@{$dups_temp->{$db_to_org->{$sdb}}->{$qprot}}, $sprot);
				#print STDERR "DUPS $qdb -- $qprot $sdb -- $sprot $p_cov_ref $p_cov_comp $all_seg_p_sim\n";
        	} else{
				# If we are here, assume the query and subject genomes are different
				#print join("\t", $db_to_org->{$qdb},$qprot,$db_to_org->{$sdb},$sprot,$all_seg_p_sim,$p_cov_ref) . "\n";
				
				# Only add subject genome into list of "good-hit" genomes for a given query gene if the subject genome is not already in the array.  
            	push (@{$hit_results{$qdb}{$qprot}}, $sdb) unless $good_hit_counts{$qdb}{$qprot}{$sdb}++ ;
        	}
        }
    }
    elsif($options{'organism_to_db_mapping'}) {
        print STDERR "Filtering out from organism_to_prefix -- query ".$db_filter->{$qdb}." or subject".$db_filter->{$sdb}."\n";
    }
    elsif($options{'db_list'}) {
        print STDERR "Filtering out from db_list -- query ".$db_filter->{$qdb}." or subject".$db_filter->{$sdb}."\n";
    }
}

sub read_db_list {
    open IN, "<".$options{'db_list'} || die "couldn't open '$options{db_list}' for reading: $!";
    while(<IN>) {
        chomp;
        $db_filter->{$_} = 1;
        $db_to_org->{$_} = $_;
    }
	close IN;
}

sub read_organism_to_db_mapping {
   open IN, "<".$options{'organism_to_db_mapping'} or die "couldn't open '$options{organism_to_db_mapping}' for reading: $!";
    while(<IN>) {
        chomp;
        $_ =~ /^(.*)\s(\S+)$/;
        my $org = $1;
        my $db = $2;
        $org =~ s/\s/_/g;
        if(!$options{db_list}) {$db_filter->{$db} = 1};
        $db_to_org->{$db} = $org;
    }
	close IN;
}
