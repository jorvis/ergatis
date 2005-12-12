#!/usr/local/bin/perl

=head1 NAME

prepare_for_geneWise.dbi

=head1 SYNOPSIS

USAGE: prepare_for_geneWise.dbi
           --Database|-D                  aa1 
           --username|-u                  access 
           --password|-p                  access
           --search_db|-s                 AllGroup.niaa
           --work_dir|-w                  working_directory
         [ --asmbl_id|-a                  24832 
           --asmbl_file|-A                list_of_asmbl_ids 
           --file_list|-f                 genewise.input_file_listing    
           --verbose|-v
           --PADDING_LENGTH               500
           --MIN_CHAIN_SCORE              50
           --MIN_PERCENT_CHAIN_ALIGN      70
           --num_tiers|n                  2
           --JUST_PRINT_BEST_LOCATIONS
         ]


=head1 OPTIONS

B<--Database,-D>
    Annotation database name

B<--username,-u>
    username to access database

B<--password,-p>
    password to access database

B<--search_db,-s>
    name of fasta file already searched and with results already loaded into 
    the evidence table of the annotation database, ie. using AAT

B<--work_dir,-w>
    working directory to post the genome sequence segments and proteins to
    be aligned using genewise

B<--asmbl_id,-a>
    asmbl_id in annotation database to process

B<--asmbl_file,-A>
    file containing a list of asmbl_ids to process, whitespace delimited

B<--file_list,-f>
    name of the output file containing the list of input sequences to genewise
    default: "genewise.input_file_listing", written to the working directory.  
    
    You can specify a full path to your required output file list, or a 
    different name for this file.

B<--help,-h>
    This help documentation

B<--verbose,-v>
    Verbose mode 

B<--PADDING_LENGTH>
    number of basepairs to extend the genome location on each end (default: 500)

B<--MIN_PERCENT_CHAIN_ALIGN>
    mimimum percent of the matching proteins length that aligns to the genome in a single alignment chain.
    default: 70%

B<--MIN_CHAIN_SCORE>
    minimum score for an alignment chain to be considered a candidate for genewise realignment.
    chain_score = sum (per_id * segment_length)
    default minimum score = 50

B<--num_tiers,-n>
    number of overlapping best location hits.
    default: 2 (only the two best hits per location is extracted)

B<--JUST_PRINT_BEST_LOCATIONS>
    The best matches are reported to stdout.  No files are written.


=head1 DESCRIPTION

The best protein match to a given genomic region is extracted along with the genomic 
sequence corresponding to that region to be realigned using the genewise program.

The algorithm for finding the best match per location is as follows:  The results of an 
AAT search of a protein database against a genome assembly should be available for querying
from the evidence table.  All alignment chains are retrieved from the database and scored 
as the sum of (per_id * length) for each alignment segment.  The chains are sorted by score
and tiled to the genome, disallowing overlap among alignment chains along the genome.  Matches
to the top and bottom strands of the genome are tiled separately.  The first tier on each
strand contains the best hits per location.  Set --num_tiers > 1 to extract multiple best hits per
location.

Each protein is retrieved from the protein database fasta file using the cdbyank utility.  Each 
genome region is extracted as a substring of the genome sequence.  Both the protein and genome
sequence substring are written as fasta files in the working directory, structured like so:

{WORKING_DIRECTORY}/{asmbl_id}/{database}.assembly.{asmbl_id}.{counter}.pep

{WORKING_DIRECTORY}/{asmbl_id}/{database}.assembly.{asmbl_id}.{counter}.fsa

for the protein and genome sequence region, respectively.  The counter is an integer incremented for
each genome location on each asmbl_id.

The header of the .fsa file is constructed like so:

{database}.assembly.{asmbl_id}.{end5}.{end3}

so that, by extracting the header for this entry, the exact genome coordinates can be inferred,
and 

A list of complete paths to each .fsa file is stored as:

{WORKING_DIRECTORY}/{file_list}  

and to be used with subsequent genewise searches as part of workflow.

=head1 CONTACT

    Brian Haas
    bhaas@tigr.org

=cut

use strict;
use warnings;
use lib '/usr/local/devel/ANNOTATION/Euk_modules/bin';
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use DBI;
use Pod::Usage;
use CdbTools;
use Data::Dumper;
use Carp;

#option processing
my ($database, $username, $password, $search_db, $work_dir,
    $asmbl_id, $asmbl_file, $file_list, $help, $verbose, $JUST_PRINT_BEST_LOCATIONS);

# default settings.
my $DEFAULT_FILE_LIST = "genewise.input_file_listing";
my $PADDING_LENGTH = 500; # sequence extended from each end of genomic location.
my $MIN_CHAIN_SCORE = 50; # min score of alignment required for possible genewise realignment.
my $num_tiers = 2; # only the two best matches per location is extracted.
my $MIN_PERCENT_CHAIN_ALIGN = 70; #minimum percent of the protein sequence found to align to the genome

&GetOptions (
             'Database|D=s' => \$database,
             'username|u=s' => \$username,
             'password|p=s' => \$password,
             'search_db|s=s' => \$search_db,
             'work_dir|w=s' => \$work_dir,
             'asmbl_id|a=s' => \$asmbl_id,
             'asmbl_file|A=s' => \$asmbl_file,
             'file_list|f=s' => \$file_list,
             'help|h' => \$help,
             'verbose|v' => \$verbose,
             'PADDING_LENGTH=i' => \$PADDING_LENGTH,
             'MIN_CHAIN_SCORE=i' => \$MIN_CHAIN_SCORE,
             'MIN_PERCENT_CHAIN_ALIGN=f' => \$MIN_PERCENT_CHAIN_ALIGN,
             'num_tiers|n=i' => \$num_tiers,

             'JUST_PRINT_BEST_LOCATIONS' => \$JUST_PRINT_BEST_LOCATIONS,

	     ) || pod2usage();

if ($help) {
    pod2usage();
}

unless ($database && $username && $password && $search_db && $work_dir 
	&& ($asmbl_id || $asmbl_file)
	) {
    carp "** Missing required options. **\n\n";
    pod2usage();
}


## establish a database connection:
my $dbproc = DBI->connect("dbi:Sybase:server=SYBTIGR; packetSize=8092",$username, $password) 
    || croak "couldn't establish a database connection\n\n";
$dbproc->{RaiseError} = 1; # don't tolerate faulty db interactions
$dbproc->do("use $database");
$dbproc->do("set textsize 50000000");

## find search_db file:
my $protein_fasta_file = &find_fasta_file ($search_db) or die "Error, cannot find fasta file for $search_db";

# get db_id:
my $db_id = &get_db_id($search_db);

# track the protein alignment lengths:
my %prot_acc_to_prot_length;


## get list of asmbls to process
my @asmbls_to_process;
if ($asmbl_id) {
    @asmbls_to_process = ($asmbl_id);
} else {
    open (my $fh, $asmbl_file) or croak ("cannot open $asmbl_file ");
    while (<$fh>) {
        # should be whitespace delimited
	    while (/(\d+)/g) {
            my $asmbl = $1;
            push (@asmbls_to_process, $asmbl);
	    }
	}
}
if (! @asmbls_to_process) {
    croak "error, do not have any asmbl_ids to process ";
}

# check for working directory:
if (! -d $work_dir) {
    croak "error, cannot find work_dir: $work_dir ";
} else {
    # change curr dir to work_dir
    chdir ($work_dir) or croak "error, cannot cd to $work_dir ";
    
}

umask(0000); #write all output files permissibly


## open file for writing input file listings
if (! $file_list) {
    $file_list = $DEFAULT_FILE_LIST;
}

my $file_list_fh;
if (! $JUST_PRINT_BEST_LOCATIONS) {
    open ($file_list_fh, ">$file_list") or croak ("cannot write to $file_list ");
}

# get the pep and genome seqs for each location on each asmbl_id
foreach my $asmbl_id (@asmbls_to_process) {
    &prepare_asmbl_data($asmbl_id);
}

#$length_fetch->finish();

if (! $JUST_PRINT_BEST_LOCATIONS) {
    close $file_list_fh;
}

exit(0);


####
sub prepare_asmbl_data {
    my ($asmbl_id) = @_;
    
    print "Processing asmbl_id: $asmbl_id\n" if $verbose;
    
    ## retrieve best hits and genome locations
    my @best_hits = &get_best_location_hits($asmbl_id);
    
    if (@best_hits) {
		
        if (!$JUST_PRINT_BEST_LOCATIONS) {
            &prepare_genewise_inputs($asmbl_id, \@best_hits);
        }
        
    } 

    else {
        carp "warning, no best hits were returned for asmbl_id: $asmbl_id\n";
    }
}


#### 
sub get_best_location_hits {
    my ($asmbl_id) = @_;
    
    ## get the alignments out of the database:
    my $query = qq{select distinct e.accession, e.end5, e.end3, e.m_lend, e.m_rend, e.chainID, e.per_id from evidence e where feat_name like "$asmbl_id.%" and db = "$db_id" and e.ev_type = 'nap' order by e.chainID};
    
    print "Query: $query\n" if $verbose;
    
    my $sth = $dbproc->prepare($query);
    $sth->execute();
    
    ## Populate data structure:
    
    my %alignments; #key on chainID
    # chain struct:
    #    m_lend, m_rend, lend, rend, score, accession, orient
    
    while (my @row = $sth->fetchrow_array()) {
        my ($acc, $end5, $end3, $m_lend, $m_rend, $chainID, $per_id) = @row;
        
        unless ($chainID) {
            carp "Error, no chainID stored for @row\n";
            next;
        }
        
        $chainID = $acc . "," . $chainID;
        
        my ($lend, $rend) = sort {$a<=>$b} ($end5, $end3);
        my $match_length = $m_rend - $m_lend + 1;
        
        if ($match_length < 50) { next;} # avoid shorties that falsely extend alignment chains.
        
        my $match_score = $match_length * $per_id/100;
        
        #print "$acc\t$end5-$end3\t$m_lend-$m_rend\t$per_id\n";
        
        if (my $chain = $alignments{$chainID}) {
            # increment score
            $chain->{score} += $match_score;
            # adjust min/max coords for feature
            if ($chain->{lend} > $lend) {
                $chain->{lend} = $lend;
            }
            if ($chain->{rend} < $rend) {
                $chain->{rend} = $rend;
            }
            
            if ($chain->{m_lend} > $m_lend) {
                $chain->{m_lend} = $m_lend;
            }
            
            if ($chain->{m_rend} < $m_rend) {
                $chain->{m_rend} = $m_rend;
            }
        } else {
            # create chain entry:
            
            if ($end5 == $end3) { next; } # don't store single coord feature
            my $orient = ($end5 < $end3) ? '+' : '-';
            
            $alignments{$chainID} = { lend => $lend,
                                      rend => $rend,
                                      m_lend => $m_lend,
                                      m_rend => $m_rend,
                                      score => $match_score,
                                      orient => $orient,
                                      accession => $acc };
        }
    }
    
    if (! %alignments) {
        return (); # no alignments found in database
    }
    
    
    ## Store only the single best non-overlapping match for each strand
    my @top_strand_tiers;
    my @bottom_strand_tiers;
    
    ## add empty tiers:
    for (1..$num_tiers) {
        push (@top_strand_tiers, []);
        push (@bottom_strand_tiers, []);
    }

    
    my @sorted_chains = reverse sort {$a->{score}<=>$b->{score}} values %alignments;
    
    
    my @structs;
    foreach my $chain_struct (@sorted_chains) {
        
        # don't bother with low scoring alignments.
        my $score = $chain_struct->{score};
        if ($score < $MIN_CHAIN_SCORE) {
            next;
        }

        ## check the alignment length
        my ($m_lend, $m_rend) = ($chain_struct->{m_lend}, $chain_struct->{m_rend});
        my $align_len = $m_rend - $m_lend + 1;
        my $prot_acc = $chain_struct->{accession};
        

        my $prot_len = -1;
        eval {
            ## if the protein database (ie. AllGroup.niaa) has been updated since the last AAT searches
            ## we might not be able to recover certain entries from the database, in which case,
            ## this part will die!
            
            $prot_len = &get_protein_length_via_accession($prot_acc);
            
        };
        if ($@) {
            print STDERR "Error trying to extract protein and determine length for entry $prot_acc\n";
            next;
        }
        
        

        if ($align_len/$prot_len * 100 < $MIN_PERCENT_CHAIN_ALIGN) {
            next; #insufficient percent alignment to genome
        }
        
        ## Tile the hit
        if ($chain_struct->{orient} eq '+') {
            &try_add_chain($chain_struct, \@top_strand_tiers);
        } else {
            &try_add_chain($chain_struct, \@bottom_strand_tiers);
        }
    }
    

    
    ## move tiers to complete lists:
    my @top_strand_feats;
    foreach my $tier (@top_strand_tiers) {
        push (@top_strand_feats, @$tier);
    }
    my @bottom_strand_feats;
    foreach my $tier (@bottom_strand_tiers) {
        push (@bottom_strand_feats, @$tier);
    }
    
    if ($verbose || $JUST_PRINT_BEST_LOCATIONS) {
        print "All alignments: \n";
        my @all_alignments = sort {$a->{lend}<=>$b->{lend}} values %alignments;
        &dump_feats(@all_alignments);
        
        
        print "\n\nTop Strand Best Hits:\n";
        &dump_feats(@top_strand_feats);
        print "\n\nBottom Strand Best Hits:\n";
        &dump_feats(@bottom_strand_feats);
    }
    
    %alignments = (); #clear
    return (@top_strand_feats, @bottom_strand_feats);
}



sub prepare_genewise_inputs {
    my ($asmbl_id, $best_hits_aref) = @_;
    
    ## prepare sequence files
    if (! -d $asmbl_id) {
        mkdir ($asmbl_id) or croak "error, couldn't mkdir $asmbl_id in $work_dir "; 
    }
    
    ## get the genome sequence:    
    my $query = qq{select sequence from assembly where asmbl_id = $asmbl_id};
    
    ## need to set textsize here
    
    my $sth = $dbproc->prepare($query);
    $sth->execute();
    my ($genome_seq) = $sth->fetchrow_array;
    $sth->finish;
    
    my $genome_seq_length = length($genome_seq);
    
    my $counter = 0;
    
    ## Prepare files: 
    foreach my $chain (@$best_hits_aref) {
        $counter++;
        
        my ($lend, $rend, $m_lend, $m_rend, $accession, $score, $orient) = ($chain->{lend},
                                                                            $chain->{rend},
                                                                            $chain->{m_lend},
                                                                            $chain->{m_rend},
                                                                            $chain->{accession},
                                                                            $chain->{score},
                                                                            $chain->{orient});
        
		my $fasta_entry;
        eval {
            $fasta_entry = cdbyank($accession, $protein_fasta_file);
        };

        if ($@) {
            print STDERR "Error, couldn't retrieve fasta entry $accession from $protein_fasta_file\n";
            next;
        }
        

        $lend -= $PADDING_LENGTH;
        if ($lend <= 0) {
            $lend = 1;
        }
        $rend += $PADDING_LENGTH;
        if ($rend > $genome_seq_length) {
            $rend = $genome_seq_length;
        }
        
        my $subseq = substr($genome_seq, $lend -1, $rend - $lend + 1);
        
        if ($orient eq '-') {
            $subseq = &reverse_complement($subseq);
            ($lend, $rend) = ($rend, $lend);
        }
        $subseq =~ s/(\w{60})/$1\n/g;
        chomp $subseq;
	
        # write genome entry
        open (my $fh, ">$asmbl_id/$database.assembly.$asmbl_id.$counter.fsa") or die $!;
        print $fh ">$database.assembly.$asmbl_id.$lend.$rend\n$subseq\n";
        close $fh;
        
        # write protein entry:
        open ($fh, ">$asmbl_id/$database.assembly.$asmbl_id.$counter.pep") or die $!;
        print $fh $fasta_entry;
        close $fh;
        
        ## add entry in file listing:
        print $file_list_fh "$work_dir/$asmbl_id/$database.assembly.$asmbl_id.$counter.fsa\n";

    }
}
    

####
sub try_add_chain {
    my ($chain_struct, $tier_list_aref) = @_;
    
    ## if chain_struct doesn't overlap any existing feature, go ahead and add it:
    
    
    my ($chain_lend, $chain_rend, $chain_acc, $chain_score) = ($chain_struct->{lend}, 
                                                               $chain_struct->{rend},
                                                               $chain_struct->{accession},
                                                               $chain_struct->{score});
    
  TIERS:
    foreach my $feature_list_aref (@$tier_list_aref) {    
        my $found_overlap_flag = 0;
      FEATURES:
        foreach my $feat (@$feature_list_aref) {
            my ($feat_lend, $feat_rend) = ($feat->{lend}, $feat->{rend});
            if ($feat_lend < $chain_rend && $feat_rend > $chain_lend) {
                # got overlap:
                # check to see if the accession is the same.  
                #Don't want the same protein aligned to the same region in different tiers
                if ($feat->{accession} eq $chain_acc) {
                    last TIERS;  
                }
                $found_overlap_flag = 1;
                last FEATURES;
            }
        }
        if (! $found_overlap_flag) {
            ## add to current tier:
            push (@$feature_list_aref, $chain_struct);
            last TIERS;
        }
    }
}


####
sub find_fasta_file {
    my ($protein_db_name) = @_;
    
    my $DB = uc $database;
    
    foreach my $loc_dir ( "/usr/local/annotation/$DB/CustomDB",
                          "/usr/local/db/common",
                          "/usr/local/db/panda/AllGroup/",
                          "/usr/local/annotation/ACLA1/CustomDB") {
        my $candidate_filename = $loc_dir . "/$protein_db_name";
        if (-s $candidate_filename) {
            return ($candidate_filename);
        }
    }
    
    return (undef);
}



####
sub get_db_id {
    my ($search_db_token) = @_;
    
    ## entry must exist in the common..search_dbs table.
    
    my $query = qq{select id from common..search_dbs where name = "$search_db_token" and iscurrent = 1};
    print "QUERY: $query\n\n" if $verbose; 
    my $sth = $dbproc->prepare($query);
    $sth->execute();
    my ($db_id) = $sth->fetchrow_array();
    $sth->finish;
    
    if (! $db_id) {
        croak "Error, couldn't find an entry in the common..search_dbs table for $search_db_token ";
    }
    return ($db_id);
}


####
sub dump_feats {
    my @feats = @_;
    
    @feats = sort {$a->{lend}<=>$b->{lend}} @feats;
    
    foreach my $feat (@feats) {
        my ($lend, $rend, $m_lend, $m_rend, $score, $accession) = ($feat->{lend},
                                                                   $feat->{rend},
                                                                   $feat->{m_lend},
                                                                   $feat->{m_rend},
                                                                   $feat->{score},
                                                                   $feat->{accession});
        
        
        $score = sprintf ("%.2f", $score);
        print "$lend-$rend\t$accession\t$score\t$m_lend-$m_rend\n";
    }
}


sub reverse_complement { # from Egc_library.pm
    my($s) = @_;
    my ($rc);
    $rc = reverse ($s);
    $rc =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
     
    return($rc);
}

sub get_protein_length_via_accession {
    my ($prot_acc) = @_;
    
    if (my $length = $prot_acc_to_prot_length{$prot_acc}) {
        return ($length);
    }

    else {
        print "-retrieving protein length for $prot_acc\n" if $verbose;
        my $fasta_entry = cdbyank($prot_acc, $protein_fasta_file);
        my ($acc, $header, $seq) = linearize($fasta_entry);
        
        my $length = length($seq);
        $prot_acc_to_prot_length{$prot_acc} = $length;
        return ($length);
    }
}

