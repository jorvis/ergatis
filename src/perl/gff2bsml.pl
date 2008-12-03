#!/usr/bin/perl

# NOTES
#  -allow multiple transcripts/gene (each gets its own BSML feature group)
#  -allow id-less CDS features (with --use_cds_parent_ids), so long as parent ids are unique
#  -removed analysis and analysis links
#  -added miRNA support (and, more generally, support for multiple transcript-level types)
#  -ditto for tRNA, ncRNA, snoRNA, snRNA, rRNA
#  -made $feat_type_map into a global variable
#  -added --default_seq_class
#  -added --peptide_fasta, --peptide_id_regex, --peptide_id_prefix, --peptide_id_suffix to handle external file of peptide seqs
#  -added --dna_fasta, --dna_id_regex, --dna_id_prefix, --dna_id_suffix to handle external file of DNA seqs
#  -removed "temp_id" strings from feature element titles
#  -added 'organism' input that gets converted into a BSML <Genome> section
#  -added --organism_genetic_code, --organism_mt_genetic_code
#   -all Sequences link to that <Genome>
#   -reference Sequence lengths are determined from the GFF, if possible
#  -modified Seq-data-imports to use absolute paths

# Current restrictions/assumptions
#  -all sequences/features belong to 1 genome, the one named by --organism param

# TODO
#  -add a strict mode in which deviations from the GFF3/canonical gene spec. are reported
#   -see http://www.sequenceontology.org/gff3.shtml
#  -if --use_cds_parent_ids is in effect the script should check that the CDS features under each mRNA are
#   nonoverlapping: if not it may be impossible to correctly reconstruct the original gene model without
#   correctly-assigned CDS ids
#  -check handling of identifiers in Seq-data-imports: can a shorter id be used without rewriting the FASTA file?
#  -add options to control the handling of '*' characters in the peptide sequences
#   -current legacy2bsml/genbank2bsml scripts appear to be inconsistent wrt whether '*' is included
#   -and if the '*' is included then it gets counted in the seqlen in chado
#  -create Research/Analysis elements for computed elements in the BSML (add gene_model analysis param?)
#  -add --ref_seq_type param? change name of --default_seq_class param
#  -need a way to specify which features get linked to Analyses in the <Research> section

=head1 NAME

gff32bsml.pl - Convert GFF(3) annotation data into BSML.

=head1 SYNOPSIS

gff32bsml.pl
         --input=/path/to/annotation.gff3
         --output=/path/to/results.bsml
         --project=rca1
         --id_repository=/usr/local/projects/rca1/workflow/project_id_repository
         --organism='Aspergillus nidulans'
        [--peptide_fasta=/path/to/peptide-seqs.fsa
         --peptide_id_regex='^>(\S+)'
         --peptide_id_prefix=''
         --peptide_id_suffix='-Protein'
         --dna_fasta=/path/to/dna-seqs.fsa
         --dna_id_regex='^>(\S+)'
         --dna_id_prefix='Chr'
         --dna_id_suffix=''
         --organism_genetic_code=1
         --organism_mt_genetic_code=4
         --use_cds_parent_ids
         --default_seq_class=assembly
         --log=/path/to/some.log
         --debug=4
         --help
         --man ]

=head1 OPTIONS

B<--input, -g>
    path to the input GFF file.

B<--peptide_fasta, -e>
    optional.  path to a multi-FASTA file containing the polypeptides for the genes in the GFF file.
    if this argument is supplied then the named file _must_ contain all of the polypeptides/proteins
    in the input GFF file.

B<--peptide_id_regex>
    optional.  regular expression used to parse peptide unique id from FASTA deflines in --peptide_fasta.
    default = ^>(\S+)

B<--peptide_id_prefix>
    optional. prefix to prepend to the id parsed by --peptide_id_regex to produce the polypeptide ID 
    used in the GFF file.

B<--peptide_id_suffix>
    optional. suffix to append to the id parsed by --peptide_id_regex to produce the polypeptide ID 
    used in the GFF file.

B<--dna_fasta>
    optional.  path to a multi-FASTA file containing the DNA sequences for the reference sequences
    in the GFF file.  if this argument is supplied then the named file _must_ contain all of the 
    reference sequences in the input GFF file.

B<--dna_id_regex>
    optional.  regular expression used to parse the sequence unique id from FASTA deflines in 
    --dna_fasta.   default = ^>(\S+)

B<--dna_id_prefix>
    optional. prefix to prepend to the id parsed by --dna_id_regex to produce the sequence ID 
    used in the GFF file.

B<--dna_id_suffix>
    optional. suffix to append to the id parsed by --dna_id_regex to produce the sequence ID 
    used in the GFF file.

B<--organism_genetic_code>
    optional.  specifies a genetic code value to appear in the BSML <Organism> element for --organism
    see http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c#SG1 for values

B<--organism_mt_genetic_code>
    optional.  specifies a mitochondrial genetic code value to appear in the BSML <Organism> element 
    for --organism.  see http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c#SG1 for values

B<--output, -o>
    output BSML file (which should not exist before running the script.)

B<--project, -p>
    project/database name used to create unique feature identifiers.

B<--id_repository, -i>
    path to the Ergatis project_id_repository for the current project.

B<--organism, -r>
    organism to place in the BSML <Genomes> section.  the first word will be used as the BSML genus,
    and the rest of the string will be placed in the BSML species attribute.

B<--use_cds_parent_ids,-c>
    optional.  whether to use the GFF parent id to identify CDS features that lack an ID of their own.

B<--default_seq_class,-e>
    optional.  default sequence class/SO type (e.g., 'assembly', 'supercontig') to use for sequences 
    for which the type cannot be parsed from the GFF file.

B<--log,-l>
    optional.  path to a log file the script should create.  will be overwritten if
    it already exists.

B<--debug,-d>
    optional.  the debug level for the logger (an integer)

B<--help,-h>
    print usage/help documentation

B<--man,-m>
    print detailed usage/help documentation

=head1 DESCRIPTION

This script parses a single GFF input file containing annotation data (for one or more
genomic sequences) and writes it out as a BSML document that can be subsequently fed 
into bsml2chado to load it into a Chado database.  This script was originally created 
by copying and then generalizing Brett Whitty's evmgff32bsml.pl utility.

=head1 INPUT

A single GFF input file.  See the GFF documentation for more details.  Note that the
script has a number of options that allow it to correctly parse/convert a variety of
GFF-format files, even those that may deviate from the official specification in one
way or another.

=head1 OUTPUT

A single BSML document.  See the BSML documentation for more details.

=head1 CONTACT

    Jonathan Crabtree
    jonathancrabtree@gmail.com

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Data::Dumper;
use BSML::BsmlBuilder;
use Ergatis::IdGenerator;
use Ergatis::Logger;
use FileHandle;
use File::Spec;

## global/default values    
my $DEFAULT_SEQ_CLASS = '';
my $GLOBAL_ID_COUNTER = 0;
my $NODES = {};

# mapping from GFF3 feature type to Chado cvterm.  if the script encounters
# a GFF type not found in this hashref it will fail
my $FEAT_TYPE_MAP = {
    'gene' => 'gene',
    'CDS' => 'CDS',
    'exon' => 'exon',
    'mRNA' => 'transcript',               
    'miRNA' => 'miRNA',
    'protein' => 'polypeptide',
    'polypeptide' => 'polypeptide',
    'tRNA' => 'tRNA',
    'ncRNA' => 'ncRNA',
    'snoRNA' => 'snoRNA',
    'snRNA' => 'snRNA',
    'rRNA' => 'rRNA',
    'three_prime_UTR' => 'three_prime_UTR',
    'five_prime_UTR' => 'five_prime_UTR',
};

# list of types that are allowed for the mRNA; must be a subset of those
# listed _as values_ in $FEAT_TYPE_MAP
my $MRNA_TYPES = {
    'transcript' => 1,
    'miRNA' => 1,
    'tRNA' => 1,
    'ncRNA' => 1,
    'snoRNA' => 1,
    'snRNA' => 1,
    'rRNA' => 1,
};

## input/options
my $options = {};
my $results = GetOptions($options, 
                         'input|g=s',
                         'peptide_fasta=s',
                         'peptide_id_regex=s',
                         'peptide_id_prefix=s',
                         'peptide_id_suffix=s',
                         'dna_fasta=s',
                         'dna_id_regex=s',
                         'dna_id_prefix=s',
                         'dna_id_suffix=s',
                         'organism_genetic_code=i',
                         'organism_mt_genetic_code=i',
                         'output|o=s',
                         'project|p=s',
                         'id_repository|i=s',
                         'organism|r=s',
                         'use_cds_parent_ids|c',
                         'default_seq_class|e=s',
                         'log|l=s',
                         'debug|d=i',
                         'help|h',
                         'man|m'
                         ) || pod2usage();

&pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT}) if ($options->{'help'});
&pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if ($options->{'man'});

&check_parameters($options);

## initialize logging
my $logfile = $options->{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile, 'LOG_LEVEL'=>$options->{'debug'});
$logger = Ergatis::Logger::get_logger();

## MAIN SECTION
my $idcreator = Ergatis::IdGenerator->new('id_repository' => $options->{'id_repository'});
my @root_nodes = ();
my($genus, $species) = ($options->{'organism'} =~ /^(\S+) (.*)$/);

## parse peptide FASTA (optional)
my $peptides = undef;

if (defined($options->{'peptide_fasta'})) {
    my $pep_prefix = $options->{'peptide_id_prefix'};
    my $pep_suffix = $options->{'peptide_id_suffix'};
    my $pep_id_fn = sub {
        my($id) = @_;
        my $new_id = (defined($pep_prefix)) ? $pep_prefix . $id : $id;
        $new_id .= $pep_suffix if (defined($pep_suffix));
        return $new_id;
    };
    $peptides = &read_sequence_lengths($options->{'peptide_fasta'}, $options->{'peptide_id_regex'}, $pep_id_fn);
}

## parse peptide FASTA (optional)
my $dna_seqs = undef;

if (defined($options->{'dna_fasta'})) {
    my $dna_prefix = $options->{'dna_id_prefix'};
    my $dna_suffix = $options->{'dna_id_suffix'};
    my $dna_id_fn = sub {
        my($id) = @_;
        my $new_id = (defined($dna_prefix)) ? $dna_prefix . $id : $id;
        $new_id .= $dna_suffix if (defined($dna_suffix));
        return $new_id;
    };
    $dna_seqs = &read_sequence_lengths($options->{'dna_fasta'}, $options->{'dna_id_regex'}, $dna_id_fn);
}

## parse GFF
open (IN, $options->{'input'}) || die "couldn't open input file '$options->{input}'";

while (<IN>) {
    chomp;
    
    if (/^\#/) {
        # process header line
    } elsif (/^>/ || /^\#.*fasta.*/i) {

        last;
    } else {
        my $record = parse_record($_);

        # if it has no parents it's a root node    
        if ((!$record->{'Parent'}) && (!$record->{'Derives_from'})) {
            push(@root_nodes, $record);
        }
    
        # store records in a set of arrays indexed by ID
        if ($record->{'ID'}) {
            push(@{$NODES->{$record->{'ID'}}->{'records'}},$record);
        }
    }
}
close IN;

# populate children for each record
foreach my $id(keys(%{$NODES})) {
    foreach my $record(@{$NODES->{$id}->{'records'}}) {
        if (defined($record->{'children'})) {
            print STDERR "WARNING - about to overwrite children field of record $record with id=$id\n";
        }
        $record->{'children'} = $NODES->{$id}->{'children'};
        delete $NODES->{$id}->{'children'};
    }
}

my $gene_nodes = [];
fetch_node_type('gene', \@root_nodes, $gene_nodes);

my $doc = new BSML::BsmlBuilder();

# add Genomes
my $genome = $doc->createAndAddGenome();
my $genome_id = $genome->{attr}->{id};
my $organism = $doc->createAndAddOrganism('genome' => $genome, 'genus' => $genus, 'species' => $species);
my $code = $options->{'organism_genetic_code'};
my $mt_code = $options->{'organism_mt_genetic_code'};

my $code_att = $doc->createAndAddBsmlAttribute($organism, 'genetic_code', $code) if (defined($code));
my $mt_code_att = $doc->createAndAddBsmlAttribute($organism, 'mt_genetic_code', $code) if (defined($mt_code));

foreach my $root(@{$gene_nodes}) {
    if ($root->{'_type'} eq 'gene') {
        my $features = {};
        process_node($root, $features);
        $features->{'_seqid'} = $root->{'_seqid'};
        gene_feature_hash_to_bsml($NODES, $features, $peptides, $dna_seqs);
    }
}

# add the analysis element
# TODO - need a way to specify which features get linked to an analysis
#$doc->createAndAddAnalysis(
#                            id => 'EVM_analysis',
#                            sourcename => $options->{'output'},
#                          );

$doc->write($options->{'output'});
exit(0);

## subroutines

sub check_parameters {
	my $options = shift;

    ## make sure required parameters were passed
    my @required = qw(input output project id_repository organism);
    for my $option ( @required ) {
        unless  ( defined $options->{$option} ) {
            die "--$option is a required option";
        }
    }

    ## default values
    $options->{'default_seq_class'} = $DEFAULT_SEQ_CLASS if (!defined($options->{'default_seq_class'}));
    $options->{'project'} =~ tr/A-Z/a-z/;

    ## check files
    my $input = $options->{'input'};
    die "input GFF file ($input) does not exist" if (!-e $input);
    die "input GFF file ($input) is not readable" if (!-r $input);

    my $peps = $options->{'peptide_fasta'};
    if (defined($peps)) {
        die "input peptide FASTA file ($peps) does not exist" if (!-e $peps);
        die "input peptide FASTA file ($peps) is not readable" if (!-r $peps);
    }

    my $dnas = $options->{'dna_fasta'};
    if (defined($dnas)) {
        die "input DNA sequence FASTA file ($dnas) does not exist" if (!-e $dnas);
        die "input DNA sequence FASTA file ($dnas) is not readable" if (!-r $dnas);
    }

    if ($options->{'organism'} !~ /^\S+\s\S+/) {
        die "--organism must specify a species and genus, e.g., --organism='Aspergillus nidulans'";
    }
}

# Parse a single GFF3 record.
#
# $line - A chomp'ed GFF3 feature line.
#
sub parse_record {
    my ($line) = @_;

    my $record = {};
    my $attrib_hash = {};
    my $record_id = '';
    my $parent_id = '';
    
    my @cols = split("\t");

    # convert to interbase coordinates:
    $cols[3]--;

    # change the + and - symbols in strand column to 0 and 1, respectively
    if ($cols[6] eq '+') {
        $cols[6] = 0;
    } elsif ($cols[6] eq '-') {
        $cols[6] = 1;
    } elsif ($cols[6] eq '.') {
        $cols[6] = 0;
        $logger->warn("converting feature with strand='.' to BSML feature with complement=false");
    } else {
        die("unknown value ($cols[6]) in strand column.  expected '+', '-', or '.'.");
    }

    $record->{'_seqid'}     = $cols[0];
    $record->{'_source'}    = $cols[1];
    $record->{'_type'}      = $cols[2];
    $record->{'_start'}     = $cols[3];
    $record->{'_end'}       = $cols[4];
    $record->{'_score'}     = $cols[5];
    $record->{'_strand'}    = $cols[6];
    $record->{'_phase'}     = $cols[7];
    
    my @attribs = split(";", $cols[8]);
    foreach my $attrib(@attribs) {
        my ($type, $val) = split("=", $attrib);
        my @vals = split(",", $val);
        $record->{$type}=\@vals;
    }

    if (!defined($record->{'ID'})) {
#        print STDERR "no ID defined for record $line\n";
    } else {
        $record->{'ID'} = $record->{'ID'}->[0];
    }

    my $parents = [];

    if (!defined($record->{'Parent'})) {
        if (defined($record->{'Derives_from'})) {
            push(@$parents, @{$record->{'Derives_from'}});
        }
#        print STDERR "no Parent attribute defined for record $line\n";
    } else {
        push(@$parents, @{$record->{'Parent'}});
    }

    foreach my $parent(@$parents) {
        # store record's hash reference as a child of Parent
        if (!defined($NODES->{$parent}->{'children'})) {
            $NODES->{$parent}->{'children'}->{$record->{'_type'}} = [];
        }
        push (@{$NODES->{$parent}->{'children'}->{$record->{'_type'}}},$record);
    }
    
    return $record;
}

# recursively traverse a node and all its children
# and process each record to create a feature hash
#
# $node - 
# $features - 
#
sub process_node {
    my ($node, $features) = @_;

    # hash of keys to ignore when extracting a record hash from a node hash
    my $ignore_keys = {
                        'children' => 1,
                      };
    
    # build a record hash from values stored in the node hash
    my $record;
    foreach my $key(keys %{$node}) {
        if (!$ignore_keys->{$key}) {
            $record->{$key} = $node->{$key};
        }
    }
    
    # process the record hash
    process_record($record, $features);
    
    # process the node's children
    foreach my $child_type(keys %{$node->{'children'}}) {
        foreach my $child_record(@{$node->{'children'}->{$child_type}}) {
                process_node($child_record, $features); 
        }
    }
    
    return; 
}

# process the records and store them in the features hash
sub process_record {
    my ($record, $features) = @_;
    
    print STDERR "process_record called on " . $record . "\n" if ($record->{'_type'} eq 'polypeptide');

    if ($record->{'_type'} eq 'CDS') {
        # CDS records can span lines and must be merged into one CDS feature
        if (!$record->{'ID'}) {
            if ($options->{'use_cds_parent_ids'}) {
                if ($record->{'Parent'}) {
                    my $trans_id = &get_parent_id_by_type({'parent' => $record->{'Parent'}}, $features, $MRNA_TYPES);
                    die "no parent id found for CDS feature" if (!defined($trans_id));
                    $record->{'ID'} = $trans_id;
                    $logger->debug("setting CDS parent id to " . $trans_id);
                } else {
                    die "CDS feature lacks ID and also has no Parent ID for --use_cds_parent_ids!";
                }
            } else {
                die "CDS feature lacks ID -> bad form!";
            }
        }
        if ($features->{'CDS'}->{$record->{'ID'}}) {
            if ($features->{'CDS'}->{$record->{'ID'}}->{'startpos'} > $record->{'_start'}) {
                $features->{'CDS'}->{$record->{'ID'}}->{'startpos'} = $record->{'_start'};
            }
            if ($features->{'CDS'}->{$record->{'ID'}}->{'endpos'} < $record->{'_end'}) {
                $features->{'CDS'}->{$record->{'ID'}}->{'endpos'} = $record->{'_end'};
            }
        } else {
            $features->{$FEAT_TYPE_MAP->{'CDS'}}->{$record->{'ID'}} = {
                                                    'complement'    => $record->{'_strand'},
                                                    'startpos'      => $record->{'_start'},
                                                    'endpos'        => $record->{'_end'},
                                                    'parent'        => $record->{'Parent'},
                                                    'type'          => $record->{'_type'},
                                                };
        }
        
    } else {
        #handle all other feature types    
        
        my $feat_type;
        
        if (!defined($FEAT_TYPE_MAP->{$record->{'_type'}})) {
            print STDERR "unexpected feature type '$record->{_type}'\n";
            $feat_type = $record->{'_type'};    
        } else {
            $feat_type = $FEAT_TYPE_MAP->{$record->{'_type'}};
        }

        my $id;
        if (!$record->{'ID'}) {
            $id = getTempId();
        } else {
            $id = $record->{'ID'};
        }
        
        my $title = undef;
        if ($record->{'Name'}) {
            my $name = shift(@{$record->{'Name'}});
            $title = $name if ($name !~ /^temp_id/);
        }

        $features->{$feat_type}->{$id} = {  
                    'complement'    => $record->{'_strand'},
                    'startpos'      => $record->{'_start'},
                    'endpos'        => $record->{'_end'},
                    'title'         => $title,
                    'parent'        => $record->{'Parent'} || $record->{'Derives_from'},
                    'type'          => $record->{'_type'},
                };
    }
    
}

# find the parent of a feature, checking only those parents in a specified
# set of parent types.  returns undef if a unique parent of one of the specified 
# types cannot be found
#
# $feat - feature whose parent to find
# $features - all features
# $types - hashref whose keys are the types of parents to check
#
sub get_parent_id_by_type {
    my($feat, $features, $types) = @_;
    my $parent_list = $feat->{'parent'};
    my $matches = [];

    foreach my $type (keys %$types) {
        my $type_list = $features->{$type};
    
        foreach my $parent_id (@$parent_list) {
            push(@$matches, $parent_id) if (defined($type_list->{$parent_id}));
        }
    }

    my $nm = scalar(@$matches);
    if ($nm == 0) {
        print Dumper $feat;
        die "no parent found for feature $feat->{title} with parents = " . join(',', @$parent_list);
    } elsif ($nm > 1) {
        print Dumper $feat;
        die "multiple parents found for feature $feat->{title} with parents = " . join(',', @$parent_list);
    } 

    return $matches->[0];
}

# convert a gene feature hash into BSML
sub gene_feature_hash_to_bsml {
    my ($nodes, $features, $peptides, $dna_seqs) = @_;
    
    my $seq_id = $features->{'_seqid'};
    delete $features->{'_seqid'};

    my $id_hash = {};
    my %id_counts;
    
    foreach my $type(keys(%{$features})) {
        $id_counts{$type} = scalar(keys(%{$features->{$type}}));
    }
    $idcreator->set_pool_size(%id_counts);
    
    my $seq;
    
    # create a sequence stub for the seq_id if it doesn't exist yet
    if (!($doc->returnBsmlSequenceByIDR($seq_id))){

        # determine length of sequence from feature in GFF, if present
        my $gff_seq = $nodes->{$seq_id};
        my $gff_seqlen = undef;

        if (defined($gff_seq)) {
            my $records = $gff_seq->{'records'};
            my $nr = scalar(@$records);
            if ($nr != 1) {
                $logger->warn("$nr GFF record(s) found for sequence with id = $seq_id");
            } else {
                # interbase coords in effect here
                $gff_seqlen = $records->[0]->{'_end'} - $records->[0]->{'_start'};
            }
        } else {
            $logger->warn("unable to determine length of Sequence $seq_id: corresponding feature not found in GFF");
        }

        $seq = $doc->createAndAddSequence( $seq_id, $seq_id, $gff_seqlen, 'dna', $options->{'default_seq_class'} );
        my $genome_link = $doc->createAndAddLink($seq, 'genome', '#' . $genome_id);

        if (defined($dna_seqs)) {
            my $dseq = $dna_seqs->{$seq_id};
            die "couldn't find DNA sequence for $seq_id" if (!defined($dseq));
            my $path = File::Spec->rel2abs($options->{'dna_fasta'});
            my $seq_data_import_elem = $doc->createAndAddSeqDataImport($seq, 'fasta', $path, undef, $dseq->{'defline'});
        }

        # TODO - need a way to specify which features get linked to an analysis
#        $seq->addBsmlLink('analysis', '#EVM', 'input_of');
    } else {
        $seq = $doc->returnBsmlSequenceByIDR($seq_id);
    }

    # each of the feature types listed in $MRNA_TYPES can be used as a transcript-level feature
    my @transcript_ids = ();
    foreach my $ttype (keys %$MRNA_TYPES) {
        my @mrna_ids = map {{'id' => $_, 'type' => $ttype}} keys(%{$features->{$ttype}});
        push(@transcript_ids, @mrna_ids);
    }
    my $feat_table = $doc->createAndAddFeatureTable($seq);

    # $feat_groups indexed by transcript id
    my $feat_groups = {};

    foreach my $trans_id(@transcript_ids) {
        my $t_id = $idcreator->next_id( 
                                        'project' => $options->{'project'}, 
                                        'type' => $trans_id->{'type'},
                                        );

        if (defined($feat_groups->{$trans_id->{'id'}})) {
            die "duplicate transcript id ($trans_id->{id}) encountered";
        } else {
            $id_hash->{$trans_id->{'id'}} = $t_id;
            $feat_groups->{$trans_id->{'id'}} = $doc->createAndAddFeatureGroup($seq, '', $t_id);
            
            # create reference for parent of transcript also
            my $trans = $features->{$trans_id->{'type'}}->{$trans_id->{'id'}};
            foreach my $parent_id (@{$trans->{'parent'}}) {
                $feat_groups->{$parent_id} = $feat_groups->{$trans_id->{'id'}};
            }
        }
    }

    foreach my $type(keys(%{$features})) {

        foreach my $feat_id(keys(%{$features->{$type}})) {
                my $id;
                if (! defined($id_hash->{$feat_id})) {
                    $id = $idcreator->next_id( 
                                               'project' => $options->{'project'}, 
                                               'type' => $type 
                                               );
                } else {
                    $id = $id_hash->{$feat_id};
                }

                my $feat_title = ($feat_id =~ /temp_id/) ? undef : $feat_id;
                my $feat = $doc->createAndAddFeature( $feat_table, $id, $feat_title, $type);

                # link polypeptides to Sequences in --peptide_fasta
                if (($type eq 'polypeptide') && defined($peptides)) {
                    my $fseq = $peptides->{$feat_id};
                    die "couldn't find peptide sequence for $feat_id" if (!defined($fseq));
                    my $pep_seq = $doc->createAndAddSequence($id . "_seq", undef, $fseq->{'seqlen'}, 'aa', 'polypeptide');
                    my $seq_genome_link = $doc->createAndAddLink($pep_seq, 'genome', '#' . $genome_id);
                    my $feat_seq_link = $doc->createAndAddLink($feat, 'sequence', '#'.$id.'_seq');

                    my $path = File::Spec->rel2abs($options->{'peptide_fasta'});
                    my $seq_data_import_elem = $doc->createAndAddSeqDataImport($pep_seq, 'fasta', $path, undef, $fseq->{'defline'});
                }

                # TODO - need a way to specify which features get linked to an analysis
#                $feat->addBsmlLink('analysis', '#EVM', 'computed_by');

                $feat->addBsmlIntervalLoc(
                                          $features->{$type}->{$feat_id}->{'startpos'},
                                          $features->{$type}->{$feat_id}->{'endpos'},
                                          $features->{$type}->{$feat_id}->{'complement'},
                                          );

                # determine which feat_group this feature belongs with
                my $feat_group = undef;
                if (defined($feat_groups->{$feat_id})) { 
                    $feat_group = $feat_groups->{$feat_id};
                } else {
                    my $trans_id = &get_parent_id_by_type($features->{$type}->{$feat_id}, $features, $MRNA_TYPES);
                    die "no parent id found for $type feature" if (!defined($trans_id));
                    $feat_group = $feat_groups->{$trans_id};
                    if (!defined($feat_group)) {
                        print Dumper $features;
                        die "no feat_group found for transcript=$trans_id";
                    }
                }
                
                $feat_group->addBsmlFeatureGroupMember( $id, $type );
            }
    }
}

# traverses an array of node references 
# and returns all nodes of the specified type
sub fetch_node_type {
    my ($type, $nodes_ref, $found_nodes) = @_;
    
    foreach my $node (@{$nodes_ref}) {
        if ($node->{'_type'} eq $type) {
            push (@{$found_nodes}, $node);
        }
        foreach my $key (keys %{$node->{'children'}}) {
            fetch_node_type($type, $node->{'children'}->{$key}, $found_nodes);
        }
    }
    
    return;
}

# read deflines and sequence lengths from a multi-FASTA file
#
# $file - 
# $id_regex - 
#
# returns a hashref mapping sequence id to { 'seqlen' => sequence length, 'defline' => sequence defline }
#
sub read_sequence_lengths {
    my($file, $id_regex, $id_fn) = @_;
    my $result = {};
    my $fh = FileHandle->new();

    my $process_seq = sub {
        my($defline, $seq, $lnum) = @_;
        my($id) = ($defline =~ /$id_regex/);
        die "unable to parse id from line $lnum of $file: $defline" if (!defined($id));
        $seq =~ s/\s+//g;
        $defline =~ s/^>//;
        $id = &$id_fn($id) if (defined($id_fn));
        die "sequence id '$id' is not unique in $file" if (defined($result->{$id}));
        $result->{$id} = { 'id' => $id, 'seqlen' => length($seq), 'defline' => $defline };
    };

    my $defline = undef;
    my $defline_lnum = undef;
    my $seq = undef;
    my $lnum = 0;

    $fh->open($file, 'r') || die "unable to read from $file";
    while (my $line = <$fh>) {
        chomp($line);
        ++$lnum;
        if ($line =~ /^>/) {
            &$process_seq($defline, $seq, $defline_lnum) if (defined($defline));
            $defline = $line;
            $defline_lnum = $lnum;
            $seq = "";
        } else {
            $seq .= $line;
        }
    }
    &$process_seq($defline, $seq, $defline_lnum) if (defined($defline));
    $fh->close();

    return $result;
}

sub getTempId {
    return "temp_id_".$GLOBAL_ID_COUNTER++; 
}
