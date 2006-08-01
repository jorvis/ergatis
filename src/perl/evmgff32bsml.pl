#!/usr/local/bin/perl

=head1 NAME

evmgff32bsml.pl

=head1 SYNOPSIS

USAGE: evmgff32bsml.pl
           --input|-i
           --output|-o
           --project|-p
           --help|-h

=head1 OPTIONS

B<--input,-i>
    Input gff3 format file created by EVM.

B<--output,-o>
    Output BSML file name.

B<--project,-p>
    Project/database name used to create feature identifiers.

B<--help,-h>
    This help documentation

=head1 DESCRIPTION

This script parses GFF output from EVM and writes it out 
as BSML suitable for import into CHADO DBs.

=head1 INPUT

Current assumption is that input files are syntactically
correct with respect to the GFF3 spec, and that they are
encoding predicted gene features.

The script could be extended to support writing BSML from
GFF3 compliant files encoding other features.

=head1 OUTPUT

Output will be BSML encoding the EVM predicted gene models.

=head1 CONTACT

    Brett Whitty
    bwhitty@tigr.org

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Data::Dumper;
use BSML::BsmlBuilder;
use Workflow::IdGenerator;

my ($input, $help, $man, $output, $project, $id_repository);

GetOptions (
             'help|h'              => \$help,
             'man|m'               => \$man,
             'output|o=s'          => \$output,
             'input|i=s'           => \$input,
             'project|p=s'         => \$project,
			 'id_repository=s'	   => \$id_repository
           ) || pod2usage();

&pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if ($man);
&pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT}) if ($help);


if (!$input){
    pod2usage("input gff_file was not defined with --input");
}

if (!$output){
    pod2usage("output BSML filename was not defined with --output");
}

if (!$project) {
        pod2usage("You must specify a project name with --project");
}

if (!$id_repository) {
        pod2usage("You must specify an id repository with --id_repository");
}

$project =~ tr/A-Z/a-z/;

my $idcreator = Workflow::IdGenerator->new('id_repository' => $id_repository);

my $global_id_counter=0;
my $nodes = {};
my @root_nodes = ();

open (IN, $input) || die "couldn't open input file '$input'";

while (<IN>) {
	chomp;
	
	if (/^#/) {
		#process header line
	} elsif (/^>/ || /^#.*fasta.*/i) {
		last;
	} else {
		my $record = parse_record($_);

		## if it has no parents it's a root node	
		if (!$record->{'Parent'}) {
			push(@root_nodes, $record);
		}
	
		## store records in a set of arrays indexed by ID
		if ($record->{'ID'}) {
			push(@{$nodes->{$record->{'ID'}}->{'records'}},$record);
		}
	}
	
}
close IN;

## populate children for each record
foreach my $id(keys(%{$nodes})) {
	foreach my $record(@{$nodes->{$id}->{'records'}}) {
		$record->{'children'} = $nodes->{$id}->{'children'};
		delete $nodes->{$id}->{'children'};
	}
}

my $gene_nodes = [];
fetch_node_type('gene', \@root_nodes, $gene_nodes);

my $doc = new BSML::BsmlBuilder();

#foreach my $root(@root_nodes) {
foreach my $root(@{$gene_nodes}) {
	if ($root->{'_type'} eq 'gene') {
		my $features = {};
		process_node($root, $features);
		$features->{'_seqid'} = $root->{'_seqid'};
		gene_feature_hash_to_bsml($features);
	}
}

# add the analysis element
$doc->createAndAddAnalysis(
                            id => 'EVM_analysis',
                           	sourcename => $output,
                          );

$doc->write($output);

exit();

sub parse_record {
	my ($line) = @_;

	my $record = {};
	my $attrib_hash = {};
	my $record_id = '';
	my $parent_id = '';
	
	my @cols = split("\t");

    ## adjust both positions so that we are numbering from zero
    $cols[3]--;
    $cols[4]--;

    ## change the + and - symbols in strand column to 0 and 1, respectively
    if ($cols[6] eq '+') {
        $cols[6] = 0;
    } elsif ($cols[6] eq '-') {
        $cols[6] = 1;
    } else {
        die("unknown value ($cols[6]) in strand column.  expected + or -.");
    }

	$record->{'_seqid'}		= $cols[0];
	$record->{'_source'}	= $cols[1];
	$record->{'_type'} 		= $cols[2];
	$record->{'_start'}		= $cols[3];
	$record->{'_end'}		= $cols[4];
	$record->{'_score'}		= $cols[5];
	$record->{'_strand'}	= $cols[6];
	$record->{'_phase'}		= $cols[7];
	
	my @attribs = split(";", $cols[8]);
	foreach my $attrib(@attribs) {
		my ($type, $val) = split("=", $attrib);
		my @vals = split(",", $val);
		$record->{$type}=\@vals;
	}

	if (!defined($record->{'ID'})) {
		#print STDERR "no ID defined for record\n";
	} else {
		$record->{'ID'} = $record->{'ID'}->[0];
	}

	if (!defined($record->{'Parent'})) {
		#print STDERR "no Parent attribute defined for record\n";
	} else {
		foreach my $parent(@{$record->{'Parent'}}) {
			## store record's hash reference as a child of Parent
			if (!defined($nodes->{$parent}->{'children'})) {
				$nodes->{$parent}->{'children'}->{$record->{'_type'}} = [];
			}
			push (@{$nodes->{$parent}->{'children'}->{$record->{'_type'}}},$record);
		}
	}

	return $record;
}


## recursively traverse a node and all its children
## and process each record to create a feature hash
sub process_node {

	my ($node, $features) = @_;

	## hash of keys to ignore when extracting a record hash from a node hash
	my $ignore_keys = {
						'children' => 1,
				      };
	
	## build a record hash from values stored in the node hash
	my $record;
	foreach my $key(keys %{$node}) {
		if (!$ignore_keys->{$key}) {
			$record->{$key} = $node->{$key};
		}
	}
	
	## process the record hash
	process_record($record, $features);
	
	## process the node's children
	foreach my $child_type(keys %{$node->{'children'}}) {
		foreach my $child_record(@{$node->{'children'}->{$child_type}}) {
				process_node($child_record, $features); 
		}
	}
	
	return;	
}

## process the records and store them in the features hash
sub process_record {
	my ($record, $features) = @_;
	
	my $feat_type_map = {
							'gene'	=>	'gene',
							'CDS'	=>	'CDS',
							'exon'	=>	'exon',
							'mRNA'	=>	'transcript',				
						};
	
	
	if ($record->{'_type'} eq 'CDS') {
		## CDS records can span lines and must be merged into one CDS feature
		if (!$record->{'ID'}) {
			die "CDS feature lacks ID -> bad form!";
		}
		if ($features->{'CDS'}->{$record->{'ID'}}) {
			if ($features->{'CDS'}->{$record->{'ID'}}->{'startpos'} > $record->{'_start'}) {
				$features->{'CDS'}->{$record->{'ID'}}->{'startpos'} = $record->{'_start'};
			}
			if ($features->{'CDS'}->{$record->{'ID'}}->{'endpos'} < $record->{'_end'}) {
				$features->{'CDS'}->{$record->{'ID'}}->{'endpos'} = $record->{'_end'};
			}
		} else {
			$features->{$feat_type_map->{'CDS'}}->{$record->{'ID'}} = {
													'complement' 	=> $record->{'_strand'},
													'startpos'		=> $record->{'_start'},
													'endpos'		=> $record->{'_end'},
																	  };
		}
		
	} else {
		##handle all other feature types	
		
		my $feat_type;
		
		if (!defined($feat_type_map->{$record->{'_type'}})) {
			print STDERR "unexpected feature type '$record->{_type}'\n";
			$feat_type = $record->{'_type'};	
		} else {
			$feat_type = $feat_type_map->{$record->{'_type'}};
		}

		my $id;
		if (!$record->{'ID'}) {
			$id = getTempId();
		} else {
			$id = $record->{'ID'};
		}
		
		my $title = undef;
		if ($record->{'Name'}) {
			$title = shift(@{$record->{'Name'}});
		}
		
		$features->{$feat_type}->{$id} = {	
					'complement' 	=> $record->{'_strand'},
					'startpos'		=> $record->{'_start'},
					'endpos'		=> $record->{'_end'},
					'title'			=> $title,
									     };
	}
	
}


## convert a gene feature hash into BSML
sub gene_feature_hash_to_bsml {
	my ($features) = @_;
	
	my $seq_id = $features->{'_seqid'};
	delete $features->{'_seqid'};

	my $id_hash = {};
	my %id_counts;
	
	foreach my $type(keys(%{$features})) {
		$id_counts{$type} = scalar(keys(%{$features->{$type}}));
	}
	$idcreator->set_pool_size(%id_counts);
	
	my $seq;
	
	## create a sequence stub for the seq_id if it doesn't exist yet
	if (!($doc->returnBsmlSequenceByIDR($seq_id))){
        $seq = $doc->createAndAddSequence( $seq_id, $seq_id, undef, 'dna', '' );
		#if ($options{'fasta_file'}) {
		#   $doc->createAndAddSeqDataImport($seq, 'fasta', $options{'fasta_file'}, '', "$features->{_seqid}");
		#}
        $seq->addBsmlLink('analysis', '#EVM', 'input_of');
    } else {
		$seq = $doc->returnBsmlSequenceByIDR($seq_id);
	}

	my @transcript_id = keys(%{$features->{'transcript'}});

	if (scalar @transcript_id > 1) {
		print Dumper $features;
		die "multiple transcripts encountered";
	}
 	my $t_id = $idcreator->next_id( 
									'project' => $project, 
									'type' => 'transcript' 
						 		  );
	
	$id_hash->{$transcript_id[0]} = $t_id;
	
	my $feat_table = $doc->createAndAddFeatureTable($seq);
	my $feat_group = $doc->createAndAddFeatureGroup($seq, '', $t_id);
    
	foreach my $type(keys(%{$features})) {
		
		foreach my $feat_id(keys(%{$features->{$type}})) {
    
			my $id;
			if (! defined($id_hash->{$feat_id})) {
		   		$id = $idcreator->next_id( 
										'project' => $project, 
										'type' => $type 
					 				);
			} else {
				$id = $id_hash->{$feat_id};
			}
			
			my $feat = $doc->createAndAddFeature( $feat_table, $id, $feat_id, $type);

			$feat->addBsmlLink('analysis', '#EVM', 'computed_by');

			$feat->addBsmlIntervalLoc(
					$features->{$type}->{$feat_id}->{'startpos'},
					$features->{$type}->{$feat_id}->{'endpos'},
					$features->{$type}->{$feat_id}->{'complement'},
									 );
   
			$feat_group->addBsmlFeatureGroupMember( $id, $type );
		}
		
	}
	
}


## traverses an array of node references 
## and returns all nodes of the specified type
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

sub getTempId {
	return "temp_id_".$global_id_counter++;	
}
