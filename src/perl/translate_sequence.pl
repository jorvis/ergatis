#!/usr/local/bin/perl

use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";
	
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Workflow::IdGenerator;
use XML::Twig;
use BSML::BsmlDoc;

my $transeq_exec;
my $transeq_flags = ' -warning 1 -error 1 -fatal 1 -die 1';

my %coords = ();
my %id2title = ();
my @genes = ();

my $bsml_sequences;
my $sequence_children;
my $exon_locs;
my $cds_locs;
my $cds_regions;
my $exon_frame;
my $cds_frame;

my %options = ();
my $results = GetOptions (\%options,
						  'transeq_bin=s',
                          'input|i=s',
						  'class|c:s',
						  'regions|r:s',
						  'frame|f:s',
						  'table|t:i',
						  'output|o=s',
						  'id_repository=s',
						  'project=s',
						  'substring',
                          'help|h') || pod2usage();

my $fasta_flag = 0;
					 
if ($options{'help'}) {
	pod2usage(verbose => 2);
}

if (!$options{'id_repository'}) {
	pod2usage("must provided --id_repository");
}
if (!$options{'project'}) {
	pod2usage("project name must be provided with --project");
}
my $id_gen = Workflow::IdGenerator->new('id_repository' => $options{'id_repository'});

if (!$options{'transeq_bin'}) {
	pod2usage("must provide path to transeq executable with --transeq_bin");
}
$transeq_exec = $options{'transeq_bin'};

unless ($options{'input'}) {
	pod2usage("fasta or bsml input must be provided with --input");
}

if ($options{'input'}) {

	unless (-e $options{'input'}) {
		pod2usage("input file '$options{input}' does not exist");
	}
	if ($options{'input'} =~ /\.fsa$/) {
		$fasta_flag = 1;
	}
	if ($fasta_flag && count_fasta_records($options{'input'}) != 1) {
		die "fasta input file must contain only one sequence record";
	}
}

if (!$options{'output'}) {
	pod2usage("must specify and output directory with --output");
}
$options{'output'} =~ s/\/$//;
unless (-d $options{'output'}) {
	pod2usage("specified output path '$options{output}' is not a directory");
}



if ($options{'table'}) {
	## translation table
	## see transeq -h for list
	$transeq_flags .= " -table $options{table}";
}


if (!$fasta_flag) {
	## if a bsml document has been provided
	## we're going to process it and perform
	## transeq on the gene models it encodes

	if ($options{'regions'}) { 
		print STDERR "WARNING: --regions flag is incompatible with BSML document input and will be ignored\n";
	}	
	if ($options{'frame'}) { 
		print STDERR "WARNING: --frame flag is incompatible with BSML document input and will be ignored\n";
	}	
	
	## scan through BSML input doc and find all assembly sequences
	## make sure sequence is present as Seq-data-import or Seq-data
	parse_bsml_sequences($options{'input'});

	## scan through BSML and pull out interval locs from feature tables	
	parse_bsml_interval_locs($options{'input'});
	
	## prepare a hash of strings for the transeq regions flag
	## for each transcript
	foreach my $transcript_id(keys(%{$exon_locs})) {
		$exon_locs->{$transcript_id} = constrain_exons_by_cds(
										$exon_locs->{$transcript_id}, 
										$cds_locs->{$transcript_id}->[0],
										$cds_locs->{$transcript_id}->[1],
							  								 );
		my @exon_ref_arr = @{$exon_locs->{$transcript_id}};
		my @cds_regions = ();
		foreach my $exon_ref(@exon_ref_arr) {
			push(@cds_regions, $exon_ref->[0]."-".$exon_ref->[1]);
		}
		$cds_regions->{$transcript_id} = join(",",@cds_regions);
	}

	my $temp_in_fsa = $options{'output'}."/temp.in.fsa";
	my $temp_out_fsa = $options{'output'}."/temp.out.fsa";
	foreach my $seq_id(keys(%{$sequence_children})) {
		my $seq = '';
		if ($options{'substring'}) {
			$seq = get_sequence($bsml_sequences->{$seq_id});
		}
		foreach my $transcript_id(@{$sequence_children->{$seq_id}}) {
			my $flags = $transeq_flags;
			my $nt_fasta = $bsml_sequences->{$seq_id};
			my $regions_string = $cds_regions->{$transcript_id};
			if ($options{'substring'}) {
				$nt_fasta = $temp_in_fsa;
		   		my $subseq = get_regions_substring($seq, $regions_string);
				write_seq_to_fasta($nt_fasta, "$seq_id", $subseq);
			} else {
				$flags .= " -regions $regions_string";
			}
			$flags .= " -sequence $nt_fasta";
			$flags .= " -outseq $temp_out_fsa";
			$flags .= " -frame $cds_frame->{$transcript_id}";
			system($transeq_exec.$flags);
			my $out_fsa = $options{'output'}."/"."$transcript_id.fsa";
			my $id_hash = replace_sequence_ids($temp_out_fsa, $out_fsa);
			unlink($temp_out_fsa);
		}
	}
	if (-e $temp_in_fsa) {unlink($temp_in_fsa);}
	
} else {
	## otherwise we're just going to run
	## transeq on the input nt sequence
	## to generate a polypeptide fasta file
	
	my $temp_out_fsa = $options{'output'}."/"."temp.polypeptides.fsa";

	if ($options{'regions'}) {
		$transeq_flags .= " -regions $options{regions}";
	}
	if (!$options{'frame'}) {
		$options{'frame'} = '1';
	} 
	$transeq_flags .= " -frame $options{frame}";
	$transeq_flags .= " -sequence $options{input}";
	$transeq_flags .= " -outseq $temp_out_fsa";
	
	system($transeq_exec.$transeq_flags);

	my $query_id = get_sequence_id($options{'input'});
	
	my $out_fsa = $options{'output'}."/"."$query_id.fsa";
	my $id_hash = replace_sequence_ids($temp_out_fsa, $out_fsa);
	
	unlink($temp_out_fsa);
}

exit();

sub replace_sequence_ids {
	my ($old_fsa_file,$new_fsa_file) = @_;
	
	my %ids = ();
	
	my $count = count_fasta_records($old_fsa_file);
	
	$id_gen->set_pool_size('polypeptide' => $count);
	
	open (IN, $old_fsa_file) || die "couldn't read fsa file";
	open (OUT, ">".$new_fsa_file) || die "couldn't write fsa file";
	while (<IN>) {
		if (/^>[^\s]+_(\d)/) {
			my $seq_id = $id_gen->next_id(
    	    		                      'project' => $options{'project'},
        	        		              'type'    => 'polypeptide'
                            		     );
			$ids{$1} = $seq_id;
			print OUT ">$seq_id\n";
		} else {
			print OUT $_;
		}
	}
	close IN;
	close OUT;

	return \%ids;
}

sub count_fasta_records {
	my ($fname) = @_;
	my $count = 0;
	open (IN, $fname) || die "couldn't open fasta file for reading";
	while (<IN>) {
		if (/^>/) {
			$count++;
		}
	}
	close IN;

	return $count;
}

sub get_sequence_id {
	my ($fname) = @_;
	my $id = '';
	open (IN, $fname) || die "couldn't open fasta file for reading";
	while (<IN>) {
		chomp;
		if (/^>([^\s]+)/) {
			$id = $1;
			last;	
		}
	}
	close IN;

	return $id;
}


sub parse_bsml_sequences {
        my ($file) = @_;
		
		my $ifh;
        if (-e $file.".gz") {
			$file .= ".gz";
		} elsif (-e $file.".gzip") {
			$file .= ".gzip";
		}

        if ($file =~ /\.(gz|gzip)$/) {
            open ($ifh, "<:gzip", $file);
        } else {
            open ($ifh, "<$file");
        }
		
		my $twig = new XML::Twig(	TwigRoots => {'Sequence' => 1, 'Feature-group' => 1},
									TwigHandlers => {'Sequence' => \&process_sequence}
								);
		
		$twig->parse($ifh);
		close $ifh;
}

sub process_sequence {
	my ($twig, $sequence) = @_;

	my $seq_id;
	my $class;
	my $molecule;
	my $seq_data_import;
	my $seq_data;
	my $source;
	my $identifier;
	my $format;
	
	$seq_id = $sequence->{'att'}->{'id'};
	$class = $sequence->{'att'}->{'class'};
	$molecule = $sequence->{'att'}->{'molecule'};
		
	if (!defined($class)) {
		print STDERR "WARNING: sequence class of '$seq_id' was not defined\n";
	}
	if (!defined($molecule) || $molecule eq 'mol-not-set') {
		print STDERR "WARNING: molecule type of '$seq_id' was not defined\n";
	}
	
	if ($seq_data_import = $sequence->first_child('Seq-data-import')) {
		
		$source = $seq_data_import->{'att'}->{'source'};
		$identifier = $seq_data_import->{'att'}->{'identifier'};
		$format = $seq_data_import->{'att'}->{'format'};
		
		unless (-e $source) {
			die "fasta file referenced in BSML Seq-data-import '$source' doesn't exist";
		}
		unless (defined($identifier)) {
			die "Seq-data-import for '$seq_id' does not have a value for identifier";
		}
		
		my $seq_count = count_fasta_records($source);
		if ($seq_count == 1) {
			unless (get_sequence_id($source) eq $identifier) {
				print STDERR "ID disagreement between BSML Seq-data-import '$identifier' and fasta file '$source'";
			}

			$bsml_sequences->{$seq_id} = $source;
			
		} elsif ($seq_count > 1) {
			my $nt_seq = get_sequence_by_id($source, $identifier);
			
			if (length($nt_seq) > 0) {
				my $sequence_file = $options{'output'}."/$seq_id.fsa";
				write_seq_to_fasta($sequence_file, $seq_id, $nt_seq);

				$bsml_sequences->{$seq_id} = $sequence_file;
				
			} else {
				die "couldn't extract sequence for '$seq_id' from Seq-data-import source '$source'";
			}
		} else {
			die "no fasta records found for BSML Seq-data-import '$source' for sequence '$seq_id'";
		}
			
		
	} elsif ($seq_data = $sequence->first_child('Seq-data')) {
		## sequence is in the BSML
		## so it will be written to a fasta file in the output dir
		
		my $sequence_file = $options{'output'}."/$seq_id.fsa";
	
		write_seq_to_fasta($sequence_file, $seq_id, $seq_data->text());
		
		$bsml_sequences->{$seq_id} = $sequence_file;
		
	} else {
		
		## there is no Seq-data or Seq-data-import for the sequence
		die "No sequence present in BSML sequence element";
	}
   
	foreach my $child ($sequence->first_child('Feature-tables')->children('Feature-group')) {
		#$feature_parent_seq->{$child->{'att'}->{'group-set'}} = $seq_id;
		push(@{$sequence_children->{$seq_id}}, $child->{'att'}->{'group-set'});
	}
		
    $twig->purge;
}

sub parse_bsml_interval_locs {
        my ($file) = @_;
		
		my $ifh;

        if (-e $file.".gz") {
			$file .= ".gz";
		} elsif (-e $file.".gzip") {
			$file .= ".gzip";
		}

        if ($file =~ /\.(gz|gzip)$/) {
            open ($ifh, "<:gzip", $file);
        } else {
            open ($ifh, "<$file");
        }

		my $twig = new XML::Twig(
								 twig_roots => {
										'Feature' 			 => \&process_feat,
                                 		'Feature-group' 	 => \&process_feat_group,
                                			   }
								);
		$twig->parse($ifh);
}


sub process_feat {
	my ($twig, $feat) = @_;
    my $id = $feat->att('id');

    if ($feat->att('class') eq 'exon') {
    	my $seq_int = $feat->first_child('Interval-loc');
    	my $complement = $seq_int->att('complement');
        my ($start_pos, $end_pos) = ($seq_int->att('startpos') + 1, $seq_int->att('endpos') + 1);
        push @{$coords{$id}}, $start_pos, $end_pos;
       	if ($complement) {
			$exon_frame->{$id} = -1;	
		} else {
			$exon_frame->{$id} = 1;
		}
    } elsif ($feat->att('class') eq 'CDS') {
    	my $seq_int = $feat->first_child('Interval-loc');
        my ($start_pos, $end_pos) = ($seq_int->att('startpos') + 1, $seq_int->att('endpos') + 1);
        push @{$coords{$id}}, $start_pos, $end_pos;
    }

    $twig->purge;
}

sub process_feat_group {
    my ($twig, $feat_group) = @_;
    my @exon_coords = ();

	my $feat_group_id = $feat_group->{'att'}->{'group-set'};
	my $count = 0;
	my $sum = 0;
    foreach my $child ($feat_group->children('Feature-group-member')) {
		if ($child->att('feature-type') eq 'exon') {
        	my $id = $child->att('featref');
            push(@exon_coords,[$coords{$id}->[0], $coords{$id}->[1]]);
			$sum += $exon_frame->{$id};
			$count++;
        } elsif ($child->att('feature-type') eq 'CDS') {
        	my $id = $child->att('featref');
			$cds_locs->{$feat_group_id} = [$coords{$id}->[0], $coords{$id}->[1]];
		}
    }
    if (scalar @exon_coords) {
    	@exon_coords = sort { $$a[0] <=> $$b[0]; } @exon_coords;
		$exon_locs->{$feat_group_id} = \@exon_coords;
		if (abs($sum/$count) != 1) {
			die "transcript '$feat_group_id' has some exons on both strands";
		} else {
			$cds_frame->{$feat_group_id} = $sum/$count;
		}
	}
	$twig->purge;
}


sub write_seq_to_fasta {
	my ($file, $header, $sequence) = @_;
	
	open (OUT, ">$file") || die "couldn't write fasta file '$file'";
		
	$sequence =~ s/\W+//g;
	$sequence =~ s/(.{1,60})/$1\n/g;
		
	print OUT ">$header\n$sequence";
	close OUT;
}

sub get_sequence_by_id {
	my ($fname, $id) = @_;
	my $seq_id = '';
	my $sequence = '';
	open (IN, $fname) || die "couldn't open fasta file for reading";
	TOP: while (<IN>) {
		chomp;
		if (/^>([^\s]+)/) {
			$seq_id = $1;
			if ($seq_id eq $id) {
				while (<IN>) {
					chomp;
					if (/^>/) {
						last TOP;
					} else {
						$sequence .= $_;
					}
				}
			}	
		}
	}
	close IN;

	return $sequence;
}


## return the sequence 
sub get_sequence {
	my ($fname) = @_; 
	
	my $sequence = '';
	my $flag = 0;
	
	open (IN, $fname) || die "couldn't open fasta file for reading";
	while (<IN>) {
		chomp;
		if (/^>/) {
			if ($flag) {
				die "Unexpectedly encountered more than one fasta record in input nt sequence file";
			}
			$flag = 1;
			next;
		} else {
			$sequence .= $_;
		}
	}
	$sequence =~ s/\W+//g;
	return $sequence;
}

## takes the string that would otherwise be provided to transeq
## with the regions flag and uses it to return the substring of
## an input sequence
sub get_regions_substring {
	my ($sequence, $regions_string)=@_;
	
	my @regions = split(",", $regions_string);

	my $subseq = '';
	
	foreach my $region(@regions) {
		my ($start, $stop) = split("-", $region);
		my $len = $stop - $start + 1;
		$start--;
		$subseq .= substr($sequence, $start, $len);
	}
	
	return $subseq;
}

## constrains the exons regions to within boundaries
## set by the CDS feature 
sub constrain_exons_by_cds {
	my ($exon_loc_ref, $cds_start, $cds_end) = @_;
	
	my @exon_locs = ();
	
	foreach my $exon_ref(@{$exon_loc_ref}) {
		if ($exon_ref->[1] < $cds_start) {
			next;
		}
		if ($exon_ref->[0] > $cds_end) {
			next;
		}
		if ($exon_ref->[0] < $cds_start) {
			$exon_ref->[0] = $cds_start;
		}
		if ($exon_ref->[1] > $cds_end) {
			$exon_ref->[1] = $cds_end;
		}
		push(@exon_locs, $exon_ref);
	}
	
	return \@exon_locs;
}
