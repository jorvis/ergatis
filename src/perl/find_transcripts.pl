#!/usr/bin/perl

=head1 NAME

find_transcripts.pl - finds transcripts from an alignment based on coverage islands.  

=head1 SYNOPSIS

 USAGE: find_transcripts.pl
	--file_type either BAM, WIG, or BigWig
	--file the BAM, WIG, or BigWig file.  
	       If stranded and using WIG or BigWig, a comma-separated pair, like <forward coverage>, <reverse coverage>
	--sizes_file needed for a WIG filetype, to indicate chromosome sizes
	--stranded= 1 or 0, indicating if the RNA-seq library was strand-specific [default 0].  
	            if stranded and using WIG or BigWig files, one must be provided for each strand
	--reverse_strand = 1 or 0, to reverse the strand of ntars, useful for wierd RNA-seq protocols [default 0]

	--min_cov = minimum coverage for a transcribed region [default 10]
	--max_intron = maximum length of an "intron", or low coverage area inside a transcript [default 50]
	--min_transcript = minimum length of a transcript [default 0]
	--min_intergenic = minimum spacing between NTARs and base (reference) genes [default 250]
        --min_ntar_intergenic = minimum spacing between successive NTARs [default 1]	
	--base_gff3= optional gff3 containing known transcripts to use for filtering for novels
	--base_feature_type= feature type to extract from base_gff3 for transcripts [default transcript]
	
	--output_stub=path and filename stub to use for output files
	--verbose output logging information [default no]
	--ntar_prefix optional prefix to add to ntar ids
    [
       --help
     ]

=head1  CONTACT

    Umar Farooq
    ufarooq@som.umaryland.edu

=cut

use strict;
use warnings;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use NGS::GFF3;
use NGS::Alignment;

my %options = ();
my $results = GetOptions (\%options, 
						  'file_type=s',
						  'file=s',
						  'sizes_file=s',
						  'stranded=s',
						  'reverse_strand=s',
						  'min_cov=s',
						  'max_intron=s',
						  'min_transcript=s',
						  'base_gff3=s',
						  'base_feature_type=s',
						  'min_intergenic=s',
			                          'min_ntar_intergenic=s',
						  'output_stub=s',
						  'verbose=s',
						  'ntar_prefix=s',
                          'help|h') || pod2usage();

# display documentation
if( $options{'help'} ) {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

# make sure parameters are correct
&check_parameters(\%options);

#a hash with start,stop, chromosome, and strand to store new transcripts
my %ntars; 

#this will assign unique gene ids
my $gene_id = 1;

#the total number of islands found
my $islands = 0;

if(!$options{stranded}){
    my $alignment;
    
    if(uc($options{file_type}) eq "BAM"){
		$alignment = NGS::Alignment->new( stranded => 0, file_type => "BAM", file => $options{file});
    }
    elsif(uc($options{file_type}) eq "WIG"){	
		$alignment = NGS::Alignment->new( stranded => 0, file_type => "WIG", file => $options{file}, sizes_file => $options{sizes_file});
    }
    elsif(uc($options{file_type}) eq "BIGWIG"){	
		$alignment = NGS::Alignment->new( stranded => 0, file_type => "BigWig", file => $options{file});
    }

    my $seq_ids = $alignment->seq_ids;

    foreach my $seq_id (@{$seq_ids}){
		print "Processing seq_id: $seq_id\n";
		
		my $coverage = $alignment->get_per_bp_coverage({seq_id => $seq_id});
		
		my $novels = find_transcripts($options{min_cov}, 
									  $options{max_intron}, 
									  $options{min_transcript}, 
									  $options{min_ntar_intergenic},
									  $coverage);
		foreach (sort keys %{$novels}){
			$ntars{$_."_".$options{ntar_prefix}} = {seq_id => $seq_id, start => $novels->{$_}{start}, end => $novels->{$_}{end}, strand => "."};
		}
    }
	print "Done.\n\n";
	
} else {

    # for stranded libraries, we treat each strand separately
    if( uc($options{file_type}) eq "BAM"){
		
		if( $options{verbose} ) {
			print "\nProcessing stranded BAM files...\n\n";
		}
		
		my $alignment = NGS::Alignment->new( stranded => 1, file_type => "BAM", file => $options{file});
		
		my $seq_ids = $alignment->seq_ids;
		
		# look for ntars on the plus strand
		
		foreach my $seq_id (@{$seq_ids}){
			print "Processing seq_id: $seq_id for + strand\n";
			
			my $coverage = $alignment->get_per_bp_coverage({seq_id => $seq_id, strand => "+"});
			
			my $novels = find_transcripts($options{min_cov}, 
										  $options{max_intron}, 
										  $options{min_transcript}, 
										  $options{min_ntar_intergenic},
										  $coverage);

			foreach (sort keys %{$novels}) {
				my $ntar_id = $_ . "_" . $options{ntar_prefix};
				$ntars{$seq_id}->{'ntar_id'}->{$ntar_id} = {start => $novels->{$_}->{start}, end => $novels->{$_}->{end}, strand => "+"};
			}
		}
		
		#look for ntars on the minus strand
		foreach my $seq_id (@{$seq_ids}){
			print "\nProcessing seq_id: $seq_id for - strand";
			
			my $coverage = $alignment->get_per_bp_coverage({seq_id => $seq_id, strand => "-"});
			
			my $novels = find_transcripts($options{min_cov}, 
										  $options{max_intron}, 
										  $options{min_transcript}, 
										  $options{min_ntar_intergenic},
										  $coverage);

			foreach (sort keys %{$novels}) {
				my $ntar_id = $_ . "_" . $options{ntar_prefix};
				$ntars{$seq_id}->{'ntar_id'}->{$ntar_id} = {start => $novels->{$_}->{start}, end => $novels->{$_}->{end}, strand => "-"};
			}
		}
		
    } elsif( (uc($options{file_type}) eq "WIG") || (uc($options{file_type}) eq "BIGWIG") ) {
		
		if( $options{verbose} ) {
			print "\nProcessing WIG files...\n\n";
		}

		# processing for WIG and BigWig file types
		my $p_alignment;
		my $m_alignment;
		
		if( uc($options{file_type}) eq "WIG" ) {
			$p_alignment = NGS::Alignment->new( stranded => 1, file_type => "WIG", file => $options{forward_coverage}, sizes_file => $options{sizes_file});
			$m_alignment = NGS::Alignment->new( stranded => 1, file_type => "WIG", file => $options{reverse_coverage}, sizes_file => $options{sizes_file});
			
		} elsif( uc($options{file_type}) eq "BIGWIG"){
			$p_alignment = NGS::Alignment->new( stranded => 1, file_type => "BigWig", file => $options{forward_coverage} );
			$m_alignment = NGS::Alignment->new( stranded => 1, file_type => "BigWig", file => $options{reverse_coverage} );
		}
		
		my $p_seq_ids = $p_alignment->seq_ids;
		my $p_sizes = $p_alignment->sizes();
		
		# look for ntars on the plus strand
		foreach my $seq_id (@{$p_seq_ids}) {
			print "Processing seq_id: $seq_id for + strand\n";
			
			my $coverage = $p_alignment->get_per_bp_coverage({seq_id => $seq_id, strand => "+"});
			
			my $novels = find_transcripts($options{min_cov}, 
										  $options{max_intron}, 
										  $options{min_transcript}, 
										  $options{min_ntar_intergenic},
										  $coverage);
			
			# add forward novel tars to ntars hash
			foreach (sort keys %{$novels}) {
				my $ntar_id = $_ . "_" . $options{ntar_prefix};
				$ntars{$seq_id}->{'ntar_id'}->{$ntar_id} = {start => $novels->{$_}->{start}, end => $novels->{$_}->{end}, strand => "+"};
			}
		}
		
		my $m_seq_ids = $m_alignment->seq_ids;
		my $m_sizes = $m_alignment->sizes();
		
		# look for ntars on the minus strand
		foreach my $seq_id (@{$m_seq_ids}) {
			print "Processing seq_id: $seq_id for - strand\n";
			
			my $coverage = $m_alignment->get_per_bp_coverage({seq_id => $seq_id, strand => "-"});

			my $novels = find_transcripts($options{min_cov}, 
										  $options{max_intron}, 
										  $options{min_transcript}, 
										  $options{min_ntar_intergenic},
										  $coverage);
			
			# add reverse novel tars to global ntars hash
			foreach (sort keys %{$novels}) {
				my $ntar_id = $_ . "_" . $options{ntar_prefix};
				$ntars{$seq_id}->{'ntar_id'}->{$ntar_id} = {start => $novels->{$_}->{start}, end => $novels->{$_}->{end}, strand => "-"};
			}
		}	
    }
}

if( $options{reverse_strand} ) {
	print "\nReversing strand on ntars";
    foreach (sort keys %ntars) {
		
		if($ntars{$_}->{strand} eq "+"){
			$ntars{$_}->{strand} = "-";
			
		} elsif($ntars{$_}->{strand} eq "-") {
			$ntars{$_}->{strand} = "+";
		}
    }
}


my $all_ntars = 0;

# create unfiltered tars gff3 file
open GFF3, ">" . $options{output_stub} . ".ntars.gff3";
print GFF3 "##gff-version   3\n";

foreach my $seq_id (sort {$a cmp $b} keys %ntars) {
	
	my $ntar_ids = $ntars{$seq_id}->{'ntar_id'};
	foreach my $ntar_id (sort {$ntar_ids->{$a}->{start} <=> $ntar_ids->{$b}->{start}} keys %$ntar_ids) {
		
		$all_ntars++;
		
		print GFF3 join("\t", $seq_id, 
						"coverage_based", 
						"transcript", 
						$ntar_ids->{$ntar_id}->{start},
						$ntar_ids->{$ntar_id}->{end},
						".",
						$ntar_ids->{$ntar_id}->{strand},
						".", 
						"ID=$ntar_id;Name=$ntar_id")."\n";
	}
}

close GFF3;


print "\n\nCoverage islands found: $islands\n";
print "Reported TARs found meeting length, distance, coverage criteria: $all_ntars\n";


# if we have an old gff3 available, use it to filter out novel transcripts that overlap or are too close to known ones
if( defined $options{base_gff3} ) {
	filter_transcripts_from_known_transcripts(\%ntars);
}

my $filtered_ntars = 0;

# create filtered tars GFF3 file
open KEPT, ">" . $options{output_stub} . ".kept_ntars.gff3";
print KEPT "##gff-version   3\n";

foreach my $seq_id (sort {$a cmp $b} keys %ntars) {
	
	my $ntar_ids = $ntars{$seq_id}->{'ntar_id'};
	foreach my $ntar_id (sort {$ntar_ids->{$a}->{start} <=> $ntar_ids->{$b}->{start}} keys %$ntar_ids) {

		$filtered_ntars++;
		
		print KEPT join("\t", $seq_id, 
						"coverage_based", 
						"transcript", 
						$ntar_ids->{$ntar_id}->{start},
						$ntar_ids->{$ntar_id}->{end},
						".",
						$ntar_ids->{$ntar_id}->{strand},
						".", 
						"ID=$ntar_id;Name=$ntar_id")."\n";
	}
}

close KEPT;

print STDERR "\nReported TARs after filtering based on overlaps: $filtered_ntars\n\n";


##############################################

# given parameters and an array of perbp coverage, returns a hash with transcript locations
# locations are converted to 1-based coordinates, to match the usual GFF3 specifications
sub find_transcripts {
    my ($min_cov, $max_intron, $min_transcript, $min_intergenic, $coverage) = @_;
	
    # initialize variables
    my $in_intergenic = 1;
    my $intron_len = 0;
    my $gene_len = 0;
	
	# the gff3 will use 1-based coordinates
    my $gstart = 0;
    my $current = 0;
	
    my %genes;

    my $prev_end = 0;

    my @tgene;

    my $estart;
    my $cov_end;

	print "\nCUTOFFS: min_cov=$min_cov, max_intron=$max_intron, min_transcript=$min_transcript, min_intergenic=$min_intergenic\n\n";

    foreach (my $i=0; $i <@{$coverage}; $i++){
		
		if( @{$coverage}[$i] >= $min_cov ) {
			
			if( $in_intergenic == 1 ) {
				# we are starting over
				if( defined $options{verbose} ) {
					print "We are in an integenic region so we are starting over...\n";
				}
				$in_intergenic = 0;
				$gstart=$current;
				$estart = $current;
				@tgene = ();
				
			} elsif( $in_intergenic == 0 && $intron_len > 0 ) {
				# we are entering an exon after leaving an intron
				if( defined $options{verbose} ) {
					print "We are entering an exon after leaving an intron....\n";
				}
				$estart = $current;
			}
			
			$gene_len++;
			$intron_len = 0;
			$cov_end = $current;
			
		} elsif( (@{$coverage}[$i] < $min_cov && $in_intergenic == 0) ) {
			
			if( $intron_len == 0 ) {
				# save the exon to our temp gene
				push @tgene, [$estart, $cov_end];
			}
			
			if( $intron_len < $max_intron ) {
				# intron less than max_intron
				$intron_len++;
				$gene_len++;		
				
			} else {
				# intron was too long, save gene if applicable, start over
				$islands++;
				
				if($gene_len >=$min_transcript && ($gstart>=($prev_end+$min_intergenic) || ($prev_end == 0) ) ){
					$genes{"ntar_".$gene_id} = {'start' => $tgene[0][0]+1, 'end' => $tgene[-1][1] + 1};
					
					for( my $i=0; $i<@tgene; $i++ ) {
						push @{$genes{"ntar_".$gene_id}{exons}}, {'start' => $tgene[$i][0], 'end' => $tgene[$i][1]};
					}
					
					$prev_end = $current;
					$gene_id++;
				}
				
				$in_intergenic = 1;
				$gene_len=0;
			}
		}
		
		# some extra checks here if we hit the end of the chromosome
		if( ($i == @{$coverage} - 1) && ($in_intergenic == 0) ) {
			
			$islands++;
			
			if( $gene_len >= $min_transcript && ($gstart >= ($prev_end + $min_intergenic) || ($prev_end == 0) ) ){
				
				if( $intron_len == 0 ) { 
					push @tgene, [$estart, $cov_end]; #save the exon to our temp gene
				}
				
				$genes{"ntar_".$gene_id} = {'start' => $tgene[0][0]+1, 'end' => $tgene[-1][1]+1};
				
				for(my $i=0; $i<@tgene; $i++){
					push @{$genes{"ntar_".$gene_id}{exons}}, {'start' => $tgene[$i][0], 'end' => $tgene[$i][1]};
				}
				
				$gene_id++;
			}
		}
		
		$current++;
	}
	
    return \%genes;
}

sub filter_transcripts_from_known_transcripts {
	my ($ntars) = @_;
	
	# gff3 file
    my $base_gff3 = NGS::GFF3->new( file => $options{base_gff3} );
	
	# get all features from the gff3 for a given feature type
    my $features = $base_gff3->get_features(feat_type => $options{base_feature_type} );
	
	# iterate through identified tars
	foreach my $seq_id (sort {$a cmp $b} keys %$ntars) {
		
		# get all novel tar ids
		my $ntar_ids = $ntars->{$seq_id}->{'ntar_id'};
		
		foreach my $ntar_id (sort {$ntar_ids->{$a}->{start} <=> $ntar_ids->{$b}->{start}} keys %$ntar_ids) {
			
			# iterate through known genes/transcripts
			foreach my $feat (@{$features}) {
				
				# skip if seq_ids do not match
				next if( $feat->seq_id ne $seq_id );
				
				# skip if not on the same strand
				next if( $options{stranded} && ($feat->{strand} ne $ntar_ids->{$ntar_id}->{strand}) );
				
				# check for overlaps
				my $overlap = $base_gff3->check_for_overlap( {start => ($feat->start - $options{min_intergenic}), 
															  end => ($feat->end + $options{min_intergenic}) },
															 {start => $ntar_ids->{$ntar_id}->{start}, 
															  end => $ntar_ids->{$ntar_id}->{end} });
				
				
				if($overlap->{action} ne "none") {
					
					#if( defined $options{verbose} ) {
					#	print "Removing $ntar_id for $seq_id: $ntar_ids->{$ntar_id}->{start}, $ntar_ids->{$ntar_id}->{end} [" . $feat->start . ", " . $feat->end . "]\n";
					#}

					# delete overalapping tar
					delete($ntar_ids->{$ntar_id});
					last;
				}
			}
		}
	}
}

sub check_parameters {
    my $options = shift;
    
    # make sure required arguments were passed
    my @required = qw( file_type file );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
	
	if( !defined $options{output_stub} ) {
		($options{output_stub}) = ($options{file} =~ /\S+\/(.*)\.\S+/);
		
		if( $options{verbose} ) {
			print "output_stub: $options{output_stub}\n";
		}
	}
	
    $options{stranded} = 0 if(!defined $options{stranded});
    $options{reverse_strand} = 0 if(!defined $options{reverse_strand});
    $options{min_cov} = 10 if(!defined $options{min_cov});
    $options{max_intron} = 50 if(!defined $options{max_intron});
    $options{min_transcript} = 0 if(!defined $options{min_transcript});
    $options{min_intergenic} = 250 if(!defined $options{min_intergenic});
    $options{min_ntar_intergenic} = 1 if(!defined $options{min_ntar_intergenic});

    $options{base_feature_type} = "transcript" if(!defined $options{base_feature_type});   
    $options{ntar_prefix} = "" if(!defined $options{ntar_prefix});

	# if the files are separated by commas, change the separator to 2 underscores (__)
	if( $options{file} =~ /.+\,.+/ ) {
		$options{file} =~ s/\,/__/;
	}
	
	if( $options{stranded} && ( (uc($options{file_type}) eq "WIG") || (uc($options{file_type}) eq "BIGWIG")) && !($options{file} =~ /(.+)__(.+)/) ) {
		die "Provide stranded coverage as a comma-separated pair to --file!\nStranded libraries require two WIG or BigWig files, one for each strand's coverage!";
    }
	
    $options{forward_coverage} = $1;
    $options{reverse_coverage} = $2;

    if( (uc($options{file_type}) eq "WIG") && !(defined $options{sizes_file}) ) {
		die "WIG files require a sizes file to indicate the size of each chromosome!";
    }
}


