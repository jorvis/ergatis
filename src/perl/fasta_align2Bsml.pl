#! /local/perl/bin/perl

=head1  NAME 

fasta_align2bsml.pl  - convert FA format alignment files into BSML documents

=head1 SYNOPSIS

USAGE:  fasta_align2bsml.pl -a ali_dir -o msf.bsml

=head1 OPTIONS

=over 4

=item *

B<--ali_dir,-a>   [REQUIRED] Dir containing FA format alignment files that ends in *.fa

=item *

B<--output,-o> [REQUIRED] output BSML file

=item *

B<--suffix,-s> suffix of the multiple sequence alignment file. Default(fa)

=item *

B<--mol_type,-m> molecule_type of alignment. Default(protein)

=item *

B<--help,-h> This help message

=back

=head1   DESCRIPTION

fasta_align2bsml.pl is designed to convert multiple sequence alignments in FA format
into BSML documents.  More specifically this script is equippd to parse the multiple 
sequence alignments of domains.  This requires the capture of the start and end 
coordinates of each "subsequence" in the alignment.  As an option, the user can
specify the suffix of the alignment files that the script should parse in the 
directory specified by --ali_dir.  The default suffix is "fa" i.e. domain_1.fa.

Samples:

1. Convert all FA formatted files that end in "aln" in ali_dir
   fasta_align2bsml.pl -a ali_dir -o alignment.bsml -s aln


NOTE:  

Calling the script name with NO flags/options or --help will display the syntax requirement.


=cut



#use lib '/export/CVS/bsml/src';
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BSML::BsmlBuilder;
use Data::Dumper;
use File::Basename;
use Pod::Usage;


my %options = ();
my $results = GetOptions (\%options, 'ali_dir|a=s', 'suffix|s=s', 'output|o=s', 'mol_type', 
                                     'program|p=s', 'verbose|v', 'help|h', 'man') || pod2usage();

###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $output    = $options{'output'};
my $ali_dir   = $options{'ali_dir'};
$ali_dir =~ s/\/+$// if($ali_dir);
my $molecule_type = $options{'mol_type'} || 'protein';
my $ali_suffix = $options{'suffix'} || "fa";  #default suffix pattern of alignment file
my $program   = $options{'program'} || "clustal";

&cmd_check();

###-------------------------------------------------------###

my $builder = new BSML::BsmlBuilder;

# add an analysis object to the document. 
my $filename = basename($output);
$builder->createAndAddAnalysis( 'sourcename' => $filename,
				'programversion' => '1.0',
				'program' => $program,
				'bsml_link_url' => 'BsmlTables',
				'bsml_link_relation' => 'MULTIPLE_ALIGNMENTS' );

my $file_found = 0;
opendir(DIR, $ali_dir) or die "Unable to access $ali_dir due to $!";
while( my $file = readdir(DIR)) {
    next if ($file =~ /^\.{1,2}$/);  #skip  "." ,  ".."
    if($file =~ /^(.+?)(\.msf\.dtab)?\.$ali_suffix$/) {
	$file_found++;
	my $accession = undef;
	my $domain_id = "domain_";
	if($file =~ /^(PF\d{5})\.$ali_suffix$/ || $file =~ /^(TIGR\d{5})\.$ali_suffix$/) { 
	    $accession = $1;
	    $domain_id .= $accession;
	}
	#my $fam_name = $1;
	my $mult_alignment = process_fasta_alignment_file("$ali_dir/$file");
	next if(keys %$mult_alignment < 1);   #skip empty alignment files
	my $table = $builder->createAndAddMultipleAlignmentTable('molecule-type' => $molecule_type);
	$table->addattr( 'id', $domain_id );

	my $summary = $builder->createAndAddAlignmentSummary( 'multipleAlignmentTable' => $table,
							      'seq-type' => $molecule_type,
							      'seq-format' => 'FastaAlign' 
                                                            );

	my $aln = $builder->createAndAddSequenceAlignment( 'multipleAlignmentTable' => $table );
	my $seqnum=0;
	my $sequences_tag;
	foreach my $seq (keys %$mult_alignment) {
	    $seqnum++;
	    my ($pro_acc, $coords) = split(/\//, $seq);
	    $pro_acc =~ s/\|\|$//;                         #remove terminal '||' if exists
	    my $alignment = $mult_alignment->{$seq}->{'alignment'};
	    my $start = $mult_alignment->{$seq}->{'start'};
	    my $length = $mult_alignment->{$seq}->{'domain_length'};
	    $builder->createAndAddAlignedSequence( 'alignmentSummary' => $summary,
						   'seqnum' => $seqnum,
						   'length' => $length,
						   'name'   => "$pro_acc:$seqnum",
                                                   'start'  => $start
                                                 );
	    
	    $builder->createAndAddSequenceData( 'sequenceAlignment' => $aln,
						'seq-name' => "$pro_acc:$seqnum",
						'seq-data' => $alignment
                                              );



	    $sequences_tag .= "$seqnum:";
	}
	$aln->addattr( 'sequences', $sequences_tag );
    }

}

if($file_found) {
    $builder->write( $output );
    chmod 0666, $output;
} else {
    print STDERR "No files that ends in \"$ali_suffix\" are found in $ali_dir\n";
}



sub process_fasta_alignment_file {

    my $file = shift;

    my $FA_alignments = {};
    open(FA, "$file") or die "Unable to open $file due to $!";
    my $line = <FA>;
    while(defined($line)) {
	my $sequence;
	if($line =~ /^>(.+)/) {
	    my $header = $1;
	    my ($pro_acc, $coords) = split(/\//, $header);
	    $pro_acc =~ s/\|\|$//;                         #remove terminal '||' if exists
	    my ($rel5, $rel3) = split(/-/, $coords);
	    if($rel3 < $rel5) {
		print STDERR "BAD coordinates for $pro_acc in $file\n";
	    }
	    my $domain_len = $rel3 - $rel5 +1;
	    while(defined($line=<FA>) and $line !~ /^>/ ) {
		next if($line =~/^\s+$/);                   #skip blank lines
		#chomp($line);
		$sequence .= $line;
	    }
	    $FA_alignments->{$header}->{'alignment'} = $sequence;
	    $FA_alignments->{$header}->{'start'} = $rel5;
	    $FA_alignments->{$header}->{'end'} = $rel3;
	    $FA_alignments->{$header}->{'domain_length'} = $domain_len;
	} else { $line = <FA> };
    }
    close FA;

    return $FA_alignments;

}


sub cmd_check {
#quality check

    if( exists($options{'man'})) {
	pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT});
    }   

    if( exists($options{'help'})) {
	pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT});
    }

    if(!$output or !$ali_dir) {
	pod2usage({-exitval => 2,  -message => "$0: All the required options are not specified", -verbose => 1, -output => \*STDERR});
    }
    
    if(! -d $ali_dir) {
	print STDERR "$ali_dir directory NOT found.  Aborting...\n";
	exit 5;
    }

}	    

