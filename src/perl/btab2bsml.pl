#!/usr/local/bin/perl

=head1  NAME 

btab2bsml.pl  - convert info stored in btab files into BSML documents

=head1 SYNOPSIS

USAGE:  btab2bsml.pl -b btab_dir -o blastp.bsml -d bsml_dir -t 1

=head1 OPTIONS

=over 4

=item *

B<--bsml_dir,-d> [REQUIRED]  Dir containing BSML documents (repository)

=item *

B<--output,-o> [REQUIRED] output BSML file containing bit_score information

=item * 

B<--btab_dir,-b> [REQUIRED] Dir containing btab files

=item *

B<--btab_type,-t> [REQUIRED] Type of btab files: 1=ber, 2=blastp

=item *

B<--help,-h> This help message

=back

=head1   DESCRIPTION

btab2bsml.pl is designed to convert information in btab files into BSML documents
There 2 types of btab files.  One is generated from allvsall.  The other is
generated from blast family of programs.  The user can specify which type of
btab files to convert by --btab_type(t) flag.  An argument of 1 means allvsall, 2 
means blastp.  

Samples:

1. convert blastp btab files in /usr/btab to blastp.bsml 

   btab2bsml.pl -d bsml_dir -b /usr/btab -o blastp.bsml -t 2 

NOTE:  

Calling the script name with NO flags/options or --help will display the syntax requirement.


=cut


use strict;
use Log::Log4perl qw(get_logger);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use English;
use BSML::BsmlBuilder;
use BSML::BsmlReader;
use BSML::BsmlParserTwig;
use File::Basename;
use File::Path;
use Pod::Usage;


my %options = ();
my $results = GetOptions (\%options, 'btab_dir|b=s', 'bsml_dir|d=s', 'btab_type|t=s', 
                                     'output|o=s', 'max_hsp_count|m=s', 'cutoff|c=s', 'verbose|v', 'help|h', 'man') || pod2usage();


###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $output     = $options{'output'};
my $output_dir = dirname($output);
my $btab_dir   = $options{'btab_dir'};
$btab_dir =~ s/\/+$//;
my $BSML_dir   = $options{'bsml_dir'};
$BSML_dir =~ s/\/+$//;        
my $btab_type     = $options{'btab_type'};   # 1 = ber , 2 = blast family 
my $verbose    = $options{'verbose'};
my $cutoff = $options{'cutoff'};
my $hspcutoff = $options{'max_hsp_count'};

&cmd_check();
###-------------------------------------------------------###

my $doc = new BSML::BsmlBuilder();


my ($gene_asmbl_id, $protID_cdsID);
my @files;
opendir(DIR, $btab_dir) or die "Unable to access $btab_dir due to $!";
while( my $file = readdir(DIR)) {
    next if ($file =~ /^\.{1,2}$/);  #skip  "." ,  ".."
    if($file =~ /(.+)\.btab$/) {
	push (@files, "$btab_dir/$file");
    }
}

#my $bsmlrepo = new BSML::BsmlRepository('REPOSITORY_ROOT'=>$BSML_dir);

#my ($files) = $bsmlrepo->list_bsml_files();

if($btab_type == 1) {
    ($gene_asmbl_id, $protID_cdsID) = build_id_lookups_ber($BSML_dir);
    $doc->makeCurrentDocument();
    parse_ber_btabs(\@files);
}elsif($btab_type == 2) {
    $gene_asmbl_id = build_id_lookups_blast($BSML_dir);
    $doc->makeCurrentDocument();
    parse_blast_btabs(\@files);
} else {
    print STDERR "Bad btab_type.  Aborting...\n";
    exit 5;
}

$doc->write($output);

sub parse_blast_btabs {

    my $btab_files = shift;

    my $num;
    foreach my $file(@$btab_files) {
	$num++;
	open (BTAB, "$file") or die "Unable to open \"$file\" due to $!";
	print STDERR "opening $file  $num using cutoff $cutoff\n" if($verbose);
	my $hsplookup;
	while(my $line = <BTAB>) {
	    chomp($line);
	    my @btab = split("\t", $line);
	    if(($btab[19] < $cutoff) && ($btab[0] ne "") && ($btab[5] ne "")){
		if(!(exists $hsplookup->{$btab[0]}->{$btab[5]})){
		    $hsplookup->{$btab[0]}->{$btab[5]} = [];
		}
		push @{$hsplookup->{$btab[0]}->{$btab[5]}},{'pvalue'=>$btab[19],
							    'line'=>$line};
	    }
	}
	close BTAB;
	foreach my $query (keys %$hsplookup){
	    foreach my $subject (keys %{$hsplookup->{$query}}){
		my @hsps = sort {$a->{'pvalue'} <=> $b->{'pvalue'}} @{$hsplookup->{$query}->{$subject}};
		my $maxhsp;
		if($hspcutoff ne ""){
		    $maxhsp = ($hspcutoff<scalar(@hsps)) ? $hspcutoff : scalar(@hsps);
		}
		else{
		    $maxhsp = scalar(@hsps);
		}
		my $queryid;
		for(my $i=0;$i<$maxhsp;$i++){
		    my $line = $hsps[$i]->{'line'};
		    my @btab = split("\t", $line);
		    print STDERR "Storing HSP $btab[0] $btab[5] $btab[19] with max hsp_count: $maxhsp\n" if($verbose);;
		    $queryid = $btab[0] if($btab[0] && (!$queryid));
		    my $align = $doc->createAndAddBtabLine(@btab);
		    my $seq = $doc->returnBsmlSequenceByIDR($btab[5]);
		    my $match_asmbl_id = $gene_asmbl_id->{$btab[5]};
		    $doc->createAndAddBsmlAttributeN('elem'=> $seq, 'key'=>'ASSEMBLY', 'value'=>$match_asmbl_id);
		}
		if($queryid){
		    my $seq = $doc->returnBsmlSequenceByIDR($queryid);
		    if($seq) {
			my $query_asmbl_id = $gene_asmbl_id->{$queryid};
			$doc->createAndAddBsmlAttributeN('elem'=> $seq, 'key'=>'ASSEMBLY', 'value'=>$query_asmbl_id);
		    }
		}
	    }
	}
    }
    $doc->createAndAddAnalysis("program" => "washu blastp", "programversion" => '2.0mp', 'sourcename' =>$output,
                               "bsml_link_relation" => 'SEQ_PAIR_ALIGNMENTS', 'bsml_link_url' => '#BsmlTables');

}
    
sub parse_ber_btabs {

    my $btab_files = shift;
    my $num;
    foreach my $file (@$btab_files) {
	$num++; 
	open (BTAB, "$file") or die "Unable to open \"$file\" due to $!";
	print STDERR "opening $file  $num\n" if($verbose);
        my (@btab, $seq, $query_name, $match_name, $query_cds_id, $query_protein_id);
	while(my $line = <BTAB>) {
	    chomp($line);
	    @btab = split("\t", $line);
	    next if($btab[20] > $cutoff );
	    next if($btab[3] ne 'praze');
	    next if( !($btab[13] > 0) || !($btab[14] > 0) );
	    next if(!$btab[0] or !$btab[5]);
	    $btab[5] =~ s/\|//g;   #get rid of trailing |
	    $query_protein_id = $btab[0];
	    $btab[0] = $protID_cdsID->{$btab[0]};
	    $query_cds_id = $btab[0];	    
	    $match_name = $btab[5];
	    splice(@btab, 19, 1);
	    
	    my $align = $doc->createAndAddBtabLine(@btab);
	    $seq = $doc->returnBsmlSequenceByIDR($match_name);
	    my $match_asmbl_id = $gene_asmbl_id->{$match_name};
	    $doc->createAndAddBsmlAttributeN('elem'=> $seq, 'key'=>'ASSEMBLY', 'value'=>"$match_asmbl_id");
	}
	close BTAB;
	next if(!defined($query_cds_id));
	$seq = $doc->returnBsmlSequenceByIDR($query_cds_id);
	if($seq) {
	    my $query_asmbl_id = $gene_asmbl_id->{$query_protein_id};
	    $doc->createAndAddBsmlAttributeN('elem'=> $seq, 'key'=>'ASSEMBLY', 'value'=>"$query_asmbl_id");
	}
    }
    $doc->createAndAddAnalysis("program" => "ber", "programversion" => '-', 'sourcename' =>$output,
                               "bsml_link_relation" => 'SEQ_PAIR_ALIGNMENTS', 'bsml_link_url' => '#BsmlTables');
}



    
sub build_id_lookups_ber {

    my $BSML_dir = shift;

    my @files;
    opendir(DIR, $BSML_dir) or die "Unable to access $BSML_dir due to $!";
    while( my $file = readdir(DIR)) {
	next if ($file =~ /^\.{1,2}$/);  #skip  "." ,  ".."
	if($file =~ /(.+)\.bsml$/) {
	    push(@files, "$BSML_dir/$file");
	}
    }
    my $parser = new BSML::BsmlParserTwig;

    my $protID_cdsID = {};
    my $gene_asmbl_id ={};
    my $reader;
    foreach my $bsml_doc (@files) {
	if (-s $bsml_doc) {
	    print STDERR "parsing $bsml_doc\n" if($verbose);
	    $reader = BSML::BsmlReader->new();
	    $parser->parse( \$reader, $bsml_doc );
	    my $hash_ref = $reader->get_all_protein_assemblyId();
	    while(my ($protein_id, $asmbl_id) = each %$hash_ref) {
		$gene_asmbl_id->{$protein_id} = $asmbl_id;
	    }
	    my $rhash = $reader->returnAllIdentifiers();
	    build_protID_cdsID_mapping($rhash, $protID_cdsID); 	
	} else {  print STDERR "Empty $bsml_doc...skipping\n" if($verbose); }
    }

    return ($gene_asmbl_id, $protID_cdsID);

}


sub build_id_lookups_blast {

    my $BSML_dir = shift;

    my @files;
    opendir(DIR, $BSML_dir) or die "Unable to access $BSML_dir due to $!";
    while( my $file = readdir(DIR)) {
	next if ($file =~ /^\.{1,2}$/);  #skip  "." ,  ".."
	if($file =~ /(.+)\.bsml$/) {
	    push(@files, "$BSML_dir/$file");
	}
    }
    my $parser = new BSML::BsmlParserTwig;

    my $gene_asmbl_id ={};
    my $reader;
    foreach my $bsml_doc (@files) {
	if (-s $bsml_doc) {
	    print STDERR "parsing $bsml_doc\n" if($verbose);
	    $reader = BSML::BsmlReader->new();
	    $parser->parse( \$reader, $bsml_doc );
	    my $hash_ref = $reader->get_all_protein_assemblyId();
	    while(my ($protein_id, $asmbl_id) = each %$hash_ref) {
		$gene_asmbl_id->{$protein_id} = $asmbl_id;
	    }
	} else {  print STDERR "Empty $bsml_doc...skipping\n" if($verbose); }
    }

    return $gene_asmbl_id;

}

sub build_protID_cdsID_mapping {
#This function builds a mapping between cdsID to proteinID. 
#The returned structure is a hash ref, where key is protID, value is cdsID

    my $rhash = shift;
    my $protID_cdsID=shift;

    foreach my $seqID (keys %$rhash) {
	foreach my $geneID (keys %{ $rhash->{$seqID} }) {
	    foreach my $transcriptID (keys %{ $rhash->{$seqID}->{$geneID} }) {
		my $cdsID = $rhash->{$seqID}->{$geneID}->{$transcriptID}->{'cdsId'};
		my $proteinID = $rhash->{$seqID}->{$geneID}->{$transcriptID}->{'proteinId'};
		$protID_cdsID->{$proteinID} = $cdsID;
	    }
	}
    }

}
    

sub cmd_check {
#quality check

    if( exists($options{'man'})) {
	pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT});
    }   
    
    if( exists($options{'help'})) {
	pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT});
    }
    
    if(!$BSML_dir or !$output or !$btab_type or !$btab_dir) {
	pod2usage({-exitval => 2,  -message => "$0: All the required options are not specified", -verbose => 1, -output => \*STDERR});    
    }

    if(!-d $BSML_dir) {
	print STDERR "bsml repository directory \"$BSML_dir\" cannot be found.  Aborting...\n";
	exit 5;
    }

    if($btab_type != '1' and $btab_type != '2') {
	pod2usage({-exitval => 2,  -message => "$0: btab_type can only be 1 or 2", -verbose => 1, -output => \*STDERR}); 
    }

    if(! -d $btab_dir) {
	print STDERR "btab directory \"btab_dir\" cannot be found.  Aborting...\n";
	exit 5;
    }

    #check for presence of output directory
    if(! -d $output_dir) {
	mkpath($output_dir) or die "Unable to create $output_dir.  Aborting...\n";
	#chmod 0777, $output_dir;
    }

}










