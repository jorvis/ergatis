#!/usr/local/bin/perl


use strict;
use Log::Log4perl qw(get_logger);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use English;
use BSML::BsmlBuilder;
use BSML::BsmlReader;
use BSML::BsmlParserTwig;
use File::Basename;



my %options = ();
my $results = GetOptions (\%options, 'btab_dir|b=s', 'bsml_dir|d=s', 'btab_type|t=s', 'output|o=s', 'verbose|v', 'help|h',);


###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $output     = $options{'output'};
my $btab_dir   = $options{'btab_dir'};
$btab_dir =~ s/\/+$//;
my $BSML_dir   = $options{'bsml_dir'};
$BSML_dir =~ s/\/+$//;         #remove terminating '/'s
my $btab_type     = $options{'btab_type'};   # 1 = allvsall , 2 = blast family 
my $verbose    = $options{'verbose'};
#Log::Log4perl->init("log.conf");
#my $logger = get_logger();



if(!$btab_dir or !$output or !$btab_type or !$BSML_dir or exists($options{'help'})) {
    #$logger->fatal("Not all of the required options have been defined.  Exiting...");
    &print_usage();
}

###-------------------------------------------------------###

if(! -d $btab_dir ) {
    print STDERR "$btab_dir does not exist! Aborting...\n";
    #$logger->fatal("The directory \"$btab_dir\" cannot be found.  Exiting...");
    #print STDERR "The directory \"$btab_dir\" cannot be found.  Exiting...\n";
    exit 10;
}

my $doc = new BSML::BsmlBuilder();


my ($gene_asmbl_id, $protID_cdsID);
my @files = <$btab_dir/*.btab>;
if($btab_type == 1) {
    ($gene_asmbl_id, $protID_cdsID) = build_id_lookups_allvsall();
    $doc->makeCurrentDocument();
    parse_allvsall_btabs(\@files);
}else {
    $gene_asmbl_id = build_id_lookups_blast();
    $doc->makeCurrentDocument();
    parse_blast_btabs(\@files);
}

$doc->write($output);
chmod 0666, $output;




sub parse_blast_btabs {

    my $btab_files = shift;

    my $num;
    foreach my $file(@$btab_files) {
	$num++;
	open (BTAB, "$file") or die "Unable to open \"$file\" due to $!";
	print STDERR "opening $file  $num\n" if($verbose);
	my (@btab, $seq);
	while(my $line = <BTAB>) {
	    chomp($line);
	    @btab = split("\t", $line);
	    next if($btab[19] > '1e-15');
	    next if(!$btab[0] or !$btab[5]);
	    my $align = $doc->createAndAddBtabLine(@btab);
	    $seq = $doc->returnBsmlSequenceByIDR($btab[5]);
	    my $match_asmbl_id = $gene_asmbl_id->{$btab[5]};
	    my $query_asmbl_id = $gene_asmbl_id->{$btab[0]};
	    $doc->createAndAddBsmlAttributeN('elem'=> $seq, 'key'=>'ASSEMBLY', 'value'=>$match_asmbl_id);
	}
	close BTAB;
	$seq = $doc->returnBsmlSequenceByIDR($btab[0]);
	if($seq) {
	    my $query_asmbl_id = $gene_asmbl_id->{$btab[0]};
	    $doc->createAndAddBsmlAttributeN('elem'=> $seq, 'key'=>'ASSEMBLY', 'value'=>$query_asmbl_id);
        }
    }
    $doc->createAndAddAnalysis("program" => "blastp 2.0mp-washu", "programversion" => '2.0', 'sourcename' =>$output,
                               "bsml_link_relation" => 'SEQ_PAIR_ALIGNMENTS', 'bsml_link_url' => '#BsmlTables');

}
    
sub parse_allvsall_btabs {

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
	    next if($btab[20] > '1e-15');
	    next if($btab[3] ne 'praze');
	    next if( !($btab[13] > 0) || !($btab[14] > 0) );
	    next if(!$btab[0] or !$btab[5]);
	    $btab[5] =~ s/\|//g;   #get rid of trailing |
	    $query_protein_id = $btab[0];
	    $btab[0] = $protID_cdsID->{$btab[0]};
	    $query_cds_id = $btab[0];	    
            #$query_name = $btab[0];
	    $match_name = $btab[5];
	    splice(@btab, 19, 1);
	    
	    my $align = $doc->createAndAddBtabLine(@btab);
	    $seq = $doc->returnBsmlSequenceByIDR($match_name);
	    my $match_asmbl_id = $gene_asmbl_id->{$match_name};
	    $doc->createAndAddBsmlAttributeN('elem'=> $seq, 'key'=>'ASSEMBLY', 'value'=>"$match_asmbl_id");
	}
	close BTAB;
	$seq = $doc->returnBsmlSequenceByIDR($query_cds_id);
	if($seq) {
	    my $query_asmbl_id = $gene_asmbl_id->{$query_protein_id};
	    $doc->createAndAddBsmlAttributeN('elem'=> $seq, 'key'=>'ASSEMBLY', 'value'=>"$query_asmbl_id");
	}
    }
    $doc->createAndAddAnalysis("program" => "allvsall", "programversion" => '1.0', 'sourcename' =>$output,
                               "bsml_link_relation" => 'SEQ_PAIR_ALIGNMENTS', 'bsml_link_url' => '#BsmlTables');
}



    
sub build_id_lookups_allvsall {

    my @files = <$BSML_dir/*.bsml>;
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

    my @files = <$BSML_dir/*.bsml>;
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
	    #my $rhash = $reader->returnAllIdentifiers();
	    #build_cdsID_protID_mapping($rhash, $cdsID_protID); 	
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
    


sub print_usage {


    print STDERR "SAMPLE USAGE:  btab2bsml.pl -b btab_dir -o output -t 1\n";
    print STDERR "  --btab_dir    = dir containing btab files\n";
    print STDERR "  --output      = file to save output to\n";
    print STDERR "  --btab_type   = type of btab files (1=allvsall, 2=blast)\n";
    print STDERR "  --help = This help message.\n";
    print STDERR "  *optional:    --bsml_dir    = dir containing BSML doc\n";
    exit 1;

}







