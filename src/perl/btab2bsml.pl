#!/usr/local/bin/perl


use lib("shared", "/usr/local/annotation/PNEUMO/clu_dir/BSML/ANNOTATION/bsml/src");
use strict;
use Log::Log4perl qw(get_logger);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use English;
use BsmlBuilder;


my %options = ();
my $results = GetOptions (\%options, 'btab_dir|b=s', 'bsml_dir', 'btab_type|t=s', 'output|o=s', 'help|h',);


###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $output     = $options{'output'};
my $btab_dir   = $options{'btab_dir'};
$btab_dir =~ s/\/+$//;
my $BSML_dir   = $options{'bsml_dir'} || "/usr/local/annotation/PNEUMO/BSML_repository";
$BSML_dir =~ s/\/+$//;         #remove terminating '/'s
my $btab_type     = $options{'btab_type'};   # 1 = allvsall , 2 = blast family 

#Log::Log4perl->init("log.conf");
#my $logger = get_logger();



if(!$btab_dir or !$output or !$btab_type or exists($options{'help'})) {
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


#Make a new bsml object
my $doc = new BsmlBuilder();

my $gene_asmbl_id = get_gene_asmbl_id();  #hash_ref stores each gene and the asmbl_id to which it belongs


my @files = <$btab_dir/*.btab>;
if($btab_type == 1) {
    parse_allvsall_btabs(\@files);
}else {
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
	print STDERR "opening $file  $num\n";
	my (@btab, $seq);
	while(my $line = <BTAB>) {
	    chomp($line);
	    @btab = split("\t", $line);
	    next if($btab[19] > '1e-10');
	    my $align = $doc->createAndAddBtabLine(@btab);
	    $seq = $doc->returnBsmlSequenceByIDR($btab[5]."_aa");
	    my $match_asmbl_id = $gene_asmbl_id->{$btab[5]};
	    my $query_asmbl_id = $gene_asmbl_id->{$btab[0]};
	    $doc->createAndAddBsmlAttributeN('elem'=> $seq, 'key'=>'ASSEMBLY', 'value'=>"PNEUMO_${match_asmbl_id}");
	}
	close BTAB;
	$seq = $doc->returnBsmlSequenceByIDR($btab[0]."_aa");
	if($seq) {
	    my $query_asmbl_id = $gene_asmbl_id->{$btab[0]};
	    $doc->createAndAddBsmlAttributeN('elem'=> $seq, 'key'=>'ASSEMBLY', 'value'=>"PNEUMO_${query_asmbl_id}");
        }
    }

}
    
sub parse_allvsall_btabs {

    my $btab_files = shift;

    my $num;
    foreach my $file (@$btab_files) {
	$num++; 
	open (BTAB, "$file") or die "Unable to open \"$file\" due to $!";
	print STDERR "opening $file  $num\n";
        my (@btab, $seq, $query_name, $match_name);
	while(my $line = <BTAB>) {
	    chomp($line);
	    @btab = split("\t", $line);
	    next if($btab[20] > '1e-10');
	    next if($btab[3] ne 'praze');
	    next if( !($btab[13] > 0) || !($btab[14] > 0) );
	    next if(!$btab[0] or !$btab[5]);
	    $btab[5] =~ s/\|//g;   #get rid of trailing |
	    $query_name = $btab[0];
	    $match_name = $btab[5];
	    splice(@btab, 19, 1);
	    my $align = $doc->createAndAddBtabLine(@btab);
	    $seq = $doc->returnBsmlSequenceByIDR($match_name."_aa");
	    my $match_asmbl_id = $gene_asmbl_id->{$match_name};
	    my $query_asmbl_id = $gene_asmbl_id->{$query_name};
	    $doc->createAndAddBsmlAttributeN('elem'=> $seq, 'key'=>'ASSEMBLY', 'value'=>"PNEUMO_${match_asmbl_id}");
	}
	close BTAB;
	$seq = $doc->returnBsmlSequenceByIDR($query_name."_aa");
	if($seq) {
	    my $query_asmbl_id = $gene_asmbl_id->{$query_name};
	    $doc->createAndAddBsmlAttributeN('elem'=> $seq, 'key'=>'ASSEMBLY', 'value'=>"PNEUMO_${query_asmbl_id}");
	}
    }
}


=hello
sub get_gene_asmbl_id2 {
#This subroutine retrieves all the genes and the asmbl_id to which they belong

    #my $logger = get_logger();
    #$logger->info("Begin retrieving mapping of genes to asmbl_ids");
    
    my ($id, $asmbl_id, $feat_name, $hash_ref);
    my $result = $CGC->fetch_all_gene_genome();
    #$logger->debug("Mapping of gene to asmbl_id contain $result->{'count'} entries");
    for(my $i=0; $i<$result->{'count'}; $i++) {
	$id = $result->{$i}->{'id'};
	$asmbl_id = $result->{$i}->{'asmbl_id'};
        $feat_name = $result->{$i}->{'feat_name'};
        $hash_ref->{$feat_name} = $asmbl_id;
    }

    #$logger->info("Finished retrieveing mapping of genes to asmbl_ids");
    return $hash_ref;

}
=cut


sub get_gene_asmbl_id {

    #my $BSML_dir = "/usr/local/annotation/PNEUMO/BSML_repository";
    #$BSML_dir =~ s/\/+$//;       #remove terminating '/'s
    
    my $gene_asmbl_id={};
    open (IN, "$BSML_dir/gene_asmbl.txt") or die "Unable to open $BSML_dir/gene_asmbl.txt due to $!";
    while(my $line = <IN>) {
	chomp($line);
	my ($feat_name, $asmbl_id) = split("\t", $line);
	$gene_asmbl_id->{$feat_name} = $asmbl_id;
    }

    return $gene_asmbl_id;

}


sub print_usage {


    print STDERR "SAMPLE USAGE:  btab2bsml.pl -b btab_dir -o output\n";
    print STDERR "  --btab_dir    = dir containing btab files\n";
    print STDERR "  --output      = file to save output to\n";
    print STDERR "  --help = This help message.\n";
    print STDERR "  *optional:    --bsml_dir    = dir containing BSML doc\n";
    exit 1;

}







