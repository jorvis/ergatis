#!/usr/local/bin/perl


use strict;
#use lib("/usr/local/annotation/PNEUMO/clu_dir/prok_CGC"); 
use lib("../..", "/usr/local/annotation/PNEUMO/clu_dir/BSML/ANNOTATION/bsml/src");
use PEffect::PEffectXML;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BsmlReader;
use BsmlParserTwig;

my %options = ();
my $results = GetOptions (\%options, 'asmbl_id|a=s', 'output|o=s', 'match_asmbl_id|b=s',  'bsml_dir|d=s', 'help|h' );

###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my @asmbls;
my $asmbl_id = $options{'asmbl_id'};
#match_asmbl_id can be specified by comma separated list or invoking -b flag multiple times
#my @match_asmbls;
#if(exists($options{'match_asmbl_ids'})) {
#    @match_asmbls= @{$options{'match_asmbl_ids'}};
#    @match_asmbls = split(/,/,join(',',@match_asmbls));
#}
my $match_asmbl_id = $options{'match_asmbl_id'};
#my $output    = $options{'output'};
#$output =~ s/\/$//;
my $BSML_dir = $options{'bsml_dir'};
$BSML_dir =~ s/\/$//;

if(!$asmbl_id or !$match_asmbl_id or !$BSML_dir or exists($options{'help'})) {
    &print_usage();
}

###-------------------------------------------------------###
my $bsml_file = "$BSML_dir/asmbl_${asmbl_id}.bsml_allvsall";
my $reader;
if (!-s $bsml_file) {
    print STDERR "The $bsml_file does not exist!  Aborting...\n";
    exit 5;
} else {
    $reader = new BsmlReader;
    my $parser = new BsmlParserTwig;
    $parser->parse( \$reader, $bsml_file );
    my $pexml = new PEffect::PEffectXML();
    addMatches($asmbl_id, $match_asmbl_id, $pexml);
    my($oref);
    $pexml->outputXML(\$oref);
    print $oref,"\n";
}



sub addMatches {

    my($asmbl_id, $match_asmbl_id, $pexml) = @_;

    my $lref;
    if($match_asmbl_id ne 'all') {
	$lref = $reader->fetch_genome_pairwise_matches( "PNEUMO_${asmbl_id}", "PNEUMO_${match_asmbl_id}" );
    } else {
	$lref = $reader->fetch_genome_pairwise_matches( "PNEUMO_${asmbl_id}", 'all');
    }
    
    foreach my $match (@$lref) {
	my $q_feat_name = $match->{'query_gene_name'};
	$q_feat_name =~ s/\_aa$//;
        my $m_feat_name = $match->{'match_gene_name'};
        $m_feat_name =~ s/\_aa$//; 
	my $per_sim     = $match->{'percent_similarity'};
	my $per_id      = $match->{'percent_identity'};
        my $pvalue      = $match->{'pval'};
	$pexml->addAlignment($q_feat_name, $m_feat_name, $per_sim, $per_id, $pvalue);
    }

}

sub print_usage {


    print STDERR "SAMPLE USAGE:  gene_match_bsml.pl -a 1 -b 3 -d bsml_dir\n";
    print STDERR "  --asmbl_ids  = assembly id\n";
    print STDERR "  --match_asmbl_ids = assembly id or 'all' \n";
    print STDERR "  --bsml_dir = directory containing the BSML documents\n";
    print STDERR "  --help = This help message.\n";
    exit 1;

}
