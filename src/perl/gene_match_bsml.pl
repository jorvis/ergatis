#!/usr/local/bin/perl


use strict;
use PEffect::PEffectXML;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BSML::BsmlReader;
use BSML::BsmlParserTwig;

my %options = ();
my $results = GetOptions (\%options, 'asmbl_id|a=s', 'output|o=s', 'verbose|v', 'bsml_btab|f=s', 'match_asmbl_id|b=s',  'bsml_dir|d=s', 'help|h' );

###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my @asmbls;
my $asmbl_id = $options{'asmbl_id'};
my $bsml_btab_file = $options{'bsml_btab'};
my $verbose = $options{'verbose'};
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

if(!$asmbl_id or !$match_asmbl_id or !$BSML_dir or !$bsml_btab_file or exists($options{'help'})) {
    &print_usage();
}

###-------------------------------------------------------###


my $cdsID_protID = build_id_lookups($BSML_dir);
my $reader;
if (!-s $bsml_btab_file) {
    print STDERR "The $bsml_btab_file does not exist!  Aborting...\n";
    exit 5;
} else {
    $reader = BSML::BsmlReader->new();
    my $parser = BSML::BsmlParserTwig->new();
    $parser->parse( \$reader, $bsml_btab_file );
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
	$lref = $reader->fetch_genome_pairwise_matches( $asmbl_id, $match_asmbl_id );
    } else {
	$lref = $reader->fetch_genome_pairwise_matches( $asmbl_id, 'all');
    }
    
    foreach my $match (@$lref) {
	my $q_feat_name = $match->{'query_gene_name'};
	$q_feat_name = $cdsID_protID->{$q_feat_name};
        my $m_feat_name = $match->{'match_gene_name'};
	my $per_sim     = $match->{'percent_similarity'};
	my $per_id      = $match->{'percent_identity'};
        my $pvalue      = $match->{'pval'};
	$pexml->addAlignment($q_feat_name, $m_feat_name, $per_sim, $per_id, $pvalue);
    }

}

sub build_cdsID_protID_mapping {
#This function builds a mapping between cdsID to proteinID. 
#The returned structure is a hash ref, where key is cdsID, value is proteinID

    my $rhash = shift;
    my $cdsID_protID=shift;

    foreach my $seqID (keys %$rhash) {
	foreach my $geneID (keys %{ $rhash->{$seqID} }) {
	    foreach my $transcriptID (keys %{ $rhash->{$seqID}->{$geneID} }) {
		my $cdsID = $rhash->{$seqID}->{$geneID}->{$transcriptID}->{'cdsId'};
		my $proteinID = $rhash->{$seqID}->{$geneID}->{$transcriptID}->{'proteinId'};
		$cdsID_protID->{$cdsID} = $proteinID;
	    }
	}
    }

}


sub build_id_lookups {

    my $BSML_dir = shift;
    my $cdsID_protID = {};

    my @files = <$BSML_dir/*.bsml>;
    my $new_parser = new BSML::BsmlParserTwig;

    foreach my $bsml_doc (@files) {
	if (-s $bsml_doc) {
	    print STDERR "parsing $bsml_doc\n" if($verbose);
	    my $reader = BSML::BsmlReader->new();
	    $new_parser->parse( \$reader, $bsml_doc );
	    my $rhash = $reader->returnAllIdentifiers();
	    build_cdsID_protID_mapping($rhash, $cdsID_protID); 
	} else {
	    print STDERR "Empty $bsml_doc...skipping\n" if($verbose);
        }
    }

    return $cdsID_protID;

}

sub print_usage {


    print STDERR "SAMPLE USAGE:  gene_match_bsml.pl -a bsp_3839_assembly -b gbs_799_assembly -d bsml_dir -f allvsall.bsml\n";
    print STDERR "  --asmbl_ids  = assembly id\n";
    print STDERR "  --match_asmbl_ids = assembly id or 'all' \n";
    print STDERR "  --bsml_dir = directory containing the BSML documents\n";
    print STDERR "  --bsml_btab = bsmldoc encoding the btab info\n";
    print STDERR "  --help = This help message.\n";
    exit 1;

}
