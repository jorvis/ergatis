#!/usr/local/bin/perl

use strict;
use PEffect::PEffectXML;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BsmlCGCReader;
use BSML::BsmlParserTwig;

my %options = ();
my $results = GetOptions (\%options, 'asmbl_id|a=s', 'output|o=s', 'verbose|v', 'bsml_btab_dir|f=s', 'gene_pos_file=s',
                                     'match_asmbl_id|b=s',  'gene_pos_check', 'bsml_dir|d=s', 'asmbl_file=s', 'help|h' );

###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $ASMBL_IDS = $options{'asmbl_id'};
my $bsml_btab_dir = $options{'bsml_btab_dir'};
$bsml_btab_dir =~ s/\/$//;
my $verbose = $options{'verbose'};
my $check_gene_pos = $options{'gene_pos_check'};
my $gene_pos_file = $options{'gene_pos_file'};
my $match_asmbl_id = $options{'match_asmbl_id'} || 'all';
#my $output    = $options{'output'};
#$output =~ s/\/$//;
my $BSML_dir = $options{'bsml_dir'};
$BSML_dir =~ s/\/$//;
my $asmbl_file      = $options{'asmbl_file'};

if($check_gene_pos and !$gene_pos_file) {
    print STDERR "Must specify the gene position xml file\n";
    &print_usage();
}

if(!$BSML_dir or !$bsml_btab_dir or exists($options{'help'})) {
    &print_usage();
}

if(!$asmbl_file and ! $ASMBL_IDS) {
    print STDERR "Either --asmbl_ids  OR --asmbl_file option is needed\n";
    &print_usage();
}

if($asmbl_file and $ASMBL_IDS) {
    print STDERR " Specify either --asmbl_ids OR --asmbl_file\n"; 
    &print_usage();
}


###-------------------------------------------------------###

my $valid_asmbl_ids = fetch_valid_asmbl_id($gene_pos_file) if($check_gene_pos);
my ($cdsID_protID, $proteinID_seqID) = build_id_lookups($BSML_dir);

my @asm_ids;
if($asmbl_file) {   #asmbl_id will be read from a flat file
    @asm_ids = read_asmbl_file($asmbl_file);
    if(!@asm_ids) {
	print STDERR "No asmbl_ids found in $asmbl_file.  Aborting...\n";
	exit 4;
    }
} else {   #asmbl_id will be read from --asmbl_ids flag
    if($ASMBL_IDS =~ /all/i) {
	my @files = <$BSML_dir/*.bsml>;
	foreach (@files) {
	    my $basename = basename($_);
	    if($basename =~ /(.+)\.bsml/) {
		push(@asm_ids, $1);
	    }
	}
    } else {
	@asm_ids = split(/,/, $ASMBL_IDS);
    }

}
my $pexml = new PEffect::PEffectXML();
my $reader;
foreach my $asmbl_id (@asm_ids) {
    my $bsml_btab_file = "$bsml_btab_dir/$asmbl_id.allvsall.bsml";
    if (!-s $bsml_btab_file) {
	print STDERR "The $bsml_btab_file does not exist!  Skipping...\n";
	next;
    } else {
	print STDERR "Processing $asmbl_id\n";
	$reader = BsmlCGCReader->new(); 
	my $parser = BSML::BsmlParserTwig->new();
	$parser->parse( \$reader, $bsml_btab_file );
	#my $pexml = new PEffect::PEffectXML();
	addMatches($asmbl_id, $match_asmbl_id, $pexml);
	#my($oref);
	#$pexml->outputXML(\$oref);
	#print $oref,"\n";
    }

}
my($oref);
$pexml->outputXML(\$oref);
print $oref,"\n";


sub addMatches {

    my($asmbl_id, $match_asmbl_id, $pexml) = @_;

    my $lref;
    if($match_asmbl_id ne 'all') {
	$lref = $reader->fetch_genome_pairwise_matches( $asmbl_id, $match_asmbl_id );
    } else {
	$lref = $reader->fetch_genome_pairwise_matches( $asmbl_id, 'all');
    }

    if($check_gene_pos) {
	if(!exists($valid_asmbl_ids->{$asmbl_id})) {
	    print STDERR "asmbl_id $asmbl_id NOT in gene position xml.  Skipping...\n";
	    return;
        }
	if($match_asmbl_id ne 'all') {
	    if(!exists($valid_asmbl_ids->{$match_asmbl_id})) {
		print STDERR "match_asmbl_id $match_asmbl_id is NOT in gene position xml. Skipping...\n";
		return;
	    }
        }
    }
	    
    
    foreach my $match (@$lref) {
	my $q_feat_name = $match->{'query_gene_name'};
	$q_feat_name = $cdsID_protID->{$q_feat_name};
        my $m_feat_name = $match->{'match_gene_name'};
	my $per_sim     = $match->{'percent_similarity'};
	my $per_id      = $match->{'percent_identity'};
        my $pvalue      = $match->{'pval'};

	#If check_gene_pos option is enabled, the asmbl_id to which each match_gene_name belongs to MUST exist in the gene position xml
	if($check_gene_pos) {              
	    my $seqID = $proteinID_seqID->{$m_feat_name};  
	    if(!exists($valid_asmbl_ids->{$seqID})) { #skip if asmbl_id NOT in gene position xml
		#print STDERR "$seqID to which $m_feat_name belongs to is NOT in gene position xml.  Skipping...\n";
                next;
	    }
	}
	$pexml->addAlignment($q_feat_name, $m_feat_name, $per_sim, $per_id, $pvalue);
    }

}

sub build_cdsID_protID_mapping {
#This function builds a mapping between cdsID to proteinID. 
#The returned structure is a hash ref, where key is cdsID, value is proteinID

    my $rhash = shift;
    my $cdsID_protID=shift;
    my $proteinID_seqID=shift;

    foreach my $seqID (keys %$rhash) {
	foreach my $geneID (keys %{ $rhash->{$seqID} }) {
	    foreach my $transcriptID (keys %{ $rhash->{$seqID}->{$geneID} }) {
		my $cdsID = $rhash->{$seqID}->{$geneID}->{$transcriptID}->{'cdsId'};
		my $proteinID = $rhash->{$seqID}->{$geneID}->{$transcriptID}->{'proteinId'};
		$cdsID_protID->{$cdsID} = $proteinID;
		$proteinID_seqID->{$proteinID} = $seqID;
	    }
	}
    }

}

sub build_id_lookups {

    my $BSML_dir = shift;
    my $cdsID_protID = {};
    my $proteinID_seqID = {};

    my @files = <$BSML_dir/*.bsml>;
    my $new_parser = new BSML::BsmlParserTwig;

    foreach my $bsml_doc (@files) {
	if (-s $bsml_doc) {
	    print STDERR "parsing $bsml_doc\n" if($verbose);
	    my $reader = BsmlCGCReader->new();
	    $new_parser->parse( \$reader, $bsml_doc );
	    my $rhash = $reader->returnAllIdentifiers();
	    build_cdsID_protID_mapping($rhash, $cdsID_protID, $proteinID_seqID); 
	} else {
	    print STDERR "Empty $bsml_doc...skipping\n" if($verbose);
        }
    }

    return ($cdsID_protID, $proteinID_seqID);

}

sub fetch_valid_asmbl_id {
#This subroutine parses the gene position xml file and grabs a list of valid asmbl_ids
#Returns a hashref with asmbl_ids as keys

    my $gene_position_xml = shift;

    open (IN, "$gene_position_xml") or die "Unable to open $gene_position_xml due to $!";
    my $line;
    my $valid_asmbl;
    while ($line = <IN>) {
	if($line =~ /<fg id=\'db_(.+)\'>/) {
	    my $asmbl_id = $1;
	    $valid_asmbl->{$asmbl_id} = 1;
	}
    }
    close IN;
    return $valid_asmbl;
}

sub read_asmbl_file {
#Reads a flat file containing a list of asmbl_ids to run on
#each asmbl_id is a on a separate line

    my $file = shift;

    my @asmbl_id_list;

    open (IN, "$file")  or die "Unable to read $file due to $!";
    my $line;
    while($line = <IN>) {
	chomp($line);
	next if($line =~ /^\s*$/);
	push(@asmbl_id_list, $line);
    }
    close IN;

    return @asmbl_id_list;

}



sub print_usage {


    print STDERR "SAMPLE USAGE:  gene_match_bsml.pl -a bsp_3839_assembly -b gbs_799_assembly -d bsml_dir -f allvsall.bsml\n";
    print STDERR "  --asmbl_ids  = assembly id\n";
    print STDERR "  --asmbl_file  = name of the file containing a list of asmbl_ids\n";
    print STDERR "  --match_asmbl_ids = assembly id or 'all' \n";
    print STDERR "  --gene_pos_check  = check to see if asmbl is already in gene_pos_file\n";
    print STDERR "  --gene_pos_file  = gene position xml file \n";
    print STDERR "  --bsml_dir = directory containing the BSML documents\n";
    print STDERR "  --bsml_btab_dir(-f) = dir containing bsml encoding btabs\n";
    print STDERR "  --help = This help message.\n";
    exit 1;

}
