#!/usr/local/bin/perl


use strict;
use PEffect::PEffectXML;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BSML::BsmlReader;
use BSML::BsmlParserTwig;


my %options = ();
my $results = GetOptions (\%options, 'asmbl_ids|a=s@', 'output|o=s', 'bsml_dir|b=s', 'help|h' );
###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my @asmbls;
#asmbl_id can be specified by comma separated list or invoking -a flag multiple times
if(exists($options{'asmbl_ids'})) {
    @asmbls= @{$options{'asmbl_ids'}};
    @asmbls = split(/,/,join(',',@asmbls));
}
my $BSML_dir = $options{'bsml_dir'};
$BSML_dir =~ s/\/$//;
my $output = $options{'output'};

if(!$BSML_dir or exists($options{'help'})) {
    &print_usage();
}

###-------------------------------------------------------###

my $parser = new BSML::BsmlParserTwig;
my $pexml = new PEffect::PEffectXML();
 
#my @asmbls;
if(!@asmbls) {
    my @files = <$BSML_dir/*.bsml>;
    foreach (@files) {
	if(/asmbl_(\d+)\.bsml/) {
	    push(@asmbls, $1);
        }
    } 
}

foreach my $asmbl_id (@asmbls){
    if($asmbl_id !~ /\d+/) {
	print STDERR "bad asmbl_id \"$asmbl_id\"\n";
        next;
    }
    addGenes($asmbl_id, "db_${asmbl_id}", $pexml);
}


#print out the entire xml to STDOUT
my($oref);
$pexml->outputXML(\$oref);
print $oref,"\n";




#----------------------------------------------------------------

sub addGenes {
#This subroutine grabs all the genes' length, orientation and starting coordinate
#for a given asmbl_id and uses PEffectXML module to make xml

    my $asmbl_id = shift;
    my $name = shift;
    my $pexml = shift;

    my $bsml_file = "$BSML_dir/asmbl_${asmbl_id}.bsml";
    if (-s $bsml_file) {
	my $reader = BsmlReader->new();
	$parser->parse( \$reader, $bsml_file );
	my $order = 0;
	my $sorted_genes = get_sorted_gene_position($asmbl_id, $reader);
	foreach my $gene (@$sorted_genes) {
	    my $attrref={};
	    my $feat_name = $gene->{'feat_name'};
	    $attrref->{'length'} = $gene->{'length'}; 
	    $attrref->{'orient'} = $gene->{'orient'}; 
	    $attrref->{'coord'}  = $gene->{'coord'}; 
	    $pexml->addFeature($feat_name, $attrref, $name, $order);
	    $order++;
	}
    } else {
	print STDERR "$bsml_file does not exist!!! skipping...\n";
    }
	
}

sub get_sorted_gene_position {

    my $asmbl_id = shift;
    my $reader   = shift;

    my $gene_pos = $reader->fetch_gene_positions("PNEUMO_${asmbl_id}");
    my $array_ref=[];
    foreach (@$gene_pos) {
	foreach my $gene (keys %$_) {
	    my $end5 = $_->{$gene}->{'startpos'};
	    my $end3 = $_->{$gene}->{'endpos'};
	    my $complement = $_->{$gene}->{'complement'} == 0 ? '-' : '+';
	    my $length = abs($end3 - $end5);
	    my $coord  = $end5 > $end3 ? $end5 : $end3;
	    push( @$array_ref, { 'feat_name' => $gene, 'length' => $length, 'orient' => $complement, 'coord' => $coord } ); 
	}
    }
    my @sorted_ref = sort { $a->{'coord'} <=> $b->{'coord'} || $a->{'feat_name'} cmp $b->{'feat_name'} } @$array_ref;
    
    return (\@sorted_ref);

}

sub print_usage {


    print STDERR "SAMPLE USAGE:  gene_pos_bsml.pl -b bsml_dir -a 1,2,3,4\n";
    print STDERR "  --bsml_dir    = dir containing BSML doc\n";
    print STDERR "  --asmbl_ids  = assembly ids (comma separated) \n";
    print STDERR "  --help = This help message.\n";
    exit 1;

}






