#!/usr/local/bin/perl


use strict;
#use Log::Log4perl qw(get_logger);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use English;
use BSML::BsmlBuilder;
use BSML::BsmlParserTwig;
use File::Basename;



my %options = ();
my $results = GetOptions (\%options, 'mummer_coords|m=s', 'bsml_dir|b=s', 'mummer_type|t=s', 'output|o=s', 'verbose|v', 'help|h',);


###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $mummer_file     = $options{'mummer_coords'};
my $output          = $options{'output'};
my $mummer_type     = $options{'mummer_type'};   # 1 = nucmer , 2 = promer
my $BSML_dir        = $options{'bsml_dir'};
$BSML_dir =~ s/\/+$//;         #remove terminating '/'s
my $verbose    = $options{'verbose'};
#Log::Log4perl->init("log.conf");
#my $logger = get_logger();

if(!$mummer_file or !$mummer_type or !$output  or exists($options{'help'})) {
    #$logger->fatal("Not all of the required options have been defined.  Exiting...");
    &print_usage();
}

###-------------------------------------------------------###

my $doc = BSML::BsmlBuilder->new();
if($mummer_type == 1) {
    parse_nucmer_coords($mummer_file);
}elsif($mummer_type == 2) {
    parse_promer_coords($mummer_file);
}else {
    print STDERR "Bad mummer_type.  Aborting...\n";
}

$doc->write($output);
chmod 0666, $output;


sub parse_promer_coords {

    my $coords_file = shift;

    my @promer;
    open (PROMER, $coords_file) or die "Unable to open \"$coords_file\" due to $!";
    while (my $line = <PROMER>) {
	chomp($line);
	@promer = split("\t", $line);	
	#In Promer ref seq refers to query
	my $ref_name    = $promer[15];
        my $qry_name    = $promer[16];
	my $ref_start   = $promer[0];
	my $ref_end     = $promer[1];
	my $ref_length  = $promer[4];
	my $qry_start   = $promer[2];
        my $qry_end     = $promer[3];
	my $qry_length  = $promer[5];
	my $percent_id  = $promer[6];
        my $percent_sim = $promer[7];
	my $percent_stop = $promer[8];
	my $qry_asmbl_length = $promer[10];
        my $ref_asmbl_length = $promer[9];
	my $frame_ref = $promer[13];
        my $frame_qry = $promer[14];

	my $aln = $doc->createAndAddSequencePairAlignment( 'refseq'        => $ref_name,
							   'compseq'       => $qry_name,
                                                           'complength'    => $qry_asmbl_length,
                                                           'reflength'     => $ref_asmbl_length
							 ); 
	my $s = $doc->createAndAddSequencePairRun( 'alignment_pair'   => $aln,
						   'refpos'           => $ref_start,
						   'runlength'        => $ref_length,
						   'comppos'          => $qry_start,
						   'comprunlength'    => $qry_length,
                                                   'percent_identity' => $percent_id,
						   'percent_similarity' => $percent_sim
                                                  );
	#additional attributes
	$s->addBsmlAttr( 'percent_stop_codon',  $percent_stop);
	$s->addBsmlAttr( 'refframe', $frame_ref );
	$s->addBsmlAttr( 'compframe', $frame_qry );



    }

}


sub parse_nucmer_coords {

    my $coords_file = shift;

    my @mummer;
    open (NUCMER, $coords_file) or die "Unable to open \"$coords_file\" due to $!";
    while (my $line = <NUCMER>) {
	chomp($line);
	@mummer = split("\t", $line);
	#In Mummer ref seq refers to query
	my $ref_name = $mummer[11];
        my $qry_name = $mummer[12];
	my $ref_start = $mummer[0];
	my $ref_end   = $mummer[1];
	my $ref_length = $mummer[4];
	my $qry_start = $mummer[2];
        my $qry_end   = $mummer[3];
	my $qry_length = $mummer[5];
	my $percent_id = $mummer[6];
	#my $percent_cov_ref = $mummer[9];
	#my $percent_cov_qry = $mummer[10];
	my $qry_asmbl_length = $mummer[8];
        my $ref_asmbl_length = $mummer[7];

	my $aln = $doc->createAndAddSequencePairAlignment( 'refseq'        => $ref_name,
							   'compseq'       => $qry_name,
                                                           'complength'    => $qry_asmbl_length,
                                                           'reflength'     => $ref_asmbl_length
							 ); 

	my $s = $doc->createAndAddSequencePairRun( 'alignment_pair'   => $aln,
						   'refpos'           => $ref_start,
						   'runlength'        => $ref_length,
						   'comppos'          => $qry_start,
						   'comprunlength'    => $qry_length,
                                                   'percent_identity' => $percent_id,
						  );
    }
}

sub print_usage {


    print STDERR "SAMPLE USAGE:  mummer2bsml.pl -m mummer_coords -b bsml_repository -t 1 -o output_file\n";
    print STDERR "  --mummer_coords     = mummer output file\n";
    print STDERR "  --bsml_dir(-b)      = BSML repository directory\n";
    print STDERR "  --output            = bsml output file\n";
    print STDERR "  --mummer_type(-t)   = type of mummer output (1=nucmer, 2=promer)\n";
    print STDERR "  --help = This help message.\n";
    exit 1;

}
