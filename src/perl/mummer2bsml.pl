#!/usr/local/bin/perl

=head1  NAME 

mummer2bsml.pl  - convert info stored in mummer coord files into BSML documents

=head1 SYNOPSIS

USAGE:  mummer2bsml.pl -m mummer_coord.txt -o nucmer.bsml -t 1

=head1 OPTIONS

=over 4

=item *

B<--mummer_coords,-m> [REQUIRED]  mummer coordinate file

=item *

B<--output,-o> [REQUIRED] output BSML file containing mummer coordinate information

=item *

B<--mummer_type,-t> [REQUIRED] Type of mummer files: 1 = nucmer, 2 = promer

=item *

B<--help,-h> This help message

=back

=head1   DESCRIPTION

mummer2bsml.pl is designed to convert information in mummer coordinate files 
into BSML documents.  There 2 types of mummer files: nucmer and promer.  The user 
can specify which type of mummer files to convert by --mummer_type(t) flag.  An 
argument of 1 means nucmer, while a 2 indicate the file is that of promer.  

Samples:

1. convert nucmer coordinate files "nucmer_coord.txt" into BSML doc 

   mummer2bsml.pl -m nucmer_coord.txt -o nucmer.bsml -t 1 

NOTE:  

Calling the script name with NO flags/options or --help will display the syntax requirement.


=cut


use strict;
#use Log::Log4perl qw(get_logger);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use English;
use BSML::BsmlBuilder;
use BSML::BsmlParserTwig;
use File::Basename;
use Pod::Usage;

umask(0000);

my %options = ();
my $results = GetOptions (\%options, 'mummer_coords|m=s', 'mummer_type|t=s', 'man', 'output|o=s', 'database|d=s', 'verbose|v', 'help|h') || pod2usage();


###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $mummer_file     = $options{'mummer_coords'};
my $output          = $options{'output'};
my $database          = $options{'database'};
my $output_dir      = dirname($output);
my $mummer_type     = $options{'mummer_type'};   # 1 = nucmer , 2 = promer
my $verbose    = $options{'verbose'};
#Log::Log4perl->init("log.conf");
#my $logger = get_logger();

&cmd_check();
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
	
	$qry_name =~ s/^$database\_//i; #strip leading database name for now
	$ref_name =~ s/^$database\_//i; #strip leading database name for now

	next if($ref_name eq $qry_name);

	

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
    $doc->createAndAddAnalysis("program" => "PROmer", "programversion" => '3.0', 'sourcename' =>$output,
                               "bsml_link_relation" => 'SEQ_PAIR_ALIGNMENTS', 'bsml_link_url' => '#BsmlTables');

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

	
	$qry_name =~ s/^$database\_//i; #strip leading database name for now
	$ref_name =~ s/^$database\_//i; #strip leading database name for now

	next if($ref_name eq $qry_name);	

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
    $doc->createAndAddAnalysis("program" => "NUCmer", "programversion" => '3.0', 'sourcename' =>$output,
                               "bsml_link_relation" => 'SEQ_PAIR_ALIGNMENTS', 'bsml_link_url' => '#BsmlTables');


}


sub cmd_check {
#quality check

    if( exists($options{'man'})) {
	pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT});
    }   
    
    if( exists($options{'help'})) {
	pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT});
    }
    
    if(!$mummer_file or !$mummer_type or !$output) {
	pod2usage({-exitval => 2,  -message => "$0: All the required options are not specified", -verbose => 1, -output => \*STDERR});    
    }


    if($mummer_type != '1' and $mummer_type != '2') {
	pod2usage({-exitval => 2,  -message => "$0: mummer_type can only be 1 or 2", -verbose => 1, -output => \*STDERR}); 
    }

    #check for presence of output directory
    if(! -d $output_dir) {
	mkpath($output_dir) or die "Unable to create $output_dir.  Aborting...\n";
	#chmod 0777, $output_dir;
    }

}

