#!/usr/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

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

B<--mummer_type,-t> [REQUIRED] Type of mummer files: 1 = full promer coords, 2 = collapsed promer coords, 3 = collapsed nucmer coords


=item *

B<--debug,-d> Debug level.  Use a large number to turn on verbose debugging. 

=item *

B<--log,-l> Log file

=item *

B<--class,-c> ref/comp sequence class type.  Default class is 'assembly'

=item *

B<--help,-h> This help message

=back

=head1   DESCRIPTION

mummer2bsml.pl is designed to convert information in mummer coordinate
files into BSML documents.  The user can specify which type of mummer
files to convert by --mummer_type(t) flag. 

Samples:

1. convert nucmer coordinate files "nucmer_coord.txt" into BSML doc 

   mummer2bsml.pl -m nucmer_coord.txt -o nucmer.bsml -t 3 

NOTE:  

Calling the script name with NO flags/options or --help will display the syntax requirement.


=cut


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use English;
use File::Basename;
use Pod::Usage;
BEGIN {
use Ergatis::Logger;
use BSML::BsmlBuilder;
use BSML::BsmlParserTwig;
}

my %options = ();
my $results = GetOptions (\%options, 
			  'mummer_coords|m=s', 
			  'mummer_type|t=s', 
			  'output|o=s', 
			  'database|d=s', 
			  'log|l=s',
			  'debug=s',
			  'class|c=s',
			  'help|h') || pod2usage();


my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});

$logger = Ergatis::Logger::get_logger();

&check_parameters(\%options);

my $class;
if (!defined($options{'class'})){
    $logger->logdie("class was not defined");
}
else{
    $class = $options{'class'};
}


my $doc = BSML::BsmlBuilder->new();
if($options{'mummer_type'} == 1) {
    parse_promer_full_coords($options{'mummer_coords'});
}elsif($options{'mummer_type'} == 2) {
    parse_mummer_collapsed_coords($options{'mummer_coords'});
}
elsif($options{'mummer_type'} == 3) {
    parse_nucmer_full_coords($options{'mummer_coords'});
}else {
    $logger->logdie("Bad options{'mummer_type'} $options{'mummer_type'}");
}

$doc->write($options{'output'});

sub parse_mummer_collapsed_coords {

    my $coords_file = shift;

    my @mummer;
    open (MUMMER, $coords_file) or die "Unable to open \"$coords_file\" due to $!";
    while (my $line = <MUMMER>) {
	chomp($line);
	@mummer = split("\t", $line);	
	#In Mummer ref seq refers to query

	my $ref_start   = $mummer[0];
	my $ref_end     = $mummer[1];
	my $qry_start   = $mummer[2];
        my $qry_end     = $mummer[3];
	my $ref_length  = $mummer[4];
	my $qry_length  = $mummer[5];
	my $ref_asmbl_length = $mummer[6];
	my $qry_asmbl_length = $mummer[7];
	my $ref_name    = $mummer[10];
        my $qry_name    = $mummer[11];

	$qry_name =~ s/^$options{'database'}\_//i; #strip leading database name for now
	$ref_name =~ s/^$options{'database'}\_//i; #strip leading database name for now

	next if($ref_name eq $qry_name);

	my $complement= ($qry_start > $qry_end) ? 1 : 0;
	if($complement){
	    ($qry_start,$qry_end) = ($qry_end,$qry_start);
	}
	if($ref_start > $ref_end){
	    $logger->logdie("Unexpected reference coordinates $ref_start $ref_end\n");
	}
	$qry_start--;
	$ref_start--;	

	my $aln = $doc->createAndAddSequencePairAlignment( 'refseq'        => $ref_name,
                                                       'compseq'       => $qry_name,
                                                       'complength'    => $qry_asmbl_length,
                                                       'reflength'     => $ref_asmbl_length,
                                                       'class'         => 'match'
                                                       ); 
	my $s = $doc->createAndAddSequencePairRun( 'alignment_pair'   => $aln,
                                               'refpos'           => $ref_start,
                                               'refcomplement'    => 0,
                                               'runlength'        => $ref_length,
                                               'comppos'          => $qry_start,
                                               'comprunlength'    => $qry_length,
                                               'compcomplement'   => $complement,
                                               'class'            => 'match_part'
                                                  );
    }
}

sub parse_promer_full_coords {
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
	
	$qry_name =~ s/^$options{'database'}\_//i; #strip leading database name for now
	$ref_name =~ s/^$options{'database'}\_//i; #strip leading database name for now

	next if($ref_name eq $qry_name);

	my $complement= ($qry_start > $qry_end) ? 1 : 0;
	if($complement){
	    ($qry_start,$qry_end) = ($qry_end,$qry_start);
	}
	if($ref_start > $ref_end){
	    $logger->logdie("Unexpected reference coordinates $ref_start $ref_end\n");
	}
	$qry_start--;
	$ref_start--;	

	my $aln = $doc->createAndAddSequencePairAlignment( 'refseq'        => $ref_name,
							   'compseq'       => $qry_name,
                                                           'complength'    => $qry_asmbl_length,
                                                           'reflength'     => $ref_asmbl_length,
							   'class'         => $class
							 ); 
	my $s = $doc->createAndAddSequencePairRun( 'alignment_pair'   => $aln,
						   'refpos'           => $ref_start,
						   'refcomplement'    => 0,
						   'runlength'        => $ref_length,
						   'comppos'          => $qry_start,
						   'comprunlength'    => $qry_length,
						   'compcomplement'   => $complement
                                                  );
	#additional attributes
	$s->addBsmlAttr( 'percent_identity',  $percent_id);
	$s->addBsmlAttr( 'percent_similarity', $percent_sim);
	$s->addBsmlAttr( 'percent_stop',  $percent_stop);
	$s->addBsmlAttr( 'ref_frame', $frame_ref );
	$s->addBsmlAttr( 'query_frame', $frame_qry );

    }
}

sub parse_nucmer_full_coords {

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
	my $qry_start = $mummer[2];
        my $qry_end   = $mummer[3];
	my $ref_length = $mummer[4];
	my $qry_length = $mummer[5];
	my $percent_id = $mummer[6];
        my $ref_asmbl_length = $mummer[7];
	my $qry_asmbl_length = $mummer[8];

	$qry_name =~ s/^$options{'database'}\_//i; #strip leading database name for now
	$ref_name =~ s/^$options{'database'}\_//i; #strip leading database name for now

	next if($ref_name eq $qry_name);	

	my $aln = $doc->createAndAddSequencePairAlignment( 'refseq'        => $ref_name,
							   'compseq'       => $qry_name,
                                                           'complength'    => $qry_asmbl_length,
                                                           'reflength'     => $ref_asmbl_length,
							   'class'         => $class
							 ); 

	my $complement= ($qry_start > $qry_end) ? 1 : 0;
	if($complement){
	    ($qry_start,$qry_end) = ($qry_end,$qry_start);
	}
	if($ref_start > $ref_end){
	    $logger->logdie("Unexpected reference coordinates $ref_start $ref_end\n");
	}
	$qry_start--;
	$ref_start--;

	my $s = $doc->createAndAddSequencePairRun( 'alignment_pair'   => $aln,
						   'refpos'           => $ref_start,
						   'refcomplement'    => 0,
 						   'runlength'        => $ref_length,
						   'comppos'          => $qry_start,
						   'comprunlength'    => $qry_length,
						   'compcomplement'   => $complement
						   );
	$s->addBsmlAttr( 'percent_identity',  $percent_id);
    }
}


sub check_parameters{
    my ($options) = @_;
    
    if(!$options{'mummer_coords'} or !$options{'mummer_type'} or !$options{'output'}) {
	pod2usage({-exitval => 2,  -message => "$0: All the required options are not specified", -verbose => 1, -output => \*STDERR});    
    }


    if($options{'mummer_type'} != '1' and $options{'mummer_type'} != '2' and $options{'mummer_type'} != '3') {
	pod2usage({-exitval => 2,  -message => "$0: mummer_type can only be 1 or 2", -verbose => 1, -output => \*STDERR}); 
    }
}

