#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use MG::ParseSeq;
use MG::RunParse;
use Ergatis::Logger;
use Bio::SeqIO;

my %options = ();
my $results = GetOptions (\%options, 
			  'output_prefix|o=s',
			  'infile|i=s', 
			  'blast_query|f=s',
			  'algo|a=s',
			  'db_file|d=s',
			  'out_format|m=s',
			  'gene_prefix|p=s',
			  'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				 'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

unless ($options{infile} && $options{blast_query}) {
    pod2usage({-exitval=>0, -verbose => 2, -output => \*STDOUT});
}

if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

my ($hitref, $seqref, %seg);
%seg = parse_search($options{algo},$options{infile},'1e-6',undef,undef,undef);

my %hregion;

## Convert blastx output to genbank format.
my $out = Bio::SeqIO->new( -file    =>  ">$options{'output_prefix'}\.gbk",
			               -format  =>  $options{'out_format'} );

## Convert blastx output to tbl format.
open TBL, ">$options{'output_prefix'}\.txt" or die $!;
my $j = 0;
foreach my $r(keys %seg) {
    my $i = 0;
    my (%raln,%fxtran);
    foreach my $g (@{$seg{$r}}) {
	$raln{$g->[0]} = 1;
    }
    my $rname = (split(/\/|\.|\|/, $r))[-1];
    create_db("$options{infile}\_$rname\.lib",$options{db_file}, keys %raln);
    create_db("$options{infile}\_$rname\.query",$options{blast_query}, $r);
    my %seqobj = runparse_fastx("$options{infile}\_$rname\.query","$options{infile}\_$rname\.lib",$rname);
    foreach my $query_name (keys %seqobj) {
	my $read = $query_name;
	my $obj = $seqobj{$query_name};
	next unless ($obj->get_SeqFeatures);
	foreach my $feat ($obj->get_SeqFeatures) {
	    $j ++;
	    my $newacc = $options{'gene_prefix'}.sprintf("%06s",$j);
	    $feat->add_tag_value('protein_id',$newacc);
	    print TBL join("\t", $read,$newacc,$feat->start,$feat->end,
			   $feat->strand,$feat->frame,'',$feat->score,'','fastx'),"\n";
	}
	$out->write_seq($obj);
    }
}

close TBL;
