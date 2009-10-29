#!/usr/bin/perl

=head1  NAME

metagene2bsml.pl  - converts raw metagene output into BSML

=head1 SYNOPSIS

USAGE:  metagene2bsml.pl -i raw_metagene_output.txt -f /path/to/fasta_input.fsa -o /path/to/outfile.bsml [-l /path/to/logfile.log] --project project_name --id_repository /path/to/id_repository

=head1 OPTIONS

B<--input_file,-i>
    The input raw output from metagene.

B<--fasta_input,-f>
    The full path to the fasta sequence file that was input to metagene.
    
B<--output,-o>
    The full path to the BSML output file to be created. 

B<--debug,-d>
    Debug level.  Use a large number to turn on verbose debugging.

B<--log,-l>
    Log file

B<--help,-h>
    This help message

=head1   DESCRIPTION

Create BSML document from metagene output. Will extract nucleotide sequences for the
predicted ORFs.

=head1  CONTACT
        Brett Whitty
        bwhitty@tigr.org

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use File::Basename;
use Chado::Gene;
use BSML::GenePredictionBsml;
use Ergatis::IdGenerator;
use Ergatis::Logger;



my %options = ();
my $results = GetOptions (\%options,
                          'input_file|i=s',
                          'fasta_input|f=s',
                          'output_directory=s',
                          'output|o=s',
                          'num_seqs|n=s',
                          'id_repository=s',
                          'analysis_id=s',
                          'project=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				 'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();
                     
if ($options{'help'}) {
    pod2usage(verbose => 2);
}

if (! $options{'analysis_id'}) {
    $options{'analysis_id'} = 'metagene_analysis';
}

if (!$options{'id_repository'}) {
    pod2usage("must provided --id_repository");
}
if (!$options{'project'}) {
    pod2usage("project name must be provided with --project");
}
my $id_gen = Ergatis::IdGenerator->new('id_repository' => $options{'id_repository'});

unless ($options{'input_file'}) {
    pod2usage("raw metagene output must be provided with --input_file");
}

unless ($options{'output'}) {
    pod2usage("must specify an output prefix with --output");
}
my $output_dir = $options{'output_directory'} || '';
if ($output_dir eq '') {
    $output_dir = '.';
}

my $output_fsa_file = $output_dir."/".basename($options{'output'}, ".bsml").".fsa";

my $orfs = {};
$options{'num_seqs'} = 1 if(!defined $options{'num_seqs'});

my $in_fh;
open $in_fh,"$options{'input_file'}";
my $orf_count = 0;
my $seq_id;
while (<$in_fh>) {
    chomp;    
    if($_ =~ /^\#\s+(\S+)/){
	if($1 ne 'gc' && $1 ne 'self:'){
	    print $1," ",scalar(keys %$orfs)," $orf_count\n";
	    $seq_id = $1;
	}
    }
    else{
	my @elts = split("\t", $_);
	if(scalar (@elts) == 11){
	    my ($geneid,
		$startpos, 
		$endpos, 
		$strand, 
		$frame, 
		$complete,
		$genescore, 
		$usedmodel,
		$rbsstart,
		$rbsend,
		$rbsscore
		) = @elts;
	    
	    ## set flags for partial ORFs
	    my ($five_prime_partial, $three_prime_partial) = split(//,$complete);
	    
	    push (@{$orfs->{$seq_id}},  {    
		'model'       =>  $usedmodel,
		'startpos'    =>  $startpos,
		'endpos'      =>  $endpos,
		'complement'  =>  ($strand eq '+') ? 0 : 1,
		'frame'       =>  $frame,
		'score'       =>  $genescore,
		'5_partial'   =>  $five_prime_partial,
		'3_partial'   =>  $three_prime_partial,
		'orf_id'      =>  $geneid
		});
	    $orf_count++;
	}        
	else{
	    print $_,"\n";
	}
    }
}

$id_gen->set_pool_size('ORF' => $orf_count, 'CDS' => $orf_count);

my $orfcount=0;

my $doc = new BSML::GenePredictionBsml( 'glimmer3', $options{'fasta_input'} );
my @bsmldocs;
my $numseqs=0;
my $outputnum=0;
foreach my $seq_id (keys(%{$orfs})) {
    if($numseqs >= $options{'num_seqs'}){
	$outputnum++;
	my $ofile = $output_dir."/".$options{'output'}.".$outputnum.bsml";
	$doc->writeBsml($ofile);
	my $doc = new BSML::GenePredictionBsml( 'glimmer3', $options{'fasta_input'} );
	$numseqs=1;
    }
    
    my $addedTo = $doc->setFasta($seq_id, $options{'fasta_input'});
    die "$seq_id was not a sequence associated with the gene" unless($addedTo);

    my @orfs;
    foreach my $orf (@{$orfs->{$seq_id}}) {
                
        my $orf_id = $options{'project'}."ORF".$orfcount++;
                     #$id_gen->next_id(
                     #   'project' => $options{'project'}, 
                     #   'type'    => 'CDS'
                     #                );

        $orf->{'id'} = $orf_id;

	#Create some genes and push them ontot he $genes array
	my $tmp = new Chado::Gene( $orfcount++,
				   $orf->{'startpos'}-1, $orf->{'endpos'}, ($orf->{'complement'} > 0) ? 0 : 1,
				   $seq_id);
	
	foreach my $type(qw(exon CDS transcript polypeptide)) {
	    $tmp->addFeature("ORF_$orfcount",
			     $orf->{'startpos'}-1, $orf->{'endpos'}, ($orf->{'complement'} > 0) ? 0 : 1,
			     $type);
	}
	
	my $count = $tmp->addToGroup($tmp->getId, { 'all' => 1 });
    }
    $numseqs++;
}

if($numseqs >= $options{'num_seqs'}){
    $outputnum++;
    my $ofile = $output_dir."/".$options{'output'}.".$outputnum.bsml";
    $doc->writeBsml($ofile);
    my $doc = new BSML::GenePredictionBsml( 'metagene', $options{'fasta_input'} );
    $numseqs=1;
}
