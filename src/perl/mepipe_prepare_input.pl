#!/usr/local/bin/perl

=head1  NAME

mepipe_prepare_input.pl - prepares input files for molecular evolution pipeline

=head1 SYNOPSIS

USAGE: mepipe_prepare_input.pl
        --input_prefix=gene_id
		--input_path=/path/to/directory/containing/input/files/
        --temp_path=/path/to/store/temporary/output
		--output_path=/path/to/store/final/output
		--emboss_bin_dir=/path/to/emboss/bin
        --debug=4

=head1 OPTIONS

B<--input_prefix,-p>
	The file prefix for the gene/cog to be analyzed.
	(NB: This was referred to as 'gene name' in the original script.)

B<--input_path,-i>
    The path to the directory containing required input files.

B<--temp_path,-t>
	The temporary path for execution steps. (eg: $;TEMP_DIR$;)
	
B<--output_path,-o>
	The path to the directory to store output.
	
B<--emboss_bin_dir, -b>
	Full path containing EMBOSS binaries (specifically tranalign).
	
B<--debug,-d>
    Debug level.  Use a large number to turn on verbose debugging.

B<--log,-l>
    Log file

B<--help,-h>
    This help message

=head1   DESCRIPTION

*** TO BE COMPLETED ***

=head1 INPUT

*** TO BE COMPLETED ***

This script takes as input a file with a list of gene names (one per line)
and, for each gene, an amino acid clustal alignment (<geneName.aa>) 
and respective nucleotide seqs (<geneName.seq>) in fasta format,
converts the aa alignment into nucleotide alignment, and creats PAUP 
and PAML input files
Runs PAUP and PAML

NOTE: for each gene, sequence names in the amino acid alignment file
and corresponding nucleotide seq file have to be identical; 
names should contain no spaces.

=head1 OUTPUT

*** TO BE COMPLETED ***

=head1 CONTACT

	Author:
	Joana Silva
	jsilva@tigr.org

	Modifications for running with workflow:
    Brett Whitty
    bwhitty@tigr.org

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use File::Copy;
BEGIN {
use Workflow::Logger;
}

umask(0000);

my %options = ();
my $results = GetOptions (\%options,
						  'input_prefix|p=s',
                          'input_path|i=s',
						  'temp_path|t=s',
                          'output_path|o=s',
#						  'paml_bin_dir=s',
#						  'paup_bin_dir=s',
						  'emboss_bin_dir|b=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || pod2usage();

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = Workflow::Logger::get_logger();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDERR} );
}

#&check_parameters(\%options);

my $aaAlign;
my $nuclSeqs;
my $nuclAlign;
my $geneName;

my $outgroup;
my $batchfile;
my @speciesNames;
my $noSpecs;
my $aaAlignLength;
my $PAMLrunmode;

foreach my $dir((
#		'paml_bin_dir', 
		'emboss_bin_dir', 
#		'paup_bin_dir',
		'temp_path',
		'input_path',
		'output_path', 
	        )) {
		if (! -d $options{$dir}) {
			die "'$dir' specified as '".$options{$dir}."' does not exist.\n";
		}
		unless ($options{$dir} =~ /\/$/) {
			$options{$dir} .= "/";
		}
}

	chdir($options{'temp_path'});
	
	$geneName = $options{'input_prefix'};
    $aaAlign = "$geneName.aln";
    $nuclSeqs = "$geneName.seq";
    $nuclAlign = "$geneName.nuc.aln";

	## Copy input files to the work directory
	copy($options{'input_path'}.$aaAlign, $options{'temp_path'}.$aaAlign); 
	copy($options{'input_path'}.$nuclSeqs, $options{'temp_path'}.$nuclSeqs); 
	
	## create the nexus formatted alignment file
	clustal2nexus($aaAlign);
	
	## read in species names from the nexus file to get a count
	## (the global variable $noSpecs is required by checkNucFastaformat())
	&readSpeciesNames($options{'temp_path'}.$aaAlign.".nex");
 	$noSpecs = scalar @speciesNames;
	
	#checks that nuc sequences are in the same order as seq in the aa alignment file
	&checkNucFastaformat($options{'temp_path'}.$nuclSeqs);     
                                         
	#tranalign part of the EMBOSS package, creates NA alignment from AA alignment
   system($options{'emboss_bin_dir'}."tranalign -asequence ".$options{'temp_path'}."$nuclSeqs -bsequence ".$options{'temp_path'}."$aaAlign -outseq ".$options{'temp_path'}."$nuclAlign 2>/dev/null");

   # convert fasta nucleotice alignment (gapped fasta) to PAML format ($geneName.nuc.aln.PAMLseq)
   &convertNuclAlignToPAML($options{'temp_path'}.$nuclAlign, $noSpecs, $aaAlignLength);

	## write species count to a file so that we can use it later	
#   open (OUT, ">".$options{'temp_path'}."$geneName.count") || die $!;
#   print OUT $noSpecs."\n";
#   close OUT;
   
exit();
						 
# end main
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

sub convertNuclAlignToPAML{  #($ nucl alignfile $noSpecs, $aaAlignLength)
	my $fastafile = "$_[0]";
	my $pamlfile = "$_[0].PAMLseq";
	my $noTaxa =  $_[1]; 
	my $nucAlignLen =  3 * $_[2]; #align length passed in $_[2] is for aa

	## THERE WILL PROBABLY BE SOME PATH PROBLEMS HERE
	open (PAML, ">$pamlfile");
	print PAML " $noTaxa $nucAlignLen\n";

	open (FASTA, "$fastafile");

	while (<FASTA>) {
    	print PAML "$_";
	}

	close FASTA;
	close PAML;
}


sub checkNucFastaformat {

    my $nucfile = shift @_;
	
    my $tempfile = "$nucfile.temp";
    my @names;
    my @nucseqs;
    my $count=0;
    my $line;
    my @dummyarr;

#read in nuc seq

	## THERE WILL PROBABLY BE SOME PATH PROBLEMS HERE
    open (NUCFILE, "$nucfile") or die $!;

    while (<NUCFILE>) {
		$line = $_;
		last if ($line =~ /^>/);
    }

    chomp($line);
    @dummyarr = split ( ' ', $line); #removes spaces in the end of the name, and any spurious comments

    $names[$count] = substr( $dummyarr[0], 1 );
      
    while (<NUCFILE>) {
		$line = $_;
		if ($line =~ /^>/) {
	    	chomp($line);
	    	@dummyarr = split ( ' ', $line); #removes spaces in the end of the name, and any spurious comments
	    	$names[++$count] = substr($dummyarr[0], 1 );
	    	next;
		}
		$nucseqs[$count] .= $line;
    }
    close NUCFILE;


#write seqs out to file in the order of the alignment file

	## THERE WILL PROBABLY BE SOME PATH PROBLEMS HERE
    open (TEMPFILE, ">$tempfile") or die $!;

    for (my $i = 0; $i <= $count ; $i++) {
		for  (my $j = 0; $j <= $count ; $j++) {

		    if ( "$speciesNames[$i]" eq "$names[$j]")  {  
				print TEMPFILE ">$names[$j]\n";
				print TEMPFILE "$nucseqs[$j]";
				last;
		    }
		}
    }

    close TEMPFILE;

	## THERE WILL PROBABLY BE SOME PATH PROBLEMS HERE
    move("$nucfile", "$nucfile.original");
	## THERE WILL PROBABLY BE SOME PATH PROBLEMS HERE
    move("$tempfile", "$nucfile");
}


sub readSpeciesNames {   #gets passed a nexus file
	
	my $infile = shift @_;

    my $count=0;
    my @dummyarr;
    my @arr;

	## THERE WILL PROBABLY BE SOME PATH PROBLEMS HERE
	open (NEXUS, "$infile") or die $!;
    while (<NEXUS>) {
	chomp;
	@arr= split(' ', $_);
	if ($arr[0]) {
	    last if ($arr[0] eq "TAXLABELS");
	}
    }
    while (<NEXUS>) {
		chomp;
		@arr= split(' ', $_);
		last if ($arr[0] eq ";");
		$speciesNames[$count++]= substr($arr[0], 1, (length($arr[0]) - 2) );
    }
    while (<NEXUS>) {
		chomp;
		@arr= split(' ', $_);
		if ($arr[1]) {
	    	if ($arr[1] =~ /^NCHAR/) {
				@dummyarr = split(/=/, $arr[1]);
				$aaAlignLength = substr($dummyarr[1], 0, (length($dummyarr[1])-1) );
				last;
	    	}
		}
    }
 
    close NEXUS;
}


sub clustal2nexus {

	my $infile = shift @_;
	
	my @seqNames= "";
	my @seqs = "";
	my $speciesNo= 0;
	my @array1="";
	my $counter = 0;
	my $lineNo = 0;
	my $len = 0;
	my $alignLen;

	my $outfile = $infile.".nex";

	open (INFILE, $infile) or die $!;

	my $line = <INFILE>;
	$lineNo++;

	if ($line !~ /^CLUSTAL/) { 
   		die("ERROR! Input file is not in CLUSTAL format\n");
	} 

	while (<INFILE>) {
    	$lineNo++;
	    chomp ($line = $_);
    	$len = length ($line);
	    last if ($len == 0);
	} #skips any extra lines in header

	while (<INFILE>) {
    	$lineNo++;
	    chomp ($line = $_);
    	last if (length ($line) != 0 );
	} #skips empty lines to the first seq

	while ( length ($line) != 0 ) {
    	last if ($line =~ /^\s/);
	    @array1 = split (' ', $line);
    	$seqNames[$speciesNo] = $array1[0];
	    $seqs[$speciesNo] = $array1[1];
    	$speciesNo++;
	    chomp ($line = <INFILE>);
    	$lineNo++;
	} #reads first piece of alignment, gets number of sequences

	while (<INFILE>) {
    	$lineNo++;
	    if ($_ =~ /^\s/ ) {
			$counter = 0; 
			next; #skips empty lines and lines starting with an empty space#
    	}
	    chomp ($line = $_);
    	@array1 = split (' ', $line);
	    $seqs[$counter] .= $array1[1];
    	$counter++;
	}
	#Finished reading Clustal alignment file
	close INFILE;

	$alignLen = length($seqs[0]);

	for (my $i = 1; $i < $speciesNo; $i++) {
    	if ( length($seqs[$i]) != $alignLen ) {
			die("ERROR: some sequences are of different length\n");
	    }
	}
	open (OUTFILE, ">$outfile");

	printf (OUTFILE "#NEXUS\n[Translated from Clustalw alignment: $infile]\n\n\n");	
	printf (OUTFILE "BEGIN TAXA;\n");	
	printf (OUTFILE "\tDIMENSIONS  NTAX=$speciesNo;\n");	
	printf (OUTFILE "\tTAXLABELS");	

	for (my $i = 0; $i < $speciesNo; $i++) {
    	printf (OUTFILE "\n\t'$seqNames[$i]'");	
	}
	printf (OUTFILE "\n;\nEND;\n\n\n");	

	printf (OUTFILE "BEGIN CHARACTERS;\n");	
	printf (OUTFILE "\tDIMENSIONS  NCHAR=$alignLen;\n");	
	printf (OUTFILE "\tFORMAT DATATYPE=PROTEIN  MISSING=? GAP=-;\n");	
	printf (OUTFILE "MATRIX\n");	

	for (my $i = 0; $i < $speciesNo; $i++) {
    	printf (OUTFILE "\n'$seqNames[$i]'    $seqs[$i]" );	
	}
	printf (OUTFILE "\n;\nEND;\n");	

	close OUTFILE;
}
