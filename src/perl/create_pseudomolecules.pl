#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
if 0; # not running under some shell

BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

###########################################################
# POD DOCUMENTATION                                       #
###########################################################
=head1 NAME

create_pseudomolecules.pl - program to generate pseudomolecules from contigs for annotation.

=head1 SYNOPSIS

    create_pseudomolecules.pl --input_file <reference accession ids file> --output_dir <outdir> --database <genBank database> --format <reference file format> --contig_file <contig file> --strain <strain name> --config_param <nucmer config params> [--log <log file> --debug <debug level> --help <usgae>]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --input_file  = /path/to/input_file. This is a tab-delimited file containing
		    genBank accession ids and a serial group number indicating which
		    reference sequences should be concatenated together for analysis .
		    e.g.CP000803.1	1
			CP000804.1	2

    --output_dir  = /path/to/output_dir. This is the directory where all output 
		    files generated during the script execution will be stored.

    --database 	  = Name of the genBank database used to download the reference 
		    genomes.e.g.nucleotide

    --format      = Format in which the reference genome sequence is downloaded
		    e.g.FASTA

    --contig_file = /path/to/contig_file. This is a fasta file containing the
		    strain genome contigs.
    
    --strain      = Name of the strain used for naming the output files.

    --config_param= Configuration parameters to execute nucmer.
		    e.g. -c 100 -maxmatch

   [--log 	  = /path/to/log_file. Log file. Optional

    --debug       = Debug level. Optional

    --help]       = Help message, script usage. Optional

=head1 DESCRIPTION

The program creates pseudomolecules from contigs for annotation purposes. Following steps are taken in pseudomolecule creation:
1. Download reference genome sequences from genBank.
2. Concatenate reference genomes based on the groups provided by the user.
3. Run nucmer to create .delta file by aligning reference genome with input contigs file
4. Run show-coords to create .coords file for the alignment.
5. Run reference_genome_match_tiler.pl to produce .map file and  for the contigs that mapped with the reference.
6. Run create_fasta_pseudomolecules.pl to produce a pseudomolecule and unmapped contigs fasta file using the .map file.
7. Run cleanFasta.pl to clean up the pseudomolecule file and generate a FASTA format file.
8. If there is more than 1 reference group then re-run all the steps above starting from step 3 with unmapped contigs file as the 
   new contig file and next group of reference genomes as the reference sequence. 
9. At the end generate the pseudomolecule of unmapped contigs using create_fasta_pseudomolecules.pl

=head1  CONTACT

Sonia Agrawal
sagrawal@som.umaryland.edu

=cut

use strict;
use File::Basename;
use FindBin qw($Bin);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use POSIX;
BEGIN {
        use Ergatis::Logger;
}

#################################################
# GLOBAL VARIABLES 				#
#################################################
my (%options, %accHash) = (); 
my ($results, $contig_file, $groupnum);
my @contents;

#################################################
# MAIN PROGRAM					#
#################################################
$results = GetOptions (\%options,
                'input_file|i=s',
		'output_dir|o=s',
		'database|d=s',
		'format|f=s',
		'contig_file|c=s',
		'strain|s=s',
		'config_param|p=s',
		'nucmer_exec|n=s',
		'coords_exec|e=s',
                'log|l=s',
                'debug|b=s',
                'help|h') || pod2usage();

## Display documentation
if( $options{'help'} ){
        pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## Getting the log file
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

## Make sure everything passed was correct
&check_parameters(\%options);

## Read the contents of the input accession ids file
@contents = &read_file($options{input_file});

foreach my $line (@contents) {
	chomp($line);
	next if($line =~ /^#/);
	my @fline = split(/\t/,$line);
## Hash keyed on group number, value is the accession number.
	if(exists($accHash{$fline[1]})) {
		$accHash{$fline[1]} = $accHash{$fline[1]}.",$fline[0]";	
	} else {
		$accHash{$fline[1]} = $fline[0];
	}
}

$contig_file = $options{contig_file};

for my $group (sort {$a<=>$b} keys %accHash) {
	chomp($accHash{$group});
	system("/usr/local/projects/ergatis/package-nightly/bin/fetch_genbank --output_dir=$options{output_dir} --database=$options{database} --query=$accHash{$group} --format=$options{format}");		       my @acc = split(/,/,$accHash{$group});
	my $ref_file = $options{output_dir}."/reference_grp".$group.".".$options{format};
	foreach my $id (@acc) {
		my $filename = $options{output_dir}."/reference_".$id.".".$options{format};
		if (-e $filename) {
			system("cat $filename >> $ref_file");
		}
	}
	my ($ref_base,$ref_dir,$ref_ext) = fileparse($ref_file,qr/\.[^.]*/);
	my ($contig_base,$contig_dir,$contig_ext) = fileparse($contig_file,qr/\.[^.]*/);
	my $nucmer_prefix = $options{output_dir}."/".$ref_base.$contig_base;
	system("$options{'nucmer_exec'} -p $nucmer_prefix $options{config_param} $ref_file $contig_file");
	my $nucout_file = "$nucmer_prefix.delta";
	my @nucout = &read_file($nucout_file);
	my $is_full = 0;
	foreach my $line (@nucout) {
		chomp($line);
		if($line =~ /^>/) {
			$is_full = 1;
			last; 
		}	
	}
	if($is_full == 0) {
		next;
	}
	my $delta_file = $options{output_dir}."/".$ref_base.$contig_base.".delta";
	my $coords_file = $options{output_dir}."/".$ref_base.$contig_base.".coords";
	system("$options{coords_exec} -T $delta_file > $coords_file");
	system("$Bin/reference_genome_match_tiler --mummer_coords_file $coords_file --min_match_length 100 --mummer_delta_file $delta_file --method nucmer --strain $options{strain} --pseudonum $group --output_dir $options{output_dir}");
	my $map_file = $options{output_dir}."/".$options{strain}."_".$group.".map";
	my $pseudo_file = $options{output_dir}."/".$options{strain}.".pseudomolecule.".$group.".fasta";
	$contig_file = $options{output_dir}."/unmapped_".$group.".fsa";
	system("$Bin/create_fasta_pseudomolecules --input_fasta_file=$options{contig_file} --map_file=$map_file --output_file=$pseudo_file  --unmapped_output=$contig_file");
	system("$Bin/clean_fasta $pseudo_file");
	$groupnum = $group;
}	

# The last unmapped set of contigs will be the last pseudomolecule
if(-e $contig_file) {
	my @unmapp = &read_file($contig_file);
	my $unmap_file = $options{output_dir}."/".$options{strain}."_unmapped.map";
	open(OUTFH, "> $unmap_file") or $logger->logdie("Could not open $unmap_file file for writing\n");
	my $molnum = $groupnum + 1;
	my $unmap_mol = "$options{strain}.pseudomolecule.$molnum";
	my $unmap_pseudo = $options{output_dir}."/".$options{strain}.".pseudomolecule.".$molnum.".fasta";
	foreach my $line (@unmapp) {
		chomp($line);
		if($line =~ /^>/) {
			$line =~ s/>//;
			print OUTFH "$unmap_mol\t" . "$line\t" . "+\n";
		}
	}
	close(OUTFH);
	system("$Bin/create_fasta_pseudomolecules --input_fasta_file=$options{contig_file} --map_file=$unmap_file --output_file=$unmap_pseudo");
	system("$Bin/clean_fasta $unmap_pseudo");
}


#########################################################
# SUBROUTINES						#
#########################################################

# Subroutine to check the supplied paramaters are correct
sub check_parameters {
        my $options = shift;

## make sure output directory, database, contig file, strain name, accession ids input file and format were passed
        unless ( $options{database} && $options{output_dir} && $options{format} && $options{input_file} && $options{contig_file} && $options{strain} ) {
		$logger->logdie("All the manadatory parameters should be passed");
	}

## make sure the output directory exists
        if (! -e "$options{output_dir}") {
		$logger->logdie("The output directory passed could not be read or does not exist");
        	
	}

## make sure the input file exists and is readable
	if (! -e "$options{input_file}") {
		$logger->logdie("The input file passed could not be read or does not exist");
	}

## make sure the contig file exists and is readable
	if (! -e "$options{contig_file}") {
		$logger->logdie("The contigs file passed could not be read or does not exist");
	}
}

## Subroutine to read files
sub read_file {
	my $filename = shift;
	my @lines;
	open(FH , "< $filename") || $logger->logdie("Could not open $filename file for reading.$!");
	@lines = <FH>;
	close(FH);
	return(@lines);
} 
