#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

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

    create_pseudomolecules.pl --output_dir <outdir>  --contig_input <contig file> --strain <strain name> [--input_file <reference accession ids file> --database <genBank database> --format <reference file format> --config_param <nucmer config params> --log <log file> --debug <debug level> --help <usgae>]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --output_dir  	= /path/to/output_dir. This is the directory where all output 
		    	  files generated during the script execution will be stored.

    --contig_input 	= /path/to/contig_file or contig_list_file. This is a file 
			  containing the strain genome contigs in fasta format 
			  or paths to contig files.
    
    --strain      	= Name of the strain used for naming the output files.

   [--config_param	= Configuration parameters to execute nucmer.
		    	  e.g. -c 100 -maxmatch

    --input_file  	= /path/to/input_file. This is a tab-delimited file containing
		    	  genBank accession ids or paths to reference genome files and 
			  a serial group number indicating which reference sequences 
			  should be concatenated together for analysis .
		    	  e.g.CP000803.1	1
			      /path/to/ref_file	1 	
			      CP000804.1	2
			      CP000806.1	2

    --database 	  	= Name of the genBank database used to download the reference 
		    	  genomes.e.g.nucleotide

    --format      	= Format in which the reference genome sequence is downloaded
		    	  e.g.FASTA
    
    --linker_sequence	= This sequence will be inserted after each sequence stitched together into
			  a pseudomolecule.  The default sequence (NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN) will be
			  inserted if this isn't specified and contains 6-frame translational stop codons.  This
			  can be an empty string.

    --log 	  = /path/to/log_file. Log file. Optional

    --debug       = Debug level. Optional

    --help]       = Help message, script usage. Optional

=head1 DESCRIPTION

The program creates pseudomolecules from contigs for annotation purposes. Following steps are taken in pseudomolecule creation:
1. If no reference file is passed then concatenate the contigs based on length - longest to shortest to create a pseudomolecule.
2. If a reference accession file containing accession ids or file paths is specified then download reference genome sequences from genBank.
   Always provide absolute paths for the reference genome files
3. Concatenate reference genomes based on the groups provided by the user.
4. Run nucmer to create .delta file by aligning reference genome with input contigs file
5. Run show-coords to create .coords file for the alignment.
6. Run reference_genome_match_tiler.pl to produce .map file and  for the contigs that mapped with the reference.
7. Run create_fasta_pseudomolecules.pl to produce a pseudomolecule and unmapped contigs fasta file using the .map file.
8. Run cleanFasta.pl to clean up the pseudomolecule file and generate a FASTA format file.
9. If there is more than 1 reference group then re-run all the steps above starting from step 4 with unmapped contigs file as the 
   new contig file and next group of reference genomes as the reference sequence. 
10. At the end generate the pseudomolecule of unmapped contigs using create_fasta_pseudomolecules.pl

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
my ($results, $contig_file, $groupnum, $orig_contig_file);
my @contents;

#################################################
# MAIN PROGRAM					#
#################################################
$results = GetOptions (\%options,
                'input_file|i=s',
		'output_dir|o=s',
		'database|d=s',
		'format|f=s',
		'contig_input|c=s',
		'strain|s=s',
		'config_param|p=s',
		'nucmer_exec|n=s',
		'coords_exec|e=s',
		'linker_sequence|k=s',
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
	next if ($line =~ /^\s*$/);
	my @fline = split(/\t/,$line);
## Hash keyed on group number, value is the accession number.
	if(exists($accHash{$fline[1]})) {
		$accHash{$fline[1]} = $accHash{$fline[1]}.",$fline[0]";	
	} else {
		$accHash{$fline[1]} = $fline[0];
	}
}

if(keys(%accHash) == 0) {
	$logger->logdie("Empty reference accession file $options{'input_file'} passed.");
} 

$contig_file = $orig_contig_file;

for my $group (sort {$a<=>$b} keys %accHash) {
	my @contig_file_size = &read_file($contig_file);
	my $file_size = @contig_file_size;
	if ($file_size > 0) {
		chomp($accHash{$group});
		system("$Bin/fetch_genbank --output_dir=$options{output_dir} --database=$options{database} --query=$accHash{$group} --format=$options{format}");		       
		my @acc = split(/,/,$accHash{$group});
		my $ref_file = $options{output_dir}."/reference_grp".$group.".".$options{format};
		foreach my $id (@acc) {
			my $filename;
			if ($id =~ /\/|\\/g) {
				$filename = $id;
			} else {
				$filename = $options{output_dir}."/reference_".$id.".".$options{format};
			}
			if (-e $filename) {
				system("cat $filename >> $ref_file");
			} else {
				$logger->logdie("Unable to download $id and create $filename or provided $filename does not exist");
			}
		}
		if (-s $ref_file) {
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
			my $unmapped_contig_file = $options{output_dir}."/unmapped_".$group.".fsa";
			system("$Bin/create_fasta_pseudomolecules --input_fasta_file=$contig_file --map_file=$map_file --output_file=$pseudo_file  --unmapped_output=$unmapped_contig_file --linker_sequence=$options{'linker_sequence'}");
			system("$Bin/clean_fasta $pseudo_file");
			system("$Bin/clean_fasta $unmapped_contig_file");
			$groupnum = $group;
			$contig_file = $unmapped_contig_file;
		} else {
			$logger->logdie("$ref_file file is empty. Some error in creating the reference genomes file");
		}
	}
}	

# The last unmapped set of contigs will be the last pseudomolecule
if(-e $contig_file) {
	my @unmapp = &read_file($contig_file);
	my $size_arr = @unmapp;
	if ($size_arr > 0) {
		my $unmap_file = $options{output_dir}."/".$options{strain}."_unmapped.map";
		open(OUTFH, "> $unmap_file")  or $logger->logdie("Could not open $unmap_file file for writing\n");
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
		system("$Bin/create_fasta_pseudomolecules --input_fasta_file=$contig_file --map_file=$unmap_file --output_file=$unmap_pseudo --linker_sequence=$options{'linker_sequence'}");
		system("$Bin/clean_fasta $unmap_pseudo");
	}
}


#########################################################
# SUBROUTINES						#
#########################################################

# Subroutine to check the supplied paramaters are correct
sub check_parameters {
        my $options = shift;

## make sure output directory, contig file and strain name were passed
        unless ($options{output_dir} && $options{contig_input} && $options{strain}) {
		$logger->logdie("All the manadatory parameters should be passed");
	}

## make sure the output directory exists
        if (! -e "$options{output_dir}") {
		$logger->logdie("The $options{output_dir} output directory passed could not be read or does not exist");
       	}
## make sure the contig file exists and is readable
	if (! -e "$options{contig_input}") {
		$logger->logdie("The $options{contig_input} contigs input file passed could not be read or does not exist");
	} else {
		my @ctg_ip = &read_file($options{'contig_input'});
		foreach my $content (@ctg_ip) {
			chomp($content);
			next if ($content =~ /^\s*$/);
			next if ($content =~ /^#/);
			if($content =~ /^>/) {
## If the contig_input is a multi fasta file containing contigs				
				$orig_contig_file = $options{'contig_input'};
				last;
			} elsif ($content =~ /\//g) {
## If the contig_input is a list file containing paths to contig files 
				&concat_contigs();
				last;
			} else {
## Else the contig_input does not contain correct data - neither fasta sequence nor list of file paths
				$logger->logdie("The $options{contig_input} file is neither a fasta file nor a list file of paths. Incorrect input");
			}
		}
## Assign default linker sequence if not supplied	
		$options{'linker_sequence'} = 'NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN' unless ( defined $options{linker_sequence});
## make sure the input file exists and is readable
		if ($options{'input_file'} && (-s $options{'input_file'})) {
			if(-e "$options{input_file}") {
				my @fetch_opts = qw(database format);
				foreach my $fetch_opt( @fetch_opts ) {
					if( !$options->{$fetch_opt} ) {
						$logger->logdie("The --input_file was passed but required $fetch_opt was missing");
					}
				}
			} else {
				$logger->logdie("The $options{input_file} file passed could not be read or does not exist");
			}
		} else {
			my $pmarks_file = "";
			$pmarks_file = $options{output_dir}."/".$options{strain}.".pseudomolecule.1.fasta.pmarks";
## If input file containing reference ids is not passed then sort the contigs file based on contigs length - longest to shortest 
			my $contig_out = $options{output_dir}."/".$options{strain}.".pseudomolecule.1.fasta";
			&sort_contigs($orig_contig_file,$contig_out, $options{'linker_sequence'}, $pmarks_file);
			system("$Bin/clean_fasta $contig_out");
			exit;
		}
	}
}

sub concat_contigs {
	$orig_contig_file = $options{output_dir}."/".$options{strain}.".multi.fasta";
	my @contig_files = &read_file($options{'contig_input'});
	foreach my $path (@contig_files) {
		chomp($path);
		next if ($path =~ /^\s*$/);
		next if ($path =~ /^#/);
		if(-e $path) {
			system("cat $path >> $orig_contig_file");			
		}
	}
	system("$Bin/clean_fasta $orig_contig_file");	
}


## Subroutine to arrange the contigs from longest to shortest in a pseudomolecule
sub sort_contigs {
	my ($in_file,$out_file,$linker,$pmarks_file) = @_;
	my $seq = "";
	my $header;
	my $flag = 1;
	my $linker_start = 0;
	my $linker_len = length($linker);
	my $linker_end = $linker_len;
	my $linker_point = $linker_end;
	my %contig_hash;
	open(OUTFH, "> $out_file")  or $logger->logdie("Could not open $out_file file for writing\n");
	open(POUT, "> $pmarks_file")  or $logger->logdie("Could not open $pmarks_file file for writing\n") if (defined $pmarks_file);
	
## Read input contig file
	my @contigs = &read_file($in_file);
	foreach my $line (@contigs) {
		chomp($line);
		if($line =~ /^>(.*)/) {
			unless($flag) {
				chomp($seq);
				my $seq_len = length($seq);
				$contig_hash{$seq_len}{$header} = $seq;
				$seq = "";
			}
			$flag = 0;
			$header = $1;
		} else {
			## skip it if it is just whitespace
			next if ($line =~ /^\s*$/);
			$seq .= $line;
		}
	}
## Concatenate the last contig sequence
	my $seq_len = length($seq);
	$contig_hash{$seq_len}{$header} = $seq if(defined($header));
## Header for the pseudomolecule	
	my ($out_base,$out_dir,$out_ext) = fileparse($out_file,qr/\.[^.]*/);
	print OUTFH ">$out_base\n";
	print POUT ">$out_base\n";
## Adding linker sequence between each contig and at the beginning and end of the pseudomolecule
	print OUTFH "$linker";
	foreach my $seqlen (sort {$b<=>$a} keys %contig_hash) {
		foreach my $head (keys %{$contig_hash{$seqlen}}) {
			print POUT "$linker_start\t$linker_end\n" if (defined $pmarks_file);
			print OUTFH "$contig_hash{$seqlen}{$head}$linker";
			$linker_start = $linker_point + $seqlen; 
        		$linker_end = $linker_start + $linker_len;
			$linker_point = $linker_end;
		}
	}
	print POUT "$linker_start\t$linker_end\n" if (defined $pmarks_file);
	close(OUTFH);
	close(POUT) if (defined $pmarks_file);
}

## Subroutine to read files
sub read_file {
	my $filename = shift;
	my @lines;
	open(FH , "< $filename")  || $logger->logdie("Could not open $filename file for reading.$!");
	@lines = <FH>;
	close(FH);
	return(@lines);
} 
