#!/usr/bin/perl -w

###########################################################
# POD DOCUMENTATION                                       #
###########################################################
=head1 NAME

insert_pmarks_in_scaffold.pl - program to replace a string of Ns within scaffolds with PMARKS.

=head1 SYNOPSIS

    insert_pmarks_in_scaffold.pl --output_dir <outdir>  --scaffold_input <scaffold file> --strain <strain name> [--linker_sequence <sequence> --log <log file> --help <usgae>]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --output_dir  	= /path/to/output_dir. This is the directory where all output 
		    	  files generated during the script execution will be stored.

    --scaffold_input 	= /path/to/scaffold_file or scaffold_list_file. This is a file 
			  containing the strain genome scaffolds in fasta format 
			  or paths to scaffold files.
    
    --strain      	= Name of the strain used for naming the output files.

   [--linker_sequence	= This sequence will replace a string of Ns within the sequencing gap of the scaffold.  
			  The default sequence (NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN) will be
			  inserted if this isn't specified and contains 6-frame translational stop codons.  This
			  can be an empty string.

    --log 	  	= /path/to/log_file. Log file. Optional

    --help]       	= Help message, script usage. Optional

=head1 DESCRIPTION

The program replaces a string of Ns that represent a sequencing gap within a scaffold with a PMARK
-If the string of N's is 70 or fewer nucleotides, then the entire string is replaced with a PMARK
-Otherwise strings of N's with greater than 70 nucleotides will have PMARKS on both sides of the string

=head1  CONTACT

Shaun Adkins
sadkins@som.umaryland.edu

=cut

use strict;
use lib ("/usr/local/projects/ergatis/package-latest/lib/perl5");
use File::Basename;
use FindBin qw($Bin);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use POSIX;
BEGIN {
        use Ergatis::Logger;
}

##### Prototypes #####

sub check_parameters($);
sub concat_files();
sub insert_pmarks($);
sub read_file($);


##### Main Program #####

my %options = ();
my $results = GetOptions (\%options,
		'output_dir|o=s',
		'scaffold_input|i=s',
		'strain|s=s',
		'linker_sequence|k=s',
                'log|l=s',
		'debug|b=s',
                'help|h') || pod2usage();

my $scaffold_file;
my $orig_scaffold_file;

## Cutoff ength of the string of N's for PMARK placement
my $length_cutoff = 70;

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
check_parameters(\%options);
$scaffold_file = $orig_scaffold_file;
insert_pmarks($scaffold_file);

##### Subroutines #####

# Subroutine to check the supplied paramaters are correct
sub check_parameters($) {
        my $options = shift;

## make sure output directory, scaffold file and strain name were passed
        unless ($options{output_dir} && $options{scaffold_input} && $options{strain}) {
		$logger->logdie("All the manadatory parameters should be passed");
	}

## make sure the output directory exists
        if (! -e "$options{output_dir}") {
		$logger->logdie("The $options{output_dir} output directory passed could not be read or does not exist");
       	}

## make sure the scaffold file exists and is readable
	if (! -e "$options{scaffold_input}") {
		$logger->logdie("The $options{scaffold_input} scaffold input file passed could not be read or does not exist");
	} else {
		my @ctg_ip = &read_file($options{'scaffold_input'});
		foreach my $content (@ctg_ip) {
			chomp($content);
			next if ($content =~ /^\s*$/);	#ignore whitespace
			next if ($content =~ /^#/);	#ignore 
			if($content =~ /^>/) {
## If the scaffold_input is a multi fasta file containing scaffolds				
				$orig_scaffold_file = $options{'scaffold_input'};
				last;
			} elsif ($content =~ /\//g) {
## If the scaffold_input is a list file containing paths to scaffold files 
				concat_files();
				last;
			} else {
## Else the scaffold_input does not contain correct data - neither fasta sequence nor list of file paths
				$logger->logdie("The $options{scaffold_input} file is neither a fasta file nor a list file of paths. Incorrect input");
			}
		}
## Assign default linker sequence if not supplied	
	$options{linker_sequence} = 'NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN' unless ( defined $options{linker_sequence});
	}
}

# Retrieves the sequence information from each path within the list file
# Next concatenates the FASTA sequences into a single multi.fasta file
sub concat_files() {
	$orig_scaffold_file = $options{output_dir}."/".$options{strain}.".multi.fasta";
	my @scaffold_files = read_file($options{'scaffold_input'});
	foreach my $path (@scaffold_files) {
		chomp($path);
		next if ($path =~ /^\s*$/);
		next if ($path =~ /^#/);
		if(-e $path) {
			system("cat $path >> $orig_scaffold_file");			
		}
	}
	system("$Bin/clean_fasta $orig_scaffold_file");	
}

# Replaces the string of N's from a scaffold sequencing gap with PMARKS
sub insert_pmarks($) {
	my $scaffold_file = shift;
	my $new_scaffold_file = $options{output_dir}."/".$options{strain}.".pmark.multi.fasta";
	my @fasta = read_file($scaffold_file);
	open (NEW, ">$new_scaffold_file"), or $logger->logdie("Could not open $new_scaffold_file file for writing: $!\n");
	my ($header, $sequence);
	my $new_sequence = "";
	my $n_counter = 0;
	my $n_start = -1;
	my $n_end = -1;
	my $append_seq = 0;

	foreach my $line (@fasta) {
		chomp $line;
		if ($line =~ /^>/) {
## 1st header/sequence combo is not defined yet
			if (defined $header) {
## Go through scaffold sequence and check for N's
				foreach my $nuc (0..(length($sequence) - 1)) {
					if (uc(substr($sequence, $nuc, 1)) eq 'N') {
						$n_counter++;
						$n_start = $nuc if ($n_counter == 1);
					} else {
						$n_counter = 0;
## If the start site has been established, declare the end site where the first non-N occurs
						if ($n_start > -1) {
							$n_end = $nuc;
## Contruct and add to new sequence with PMARKS based on the number of consecutive N's
							if (($n_end - $n_start) <= $length_cutoff) {
								if (length($new_sequence) == 0) {
									$new_sequence = substr($sequence, 0, $n_start);
								} else {
									$new_sequence .= substr($sequence, $append_seq, $n_start-$append_seq);
								}
								$new_sequence .= $options{linker_sequence};
							} else {
								if (length($new_sequence) == 0) {
									$new_sequence = substr($sequence, 0, $n_start);
								} else {
									$new_sequence .= substr($sequence, $append_seq, $n_start-$append_seq);
								}
								$new_sequence .= $options{linker_sequence};
								for (my $i = 0; $i < ($n_end - $n_start); $i++) {
									$new_sequence .= 'N';
								}
								$new_sequence .= $options{linker_sequence};
							}
## Reset the start site to mark the site of new N-strings later in the same scaffold
## No need to reset the end site as the start site dictates when it is updated
## The final instance of the end site will help to finish appending the sequence
							$n_start = -1;
							$append_seq = $n_end;
						}
					}
				}
				$new_sequence .= substr($sequence, $append_seq);
## prints previously saved header/new_sequence combo
				print NEW ">".$header, "\n";
				print NEW $new_sequence, "\n";
			}
## Save header of next sequence and initialize next sequence, counter, start and end sites
			$header = substr($line,1);
			$sequence = "";	
			$new_sequence = "";
			$n_counter = 0;
			$n_start = -1;
			$n_end = 0;
			$append_seq = 0;
		} else {
			next if ($line =~ /^\s*$/);	#ignore whitespace
			$sequence .= uc($line);
		}
	}	

## Print last set of header and new sequence
	foreach my $nuc (0..(length($sequence) - 1)) {
		if (uc(substr($sequence, $nuc, 1)) eq 'N') {
			$n_counter++;
			$n_start = $nuc if ($n_counter == 1);
		} else {
			$n_counter = 0;
			if ($n_start > -1) {
				$n_end = $nuc;
				if (($n_end - $n_start) <= $length_cutoff) {
					if (length($new_sequence) == 0) {
						$new_sequence = substr($sequence, 0, $n_start);
					} else {						
						$new_sequence .= substr($sequence, $append_seq, $n_start-$append_seq);
					}
					$new_sequence .= $options{linker_sequence};
				} else {
					if (length($new_sequence) == 0) {
						$new_sequence = substr($sequence, 0, $n_start);
					} else {
						$new_sequence .= substr($sequence, $append_seq, $n_start-$append_seq);
					}
					$new_sequence .= $options{linker_sequence};
					for (my $i = 0; $i < ($n_end - $n_start); $i++) {
						$new_sequence .= 'N';
					}
					$new_sequence .= $options{linker_sequence};
				}
				$n_start = -1;
				$append_seq = $n_end;
			}
		}
	}
	$new_sequence .= substr($sequence,- $append_seq);
	print NEW ">".$header, "\n";
	print NEW $new_sequence, "\n";
	close NEW;
}

# Subroutine to read files
sub read_file($) {
	my $filename = shift;
	my @lines;
	open(FH , "<$filename")  || $logger->logdie("Could not open $filename file for reading.$!");
	@lines = <FH>;
	close(FH);
	return(@lines);
}
