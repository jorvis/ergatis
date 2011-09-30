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
use FindBin qw($Bin);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Ergatis::Logger;

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
	my $new_scaffold_file = $options{output_dir}."/".$options{strain}.".fasta";
	my $coord_file = $options{output_dir}."/".$options{strain}.".fasta.pmarks";
	my @fasta = read_file($scaffold_file);
	open (NEW, ">$new_scaffold_file"), or $logger->logdie("Could not open $new_scaffold_file file for writing: $!\n");
	open (COORD, ">$coord_file"), or $logger->logdie("Could not open $coord_file file for writing: $!\n");
        my($header, $sequence) = (undef, '');
       
## Subroutine to insert pmarks where any number of N's are present in the sequence
        my $ins_pmarks = sub {
	    my $old_seq_length = length($sequence);
	    my ($pmark_added, $n_removed) = (0,0);
	    my $pmark_len = length($options{'linker_sequence'});
            
	    return if (!defined($header));
	    my @n_array = ($sequence =~ /(n+)/gi);
	    foreach my $string (@n_array) {
		if (length($string) <= $length_cutoff) {
		    $n_removed += length($string);
		    $pmark_added ++;    
		} else {
		    $pmark_added += 2;
		}
	    }
	    my $new_seq_length = $old_seq_length + ($pmark_added * $pmark_len) - $n_removed;

	    $sequence =~ s/(n+)/(length($1) <= $length_cutoff) ? $options{'linker_sequence'} : $options{'linker_sequence'} . $1 . $options{'linker_sequence'}/gie;  
	    my $actual_seq_length = length($sequence);

            print NEW ">".$header, "\n";
	    print NEW "$1\n" while( $sequence =~ /(\w{1,60})/g );
	    #print NEW $sequence, "\n";	 

	    print COORD ">".$header, "\n";   
   	    while ($sequence =~ /$options{'linker_sequence'}/gi) {
		my $l_end = pos($sequence);
		my $l_start = $l_end - length($options{'linker_sequence'});
		print COORD $l_start, "\t", $l_end,"\n"
	    }
	    if ($new_seq_length != $actual_seq_length) {
	    	$logger->logdie("Projected sequence length ($new_seq_length) does not match actual sequence length ($actual_seq_length).  Please report this bug.");
	    }	
        };

## Go through line by line and run the ins_pmarks subroutine
        foreach my $line (@fasta) {
            chomp $line;
            if ($line =~ /^>/) {
                my $next_header = substr($line, 1);
                &$ins_pmarks();
                $header = $next_header;
                $sequence = '';
            } else {
                next if ($line =~ /^\s*$/);  #ignore whitespace
                $sequence .= uc($line);
            }
        }
## Run ins_pmarks on the final sequence
        &$ins_pmarks();
        close NEW;
	close COORD;
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
