#!/usr/bin/env perl

=head1 NAME

get_pmark_locations.pl - program to get pmark coordinates

=head1 SYNOPSIS

    create_pmark_locations.pl --output_file <outfile.pmarks>  --input_file <fasta file> [ --linker_seq <NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN> --help ]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1  CONTACT
    Shaun Adkins
    sadkins@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

#################################################
# MAIN PROGRAM					#
#################################################
my %options;
my $results = GetOptions (\%options,
                'input_file|i=s',
		        'output_file|o=s',
		        'linker_seq|k=s',
                'help|h') || pod2usage();

## Display documentation
if( $options{'help'} ){
        pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## Make sure everything passed was correct
check_parameters(\%options);
determine_pmark_coords($options{'input_file'}, $options{'linker_seq'}, $options{'out_file'});

# Subroutine to check the supplied paramaters are correct
sub check_parameters {
    my $options = shift;

    ## make sure output directory, contig file and strain name were passed
     unless ($options{output_file} && $options{input_file}) {
		die("All the manadatory parameters should be passed");
	}
	$options{'linker_sequence'} = 'NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN' unless (defined $options{linker_sequence});
}

## Subroutine to read files
sub read_file {
	my $filename = shift;
	my @lines;
	open(FH , "< $filename")  || die("Could not open $filename file for reading.$!");
	@lines = <FH>;
	close(FH);
	return(@lines);
}

sub determine_pmark_coords {
	my ($in_file,$linker,$pmarks_file) = @_;
	my $seq = "";
	my $header;
	my $flag = 1;
	my $linker_start = 0;
	my $linker_len = length($linker);
	my $linker_end = $linker_len;
	my %contig_hash;
	open(POUT, "> $pmarks_file")  or die("Could not open $pmarks_file file for writing\n");
	
    ## Read input contig file
	my @contigs = &read_file($in_file);
	foreach my $line (@contigs) {
		chomp($line);
		if($line =~ /^>(.*)/) {
			unless($flag) {
				chomp($seq);
				$contig_hash{$header} = $seq;
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
	$contig_hash{$header} = $seq if(defined($header));

    ## Header for the pseudomolecule	
    foreach my $head (keys %contig_hash) {
        my $offset = 0;
        my $result = index($contig_hash{$head}, $linker, $offset);
        while ($result != -1) {
	        print POUT ">$head\n";
            $linker_start = $result;
            $linker_end = $linker_start + $linker_len;
            print POUT "$linker_start\t$linker_end\n";
            $offset = $result + 1;
            $result = index($contig_hash{$head}, $linker, $offset);
        }
    }
	close(POUT);
}