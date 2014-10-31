#!/usr/bin/env perl

=head1 NAME

pseudomolecule2glimmer3.pl - Get coordinates from genes that were put into a pseudomolecule in glimmer3 form

=head1 SYNOPSIS

 USAGE: perl_template.pl
       --input_file=/path/to/some/input.file
       --output=/path/to/transterm.file
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--pmarks_bsml_file,-b>

B<--input_fsa_file, -f>

B<--output_file,-o>

B<--contig_name, -c>
	Name of contig to be used.  Will append count numbers at the end of the name

B<--linker_seq, -l>
	Pmarks spacer sequence
	
B<-pseudomolecule_file, -p>

B<--log,-l>
    Logfile.

B<--debug,-d>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

 DESCRIPTION
 
=head1  INPUT
	A pmarks2bsml BSML output file.  This file will output coordinates for linker pmark sequences placed in a pseudomolecule.  This allows for us to precisely determine where the start and stop coords for each gene is.
	This file can be either a list file or an individual fasta file
	
	Alternatively just read in the multifasta file for our genes, and compare to the pseudomolecule.  This list should only have 1 pseudomocule file, since this script is intended for a pseudomolecule created by genecalls.

=head1 OUTPUT
	The glimmer3 tab-delimited file should have four columns (seqid, start, end, strand) as shown below, where "start" is always less than "end".  There is also a ">" in front of the seqid.  
	>GNMG2136.assembly.1	7	498	+
	
	This output can be used in genecalls2bsml.pl to create bsml to input to the translate_sequence.translate_prediction, bsml2fsa.prediction_CDS, and promote_gene_prediction.promote_prediction components in the Prokaryotic Annotation Pipeline

=head1 FUTURE CHANGES
	1. Add option to consider linker sequence
	2. Add ability to choose b/t pmarks2bsml output or straight read multifasta file 

=head1  CONTACT
    Shaun Adkins
    sadkins@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
my $contig_name = 'orf';
####################################################

my %options;

my $results = GetOptions (\%options,
                         "pmarks_bsml_file|b=s",
                         "input_fsa_file|f=s",
                         "pseudomolecule_file|p=s",
                         "linker_seq|l=s",
                         "output_file|o=s",
                         "contig_name|c=s",
                         "log|l=s",
                         "debug|d=s",
                         "help|h"
                          );

&check_options(\%options);

open my $out, '>'. $options{'output_file'} or die "Cannot open output file for writing: $!\n";

if (defined $options{'pmarks_bsml_file'}) {
	&read_bsml($options{'pmarks_bsml_file'}, $options{'contig_name'}, $out);
	
} else {
	&read_fasta($options{'input_fsa_file'});
}

close $out;
exit(0);

sub read_bsml {
	my $bsml_file = shift;
	my $header = shift;
	my $outfh = shift;
	
	my $start = 1;
	my $end = "";
	my $count = 0;
	
	open BSML, $bsml_file or die "Cannot open $bsml_file for reading: $!\n";
	while (<BSML>) {
		my $line = $_;
		chomp $line;
		
		if ($line =~ /class=\"assembly\" id=\"([^\"]+)\" molecule=/) {
			print $outfh ">$1\n";
		} 
		# Grab the start and end pmark linker coords, and use to calculate gene coords
		if ($line =~ /startpos=\"(\d+)\" endpos=\"(\d+)\">/ ) {
			$end = $1;
			next if ($end == 0);	#Pmark spacer at beginning of file.
			write_tab($outfh, $start, $end, $header, $count++);
			$start = $2 + 1;
		}
	}
	
	close BSML;
	return;
}

sub read_fasta {
	my $fasta_file = shift;
	print "Currently not supporting this command... check back later\n";
	return;
}

sub write_tab {
	my ($outfh, $s, $e, $h, $c) = @_;
	#print "$h$c\t$s\t$e\t+1\n";
	print $outfh "$h$c\t$s\t$e\t+1\n";
	return;
}

## Subroutine to read files
sub read_file {
	my $filename = shift;
	my @lines;
	open(FH , "< $filename")  || &_log($ERROR, "Could not open $filename file for reading.$!");
	@lines = <FH>;
	close(FH);
	return(@lines);
} 

sub get_bsml {
	my @bsml_files = &read_file($options{'pmarks_bsml_file'});
	&_log($ERROR, "Currently only supporting a list file of only 1 BSML file path...check back later\n") if (scalar @bsml_files > 1);
	$options{'pmarks_bsml_file'} = shift(@bsml_files);
}

sub check_options {
   my $opts = shift;
   if( $opts->{'help'} ) {
       &_pod;
   }

   if( $opts->{'log'} ) {
       open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
   }

   $debug = $opts->{'debug'} if( $opts->{'debug'} );

   foreach my $req ( qw(pmarks_bsml_file output_file) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }
   
   	my @ctg_ip = &read_file($options{'pmarks_bsml_file'});
	foreach my $content (@ctg_ip) {
		chomp($content);
		next if ($content =~ /^\s*$/);
		next if ($content =~ /^#/);
		if($content =~ /^</) {
## If the pmarks_bsml is a bsml file... just exit
			last;
		} elsif ($content =~ /\//g) {
## If the pmarks_bsml is a list file containing paths to bsml files 
			&get_bsml();
			last;
		} else {
## Else the plmarks_bsml inputdoes not contain correct data - neither bsml sequence nor list of file paths
			&_log($ERROR, "The $options{contig_input} file is neither a fasta file nor a list file of paths. Incorrect input");
		}
	}
   
 	$options{'contig_name'} = $contig_name unless ( defined $options{'contig_name'});  
}

sub _log {
   my ($level, $msg) = @_;
   if( $level <= $debug ) {
      print STDOUT "$msg\n";
   }
   print $logfh "$msg\n" if( defined( $logfh ) );
   exit(1) if( $level == $ERROR );
}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
