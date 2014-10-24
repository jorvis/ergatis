#!/usr/bin/env perl
package main;    # Change this to reflect name of script.  If you want to use as module, make sure to save as .pm instead of .pl

=head1 NAME

map_header_to_transcript_in_tab.pl - Description

=head1 SYNOPSIS

 USAGE: map_header_to_transcript_in_tab.pl
 	   --tabfile=/path/to/annotated.tab
       --feature_relationship_file=/path/to/some/bsml2featurerelaionship.mapping.txt
       --fasta_file=/path/to/multifasta.fa
       --output_file=/path/to/transterm.file
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--tabfile,-t> Annotation tabfile.

B<--feature_relations_file, -r> A bsml2featurerelationship mapping file

B<--fasta_file, -f>	Multifasta file used as initial input for the pipeline

B<--output_file,-o> Will be new tabfile with headers instead of transcript IDs at the start

B<--log,-l>
    Logfile.

B<--debug,-d>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

This program will take headers from the initial multifasta file of a pipeline and used them as substitutes for the transcript IDs in a given annotation tabfile.  In the gene calls pipeline, each sequence is processed in order of appearance in the multifasta file to create a pseudomolecule.  The coordinates for each sequence appear in the bsml2featurerelationships mapping file and since the order of the headers is known, it is easy to map the transcript ID from there.  To avoid any potential parsing issues in the tabfile, the header will be parsed up to the first space character (or newline if the header is space-free)
 
=head1  INPUT
    Describe the input

=head1 OUTPUT
    Describe the output

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

my $tabfile;
my $frfile;
my $fastafile;
my $outfile;

my $headers;	#array_ref
my $transcripts;	#array_ref
####################################################

my %options;

# Allow program to run as module for unit testing if necessary
main() unless caller();

sub main {
	my $results = GetOptions (\%options,
                         "tabfile|t=s",
                         "feature_relationship_file|r=s",
                         "fasta_file|f=s",
                         "output_file|o=s",
                         "log|l=s",
                         "debug|d=s",
                         "help|h"
                          );

    &check_options(\%options);
    
    $tabfile = $options{'tabfile'};
    $frfile = $options{'feature_relationship_file'};
    $fastafile = $options{'fasta_file'};
    $outfile = $options{'output_file'};
    
    $headers = get_ordered_headers($fastafile);
    $transcripts = get_ordered_transcripts($frfile);
    #print "scalar headers: " . scalar(@$headers) . "\n";
    #print "scalar transcripts: " . scalar(@$transcripts) . "\n";
    &_log($ERROR, "Incorrect number of headers to transcript IDs") if (scalar @$headers != scalar @$transcripts);
    sub_in_headers_in_tabfile($tabfile, $outfile, $headers, $transcripts);
    exit(0);
}

sub get_ordered_headers {
	my $fasta = shift;
	my @headers;
	open FASTA, $fasta or die "Cannot open fasta file $fasta for reading: $!\n";
	while (my $line = <FASTA>) {
		chomp $line;
		push @headers, $1 if ($line =~ /^>(\S+)/);
	}
	close FASTA;
	return \@headers;
}

sub get_ordered_transcripts {
	my ($frfile) = @_;
	my @transcripts;
	my %data;
	open FR, $frfile or die "Cannot open feature_relationship file $frfile for reading: $!\n";
	while (my $line = <FR>) {
		chomp $line;
		#polypeptide CDS transcript pseudomolecule start:end
		if ($line =~ /^\S+\s+\S+\s+(\S+)\s+\S+\s+(\d+):\d+\/0$/) {
			my $key = $2;
			my $transcript = $1;
			&_log($ERROR, "Multiple entries found for starting coordinate $key") if (exists $data{$key});
			$data{$key} = $transcript;	
		}
	}
	
	#get numerically sorted version of transcript IDs
	foreach my $k (sort {$a <=> $b} keys %data) {
		#print $k, "\t", $data{$k}, "\n";
		push @transcripts, $data{$k};
	}
	
	close FR;
	return \@transcripts;
}

sub sub_in_headers_in_tabfile {
	my ($tab, $out, $headers, $transcripts) = @_;
	open TAB, $tab or die "Cannot open tabfile $tab for reading: $!\n";
	open OUT, ">$out" or die "Cannot open outfile $out for writing: $!\n";
	while (my $line = <TAB>) {
		chomp $line;
		my $t = $1 if ($line =~ /^(\S+)\s/);
		#The last number in the ID may be different (example -> .1 and .2 for pre/post overlap analysis)		
		(my $t2 = $t) =~ s/\.\d+$//;
		my $h = find_header($t2, $headers, $transcripts);
		$line =~ s/$t/$h/;
		print OUT $line . "\n";
	}	
	close TAB;
	close OUT;
	return;
}

sub find_header {
	my ($t, $headers, $transcripts) = @_;
	#Grep for the chopped transcript ID, and get index number
	my ($index) = grep {$$transcripts[$_] =~ /$t/} 0..(scalar @$transcripts - 1);
	return $$headers[$index];
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

   foreach my $req ( qw(tabfile feature_relationship_file fasta_file output_file) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }
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
