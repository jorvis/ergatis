#!/usr/bin/env perl
package main;    # Change this to reflect name of script.  If you want to use as module, make sure to save as .pm instead of .pl

=head1 NAME

map_header_to_transcript_in_tab.pl - Description

=head1 SYNOPSIS

 USAGE: map_header_to_transcript_in_tab.pl
 	   --tabfile=/path/to/annotated.tab
       --feature_relationship_file=/path/to/some/bsml2featurerelaionship.mapping.txt
       --fasta_list=/path/to/multifasta.list
       --output_file=/path/to/transterm.file
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--tabfile,-t> Annotation tabfile.

B<--feature_relations_file, -r> A bsml2featurerelationship mapping file

B<--pseudomolecule_list, -p> List file pointing to a pseudomolecule fasta file.  The filename will be used to grab the contig order file that was also outputted from the create_pseudomolecule script.  Needs to only be 1 file in list

B<--output_file, -o> Will be new tabfile with headers instead of transcript IDs at the start

B<--print_mapping, -m> Print a mapping file showing the header ID, the transcript ID, and the transcript coords

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
my $pseudofile;
my $order_file;
my $outfile;

my $mapping_path = undef;

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
                         "pseudomolecule_list|p=s",
                         "output_file|o=s",
                         "print_mapping|m=s",
                         "log|l=s",
                         "debug|d=s",
                         "help|h"
                          );

    &check_options(\%options);

    $tabfile = $options{'tabfile'};
    $frfile = $options{'feature_relationship_file'};

    my @pseudo_files= &read_file($options{'pseudomolecule_list'});
	&_log($ERROR, "Currently only supporting a list file of only 1 pseudomolecule file path...check back later\n") if (scalar @pseudo_files > 1);
	$pseudofile = shift(@pseudo_files);
	chomp $pseudofile;	# Wasn't chomped when pushed into array
	$order_file = $pseudofile . ".order";
    $outfile = $options{'output_file'};

    $headers = get_ordered_headers($order_file);
    $transcripts = get_ordered_transcripts($frfile);
    #print "scalar headers: " . scalar(@$headers) . "\n";
    #print "scalar transcripts: " . scalar(@$transcripts) . "\n";
    &_log($ERROR, "Incorrect number of headers to transcript IDs\n". scalar @$headers. " headers to ". scalar @$transcripts. " transcipts\n") if (scalar @$headers != scalar @$transcripts);
    sub_in_headers_in_tabfile($tabfile, $outfile, $headers, $transcripts);
    print_mapped_header_transcript($headers, $transcripts, $frfile, $mapping_path) if (defined $mapping_path);
    exit(0);
}

# Retrieve the list of headers sorted by appearance in pseudomolecule
sub get_ordered_headers {
	my $order = shift;
	my @headers;
	open ORDER, $order or die "Cannot open pseudomolecule ordering file $order for reading: $!\n";
	while (my $line = <ORDER>) {
		chomp $line;
		my @split = split(/\t/, $line);
		push @headers, $split[0];
	}
	close ORDER;
	return \@headers;
}

# Get list of transcript IDs sorted by coordinate position
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

# Replace transcript ID with mapped header
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

# Search for the header in its array by determining the index position of the current transcript
sub find_header {
	my ($t, $headers, $transcripts) = @_;
	#Grep for the chopped transcript ID, and get index number
	my ($index) = grep {$$transcripts[$_] =~ /$t/} 0..(scalar @$transcripts - 1);
	return $$headers[$index];
}

# Print a tab-delimited file showing genecalls header, transcript ID and coordinates.
# This subroutine is rough around the edges since I just need it for a specific debugging purpose
sub print_mapped_header_transcript {
	my ($headers, $transcripts, $frfile, $mapfile) = @_;
	my %data;
	open FR, $frfile or die "Cannot open $frfile for reading: $!\n";
	while (my $line = <FR>) {
		chomp $line;
		#polypeptide CDS transcript pseudomolecule start:end
		if ($line =~ /^\S+\s+\S+\s+(\S+)\s+\S+\s+\d+:\d+\/0$/) {
			my $transcript = $1;
			$data{$transcript} = $line;
		}
	}
	close FR;
	open MAP, ">".$mapfile or die "Cannot open $mapfile for writing: $!\n";
	for (my $i = 0; $i <= $#{$headers}; $i++){
		print MAP $$headers[$i], "\t", $data{$$transcripts[$i]}, "\n";
	}
	close MAP;
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

# Make sure options are handled correctly
sub check_options {
   my $opts = shift;
   if( $opts->{'help'} ) {
       &_pod;
   }

   if( $opts->{'log'} ) {
       open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
   }

   $debug = $opts->{'debug'} if( $opts->{'debug'} );

   # If feature_relationhip file is not present, then we assume mapping does not need to take place
   if (! $opts->{'feature_relationship_file'}){
		print STDERR "No feature_relationship file found so skipping rest of mapping program.\n";
		exit(0);
	}


   foreach my $req ( qw(tabfile feature_relationship_file pseudomolecule_list output_file) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }

   $mapping_path = $opts->{'print_mapping'} if ($opts->{'print_mapping'});
}

# Private method for printing logging messages
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
