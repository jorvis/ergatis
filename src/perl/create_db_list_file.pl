#!/usr/bin/env perl

=head1 NAME

create_db_list_file.pl - Converts a csv or tab delimited file with database properties into a list
	of individual file paths for each individual entry

=head1 SYNOPSIS

 USAGE: perl_template.pl
       --input_file=/path/to/some/db_info.tab
       --output_dir=/path/to/output_dir
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>

B<--output_dir,-o>

B<--log,-l>
    Logfile.

B<--debug,-d>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

 DESCRIPTION
 
=head1  INPUT
    A tab-separated file containing the following columns:
 1. Name of the manatee database used while running the prokaryotic annotation pipeline
 2. NCBI_locus_tag prefixto be assigned to the genes.
 3. Absolute path to the tab-delimited rules file for fixing gene product names in the database.
 4. Absolute path to the tab-delimited gene symbol rules file for fixing gene symbols in the database.
 5. Bioproject Id beginning with PRJ [for .sbt file creation]
 6. Organism name [for tbl2asn]
 7. Strain name [for tbl2asn]
 8. Organism serotype [for tbl2asn]
 9. Name of the host from which the bacterial sample was isolated [for tbl2asn]
 10. Date on which bacterial sample was collected (mostly year only) [for tbl2asn]
 11. Country from which the bacterial sample was isolated [for tbl2asn]
 12. Assembly method used to assemble the genome sequence [for .cmt file creation]
 13. Genome Coverage [for .cmt file creation]
 14. Sequencing platform [for .cmt file creation]
 15. Contact person's name for Genbank submission. Format : Last name\sFirst name [for .sbt file creation]
 16. Contact person's email address [for .sbt file creation]
 17. Comma-separated author list for Genbank submission. Format : Each name should be Last name\sFirst name\sMiddle Initial [for .sbt file creation]
 18. Title of the publication [for .sbt file creation]
 19. Absolute path to list of deprecated or bad EC numbers, one per row (used in update_ec_numbers)
 20. Isolation source of the bacterial sample.
 
 A comma-separated file can be used but it is not recommended unless you plan on entering only the first 4 columsn or less of data

=head1 OUTPUT
	A list file containing paths to individual comma-separated or tab-delimited files.
	Each file will be in the form of:
	dbname	id_prefix	/path/to/rules.txt	/path/to/gene_syms.txt	Other_db_info	...
	
	
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

my $input_file;
my $output_dir;
my %seen_db;
my @files;
####################################################

my %options;
my $results = GetOptions (\%options,
                         "input_file|i=s",
                         "output_dir|o=s",
                         "log|l=s",
                         "debug|d=s",
                         "help|h"
                          );

&check_options(\%options);

parse_delimited_file($input_file);
write_list_file(\@files);
exit(0);

# Parse the input file and write each row to individual files
sub parse_delimited_file {
    my $in = shift;
    
    open IN, $in or die "Cannot open input file $in for reading: $!\n";
    while (<IN>) {
        my $line = $_;
        chomp $line;
        next if $line =~ /^\s*$/;	# Skip blank lines
        my ($db, $rest) = split (/,|\t/, $line, 2);
       # print $db, "\n";
        if ($seen_db{$db}++) {
        	&_log($WARN, "$db has already been encountered in this file.  Skipping this instance...");
        	next;
        }
        
        # Write line to output file
        my $output_file = $output_dir . "/$db.txt";
        #print $output_file, "\n";
        open OUT, ">$output_file" or die "Cannot open output_file $output_file for writing: $!\n";
        print OUT $line, "\n";
        close OUT;
        push @files, $output_file;
    }
    
    close IN;
    return
}

# Create a list file from each of the individual files
sub write_list_file {
	my $outputs = shift;
	my $list_file = $output_dir . "/db.list";
	open LIST, ">$list_file" or die "Cannot open list file $list_file for writing: $!\n";
	print LIST $_ . "\n" foreach @$outputs;
	close LIST;
	return;
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

   foreach my $req ( qw(input_file output_dir) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }
   
   $input_file = $opts->{'input_file'};
   $output_dir = $opts->{'output_dir'};
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
