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
    Either a comma-separated or tab-separated file containing the following columns:
    1) Name of Database
    2) Locus_tag ID prefix
    3) Path to a curate_common_names rules file
    4+) Any other DB related information (to come later)

=head1 OUTPUT
	A list file containing paths to individual comma-separated files.
	Each file will be in the form of:
	dbname,id_prefix,/path/to/rules.txt,Other_db_info
	
	
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
	my $list_file = $output_dir . "db.list";
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
