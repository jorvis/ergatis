#!/usr/bin/env perl

=head1 NAME

RAPSearch2btab.pl - Convert RAPSearch .m8 to BLAST .btab format

=head1 SYNOPSIS

 USAGE: RAPSearch2btab.pl
       --input_file=/path/to/some/results.m8
       --output=/path/to/results.btab
	   --uniprot_db=/path/to/uniprot.sqlite
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>

B<--output_file,-o>

B<--uniprot_db,-u>
	This is the path to the UniProt SQLite database.
	Note this is not the pre-RAPSearch Uniref database

B<--log,-l>
    Logfile.

B<--debug,-d>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

 DESCRIPTION
 
=head1  INPUT

    Describe the input

=head1 OUTPUT

The output is defined using the --output option.  The file created is tab-delimited and
composed of the following columns.

    1   query_name
    2   date
    3   query_length
    4   algorithm
    5   database_name
    6   hit_name
    7   qry_start
    8   qry_end
    9   hit_start
    10  hit_end
    11  percent_identity
    12  percent_similarity
    13  raw_score
    14  bit_score
    15  NULL
    16  hit_description
    17  blast_frame
    18  qry_strand (Plus | Minus)
    19  hit_length
    20  e_value
    21  p_value

=head1	CONTACT

    Shaun Adkins
    sadkins@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use bigint;
use DBI;
use Bio::DB::Fasta;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
####################################################

my %options;
my $ifh;
my $currofh;

my $results = GetOptions (\%options,
                         "input_file|i=s",
                         "output_file|o=s",
						 "uniprot_db|u=s",
                         "log|l=s",
                         "debug|d=s",
                         "help|h"
                          );

&check_options(\%options);
# Connect to the SQLite db
my $dbh = &_connect($options{'uniprot_db'});

## get a filehandle on the input
if ($options{'input_file'} =~ /\.(gz|gzip)$/) {
    open ($ifh, "<:gzip", $options{'input'})
      || &_log($ERROR, "can't open input file:\n$!");
} else {
    open($ifh, "<$options{input_file}") || &_log($ERROR, "can't read the input file: $!");
}

open ($currofh, ">$options{output_file}") || &_log($ERROR, "can't create output file for BLAST report: $!");

read_input($ifh, $currofh);
close $ifh;
close $currofh;
exit(0);

sub read_input {
	my $ifh = shift;
	my $ofh = shift;

	my $query_file;
	my $subject_file;
	my $btab_arr;
	my $bio_db;

	my $line = <$ifh>;	# RAPSearch line
	$line = <$ifh>;		# Job submission line
	$line = <$ifh>;		# Query line
	chomp $line;
	if ($line =~ /^#\sQuery\s:\s+(.+)/){
		$query_file = $1;
		# Need to do this since Bio::DB fails if blank line is present
		remove_blank_lines($query_file);
		# Index seqs for utility reasons
		$bio_db = Bio::DB::Fasta->new($query_file);
	}
	$line = <$ifh>;		# Subject line
	chomp $line;
	if ($line =~ /^#\sSubject\s:\s+(.+)/){
		$subject_file = $1;
	}
	$line = <$ifh>;		# Tabbed headers line

	# Now we get to the "real" parsing
	while ($line = <$ifh>){
		chomp $line;
		$btab_arr = parse_line($line);
		# More btab additions
		$btab_arr->[1] = '';	# Never actually seen the date added in past btabs
		$btab_arr->[2] = get_query_len($btab_arr->[0], $bio_db);
		$btab_arr->[3] = "RAPSearch2";
		$btab_arr->[4] = $subject_file;
		$btab_arr->[11] = 0;	# setting %sim = 0 for now
		$btab_arr->[12] = 0;	# setting raw = 0 for now
		$btab_arr->[14] = '';	# NULL field
		$btab_arr->[16] = 1;	# Strand and frame are only present in XML output
		$btab_arr->[17] = 1;
		$btab_arr->[18] = get_hit_len($btab_arr->[5]);	
		$btab_arr->[20] = calculate_pvalue($btab_arr->[19]);
		#&_log($DEBUG, join("\t", @$btab_arr));
		print $ofh join("\t", @$btab_arr), "\n";
	}
	return;
}

# Parse m8 line to get basic btab information
sub parse_line {
	my $line = shift;
	my @btab_arr;
	$#btab_arr = 20;	# Pre-sizing the btab array
	my ($query, $subj, $identity, $aln_len, $mismatch, $gap_openings, $q_start, $q_end, $s_start, $s_end, $e_val, $bit) = split(/\t/, $line);

	$btab_arr[0] = $query;
	my ($hit, $desc);	#For whatever reason, 'split (/\s+/, $subj, 2) ignored the limit part
	if ($subj =~ /^(UniRef100_\w+)\s+(.+)/){
		$hit = $1;
		$desc = $2;
	}
	$btab_arr[5] = $hit;
	$btab_arr[6] = $q_start;
	$btab_arr[7] = $q_end;
	$btab_arr[8] = $s_start;
	$btab_arr[9] = $s_end;
	$btab_arr[10] = $identity;
	$btab_arr[13] = $bit;
	$btab_arr[15] = $desc;
	$btab_arr[19] = $e_val;

	return \@btab_arr;
}
# Remove blank lines in fasta file.  Typically at the end of file.
sub remove_blank_lines {
	my $fasta = shift;
	my $cmd = "perl -i.orig -lne 'print if length' $fasta";
	system($cmd);
	return;
}

# Get query length of given query from query file
sub get_query_len {
	my ($q, $db) = @_;
	return $db->length($q);
}

# Get Uniref hit length from SQLite database
sub get_hit_len {
	my $hit = shift;
	my $acc;
	if ($hit =~/UniRef100_(\w+)/){
		$acc = $1;
	}	
	
	my $sth = $dbh->prepare("SELECT res_length FROM uniref_acc WHERE accession = ?");
	$sth->execute($acc);
	
	# Should only have 1 row and 1 value;
	my @row = $sth->fetchrow_array();
	my $hit_len = $row[0];
	
	&_log($WARN, "Uniref Accession $acc did not have a residue length in the SQLite3 table.  Ensure the accession is correct or that the entry has a defined 'res_length' field") unless (defined $hit_len);

	# Joshua Orvis created the SQLite table and said some fields may be null,
	# so adding this dummy value in case he needs to rebuild the database
	# Just comment out the log message above if this occurs or switch between WARN and ERROR.
	if (! defined ($hit_len) ){
		$hit_len = 1;
	}

	return $hit_len;
}

#See http://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html
sub calculate_pvalue {
    my $evalue = shift;

    my $estimate = 0.57721566490153;

    #my $p = 1 - (bexp( (-1*$evalue), 4 ) );
	#( 1 - ( $estimate**(-1*$evalue) ) );
    return $evalue;

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

   foreach my $req ( qw(input_file output_file uniprot_db) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }
}

sub _connect {
    my ($db) = @_; 

    my $dbh;
    eval {
        $dbh = DBI->connect("DBI:SQLite:$db", "", "",
                {
                    'RaiseError' => 1,
                }   );  
    };  
    if ($@) {
        &_log($ERROR, "Could not connect to database " . $@);
    }   
    return $dbh;
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
