#!/usr/local/bin/perl

=head1  NAME 

predotar2legacydb.pl - load predotar data into the legacy database schema.

=head1 SYNOPSIS

USAGE: predotar2legacydb.pl 
        --input_list=/path/to/somefile.raw.list
        --database=sma1
      [ --log=/path/to/some.log ]

=head1 OPTIONS

B<--input_list,-i> 
    raw list file from a predotar workflow component run.

B<--database,-d> 
    Sybase project database ID.

B<--log,-l> 
    optional.  will create a log file with summaries of all actions performed.

B<--help,-h> 
    This help message/documentation.

=head1   DESCRIPTION

This script is used to load predotar evidence into the legacy database schema.
Only raw output files from running predotar on a single input sequence are supported.
Both plant and animal/yeast output formats should be supported.

=head1 INPUT

The input should be the '*.predotar.raw.list' list file from a predotar workflow component run,
or a similarly formatted list of predotar raw output files with filenames matching the following
regex:
		/(\d+)_(\d+)\.predotar\.raw/
that allows translation to a corresponding feat_name of:
		$1.m$2

=head1 OUTPUT

This is a database loading script.  There is no other output unless you use 
the --log option.

=head1 CONTACT

    Brett Whitty
    bwhitty@tigr.org

=cut

use strict;
use warnings;
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
			  'input_list|i=s',
              'database|d=s',
              'log|l=s',
			  'help|h') || pod2usage();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

## make sure all passed options are peachy
&check_parameters(\%options);

## get the current user
my ($user) = getpwuid($<);

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">$options{log}") || die "can't create log file: $!";
}

## create the database connection
my $dbh = DBI->connect("dbi:Sybase:server=SYBTIGR; packetSize=8092", 'egc', 'egcpwd', {RaiseError => 1, PrintError => 1});
my $result = $dbh->do("use $options{database}");

## prepare a statement for inserting into ORF_attribute
my $qry = "INSERT INTO ORF_attribute (feat_name, att_type, curated, method, date, assignby, " .
       "score, score_desc, score2, score2_desc, score3, score3_desc, score4, score4_desc, score5, score5_desc) " .
       "VALUES ( ?, 'predotar', 0, 'predotar', GETDATE(), '$user', ?, 'prediction', ?, 'Mit', ?, 'Plast', ?, 'ER', ?, 'None' )";
my $orf_attribute_insert = $dbh->prepare($qry);

## prepare a statement for deleting old predotar 
$qry = 'DELETE FROM ORF_attribute ' .
       "WHERE att_type in ('predotar') " .
       ' and feat_name = ?';
my $delete = $dbh->prepare($qry);


## load the file list
my @files;
open(my $listfh, "<$options{input_list}") || die "can't open file list: $!";
while (<$listfh>) {
	chomp;
    next if (/^\s*$/);
    _log("adding $_ to the process list");
    push @files, $_;
}
close $listfh;


foreach my $file (@files) {

    _log("processing file $file");

    ## extract the feat_name from the file name:
    ##  bma1.model.14990_07908.iprscan.raw -> 14990.m07908
    my $feat_name;
    if ( $file =~ /(\d+)_(\d+)\.predotar\.raw/ ) {
        $feat_name = "$1.m$2";
    } else {
        die "improperly named input file.  couldn't extract feat name from $file";
    }

   
	## parse the input file
	my $predotar_hash_ref = parse_predotar_raw_file($file);
	
	my $key_count = scalar(keys(%{$predotar_hash_ref}));

	## if the input sequence was discarded we will have # of hash keys = 0
	## otherwise there will be 4 or 5 depending on whether input setting was plant or animal/yeast
	if ($key_count < 4) {
		_log("no data added for $feat_name");
	} else {
		## delete the existing predotar evidence for this feat_name
	    _log("deleting all predotar evidence for $feat_name");
		$delete->execute($feat_name);
		## store the data in the db
	    _log("orf_attribute_insert->execute($feat_name," . $predotar_hash_ref->{'Prediction'}.
										 ",".$predotar_hash_ref->{'Mit'}.",".$predotar_hash_ref->{'Plast'}.
										 ",".$predotar_hash_ref->{'ER'}.",".$predotar_hash_ref->{'None'}.");");	
  	  $orf_attribute_insert->execute( 
										$feat_name,
										$predotar_hash_ref->{'Prediction'},
										$predotar_hash_ref->{'Mit'},
										$predotar_hash_ref->{'Plast'},
										$predotar_hash_ref->{'ER'},
										$predotar_hash_ref->{'None'},
								  	);
	}
}
$delete->finish();
$orf_attribute_insert->finish();
$dbh->disconnect();

exit();

sub parse_predotar_raw_file {
	
	my $raw_file = shift @_;

	open (IN, $raw_file) || die("Couldn't open '$raw_file' for reading.");
	
	my %results;
	my @result_fields;
	my @field_ids;
	while (<IN>) {
    	chomp;
    	
		next if ( /^\s*$/ );

	    ##recognize header fields
    	if (/^Seq\t/) {
			
			##save the field ids for verifying data later
			@field_ids = split("\t");

			my $result_line = <IN>;
			chomp $result_line;
			if ($result_line =~ /Discarding/) {
				_log("Input sequence discarded by predotar in '$raw_file'");
				last;
			}
			
			my @result_fields = split("\t", $result_line);
			
			my $col = scalar @field_ids;
			for (my $i=0; $i < $col; $i++) {
				$results{$field_ids[$i]} = $result_fields[$i];
			}
			last;
		}
    }
	close IN;
	
	return \%results;
}

sub _log {
    my $msg = shift;
    
    print $logfh "$msg\n" if $logfh;
}

sub check_parameters {
    
    ## database and input_list are required
    unless ( defined $options{database} && $options{input_list} ) {
        print STDERR "database and input_list options are required\n\n";
        pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
    }
    
    ## make sure input list exists
    unless ( -e $options{input_list} ) {
        print STDERR "\n\ninput_list $options{input_list} does not exist\n\n";
        exit(1);
    }
}

