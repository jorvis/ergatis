#!/usr/local/bin/perl

=head1  NAME 

tmhmm2legacydb.pl - load tmhmm data into the legacy database schema.

=head1 SYNOPSIS

USAGE: tmhmm2legacydb.pl 
        --input_list=/path/to/somefile.bsml.list
        --database=sma1 [--server SOME_SERVER]
      [ --log=/path/to/some.log ]

=head1 OPTIONS

B<--input_list,-i> 
    BSML list file from a tmhmm workflow component run.

B<--database,-d> 
    Sybase project database ID.

B<--server,-s>
	Server hosting database (optional, defaults to SYBTIGR).

B<--log,-l> 
    optional.  will create a log file with summaries of all actions performed.

B<--help,-h> 
    This help message/documentation.

=head1   DESCRIPTION

This script is used to load tmhmm evidence into the legacy database schema.
Only BSML output files from running tmhmm on a single input sequence are supported.

=head1 INPUT

The input should be the '*.tmhmm.bsml.list' list file from a tmhmm workflow component run,
or a similarly formatted list of tmhmm bsml output files with filenames matching the following
regex:
		/(\d+)_(\d+)\.tmhmm\.bsml/
that allows translation to a corresponding feat_name of:
		$1.m$2

=head1 OUTPUT

This is a database loading script.  There is no other output unless you use 
the --log option.

=head1 CONTACT

    Brett Whitty
    bwhitty@tigr.org

=cut

use warnings;
use strict;
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use lib '/usr/local/devel/ANNOTATION/ard/current/lib/site_perl/5.8.5';
use BSML::BsmlReader;
use BSML::BsmlParserSerialSearch;

my @pos = ();
my %tmhmm = ();
my $seq_id = '';
my %signalps = ();

my %options = ();
my $results = GetOptions (\%options, 
			  'input_list|i=s',
              'database|d=s',
			  'server|s:s',
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

## allow server to be defined
if (!$options{'server'}) {
	$options{'server'} = 'SYBTIGR';
}

## create a new BSML parser
my $parser = new BSML::BsmlParserSerialSearch( SequenceCallBack => \&sequenceHandler,
											   ReadFeatureTables => 1
										     );

## create the database connection
my $dbh = DBI->connect("dbi:Sybase:server=$options{server}; packetSize=8092", 'egc', 'egcpwd', {RaiseError => 1, PrintError => 1});
my $result = $dbh->do("use $options{database}");

## prepare a statement for inserting into ORF_attribute
my $qry = "INSERT INTO ORF_attribute (feat_name, att_type, curated, method, date, assignby, " .
          "score, score_desc, score2, score2_desc, score3, score3_desc, score4, score4_desc, score5, score5_desc) " .
          "VALUES ( ?, 'GES', 0, 'TMHMM2.0', GETDATE(), '$user', ?, 'PredHel', ?, 'Topology', ?, 'ExpAA', ?, 'First60', ?, 'ProbNin' )";
my $orf_attribute_insert = $dbh->prepare($qry);

## prepare a statement for deleting old tmhmm 
$qry = 'DELETE FROM ORF_attribute ' .
       "WHERE att_type = 'GES'" .
       " and method = 'TMHMM2.0'" .
	   " and feat_name = ?";

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
	@pos = ();
	%tmhmm = ();

    _log("processing file $file");

    ## extract the feat_name from the file name:
    ##  bma1.model.14990_07908.tmhmm.bsml -> 14990.m07908
    my $feat_name;
    if ( $file =~ /(\d+)_(\d+)\.tmhmm\.bsml/ ) {
        $feat_name = "$1.m$2";
    } else {
        warn "Improperly named input file, or a proper file without results.  Couldn't extract feat name from $file";
	next;
    }

	## parse the input file
	$parser->parse($file);

	## format Topology string
	$tmhmm{'score2'} = join(":", @pos);

	## if the input sequence was discarded we will have # of hash keys = 0
	## otherwise there will be 4 or 5 depending on whether input setting was plant or animal/yeast
	if (scalar(@pos) == 0) {
		_log("no data added for $feat_name");
	} else {
		## delete the existing tmhmm evidence for this feat_name
	    _log("deleting all tmhmm evidence for $feat_name");
		$delete->execute($feat_name);

		## store the data in the db
	    _log("orf_attribute_insert->execute($feat_name, $tmhmm{score}, $tmhmm{score2}, $tmhmm{score3}, $tmhmm{score4}, $tmhmm{score5});");	
									 
  	  $orf_attribute_insert->execute( 
										$feat_name,
										$tmhmm{'score'},
										$tmhmm{'score2'},
										$tmhmm{'score3'},
										$tmhmm{'score4'},
										$tmhmm{'score5'},
								  	);
	}
}
$delete->finish();
$orf_attribute_insert->finish();
$dbh->disconnect();

exit(0);

sub sequenceHandler {
	my $sequence = shift;
	print $sequence->returnattr('id')."\n";
	## score
	$tmhmm{'score'} = $sequence->{'BsmlAttr'}->{'tmh_count'};
	## score3
	$tmhmm{'score3'} = $sequence->{'BsmlAttr'}->{'exp_aa_in_tmh'};
	## score4
	$tmhmm{'score4'} = $sequence->{'BsmlAttr'}->{'exp_first_60'};
	## score5 
	$tmhmm{'score5'} = $sequence->{'BsmlAttr'}->{'prob_n_in'};
	
	my $ftbl_list_arr_ref = $sequence->returnBsmlFeatureTableListR();
	foreach my $ftbl_ref(@{$ftbl_list_arr_ref}) {
		my $feat_list_arr_ref = $ftbl_ref->returnBsmlFeatureListR();
		foreach my $feat_ref(@{$feat_list_arr_ref}) {
			my $link_arr_ref = $feat_ref->returnBsmlLinkListR();
			#foreach my $link_ref(@{$link_arr_ref}) {
				#	print $link_ref->{'href'}."\n";
			#}
			if ($feat_ref->returnattr('class') eq 'located_sequence_feature') {
				#	print $feat_ref->returnattr('class')." ";
				#	print $feat_ref->returnattr('title')."\n";
				if ($feat_ref->returnattr('title') eq 'TMhelix') {
					my $interval_loc_list_ref = $feat_ref->returnBsmlIntervalLocListR();
					foreach my $interval_loc_ref( @{$interval_loc_list_ref}) {
						push (@pos, $interval_loc_ref->{'startpos'}."-".$interval_loc_ref->{'endpos'});
					}
				}
			}
			#my $attr_hash_ref = $feat_ref->returnBsmlAttrHashR();
			#foreach my $attr(keys(%{$attr_hash_ref})) {
				#	print $attr."=".$attr_hash_ref->{$attr}."\n";
			#}
		}
	}
	return;
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
