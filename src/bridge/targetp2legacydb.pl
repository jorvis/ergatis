#!/usr/local/bin/perl

=head1  NAME 

targetp2legacydb.pl - load targetp data into the legacy database schema.

=head1 SYNOPSIS

USAGE: targetp2legacydb.pl 
        --input_list=/path/to/somefile.raw.list
        --database=sma1
      [ --log=/path/to/some.log ]

=head1 OPTIONS

B<--input_list,-i> 
    raw list file from a targetp workflow component run.

B<--database,-d> 
    Sybase project database ID.

B<--log,-l> 
    optional.  will create a log file with summaries of all actions performed.

B<--help,-h> 
    This help message/documentation.

=head1   DESCRIPTION

This script is used to load targetp evidence into the legacy database schema.
Only raw output files from running targetp on a single input sequence are supported.
Both plant and non-plant output formats should be supported.

=head1 INPUT

The input should be the '*.targetp.raw.list' list file from a targetp workflow component run,
or a similarly formatted list of targetp raw output files with filenames matching the following
regex:
		/(\d+)_(\d+)\.targetp\.raw/
that allows translation to a corresponding feat_name of:
		$1.m$2

=head1 OUTPUT

This is a database loading script. Status messages will be printed to STDERR unless the
log option is used.   

=head1 CONTACT

    Brett Whitty
    bwhitty@tigr.org

=cut

use strict;
use warnings;
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use File::Basename;

my %options = ();
my $results = GetOptions (\%options, 
			  			  'input_list|i=s',
						  'database|d=s',
						  'log|l=s',
						  'help|h',
		  				 ) || pod2usage();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

## make sure all passed options are peachy
check_parameters(\%options);

my $dbname = $options{'database'};
$dbname =~ tr/A-Z/a-z/;

## get the current user
my ($user) = getpwuid($<);

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">$options{log}") || die "Can't create log file: $!";
}

## create the database connection
my $dbh = DBI->connect("dbi:Sybase:server=SYBTIGR; packetSize=8092", 'egc', 'egcpwd', {RaiseError => 1, PrintError => 1});
my $result = $dbh->do("use $options{database}");
	 
my ($time, $date) = time_date();
_log("\t\tSTART: $time  $date\n");

my @model_list;
my @all_model_featnames; #store all models to be processed

open (IN, $options{'input_list'}) || die "Couldn't open input list";
my @raw_list = <IN>;
foreach my $raw_file(@raw_list) {
	chomp $raw_file;
	my $basename = basename($raw_file, ".raw");
	_log("Adding file '$raw_file' for processing\n");
	if ($basename =~ /$dbname\.model\.(\d+)_(\d+)\.targetp/) {
		push(@model_list, [$1.".m".$2, $raw_file]);
	    push(@all_model_featnames, [$1.".m".$2, $1]);
	} else {
		die "Couldn't parse model name from '$basename'\n";
	}
}

my %curated_info;

foreach my $model_info (@all_model_featnames) {
    my ($model_feat_name, $asmbl_id) = @$model_info;
    _log("Deleting old targetp data from '$model_feat_name'\n");
	
	check_curated_info($dbh, $model_feat_name);
	remove_targetP_from_ORF_attribute($dbh,$model_feat_name);

}

my $query = "insert ORF_attribute ("
		  . "feat_name, att_type, curated, method, date, assignby, "
		  . "score, score_desc, score2, score2_desc, score3, "
		  . "score3_desc, score4, score4_desc, score5, score5_desc"
		  . ") values ("
		  . "?, \"targetP\", ?, \"targetP_update.dbi\", getdate(), ?, "
	      . "?, \"Location\", ?, \"c:m:s:o scores\", ?, "
	      . "\"RC-value\", ?, \"network:peplen\", ?, \"c:m:s:o cutoffs\""
		  . ")";   
		
my $db_query = $dbh->prepare($query);

foreach my $model_ref(@model_list) {
	my ($model, $raw_file) = ($model_ref->[0], $model_ref->[1]);
    
	my $tp_ref = parse_targetp($raw_file);
	my %targetp = %{$tp_ref};
	
	my $curated_flag = $curated_info{$model}->{$targetp{'loc'}} || 0;
	
	$db_query->execute(
		$model, 
		$curated_flag, 
		$user, 
		$targetp{'loc'}, 
		"$targetp{ctp}:$targetp{mtp}:$targetp{sp}:$targetp{other}", 
		$targetp{'rc'},
		"$targetp{network}:$targetp{len}", 
		"$targetp{ctp_cutoff}:$targetp{mtp_cutoff}:$targetp{sp_cutoff}:$targetp{other_cutoff}"
					  );
	
	_log("inserting -> ('"
		. join( "', '", (
				$model, 
				$curated_flag, 
				$user, 
				$targetp{'loc'}, 
				"$targetp{ctp}:$targetp{mtp}:$targetp{sp}:$targetp{other}", 
				$targetp{'rc'},
				"$targetp{network}:$targetp{len}", 
				"$targetp{ctp_cutoff}:$targetp{mtp_cutoff}:$targetp{sp_cutoff}:$targetp{other_cutoff}",
			            )
			  )
		. "')\n"
	    );				  
	
	$db_query->finish();
}

($time, $date) = time_date();
_log("\t\tFINISH: $time  $date\n");

$dbh->disconnect();

exit(0);

####
sub parse_targetp {
	my ($raw_file) = @_;

	my %tp;
	my %loc_table = (
						'_' => 'unknown',
						'C' => 'chloroplast',
						'M' => 'mitochondria',
						'S' => 'secreted',
				 	);
	
	open (CMD, $raw_file) || die "Failed opening raw file '$raw_file'";
				
	_log("Parsing '$raw_file'\n");

	while (<CMD>) {
		if (/Using ([^\s]+) networks/){
			$tp{'network'} = $1;
			$tp{'network'} =~ tr/A-Z/a-z/;
		} elsif (/^(?:$dbname\.model\.[^\s]+)	# model_feat_name
		    \s+(\d+)       						# peptide length  $1
		    \s*(\d+\.\d+)? 						# cTP score       $2
		    \s+(\d+\.\d+)  						# mTP score       $3
		    \s+(\d+\.\d+)  						# SP score        $4
		    \s+(\d+\.\d+)  						# other score     $5
		    \s+([_CMS])      					# Location        $6
		    \s+(\d)/x      						# RC-value        $7
		   ) {

		    $tp{'len'} 		= $1;
		    $tp{'ctp'} 		= $2;
		    $tp{'mtp'} 		= $3;
		    $tp{'sp'} 		= $4;
		    $tp{'other'} 	= $5;
		    $tp{'loc'} 		= $6;
		    $tp{'rc'}		= $7;
		
			$tp{'loc'} = $loc_table{$tp{'loc'}};
		} elsif (/^cutoff
			\s*(\d+\.\d+)?
			\s+(\d+\.\d+)
			\s+(\d+\.\d+)
			\s+(\d+\.\d+)/x
			) {
			$tp{'ctp_cutoff'} 	= $1;
			$tp{'mtp_cutoff'} 	= $2;
			$tp{'sp_cutoff'} 	= $3;
			$tp{'other_cutoff'} = $4;
		}	
    }
    close CMD;
	
	return \%tp;
}

####
sub remove_targetP_from_ORF_attribute {
    my ($dbh,$model_feat_name) = @_;
    my $query = "delete ORF_attribute where att_type = \"targetP\" and feat_name = ?";
	my $db_query = $dbh->prepare($query);
	
	$db_query->execute($model_feat_name);
	
	$db_query->finish();
}

####
sub check_curated_info {
    my ($dbh, $model_feat_name) = @_;

	my $query = "select score, curated from ORF_attribute where att_type = \"targetP\" and feat_name = ?";
	my $db_query = $dbh->prepare($query);
	
	$db_query->execute($model_feat_name);
	
    my $row = $db_query->fetchrow_arrayref();
    my $score = $row->[0];
    my $curated = $row->[1];
    
	$db_query->finish();

    if (defined($score)) {
		$curated_info{$model_feat_name}->{$score} = $curated;
    }
}

sub _log {
    my $msg = shift;
    
    if ($logfh) {
		print $logfh "$msg";
	} else {
		print STDERR "$msg";
	}
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

sub time_date {
    my ($sec, $min, $hr, $mday, $m, $year, $wday);
    my (%month, %day);
    ($sec, $min, $hr, $mday, $m, $year, $wday) = localtime(time());
    $month{0}="January"; $month{1}="February"; $month{2}="March";
    $month{3}="April"; $month{4}="May"; $month{5}="June";
    $month{6}="July"; $month{7}="August"; $month{8}="September";
    $month{9}="October"; $month{10}="November"; $month{11}="December";
    $day{0}="Sunday"; $day{1}="Monday"; $day{2}="Tuesday";
    $day{3}="Wednesday"; $day{4}="Thursday"; $day{5}="Friday";
    $day{6}="Saturday";
    
    $time = "$hr:$min:$sec";
    $date = "$day{$wday}, $month{$m} $mday";
    return ($time, $date);
}
