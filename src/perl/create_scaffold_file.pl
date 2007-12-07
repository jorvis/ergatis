#!/usr/bin/perl

=head1  NAME

create_scaffold_file.pl - Converts a set of scaffolds into a tab-delimited scaffold description file.

=head1 SYNOPSIS

USAGE: create_scaffold_file.pl
      --input_file=pva1-scaffolds.txt
      --scaffold_file=pva1.scf
      --db=pva1 
      --username=user1 
      --password=password1 
    [ --server=SYBTIGR ]
    [ --fixed_gap_size=100 ]
    [ --log_file=/path/to/file.log ]
    [ --log_level=4 ]
    [ --help ]
    [ --man ]

=head1 OPTIONS

--input_file
    A plain text description of the scaffolds to create.  See INPUT.

--scaffold_file
    The output file to create; this will be a tab-delimited scaffold description file 
    suitable for input into create_bsml_scaffolds.pl

--db
    Annotation database from which the lengths of the assemblies mentioned in 
    --input_file can be read.

--username
    The name of a user with SELECT permissions on the database specified by --db

--password
    The database password for --username.

--server
    The name of the Sybase server that hosts --annot_db

--fixed_gap_size
    Optional.  If defined this value will override the gap sizes specified in the input
    and a gap of this length (in base pairs) will be placed between each pair of adjacent
    contigs/assemblies in the output scaffold file.

--log_file
    Path to log file to be used by Ergatis::Logger;

--log_level
    Ergatis::Logger log level.

--help
    Displays the usage statement.   

--man
    Displays this man page.


=head1   DESCRIPTION

Converts a set of scaffolds into a tab-delimited scaffold description file.

=head1  INPUT

Scaffold descriptions are often e-mailed around in the following form:

chromo 1:     B396E...4kb...B399E
chromo 2:     B402E...3kb...B392E
chromo 3:     E2922B...6kb...E401B...1.5kb...B408E..7kb...B2920E...4kb...B412E...0kb...B422E
chromo 4:     B414E...1.5kb...B2924E...1.8kb.....E2926B
chromo 5:     B2927E
chromo 6:     B2923E...3kb...E2918B
chromo 7:    E404B...B2919E 
chromo 8:    E400B....0.2kb...B429E                        
chromo 9:     B398E
chromo 10:    B391E....2.8kb...E403B
chromo 11:    E426B....4kb...E434B
chromo 12:    E393B...6kb...B428E
chromo 13:    B394E
chromo 14:    B2929E..3kb...B405E 

"B" and "E" indicate the beginning and end of each contig and the number in the
middle is the annotation database asmbl_id.  The size of aps may or may not be 
indicated.

=head1  OUTPUT

A tab-delimited scaffold description file suitable for input to create_bsml_scaffolds.pl
For example, here is the output that corresponds to INPUT:

SCF	396	1	1	F	0	0	565852	0	565852
SCF	399	1	2	F	0	0	260516	569852	830368
SCF	402	2	1	F	0	0	162059	0	162059
SCF	392	2	2	F	0	0	589976	165059	755035
SCF	2922	3	1	R	0	0	533272	533272	0
SCF	401	3	2	R	0	0	388176	927448	539272
SCF	408	3	3	F	0	0	13965	928948	942913
SCF	2920	3	4	F	0	0	35881	949913	985794
SCF	412	3	5	F	0	0	10181	989794	999975
SCF	422	3	6	F	0	0	10652	999975	1010627
SCF	414	4	1	F	0	0	15769	0	15769
SCF	2924	4	2	F	0	0	413554	17269	430823
SCF	2926	4	3	R	0	0	444168	876791	432623
SCF	2927	5	1	F	0	0	1370936	0	1370936
SCF	2923	6	1	F	0	0	329199	0	329199
SCF	2918	6	2	R	0	0	701189	1033388	332199
SCF	404	7	1	R	0	0	1198945	1198945	0
SCF	2919	7	2	F	0	0	298774	1198945	1497719
SCF	400	8	1	R	0	0	1165049	1165049	0
SCF	429	8	2	F	0	0	513347	1165249	1678596
SCF	398	9	1	F	0	0	1923364	0	1923364
SCF	391	10	1	F	0	0	895497	0	895497
SCF	403	10	2	R	0	0	521442	1419739	898297
SCF	426	11	1	R	0	0	2021996	2021996	0
SCF	434	11	2	R	0	0	41358	2067354	2025996
SCF	393	12	1	R	0	0	1012632	1012632	0
SCF	428	12	2	F	0	0	1986252	1018632	3004884
SCF	394	13	1	F	0	0	2031768	0	2031768
SCF	2929	14	1	F	0	0	2132794	0	2132794
SCF	405	14	2	F	0	0	984623	2135794	3120417

See create_bsml_scaffolds.pl for details.

=head1  CONTACT

    Jonathan Crabtree
    crabtree@tigr.org

=cut

# standard ergatis script 
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

use strict;
use warnings;

use DBI;
use FileHandle;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Ergatis::Logger;
use Ergatis::IdGenerator;

# ------------------------------------------------------------------
# Input
# ------------------------------------------------------------------

my($input_file,
   $scaffold_file,
   $db,
   $username,
   $password,
   $server,
   $fixed_gap_size,
   $log_file,
   $log_level,
   $help,
   $man
   );

&GetOptions("input_file=s" => \$input_file, 
	    "scaffold_file=s" => \$scaffold_file,
	    "db=s" => \$db,
	    "server=s" => \$server,
	    "username=s" => \$username,
	    "password=s" => \$password,
	    "server=s" => \$server,
	    "fixed_gap_size=i" => \$fixed_gap_size,
	    "log_file=s" => \$log_file,
	    "log_level=s" => \$log_level,
	    "help" => \$help,
	    "man" => \$man,
	    );

# print usage/documentation
pod2usage(1) if $help;
pod2usage({-verbose => 2}) if $man;

if (!$input_file || !(-e $input_file) || !(-r $input_file)) {
    pod2usage({-message => "Error:\n     --input_file is not readable\n", -exitstatus => 1, -verbose => 0});
}

pod2usage({-message => "Error:\n     --db must be specified\n", -exitstatus => 1, -verbose => 0}) if (!$db);
pod2usage({-message => "Error:\n     --username must be specified\n", -exitstatus => 1, -verbose => 0}) if (!$username);
pod2usage({-message => "Error:\n     --password must be specified\n", -exitstatus => 1, -verbose => 0}) if (!$password);

# defaults
$server = 'SYBTIGR' if (!defined($server));

# logger
$log_file = Ergatis::Logger::get_default_logfilename() if (!defined($log_file));
my $elogger = new Ergatis::Logger('LOG_FILE'=>$log_file, 'LOG_LEVEL'=>$log_level);
my $logger = $elogger->get_logger(__PACKAGE__);

# ------------------------------------------------------------------
# Main program
# ------------------------------------------------------------------

# read assembly lengths from the annotation database
my $dbi_dsn = "DBI:Sybase:server=$server;packetSize=8192";
my $dbh = DBI->connect($dbi_dsn, $username, $password);
my $asmblId2Length = &readAssemblyLengths($dbh, $db);

# parse the input file
my $ifh = FileHandle->new();
$ifh->open($input_file) || die "unable to read from $input_file";
my $lnum = 0;

# write output to scaffold file
my $ofh = FileHandle->new();
$ofh->open(">$scaffold_file") || die "unable to write to $scaffold_file";

while (my $line = <$ifh>) {
    chomp($line);
    ++$lnum;

    if ($line =~ /^chromo (\d+):\s+(\S+)/) {
	my $chrom = $1;
	my $content = $2;
	my @elts = split(/\.\.+/, $content);
	$logger->debug("chrom=$chrom content=$content elts=" . join(',', @elts));

	my $index = 1;
	my $last_fmax = 0;
	my $ne = scalar(@elts);
	
	for (my $i = 0;$i < $ne;++$i) {
	    my $elt = $elts[$i];
	    my $is_first = ($i == 0) ? 1 : 0;
	    my $is_last = ($i == ($ne-1)) ? 1 : 0;
	    $logger->debug("elt=$elt i=$i is_first=$is_first is_last=$is_last");
	    
	    # gap
	    if ($elt =~ /^([\d\.]+)kb$/) {
		# only use gap size from file if no fixed gap size defined
		if (!defined($fixed_gap_size)) {
		    my $gapSize = $1 * 1000;
		    $logger->debug(" gap=$gapSize");
		    $last_fmax += $gapSize;
		}
	    }
	    # contig
	    else {
		my $asmbl_id;
		my $is_rev = 0;
		
		# fixed gap size
		if (($i > 1) && (defined($fixed_gap_size))) {
		    $logger->debug(" fixed gap=$fixed_gap_size");
		    $last_fmax += $fixed_gap_size;
		}
		
		# one way to specify the orientation of a single contig: "408B--408E"
		if ($elt =~ /^([BE])?(\d+)([BE])\-\-([BE])?(\d+)([BE])$/) {
		    my ($e1o,$a1,$e1,$e2o,$a2,$e2) = ($1,$2,$3,$4);
		    die "asmbl_id mismatch in $elt" if ($a1 != $a2);
		    die "contig end mismatch in $elt" if ($e1 eq $e2);
		    die "nonsense end specification in $elt" if (defined($e1o) && ($e1o eq $e1));
		    die "nonsense end specification in $elt" if (defined($e2o) && ($e2o eq $e2));
		    $is_rev = 1 if ($e1 eq 'E');
		    $asmbl_id = $a1;
		} 

		# an alternative: "B408E" or "E408B"
		elsif ($elt =~ /([BE])?(\d+)([BE])$/) {
		    my $e1o = $1;
		    $asmbl_id = $2;
		    my $e1 = $3;
		    $is_rev = 1 if ($e1 eq 'B');
		    die "nonsense end specification in $elt" if (defined($e1o) && ($e1o eq $e1));
		}
		
		# "408B"
		elsif ($elt =~ /^(\d+)([BE])$/) {
		    $asmbl_id = $1;
		    my $e1 = $2;
		    $is_rev = 1 if ($is_first && ($e1 eq 'B'));
		    $is_rev = 1 if ($is_last && ($e1 eq 'E'));
		    die "single element scaffold spec ('$elt') at internal node" if (!$is_first && !$is_last);
		}

		# "408"
		elsif ($elt =~ /^(\d+)$/) {
		    $asmbl_id = $1;
		    # not specifying the orientation is only allowed if this is a singleton scaffold, in 
		    # which case we assume that the contig is in the forward orientation
		    die "bare element scaffold spec in non-trivial scaffold" if (!$is_first || !$is_last);
		}

		else {
		    die "unable to parse scaffold element '$elt'";
		}
		
		my $seqlen = $asmblId2Length->{$asmbl_id} || die "unable to get sequence length for asmbl $asmbl_id";
		my $sc = $is_rev ? 'R' : 'F';
		my $fmin = $last_fmax;
		my $fmax = $fmin + $seqlen;

		if ($is_rev) {
		    $ofh->print(join("\t", ('SCF', $asmbl_id, $chrom, $index, $sc, '0', '0', $seqlen, $fmax, $fmin)) . "\n");
		} else {
		    $ofh->print(join("\t", ('SCF', $asmbl_id, $chrom, $index, $sc, '0', '0', $seqlen, $fmin, $fmax)) . "\n") ;
		}

		$last_fmax = $fmax;
		++$index;
	    }
	}
    } else {
	die "unable to parse line $lnum of $input_file: '$line'";
    }
}

$ifh->close();
$ofh->close();

exit(0);

# ------------------------------------------------------------------
# Subroutines
# ------------------------------------------------------------------

sub readAssemblyLengths {
    my($dbh, $db) = @_;
    my $sql = "SELECT asmbl_id, char_length(sequence) " .
	"FROM ${db}..assembly ";

    my $sth = $dbh->prepare($sql);
    $sth->execute();
    my $id2len = {};

    while(my($asmblId, $length) = $sth->fetchrow_array()) {
	$id2len->{$asmblId} = $length;
    }

    $sth->finish();
    return $id2len;
}

__END__
