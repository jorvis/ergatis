#!/usr/bin/perl

=head1  NAME

create_bsml_scaffolds.pl - Convert a tab-delimited scaffold description file into a set of BSML scaffold documents.

=head1 SYNOPSIS

USAGE: create_bsml_scaffolds.pl
      --project=apx3
      --annot_db=pva1 
      --scaffold_file=pva1.scf
      --id_repository=/usr/local/annotation/APX3/workflow/project_id_repository
      --assembly_bsml_list=/usr/local/annotation/APX3/output_repository/legacy2bsml/7881_default/legacy2bsml.bsml.list
      --organism_name="Plasmodium vivax"
      --organism_type=euk
    [ --output_dir=. ]
    [ --parse_gaps ]
    [ --parse_only ]
    [ --using_contig_ids ]
    [ --scaffold_id_is_chromosome ]
    [ --annot_server=SYBTIGR ]
    [ --annot_username=user1 ]
    [ --annot_password=password1 ]
    [ --log_file=/path/to/file.log ]
    [ --log_level=4 ]
    [ --help ]
    [ --man ]

=head1 OPTIONS

--project
    The name of the project for which to create the BSML scaffold documents.

--scaffold_file
   The tab-delimited flat file that contains the descriptions of the scaffolds to be created.

--id_repository
   A valid id repository.  See Ergatis::IdGenerator for details.

--assembly_bsml_list
   The full path to the file that contains a list of all the assembly BSML files (e.g., the
   legacy2bsml.bsml.list file produced by the legacy2bsml process used to populate the chado
   comparative database that corresponds to --project.)

--organism_name
   The full name ("Genus species") of the organism to which the scaffolds should be assigned.

--organism_type
   Either 'euk' or 'prok'; the type of TIGR annotation database used to store --organism_name.

--output_dir
   Where to write the scaffold BSML documents.

--parse_gaps
   Use this option to specify that the gaps between contigs should be determined by the coordinates
   supplied in the --scaffold_file.  If this option is not set then a uniform gap size of 100bp will
   be used to separate each pair of adjacent contigs/assemblies.  Note that the script will complain
   if the --parse_gaps option is not given and the coordinates in the file do not match the 100bp
   gap assumption.

--parse_only
   Parse the input file(s) and report any errors but don't create any BSML documents.

--using_contig_ids
   Whether the scaffold file uses contig ids (clone_info.clone_name) instead of asmbl_ids.  
   If set to true then --annot_username and --annot_password must also be supplied.

--scaffold_id_is_chromosome
   Use this option if the scaffold ids in the --scaffold_file should be saved as the scaffold's
   chromosome in the BSML file.

--annot_db
    The name of the TIGR annotation database from which to read the asmbl_id -> clone_name mapping.

--annot_server
    The name of the Sybase server that hosts --annot_db

--annot_username
    The name of a user with SELECT permissions on the database specified by --annot_db

--annot_password
    The database password for --annot_username.

--log_file
    Path to log file to be used by Ergatis::Logger;

--log_level
    Ergatis::Logger log level.

--help
    Displays the usage statement.   

--man
    Displays this man page.

=head1   DESCRIPTION

Convert a tab-delimited scaffold description file into a set of BSML scaffold documents.

=head1  INPUT

The input --scaffold_file should look like this:

SCF 396	1	1	F	0	0	565852	0	565852
SCF	399	1	2	F	0	0	260516	569852	830368
SCF	402	2	1	F	0	0	162059	0	162059
SCF	392	2	2	F	0	0	589976	165059	755035
SCF	2922	3	1	R	0	0	533272	533272	0
SCF	401	3	2	R	0	0	388176	927448	539272
SCF	408	3	3	F	0	0	13965	928948	942913
SCF	2920	3	4	F	0	0	35881	949913	985794
SCF	412	3	5	F	0	0	10181	989794	999975
SCF	422	3	6	F	0	0	10652	999975	1010627

The columns are tab-delimited (actually whitespace-delimited) and they are interpreted as follows:

1. SCF - This column always contains the string "SCF".  This is a historical artifact.
2. contig_id - The asmbl_id or clone_info.clone_name (if --using_contig_ids) of an 
               assembly in the scaffold.
3. scaff_id - The id of the scaffold to which this contig belongs.
4. scaff_order - An integer that indicates the order of the contig in the scaffold, counting
                 from 1 at the 5' end of the scaffold sequence.
5. direction - The orientation of the contig - either 'F' for forward or 'R' for reverse.
6. gap_mean - Mean size of the gap between this contig and the next; not currently used.
7. gap_std - Standard deviation of the gap between this contig and the next; not currently used.
8. contig_len - Length of the contig.
9. contig_start - Location of the 5' end of the contig on the scaffold in chado interbase coordinates
                  (numbering from 0)
10. contig_end - Location of the 3' end of the contig on the scaffold in chado interbase coordinates
                  (numbering from 0)

=head1  OUTPUT

The script outputs a BSML-format scaffold file for each scaff_id in the INPUT.

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

use BSML::BsmlBuilder;
use BSML::BsmlReader;
use BSML::BsmlParserTwig;
use BSML::BsmlScaffolds;

# ------------------------------------------------------------------
# Input
# ------------------------------------------------------------------

my($project,
   $scaffold_file,
   $id_repository,
   $assembly_bsml_list,
   $organism_name,
   $organism_type,
   $output_dir,
   $parse_gaps,
   $parse_only,
   $using_contig_ids,
   $scaffold_id_is_chromosome,
   $annot_db,
   $annot_server,
   $annot_username,
   $annot_password,
   $log_file,
   $log_level,
   $help,
   $man);

&GetOptions("project=s" => \$project, 
	    "scaffold_file=s" => \$scaffold_file,
	    "id_repository=s" => \$id_repository,
	    "assembly_bsml_list=s" => \$assembly_bsml_list,
	    "organism_name=s" => \$organism_name,
	    "organism_type=s" => \$organism_type,
	    "output_dir=s" => \$output_dir,
	    "parse_gaps!" => \$parse_gaps,
	    "parse_only!" => \$parse_only,
	    "using_contig_ids!" => \$using_contig_ids,
	    "scaffold_id_is_chromosome!" => \$scaffold_id_is_chromosome,
	    "annot_db=s" => \$annot_db,
	    "annot_server=s" => \$annot_server,
	    "annot_username=s" => \$annot_username,
	    "annot_password=s" => \$annot_password,
	    "log_file=s" => \$log_file,
	    "log_level=s" => \$log_level,
	    "help" => \$help,
	    "man" => \$man,
	    );

# print usage/documentation
pod2usage(1) if $help;
pod2usage({-verbose => 2}) if $man;

if (!$id_repository || !(-d $id_repository) || !(-r $id_repository)) {
    pod2usage({-message => "Error:\n     --id_repository=$id_repository is not readable\n", -exitstatus => 1, -verbose => 0});
}
if (!$assembly_bsml_list || !(-r $assembly_bsml_list)) {
    pod2usage({-message => "Error:\n     --assembly_bsml_list=$assembly_bsml_list is not readable\n", -exitstatus => 1, -verbose => 0});
}
unless ($organism_type =~ /^euk$/i || $organism_type =~ /^prok$/) {
    pod2usage({-message => "Error:\n     --organism_type must be 'euk' or 'prok'\n", -exitstatus => 1, -verbose => 0});
}
if (!$scaffold_file || !(-e $scaffold_file) || !(-r $scaffold_file)) {
    pod2usage({-message => "Error:\n     --scaffold_file is not readable\n", -exitstatus => 1, -verbose => 0});
}
pod2usage({-message => "Error:\n     --project must be specified\n", -exitstatus => 1, -verbose => 0}) if (!$project);
pod2usage({-message => "Error:\n     --annot_db must be specified\n", -exitstatus => 1, -verbose => 0}) if (!$annot_db);
pod2usage({-message => "Error:\n     --organism_name must be specified\n", -exitstatus => 1, -verbose => 0}) if (!$organism_name);

# only need database login to convert contig_ids to asmbl_ids
if ($using_contig_ids) {
    pod2usage({-message => "Error:\n     --annot_username must be specified\n", -exitstatus => 1, -verbose => 0}) if (!$annot_username);
    pod2usage({-message => "Error:\n     --annot_password must be specified\n", -exitstatus => 1, -verbose => 0}) if (!$annot_password);
    pod2usage({-message => "Error:\n     --annot_server must be specified\n", -exitstatus => 1, -verbose => 0}) if (!$annot_server);
}

# parse $organism_name into chado organism.genus, species
$organism_name =~ s/^\s+|\s+$//g;
my($genus, $species) = split(" ", $organism_name, 2);

if ($genus eq '' || $species eq '') {
    pod2usage({-message => "Error:\n     --organism_name must be specified as 'genus species'\n", -exitstatus => 1, -verbose => 0});
}

# logger
$log_file = Ergatis::Logger::get_default_logfilename() if (!defined($log_file));
my $elogger = new Ergatis::Logger('LOG_FILE'=>$log_file, 'LOG_LEVEL'=>$log_level);
my $logger = $elogger->get_logger(__PACKAGE__);

# default output dir
$output_dir = "." if (!defined($output_dir));

# ------------------------------------------------------------------
# Main program
# ------------------------------------------------------------------

$| = 1;

# read scaffolds from flat file
my $scaffolds = &readScaffoldsFromFlatFile($scaffold_file);
my $nContigs = 0;
my $chash = {};
map { $nContigs += scalar(@{$_->{'data'}}); } @$scaffolds;
$logger->info("read " . scalar(@$scaffolds) . " scaffold(s) and $nContigs contig(s) from $scaffold_file");

# using contig ids: read contig_id -> asmbl_id mapping from --annot_db
if ($using_contig_ids) {
    my $dbiDsn = "DBI:Sybase:server=$annot_server;packetSize=8192";
    my $contig2Asmbl = &readContigToAsmblMapping($dbiDsn, $annot_username, $annot_password, $annot_db);
    $logger->info("read " . scalar(keys(%$contig2Asmbl)) . " asmbl_ids from $annot_db database");

    # map contig_ids to asmbl_ids
    foreach my $s (@$scaffolds) {
	foreach my $c (@{$s->{'data'}}) {
	    my $asmblId = $contig2Asmbl->{$c->{'contig_id'}};
	    if (!defined($asmblId)) {
		die "unable to map contig_id $c->{'contig_id'} to asmbl_id";
	    } else {
		$c->{'asmbl_id'} = $asmblId;
	    }
	}
    }
}

# otherwise the contig_id is the same as the asmbl_id
else {
    foreach my $s (@$scaffolds) {
	foreach my $c (@{$s->{'data'}}) {
	    $c->{'asmbl_id'} = $c->{'contig_id'};
	}
    }
}

# use IdGenerator to generate new scaffold ids
my $id_gen = Ergatis::IdGenerator->new( 'id_repository' => $id_repository, logging => 0);
$id_gen->set_pool_size( 'supercontig' => 1 );
my $new_id_fn = sub {
    # keep --annot_db as the prefix of the new ids (for convenience)
    return $id_gen->next_id('type' => 'supercontig', project => lc($annot_db) );
};

if ($scaffold_id_is_chromosome) {
    map { $_->{'chromosome'} = $_->{'scaff_id'}; } @$scaffolds;
}

# write BSML files
&BSML::BsmlScaffolds::writeBsmlScaffoldFiles($scaffolds, $annot_db, $organism_type, $project, $assembly_bsml_list,
					     $output_dir, $genus, $species, $parse_gaps, $parse_only, $new_id_fn, $elogger);

# all done
exit(0);

# ------------------------------------------------------------------
# Subroutines
# ------------------------------------------------------------------

# Returns a listref of scaffolds read from the tab-delimited flat file.
# 
# $filename - path to the scaffold file
#
sub readScaffoldsFromFlatFile {
    my($filename) = @_;
    my $scaffIdToContigs = {};
    my $contigIds = {}; # used to check that contig ids are unique

    my $fh = FileHandle->new();
    $fh->open($filename, 'r') || return undef;
    my $lnum = 0;
    while (my $line = <$fh>) {
	++$lnum; chomp($line);
	# skip whitespace
	next if ($line =~ /^\s*$/);
	# could just split on /\s+/, but this adds some error-checking
	if ($line =~ /^SCF\s+(\d+)\s+(\S+)\s+(\d+)\s+(F|R)\s+(-?\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)$/) {
	    my $data = {
		'contig_id' => $1,
		'scaff_id' => $2,
		'scaff_order' => $3,
		'direction' => $4,
		'gap_mean' => $5, # not used
		'gap_std' => $6, # not used
		'contig_len' => $7,
		'contig_start' => $8,
		'contig_end' => $9,
	    };
	    my $scaffId = $2;
	    my $contigId = $1;

	    my $list = $scaffIdToContigs->{$scaffId};
	    if (!defined($list)) { $list = $scaffIdToContigs->{$scaffId} = []; }
	    push(@$list, $data);
	    if (defined($contigIds->{$contigId})) {
		die "contig_id $contigId is duplicated at line $lnum";
	    } else {
		$contigIds->{$contigId} = 1;
	    }
	} else {
	    die "unable to parse line $lnum of $filename: $line";
	}
    }
    $fh->close();
    $logger->info("read $lnum lines from $filename");

    # process scaffolds
    my $scaffolds = [];

    foreach my $scaffId (keys %$scaffIdToContigs) {
	my $contigs = $scaffIdToContigs->{$scaffId};
	my @sorted = sort { $a->{'scaff_order'} <=> $b->{'scaff_order'}} @$contigs;
	push(@$scaffolds, { 'scaff_id' => $scaffId, 'data' => \@sorted });
    }
    
    return $scaffolds;
}

# Read contig_id -> asmbl_id mapping from annot_db.
#
sub readContigToAsmblMapping {
    my($dbi_dsn, $annot_username, $annot_password, $annot_db) = @_;
    my $mapping = {};
    my $sql = "select distinct clone_name, asmbl_id from ${annot_db}..clone_info where asmbl_id is not null";

    my $dbh = DBI->connect($dbi_dsn, $annot_username, $annot_password);
    my $sth = $dbh->prepare($sql);
    $sth->execute();
    while(my($contigId, $asmblId) = $sth->fetchrow_array()) {
	my $currVal = $mapping->{$contigId};
	if (defined($currVal) && ($currVal ne $asmblId)) {
	    die "clone_name $contigId maps to multiple distinct asmbl_ids in $annot_db";
	} else {
	    $mapping->{$contigId} = $asmblId;
	}
    }
    $sth->finish();
    $dbh->disconnect();

    return $mapping;
}

__END__
