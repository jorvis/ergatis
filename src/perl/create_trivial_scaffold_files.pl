#!/usr/local/bin/perl

=head1  NAME

create_trivial_scaffold_files.pl - Create a tab-delimited scaffold description file for each set of assemblies not already in a scaffold.

=head1 SYNOPSIS

USAGE: create_trivial_scaffold_files.pl
      --project=apx3
      --scaffold_dir=/usr/local/annotation/APX3/scaffolds
      --chado_db=apx3
      --username=user1 
      --password=password1 
    [ --output_dir=. ]
    [ --server=SYBTIGR ]
    [ --log_file=/path/to/file.log ]
    [ --log_level=4 ]
    [ --help ]
    [ --man ]

=head1 OPTIONS

--project
    The name of the project for which to create the scaffold files.

--scaffold_dir
   The full path to the top-level directory that contains all the BSML scaffold documents for --project.

--chado_db
    The name of the chado comparative database to query for the list of assembly/contig sequences.

--server
    The name of the Sybase server that hosts --chado_db

--username
    The name of a user with SELECT permissions on the database specified by --chado_db

--password
    The database password for --username.

--output_dir
    The directory where the scaffold files should be written.

--log_file
    Path to log file to be used by Ergatis::Logger.

--log_level
    Ergatis::Logger log level.

--help
    Displays the usage statement.   

--man
    Displays this man page.

=head1   DESCRIPTION

Create a tab-delimited scaffold description file for each organism in --chado_db that has 
assembly/contig sequences that are not already part of a scaffold in --scaffold_dir.
Each such assembly will be assigned to a unique (and hence "trivial") scaffold.

=head1  OUTPUT

A set of tab-delimited scaffold description files, one for each organism in --chado_db
that has assembly sequences not currently assigned to a scaffold by a BSML scaffold
document in --scaffold_dir.  These scaffold description files may be used as the
input to create_bsml_scaffolds.pl in order to create a BSML document for each trivial
scaffold.

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
use File::Spec;
use Getopt::Long;
use Pod::Usage;
use Ergatis::Logger;

# ------------------------------------------------------------------
# Input
# ------------------------------------------------------------------

my($project,
   $scaffold_dir,
   $chado_db,
   $username,
   $password,
   $output_dir,
   $server,
   $log_file,
   $log_level,
   $help,
   $man);

&GetOptions("project=s" => \$project,
	    "scaffold_dir=s" => \$scaffold_dir,
	    "chado_db=s" => \$chado_db, 
	    "username=s" => \$username,
	    "password=s" => \$password,
	    "output_dir=s" => \$output_dir,
	    "server=s" => \$server,
	    "log_file=s" => \$log_file,
	    "log_level=s" => \$log_level,
	    "help" => \$help,
	    "man" => \$man,
	    );

pod2usage(1) if $help;
pod2usage({-verbose => 2}) if $man;

if (!$scaffold_dir || !(-d $scaffold_dir) || !(-r $scaffold_dir)) {
    pod2usage({-message => "Error:\n     --scaffold_dir=$scaffold_dir is not readable\n", -exitstatus => 1, -verbose => 0});
}

pod2usage({-message => "Error:\n     --project must be specified\n", -exitstatus => 1, -verbose => 0}) if (!$project);
pod2usage({-message => "Error:\n     --chado_db must be specified\n", -exitstatus => 1, -verbose => 0}) if (!$chado_db);
pod2usage({-message => "Error:\n     --username must be specified\n", -exitstatus => 1, -verbose => 0}) if (!$username);
pod2usage({-message => "Error:\n     --password must be specified\n", -exitstatus => 1, -verbose => 0}) if (!$password);

# default values
# server
$server = 'SYBTIGR' if (!defined($server));

# logger
$log_file = Ergatis::Logger::get_default_logfilename() if (!defined($log_file));
my $elogger = new Ergatis::Logger('LOG_FILE'=>$log_file, 'LOG_LEVEL'=>$log_level);
my $logger = $elogger->get_logger(__PACKAGE__);

# output_dir
$output_dir = "." if (!defined($output_dir));

# ------------------------------------------------------------------
# Main program
# ------------------------------------------------------------------

# read assembly names from the existing scaffold files
my $assemblyIds = &readAssemblyIdsFromScaffoldFiles($scaffold_dir);
$logger->info("read " . scalar(@$assemblyIds) .  " assembly ids from $scaffold_dir");

# read assembly names and lengths from the database
my $dbiDsn = "DBI:Sybase:server=$server;packetSize=8192";
my $chadoAssems = &readAssembliesFromChadoDb($dbiDsn, $username, $password, $chado_db);
$logger->info("read " . scalar(@$chadoAssems) . " assemblies from $chado_db");

# get list of assemblies that are not already in a scaffold
my $ids = {};
map { $ids->{$_} = 1; } @$assemblyIds;
my $assemblies = [];
map { push(@$assemblies, $_) if (!defined($ids->{$_->{'uniquename'}})); } @$chadoAssems;
$logger->info("generating trivial scaffolds for " . scalar(@$assemblies) . " assemblies:");

# group by organism/database
my $assemsByDb = {};
foreach my $assem (@$assemblies) {
    my $db = $assem->{'annot_db'};
    my $list = $assemsByDb->{$db};
    $list = $assemsByDb->{$db} = [] if (!defined($list));
    push(@$list, $assem);
}

foreach my $db (keys %$assemsByDb) {
    my $assems = $assemsByDb->{$db};
    my $na = scalar(@$assems);
    $logger->info(" $db: $na assemblies");
}

# generate a scaffold file for each source database/genome
foreach my $db (keys %$assemsByDb) {
    my $assems = $assemsByDb->{$db};
    my $scaffoldFile = File::Spec->catfile($output_dir, "${db}-trivial.scf");
    $logger->info("generating $scaffoldFile with " . scalar(@$assems) . " $db assemblies");
    my $fh = FileHandle->new();
    $fh->open(">$scaffoldFile") || die "unable to write to $scaffoldFile";
    
    foreach my $assem (@$assems) {
	my $asmblId = $assem->{'asmbl_id'};
	my $asmblLen = $assem->{'length'};
	$fh->printf("SCF %10s %10s %5s F     0     0 %10s %10s %10s\n", $asmblId, $asmblId, '1', $asmblLen, 0, $asmblLen);
    }

    $fh->close();
}

# all done
exit(0);

# ------------------------------------------------------------------
# Subroutines
# ------------------------------------------------------------------

# Returns a listref of those assemblies (by assembly name) already included in a scaffold.
#
sub readAssemblyIdsFromScaffoldFiles {
    my($scaffold_dir) = @_;
    my $names = [];

    my $cmd = "find ${scaffold_dir} -name '*.bsml' -exec egrep '\.assembly' '{}' ';' |";
    my $fh = FileHandle->new();
    $fh->open($cmd);
    while (my $line = <$fh>) {
	if ($line =~ /id=\"(\S+\.assembly\.\d+\.\d+)\"/) {
	    push(@$names, $1);
	} else {
	    die "could not parse output of find command: '$line'";
	}
    }
    $fh->close();
    return $names;
}

# Returns a listref of assembly sequence names in $chado_db.
#
sub readAssembliesFromChadoDb {
    my($dbiDsn, $username, $password, $chado_db) = @_;

    my $sql = "select f.uniquename, f.seqlen, dbx.accession, db.name " .
	"from ${chado_db}..feature f, ${chado_db}..cvterm cv, ${chado_db}..dbxref dbx, ${chado_db}..db db " .
	"where f.type_id = cv.cvterm_id " .
	"and cv.name = 'assembly' " .
	"and f.dbxref_id = dbx.dbxref_id " .
	"and dbx.db_id = db.db_id " .
	"and db.name like 'TIGR_%' ";

    my $assemblies = [];
    my $dbh = DBI->connect($dbiDsn, $username, $password);
    my $sth = $dbh->prepare($sql);
    $sth->execute();

    while(my($name, $seqlen, $asmblId, $db) = $sth->fetchrow_array()) {
	my($annotDb) = ($db =~ /TIGR_(\S+)$/);
	push(@$assemblies, {
	    'uniquename' => $name, 
	    'length' => $seqlen, 
	    'asmbl_id' => $asmblId, 
	    'annot_db' => lc($annotDb),
	}); 
    }

    $sth->finish();
    $dbh->disconnect();
    return $assemblies;
}

__END__
