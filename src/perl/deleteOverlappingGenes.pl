#!/usr/bin/env perl

=head1  NAME

deleteOverlappingGenes.pl - Parse RNA_CDS_OVERLAP errors from the tbl2asn discrepancy file and use the Manatee delete_gene.cgi script to delete the corresponding genes.

=head1 SYNOPSIS

deleteOverlappingGenes.pl
         --db_mapping_file=/path/to/db_map.txt
         --discrep_file=/path/to/discrep.txt
         --user=manatee_username
         --password=manatee_password
         --server=jabba.igs.umaryland.edu
         --delete_url=http://
        [--log=/path/to/debug/logfile.txt
         --no_deletes
         --help
         --man]

=head1 OPTIONS

B<--db_mapping_file,-d>
    File that maps databases to other various fields, Must be a comma or tab-delimited file and the first field must be the db name

    Example:
    database1,id_prefix,/path/to/rules/txt

B<--discrep_file,-f>
    Path to a particular tbl2asn discrepancy file.

B<--user,-u>
    Manatee username.

B<--password,-p>
    Manatee password.

B<--server,-s>
    Manatee database server.

B<--delete_url,-U>
    URL of the Manatee delete_gene.cgi script that can be used to delete the offending gene(s)

B<--log,-l>
    Optional. Path to a log file into which to write detailed (DEBUG-level) logging info.

B<--no_deletes,-n>
    Optional.  Run without executing the delete commands.

B<--help,-h>
    Display the script documentation.

B<--man,-m>
    Display the script documentation.

=head1 DESCRIPTION

Parse RNA_CDS_OVERLAP errors from the tbl2asn discrepancy file and use the Manatee delete_gene.cgi script to delete the corresponding genes.

=head1 INPUT

A discrepancy file from tbl2asn that contains one or more errors that match 'DiscRep_SUB:RNA_CDS_OVERLAP'.

=head1 OUTPUT

Prints the number of genes successfully deleted.

=head1 CONTACT

    Jonathan Crabtree
    jonathancrabtree@gmail.com

=cut

use strict;
use FileHandle;
use Carp;
use Data::Dumper;
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Log::Log4perl qw(:easy);
use LWP::UserAgent;

## globals
my $TRANSCRIPT_ID_QUERY =
   "SELECT t.uniquename " .
   "FROM dbxref dbx, feature_dbxref fdbx, feature_relationship fr, feature t, cvterm cvt " .
   "WHERE dbx.accession = ? " .
   "AND dbx.version = 'public_locus' " .
   "AND dbx.dbxref_id = fdbx.dbxref_id " .
   "AND fdbx.feature_id = fr.object_id " .
   "AND t.feature_id = fr.subject_id " .
   "AND t.type_id = cvt.cvterm_id " .
   "AND cvt.name = 'transcript' ";

my $RNA_ID_QUERY =
   "SELECT t.uniquename " .
   "FROM dbxref dbx, feature_dbxref fdbx, feature_relationship fr, feature t, cvterm cvt " .
   "WHERE dbx.accession = ? " .
   "AND dbx.version = 'public_locus' " .
   "AND dbx.dbxref_id = fdbx.dbxref_id " .
   "AND fdbx.feature_id = fr.object_id " .
   "AND t.feature_id = fr.subject_id " .
   "AND t.type_id = cvt.cvterm_id " .
   "AND cvt.name IN ('rRNA', 'tRNA') ";

## input
my $options = {};

&GetOptions($options,
            "db_mapping_file|d=s",
            "discrep_file|f=s",
            "user|u=s",
            "password|p=s",
            "password_file|P=s",
            "server|s=s",
            "delete_url|U=s",
            "no_deletes!",
            "log|l=s",
            "help|h",
            "man|m",
           ) || pod2usage();

## display documentation
if ( $options->{'help'} || $options->{'man'} ) {
  pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

my $password;
my $overlap_flag =0;
my %db_map;	#maps database to discrepancy file path

my $cds_ids;
my $genes_to_delete;
my $rnas_to_delete;

&check_parameters($options);

## logging
my $log4perl_conf = q(
                      log4perl.rootLogger                = INFO, screen
                      log4perl.appender.screen           = Log::Log4perl::Appender::Screen
                      log4perl.appender.screen.stderr    = 1
                      log4perl.appender.screen.layout    = Log::Log4perl::Layout::SimpleLayout
                      );

Log::Log4perl::init( \$log4perl_conf );
my $logger = Log::Log4perl->get_logger($0);

if (defined($options->{'log'})) {
  my $layout = Log::Log4perl::Layout::SimpleLayout->new();
  my $appender = Log::Log4perl::Appender->new("Log::Log4perl::Appender::File", name => "logfile", filename => $options->{log});
  $appender->layout($layout);
  $logger->add_appender($appender);
  $logger->level($DEBUG);
}
$logger->info("logging initialized");

## main program
my($dfile, $server, $user, $no_deletes, $delete_url) = map {$options->{$_}} ('discrep_file', 'server', 'user', 'no_deletes', 'delete_url');

my $db_file = $options->{'db_mapping_file'};
open MAP, $db_file or die "Cannot open mapping file $db_file for reading: $!\n";
my $line = <MAP>;
chomp $line;
my ($db, $rest) = split(/,|\t/, $line, 2);
close MAP;

my $genes = &read_genes_to_delete($dfile);
my $ng = scalar(@$genes);
$logger->info("read $ng unique gene(s) to delete from $dfile");

# connect to server
my $dbi_str = "DBI:mysql:$db:$server";
print STDERR "dbi str = '$dbi_str'\n";
my $dbh = DBI->connect($dbi_str, $user, $password, {'RaiseError' => 1,'AutoCommit' => 0});
die("Could not connect to database ".DBI->errstr) if ($@);
my $get_transcript_sth = $dbh->prepare($TRANSCRIPT_ID_QUERY);
my $get_rna_sth = $dbh->prepare($RNA_ID_QUERY);

# map all the locus ids to internal transcript uniquenames
$cds_ids = {};
$genes_to_delete = [];
$rnas_to_delete = [];
foreach my $gene (@$genes) {
  	my($cds, $rna, $name, $lines, $hyp) = map {$gene->{$_}} ('cds', 'rna', 'name', 'lines', 'hypo');
  	$logger->debug("CDS $cds overlaps with RNA $rna");
  	# duplicates should have been eliminated already:
  	if (defined($cds_ids->{$cds})) {
  		$logger->logdie("internal error: found duplicate gene $cds");
  	}
    # Analyst will do manual review of 5S or 16S rRNA overlaps
    next if $name =~ /5S|16S/;

    if ($hyp or $name =~ /23S/) {
        # If gene is hypothetical or rRNA is 23S, the delete the gene and keep RNA
  		my $transcript_id = &get_transcript_id($get_transcript_sth, $cds);
  		if (!defined($transcript_id)) {
    		$logger->logdie("unable to map cds id $cds to transcript_id");
  		}
  		$logger->debug("mapped cds id $cds to $transcript_id");
  		$gene->{'transcript_id'} = $transcript_id;
  		push(@$genes_to_delete, $gene);
  		$logger->debug("DELETE CDS\t" . $gene->{'lines'}->[0]);
  	} else {
  		my $transcript_id = &get_transcript_id($get_rna_sth, $rna);
  		if (!defined($transcript_id)) {
    		$logger->logdie("unable to map RNA $rna to transcript_id");
  		}
  		$logger->debug("mapped RNA $rna to $transcript_id");
  		$gene->{'transcript_id'} = $transcript_id;
  		push(@$rnas_to_delete, $gene);
  		$logger->debug("DELETE RNA\t" . $gene->{'lines'}->[1]);
  	}
  	$cds_ids->{$cds} = 1;
}

# TODO - do the deletes, unless --no_deletes in effect
my $num_to_delete = scalar(@$genes_to_delete);
my $num_deleted = 0;
foreach my $gene (@$genes_to_delete) {
	  my($cds, $rna, $lines, $transcript_id) = map {$gene->{$_}} ('cds', 'rna', 'lines', 'transcript_id');
	  $logger->debug("deleting $cds / $transcript_id, which overlaps with RNA $rna");
	  if (!$no_deletes) {
	    	$num_deleted += &delete_gene($delete_url, $user, $password, $transcript_id, $db);
	  }
}
$logger->info("deleted $num_deleted/$num_to_delete gene(s) from $db that have bad CDS features");
my $cds_delete = $num_to_delete;

# TODO - do the deletes, unless --no_deletes in effect
$num_to_delete = scalar(@$rnas_to_delete);
$num_deleted = 0;
foreach my $rnas (@$rnas_to_delete) {
	  my($cds, $rna, $lines, $transcript_id) = map {$rnas->{$_}} ('cds', 'rna', 'lines', 'transcript_id');
	  $logger->debug("deleting $rna / $transcript_id, which overlaps with CDS $cds");
	  if (!$no_deletes) {
	    	$num_deleted += &delete_gene($delete_url, $user, $password, $transcript_id, $db);
	  }
}
$logger->info("deleted $num_deleted/$num_to_delete gene(s) from $db that have bad RNA features");
my $rna_delete = $num_to_delete;
my $total = $rna_delete + $cds_delete;

$logger->warn("number of deleted genes (" . $total . ") does not match number of expected unique genes to delete ($ng)") if ($total ne $ng);

$dbh->disconnect();
exit(0);

## subroutines

sub check_parameters {
  my $options = shift;

  ## make sure required parameters were passed
  my @required = qw(discrep_file db_mapping_file user server delete_url);
  for my $option ( @required ) {
    unless ( defined $options->{$option} ) {
      die("--$option is a required option");
    }
  }
    #Assign password to be read from the file if it exists.
    if (defined ($options->{'password_file'} && -s $options->{'password_file'} ) ) {
	open PASS, $options->{'password_file'} or die ("Cannot open password file ". $options->{'password_file'} ." : $!\n");
	print STDERR ("Password from file will take priority over --password option\n") if (defined $options->{'password'});
	$password= <PASS>;
	chomp $password;
	close PASS;
    } elsif (defined ($options->{'password'}) ){
	$password = $options->{'password'};
    } else {
	die("Neither a password or a password file were supplied.  Please supply one or the other");
    }

  ## defaults

  ## additional parameter checking?
}

sub read_genes_to_delete {
	my($dfile) = @_;
  	my $genes = [];
  	my $start_index = undef;
  	my $num_expected = undef;
  	my $num_skipped = undef;

  	# $pattern[0] = completely overlapping, $pattern[1] = partially overlapping
  	my $pattern = ['FATAL\:\sDiscRep_SUB:RNA_CDS_OVERLAP', 'DiscRep_SUB:RNA_CDS_OVERLAP'];

  	my $fh = FileHandle->new();
  	my $lnum = 0;
  	my $reading_list = 0;
  	$fh->open($dfile) || die "unable to read from $dfile";
  	while (my $line = <$fh>) {
		chomp($line);
    	++$lnum;
    	if (($line =~ /^$pattern->[0]/) || ($line =~ /^$pattern->[1]/)) {
        	$reading_list = 1;
        	if ($line =~ /RNA_CDS_OVERLAP\:\:(\d+) coding regions/) {
         		$num_expected = $1;
          		$num_skipped = 0;
          		$start_index = scalar(@$genes);
        	}
      	} elsif ($reading_list) {

        	# end of current block of problem genes
        	if ($line =~ /^\s*$/) {
     			$reading_list = 0;
          		# check count
          		my $num_found = scalar(@$genes) - $start_index;
          		die "read $num_found instead of the expected $num_expected problem genes at line $lnum of $dfile" if (($num_found + $num_skipped) != $num_expected);
          		$start_index = $num_expected = undef;
        	}
       		# still in a block, so the features should be listed in pairs
        	else {
          		my $next_line = <$fh>;
          		chomp($next_line);
          		++$lnum;
          		my $cds = undef;
          		my $rna = undef;
                my $rna_name = undef;
          		my $hypo = 0;

          		# gecDEC10A:CDS	ybl206 protein	lcl|gecDEC10A.contig.0:c350798-350679	ECDEC10A_0352
          		# gecDEC10A:rRNA	23S ribosomal RNA	lcl|gecDEC10A.contig.0:350459-353357	ECDEC10A_0351
          		if ($line =~ /^[^\:]+\:CDS.*/) {
            		my($type, $name, $coords, $locus) = split(/\t/, $line);
            		# Do not add genes to array unless they are conserved hypothetical proteins or hypothetical proteins
            		if (($name eq 'conserved hypothetical protein') || ($name eq 'hypothetical protein') || ($name =~ /^putative/)) {
            		    $hypo = 1;
          			}
            		die "failed to parse CDS id from $line" unless defined($locus);
            		$cds = $locus;
          		} else {
            		$logger->logdie("failed to parse CDS id from first line of CDS/rRNA pair at line " . ($lnum-1) . ": $line");
          		}
          		if ($next_line =~ /^[^\:]+\:\SRNA.*/) {
            		my($type, $name, $coords, $locus) = split(/\t/, $next_line);
            		die "failed to parse RNA id from $next_line" unless defined($locus);
            		$rna = $locus;
                    $rna_name = $name;
          		} else {
           			$logger->logdie("failed to parse RNA id from second line of CDS/rRNA pair at line $lnum: $next_line");
          		}
          		$logger->logdie("discrepancy file reports CDS $cds overlapping with RNA $rna with the same id") if ($cds eq $rna);
          		if (defined($cds_ids->{$cds})) {
            		$logger->warn("ignoring duplicate gene $cds");
            		++$num_skipped;
          		} else {
          			push(@$genes, {'cds' => $cds, 'rna' => $rna, 'name' => $rna_name, 'lines' => [$line, $next_line], 'hypo' => $hypo});
          			$cds_ids->{$cds} = 1;
          		}
      		}
    	}
  	}
  $fh->close();
  return $genes;
}

sub get_transcript_id {
  my($sth, $cds) = @_;
  $sth->execute($cds);
  my $r = $sth->fetchall_arrayref();
  return $r->[0]->[0];
}

sub delete_gene {
  my($delete_url, $user, $password, $transcript_id, $db) = @_;
  my $post_data =
    {
     'user' => $user,
     'password' => $password,
     'db' => $db,
     'id' => $transcript_id,
     # gcp = redirect to Gene Curation Page when done
     'flag' => 'gcp',
     'feature_deleter' => ''
    };

  my $ua = LWP::UserAgent->new();
  my $res = $ua->post($delete_url, $post_data);

  $logger->debug("POSTing to $delete_url to delete $transcript_id");
  if (!$res->is_success) {
    # we expect to get a 302 Moved response
    if ($res->status_line =~ /302 Moved/) {
      my $location = $res->header('Location');
      my $tidrw = $transcript_id;
      $tidrw =~ s/\./\\./;
      if ($location =~ /ORF_infopage\.cgi\?orf=$tidrw/) {
        $logger->info("deleted $transcript_id successfully");
        $logger->debug("delete response=" . $res->content());
        $logger->debug("302 Location=" . $location);
        return 1;
      }
    }
  }

  $logger->debug("delete status_line=" . $res->status_line());
  $logger->debug("delete response=" . $res->content());
  return 0;
}
