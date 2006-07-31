#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

prepare_ber_extended_nt_db.pl

=head1 SYNOPSIS

USAGE: prepare_ber_extended_nt_db.pl
           --input_list|-i
           --input_file|-f
           --input_dir|-D
           --input_extension|-x           .fsa
           --output_dir|-o
           --database|-d
           --server|-s                    SYBTIGR
           --username|-u                  access
           --password|-p                  access  
           --extend_by|-e                 300
           --help|-h

=head1 OPTIONS

B<--input_list,-i>
    Optional list of fasta sequence files containing polypeptide sequence ids.

B<--input_file,-f>
    Optional fasta sequence file containing polypeptide sequence ids.

B<--input_dir,-D>
    Optional directory containing fasta sequence files with polypeptide sequence ids.

B<--input_extension,-x>
    Optional extension used for file find when input_dir is specified.

B<--output_dir,-o>
    Output directory.

B<--extend_by,-e>
    Number of bases to extend the nucleotide sequences by (default 300).

B<--database,-d>
    Name of database from which to extract sequences.

B<--is_chado>
    Database schema is chado (1), or legacy db (0).

B<--server,-s>
    Server hosting the database.

B<--username,-u>
    Username to access database

B<--password,-p>
    Password to access database

B<--help,-h>
    This help documentation

=head1 DESCRIPTION

Prepares a database of extended NT sequences that correspond to
a list of query polypeptide sequences. The start and end positions
of the NT sequences are extended by the specified number of bases.

=head1 OUTPUT

The script creates two output files in the specified directory:

ber.cds.fsa
  FASTA sequence file containing extended NT sequences.
  
ber.mapping.list
  Mapping list that maps the query polypeptide IDs to the
  NT sequence IDs.

=head1 CONTACT

    Brett Whitty
    bwhitty@tigr.org

=cut

$|++;
use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Find;
use DBI;
use Data::Dumper;

my %opt;
GetOptions ( \%opt,
             'input_list|i=s',
             'input_file|f=s',
             'input_dir|D=s',
			 'input_extension|x=s',
             'output_dir|o=s',
			 'database|d=s',
             'server|s=s',
			 'is_chado=i',
             'username|u=s',
             'password|p=s',
			 'extend_by|e=i',
			 'help|h',
			 'man|m',
             ) || pod2usage();

if ($opt{'help'} || $opt{'man'}) {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDOUT} );
}
if (!$opt{'input_list'} && !$opt{'input_file'} && !$opt{'input_dir'}) {
	pod2usage("one of --input_list, --input_file, or --input_dir must be specified");
}
if (!$opt{'output_dir'}) {
	pod2usage("output directory must be specified using --output_dir flag");
}
if (!$opt{'database'}) {
	pod2usage("source database must be provided with --database flag");
}
if (!defined($opt{'is_chado'})) {
	pod2usage("must specify if database is chado or legacy db with --is_chado flag");
}
if (!$opt{'extend_by'}) {
	$opt{'extend_by'}=300;
}
if (!$opt{'server'}) {
	$opt{'server'} = 'SYBTIGR';
}
if (!$opt{'username'} || !$opt{'password'}) {
	$opt{'username'} = 'access';
	$opt{'password'} = 'access';
}
$opt{'output_dir'} =~ s/\/$//;

my @infiles = ();

if ($opt{'input_list'}) {
	unless (-e $opt{'input_list'}) {
		die "specified input list '$opt{input_list}' does not exist";
	}
	open (IN, $opt{'input_list'}) || die "couldn't open input list '$opt{input_list}' for reading";
	while (<IN>) {
		chomp;
		if (-e $_) {
			push (@infiles,$_);
		} else {
			print STDERR "ERROR: file '$_' specified in input list '$opt{input_list}' doesn't exist.";
		}
	}
	close IN;
}
if ($opt{'input_dir'}) {
	unless (-d $opt{'input_dir'}) {
		die "specified input dir '$opt{input_dir}' is not a directory";
	}
	if (!$opt{'input_suffix'}) {
		$opt{'input_suffix'} = '.fsa';
	}
    find({wanted => sub{ if(/$opt{input_suffix}$/) {push(@infiles,$_)}}, no_chdir => 1}, ($opt{'input_dir'}));
}
if ($opt{'input_file'}) {
		if (-e $opt{'input_file'}) {
			push (@infiles,$opt{'input_file'});
		} else {
			print STDERR "ERROR: specified input file '$opt{input_file}' doesn't exist.";
		}
}

my @sequence_ids;
foreach my $file (@infiles) {
	open (IN, $file) || die "couldn't open file '$file' for reading";
	while (<IN>) {
		chomp;
		if (/^>([^\s]+)/) {
			push (@sequence_ids, $1);
		}
	}
}

my $dbh = DBI->connect(
	"dbi:Sybase:server=$opt{server}; packetSize=8092", 
	$opt{username}, 
	$opt{password}, 
	{PrintError=>1, RaiseError=>1}
);

$dbh->do("use $opt{database}");


if ($opt{'is_chado'}) {
	
## get the maximum length of an assembly
my $max_feature_length = 0;
my $query = "SELECT max(f.seqlen), count(f.feature_id) " .
          "FROM feature f, cvterm c " .
          "WHERE f.type_id=c.cvterm_id " .
          "  AND c.name = 'assembly' " .
          "  AND f.is_analysis = 0 " .
          "  AND f.is_obsolete = 0 ";

my $feature_length_selector = $dbh->prepare( $query );
   $feature_length_selector->execute();
   
$max_feature_length = ( $feature_length_selector->fetchrow_array )[0] || 0;

$feature_length_selector->finish();

## set the textsize to support this max length
$dbh->do("set textsize $max_feature_length");

## quick length check
if (! $max_feature_length ) {
    die("no features of type 'assembly' were found with residue length > 0");
}

## select the CDS->polypeptide mapping table
$query = "select g.uniquename"
		 . " from feature_relationship r, feature f, feature g, cvterm c, cvterm v, cvterm t"
		 . " where r.subject_id = f.feature_id"
		 . " and r.object_id = g.feature_id"
	     . " and f.type_id = t.cvterm_id"
	     . " and g.type_id = c.cvterm_id"
	     . " and r.type_id = v.cvterm_id"
		 . " and c.name = 'CDS'"
	     . " and t.name = 'polypeptide'"
	     . " and v.name = 'derives_from'"
		 . " and f.uniquename = ?"
 		 . " AND f.is_analysis = 0" 
		 . " AND f.is_obsolete = 0";

my $find_cds = $dbh->prepare($query);

my @cds_names = ();
my $mapping_path = $opt{'output_dir'}."/ber.mapping.list";
open (OUT, ">$mapping_path") || die "couldn't write mapping file to '$mapping_path'";
foreach my $polypeptide (@sequence_ids) {
	$find_cds->execute($polypeptide);
	my $cds = ( $find_cds->fetchrow_array )[0];
	print OUT $polypeptide."\t".$cds."\n";
	push (@cds_names, $cds);
}
close OUT;
$find_cds->finish();

## query to locate the CDS on the assembly sequences
$query = "SELECT f.uniquename, l.srcfeature_id, l.fmin, l.fmax, l.strand"
		 . " FROM feature f, featureloc l"
		 . " WHERE f.feature_id=l.feature_id"
		 . " AND f.uniquename=?"
		 . " AND f.is_analysis = 0"
		 . " AND f.is_obsolete = 0";

my $cds_features = $dbh->prepare($query);

my $cds_loc = {};

foreach my $cds_name (@cds_names) {
	$cds_features->execute($cds_name);
	while ( my $row = $cds_features->fetchrow_hashref ) {
		$cds_loc->{$row->{'srcfeature_id'}}->{$row->{'uniquename'}}->{'fmin'} = $row->{'fmin'};
		$cds_loc->{$row->{'srcfeature_id'}}->{$row->{'uniquename'}}->{'fmax'} = $row->{'fmax'};
		$cds_loc->{$row->{'srcfeature_id'}}->{$row->{'uniquename'}}->{'strand'} = $row->{'strand'};
	}
}
$cds_features->finish();		 

$query = "SELECT residues " .
       	 "FROM feature f " .
         "WHERE feature_id = ?";

my $assembly_sequence = $dbh->prepare($query);

my $cds_path = "$opt{output_dir}/ber.cds.fsa";

open(OUT, ">$cds_path") || die "can't create $cds_path: $!";

foreach my $assembly_id(keys(%{$cds_loc})) {
	
	$assembly_sequence->execute($assembly_id);
	
	my $assembly = ( $assembly_sequence->fetchrow_array )[0];
	
	foreach my $cds(keys(%{$cds_loc->{$assembly_id}})) {
		my $start = $cds_loc->{$assembly_id}->{$cds}->{'fmin'};
		my $end = $cds_loc->{$assembly_id}->{$cds}->{'fmax'};
		
		$start = $start <= $opt{'extend_by'} 
				 ? 0 : $start - $opt{'extend_by'} - 1; 
		
		my $len = 1 + $end + $opt{'extend_by'} - $start;
        
		my $ntseq = substr($assembly, $start, $len);
		
		## reverse complement if necessary
        if ($cds_loc->{$assembly_id}->{$cds}->{'strand'} == -1) {
			reverseComplementDNA(\$ntseq);
		}
		
       	$ntseq =~ s/\W+//g;
        $ntseq =~ s/(.{1,60})/$1\n/g;
		my $header = "$cds";
                     
        print OUT ">$header\n$ntseq";
	}
}
close OUT;
$assembly_sequence->finish();

} else {

	## set textsize to maximum length of assemblies
	my $max_ln_txt;
	my $query = "SELECT MAX(length) FROM clone_info WHERE is_public = 1";
	my $sth = $dbh->prepare($query);
	$sth->execute();
	$sth->bind_columns(undef, \$max_ln_txt);
	$sth->fetch();
	$sth->finish();
	
	die "\n\nProblem fetching max sequence length from legacy db\n\n" unless defined $max_ln_txt;
	
	$max_ln_txt += 10;
	$dbh->do("SET TEXTSIZE $max_ln_txt");

	my $models = {};  
	
	my $get_models_qry = "SELECT a.end5, a.end3 FROM asm_feature a WHERE a.feat_name =  ?";
    my $get_models = $dbh->prepare($get_models_qry);

	foreach my $seq_id (@sequence_ids) {
		if ($seq_id =~ /^[^\.]+\.model\.(\d+)_(\d+)/) {
			my $model_name = "$1.m$2";
			my $asmbl_id = $1;
			my ($end5, $end3);
			$get_models->execute($model_name);
			$get_models->bind_columns(undef, \$end5, \$end3);
			$get_models->fetch();
			$models->{$asmbl_id}->{$model_name}->{'end5'} = $end5;
			$models->{$asmbl_id}->{$model_name}->{'end3'} = $end3;
			$models->{$asmbl_id}->{$model_name}->{'seq_id'} = $seq_id;
		} else {
			die "sequence id '$seq_id' doesn't match expected format for legacy db";
		}
	}
	$get_models->finish();

	my $get_asmbl_seq_qry = "SELECT sequence FROM assembly WHERE asmbl_id = ?";
	my $get_asmbl_seq = $dbh->prepare($get_asmbl_seq_qry);

	my $cds_path = "$opt{output_dir}/ber.cds.fsa";
	open(CDS, ">$cds_path") || die "can't create $cds_path: $!";

	my $mapping_path = $opt{'output_dir'}."/ber.mapping.list";
	open (MAPPING, ">$mapping_path") || die "couldn't write mapping file to '$mapping_path'";
	
	foreach my $asmbl_id (keys(%{$models})) {
		my $asmbl_seq;
		$get_asmbl_seq->execute($asmbl_id);
		$get_asmbl_seq->bind_columns(undef, \$asmbl_seq);
		$get_asmbl_seq->fetch();
		$get_asmbl_seq->finish();
		
		die "Sequence for assembly $asmbl_id is not in the database\n\n" unless defined $asmbl_seq;
		
		foreach my $model_name (keys(%{$models->{$asmbl_id}})){
			my $model = $models->{$asmbl_id}->{$model_name};
			
			my ($start, $end, $minusstrand) = $model->{'end5'} < $model->{'end3'} 
				? ($model->{'end5'}, $model->{'end3'}, 0) : ($model->{'end3'}, $model->{'end5'}, 1);
			
			$model_name =~ s/\.m/_/;
			
			my $header = "$opt{database}.model.$model_name";
			
			print MAPPING "$model->{seq_id}\t$header\n";
			
            $start = $start <= $opt{'extend_by'} ? 0 : $start - $opt{'extend_by'} - 1;
            my $ln = 1 + $end + $opt{'extend_by'} - $start;
			
            my $nucl = substr($asmbl_seq, $start, $ln);
			
			if ($minusstrand) {
				reverseComplementDNA(\$nucl);
			}
			
			$nucl =~ s/\W+//g;
			$nucl =~ s/(.{1,60})/$1\n/g;
			
			print CDS ">$header\n$nucl";
		}
	}
	close CDS;
	close MAPPING;
}


sub reverseComplementDNA {
        my ($r_seq) = @_;
        ${$r_seq} =~ tr/AaCcGgTtMmRrWwSsYyKkVvHhDdBb/TtGgCcAaKkYyWwSsRrMmBbDdHhVv/;
        ${$r_seq} = reverse(${$r_seq});
        return (length(${$r_seq}));
}
