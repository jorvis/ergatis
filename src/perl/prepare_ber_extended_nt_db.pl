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
           --bsml_source_list
           --bsml_source_file
           --bsml_source_dir
           --output_dir|-o
           --lookup_db
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

B<--temp_dir,-t>
    Directory for writing temporary files.
    
B<--extend_by,-e>
    Number of bases to extend the nucleotide sequences by (default 300).

B<--lookup_db>
    Where the id lookup file should be.

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

Also creates an id lookup for later stages (bsml creation).  This
link the predicted features back to the parent sequence. 

=head1 OUTPUT

The script creates two output files in the specified directory:

ber.nt.fsa
  FASTA sequence file containing extended NT sequences.
  
ber.mapping.list
  Mapping list that maps the query polypeptide IDs to the
  NT sequence IDs.

lookup_db
    A persistent hash created using the DB_File module. 
    Can be accessed using some polypeptide id with the
    value being the parent sequence id.

=head1 CONTACT

    Brett Whitty
    bwhitty@tigr.org

=cut

$|=1;
use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Workflow::Logger;
use DBI;
use File::Find;
use XML::Twig;
use DB_File;
use Data::Dumper;

my %opt;
GetOptions ( \%opt,
             'input_list|i=s',
             'input_file|f=s',
             'input_dir|D=s',
             'input_extension|x=s',
             'output_dir|o=s',
             'temp_dir|t=s',
             'database|d:s',
             'server|s:s',
             'is_chado:i',
             'username|u:s',
             'password|p:s',
             'bsml_source_list:s',
             'bsml_source_file:s',
             'bsml_source_dir:s',
             'extend_by|e=i',
             'lookup_db=s',
             'help|h',
             'man|m',
             ) || pod2usage();

if ($opt{'help'} || $opt{'man'}) {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDOUT} );
}

my $logfile = $opt{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$opt{'debug'});
$logger = $logger->get_logger();

my $bsml_input_flag = ($opt{'bsml_source_list'} || $opt{'bsml_source_file'} || $opt{'bsml_source_dir'})
                    ? 1 : 0;

if (!$opt{'input_list'} && !$opt{'input_file'} && !$opt{'input_dir'}) {
    pod2usage("one of --input_list, --input_file, or --input_dir must be specified");
}
if (!$opt{'output_dir'}) {
    pod2usage("output directory must be specified using --output_dir flag");
}
if (!$opt{'temp_dir'}) {
    pod2usage("output directory must be specified using --temp_dir flag");
}
if (!$opt{'database'} && !$bsml_input_flag) {
    pod2usage("source database must be provided with --database flag");
}
if (!defined($opt{'is_chado'}) && !$bsml_input_flag) {
    pod2usage("must specify if database is chado or legacy db with --is_chado flag");
}
if (!$opt{'extend_by'}) {
    $opt{'extend_by'}=300;
}
if (!$opt{'lookup_db'}) {
    pod2usage("must specifiy lookup_db output with --lookup_db");
}
if (!$opt{'server'}) {
    $opt{'server'} = 'SYBTIGR';
}
if (!$opt{'username'} || !$opt{'password'}) {
    $opt{'username'} = 'access';
    $opt{'password'} = 'access';
}
$opt{'output_dir'} =~ s/\/$//;
$opt{'temp_dir'} =~ s/\/$//;

my $cds_path = "$opt{output_dir}/ber.nt.fsa";
my $mapping_path = "$opt{output_dir}/ber.mapping.list";

## bsml-parsing-specific globals
my $bsml_cds_parent;
my $bsml_cds_locs;
my $bsml_cds_strand;
my $bsml_sequences;
my $bsml_mapping;
## ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

## Used to make id lookup table
my %idLookup;
tie(%idLookup, 'DB_File', $opt{'lookup_db'}) or
    $logger->logdie("Could not open $opt{lookup_db} ($!)");
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^

my @infiles = get_input_files();

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

my @bsmlfiles = ();
if ($bsml_input_flag) {
    ## do bsml stuff    
    my $cds_fh;
    my $mapping_fh;
    
    my @bsmlfiles = get_bsml_sources();
    
    open ($cds_fh,    ">$cds_path")     || die "couldn't open cds fsa file '$cds_path' for writing: $!";
    open ($mapping_fh, ">$mapping_path") || die "couldn't open mapping file '$mapping_path' for writing: $!";
    
    foreach my $bsmlfile(@bsmlfiles) {
        do_parse_bsml($bsmlfile);
    }
    do_extract_bsml_sequences($cds_fh, $mapping_fh, @sequence_ids);
    
    close $cds_fh;
    close $mapping_fh;
    
} else {
    #do database stuff
    my $dbh = DBI->connect(
        "dbi:Sybase:server=$opt{server}; packetSize=8092", 
        $opt{username}, 
        $opt{password}, 
        {PrintError=>1, RaiseError=>1}
    );
    
    $dbh->do("use $opt{database}");

    if ($opt{'is_chado'}) {
        do_extract_from_chado($dbh, @sequence_ids);
    } else {
        do_extract_from_legacy($dbh, @sequence_ids);
    }
}

untie %idLookup;

exit();

sub do_extract_from_chado {

    my ($dbh, @sequence_ids) = @_;
    
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
   
    $max_feature_length = ($feature_length_selector->fetchrow_array)[0] || 0;

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
            $idLookup{$row->{'uniquename'}} = $row->{'srcfeature_id'};
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
            $ntseq =~ tr/a-z/A-Z/;
            my $header = "$cds";
                     
            print OUT ">$header\n$ntseq";
        }
    }
    close OUT;
    $assembly_sequence->finish();
}

sub do_extract_from_legacy {
    my ($dbh, @sequence_ids) = @_;

    ## set textsize to maximum length of assemblies
    my $max_ln_txt;
    my $query = "SELECT MAX(DATALENGTH(sequence)) FROM assembly";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    $sth->bind_columns(undef, \$max_ln_txt);
    $sth->fetch();
    $sth->finish();

    #$max_ln_txt = 30000000;
    
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

            #Fill the idLookup database
            $idLookup{$seq_id} = "$1.assembly.$2" if($seq_id =~ /^([^\.]+)\.model\.(\d+)_\d+/);

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

    open(CDS, ">$cds_path") || die "can't create $cds_path: $!";

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

            if (length($nucl) == 0) {
                die("'$header' had zero length sequence");
            }
            
            $nucl =~ s/\W+//g;
            $nucl =~ s/(.{1,60})/$1\n/g;
            $nucl =~ tr/a-z/A-Z/;
            
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

sub do_parse_bsml {
    my ($bsmlfile) = @_;

    print STDERR "parsing '$bsmlfile'\n";
    parse_bsml_sequences($bsmlfile);
    parse_bsml_interval_locs($bsmlfile);
#print Dumper $bsml_cds_parent;
#print Dumper $bsml_cds_locs;
#print Dumper $bsml_cds_strand;
#print Dumper $bsml_sequences;
#print Dumper $bsml_mapping;
}

sub do_extract_bsml_sequences {
    my ($cds_fh, $mapping_fh, @sequence_ids) = @_;

    my $extract_hash;
    
    foreach my $pep_id(@sequence_ids) {
        if (!defined($bsml_mapping->{$pep_id})) {
            die "no mapping for $pep_id to CDS";
        }
        $extract_hash->{$bsml_cds_parent->{$bsml_mapping->{$pep_id}}}->{$bsml_mapping->{$pep_id}} = 1;
        print $mapping_fh "$pep_id\t$bsml_mapping->{$pep_id}\n";
    }

    foreach my $asmbl_id(keys(%{$extract_hash})) {
        my $asmbl_seq = get_sequence_by_id(
                                            $bsml_sequences->{$asmbl_id}->{'source'}, 
                                            $bsml_sequences->{$asmbl_id}->{'id'}
                                          );
        
        if (length($asmbl_seq) == 0) {
            die "couldn't fetch sequence for assembly '$asmbl_id' from '$bsml_sequences->{$asmbl_id}'";
        }
        
        foreach my $cds_id(keys(%{$extract_hash->{$asmbl_id}})) {
            $idLookup{$cds_id} = $asmbl_id;
            my $extended_cds = get_extended_sequence(
                                                        $asmbl_seq, 
                                                        $bsml_cds_locs->{$cds_id}->[0], 
                                                        $bsml_cds_locs->{$cds_id}->[1], 
                                                        $opt{'extend_by'},
                                                    );
            if ($bsml_cds_strand->{$cds_id}) {
                reverseComplementDNA(\$extended_cds);
            }
            write_seq_to_fasta_fh($cds_fh, $cds_id, $extended_cds);     
        }
    }
}


sub parse_bsml_sequences {
        my ($file) = @_;
        
        my $ifh;
        if (-e $file.".gz") {
            $file .= ".gz";
        } elsif (-e $file.".gzip") {
            $file .= ".gzip";
        }

        if ($file =~ /\.(gz|gzip)$/) {
            open ($ifh, "<:gzip", $file);
        } else {
            open ($ifh, "<$file");
        }
        
        my $twig = new XML::Twig(   TwigRoots => {'Sequence' => 1, 'Feature-group' => 1},
                                    TwigHandlers => {'Sequence' => \&process_sequence}
                                );
        
        $twig->parse($ifh);
        close $ifh;
}

sub process_sequence {
    my ($twig, $sequence) = @_;

    my $seq_id;
    my $class;
    my $molecule;
    my $seq_data_import;
    my $seq_data;
    my $source;
    my $identifier;
    my $format;
    
    $seq_id = $sequence->{'att'}->{'id'};
    $class = $sequence->{'att'}->{'class'};
    $molecule = $sequence->{'att'}->{'molecule'};

    if ($seq_data_import = $sequence->first_child('Seq-data-import')) {
        
        $source = $seq_data_import->{'att'}->{'source'};
        $identifier = $seq_data_import->{'att'}->{'identifier'};
        $format = $seq_data_import->{'att'}->{'format'};
        
        unless (-e $source) {
            die "fasta file referenced in BSML Seq-data-import '$source' doesn't exist";
        }
        unless (defined($identifier)) {
            die "Seq-data-import for '$seq_id' does not have a value for identifier";
        }
        
        $bsml_sequences->{$seq_id}->{'source'} = $source;
        $bsml_sequences->{$seq_id}->{'id'} = $identifier;
        
    } elsif ($seq_data = $sequence->first_child('Seq-data')) {
        ## sequence is in the BSML
        ## so it will be written to a fasta file in the output dir
        
        my $sequence_file = $opt{'temp_dir'}."/$seq_id.fsa";
    
        write_seq_to_fasta($sequence_file, $seq_id, $seq_data->text());
        
        $bsml_sequences->{$seq_id}->{'source'} = $sequence_file;
        $bsml_sequences->{$seq_id}->{'id'} = $seq_id;
        
    } else {
        ## there is no Seq-data or Seq-data-import for the sequence
        die "No sequence present in BSML sequence element";
    }
  
    if ($sequence->first_child('Feature-tables')) {
        foreach my $child ($sequence->first_child('Feature-tables')->children('Feature-group')) {
            foreach my $grandchild($child->children('Feature-group-member')) {
                if ($grandchild->{'att'}->{'feature-type'} eq 'CDS') {
                    $bsml_cds_parent->{$grandchild->{'att'}->{'featref'}} = $seq_id;
                }
            }
        }
    }
        
    $twig->purge;
}

sub write_seq_to_fasta {
    my ($file, $header, $sequence) = @_;
    
    open (my $fh, ">$file") || die "couldn't write fasta file '$file'";
        
    write_seq_to_fasta_fh($fh, $header, $sequence); 
    
    close $fh;
}

sub write_seq_to_fasta_fh {
    my ($file_handle, $header, $sequence) = @_;
    
    $sequence =~ s/\W+//g;
    $sequence =~ s/(.{1,60})/$1\n/g;
    $sequence =~ tr/a-z/A-Z/;
        
    print $file_handle ">$header\n$sequence";
}

sub get_sequence_by_id {
    my ($fname, $id) = @_;
    my $seq_id = '';
    my $sequence = '';
    open (IN, $fname) || die "couldn't open fasta file '$fname' for reading";
    TOP: while (<IN>) {
        chomp;
        if (/^>([^\s]+)/) {
            $seq_id = $1;
            if ($seq_id eq $id) {
                while (<IN>) {
                    chomp;
                    if (/^>/) {
                        last TOP;
                    } else {
                        $sequence .= $_;
                    }
                }
            }   
        }
    }
    close IN;

    return $sequence;
}

sub parse_bsml_interval_locs {
        my ($file) = @_;
        
        my $ifh;

        if (-e $file.".gz") {
            $file .= ".gz";
        } elsif (-e $file.".gzip") {
            $file .= ".gzip";
        }

        if ($file =~ /\.(gz|gzip)$/) {
            open ($ifh, "<:gzip", $file);
        } else {
            open ($ifh, "<$file");
        }

        my $twig = new XML::Twig(
                                 twig_roots => {
                                        'Feature'            => \&process_feat,
                                        'Feature-group'      => \&process_feat_group,
                                               }
                                );
        $twig->parse($ifh);
}

sub process_feat {
    my ($twig, $feat) = @_;
    my $id = $feat->{'att'}->{'id'};

    if ($feat->{'att'}->{'class'} eq 'CDS') {
        my $seq_int = $feat->first_child('Interval-loc');
        my $complement = $seq_int->{'att'}->{'complement'};
        $bsml_cds_strand->{$id} = $complement;
        $bsml_cds_locs->{$id} = [$seq_int->{'att'}->{'startpos'}, $seq_int->{'att'}->{'endpos'}];
    }

    $twig->purge;
}

sub process_feat_group {
    my ($twig, $feat_group) = @_;

    my $cds_id = '';
    my $polypeptide_id = '';
    
    my $feat_group_id = $feat_group->{'att'}->{'group-set'};
    
    foreach my $child ($feat_group->children('Feature-group-member')) {
        if ($child->{'att'}->{'feature-type'} eq 'CDS') {
            if ($cds_id ne '') {
                die "Feature-group '$feat_group_id' contains more than one cds feature";
            }
            $cds_id = $child->{'att'}->{'featref'};
        } elsif ($child->{'att'}->{'feature-type'} eq 'polypeptide') {
            if ($polypeptide_id ne '') {
                die "Feature-group '$feat_group_id' contains more than one polypeptide feature";
            }
            $polypeptide_id = $child->{'att'}->{'featref'};
        }
    }
#   print STDERR "$polypeptide_id $cds_id\n";
    if ($cds_id ne '' && $polypeptide_id ne '') {
        $bsml_mapping->{$polypeptide_id} = $cds_id;
    }
    
    $twig->purge;
}

sub get_bsml_sources {
    my @bsmlfiles;
    
    if ($opt{'bsml_source_list'}) {
        unless (-e $opt{'bsml_source_list'}) {
            die "specified bsml source list '$opt{bsml_source_list}' does not exist";
        }
        open (IN, $opt{'bsml_source_list'}) || die "couldn't open bsml source list '$opt{bsml_source_list}' for reading";
        while (<IN>) {
            chomp;
            if (-e $_) {
                push (@bsmlfiles,$_);
            } else {
                print STDERR "ERROR: file '$_' specified in input list '$opt{bsml_source_list}' doesn't exist.";
            }
        }
        close IN;
    }
    if ($opt{'bsml_source_dir'}) {
        unless (-d $opt{'bsml_source_dir'}) {
            die "specified bsml source dir '$opt{bsml_source_dir}' is not a directory";
        }
        find({wanted => sub{ if(/\.bsml$|\.bsml\.gz$/) {push(@infiles,$_)}}, no_chdir => 1}, ($opt{'bsml_source_dir'}));
    }
    if ($opt{'bsml_source_file'}) {
        if (-e $opt{'bsml_source_file'}) {
            push (@bsmlfiles,$opt{'bsml_source_file'});
        } else {
            print STDERR "ERROR: specified input file '$opt{bsml_source_file}' doesn't exist.";
        }
    }
    
    return @bsmlfiles;  
}

## get extended nt sequence
sub get_extended_sequence {
        my ($sequence, $startpos, $endpos, $extend) = @_;

        my $subseq = '';

        $startpos -= $extend;
        $endpos += $extend;
        
        if ($startpos > $endpos) {
            ## at this point we should have startpos < endpos
            ## if not, this should be resolved at the bsml parsing stage
            $logger->logdie("BSML contains an exon with startpos > endpos");
        }

        my $len = $endpos - $startpos;

        if ($startpos < 0) {
        ## for circular molecules
            $subseq .= substr($sequence, $startpos, abs($startpos));
            $subseq .= substr($sequence, 0, $endpos);
        } else {
        ## otherwise
            $subseq .= substr($sequence, $startpos, $len);
        }

        return $subseq;
}

#sub get_sequence {
#        my ($fname) = @_;
#
#        my $sequence = '';
#        my $flag = 0;
#
#        open (IN, $fname) || $logger->logdie("couldn't open fasta file '$fname' for reading");
#        while (<IN>) {
#                chomp;
#                if (/^>/) {
#                        if ($flag) {
#                                $logger->logdie("Unexpectedly encountered more than one fasta record in input nt sequence file");
#                        }
#                        $flag = 1;
#                        next;
#                } else {
#                        $sequence .= $_;
#                }
#        }
#        $sequence =~ s/\W+//g;
#        return $sequence;
#}

#sub get_sequence_id {
#        my ($fname) = @_;
#        my $id = '';
#        open (IN, $fname) || die "couldn't open fasta file '$fname' for reading";
#        while (<IN>) {
#                chomp;
#                if (/^>([^\s]+)/) {
#                        $id = $1;
#                        last;
#                }
#        }
#        close IN;
#
#        return $id;
#}

sub get_input_files {
    my @infiles;

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

    return @infiles;
}
