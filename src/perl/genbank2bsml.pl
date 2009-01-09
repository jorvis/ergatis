#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use File::Basename;
use FileHandle;
use Bio::SeqIO;
use Bio::Tools::CodonTable; # for creating CDSs from mRNAs, bug #3300
use BSML::BsmlBuilder;
use Ergatis::IdGenerator;
use Data::Dumper;
umask(0000);

# optional mapping parsed from --organism_to_prefix_mapping argument
my $organism_to_prefix_map = undef;
my $current_prefix = undef;

use Data::Dumper;

my %options = &parse_options();
# output options
my $ifile = $options{input_file};
my $odir = $options{output_dir};
my $ofile = $options{output_bsml};
my %fsa_files; # file names and handlers and for each fasta output file
my $idcreator;  # Ergatis::IdGenerator object
my $feature_id_lookup = {};

my %TagCount; #NOTE: rename this or something
print "Parsing $ifile\n";
my %gbr = %{parse_genbank_file($ifile)};
my $doc = new BSML::BsmlBuilder();
print "Converting $ifile to bsml \n";
&to_bsml(\%gbr, $doc);

#output the BSML
#my $bsmlfile = "$odir/".$gbr{'gi'}.".bsml";
print "Outputting bsml to $ofile\n";
$doc->createAndAddAnalysis(
                            id => $options{'analysis_id'},
                            sourcename => $options{'output_bsml'},
                            algorithm => 'genbank2bsml',
                            program => 'genbank2bsml',
                            version => 'current',
                          );
# gzip output
$doc->write($ofile, '', 1);
# $doc->write($ofile);
chmod (0777, $ofile);

# close file handlers
foreach my $class (keys %fsa_files) {
  close($fsa_files{$class}{fh});  
}

%gbr = ();

#tag count info
print "\nTag Counts:\n";
foreach (sort { $TagCount{$b} <=> $TagCount{$a} } keys %TagCount) {
    print "  $_: ".$TagCount{$_}."\n";
}

#
# Subroutines follow
#
sub print_usage {
        my $byebye = shift;
        my $progname = $0;
        die << "END";
$byebye
usage: $progname --input_file|i <file> --output_bsml|b <file> --output_dir|o <dir>
END
}


sub parse_options {
    my %options = ();
    GetOptions( \%options,
        'input_file|i=s',
        'output_dir|o=s',
        'output_bsml|b=s',  #output bsml file
        'id_repository=s',
        'project=s',
        'organism_to_prefix_mapping=s',  #tab-delimited file that maps organism name to id prefix
        'generate_new_seq_ids=s',        #whether to replace existing sequence ids with ergatis ones
        'analysis_id|a=s',
	'skip_unknown_dbxref=s',
        ) || &print_usage("Unprocessable option");

    # check for required parameters

    unless($options{'analysis_id'}) {
        $options{'analysis_id'} = 'genbank2bsml_analysis';
    }
    (defined($options{input_file})) || die "input_file is required parameter";
    # we may be passed a file name, but the .gz is what is actually there
    unless (-e $options{input_file}) {
    $options{input_file} .= ".gz";
    }
    (-r $options{input_file}) || die "input_file ($options{input_file}) not readable: $!";

    (-w $options{output_dir} && -d $options{output_dir}) 
    || die "output_dir directory ($options{output_dir}) not writeable: $!";

    (defined($options{output_bsml})) || die "output_bsml is a required option";

    if (defined($options{organism_to_prefix_mapping})) {
	(-r $options{organism_to_prefix_mapping}) 
	    || die "organism_to_prefix_mapping ($options{organism_to_prefix_mapping}) not readable: $!";
	$organism_to_prefix_map = &parse_organism_to_prefix_mapping($options{organism_to_prefix_mapping});
    }

    # default of skip_unknown_dbxref is not to skip
    unless (defined($options{output_bsml})) {
      $options{skip_unknown_dbxref} = 0;
    }

    ## Now set up the id generator stuff
    if ($options{'id_repository'}) {

        # we're going to generate ids
        $idcreator = new Ergatis::IdGenerator('id_repository' => $options{'id_repository'});
        $idcreator->set_pool_size('exon'        => 30,
                                  'polypeptide' => 30,
                                  'transcript'  => 30,
                                  'CDS'         => 30,
                                  'gene'        => 30);
    } else {
        die "--id_repository is a required option";
    }

    print "Executing $0 with options\n";
    foreach (keys %options) { print "  $_: $options{$_}\n";}

    return %options;
}

# in: genbank/refseq file
#out: hashref of all the info
sub parse_genbank_file {
    my $gb_file = shift;
    my %gbr;

    # support for compressed input
    # example taken from http://bioperl.org/wiki/HOWTO:SeqIO
    my $seq_obj;
    if ($gb_file =~ /\.(gz|gzip)$/) {
    $seq_obj = Bio::SeqIO->new(-file   => "gunzip -c $gb_file |",
                   -format => 'genbank');
    }
    else {   
    $seq_obj = Bio::SeqIO->new(-file   => "<$gb_file",
                   -format => 'genbank');
    }

    # force die on multi-sequence genbank file
    # see http://jorvis-lx:8080/bugzilla/show_bug.cgi?id=3316
    my $n_seqs = 0;
    while (my $seq = $seq_obj->next_seq) {
    ++$n_seqs;
    die "Unsupported: multiple records ($n_seqs) in $gb_file." if ($n_seqs > 1);
    $gbr{'component'} = 'chromosome'; #default;
    $gbr{'strain'} = '';

    $gbr{'accession'} = $seq->accession_number;
    $gbr{'taxon_id'} = $seq->species->ncbi_taxid();
    $gbr{'gi'} = $seq->primary_id();
    defined($gbr{'gi'}) || die "No gi in $gb_file";

    #first word is genus, all follows is species (a workaround to encode the entire scientific name in BSML/chado)
    #this pulls the ORGANISM field up to the taxonomic lineage
    $gbr{'organism'} = $seq->species->scientific_name;

    if (defined($organism_to_prefix_map)) {
	$current_prefix = $organism_to_prefix_map->{$gbr{'organism'}};
	if (!defined($current_prefix) || ($current_prefix =~ /^\s*$/)) {
	    # TODO - this should result in a warning, but logging doesn't appear to be enabled
	    warn "No prefix defined in $options{organism_to_prefix_mapping} for organism scientific name ($gbr{organism}) "; #AARON
	}
    }

    $gbr{'organism'} =~ m/(\S+)\s+(.*)/;
    $gbr{'genus'} = $1;
    $gbr{'species'} = $2;

    ( ($gbr{'genus'}) && ($gbr{'species'}) ) || die "Unable to parse genus/species from organism ($gbr{organism})";

    $gbr{'sequence'} = $seq->seq();
    $gbr{'seq_len'} = $seq->length();

    # $gbr{'molecule'} = $seq->molecule();
    # Placed in Sequence@molecule
    # valid values are: mol-not-set|dna|rna|aa|na|other-mol
    if ($seq->molecule() =~ /dna/i) {
        $gbr{'molecule'} = 'dna';
    }
    elsif ($seq->molecule() =~ /rna/i) {
        $gbr{'molecule'} = 'rna';
    }
    elsif ($seq->molecule() =~ /aa/i) {
        die "molecule=aa; dealing with an amino acid record is so vastly unsupported script is killing itself to preserve its dignity."
    }
    else {
        die "Unknown molecule: ".$seq->molecule();
    }

    #Making molecule name unique yet descriptive
    $gbr{'molecule_name'} = "$gbr{'organism'} $gbr{'accession'}";

    # Get expanded "polymer_type" as entered in ard.obo
    # default is $seq->molecule()
    # change it if we know by lineage it's something specific (some older records would just have RNA)
    # cross reference this against $gbr{molecule} for consistency
    # see bug #3251: http://jorvis-lx:8080/bugzilla/show_bug.cgi?id=3251
    # taxonomy names taken from http://www.ncbi.nlm.nih.gov/genomes/VIRUSES/vifam.html
    # polymer types adapted from genbank list: http://www.bio.net/bionet/mm/genbankb/2001-October/000103.html
    # 34-36      Blank, ss- (single-stranded), ds- (double-stranded), or
    #            ms- (mixed-stranded)
    # 37-42      Blank, DNA, RNA, tRNA (transfer RNA), rRNA (ribosomal RNA), 
    #            mRNA (messenger RNA), uRNA (small nuclear RNA), snRNA
    $gbr{'polymer_type'} = $seq->molecule(); #default
    my $classification = join(" ",$seq->species->classification);
    if ($classification =~ /dsDNA viruses/) {
        # dsDNA = DNA
        $gbr{'polymer_type'} = "DNA";
    }
    elsif ($classification =~ /dsRNA viruses/) {
        $gbr{'polymer_type'} = "ds-RNA";
    }
    elsif ($classification =~ /ssDNA viruses/) {
        $gbr{'polymer_type'} = "ss-DNA";
    }
    elsif ($classification =~ /ssRNA negative-strand viruses/) {
        $gbr{'polymer_type'} = "ss-RNA-";
    }
    elsif ($classification =~ /ssRNA positive-strand viruses/) {
        $gbr{'polymer_type'} = "ss-RNA+";
    }
    # these shouldn't be case-insensitive since they pop up from time to time in words
    # taking this out because most of the time it fails it's nothing
#   elsif ($classification =~ /RNA/) {
#       die "Unable to set polymer_type, unknown RNA in classification ($classification)";
#   }
#   elsif ($classification =~ /DNA/) {
#       die "Unable to set polymer_type, unknown DNA in classification ($classification)";
#   }

    # check for potential misparse
    # changing these to warnings, and having taxonomic classification take priority
    # see bug #3251 http://jorvis-lx:8080/bugzilla/show_bug.cgi?id=3251
    if ( ($gbr{'molecule'} eq 'dna') && !($gbr{'polymer_type'} =~ /DNA/) ) {
        warn "Mismatch between molecule ($gbr{molecule}, parsed from ".$seq->molecule().") and polymer_type ($gbr{polymer_type})";
        if ($gbr{'polymer_type'} =~ /RNA/) {
        $gbr{'molecule'} = 'rna';
        }
        else {
        die "This is an unreconcialable mismatch";
        }
    }
    if ( ($gbr{'molecule'} eq 'rna') && !($gbr{'polymer_type'} =~ /RNA/) ) {
        warn "Mismatch between molecule ($gbr{molecule}, parsed from ".$seq->molecule().") and polymer_type ($gbr{polymer_type})";
        if ($gbr{'polymer_type'} =~ /DNA/) {
        $gbr{'molecule'} = 'dna';
        }
        else {
        die "This is an unreconcialable mismatch";
        }
    }

    # Add strand attribute (this may not be used)
    # direct creation supported via createAndAddExtendedSequence, but this is simpler
    # valid values are: std-not-set|ss|ds|mixed|std-other
    if ($gbr{'polymer_type'} eq 'DNA') {
        $gbr{'strand'} = "ds";
    }
    elsif ($gbr{'polymer_type'} =~ /ss/ ) {
        $gbr{'strand'} = "ss";
    }
    elsif ($gbr{'polymer_type'} =~ /ds/ ) {
        $gbr{'strand'} = "ds";      
    }

    # correct polymer_type (this could be earlier, but then
    # have to change strand setting)
    if ($gbr{polymer_type} eq 'ss-RNA') {
      $gbr{polymer_type} = 'ss_RNA_viral_sequence';
    } elsif ($gbr{polymer_type} eq 'ss-RNA+') {
      $gbr{polymer_type} = 'postive_sense_ssRNA_viral_sequence';
    } elsif ($gbr{polymer_type} eq 'ss-RNA-') {
      $gbr{polymer_type} = 'negative_sense_ssRNA_viral_sequence';
    } elsif ($gbr{polymer_type} eq 'ds-RNA') {
      $gbr{polymer_type} = 'ds_RNA_viral_sequence';
    } elsif ($gbr{polymer_type} eq 'ss-DNA') {
      $gbr{polymer_type} = 'single_stranded_DNA_chromosome';
    }

    #currently unused
    #$gbr{'sub_species'} = $seq->species->sub_species();
    #$gbr{'variant'} = $seq->species->variant();    
    #$gbr{'genus'} = $seq->species->genus();   #changed to organism parsing
    #$gbr{'species'} = $seq->species->species(); #changed to organism parsing
    #$gbr{'strain'} = $seq->species->strain();
    #$gbr{'desc'} = $seq->desc();
    #$gbr{'organism'} = $seq->species->common_name;

    #$seq->moltype() will say dna if the sequence is dna, this is a problem for genbank records
    #where an rna sequence will be represented as dna
    
    my $anno_collection = $seq->annotation;
    my @annotations = $anno_collection->get_Annotations(); #or Annotations('comment')

    #search FEATURE.source for:
    #  -/plasmid
    #  -/strain
    my $ugc = 0; #unknown group count
    my %FeatureCount = ();
    my %transl_tables; #track of translation tables in CDSs
    for my $feat_object ($seq->get_SeqFeatures()) {
        my $primary_tag = $feat_object->primary_tag;
        ++$FeatureCount{$primary_tag};
        
        if ($primary_tag eq 'source') {
        for my $tag ($feat_object->get_all_tags) {
            if ($tag eq 'plasmid') {
            $gbr{'component'} = "plasmid";
            }
            elsif ($tag eq 'strain') {
            $gbr{'strain'} = join(' ', $feat_object->get_tag_values($tag));
            }
            else {
            # record unknown tag
            ++$TagCount{"unknown_source_tag_".$tag};
            }
        }
        }
        else { #!= source but is something we care about

        #this test junk is examing the types of inexact locations
        #currently counting it on all features
        my $loc_thing = $feat_object->location->location_type()."_".
            $feat_object->location->start_pos_type."_".
            $feat_object->location->end_pos_type;
        
        ++$TagCount{$loc_thing};
        unless ($loc_thing =~ /^EXACT/) {
            print "\n$loc_thing: ".$gbr{'accession'};
        }
        #/test junk
        
        # random testing
        die "Hey, there are transcript features" if $primary_tag eq 'transcript';
        #/testing


        if ($primary_tag eq 'CDS'  ||
            $primary_tag eq 'gene' || 
            #$primary_tag eq 'misc_feature' ||  #not supported by chado (not in cvterm)
            $primary_tag eq 'promoter' ||
            $primary_tag eq 'exon' ||
            #$primary_tag eq 'variation' ||  #this is an example of an IN-BETWEEN tag
            $primary_tag eq 'intron' ||
            $primary_tag eq 'mRNA' ||
	    $primary_tag eq 'tRNA' || # SO:0000253
	    $primary_tag eq 'rRNA'    # SO:0000252
            ) {
            
            # obtain feature group
            my $feature_group = '';
            my $feature_tag = '';
            my $feature_value = '';
            # could be multiples of same gene, so look for locus_tag first
            if ($feat_object->has_tag("locus_tag")) {
            $feature_tag = 'locus_tag';
            $feature_value = join('_',$feat_object->get_tag_values("locus_tag"));
            }
            elsif ($feat_object->has_tag("gene")) {
            $feature_tag = 'gene';
            $feature_value = join('_',$feat_object->get_tag_values("gene"));
            }
            elsif ($feat_object->has_tag("protein_id")) {
            $feature_tag = 'protein_id';
            $feature_value = join('_',$feat_object->get_tag_values("protein_id"));
            }
            else {
            #no way to group feature, so put it in its own feature_group
            $feature_tag = 'unknown';
            $feature_value = $ugc;              
            ++$ugc;
            }
            $feature_group = $feature_tag."_".$feature_value;
            
            # id is unique across all bsml documents (use accession instead?  can't use taxon id)
            # //Feature/@id has to begin with alpha character?
            # old way: my $feature_id = "gi_".$gbr{'gi'}.".".$feature_group.".".$primary_tag."_".$FeatureCount{$primary_tag};
            my $feature_id = "gi_".$gbr{'gi'}.".".$primary_tag."_".$FeatureCount{$primary_tag};
            $feature_id =~ s/\///; #/'s in this will cause trouble down the road
            if (exists $gbr{'Features'}->{$feature_group}->{$feature_id}) { die "feature id is not unique!  feature_id: $feature_id"; }
            
            my $is_complement;
            if ($feat_object->strand == 1) {
            $is_complement = 0;
            }
            elsif ($feat_object->strand == -1) { 
            $is_complement = 1;
            }
            else {
            die "unknown feature strand in $feature_group $feature_id: (".$feat_object->strand.")"; 
            }

            #store information on the feature in the hashref
            $gbr{'Features'}->{$feature_group}->{$feature_id} = {
                                      class => $primary_tag,
                                      is_complement => $is_complement,
                                      start => $feat_object->start,
                                      end => $feat_object->end,
                                      location_type => $feat_object->location->location_type,
                                      start_type => $feat_object->location->start_pos_type,
                                      end_type => $feat_object->location->end_pos_type,
                                      feature_tag => $feature_tag,
                                      feature_value => $feature_value
                                    };

            # store locations as an array because there might be more than one of them and how we deal with them
            # will be determined by the class and the other members of the feature_group
            my @locations = $feat_object->location->each_Location();
            foreach my $loc (@locations) {
            #print $loc->start."\t".$loc->start_pos_type."\t".$loc->end."\t".$loc->end_pos_type."\n";
            push(@{$gbr{'Features'}->{$feature_group}->{$feature_id}->{locations}},
                 {start => $loc->start, start_pos_type => $loc->start_pos_type, 
                  end => $loc->end, end_pos_type => $loc->end_pos_type});
            #print $gbr{'Features'}{$feature_group}{$feature_id}{locations}[0]{start}."\n";
            }

            # make a fake translation for mRNAs later (see bug #3300 http://jorvis-lx:8080/bugzilla/show_bug.cgi?id=3300)
            if ($primary_tag eq 'CDS' || $primary_tag eq 'mRNA') {
	      #store (possibly joined) nucleotide sequence using spliced_seq()
	      $gbr{'Features'}->{$feature_group}->{$feature_id}->{'spliced_seq'} = $feat_object->spliced_seq->seq();
	      $gbr{'Features'}->{$feature_group}->{$feature_id}->{'feature_len'} = $feat_object->spliced_seq->length();
	      $gbr{'Features'}->{$feature_group}->{$feature_id}->{'translation'} = 0; #default null value
            }

            #See bug 2414
            #Can assume all CDSs have transl_table?  No.
            #Will a document have uniform transl_tables?  Dunno.
            #See http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c#SG1
            #  By default all transl_table in GenBank flatfiles are equal to id 1, and this is not shown. 
            #  When transl_table is not equal to id 1, it is shown as a qualifier on the CDS feature.
            if ($feat_object->has_tag("transl_table")) {                
                ++$TagCount{join('',$feat_object->get_tag_values("transl_table"))};             
                ++$transl_tables{join('',$feat_object->get_tag_values("transl_table"))};
            }

        # check for db_xrefs, otherwise store it as an unknown tag
        foreach my $tag ($feat_object->get_all_tags() ) {
        if ($tag eq "db_xref") {
            # Assuming the possiblity of multiple tags of the same database
            foreach ($feat_object->get_tag_values($tag)) {
            push(@{$gbr{'Features'}->{$feature_group}->{$feature_id}->{db_xrefs}}, $_);
            }
        }
	# bug 5338 add gene_product_name and comments as shared or individual
	elsif ($tag eq 'translation') {
	  $gbr{'Features'}->{$feature_group}->{$feature_id}->{$tag}=join('',$feat_object->get_tag_values($tag));
	  #$gbr{'Features'}->{$feature_group}->{$feature_id}->{$tag} = [$feat_object->get_tag_values($tag)]
	}
#	elsif ($tag eq 'EC_number' || $tag eq 'product' || $tag eq 'note') {  # use this to treat it as an arrayref
	elsif ($tag eq 'EC_number' || $tag eq 'product' || $tag eq 'note' ||
	       $tag eq "gene" || $tag eq "locus_tag" || $tag eq "protein_id") {  # use this to treat it as an arrayref
	  $gbr{'Features'}->{$feature_group}->{$feature_id}->{$tag} = [$feat_object->get_tag_values($tag)]
	}	
 #       elsif ($tag eq "gene" || $tag eq "locus_tag" || $tag eq "protein_id" || 
	elsif ( $tag eq "transl_table") {
            # do nothing
        }
        else {
            ++$TagCount{"unknown_feature_tag_".$tag};
        }
        }
        } #if big || statement
        else {
            # record unknown primary tag
            ++$TagCount{"unknown_primary_tag_".$primary_tag};
        }
        } #else != source
    } #for all the SeqFeatures

    #Most likely no /transl_table so use default, 
    #otherwise check make sure only one non-standard transl_table appeared
    my @transl_tables = (keys %transl_tables);
    if (@transl_tables == 0) {
        print "Normal transl_table: 1\n";
        $gbr{'transl_table'} = 1;
    }
    elsif (@transl_tables == 1) {
        print "Using a different transl_table: ".$transl_tables[0]."\n";
        $gbr{'transl_table'} = $transl_tables[0];
    }
    else {
        die "Conflicting /transl_table values for CDSs in $gb_file";
    }

    #if there wasn't a valid strain provided, see if one can be grabbed from organism name
    unless ($gbr{'strain'}) {
        if ($gbr{'organism'} =~ m/((sp\.)|(str\.)|(subsp\.))\s*(.*)/) {
        $gbr{'strain'} = $+; #last match
        }
    }   
    }

    # die if we didn't actually process any records
    # see bug #3315 http://jorvis-lx:8080/bugzilla/show_bug.cgi?id=3315
    die "No genbank records in $gb_file" if ($n_seqs == 0);

    undef($seq_obj);
    return \%gbr;
}

sub parse_organism_to_prefix_mapping {
    my $map_file = shift;
    my $map = {};
    my $fh = FileHandle->new();
    my $lnum = 0;
    $fh->open($map_file, 'r') || die "unable to read from $map_file";

    while (my $line = <$fh>) {
	chomp($line);
	++$lnum;
	if ($line =~ /^\s*$|^[\#]|^\/\//) {
	    next;
	} elsif ($line =~ /^\s*(\S.*)\s+(\S+)$/) {
	    my($organism, $prefix) = ($1, $2);
	    $map->{$organism} = $prefix;
	} else {
	    die "unable to parse line $lnum of $map_file: $line";
	}
    }
    $fh->close();

    return $map;
}

sub to_bsml {
    my %gbr = %{$_[0]};
    my $doc = $_[1];
    my $genome = $doc->createAndAddGenome();
    my $genome_id = $genome->{attr}->{id}; # have to track the active genome and propagate through to Link Sequences

    # modified Cross-reference@database to exclude identifier
    # see bug #3278 http://jorvis-lx:8080/bugzilla/show_bug.cgi?id=3278
    my %xref;
    $xref{'taxon_id'} = $doc->createAndAddCrossReference(
                        'parent'          => $genome, 
                        'database'        => 'taxon',           # //Genome/Cross-reference/@database
                        'identifier'      => $gbr{'taxon_id'},  # //Genome/Cross-reference/@identifier
                        'identifier-type' => 'taxon_id'         # //Genome/Cross-reference/@identifier-type
                        );

    # add Cross-references
    # use GO standard database names: http://www.geneontology.org/doc/GO.xrf_abbs
    $xref{'accession'} = $doc->createAndAddCrossReference(
                        'parent'          => $genome, 
                        'database'        => 'GenBank',          # //Genome/Cross-reference/@database
                        'identifier'      => $gbr{'accession'},  # //Genome/Cross-reference/@identifier
                        'identifier-type' => 'genbank flat file' # //Genome/Cross-reference/@identifier-type
                        );

    $xref{'gi'} = $doc->createAndAddCrossReference(
                        'parent'          => $genome, 
                        'database'        => 'NCBI_gi',  # //Genome/Cross-reference/@database
                        'identifier'      => $gbr{'gi'}, # //Genome/Cross-reference/@identifier
                        'identifier-type' => 'gi'        # //Genome/Cross-reference/@identifier-type
                        );

    my $organism = $doc->createAndAddOrganism( 
                           'genome'  => $genome,
                           'genus'   => $gbr{'genus'},   # //Genome/Organism/@genus
                           'species' => $gbr{'species'}, # //Genome/Organism/@species
                           );

    # add genetic code / translation table / transl_table
    # values of 2, 3, 4, 5, 9, 13, 14, 16, 21, 22, 23 are also mitochondrial codes
    # see http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c#SG1
    # see bugzilla 2414
    my %mt_codes = (2 => 1, 3 => 1, 4 => 1, 5 => 1, 9 => 1, 13 => 1, 14 => 1, 16 => 1, 21 => 1, 22 => 1, 23 => 1);
    my $transl_table;
    if ($mt_codes{$gbr{'transl_table'}}) {
	$transl_table = $doc->createAndAddBsmlAttribute($organism, 'mt_genetic_code', $gbr{'transl_table'});
    }
    else {
	$transl_table = $doc->createAndAddBsmlAttribute($organism, 'genetic_code', $gbr{'transl_table'});
    }

    if($gbr{'molecule_name'}){
	$doc->createAndAddBsmlAttribute($organism, 'molecule_name', $gbr{'molecule_name'});
    }

    # only add a strain element if there is valid strain info
    # see bug #3518 http://jorvis-lx:8080/bugzilla/show_bug.cgi?id=3518
    die "Strain of just one space!!!" if ($gbr{'strain'} eq ' ');
    $gbr{'strain'} =~ s/\s+$//;
    if ($gbr{'strain'}) {
    my $strain_elem = $doc->createAndAddStrain(
					       'name'     => $gbr{'strain'},
					       'organism' => $organism
					       );
    }

    # class is always assembly
    # chromosome|plasmid will be denoted in the "secondary type" ie an Attribute-list
    # if seq->molecule (taken from LOCUS line) matches dna, molecule=dna
    # if it matches rna, molecule=rna
    # otherwise it dies

    my $seqid = $gbr{'accession'};
    my $fastafile = "$odir/$gbr{'gi'}.assembly.fsa";
    my $fastaname = "gi|".$gbr{'gi'}."|ref|".$gbr{'accession'}."|".$gbr{'organism'};
    # in order for this to get parsed correctly by bsml2chado's BSML::Indexer::Fasta component, 
    # Seq-data-import/@identifier must equal the fasta header up to the first space
    # currently dealing with this by:
    $fastaname =~ s/\s+/_/g;
    my $seq_type = 'assembly';

    # generate a new ergatis-compliant sequence id for the DNA sequence
    # (for the benefit of other components that fail to work with non-ergatis-compliant ids)
    if ($options{'generate_new_seq_ids'}) {
	my $project = $current_prefix || $options{'project'};
	$seqid = $idcreator->next_id('type' => $seq_type, 'project' => $project);
	$fastafile = "${odir}/${seqid}.fsa";
	# can't use ($seqid . ' ' . $fastaname) because other components 
	# follow different rules for identifying sequences (e.g., bsml2fasta.pl
	# expects the id to be the entire FASTA defline, not just the part up to
	# the first space)
	$fastaname = $seqid;
    }

    my $sequence = $doc->createAndAddSequence(
  		          $seqid,            #id
                          undef,             #title
                          $gbr{'seq_len'},   #length
                          $gbr{'molecule'},  #molecule
                          $seq_type,         #class
                          );

    # store the main genomic sequence
    # store everything in odir

    my $seq_obj = $doc->createAndAddSeqDataImport(
                          $sequence,                     # Sequence element object reference
                          'fasta',                       # //Seq-data-import/@format
                          $fastafile,                    # //Seq-data-import/@source
                          undef,                         # //Seq-data-import/@id
                          $fastaname                     # //Seq-data-import/@identifier
                          );

    # bug #4005 add //Sequence/Attribute/@name="defline", @content=//Seq-data-import/@identifier
    $doc->createAndAddBsmlAttribute($seq_obj, 'defline', $fastaname);

    open (my $FFILE, ">$fastafile") || die "Unable to open for writing $fastafile";
    print {$FFILE} fasta_out($fastaname, $gbr{'sequence'});
    close($FFILE);
    chmod(0666, $fastafile);

    # add the "secondary types" as  //Attribute-list/Attribute[@name="SO", content=<class>]
    # presumably class is from SO (not really), content=chromosome|plasmid
    $sequence->addBsmlAttributeList([{name => 'SO', content=> $gbr{'component'}}]);

    # add expanded ard.obo polymer_type (or SO term, if DNA|RNA)
    # see bug #3251 http://jorvis-lx:8080/bugzilla/show_bug.cgi?id=3251
    if ($gbr{polymer_type} eq "DNA" || 
	$gbr{polymer_type} eq "RNA" ||
	$gbr{polymer_type} eq "mRNA" ||
	$gbr{polymer_type} eq "ss_RNA_viral_sequence" ||
	$gbr{polymer_type} eq "postive_sense_ssRNA_viral_sequence" ||
	$gbr{polymer_type} eq "negative_sense_ssRNA_viral_sequence" ||
	$gbr{polymer_type} eq "ds_RNA_viral_sequence" || 
	$gbr{polymer_type} eq "single_stranded_DNA_chromosome"
       ) {
       $sequence->addBsmlAttributeList([{name => 'SO', content=> $gbr{polymer_type}}]);
    }
    else {
       $sequence->addBsmlAttributeList([{name => 'ARD', content=> $gbr{polymer_type}}]);
    }

    # Add strand attribute (this may not be used)
    # direct creation supported via createAndAddExtendedSequence, but this is simpler
    # valid values are: std-not-set|ss|ds|mixed|std-other
    if (defined($gbr{strand})) {
    $sequence->setattr( 'strand', $gbr{strand});
    }

    my $link_elem = $doc->createAndAddLink(
                       $sequence,
                       'genome',        # rel
                       '#'.$genome_id    # href
                       );

    #
    # features, feature-groups
    #
    #create the feature table
    my $feature_table_elem = $doc->createAndAddFeatureTable($sequence);

    foreach my $feature_group (keys %{$gbr{'Features'}}) {
    # testing
    next if $feature_group =~ /unknown/;
    #/testing
    my $feature_group_elem = $doc->createAndAddFeatureGroup(
                                $sequence, #the <Sequence> element containing the feature group
                                undef,                    # id
                                "$feature_group"             # groupset
                                );
    die ("Could not create <Feature-group> element for uniquename '$sequence'") if (!defined($feature_group_elem));

    #count the number of each parsed feature in each feature group
    # really this ought to be rewritten with a master feature type list
    my %feature_type = (gene => [], CDS => [], promoter => [], exon => [], intron => [], mRNA => [], tRNA => [], rRNA => []); #hash of array

    # map feature tags to their named versions
    # if these are universal then keep them all
    my $shared_gene_product_name = get_shared_feature_tag($gbr{'Features'}->{$feature_group}, 'product');
    my $shared_comment = get_shared_feature_tag($gbr{'Features'}->{$feature_group}, 'note');
    my $ec_numbers;
    
    foreach my $feature (keys %{$gbr{'Features'}->{$feature_group}}) {
        if(ref $gbr{'Features'}->{$feature_group}->{$feature}){
        $gbr{'Features'}->{$feature_group}->{$feature}->{id} = $feature;
        }
        
	# mapping from genbank feature tag to attributes for potentially shared attributes
	if ($shared_gene_product_name) {
	  $gbr{'Features'}->{$feature_group}->{$feature}->{attributes}->{'gene_product_name'} = $shared_gene_product_name;
	}	
        elsif (exists $gbr{'Features'}->{$feature_group}->{$feature}->{'product'}){
	  $gbr{'Features'}->{$feature_group}->{$feature}->{attributes}->{'gene_product_name'} = $gbr{'Features'}->{$feature_group}->{$feature}->{'product'};
	}

	if ($shared_comment) {
	  $gbr{'Features'}->{$feature_group}->{$feature}->{attributes}->{'comment'} = $shared_comment;
	}	
        elsif (exists $gbr{'Features'}->{$feature_group}->{$feature}->{'note'}){
	  $gbr{'Features'}->{$feature_group}->{$feature}->{attributes}->{'comment'} = $gbr{'Features'}->{$feature_group}->{$feature}->{'note'};
	}
	
	# TODO: Handle ec_number similar to above
        if (exists $gbr{'Features'}->{$feature_group}->{$feature}->{'EC_number'}){
	  $ec_numbers = $gbr{'Features'}->{$feature_group}->{$feature}->{'EC_number'};
        }

	# mapping from genbank feature tag to attributes
	if (exists $gbr{'Features'}->{$feature_group}->{$feature}->{'gene'}){
	  $gbr{'Features'}->{$feature_group}->{$feature}->{attributes}->{'gene'} = $gbr{'Features'}->{$feature_group}->{$feature}->{'gene'};
	}
	if (exists $gbr{'Features'}->{$feature_group}->{$feature}->{'protein_id'}){
	  $gbr{'Features'}->{$feature_group}->{$feature}->{attributes}->{'protein'} = $gbr{'Features'}->{$feature_group}->{$feature}->{'prodtein_id'};
	}

	# store locus_tag as Cross-reference not Attribute
	# database=NCBILocus
	if (exists $gbr{'Features'}->{$feature_group}->{$feature}->{'locus_tag'}){
#	  $gbr{'Features'}->{$feature_group}->{$feature}->{attributes}->{'locus'} = $gbr{'Features'}->{$feature_group}->{$feature}->{'locus_tag'};
	  foreach my $locus_dbxref (@{$gbr{'Features'}->{$feature_group}->{$feature}->{'locus_tag'}}) {
	    push(@{$gbr{'Features'}->{$feature_group}->{$feature}->{db_xrefs}}, "NCBILocus:".$locus_dbxref);
	  }
	}

        if ($feature =~ /gene/) {
        push(@{$feature_type{gene}}, $feature);
        }
        elsif ($feature =~ /CDS/) {
        push(@{$feature_type{CDS}}, $feature);
        }
        elsif ($feature =~ /promoter/) {
        push(@{$feature_type{promoter}}, $feature);
        }
        elsif ($feature =~ /exon/) {
        push(@{$feature_type{exon}}, $feature);
        }
        elsif ($feature =~ /intron/) {
        push(@{$feature_type{intron}}, $feature);
        }
        elsif ($feature =~ /mRNA/) {
        push(@{$feature_type{mRNA}}, $feature);
        }
        elsif ($feature =~ /tRNA/) {
        push(@{$feature_type{tRNA}}, $feature);
        }
        elsif ($feature =~ /rRNA/) {
        push(@{$feature_type{rRNA}}, $feature);
        }
        else {
        die "Unexpected feature ($feature) in $feature_group";
        }
    }
    
    #
    # add gene
    #
    if (@{$feature_type{gene}} == 1) {
        &addFeature($gbr{'Features'}->{$feature_group}->{$feature_type{gene}->[0]}, $doc, $genome_id, $feature_table_elem, $feature_group_elem);
        # support for feature_groups of a just a gene, gene+tRNA, gene+rRNA
        # see bug #3298 http://jorvis-lx:8080/bugzilla/show_bug.cgi?id=3298, bug #5328
        if (scalar(keys %{$gbr{'Features'}->{$feature_group}}) == 1) {
        next; # goto next feature_group if this is the only thing (ie don't die)
        }
        elsif ( (scalar(keys %{$gbr{'Features'}->{$feature_group}}) == 2) && (@{$feature_type{tRNA}} == 1 )) {
        &addFeature($gbr{'Features'}->{$feature_group}->{$feature_type{tRNA}->[0]}, $doc, $genome_id, $feature_table_elem, $feature_group_elem);
	# see bug 5328: need an exon in the feature group
	derive_and_add_exons_from_Feature($feature_type{tRNA}, $feature_group, $doc, $genome_id, $feature_table_elem, $feature_group_elem);
        next; # next feature_group (no die)
        }
        elsif ( (scalar(keys %{$gbr{'Features'}->{$feature_group}}) == 2) && (@{$feature_type{rRNA}} == 1 )) {
        &addFeature($gbr{'Features'}->{$feature_group}->{$feature_type{rRNA}->[0]}, $doc, $genome_id, $feature_table_elem, $feature_group_elem);
	derive_and_add_exons_from_Feature($feature_type{rRNA}, $feature_group, $doc, $genome_id, $feature_table_elem, $feature_group_elem);
        next; # next feature-group, no die
        }
    }
    elsif (@{$feature_type{gene}} > 1) {
        die "Multiple gene tags (".@{$feature_type{gene}}.") in feature group $feature_group";
    }
    # support for multiple CDSs bug #3299 http://jorvis-lx:8080/bugzilla/show_bug.cgi?id=3299
    elsif (@{$feature_type{CDS}} >= 1) { # use CDSs
        my $gene_featref = &copy_featref($gbr{'Features'}{$feature_group}{$feature_type{CDS}->[0]}, 'gene');
        $gene_featref->{id} =~ s/CDS/gene_from_CDS/;
        # obtain max span
        foreach my $cds (@{$feature_type{CDS}}) {
        if ($gbr{'Features'}{$feature_group}{$cds}->{start} < $gene_featref->{start}) {
            $gene_featref->{start} = $gbr{'Features'}{$feature_group}{$cds}->{start};
            $gene_featref->{start_type} = $gbr{'Features'}{$feature_group}{$cds}->{start_type};
        }
        if ($gbr{'Features'}{$feature_group}{$cds}->{end} > $gene_featref->{end}) {
            $gene_featref->{end} = $gbr{'Features'}{$feature_group}{$cds}->{end};
            $gene_featref->{end_type} = $gbr{'Features'}{$feature_group}{$cds}->{end_type};
        }
        }
        my $gene_elem = &addFeature($gene_featref, $doc, $genome_id, $feature_table_elem, $feature_group_elem);
    }
    else {
        die "Unable to create gene object in $feature_group";
    }

    #
    # add transcript
    #
    # adding support for multiple mRNAs.  See bug #3308 http://jorvis-lx:8080/bugzilla/show_bug.cgi?id=3308
    if (@{$feature_type{mRNA}} >= 1) {
        my $trans_featref = &copy_featref($gbr{'Features'}{$feature_group}{$feature_type{mRNA}->[0]}, 'transcript');
        $trans_featref->{id} =~ s/mRNA/transcript_from_mRNA/; #no.
        # obtain max span
        foreach my $mrna (@{$feature_type{mRNA}}) {
        if ($gbr{'Features'}{$feature_group}{$mrna}->{start} < $trans_featref->{start}) {
            $trans_featref->{start} = $gbr{'Features'}{$feature_group}{$mrna}->{start};
            $trans_featref->{start_type} = $gbr{'Features'}{$feature_group}{$mrna}->{start_type};
        }
        if ($gbr{'Features'}{$feature_group}{$mrna}->{end} > $trans_featref->{end}) {
            $trans_featref->{end} = $gbr{'Features'}{$feature_group}{$mrna}->{end};
            $trans_featref->{end_type} = $gbr{'Features'}{$feature_group}{$mrna}->{end_type};
        }
        }
        my $trans_elem = &addFeature($trans_featref, $doc, $genome_id, $feature_table_elem, $feature_group_elem);
        if($ec_numbers){
        foreach my $ec_number (@$ec_numbers){
            $trans_elem->addBsmlAttributeList([{name => 'EC', content=> $ec_number}]);
        }
        }
    }
    elsif (@{$feature_type{gene}} == 1) {
        
        my $trans_featref = &copy_featref($gbr{'Features'}{$feature_group}{$feature_type{gene}->[0]}, 'transcript');
        $trans_featref->{id} =~ s/gene/transcript_from_gene/;
        my $trans_elem = &addFeature($trans_featref, $doc, $genome_id, $feature_table_elem, $feature_group_elem);
        if ($ec_numbers){
            foreach my $ec_number (@$ec_numbers){
                $trans_elem->addBsmlAttributeList([{name => 'EC', content=> $ec_number}]);
            }
        }
    }
    elsif (@{$feature_type{CDS}} >= 1) { # use CDSs
        my $trans_featref = &copy_featref($gbr{'Features'}{$feature_group}{$feature_type{CDS}->[0]}, 'transcript');
        $trans_featref->{id} =~ s/CDS/transcript_from_CDS/;
        # obtain max span
        foreach my $cds (@{$feature_type{CDS}}) {
        if ($gbr{'Features'}{$feature_group}{$cds}->{start} < $trans_featref->{start}) {
            $trans_featref->{start} = $gbr{'Features'}{$feature_group}{$cds}->{start};
            $trans_featref->{start_type} = $gbr{'Features'}{$feature_group}{$cds}->{start_type};
        }
        if ($gbr{'Features'}{$feature_group}{$cds}->{end} > $trans_featref->{end}) {
            $trans_featref->{end} = $gbr{'Features'}{$feature_group}{$cds}->{end};
            $trans_featref->{end_type} = $gbr{'Features'}{$feature_group}{$cds}->{end_type};
        }
        }
        my $trans_elem = &addFeature($trans_featref, $doc, $genome_id, $feature_table_elem, $feature_group_elem);
        if($ec_numbers){
        foreach my $ec_number (@$ec_numbers){
            $trans_elem->addBsmlAttributeList([{name => 'EC', content=> $ec_number}]);
        }
        }
    }
    else {
        die "Unable to create transcript object in $feature_group";
    }

    #
    #  add CDS
    #
    # known bugs: 
    #  -multiple CDSs are straight up just added
    #  -CDSs might be joined segments and this is kind of ignored
    # see bug #3299 http://jorvis-lx:8080/bugzilla/show_bug.cgi?id=3299 for discussion
    if (@{$feature_type{CDS}} >= 1) { # use multiple CDSs
        foreach my $cds (@{$feature_type{CDS}}) {
        # create translation if there isn't one
        unless ($gbr{'Features'}{$feature_group}{$cds}->{translation}) {
            my $codon_table  = Bio::Tools::CodonTable -> new ( -id => $gbr{transl_table} );
            $gbr{'Features'}{$feature_group}{$cds}->{translation} = $codon_table->translate($gbr{'Features'}{$feature_group}{$cds}->{spliced_seq});
        }
        &addFeature($gbr{'Features'}->{$feature_group}->{$cds}, $doc, $genome_id, $feature_table_elem, $feature_group_elem);
        }
    }
    # create CDS from mRNA if present
    # see bug #3300 http://jorvis-lx:8080/bugzilla/show_bug.cgi?id=3300
    # revised in bug #3308 to create one CDS per mRNA
    elsif (@{$feature_type{mRNA}} >= 1) {
        foreach my $mrna (@{$feature_type{mRNA}}) {
        my $cds_featref = &copy_featref($gbr{'Features'}{$feature_group}{$mrna}, 'CDS');
        $cds_featref->{id} =~ s/mRNA/CDS_from_mRNA/;
        # create translation from known transl_table value ($gbr{'transl_table'})
        # doc here: http://search.cpan.org/~birney/bioperl-1.2.3/Bio/Tools/CodonTable.pm
        my $codon_table  = Bio::Tools::CodonTable -> new ( -id => $gbr{transl_table} );
        $cds_featref->{translation} = $codon_table->translate($cds_featref->{spliced_seq});

        my $cds_elem = &addFeature($cds_featref, $doc, $genome_id, $feature_table_elem, $feature_group_elem);
#       $doc->createAndAddBsmlAttribute($cds_elem, "comment", "Derived from mRNA tag");
        }
    }
    else {
      print Dumper($gbr{'Features'}{$feature_group})."\n\n";
      
        die "Unable to create CDS object in CDS: @{$feature_type{CDS}}, mRNA: @{$feature_type{mRNA}}, feature group: $feature_group";
    }


    #
    # add exon
    #
    if (@{$feature_type{exon}} > 0) { # use exons if available
        foreach my $exon (@{$feature_type{exon}}) {
        #check if joined locations?
        &addFeature($gbr{'Features'}->{$feature_group}->{$exon}, $doc, $genome_id, $feature_table_elem, $feature_group_elem);
        }
    }
    # derive from mRNA if present
    elsif (@{$feature_type{mRNA}} >= 1) {
      derive_and_add_exons_from_Feature($feature_type{mRNA}, $feature_group, $doc, $genome_id, $feature_table_elem, $feature_group_elem);      	
    }
    # See bug #3299 http://jorvis-lx:8080/bugzilla/show_bug.cgi?id=3299 for discussion of multiple CDSs
    # Regardless of the number of CDSs, each segment is used as an exon
    # see bug $3305 http://jorvis-lx:8080/bugzilla/show_bug.cgi?id=3305
    elsif (@{$feature_type{CDS}} >= 1) { # otherwise one exon for each CDS fragment
      derive_and_add_exons_from_Feature($feature_type{CDS}, $feature_group, $doc, $genome_id, $feature_table_elem, $feature_group_elem);      	
    }
    # what if the CDS was derived from an MRNA, see bug #3300 http://jorvis-lx:8080/bugzilla/show_bug.cgi?id=3300
    elsif (@{$feature_type{mRNA}} == 1) {
      derive_and_add_exons_from_Feature($feature_type{mRNA}, $feature_group, $doc, $genome_id, $feature_table_elem, $feature_group_elem);
    }
    else {
        die "Unable to create exon object in $feature_group";
    }

    #
    # add promoters and introns
    #
    if (@{$feature_type{promoter}} > 0) { # use promoters if available
        foreach my $promoter (@{$feature_type{promoter}}) {
        #check if joined locations?
        &addFeature($gbr{'Features'}->{$feature_group}->{$promoter}, $doc, $genome_id, $feature_table_elem, $feature_group_elem);
        }
    }
    if (@{$feature_type{intron}} > 0) { # use introns if available
        foreach my $intron (@{$feature_type{intron}}) {
        #check if joined locations?
        &addFeature($gbr{'Features'}->{$feature_group}->{$intron}, $doc, $genome_id, $feature_table_elem, $feature_group_elem);
        }
    }
    } #foreach feature_group
} #/to_bsml


# create and add an exon to the $doc for each
# type of feature and each location of that feature
sub derive_and_add_exons_from_Feature {
  my ($feature_type, $feature_group, $doc, $genome_id, $feature_table_elem, $feature_group_elem) = @_;
  
  my $numexons = 0;
  foreach my $a_feature (@{$feature_type}) {
    foreach my $loc (@{$gbr{'Features'}{$feature_group}{$a_feature}{locations}}) {
      my $exon_featref = &copy_featref($gbr{'Features'}{$feature_group}{$a_feature}, 'exon');
      $exon_featref->{id} = "exon_".$numexons."_from_".$exon_featref->{id};
      $exon_featref->{start} = $loc->{start};
      $exon_featref->{start_type} = $loc->{start_pos_type};
      $exon_featref->{end} = $loc->{end};
      $exon_featref->{end_type} = $loc->{end_pos_type};
      my $exon_elem = &addFeature($exon_featref, $doc, $genome_id, $feature_table_elem, $feature_group_elem);
      ++$numexons;
    }
  } 
}

# check a feature group for a unique value of some feature tag
# if it differs among any two features, return 0
# if it is the same or absent in all features, return it
sub get_shared_feature_tag {
  my $feature_group = shift;
  my $tag = shift;
  
  my $tagref = undef;
  
  foreach my $feature (keys %$feature_group) {
    my $fref = $feature_group->{$feature};
    next unless defined($fref->{$tag});  # skip if this isn't defined, no penalty
    unless (defined($tagref)) {
      $tagref = $fref->{$tag};
      next;
    }
    # otherwise, compare
    for (my $i = 0; $i < @{$fref->{$tag}}; ++$i) {
      unless ($tagref->[$i] eq $fref->{$tag}->[$i]) {
	# warn "Mismatch on $tag between (".$tagref->[$i].") and (".$fref->{$tag}->[$i]."): no shared tag $tag for feature_group ($feature_group)";
	return 0;
      }
    }
    # extra in current
    if (@$tagref > @{$fref->{$tag}}) {
      warn "Tag count mismatch, no shared tag ($tag) for feature_group ($feature_group)";      
      return 0;
    }
  }
  return $tagref;
}



#this doesn't conflict w/ BSML namespace
sub addFeature {
    my $featref = shift;
    my $doc = shift;
    my $genome_id = shift;
    my $feature_table_elem = shift;
    my $feature_group_elem = shift;

    my $old_id = $featref->{'id'};
    
    my $class = $featref->{'class'};
    my $is_complement = $featref->{'is_complement'};
    my $start = $featref->{'start'} - 1; #convert to space based coords
    my $end = $featref->{'end'}; #this should NOT be +1, right?
    my $feature_tag = $featref->{'feature_tag'};
    my $feature_value = $featref->{'feature_value'};

    my $location_type = $featref->{'location_type'};
    my $start_type = $featref->{'start_type'};
    my $end_type = $featref->{'end_type'};
    
    my $project = $current_prefix || $options{'project'};
    my $id = $idcreator->next_id('type' => $class, 'project' => $project);

    $feature_id_lookup->{$old_id} = $id;
    
    #Add //Feature
    my $feature_elem = $doc->createAndAddFeature(
                         $feature_table_elem,  # <Feature-table> element object reference
                         "$id",                # id
                         undef,                # title
                         $class,               # class
                         undef,                # comment
                         undef,                # displayAuto
                         );
    if (!defined($feature_elem)) {
    die("Could not create feature element ($id)");
    }

    #Add the location element
    if ($location_type eq "EXACT") {
    if ($start_type eq "EXACT" && $end_type eq "EXACT") {
        my $loc_loc = $doc->createAndAddIntervalLoc($feature_elem, $start, $end, $is_complement);
    }
    elsif ($start_type eq "BEFORE" && $end_type eq "EXACT") {
        my $loc_loc = $doc->createAndAddIntervalLoc($feature_elem, $start, $end, $is_complement, 1, 0);
    }
    elsif ($start_type eq "EXACT" && $end_type eq "AFTER") {
        my $loc_loc = $doc->createAndAddIntervalLoc($feature_elem, $start, $end, $is_complement, 0, 1);
    }
    elsif ($start_type eq "BEFORE" && $end_type eq "AFTER") {
        my $loc_loc = $doc->createAndAddIntervalLoc($feature_elem, $start, $end, $is_complement, 1, 1);
    }
    else {
        die "Unsupported start_type ($start_type) or end_type ($end_type)";
    }
    }
    elsif ($location_type eq "IN-BETWEEN") {
    #Check that it's a reasonable conversion of base-coordinates into interval coordinates describing an interbase position
    ($start == $end - 2) || die "Incongruous site-loc: start=$start, end=$end";
    print "Creating site: $start $end\n";
    #Add a Site-Loc (since in interbase coords, site is one before the end)
    my $site_loc = $doc->createAndAddSiteLoc($feature_elem, $end-1,$is_complement, $class );
    }
    else {
    die "Unsupported location_type: $location_type";
    }

    #Add a reference to the new feature to the //Feature-group
    my $feature_group_member_elem = $doc->createAndAddFeatureGroupMember( $feature_group_elem,  # <Feature-group> element object reference
                                      $id,                  # featref
                                      $class,               # feattype
                                      undef,                # grouptype
                                      undef,                # cdata
                                    ); 
    if (!defined($feature_group_member_elem)){
    die("Could not create <Feature-group-member> element object reference for gene model subfeature '$id'");
    }

    # add Attribute for which tag was used to create feature-group
    # See http://jorvis-lx:8080/bugzilla/show_bug.cgi?id=2506#c3 for discussion
    if (my $attribute_feature_value = $feature_group_elem->returnBsmlAttr($feature_tag)) {
    unless ($feature_value eq $attribute_feature_value->[0]) {  #check it's the same
        die("conflict in feature_tag '$feature_tag', existing value '$attribute_feature_value->[0]' conflicts with '$feature_value'");
    }
    }
    else {
    my $attribute_feature = $doc->createAndAddBsmlAttribute($feature_group_elem, $feature_tag, $feature_value);
    if (!defined($attribute_feature)){
        die("Could not create <Attribute> element ");
    }
    }

    # TODO: fix up above so it's like this
    # add gene_product_name, comment (others) (which may have been derived from feature_group: bug 5338)
    if (defined($featref->{attributes})) {
      foreach my $attribute (keys %{$featref->{attributes}}) {	
	foreach (@{$featref->{attributes}->{$attribute}}) {
	  $doc->createAndAddBsmlAttribute($feature_elem, $attribute, $_)
	}  
      }      
    }

    # add any db_xrefs associated with the feature as Cross-references
    # list of database names taken from http://www.geneontology.org/doc/GO.xrf_abbs
    # list of databases in genbank records: http://www.ncbi.nlm.nih.gov/collab/db_xref.html
    if (defined($featref->{db_xrefs})) {
    foreach (@{$featref->{db_xrefs}}) {
        # print "\tdb_xref\t$_\n";
        # new usage convention is everything up to the first : is db, everything after (including :) is id
        # see bug #3278 comment #3 http://jorvis-lx:8080/bugzilla/show_bug.cgi?id=3278#c3
        ($_ =~ /([^:]*):(.*)/) || die "Unable to parse database and identifier from db_xref ($_)";
        my $database = $1;
        my $identifier = $2;

        # die if it's some new kind of database I never saw before
        my %known_dbxrefs = ( GI => 1, GeneID => 1, CDD => 1, ATCC => 1, Interpro => 1, UniProtKB => 1, GOA => 1,
                  HSSP => 1, PSEUDO => 1, DDBJ => 1, COG => 1, ECOCYC => 1, ASAP => 1, ISFinder => 1,
                  EMBL => 1, GenBank => 1, InterPro => 1, 'UniProtKB/TrEMBL' => 1, 'UniProtKB/Swiss-Prot' => 1,
                  dictyBase => 1, FlyBase => 1, VectorBase => 1, SGD => 1, SGDID => 1, NCBILocus => 1);
        unless (defined($known_dbxrefs{$database})) {
	  warn "Unknown database in dbxref ($database)";
	  if ($options{skip_unknown_dbxref}) {
	    warn "Skipping feature with dbxref $database:$identifier";
	    next;
	  }
	}
        
        # mod database to GO xref standard as neccessary http://www.geneontology.org/doc/GO.xrf_abbs
	# VectorBase from http://neuron.cse.nd.edu/vectorbase/index.php/GenBank_submission#db_xref
        if ($database eq 'GI') {
        $database = 'NCBI_gi';
        }
        elsif ($database eq 'Interpro') {
        $database = 'InterPro';
        }
        elsif ($database eq 'UniProtKB') {
        $database = 'UniProt';
        }
        elsif ($database eq 'ECOCYC') {
        $database = 'EcoCyc';
        }
        elsif ($database eq 'dictyBase') {
        $database = 'DDB';
        }
        elsif ($database eq 'FlyBase') {
        $database = 'FB';
        }
        elsif ($database eq 'SGDID') {
        $database = 'SGD';
        }
        elsif ($database eq 'COG') {
        if ($identifier =~ /^\s+COG/) {
            $database = 'COG_Cluster';
        }
        elsif ($identifier =~ /^\d$/) { # single digit
            $database = 'COG_Pathway';
        }
        elsif ($identifier =~ /^\w$/) { # single letter (since not digit)
            $database = 'COG_Function';
        }
        else {
            die "Unable to parse COG database from identifier ($identifier)";
        }
        }

	# add in identifier-type
	$doc->createAndAddCrossReference(
					 'parent'          => $feature_elem, 
					 'database'        => $database,          # //Genome/Cross-reference/@database
					 'identifier'      => $identifier,        # //Genome/Cross-reference/@identifier
					 # identifier-type stored as dbxref.version
					 #'identifier-type' => 'genbank flat file' # //Genome/Cross-reference/@identifier-type
					);
    }
    }

    # if the feature type is CDS
    # - add dna sequence (*or* rna)
    # - add polypeptide feature
    # - add aa sequence for polypeptide
    if ($class eq 'CDS') {
    &addSeqToFeature($featref, $doc, $genome_id, $feature_table_elem, $feature_group_elem, $feature_elem, $gbr{'molecule'});
    if ($featref->{'translation'}) {
        # create dummy featref for polypeptide
        # use copy_featref?
        my $poly_featref = {
                 'id' => $featref->{'id'},
                 'class' => 'polypeptide',
                 'is_complement' => $featref->{'is_complement'},
                 'start' => $featref->{'start'},
                 'end' => $featref->{'end'}, #this should NOT be +1, right?
                 'feature_tag' => $featref->{'feature_tag'},
                 'feature_value' => $featref->{'feature_value'},                 
                 'location_type' => $featref->{'location_type'},
                 'start_type' => $featref->{'start_type'},
                 'end_type' => $featref->{'end_type'},
                 'feature_len' => length($featref->{'translation'}),
                 'spliced_seq' => $featref->{'translation'}
                 };
        $poly_featref->{id} =~ s/CDS/polypeptide/; # change id
        my $poly_feat_elem = &addFeature($poly_featref, $doc, $genome_id, $feature_table_elem, $feature_group_elem);
        &addSeqToFeature($poly_featref, $doc, $genome_id, $feature_table_elem, $feature_group_elem, $poly_feat_elem, 'aa');
    } # if translation
    }# if CDS

    # return the feature that was created
    return $feature_elem;

} #end addFeature



sub addSeqToFeature {
    # add Sequence object
    # link sequence to genome
    # add sequence link to feature
    my $featref = shift;
    my $doc = shift;
    my $genome_id = shift;
    my $feature_table_elem = shift;
    my $feature_group_elem = shift;
    my $feature_elem = shift;
    my $molecule_type= shift;

    my $old_id = $featref->{'id'};

    my $id = '';
    if (defined($feature_id_lookup->{$old_id})) {
        $id = $feature_id_lookup->{$old_id};
    } else {
        die "Couldn't find feature id in lookup for '$old_id'";
    }
    
    my $class = $featref->{'class'};
    my $feature_len = $featref->{'feature_len'};
    my $spliced_seq = $featref->{'spliced_seq'};

    my $seq_name = $id."_seq";

    # create sequence object
    my $feature_seq = $doc->createAndAddSequence(
                         $seq_name,    # id
                         undef,        # title
                         $feature_len, # $gbr{'seq_len'},   #length
                         $molecule_type,        # molecule
                         $class  #class
                         );
    # write sequence object to file
    unless (defined($fsa_files{$class})) {
	$fsa_files{$class}{name} = "$odir/$gbr{gi}.$class.fsa";
	open($fsa_files{$class}{fh}, ">".$fsa_files{$class}{name});
    }

    my $fastaname = $id;
    #Seq-data-import/@identifier must equal the fasta header up to the first space  
    my $feature_seq_obj = $doc->createAndAddSeqDataImport(
                              $feature_seq,                  # Sequence element object reference
                              'fasta',                       # //Seq-data-import/@format
                              $fsa_files{$class}{name},      # //Seq-data-import/@source
                              undef,                         # //Seq-data-import/@id
                              $fastaname                     # //Seq-data-import/@identifier
                              );

    # bug #4005 add //Sequence/Attribute/@name="defline", @content=//Seq-data-import/@identifier
    $doc->createAndAddBsmlAttribute($feature_seq, 'defline', $fastaname);
    
    print {$fsa_files{$class}{fh}} fasta_out($fastaname, $spliced_seq);

    # link sequence to genome
    # I am reasonably certain that not including the genome@id value for the href is okay
    # the line above is a ***lie***
    my $feature_link_elem = $doc->createAndAddLink(
                           $feature_seq,
                           'genome',        # rel
#                          '#'.$genome->{'attr'}->{'id'}    # href
                           '#'.$genome_id
                           );


    # link feature to sequence
    my $feature_seq_link = $doc->createAndAddLink(
                          $feature_elem,
                          'sequence',      # rel
                          '#'.$seq_name    # href (same as Sequence@id)
                          );
    return $feature_seq;
}



# create a copy of a featref
# inputs: featref, newfeatref, [new class]
sub copy_featref {
    my $newfeat;
    my $featref = shift;
    my $class = shift;
    foreach (keys %{$featref}) {
    $newfeat->{$_} = $featref->{$_};
    }

    if (defined($class)) {
    $newfeat->{class} = $class;
    }

    return $newfeat;

}

#-------------------------------------------------------------------------
# fasta_out()
# taken without remorse from legacy2bsml.pl
#-------------------------------------------------------------------------
sub fasta_out {
    #This subroutine takes a sequence name and its sequence and
    #outputs a correctly formatted single fasta entry (including newlines).  

    #$logger->debug("Entered fasta_out") if $logger->is_debug;

    my ($seq_name, $seq) = @_;

    die "Empty sequence ($seq_name) passed to fasta_out!" unless defined($seq);

    my $fasta=">"."$seq_name"."\n";
    $seq =~ s/\s+//g;
    for(my $i=0; $i < length($seq); $i+=60){
    my $seq_fragment = substr($seq, $i, 60);
    $fasta .= "$seq_fragment"."\n";
    }
    return $fasta;

}

