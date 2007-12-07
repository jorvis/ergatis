#!/usr/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 SYNOPSIS

    gff3Valid.pl some_file.gff3

=head1 DESCRIPTION

This is a GFF3 validator for the BRC Interoperability Working Group (iowg.brcdevel.org).  This is
based on Lincoln Stein's GFF3 1.01 spec but is more strict.  Our usage is defined here:

    http://iowg.brcdevel.org/gff3.html
    
This first version of the validator does not yet check many of the rules given in the document
above.  The VALIDATION section below describes each of the current validation checks.

Ontology directive targets may currently be local files or URLs (http or ftp).  URNs are not
currently supported.

=head1 VALIDATION

The list below represents the current set of checks performed by this validator.

=over 4

=item *
note: general character class: [a-zA-Z0-9.:_-]

=item *
the gff3 file passed must exist

=item *
the first line in the file must be the version definition

=item *
implied fasta directive usage not encouraged (warning)

=item *
feature-ontology directive should be passed (warning)

=item *
feature-ontology directive URI, if passed, must be accessible

=item *
attribute-ontology directive should be passed (warning)

=item *
attribute-ontology directive URI, if passed, must be accessible

=item *
format of ##sequence-region directive must be: ##sequence-region seqid start end

=item *
directives not in the spec are not allowed

=item *
data lines must have 9 tab-delimited columns

=item *
column1, seqid, must match the general character class

=item *
column2, source, cannot be whitespace or empty (. for none)

=item *
column3, type, must have an entry in the feature ontology

=item *
column4, start, must be an integer

=item *
column5, end, must be an integer

=item *
start must be <= end

=item *
column6, score, must be a floating point number such as 5e-3 -5e+10 5E-3 5E+20 but NOT e-10

=item *
column7, strand, must be one of [-+.?]

=item *
column8, phase, must be one of [012.]

=item *
column9, attributes, must be in key=value format, with multiple values within a key separated by commas and multiple key=value pairs separated by a semi-colon

=item *
each attribute key must be defined in the attribute ontology

=item *
Dbxref attributes must be in the form db:id

=item *
Ontology_term attributes must be in the form AA:0000000

=item *
ec_number attributes must be in the form of N.N.N.N

=item *
molecule_type attributes must be either dsDNA | ssDNA | dsRNA | ssRNA

=item *
topology attributes must be either linear | circular

=item *
localization attributes must be either mitochondrion | plastid | episome | plasmid | nuclear | chromosomal

=item *
if the feature type (column 3) is mRNA or SO:0000234, the following attributes are required: ID, Name, Dbxref, description

=item *
if the feature type (column 3) is contig, SO:0000149, supercontig or SO:0000148, the following attributes are required: ID, Name, molecule_type, organism_name, translation_table, Dbxref.  Also, Dbxref must point to a taxon id.

=item *
each seqid referenced in column 1 must also have a contig or supercontig definition within the GFF3 file.

=item *
each seqid referenced in column 1 must also have a FASTA sequence within the GFF3 file.

=item *
there must be at least one row of type mRNA cds and exon

=back

=cut

use strict;
use GO::Parser;
use LWP::Simple;

#url of website with documentation info
#my $url = "file:///usr/local/devel/ANNOTATION/agussman/brcdevel.org/public_html/iowg/gff3.html"; #for test
my $url = "http://iowg.brcdevel.org/gff3.html"; #for real

my $file = shift || die "\nUSAGE: $0 file.gff\n\n";

## open the GFF file
open(my $gff3_fh, "<$file") || die "can't read $file : $!";

## this counts the number of errors (not warnings) found in the document
##  it is also the value returned by the script.
my $error_count = 0;

## we are using a more restrictive general character class than
##  that given in the GFF3 spec.
##  spec allows: [a-zA-Z0-9. :^*$@!+_?-]
##  we   allow : [a-zA-Z0-9.:_-]
my $gencc = '[a-zA-Z0-9.:_-]';

## the first line in the file must be the version definition
my $first_line = readline $gff3_fh;
if ($first_line !~ /^##gff-version\s+3/) {
    record_error("First line of file $file must be: \#\#gff-version 3 not: $first_line", "a_gff3d");
}

my ($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes);
my $checking_feature_ontology = 0;
my $checking_attribute_ontology = 0;
my ($feat_ont, $att_ont);
my $within_fasta = 0;

## holds a list of all the seqids defined in the document (reference seqs).  Later
##  we'll check to make sure all the seqids referenced were also defined.
my %seqids_defined;

## holds a list (and counts) of all the seqids referenced in the document.  These
##  should each have a FASTA entry of this same id (checked later)
my %seqids_referenced;

## holds the FASTA seqs found within the document (just ids so far)
my %seqs;

## count the number of mrna, cds, and exons
## later, ensure that we saw one of each
#my %type_count = ( mrna => 0, cds => 0, exon => 0);
my %type_count = ( gene => 0, cds => 0);

while (<$gff3_fh>) {
    chomp;

    ## skip comments
    next if (/^#[^#]/);

    ## skip blank lines
    next if (/^\s*$/);
    
    if ($within_fasta) {
        if (/^\>(\S+)/) {
            $seqs{$1} = 1;
        }
        
        next;
    }
    
    ## fasta directive
    if (/^##FASTA/) {
        $within_fasta = 1;
        next;
    }
    
    ## warn if the user uses an implied fasta directive
    if (/^\>/) {
        record_error("Implied FASTA directive found on line $. but no ##FASTA directive encountered","a_fasta");
        $within_fasta = 1;
        next;
    }
    
    ## references resolved?
    if (/^###/) {
        ## we can do some memory clearing here later
        
        next;
    }
    
    ##feature-ontology URI
    if (/^##feature-ontology\s+(.+?)\s*$/) {
        ($feat_ont, $checking_feature_ontology) = load_ontology($1);
        next;
    }

    ##attribute-ontology URI    
    if (/^##attribute-ontology\s+(.+?)\s*$/) {
        ($att_ont, $checking_attribute_ontology) = load_ontology($1);
        next;
    }

    ##sequence-region seqid start end
    if (/^##sequence-region/) {
        ## make sure it is in the format: seqid start end
        unless (/^##sequence-region\s+(.+?)\s+(\d+)\s+(\d+)/) {
            record_error("Sequence-region directive improperly formatted on line $.  should be ##sequence-region seqid start end", "a_gff_line1");
        }
        next;
    }

    ##source-ontology seqid start end
    if (/^##source-ontology/) {
        ## TEMP:  not doing anything with this yet
        next;
    }

    ## TEMP
    ## aren't handling other directives yet    
    if (/^##.+/) {
        record_error("Unknown directive ($_) on line $.", "a_gff_directives");
        next;
    };
    
    ###################
    ## syntax checks ##
    ###################
    
    my @cols = split(/\t/);
    
    ## we should have 9 columns.
    if ($#cols != 8) {
        record_error("Incorrect column count, line $.", "a_gff3d");
        next;
    }
    
    ## more transparent
    ($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = @cols;
    undef @cols;

    $type = lc($type); #because we don't care
    
    ## column 1: "seqid"
    #  must match the general character class
    if ($seqid !~ /^${gencc}+$/) {
        record_error("Seqid ($seqid) on line $. contains invalid characters.  Valid characters are $gencc", "a_gff3d");
    }

    ## remember this reference so we can make sure a FASTA sequence is defined later
    $seqids_referenced{$seqid}++;
    
    ## column 2: "source"
    #  no character class restrictions.
    #  cannot be all whitespace or null
    if ($source =~ /^\s*$/) {
        record_error("Source on line $. seems to be empty.  (use . for null fields)", "a_gff3d");
    }
    
    ## column 3: "type"
    ## are we doing ontology checking?
    if ($checking_feature_ontology) {
        my $term = fetch_term($type, $feat_ont);
        if (! $term) {
            record_error("Type ($type) not defined in feature ontology on line $.", "a_canon");
        }
    }
    ## increment count of type
    ++$type_count{$type};

    ## column 4 & 5: "start" and "end"
    #  start must be a positive integer
    if ($start !~ /^[0-9]+$/) {
        record_error("Start ($start) on line $. not defined as an integer", "a_gff_startend");
    }
    
    # end must be a positive integer
    if ($end !~ /^[0-9]+$/) {
        record_error("End ($end) on line $. not defined as an integer", "a_gff_startend");
    }
    
    # start must be less than or = end
    if ($start > $end) {
        record_error("Start ($start) not <= end ($end) on line $.","a_start_end");
    }
    
    ## column 6: "score"
    #  must be a floating point number like 5 5.02 0.0005 5e-3 5e+10 5E-3 5E+20
    #  check for scientific notation first
    if ($score =~ /e/i) {
        # valid notations are like 5e-3 -5e+10 5E-3 5E+20 but NOT e-10
        if ($score !~ /^[\-0-9.]+e[+\-][0-9]+$/i) {
            record_error("Score ($score) given in unrecognized scientific notation on line $.", "a_gff3d");
        }

    #  else check for regular floating point numbers
    } elsif ($score !~ /^[\-0-9.]+$/) {
        record_error("Score ($score) on line $. not recognized as a floating point number", "a_gff_score");
    }
    
    ## column 7: "strand"
    #  legal values for strand are - + . ?
    if ($strand !~ /^[\-+.?]$/) {
        record_error("Strand value ($strand) on line $. not a legal value. [-+.?]", "a_gff_strand");
    }
    
    ## column 8: "phase"
    #  legal values are 0 1 2 .
    if ($phase !~ /^[0-2.]$/) {
        record_error("Phase value ($phase) on line $. not a legal value. [012.]", "a_gff_phase");
    }
    
    ## column 9: "attributes"
    my %atth;
    if ($attributes !~ /^\.$/) {
        ## if attributes ends with a semi-colon (and spaces), remove it
        if ($attributes =~ /^(.+)\;\s*$/) {
            $attributes = $1;
        }
    
        ## build a hash of the attributes
        for my $pair ( split(/\;/, $attributes) ) {
            ## we have key=value pairs
            if ($pair =~ /^\s*(.+)=(.+)\s*$/) {
                ## load the attributes hash with these values
                push @{$atth{$1}},split(',', $2)
                
            } else {
                record_error("Attribute ($pair) in column 9, line $. does not seem to be in key=value format", "a_doa");
            }
        }
        
        ## check these attributes
        for my $att (keys %atth) {
            ## make sure this att is within our attribute ontology
            if ($checking_attribute_ontology) {
                my $term = fetch_term($att, $att_ont);

                if (! $term) {
                    record_error("Term ($att) not defined in attribute ontology on line $.", "a_gff_ont", "a_gff_ont");
                }
            }
            
            ## here we can check the values of specific attributes
            if ($att eq 'ID') {
                ## if this is a reference sequence definition, remember it
                if ($type =~ /contig|SO:0000149|supercontig|SO:0000148/) {
                    $seqids_defined{ $atth{$att}[0] } = 1;
                }
            
            } elsif ($att eq 'Name') {
                ## a check here later?
            } elsif ($att eq 'Alias') {
                ## a check here later?
            } elsif ($att eq 'Parent') {
                ## are we checking ontology?
                if ($checking_feature_ontology) {
                    ## throw an error if the Parent attribute does not reference the
                    #   ID attribute of a valid ancestor feature.
                    ## TODO
                }

            } elsif ($att eq 'Target') {
                ## a check here later?
            } elsif ($att eq 'Gap') {
                ## a check here later?
            } elsif ($att eq 'Note') {
                ## a check here later?
            } elsif ($att eq 'Dbxref') {
                ## dbxrefs are in the format of db:id
                for my $val ( @{$atth{$att}} ) {
                    if ($val !~ /.+\:.+/) {
                        if ($val !~ /\%3A/) {
                          record_error("Dbxref ($val) does not conform to the format db:id", "a_dbxref");
                        }
                        else {
                          record_error("Dbxref ($val) in incorrect format, apparently : is replaced with \%3A", "a_dbxref"); 
                        }
                    }
                }
            
            } elsif ($att eq 'Ontology_term') {
                ## make sure the ontology term is given in a valid format ( AA:0000000 )
                for my $val ( @{$atth{$att}} ) {
                    if ($val !~ /^[A-Z]{2}\:[0-9]{7}$/) {
                        record_error("Ontology_term ($val) on line $. does not conform to the format AA:0000000", "a_ontterm");
                    }
                }

            } elsif ($att eq 'ec_number') {
                ## required format for ec_number is N.N.N.N
                for my $val ( @{$atth{$att}} ) {
                    if ($val !~ /.+\..+\..+\..+/) {
                        record_error("ec_number ($val) on line $. does not conform to the format N.N.N.N", "a_ecnum");
                    }
                }
            
            } elsif ($att eq 'molecule_type') {
                ## legal values are dsDNA | ssDNA | dsRNA | ssRNA
                for my $val ( @{$atth{$att}} ) {
                    if ($val !~ /dsDNA|ssDNA|dsRNA|ssRNA/) {
                        record_error("Invalid molecule_type ($val) on line $..  (dsDNA | ssDNA | dsRNA | ssRNA)", "a_moltype");
                    }
                }
            
            } elsif ($att eq 'topology') {
                ## legal values are linear | circular
                for my $val ( @{$atth{$att}} ) {
                    if ($val !~ /linear|circular/) {
                        record_error("Invalid topology ($val) on line $..  (linear | circular)", "a_topology");
                    }
                }

            } elsif ($att eq 'localization') {
                ## legal values are mitochondrion | plastid | episome | plasmid | nuclear | chromosomal
                for my $val ( @{$atth{$att}} ) {
                    if ($val !~ /mitochondrion|plastid|episome|plasmid|nuclear|chromosomal/) {
                        record_error("Invalid localization ($val) on line $..  use (mitochondrion | plastid | episome | plasmid | nuclear | chromosomal)", "a_localization");
                    }
                }
            }
        }
        
        ## if the feature type is mRNA, check the required attributes
        if ($type eq 'mrna' || $type eq 'SO:0000234') {
            for ( 'ID', 'Name', 'Dbxref', 'description' ) {
                if (! defined $atth{$_}) {
                    record_error("$_ is a required attribute for mRNA feature types (line $.)", "a_transcript");
                }
            }
        }
        
        ## if the feature type is contig or supercontig, check the required attributes
        if ($type =~ /contig|SO:0000149|supercontig|SO:0000148/) {
            ## make sure these are defined
            for ( 'ID', 'Name', 'molecule_type', 'organism_name', 'translation_table', 'Dbxref' ) {
                if (! defined $atth{$_}) {
                    record_error("$_ is a required attribute for reference sequence types (line $.)", "a_sod");
                }
            }
            
            ## make sure Dbxref points to a taxon id
            my $taxonid_found = 0;
            for my $val ( @{$atth{'Dbxref'}} ) {
                if ( $val =~ /taxon\:/ ) {
                    $taxonid_found = 1;
                }
            }
            
            if (! $taxonid_found) {
                record_error("Dbxref must include a taxon reference when defining a reference sequences on line $.", "a_sod");
            }
        }
    }
}

## make sure each of the seqids referenced have a corresponding FASTA sequence and that
## they were defined in the document
for my $seqid (keys %seqids_referenced) {
    if (! defined $seqs{$seqid}) {
        record_error("A seqid ($seqid) is referenced in file but no corresponding FASTA sequence is provided", "a_fasta");
    }
    
    if (! defined $seqids_defined{$seqid} ) {
        #this links to section in The Canonical Gene, also check for at least one mRNA, exon, CDS
        record_error("A seqid ($seqid) is referenced but not defined in this file", "a_sod");
    }
}


## warn if either attribute or feature ontologies were not given
unless ($checking_feature_ontology) {
    print "<a href='$url#a_gff_ont'>WARNING:</a> feature ontologies were not checked.\n";
}

if (! $checking_attribute_ontology) {
    print "<a href='$url#a_gff_ont'>WARNING:</a> attribute ontologies were not checked.\n";
}

## check to ensure that we saw the required types: mRNA, cds, exon
unless ($type_count{'gene'} && $type_count{'cds'}) {
    my @missing_types;
    foreach ('gene', 'cds') {
      unless ($type_count{$_}) {
        push(@missing_types, $_);
      }
    }
    record_error("Required type(s) not present in the file: @missing_types", "a_canon");
}


## return code is the number of errors found in the document (not warnings)
exit($error_count);

sub fetch_term {
    my ($id, $ont) = @_;
    my $term = 0;
   
    ## is this an id like GO:0000245 ?
    if ($id =~ /^[A-Z]{2}\:[0-9]{7}$/) {
        $term = $ont->get_term($id);
    } else {
        ## try the name
        $term = $ont->get_term_by_name($id);
    }
    
    return $term;
}

sub load_ontology {
    my $uri = shift;
    my $ont;
    my $status = 0;

    my $parser = new GO::Parser( {handler => 'obj' } );
    
    ## is it a urn?
    if ($uri =~ /^urn/i) {
        print "WARNING: I can't parse URN's yet.  Skipping ontology checking";

    ## is it a local file?
    } elsif (-f $uri) {
        $parser->parse($uri);
        $ont = $parser->handler->graph;
        $status = 1;
    
    ## is it a url?  (will fetch FTP or HTTP)
    } elsif ( $uri =~ m|\://| ) {
        my $returncode = getstore($uri, "/tmp/$$.gff3Valid.temp.obo");
        
        ## return code will be either 500 (http) or 200 (ftp) upon success
        if ($returncode == 500 || $returncode == 200) {
            $parser->parse("/tmp/$$.gff3Valid.temp.obo");
            $ont = $parser->handler->graph;
            $status = 1;
        }
        
    } else {
        print "WARNING: Don't know what to do with the following URI: $uri\n";
    }
    
    return ($ont, $status);
}

sub record_error {
    my $msg = shift;
    my $anchor = shift;

    ++$error_count;
    
    if ($anchor) {
	print "<a href='$url#$anchor'>ERROR:</a> $msg\n";
    }
    else {
	print "ERROR: $msg\n";
    }
}
