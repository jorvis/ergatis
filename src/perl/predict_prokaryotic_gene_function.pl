#!/usr/bin/perl

=head1 NAME

predict_gene_function.pl - Reads annotation evidence for a gene and generates a functional
prediction with product name, gene symbol, EC number and GO terms where possible.

=head1 SYNOPSIS

USAGE: predict_gene_function.pl 
            --input_list=/path/to/bsml.list
            --output_directory=/some/dir
          [ --hmm_list=/path/to/hmm_bsml.list
            --hmm_info=/path/to/hmminfo.db
            --ber_list=/path/to/ber_bsml.list
            --ber_info=/path/to/berinfo.db
            --annotate_on=transcript ]

=head1 OPTIONS

B<--input_list>
    BSML input list, usually the output of the auto_gene_curation Ergatis component.

B<--output_directory>
    This is where the output BSML containing functional predictions will be written.
    
B<--hmm_list>
    Optional. Ouput list of HMM evidence in BSML format.

B<--hmm_info>
    Optional (unless --hmm_list defined) MLDBM file (tied hash) on disk that stores 
    attributes of HMM entries in your search database.  See the hmmlib_to_mldbm.pl 
    script for more info.

B<--ber_list>
    Optional. Output list of BER evidence in BSML format.

B<--ber_info>
    Optional (unless --ber_list is defined).  MLDBM file (tied hash) on disk that stores 
    attributes of matched BER accessions from your search.  See the tchar_to_mldbm.pl
    script for more info.

B<--annotate_on>
    Optional. On which feature type should annotation be assigned?  (default = transcript)
    
B<--log,-l> 
    Log file of all steps in the prediction process.  You should probably use this.

B<--help,-h>
    This help message

=head1  DESCRIPTION

A relatively simple approach to functional prediction, this script considers evidence in
the following order:

   1. Equivalog HMM above trusted cutoff
   2. Characterized BER match with at least 35% identity (configurable) and over 80% of the query length (configurable)
   3. Full-length subfamily HMM (will append 'family protein')
   4. Equivalog-domain HMM
   5. Superfamily HMM (will append 'family protein')
   6. Subfamily domain (will append 'domain protein')
   7. Domain HMM (will append 'domain protein')
   8. PFAM isology (will append 'family protein')
   9. consider tmhmm - if 50% (configurable) of protein is involved in tmhmm matches gene product will be 'putative membrane protein'
  10. consider lipoprotein - If matches, gene product will be 'putative lipoprotein'
  11. Hypothetical equivalog HMM (will name 'conserved hypothetical protein')

=head1  INPUT

Throughout each of these BSML description sections I've removed elements from the examples
that aren't required for this script.  Other elements are optional and should be ignored

=head2 Input BSML format

The input BSML is first scanned so that relationships between features can be defined.  Most
of the computes are run on polypeptides, so if we're attaching annotations to CDS we need
to be able to define which polypeptide corresponds to which CDS.  The section of the BSML
that defines these relationships looks like this:

    <Feature-group group-set="ctha1.gene.16128.1" id="Bsml958">
        <Feature-group-member featref="ctha1.gene.18117.1" feature-type="gene"></Feature-group-member>
        <Feature-group-member featref="ctha1.polypeptide.17997.1" feature-type="polypeptide"></Feature-group-member>
        <Feature-group-member featref="ctha1.CDS.18117.1" feature-type="CDS"></Feature-group-member>
        <Feature-group-member featref="ctha1.transcript.17997.1" feature-type="transcript"></Feature-group-member>
        <Feature-group-member featref="ctha1.exon.18117.1" feature-type="exon"></Feature-group-member>
    </Feature-group>

Once these and the evidence are processed, annotation is added to the <Feature> elements.  (See OUTPUT section)

=head2 HMM BSML format

The only section of the HMM BSML input considered are the alignments.  Each looks like this:

    <Seq-pair-alignment refseq="ctha1.polypeptide.20147.1" refstart="0" compseq="TIGR00633" class="match">
        <Attribute name="total_e_value" content="4.3e-53"></Attribute>
        <Attribute name="total_score" content="187.6"></Attribute>
        <Seq-pair-run compcomplement="0" runlength="238" refpos="0" comprunlength="279" refcomplement="0" runscore="187.6" comppos="0">
            <Attribute name="class" content="match_part"></Attribute>
            <Attribute name="domain_num" content="1"></Attribute>
            <Attribute name="domain_of" content="1"></Attribute>
            <Attribute name="e_value" content="4.3e-53"></Attribute>
        </Seq-pair-run>
        <Link rel="analysis" href="#hmmpfam_analysis" role="computed_by"></Link>
    </Seq-pair-alignment>

=head2 BER BSML format

Not yet written

=head1  OUTPUT

The output of the script is a copy of the input BSML but with additional function definitions
added.  What was this, unannotated:

    <Feature class="transcript" title="ctha1.transcript.1397901.1" id="ctha1.transcript.1397901.1">
      <Interval-loc complement="0" endpos="1378390" startpos="1378195"></Interval-loc>
    </Feature>

Will become something like this, depending on how much evidence is present:

    <Feature class="transcript" title="ctha1.transcript.1397901.1" id="ctha1.transcript.1397901.1">
      <Interval-loc complement="0" endpos="1378390" startpos="1378195"></Interval-loc>
      <Attribute name="gene_product_name" content="polyphosphate kinase"></Attribute>
      <Attribute name="gene_product_name_source" content="PF02503"></Attribute>
      <Attribute name="gene_symbol" content="ppk"></Attribute>
      <Attribute name="gene_symbol_source" content="PF02503"></Attribute>
      <Link rel="analysis" href="#predict_prokaryotic_gene_function_analysis" role="computed_by"></Link>
      <Attribute-list>
        <Attribute name="GO" content="GO:0006799"></Attribute>
        <Attribute name="IEA" content="PF02503"></Attribute>
      </Attribute-list>
      <Attribute-list>
        <Attribute name="GO" content="GO:0008976"></Attribute>
        <Attribute name="IEA" content="PF02503"></Attribute>
      </Attribute-list>
      <Attribute-list>
        <Attribute name="GO" content="GO:0009358"></Attribute>
        <Attribute name="IEA" content="PF02503"></Attribute>
      </Attribute-list>
      <Attribute-list>
        <Attribute name="EC" content="2.7.4.1"></Attribute>
        <Attribute name="IEA" content="PF02503"></Attribute>
      </Attribute-list>
    </Feature>

The log file, if generated with the --log option, contains a great deal of information about the
decision process for functional curation of each gene.

=head1 TESTING

This is my command for testing.  This should be removed eventually, but I'm putting it here
because it takes a while to build the command.  :)

perl -I /usr/local/projects/ergatis/package-devel/lib/perl5/site_perl/5.8.8/ predict_prokaryotic_gene_function.pl --input_list=/usr/local/projects/aengine/output_repository/auto_gene_curation/537_overlap_analysis/auto_gene_curation.bsml.list.all --output_directory=/tmp/annotate --hmm_list=/usr/local/projects/aengine/output_repository/hmmpfam/537_pre_overlap_analysis/hmmpfam.bsml.list --hmm_info=/usr/local/projects/db/coding_hmm/coding_hmm.lib.db --ber_list=/usr/local/projects/aengine/output_repository/ber/537_pre_overlap_analysis/ber.bsml.list --ber_info=/usr/local/projects/db/tchar/tchar.db --log=/tmp/annotate/annotation.log

=head1 DEVEL NOTES

Do we have a problem when features are versioned up and we point at BLAST or HMM results from
the previous version?  This needs to be considered, but since a version change usually implies
a sequence change the old result files would be stale.

=head1  CONTACT

    Joshua Orvis
    jorvis@users.sf.net

=cut

use strict;
use File::Basename;
use Fcntl qw( O_RDONLY );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use MLDBM 'DB_File';
use Pod::Usage;
use XML::Twig;

my %options = ();
my $results = GetOptions (\%options, 
                          'input_list=s',
                          'output_directory=s',
                          'hmm_list=s',
                          'hmm_info=s',
                          'ber_list=s',
                          'ber_info=s',
                          'annotate_on=s',
                          'log|l=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">$options{log}") || die "can't create log file: $!";
}

## make sure we can tie the hmminfo file
my %hmm_info;
if ( $options{hmm_info} ) {
    _log("INFO: attempting to read tied hmm_info hash");
    tie(%hmm_info, 'MLDBM', $options{hmm_info}, O_RDONLY ) or die("Unable to tie hash to $options{hmm_info}");
    _log("INFO: successfully tied hmm_info hash");
}

## make sure we can tie the berinfo file
my %ber_info;
if ( $options{ber_info} ) {
    _log("INFO: attempting to read tied ber_info hash");
    tie(%ber_info, 'MLDBM', $options{ber_info}, O_RDONLY ) or die("Unable to tie hash to $options{ber_info}");
    _log("INFO: successfully tied ber_info hash");
}

my $initial_product_name = 'hypothetical protein';

## priority system for the annotation rules we want to apply
##  you can simply reorder the values of these to reprioritize, but
##  adding (or removing) priority types requires code changes.
my $annot_priorities = {
                            'hmm_equivalog' => 1,
                            'ber_characterized' => 2,
                            'hmm_subfamily' => 3,
                            'hmm_equivalog_domain' => 4,
                            'hmm_superfamily' => 5,
                            'hmm_subfamily' => 6,
                            'hmm_domain' => 7,
                            'hmm_pfam' => 8,
                            'tmhmm' => 9,                   ## not currently used
                            'lipoprotein' => 10,            ## not currently used
                            'hmm_hypoth_equivalog' => 11,
                       };

## this could be a tied hash later.  the _score names refer to the value in
##  %annot_priorities that led to the assignment.  this is how we know whether
##  to override with new evidence.
## structure like:
##      $h{ $polypeptide_id } => {
##                                  product => $initial_product_name,
##                                  product_score => 2,
##                                  product_source => 'PF02503',
##                                  gene_sym => $gene_symbol,
##                                  gene_sym_score => 2,
##                                  gene_sym_source => 'PF02503',
##                                  ec_num => '1.2.3.4',
##                                  ec_num_score => 2,
##                                  ec_num_source => 'PF02503',
##                                  go => [ { term => 'GO:0003735', src  => 'PF02503' }, ... ]
##                               }
my %features;

## indexed on annotated feature IDs, this gives the corresponding polypeptide ID
## structure like:
##      $h{ $annotate_on_id } => $polypeptide_id;
my %peplookup;

## index on CDS feature IDs, this gives the corresponding polypeptide ID
## structure like:
##      $h{ $annotate_on_id } => $polypeptide_id;
my %cdslookup;

my $input_files = parse_multi_list( $options{input_list} );
_log("INFO: found " . scalar(@$input_files) . " input files to process");

## process the feature-groups to get feature relationships
&process_feature_relationships( $input_files );

## apply hmm evidence
my $hmm_files = parse_multi_list( $options{hmm_list} );
_log("INFO: found " . scalar(@$hmm_files) . " HMM files to process");
&process_hmm_results( $hmm_files );

## apply ber evidence
my $ber_files = parse_multi_list( $options{ber_list} );
_log("INFO: found " . scalar(@$hmm_files) . " BER files to process");
&process_ber_results( $ber_files );

## perform post-processing on all predicted names
&check_gene_product_names( \%features );

## write out result BSML
&output_annotated_bsml( $input_files );

untie(%hmm_info) if $options{hmm_info};
untie(%ber_info) if $options{ber_info};

exit(0);


## this is called when we finally start adding annotation to the output BSML.
#   it reads a single Feature element and, if it is of the correct type ($annotate_on)
#   child elements representing the functional annotation are added.
sub annotate_feature_element {
    my ($t, $feat, $ofh) = @_;
    
    ## are we annotating on this feature type?
    if ( $feat->{att}->{class} eq $options{annotate_on} ) {
        _log("INFO: exporting annotation to feature $$feat{att}{id}");
        
        if ( exists $peplookup{ $$feat{att}{id} } ) {
            my $polypeptide_id = $peplookup{ $$feat{att}{id} };
            
            ## add the common name
            my $gene_product_name_elt = XML::Twig::Elt->new( Attribute => {  
                                                                    name => 'gene_product_name',
                                                                    content => $features{ $polypeptide_id }{product},
                                                             } );
            $gene_product_name_elt->paste( last_child => $feat );
            
            ## attach product name source, if available
            if ( $features{ $polypeptide_id }{product_source} ) {
                my $gene_product_name_source_elt = XML::Twig::Elt->new( Attribute => {  
                                                                        name => 'gene_product_name_source',
                                                                        content => $features{ $polypeptide_id }{product_source},
                                                                 } );
                $gene_product_name_source_elt->paste( last_child => $feat );
            }
            
            ## whitespace seems to creep in from everywhere here ... make sure we don't propagate it.
            
            my $gene_symbol = $features{  $polypeptide_id  }{gene_sym};
            if ( $gene_symbol && $gene_symbol =~ /\S/ ) {
                my $gene_symbol_elt = XML::Twig::Elt->new( Attribute => {  
                                                             name => 'gene_symbol',
                                                             content => $gene_symbol,
                                                      } );
                $gene_symbol_elt->paste( last_child => $feat );
                
                my $gene_symbol_source_elt = XML::Twig::Elt->new( Attribute => {  
                                                             name => 'gene_symbol_source',
                                                             content => $features{  $polypeptide_id  }{gene_sym_source},
                                                      } );
                $gene_symbol_source_elt->paste( last_child => $feat );
            }
            
            ## EC numbers are handled a little differently.  They are represented as an Attribute-list
            my $ec_num = $features{ $polypeptide_id }{ec_num};
            
            ## sometimes the EC number contains labels like 'EC ' ... remove this
            if ( $ec_num =~ /EC\s*(.+)/i ) {
                $ec_num = $1;
            }
            
            if ( $ec_num && $ec_num =~ /\S/ ) {
                ## create the Attribute list
                my $att_list = XML::Twig::Elt->new( 'Attribute-list' );
                
                my $ec_att = XML::Twig::Elt->new( Attribute => { name => 'EC', content => $ec_num } );
                   $ec_att->paste( last_child => $att_list );
                   
                my $method_att = XML::Twig::Elt->new( Attribute => { name => 'IEA', 
                                                                     content => $features{  $peplookup{ $$feat{att}{id} }  }{ec_num_source} } );
                   $method_att->paste( last_child => $att_list );
                
                $att_list->paste( last_child => $feat );
            }
            
            ## any GO terms to add?
            if ( scalar @{ $features{$polypeptide_id}{go} } ) {
                for my $go ( @{ $features{$polypeptide_id}{go} } ) {
                    
                    ## create the Attribute list
                    my $att_list = XML::Twig::Elt->new( 'Attribute-list' );

                    my $go_att = XML::Twig::Elt->new( Attribute => { name => 'GO', content => $go->{term} } );
                       $go_att->paste( last_child => $att_list );

                    my $ev_att = XML::Twig::Elt->new( Attribute => { name => 'IEA', 
                                                                     content => $go->{src} },
                                                    );
                       $ev_att->paste( last_child => $att_list );

                    $att_list->paste( last_child => $feat );
                    
                }
            }
            
            ## add a link to the analysis (this one).  the link may not resolve yet, but it will after
            #  the next step in the pipeline, store_config_params, is run.  this should be part of the
            #  Ergatis component.
            my $link = XML::Twig::Elt->new( Link => { rel => 'analysis', 
                                                      href => '#predict_prokaryotic_gene_function_analysis',
                                                      role => 'computed_by',
                                                    } );
            $link->paste( last_child => $feat );
        }
        
    } else {
        _log("INFO: we're not annotating class $$feat{att}{class} features, skipping $$feat{att}{id}");
    }
    
    $feat->flush($ofh, pretty_print => 'indented' );
}


## applies the list of name change rules defined by manual curation team
#   to all gene product names
sub check_gene_product_names {
    my $feats = shift;
    
    for my $feature (keys %$feats) {
        if ( defined $$feats{$feature}{product} ) {
            
            ## names with 'family, family' or 'family family' truncated to just 'family'
            #   ex: AcrB/AcrD/AcrF family family protein -> AcrB/AcrD/AcrF family protein
            $$feats{$feature}{product} =~ s|family[, ]+family|family|ig;
            
            ## names with 'domain, domain' or 'domain domain' truncated to just 'domain'
            $$feats{$feature}{product} =~ s|domain[, ]+domain|domain|ig;
            
            ## Names with family and domain should just be domain protein:
            #   ex. PHP domain family protein -> PHP domain protein. 
            $$feats{$feature}{product} =~ s|domain family protein|domain protein|ig;
            
            ## if 'putative' and 'family' are in the same name, remove the putative
            #   ex. putative methyltransferase family protein -> methyltransferase family protein
            if ( $$feats{$feature}{product} =~ m|putative|i && $$feats{$feature}{product} =~ m|family|i ) {
                $$feats{$feature}{product} =~ s|putative\s*||ig;
            }
            
            ## Remove The from the beginning of all names. 
            #   ex. The GLUG motif family protein -> GLUG motif family protein
            $$feats{$feature}{product} =~ s|^the\s+||ig;
            
            ## PAS domain proteins should just be called 'sensory box protein'
            if ( $$feats{$feature}{product} =~ m|PAS domain|i ) {
                $$feats{$feature}{product} = 'sensory box protein';
            }
            
            ## S-box domain proteins should just be called 'sensory box protein'
            if ( $$feats{$feature}{product} =~ m|S\-box domain protein|i ) {
                $$feats{$feature}{product} = 'sensory box protein';
            }

            ## if 'response regulator' is anywhere within the name, that should be the entire name
            if ( $$feats{$feature}{product} =~ m|response regulator|i ) {
                $$feats{$feature}{product} = 'response regulator';
            }            
            
            ## any name with conserved hypothetical, DUF family, or protein of unknown function should be 
            #   conserved hypothetical protein 
            #   exs conserved hypothetical family protein; Protein of unknown function (DUF454) family protein -> conserved hypothetical protein
            if ( $$feats{$feature}{product} =~ m|conserved hypothetical|i || 
                 $$feats{$feature}{product} =~ m|DUF.*protein| ||
                 $$feats{$feature}{product} =~ m|protein.*unknown function|i   ) {
                $$feats{$feature}{product} = 'conserved hypothetical protein';
            }
            
            ## anything with 'uncharacterized ACR' should be changed to conserved hypothetical
            if ( $$feats{$feature}{product} =~ m|ncharacterized ACR| ) {
                $$feats{$feature}{product} = 'conserved hypothetical protein';
            }
            
            ## any gene symbol family starting with Y should be a conserved hypothetical
            #   ex. YfiH family COG1496 family protein -> conserved hypothetical protein
            if ( $$feats{$feature}{product} =~ m|Y\S*[A-Z] family| ) {
                $$feats{$feature}{product} = 'conserved hypothetical protein';
            }
            
            ## names with superfamily and family should just be family
            #   ex. PAP2 superfamily family protein -> PAP2 family protein
            $$feats{$feature}{product} =~ s|superfamily family|family|ig;
            
            ## Truncate "family protein" when if follows subunit. 
            #   ex. electron transport complex, RnfABCDGE type, D subunit family protein -> electron transport complex, RnfABCDGE type, D subunit
            $$feats{$feature}{product} =~ s|subunit family protein|subunit|i;
            
            ## anytime there's a protein something protein, take off the first instance of protein
            $$feats{$feature}{product} =~ s|protein (\S+ protein)|$1|i;
            
            ## replace 'or' with /. 
            #   ex. succinate dehydrogenase or fumarate reductase, flavoprotein subunit family protein -> succinate dehydrogenase/fumarate reductase, flavoprotein subunit
            $$feats{$feature}{product} =~ s| or |\/|i;
            
            ## remove leading predicted from names
            #   ex. Predicted permease family protein -> permease family protein
            $$feats{$feature}{product} =~ s|^predicted\s+||i;
        }
    }
}

## reads through input files and, for each, creates a new one in the output_directory
##  with the same base name but with annotations attached to $annotate_on features.
sub output_annotated_bsml {
    my $infiles = shift;
    
    for my $infile ( @$infiles ) {
        my $base = fileparse( $infile, '.bsml' );
        
        ## create an output file for this
        open(my $ofh, ">$options{output_directory}/$base.bsml") || die "failed to create output file: $!";
        
        my $t = XML::Twig->new(
                    twig_handlers => {
                        Feature => sub { annotate_feature_element($_[0], $_[1], $ofh) },
                    },
                );
                
        $t->parsefile($infile);
        $t->flush($ofh, pretty_print => 'indented' );
    }
}

## pulls the <Feature-group> elements from the input BSML files and passes
#  them to &process_input_feature_group
sub process_feature_relationships {
    my $files = shift;
    
    for my $file ( @$files ) {
        _log("INFO: processing feature relationships in file: $file");
        
        my $t = XML::Twig->new(
                    twig_roots => {
                        'Feature-group' => \&process_input_feature_group,
                    }
                );
        $t->parsefile( $file );
    }
}

sub process_ber_alignment {
    my ($t, $spa, $id_lookup, $input_lengths) = @_;
    
    my $ref_id = $spa->{att}->{refseq};
    my $comp_id = $spa->{att}->{compseq};
    
    ## annotation is mapped internally on the polypeptide, so we need that ID
    my $polypeptide_id;
    if ( exists $cdslookup{$ref_id} ) {
        $polypeptide_id = $cdslookup{$ref_id};
    } else {
        _log("WARN: CDS -> polypeptide lookup failed for ref id $ref_id, skipping this BER hit");
        return;
    }
    
    ## is it in the characterized db?
    if ( exists $ber_info{ $$id_lookup{$comp_id} } ) {
        
        ## we got here, so this is a characterized BER match.
        
        ## we can just compare with the product_score here since all scores are currently
        #   set together.  in the future this could change.
        if ( $$annot_priorities{"ber_characterized"} < $features{$polypeptide_id}{product_score} ) {
            _log("INFO: $polypeptide_id: setting product, gene_sym and ec_num based on $ref_id BER hit to $$id_lookup{$comp_id}");
            
            $features{$polypeptide_id}{gene_sym} = $ber_info{ $$id_lookup{$comp_id} }{gene_symbol} || '';
            $features{$polypeptide_id}{gene_sym_score} = $$annot_priorities{"ber_characterized"};
            $features{$polypeptide_id}{gene_sym_source} = $$id_lookup{$comp_id};

            $features{$polypeptide_id}{ec_num} = $ber_info{ $$id_lookup{$comp_id} }{ec_num} || '';
            $features{$polypeptide_id}{ec_num_score} = $$annot_priorities{"ber_characterized"};
            $features{$polypeptide_id}{ec_num_source} = $$id_lookup{$comp_id};

            $features{$polypeptide_id}{product} = $ber_info{ $$id_lookup{$comp_id} }{com_name};
            $features{$polypeptide_id}{product_score} = $$annot_priorities{"ber_characterized"};
            $features{$polypeptide_id}{product_source} = $$id_lookup{$comp_id};
            
            ## handle any GO terms on the accession
            if ( scalar @{ $ber_info{$$id_lookup{$comp_id}}{go} } ) {
                for my $term ( @{ $ber_info{$$id_lookup{$comp_id}}{go} } ) {
                    push @{ $features{$polypeptide_id}{go} }, { term => $term, src => $$id_lookup{$comp_id} };
                }

                _log("INFO: $ref_id: set " . scalar(@{ $ber_info{$$id_lookup{$comp_id}}{go} }) . " GO terms based on hit to $$id_lookup{$comp_id}");
            }
            
        } else {
            _log("INFO: $polypeptide_id: found characterized $ref_id BER hit to $$id_lookup{$comp_id}, but higher priority annotation exists ($features{$ref_id}{product_score})");
        }
        
    } else {
        _log("INFO: $polypeptide_id: skipping $ref_id hit to $$id_lookup{$comp_id} - not in characterized database");
    }
}

sub process_ber_results {
    my $files = shift;
    
    ## within each BER result file the polypeptide Sequence elements have IDs like this:
    #
    #       RF_NP_343702.1_15899097_NC_002754
    #
    #   but also has elements like this:
    #
    #       <Attribute name="defline" content="RF|NP_343702.1|15899097|NC_002754"></Attribute>
    #
    #   It's these second values that will be contained within the characterized lookup.  Pass
    #   through the file building the lookup and go through the alignments
    
    for my $file ( @$files ) {
        _log("INFO: processing BER results in $file");
        
        ## structure like:
        #   $h->{RF_NP_343702.1_15899097_NC_002754} = 'RF|NP_343702.1|15899097|NC_002754'
        my $id_lookup = {};
        
        ## holds the lengths of the input sequences (usually only one, but a hash here just in case)
        my $input_lengths = {};
        
        my $t = XML::Twig->new(
            twig_roots => {
                ## we need to parse each of the sequences here to populate $id_lookup
                ##
                ##  each sequence looks like:
                #  <Sequence class="polypeptide" title="holliday junction resolvase ..." id="RF_NP_147250.1_14600729_NC_000854" molecule="aa">
                #    <Attribute name="defline" content="RF|NP_147250.1|14600729|NC_000854"></Attribute>
                #  </Sequence>
                
                'Sequence' => sub {
                    my ($t, $elt) = @_;
                    
                    my $id = $elt->{att}->{id} || die "failed to get ID attribute of Sequence in BER alignment file $file";
                    my $defline;
                    
                    ## check for input sequences.  They will be a CDS and have an input_of role
                    if ( $elt->{att}->{class} eq 'CDS' ) {
                        for my $link ( $elt->children('Link') ) {
                            if ( $link->{att}->{role} eq 'input_of' ) {
                                $$input_lengths{$id} = $elt->{att}->{length} || die "failed to determine sequence length for $id, an input to BER within file $file";
                            }

                            ## we can return here, since these won't have any attributes to save.
                            return;
                        }
                    }
                    
                    for my $att ( $elt->children('Attribute') ) {
                        if ( $att->{att}->{name} eq 'defline' ) {
                        
                            ## if the comp ID looks like this:
                            #       RF|YP_235471.1|66045630|NC_007005
                            ## we only need to look at the first bit, since that's how accessions are stored in the 
                            ## ber_info file:
                            #       RF|YP_235471.1
                            if ( $att->{att}->{content} =~ /^([A-Z]+\|[A-Z0-9_\.]+)/ ) {
                                $$id_lookup{$id} = $1;
                            } else {
                                $$id_lookup{$id} = $att->{att}->{content};
                            }
                            
                            _log("");
                        }
                    }
                },
                'Seq-pair-alignment' => sub {
                    my ($t, $spa) = @_;
                    &process_ber_alignment($t, $spa, $id_lookup, $input_lengths),
                },
            },
        );
        $t->parsefile( $file );
    }
    
}

## pulls the <Seq-pair-alignment> elements from the input HMM files and passes
#   them to process_hmm_coding_alignment, where the actual works happens
sub process_hmm_results {
    my $files = shift;
    
    for my $file ( @$files ) {
        _log("INFO: processing HMM results in $file");
        
        ## these files can have results from multiple searches in them, so we can't
        #  rely on file names or pulling the first Sequence element ID.  We really
        #  only need to look at refseq values within the alignments.
        my $t = XML::Twig->new(
                    twig_roots => {
                        'Seq-pair-alignment' => \&process_hmm_coding_alignment,
                    },
                );
        $t->parsefile( $file );
    }
}

sub process_hmm_coding_alignment {
    my ($t, $spa) = @_;
    
    my $ref_id = $spa->{att}->{refseq};
    my $comp_id = $spa->{att}->{compseq};
    my $total_score;
    
    ## this ref_id needs to be part of our input set
    if ( exists $features{$ref_id} ) {
        _log("INFO: $ref_id: processing HMM match to $comp_id");
    } else {
        _log("WARN: $ref_id: <- not recognized as part of input set to be annotated.  skipping.");
        return;
    }
    
    for my $att ( $spa->children('Attribute') ) {
        if ( $att->{att}->{name} eq 'total_score' ) {
            $total_score = $att->{att}->{content};
        }
    }
    
    ## we need both IDs and the score to continue
    unless ( $ref_id && $comp_id && defined $total_score ) {
        die "failed to get ref_id, comp_id and total_score from Seq-pair-alignment";
    }
    
    ## make sure we have info on this HMM
    if ( exists $hmm_info{$comp_id} ) {
        
        if ( $total_score >= $hmm_info{$comp_id}{trusted_cutoff} ) {
            
            ## what is the priority of this match, based on its type?
            my $isotype = $hmm_info{$comp_id}{isotype};
            if ( exists $$annot_priorities{ "hmm_$isotype" } ) {
                
                ## we can just compare with the product_score here since all scores are currently
                #   set together.  in the future this could change.
                if ( $$annot_priorities{"hmm_$isotype"} < $features{$ref_id}{product_score} ) {
                    _log("INFO: $ref_id: setting product, gene_sym and ec_num based on $isotype HMM hit to $comp_id ($isotype)");
                    
                    $features{$ref_id}{gene_sym} = $hmm_info{$comp_id}{gene_symbol};
                    $features{$ref_id}{gene_sym_score} = $$annot_priorities{"hmm_$isotype"};
                    $features{$ref_id}{gene_sym_source} = $comp_id;
                    
                    $features{$ref_id}{ec_num} = $hmm_info{$comp_id}{ec_num};
                    $features{$ref_id}{ec_num_score} = $$annot_priorities{"hmm_$isotype"};
                    $features{$ref_id}{ec_num_source} = $comp_id;
                    
                    $features{$ref_id}{product_score} = $$annot_priorities{"hmm_$isotype"};
                    $features{$ref_id}{product_source} = $comp_id;
                    
                    ## what we do now depends on the isotype
                    if ( $isotype eq 'equivalog' ) {
                        $features{$ref_id}{product} = $hmm_info{$comp_id}{hmm_com_name};
                    
                    } elsif ( $isotype eq 'subfamily' ) {
                        $features{$ref_id}{product} = $hmm_info{$comp_id}{hmm_com_name} . ' family protein';
                    
                    } elsif ( $isotype eq 'equivalog_domain' ) {
                        $features{$ref_id}{product} = $hmm_info{$comp_id}{hmm_com_name};
                    
                    } elsif ( $isotype eq 'superfamily' ) {
                        $features{$ref_id}{product} = $hmm_info{$comp_id}{hmm_com_name} . ' family protein';
                    
                    } elsif ( $isotype eq 'subfamily' ) {
                        $features{$ref_id}{product} = $hmm_info{$comp_id}{hmm_com_name} . ' domain protein';
                    
                    } elsif ( $isotype eq 'domain' ) {
                        $features{$ref_id}{product} = $hmm_info{$comp_id}{hmm_com_name} . ' domain protein';
                    
                    } elsif ( $isotype eq 'pfam' ) {
                        $features{$ref_id}{product} = $hmm_info{$comp_id}{hmm_com_name} . ' family protein';
                    
                    } elsif ( $isotype eq 'hypoth_equivalog' ) {
                        $features{$ref_id}{product} = 'conserved hypothetical protein';
                    
                    }
                    
                    ## handle any GO terms on the accession
                    if ( scalar @{ $hmm_info{$comp_id}{go} } ) {
                        for my $term ( @{ $hmm_info{$comp_id}{go} } ) {
                            push @{ $features{$ref_id}{go} }, { term => $term, src => $comp_id };
                        }
                    
                        _log("INFO: $ref_id: set " . scalar(@{ $hmm_info{$comp_id}{go} }) . " GO terms based on $isotype HMM hit to $comp_id");
                    }
                    
                } else {
                    _log("INFO: $ref_id: $isotype HMM match to $comp_id not better priority than current match.  skipping");
                }
                
            } else {
                _log("WARN: $ref_id: ignoring trusted HMM hit to $comp_id - can't tell which isotype $comp_id is");
            }
        } else {
            _log("INFO: $ref_id: ignoring HMM hit to $comp_id - score ($total_score) >= trusted cutoff ($hmm_info{$comp_id}{trusted_cutoff}) ");
        }
        
    } else {
        _log("WARN: $ref_id: no entry for HMM hit $comp_id found in hmm_info file.  ignored.");
    }
}


## given a Feature-group element this searches the Feature-group-member elements to find
#   the ID relationship between the polypeptide and the feature being annotated
sub process_input_feature_group {
    my ($t, $feat_group) = @_;
    
    ## the gene organizes the feature group
    my $gene_name = $feat_group->{att}->{'group-set'};
    _log("INFO: processing features within gene: $gene_name");
    
    ## feat_group here is an XML::Twig::Elt
    my @feat_group_members = $feat_group->getElementsByTagName('Feature-group-member');
    my ($polypeptide_id, $cds_id, $annot_feat_id);
    
    for my $gm ( @feat_group_members ) {
        if ( $gm->{att}->{'feature-type'} eq 'polypeptide' ) {
            $polypeptide_id = $gm->{att}->{featref};
        }
        
        if ( $gm->{att}->{'feature-type'} eq 'CDS' ) {
            $cds_id = $gm->{att}->{featref};
        }

        if ( $gm->{att}->{'feature-type'} eq $options{annotate_on} ) {
            $annot_feat_id = $gm->{att}->{featref};
        }
    }
    
    ## both must have been found
    if ( defined $polypeptide_id && defined $annot_feat_id ) {
        _log("INFO: annotating $polypeptide_id evidence from gene $gene_name on $annot_feat_id");
    } else {
        die "FATAL: failed to find polypeptide and $options{annotate_on} feature types within gene $gene_name Feature-group";
    }
    
    ## remember the relationship
    $peplookup{$annot_feat_id} = $polypeptide_id;
    
    _log("INFO: $cds_id: storing relationship to polypeptide ID: $polypeptide_id");
    $cdslookup{$cds_id} = $polypeptide_id;
    
    
    ## initialize it, with crappy score values
    $features{$polypeptide_id} = {
        product => $initial_product_name,
        product_score => 1000,
        gene_sym => '',
        gene_sym_score => 1000,
        ec_num => '',
        ec_num_score => 1000,
        ec_num_source => '',
        go => [],
    };
}

## several of the options can be comma-separated input lists.  this
##  resolves them to an arrayref of paths.
sub parse_multi_list {
    my $list = shift;
    my $files = [];
    
    my @list_files = split(',', $list);
    
    return $files unless scalar @list_files;
    
    for my $list_file ( @list_files ) {
        open(my $ifh, "<$list_file") || die "failed to open list file ($list_file): $!";
        
        while ( <$ifh> ) {
            chomp;
            next if /^\s*$/;
            
            _log("INFO: adding $_ to input file processing list");
            push @$files, $_;
        }
    }
    
    return $files;
}

sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
}


sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( input_list output_directory );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
    
    ##
    ## you can do other things here, such as checking that files exist, etc.
    ##
    
    ## handle some defaults
    $options{annotate_on} = 'transcript' unless $options{annotate_on};
    $options{hmm_list} = '' unless $options{hmm_list};
    $options{ber_list} = '' unless $options{ber_list};
    $options{hmm_info} = '' unless $options{hmm_info};
    
    ## annotate_on must be either transcript, CDS or polypeptide
    if ( $options{annotate_on} ne 'transcript' && 
         $options{annotate_on} ne 'CDS' &&
         $options{annotate_on} ne 'polypeptide' ) {
    
         die "--annotate_on value must be either transcript, CDS or polypeptide";
    }
    
    ## if hmm_list was defined, hmm_info needs to be
    if ( $options{hmm_list} && ! $options{hmm_info} ) {
        die "--hmm_info required if passing --hmm_list";
    }

    ## if ber_list was defined, ber_info needs to be
    if ( $options{ber_list} && ! $options{ber_info} ) {
        die "--ber_info required if passing --ber_list";
    }
}
