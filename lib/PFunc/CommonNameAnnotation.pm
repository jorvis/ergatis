package PFunc::CommonNameAnnotation;

use strict;
use warnings;
use Exporter 'import';
use vars qw(@EXPORT_OK);

@EXPORT_OK = qw(clean_common_name);

sub clean_common_name {
    my ($new_product_name) = @_;
    
    ## names with 'family, family' or 'family family' truncated to just 'family'
    #   ex: AcrB/AcrD/AcrF family family protein -> AcrB/AcrD/AcrF family protein
    $new_product_name =~ s|family[, ]+family|family|ig;
    
    ## names with 'domain, domain' or 'domain domain' truncated to just 'domain'
    $new_product_name =~ s|domain[, ]+domain|domain|ig;
    
    ## Names with family and domain should just be domain protein:
    #   ex. PHP domain family protein -> PHP domain protein. 
    $new_product_name =~ s|domain family protein|domain protein|ig;
    
    ## if 'putative' and 'family' are in the same name, remove the putative
    #   ex. putative methyltransferase family protein -> methyltransferase family protein
    if ( $new_product_name =~ m|putative|i && $new_product_name =~ m|family|i ) {
        $new_product_name =~ s|putative\s*||ig;
    }
    
    ## Remove The from the beginning of all names. 
    #   ex. The GLUG motif family protein -> GLUG motif family protein
    $new_product_name =~ s|^the\s+||ig;
    
    ## PAS domain proteins should just be called 'sensory box protein'
    if ( $new_product_name =~ m|PAS domain|i ) {
        $new_product_name = 'sensory box protein';
    }
    
    ## S-box domain proteins should just be called 'sensory box protein'
    if ( $new_product_name =~ m|S\-box domain protein|i ) {
        $new_product_name = 'sensory box protein';
    }

    ## if 'response regulator' is anywhere within the name, that should be the entire name
    if ( $new_product_name =~ m|response regulator|i ) {
        $new_product_name = 'response regulator';
    }            
    
    ## any name with conserved hypothetical, DUF family, or protein of unknown function should be 
    #   conserved hypothetical protein 
    #   exs conserved hypothetical family protein; Protein of unknown function (DUF454) family protein -> 
    #   conserved hypothetical protein
    if ( $new_product_name =~ m|conserved hypothetical|i || 
         $new_product_name =~ m|DUF.*protein| ||
         $new_product_name =~ m|protein.*unknown function|i   ) {
        $new_product_name = 'conserved hypothetical protein';
    }
    
    ## anything with 'uncharacterized ACR' should be changed to conserved hypothetical
    if ( $new_product_name =~ m|ncharacterized ACR| ) {
        $new_product_name = 'conserved hypothetical protein';
    }
    
    ## any gene symbol family starting with Y should be a conserved hypothetical
    #   ex. YfiH family COG1496 family protein -> conserved hypothetical protein
    if ( $new_product_name =~ m|Y\S*[A-Z] family| ) {
        $new_product_name = 'conserved hypothetical protein';
    }
    
    ## names with superfamily and family should just be family
    #   ex. PAP2 superfamily family protein -> PAP2 family protein
    $new_product_name =~ s|superfamily family|family|ig;
    
    ## Truncate "family protein" when if follows subunit. 
    #   ex. electron transport complex, RnfABCDGE type, D subunit family protein -> 
    #   electron transport complex, RnfABCDGE type, D subunit
    $new_product_name =~ s|subunit family protein|subunit|i;
    
    ## anytime there's a protein something protein, take off the first instance of protein
    $new_product_name =~ s|protein (\S+ protein)|$1|i;
    
    ## replace 'or' with /. 
    #   ex. succinate dehydrogenase or fumarate reductase, flavoprotein subunit family protein -> 
    #   succinate dehydrogenase/fumarate reductase, flavoprotein subunit
    $new_product_name =~ s| or |\/|i;
    
    ## remove leading predicted from names
    #   ex. Predicted permease family protein -> permease family protein
    $new_product_name =~ s|^predicted\s+||i;
    
    ## remove 'precursor' from the protein name
    #   ex. halocyanin precursor -> halocyanin
    $new_product_name =~ s|precursor||i;
    
    ## rename 'probable' to 'putative'
    $new_product_name =~ s|probable|putative|gi;
    
    ## remove leading whitespace
    if ( $new_product_name =~ /^\s*(.+)$/ ) {
        $new_product_name = $1;
    }
    
    ## remove trailing whitespace
    if ( $new_product_name =~ /^(.+?)\s+$/ ) {
        $new_product_name = $1;
    }
    
    ## remove trailing periods
    if ( $new_product_name =~ /^(.+?)\.$/ ) {
        $new_product_name = $1;
    }
    
    ## remove EC numbers from names. this will catch the following case-insensitive convention:
    #   ex. Aspartate aminotransferase (EC 2.6.1.1) -> Aspartate aminotransferase
    $new_product_name =~ s|\(ec .+?\)||ig;
    
    ## The beginnings of all names should be lowercase unless it is an abbreviation. 
    #   ex. Alkaline phosphatase family protein -- alkaline phosphatase family protein
    if ( $new_product_name =~ /^([A-Z])([a-z].+)/ ) {
        $new_product_name = lc($1) . $2;
    }
    
    ## multiple whitespace should be collapsed
    $new_product_name =~ s|\s{2,}| |g;

    return $new_product_name;
}

1;
