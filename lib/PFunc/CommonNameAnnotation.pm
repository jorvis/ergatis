package PFunc::CommonNameAnnotation;

use strict;
use warnings;
use Exporter 'import';
use vars qw(@EXPORT_OK);

@EXPORT_OK = qw(clean_common_name clean_gene_symbol);

sub clean_common_name {
    my ($new_product_name, $chp) = @_;
    
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

    if( $new_product_name eq 'transposase and inactivated derivative' || $new_product_name =~ /^IS/ ) {
        $new_product_name = 'putative transposase';
    }
    
    if( $new_product_name =~ /phospholipase d active site motif family protein/i ) {
        $new_product_name = 'phospholipase D family protein';
    }
    
    if( ($new_product_name =~ /_/) && ($new_product_name eq 'menC_gamma/gm+: o-succinylbenzoic acid (OSB) synthetase')) {
        $new_product_name = 'o-succinylbenzoic acid (OSB) synthetase';
    }

    ## Misc name changes
    if( $new_product_name eq "transcriptional regulatory protein, C terminal family protein" ) {
        $new_product_name = 'putative transcriptional regulator';
    } elsif( $new_product_name eq "bacterial regulatory proteins, luxR family protein" ) {
        $new_product_name = "transcriptional regulator, LuxR family";
    } elsif( $new_product_name eq "bacterial regulatory helix-turn-helix proteins, AraC family protein" ) {
        $new_product_name = "transcriptional regulator, AraC family";
    } elsif( $new_product_name eq "bacterial regulatory proteins, tetR family protein" ) {
        $new_product_name = "transcriptional regulator, TetR family";
    } elsif( $new_product_name eq "bacterial regulatory proteins, gntR family protein" ) {
        $new_product_name = "transcriptional regulator, GntR family";
    } elsif( $new_product_name eq "bacterial regulatory proteins, lacI family protein" ) {
        $new_product_name = "transcriptional regulator, LacI family";
    } elsif( $new_product_name eq "FGGY family of carbohydrate kinases, C-terminal domain protein" ) {
        $new_product_name = "carbohydrate kinase, FGGY family";
    } elsif( $new_product_name eq "tripartite ATP-independent periplasmic transporters, DctQ component family protein" ) {
        $new_product_name = "tripartite ATP-independent periplasmic transporter, DctQ family";
    } elsif( $new_product_name eq "thiamin/thiamin pyrophosphate ABC transporter, thiamin/thiamin pyrophospate-binding protein" ) {
        $new_product_name = "thiamin/thiamine pyrophosphate ABC transporter, thiamin/thiamine pyrophospate-binding protein";
    } elsif( $new_product_name eq "GSPII_E N-terminal domain protein" ) {
        $new_product_name = "bacteriophage N4 adsorption protein B";
    } elsif( $new_product_name eq "transcriptional activator of defense systems" ) {
        $new_product_name = "multiple antibiotic resistance protein MarA";
    } elsif($new_product_name eq "type IV secretory pathway VirD2 components") {
	$new_product_name = "type IV secretory pathway protein";
    } elsif($new_product_name eq "hydro-lases, Fe-S type, tartrate/fumarate subfamily, beta region domain protein") {
	$new_product_name = "fumarate hydratase family protein";
    } elsif($new_product_name eq "glycogen/starch synthases, ADP-glucose type family protein") {
	$new_product_name = "glycogen/starch synthase";
    } elsif($new_product_name eq "glutamate synthases, NADH/NADPH, small subunit domain protein") {
	$new_product_name = "glutamate synthase, NADH/NADPH, small subunit";
    } elsif($new_product_name eq "K+ potassium transporter family protein") {
	$new_product_name = "potassium uptake protein";
    } elsif($new_product_name eq "domain related to MnhB subunit of Na+/H+ antiporter family protein") {
	$new_product_name = "Na+/H+ antiporter family protein";	
    } elsif($new_product_name eq "arginine-tRNA-transferase, C terminus family protein") {
	$new_product_name = "putative arginine-tRNA-transferase";
    } elsif($new_product_name eq "cytochrome b(C-terminal)/b6/petD family protein") {
	$new_product_name = "cytochrome b family protein";
    } elsif($new_product_name eq "traG-like , N-terminal region family protein") {
	$new_product_name = "putative traG protein";
    } elsif($new_product_name eq "cyclic di-GMP binding protein VCA0042") {
	$new_product_name = "cyclic di-GMP binding protein";
    } elsif($new_product_name eq "alr5027 protein") {
	$new_product_name = "heme-binding protein HutZ";
    } elsif($new_product_name eq "putative 2-hydroxyacid dehydrogenase HI_1556") {
	$new_product_name = "putative 2-hydroxyacid dehydrogenase";
    } elsif($new_product_name eq "SULFATE TRANSPORTER SULFATE TRANSPORTER FAMILY PROTEIN") {
	$new_product_name = "sulfate permease family protein";
    } elsif($new_product_name eq "conserved protein with nucleoside triphosphate hydrolase domain") {
	$new_product_name = "putative ATP-dependent endonuclease";
    } elsif($new_product_name eq "gene 25-like lysozyme family protein") {
	$new_product_name = "lysozyme family protein";
    } elsif($new_product_name eq "bordetella uptake gene (bug) product family protein") {
	$new_product_name = "bug family protein";
    } elsif($new_product_name eq "phage/plasmid replication , gene II/X family protein") {
	$new_product_name = "phage/plasmid replication protein, gene II/X family";
    } elsif($new_product_name eq "invasion gene expression up-regulator, SirB family protein") {
	$new_product_name = "invasion gene expression up-regulator";
    } elsif(($new_product_name eq "PIII") || ($new_product_name eq "zn-dependent hydrolase of the beta-lactamase fold")) {
	$new_product_name = $chp;
    }

    ## If protein name begins with 'orf' or 'residues' or 'ttg start' or contains 'similar to' 
    #  or has no significant matches or unnamed or has only 'protein' 
    #  or has short bogus name or has 'gene' number 'protein'
    #  then change to conserved hypothetical protein
    if( $new_product_name =~ /^(orf)/i ||
	$new_product_name =~ /similar to/ ||
	$new_product_name =~ /^(cons|no significant matches)$/i ||
 	$new_product_name =~ /\bunnamed\b/i ||
	$new_product_name =~ /^\s*protein\s*$/i ||
        $new_product_name =~ /^residues\b/i ||
        $new_product_name =~ /^[A-Za-z]{3,}\d{4,}\s+/ ||
        $new_product_name =~ /^\w{1,2}\d{1,3}$/ ||
        $new_product_name =~ /^ttg start/i ||
	$new_product_name =~ /gene \d+ protein/ ||
	$new_product_name =~ /conserved domain protein/ ) {
	$new_product_name = $chp;
    }

    if( $new_product_name =~ /homolog/ ) {
        if( $new_product_name =~ /shiA homolog/ ) {
            $new_product_name = 'shiA protein';
        } elsif( $new_product_name =~ /virulence factor mviM homolog/ ) {
            $new_product_name = 'virulence factor mviM';
        } elsif( $new_product_name =~ /protein phnA homolog/ ) {
            $new_product_name = 'phnA protein';
        } elsif( $new_product_name =~ /protein seqA homolog/ ) {
            $new_product_name = 'seqA protein';
        } else {
            $new_product_name = $chp;
        }
    }

    ## any name with conserved hypothetical, DUF or UPF family, or protein of unknown function should be 
    #   conserved hypothetical protein 
    #   exs conserved hypothetical family protein; Protein of unknown function (DUF454) family protein -> 
    #   conserved hypothetical protein
    if ( 
#	 $new_product_name =~ m|conserved hypothetical|i ||
	 $new_product_name =~ /conserved hypothetica\b/i || 
         $new_product_name =~ m|DUF.*protein| ||
	 $new_product_name =~ /\b(DUF|UPF)/ ||
	 $new_product_name =~ /golgi/i ||
         $new_product_name =~ m|protein.*unknown function|i ) {
         $new_product_name = $chp;
    } elsif ($new_product_name =~ /hypothetica\b/) {
	$new_product_name = 'hypothetical protein';
    } 

    ## anything with 'uncharacterized ACR' should be changed to conserved hypothetical
    if ( $new_product_name =~ /\buncharacteri(z|s)ed\b/i || $new_product_name =~ m/ncharacteri(z|s)ed ACR/ ) {
        $new_product_name = $chp;
    }
    
    ## any gene symbol family starting with Y or beginning with ZB locus should be a conserved hypothetical
    #   ex. YfiH family COG1496 family protein -> conserved hypothetical protein
#    if ( $new_product_name =~ m|Y\S*[A-Z] family| || $new_product_name =~ /^\s*(Z|B)\d+\s+gene\s+product/ ) {
     if ($new_product_name =~ /^\s*(Z|B)\d+\s+gene\s+product/) {
        $new_product_name = $chp;
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
    
    ## Americanize some words
    $new_product_name =~ s/utilisation/utilization/g;
    $new_product_name =~ s/utilising/utilizing/g;
    $new_product_name =~ s/mobilisation/mobilization/g;
    $new_product_name =~ s/dimerisation/dimerization/g;
    $new_product_name =~ s/disulphide/disulfide/g;
    $new_product_name =~ s/sulphur/sulfur/g;
    if( $new_product_name =~ /\bhaem(\w+)/ ) {
        my $p = $1;
        $new_product_name =~ s/haem$p/hem$p/;
    }
    
    $new_product_name =~ s/\d*\s*C-terminus//;
    $new_product_name =~ s/\d*\s*N-terminus//;

    ## Make the words singular
    $new_product_name =~ s/desulfurases/desulfurase/g;
    $new_product_name =~ s/synthases/synthase/g;
    $new_product_name =~ s/kinases/kinase/g;
    $new_product_name =~ s/decarboxylases/decarboxylase/g;
    $new_product_name =~ s/oxidases/oxidase/g;
    $new_product_name =~ s/\bgenes\b/gene/g;
    ## Spelling corrections
    $new_product_name =~ s/hypotheical|hypothetic\b/hypothetical/g;
    
    ## replace 'or' with /. 
    #   ex. succinate dehydrogenase or fumarate reductase, flavoprotein subunit family protein -> 
    #   succinate dehydrogenase/fumarate reductase, flavoprotein subunit
    $new_product_name =~ s| or |\/|i;
    
    ## Change leading predicted, possible, potential, probable from names to putative
    #   ex. Predicted permease family protein -> putative permease family protein
    $new_product_name =~ s/^(predicted|possible|potential|probable)/putative/i;
    
    ## remove 'precursor' from the protein name
    #   ex. halocyanin precursor -> halocyanin
    $new_product_name =~ s|precursor||i;
    
    ## rename 'probable' or 'putaive' to 'putative'
    $new_product_name =~ s/probable|putaive/putative/gi;
    
    ## Correct asparate to aspartate
    $new_product_name =~ s/\basparate\b/aspartate/;

    ## remove leading whitespace
    if ( $new_product_name =~ /^\s*(.+)$/ ) {
        $new_product_name = $1;
    }
    
    ## remove trailing whitespace
    if ( $new_product_name =~ /^(.+?)\s+$/ ) {
        $new_product_name = $1;
    }
    
    ## remove trailing periods, comma, hyphen, underscore, colon or forward slash
    if ( $new_product_name =~ /^(.+?)(\.|\,|\-|\_|\:|\/)$/ ) {
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
    
    
    ## if we have proteins...protein then remove 'proteins'
    #	ex. bacterial regulatory proteins, luxR family protein -- bacterial regulatory, luxR family protein
    if ($new_product_name =~ /proteins(.*protein)/ig) {
	$new_product_name =~ s/proteins(.*protein)/$1/ig;
	#some instances of 'proteins' are followed with a comma, so we need to remove the space before the comma
	$new_product_name =~ s/\s+(,\s+)/$1/;	
    }
    
    ## make sure we don't have some words more than once.
    ## something like abc domain domain protein
    while( $new_product_name =~ /(protein|domain|family|possible).*\1/ ) {
    	last if ($new_product_name =~/protein-\w+.*protein/);
        $new_product_name =~ s/$1\s*//;
    }

    ## multiple whitespace should be collapsed
    $new_product_name =~ s|\s{2,}| |g;

    ## remove periods that that come before 'domain protein'
    if( $new_product_name =~ /\.\s*(domain protein|family protein)?$/ ) {
        my $dp = $1 if( $1 );
        $new_product_name =~ s/\s*\.\s*(domain protein|family protein)?$//;
        $new_product_name .= " $dp" if( $dp );
    }
    
    ## if protein name ends with "binding" or "like", add "protein" to the end
    ## and remove "protein from the beginning"
    if ($new_product_name =~ /(binding|like)\s*$/) {
	$new_product_name =~ s/^(protein)\s//;
	$new_product_name .= ' protein';
    }
    ## remove strings with trailing parens
    if( $new_product_name =~ /\(.*\)\s*(domain protein|family protein)?$/ ) {
        my $dp = $1 if( $1 );
        $new_product_name =~ s/\s*\(.*\)\s*(domain protein|family protein)?$//;
        $new_product_name .= " $dp" if( $dp );
    }

    ## if the term is long (greater than 128 chars), try to remove extra info
    ## from parens or braces.
    if( length( $new_product_name ) > 128 ) {
        $new_product_name =~ s/[\[\(].*[\]\)]\.?//;
    }

    return $new_product_name;
}

# Subroutine to clean gene symbol
sub clean_gene_symbol {
	my ($new_gene_symbol) = @_;
	
	## gene sysmbol should not start with orf. It cannot have 4 consecutive numbers
	if(($new_gene_symbol =~ /^orf/i) || ($new_gene_symbol =~ /\d{4}/)) {
		$new_gene_symbol = "";
	}
	return $new_gene_symbol;
}
1;
