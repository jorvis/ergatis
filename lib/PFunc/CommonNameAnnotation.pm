package PFunc::CommonNameAnnotation;

use strict;
use warnings;
use Exporter 'import';
use vars qw(@EXPORT_OK);

@EXPORT_OK = qw(clean_common_name);

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
        $new_product_name =~ /^\w{1,2}\d{1,3}$/ ||
        $new_product_name =~ /^ttg start/i ||
	$new_product_name =~ /gene \d+ protein/ ) {
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
    if ( $new_product_name =~ m|conserved hypothetical|i ||
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
    if ( $new_product_name =~ m|Y\S*[A-Z] family| || $new_product_name =~ /^\s*(Z|B)\d+\s+gene\s+product/ ) {
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
    $new_product_name =~ s/dimerisation/dimerization/g;
    $new_product_name =~ s/disulphide/disulfide/g;
    $new_product_name =~ s/sulphur/sulfur/g;
    if( $new_product_name =~ /\bhaem(\w+)/ ) {
        my $p = $1;
        $new_product_name =~ s/haem$p/hem$p/;
    }

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
    
    ## make sure we don't have some words more than once.
    ## something like abc domain domain protein
    while( $new_product_name =~ /(protein|domain|family|possible).*\1/ ) {
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

1;
