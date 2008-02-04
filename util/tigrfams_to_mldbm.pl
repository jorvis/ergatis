#!/usr/bin/perl

=head1 NAME

tigrfams_to_mldbm.pl - Reads an HMM file and creates an MLDBM database of commonly-accessed 
attributes for each HMM.

=head1 SYNOPSIS

USAGE: tigrfams_to_mldbm.pl 
            --hmm_file=/path/to/some_file.lib
            --output_file=/path/to/some.db

=head1 OPTIONS

B<--hmm_file>
    Input file can be a single HMM or collection.  See INPUT section below for format example.

B<--output_file>
    Ouput MLDBM file to be created from the parsed input file.  See OUTPUT section for data
    structure definition.

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

Several scripts perform lookups on HMM attributes such as annotation, trusted cutoff, etc and need
a rapid way to access this information.  This script reads an HMM library and generates an MLDBM
file of the most commonly accessed attributes.

=head1  INPUT

As described above, the input can be a library of HMMs or a single HMM file.  Records are separated
with a // symbol.  An example record:

    HMMER2.0
    NAME  S16
    ACC   TIGR00002
    DESC  S16: ribosomal protein S16
    LENG  81
    ALPH  Amino
    RF    no
    CS    no
    MAP   yes
    COM   hmmbuild --pam blosum62 --pamwgt 20 ogc_results//853/17288/17288_1_3.HMM ogc_results//853/17288/17288_1.msf
    COM   hmmcalibrate --cpu 1 --num 5000 --histfile ogc_results//853/17288/17288_1_3.hist ogc_results//853/17288/17288_1_3.HMM
    NSEQ  20
    DATE  Mon Feb 11 17:46:41 2002
    CKSUM 95
    GA    35.00 35.00
    TC    35.00 35.00
    NC    10.00 10.00
    XT      -8455     -4  -1000  -1000  -8455     -4  -8455     -4 
    NULT      -4  -8455
    NULE     595  -1558     85    338   -294    453  -1158    197    249    902  -1085   -142    -21   -313     45    531    201    384  -1998   -644 
    EVD   -40.960884   0.269394
    HMM        A      C      D      E      F      G      H      I      K      L      M      N      P      Q      R      S      T      V      W      Y    
             m->m   m->i   m->d   i->m   i->i   d->m   d->d   b->m   m->e
              -71      *  -4389
         1  -1123  -1482  -2599  -2123  -1302  -2623  -2482    844  -2004   1300   -331  -2482  -2123  -2004  -2331  -1982  -1123   2752  -2331  -1482     1
         -   -149   -500    233     43   -381    399    106   -626    210   -466   -720    275    394     45     96    359    117   -369   -294   -249 
         -     -3  -9476 -10518   -894  -1115   -701  -1378    -71      *

        ... [ and so on for all positions ] ...

             81   -439  -1642  -2419  -2205   2226  -2269  -1903   1024  -2082    446   1484  -2214  -2388  -1909  -2082  -1628    553   -751   3736   -339    91
         -      *      *      *      *      *      *      *      *      *      *      *      *      *      *      *      *      *      *      *      * 
         -      *      *      *      *      *      *      *      *      0 
//


=head1  OUTPUT

The output file is created using the MLDBM library, which can be accessed as a tied hash with
this structure:

    $h->{$accession} = {
                            hmm_com_name => ?,
                            hmm_len => ?,
                            trusted_cutoff => ?,
                            noise_cutoff => ?,
                            trusted_cutoff2 => ?,
                            noise_cutoff2 => ?,
                       };

Future attributes that could be added:

iso_type, ec_num, gene_sym, role_ids, go_terms

=head1  CONTACT

    Joshua Orvis
    jorvis@users.sf.net

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use MLDBM 'DB_File';
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
                          'hmm_file=s',
                          'output_file=s',
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


open(my $ifh, "<$options{hmm_file}") || die "can't read input file: $!";

my %info;

## create the tied hash
tie(%info, 'MLDBM', $options{'output_file'});

## store attributes for each HMM as we encounter them.  This is purged after each
##  HMM in the file
my $att;

while ( my $line = <$ifh> ) {
    
    ## are we at the end of an entry?
    if ( $line =~ m|^//| && defined $att ) {
        
        my ($trusted_global, $trusted_domain, $noise_global, $noise_domain);
        
        if ( $$att{TC} =~ /([0-9\.\-]+)\s+([0-9\.\-]+)/ ) {
            ($trusted_global, $trusted_domain) = ($1, $2);
            
        } else {
            die "failed to parse trusted global and trusted domains for $$att{ACC}\n";
        }
        
        if ( $$att{NC} =~ /([0-9\.\-]+)\s+([0-9\.\-]+)/ ) {
            ($noise_global, $noise_domain) = ($1, $2);
            
        } else {
            die "failed to parse noise global and noise domains for $$att{ACC}\n";
        }
        
        ## handle the product name.  if it contains the NAME, remove that
        my $product_name = $$att{DESC};
        
        if ( $product_name =~ /^$$att{NAME}: (.+)/ ) {
            $product_name = $1;
        }
        
        ## add entry to MLDBM
        $info{$$att{ACC}} = {
                                hmm_com_name => $product_name,
                                hmm_len => $$att{LENG},
                                trusted_cutoff => $trusted_global,
                                noise_cutoff => $noise_global,
                                trusted_cutoff2 => $trusted_domain,
                                noise_cutoff2 => $noise_domain,
                            };
        
        undef $att;
        next;
    }
    
    if ( $line =~ /^([A-Z]+)\s+(.*)$/ ) {
        $att->{$1} = $2;
    }
}

untie(%info);

exit(0);


sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
}


sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( hmm_file output_file );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
    
    ##
    ## you can do other things here, such as checking that files exist, etc.
    ##
    
    ## handle some defaults
    #$options{optional_argument2}   = 'foo'  unless ($options{optional_argument2});
}
