#!/usr/bin/perl

=head1 NAME

tigrfam_info_to_mldbm.pl - Reads a directory of TIGR*.INFO files and adds attributes to an 
existing HmmInfo MLDBM file.

=head1 SYNOPSIS

USAGE: tigrfam_info_to_mldbm.pl 
            --mldbm_file=/path/to/some_file.db
            --info_dir=/path/to/some_dir

=head1 OPTIONS

B<--mldbm_file>
    This should be an MLDBM file created by the hmmlib_to_mldbm script.

B<--info_dir>
    A directory of files named like TIGR00572.INFO with one 

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

The script hmmlib_to_mldbm.pl reads an HMM library and creates a perl tied hash to disk
for the most commonly accessed attribute of each HMM.  This script adds additional annotations
to that data when they're available.

=head1  INPUT

The --info_dir is read for .INFO files, each with a format like this:

    ID  purD
    AC  TIGR00877
    DE  phosphoribosylamine--glycine ligase
    AU  Haft DH
    TC  200.00 200.00
    NC  -50.00 -50.00
    AL  clustalw_manual
    IT  equivalog_domain
    EN  phosphoribosylamine--glycine ligase
    GS  purD
    EC  6.3.4.13
    TP  TIGRFAMs
    CC  Alternate name: glycinamide ribonucleotide synthetase (GARS).
    CC  This enzyme appears as a monofunctional protein in prokaryotes but as part of a larger, multidomain protein in eukaryotes.
    DR  EXPERIMENTAL; EGAD|7937|EC4005; Escherichia coli
    DR  ECOCYC; EG10792; purD

The only fields I currently look for are these:

    IT  equivalog_domain
    GS  purD
    EC  6.3.4.13

=head1  OUTPUT

The following are added to the existing tied hash.  Each accession entry should exist alraedy.

    $h->{$accession} = {
                            isotype => 'equivalog_domain'
                            ec_num => '6.3.4.13',
                            gene_sym => 'purD',
                       };

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
                          'mldbm_file=s',
                          'info_dir=s',
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


opendir(my $idh, "$options{info_dir}") || die "can't read info dir: $!";

my %info;

## create the tied hash
tie(%info, 'MLDBM', $options{'mldbm_file'});

_log("INFO: scanning $options{info_dir} for .INFO files");
while ( my $file = readdir($idh) ) {
    ## make sure it's a file
    next unless ( -f "$options{info_dir}/$file");

    _log("INFO: processing file $options{info_dir}/$file");

    my $accession = '';

    ## it needs to match our pattern
    if ( $file =~ /^(TIGR\d+)\.INFO$/ ) {
        $accession = $1;
    } else {
        _log("WARN: skipping $file because it doesn't match TIGRNNNNN.INFO");
        next;
    }
    
    ## make sure this accession even exists in hmmInfo
    if (! exists $info{$accession} ) {
        _log("WARN: skipping accession $accession because it doesn't exist in the MLDBM");
        next;
    }
    
    open(my $ifh, "<$options{info_dir}/$file") || die "can't read file $options{info_dir}/$file: $!";
    
    my $ec_num = $info{$accession}{ec_num} || '';
    my $gene_symbol = $info{$accession}{gene_symbol} || '';
    my $isotype = $info{$accession}{isotype} || '';
    
    while ( my $line = <$ifh> ) {
        chomp $line;
        
        ## look for isotype lines
        if ( $line =~ /^IT\s+(.+)$/ ) {
            _log("INFO: $accession isotype was [$isotype], changing to [$1]");
            $info{$accession}{isotype} = $1;
        
        ## look for ec num lines
        } elsif ( $line =~ /^EC\s+(.+)$/ ) {
            _log("INFO: $accession ec_num was [$ec_num], changing to [$1]");
            $info{$accession}{ec_num} = $1;
            
        } elsif ( $line =~ /^GS\s+(.+)$/ ) {
            _log("INFO: $accession gene_symbol was [$gene_symbol], changing to [$1]");
            $info{$accession}{gene_symbol} = $1;
        }
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
    my @required = qw( mldbm_file info_dir );
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
