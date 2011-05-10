#!/usr/bin/perl

=head1 NAME

tigrfam_info_to_mldbm.pl - Reads a directory of TIGR*.INFO files and adds attributes to an 
existing HmmInfo MLDBM file.

=head1 SYNOPSIS

USAGE: tigrfam_info_to_mldbm.pl 
            --mldbm_file=/path/to/some_file.db
            --info_dir=/path/to/some_dir
          [ --go_link=/path/to/TIGRFAMS_GO_LINK 
	    --role_link=/path/to/TIGRFAMS_ROLE_LINK
	  ]

=head1 OPTIONS

B<--mldbm_file>
    This should be an MLDBM file created by the hmmlib_to_mldbm script.

B<--info_dir>
    A directory of files named like TIGR00572.INFO with one 

B<--go_link>
    A tab-delimited file containing relationships of TIGRFAMs to GO terms
    
B<--role_link>
    A tab-delimited file containing relationships of TIGRFAMs to TIGR role ids

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
    AC  TIGR00877
    TC  200.00 200.00
    NC  -50.00 -50.00
    CC  This enzyme appears as a monofunctional protein in prokaryotes but as part of a larger, multidomain protein in eukaryotes.  

The --go_link is a tab-delimited text file with this format:

    TIGR00025       GO:0005887      NULL
    TIGR00025       GO:0006810      NULL
    TIGR00025       GO:0043190      NULL
    TIGR00025       GO:0042626      contributes_to
    TIGR00029       GO:0003735      NULL
    TIGR00029       GO:0006412      NULL
    TIGR00029       GO:0000312      NULL

Currently, this script only considers the first two columns.  This file is optional.

The --role_link is a tab-delimited test file with this format:

    TIGR00001	158
    TIGR00002	158
    TIGR00003	145
    TIGR00004	134
    TIGR00005	168
    TIGR00006	89
    TIGR00007	161
    TIGR00008	169
    TIGR00009	158


=head1  OUTPUT

The following are added to the existing tied hash.  Each accession entry should exist alraedy.

    $h->{$accession} = {
                            isotype => 'equivalog_domain'
                            ec_num => '6.3.4.13',
                            gene_sym => 'purD',
                            go => [ 'GO:0003735', 'GO:0006412', ... ]
			    roleid => 158
                       };

=head1  CONTACT

    Joshua Orvis
    jorvis@users.sf.net

=cut

use strict;
use Fcntl qw( O_RDWR );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use MLDBM 'DB_File';
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
                          'mldbm_file=s',
                          'info_dir=s',
                          'go_link=s',
			  'role_link=s', 
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

## read in all the GO associations first, if there are any
my $go_terms = {};
if ( $options{go_link} ) {
    _log("INFO: reading go_link file $options{go_link}");
    open(my $go_fh, "<$options{go_link}") || die "couldn't read --go_link file: $!";
    
    while ( my $line = <$go_fh> ) {
        chomp $line;
        next if $line =~ /^\s*$/;
        
        my @cols = split(/\t/, $line);
        push @{$$go_terms{$cols[0]}}, $cols[1];
    }
}

## read all role ids for TIGRFAMs
my $role_ids = {};
if ($options{role_link}) {
    _log("INFO: reading role_link file $options{role_link}");
    open(my $role_fh, "<$options{role_link}") || die "couldn't read --role_link file: $!";
    while ( my $line = <$role_fh> ) {
         chomp $line;
         next if $line =~ /^\s*$/;
         my @cols = split(/\t/, $line);
	 push @{$$role_ids{$cols[0]}}, $cols[1];
    }
}

opendir(my $idh, "$options{info_dir}") || die "can't read info dir: $!";

my %info;

## open the tied hash
tie(%info, 'MLDBM', $options{'mldbm_file'}, O_RDWR );

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
    my $comment = $info{$accession}{hmm_comment} || '';
    my $trusted_cutoff = $info{$accession}{trusted_cutoff} || '';
    my $trusted_cutoff2 = $info{$accession}{trusted_cutoff2} || '';
    my $noise_cutoff = $info{$accession}{noise_cutoff} || '';
    my $noise_cutoff2 = $info{$accession}{noise_cutoff2} || '';
    my $gathering_cutoff = $info{$accession}{gathering_cutoff} || '';
    my $gathering_cutoff2 = $info{$accession}{gathering_cutoff2} || '';

    ## we have to store this into a variable first, make modifications, then write it
    ##  back.  this is dumb, but see the BUGS section of the MLDBM perldoc for the
    ##  explanation.
    my $tied_acc = $info{$accession};

    while ( my $line = <$ifh> ) {
        chomp $line;
        
        ## look for isotype lines
        if ( $line =~ /^IT\s+(.+)$/ ) {
            _log("INFO: $accession isotype was [$isotype], changing to [$1]");
            $$tied_acc{isotype} = $1;
        
        ## look for ec num lines
        } elsif ( $line =~ /^EC\s+(.+)$/ ) {
            _log("INFO: $accession ec_num was [$ec_num], changing to [$1]");
            $$tied_acc{ec_num} = $1;
        } elsif ( $line =~ /^GS\s+(.+)$/ ) {
            _log("INFO: $accession gene_symbol was [$gene_symbol], changing to [$1]");
            $$tied_acc{gene_symbol} = $1;
        } elsif ( $line =~ /^CC\s+(.+)$/ ) {
            _log("INFO: $accession hmm_comment was [$comment], changing to [$1]");
            $$tied_acc{hmm_comment} = $1;
        } elsif ( $line =~ /^TC\s+([0-9\.\-]*)\s+([0-9\.\-]*)$/ ) {
	    _log("INFO: $accession trusted global cutoff was [$trusted_cutoff], changing to [$1] and trusted domain cutoff was [$trusted_cutoff2], changing to [$2]");
	    $$tied_acc{trusted_cutoff} = $1;
            $$tied_acc{trusted_cutoff2} = $2;
	} elsif ( $line =~ /^NC\s+([0-9\.\-]*)\s+([0-9\.\-]*)$/ ) { 
            _log("INFO: $accession noise global cutoff was [$noise_cutoff], changing to [$1] and noise domain cutoff was [$noise_cutoff2], changing to [$2]");
            $$tied_acc{noise_cutoff} = $1; 
            $$tied_acc{noise_cutoff2} = $2; 
        } elsif ( $line =~ /^GA\s+([0-9\.\-]*)\s+([0-9\.\-]*)$/ ) { 
            _log("INFO: $accession gathering global cutoff was [$gathering_cutoff], changing to [$1] and gathering domain cutoff was [$gathering_cutoff2], changing to [$2]");
            $$tied_acc{gathering_cutoff} = $1; 
            $$tied_acc{gathering_cutoff2} = $2; 
        }
    }
    
    ## now add any GO terms
    if ( exists $$go_terms{$accession} ) {
        _log("INFO: attaching " . scalar( @{$$go_terms{$accession}} ) . " to accession $accession");
        $$tied_acc{go} = $$go_terms{$accession};
    } else {
        _log("INFO: no GO terms found for accession $accession");
    }
        
    ## now add any role_ids
    if ( exists $$role_ids{$accession} ) { 
        _log("INFO: attaching " . scalar(@{$$role_ids{$accession}}) . " to accession $accession");
        $$tied_acc{role_id} = $$role_ids{$accession};
    } else {
        _log("INFO: no role_id found for accession $accession");
    } 

    $info{$accession} = $tied_acc;
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
    $options{go_link} = '' unless ($options{go_link});
    $options{role_link} = '' unless ($options{role_link});
}
