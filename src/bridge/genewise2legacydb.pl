#!/usr/local/bin/perl

=head1 NAME

genewise2legacydb.pl - uses Brian's geneWiseSearch.dbi to load genewise data
into the legacy database.

=head1 SYNOPSIS

    USAGE: genewise2legacydb.pl 
                --database=sma1
                --pipeline_config=/path/to/some/component/pipeline.config
              [ --log=/path/to/some/optional/file.log ]

=head1 OPTIONS

B<--database,-d>
    The Sybase database to load.

B<--pipeline_config,-p>
    The path to some component's pipeline.config file.  This will be created when running a 
    genewise component in Workflow and should have the following variables within:
    
        $;SEARCH_DB$;
        $;ASMBL_LIST_FILE$; (first check)
        or
        $;ASMBL_ID$;

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

This program assumes you've already run move_genewise_results.pl, which places the output
files into a directory convention like:

    $repository_root/asmbls/$asmbl_id/genewise/$database/
    
and the files are named like:

    /$asmbl_id.$fragment.fsa
    /$asmbl_id.$fragment.pep
    /$asmbl_id.$fragment.genewise

=head1  CONTACT

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;
use Config::IniFiles;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Copy;
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
                          'database|d=s',
                          'pipeline_config|p=s',
                          'log|l=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

my $loader = '/usr/local/devel/ANNOTATION/euk_genome_control/bin/geneWiseSearch.dbi';

## open the log file if needed:
my $logfh;

if ($options{log}) {
    open($logfh, ">$options{log}") || die "can't open log file: $!";
}

my $cfg = new Config::IniFiles( -file => $options{pipeline_config} );

## get the variables from the conf file
my $search_db       = $cfg->val( 'parameters', '$;SEARCH_DB$;' ) || $cfg->val( 'parameters genewise_best_loc', '$;SEARCH_DB$;' )  || die "couldn't find SEARCH_DB in config file";

## holds the assembly ids to process
my @asmbl_ids;

## did the user specify a list of assembly ids, or just a single one?
if ( $cfg->val( 'input', '$;ASMBL_LIST_FILE$;' ) || $cfg->val( 'input genewise_best_loc', '$;INPUT_FILE_LIST$;' ) ) {
    my $asmbl_list_file = $cfg->val( 'input', '$;ASMBL_LIST_FILE$;' ) || $cfg->val( 'input genewise_best_loc', '$;ASMBL_FILE_LIST$;' )  || die "couldn't find ASMBL_FILE_LIST in config file";
    
    print $logfh "processing asmbl list file $asmbl_list_file\n" if $logfh;
    
    ## open the list file and get each asmbl id
    open (my $lfh, $asmbl_list_file) || die "couldn't read input ASMBL_LIST_FILE ($asmbl_list_file): $!";
    while (<$lfh>) {
        chomp;
        if ( /^\s*(\d+)\s*$/ ) {
            push @asmbl_ids, $1;
        }
    }
    
} elsif ( $cfg->val( 'input', '$;ASMBL_ID$;' ) || $cfg->val( 'input genewise_best_loc', '$;ASMBL_ID$;' )  ) {
    my $aId = $cfg->val( 'input', '$;ASMBL_ID$;' ) || $cfg->val( 'input genewise_best_loc', '$;ASMBL_ID$;'  ) || die "couldn't read input ASMBL_ID";
    push @asmbl_ids, $aId;
    
} else {
    die "ASMBL_LIST_FILE or ASMBL_ID not defined in pipeline.config";
}

## use brian's script to load the results of each assembly
my $cmd;
foreach my $asmbl_id ( @asmbl_ids ) {
    print $logfh "loading asmbl_id $asmbl_id\n" if $logfh;
    $cmd = "$loader -D $options{database} -p /home/jorvis/.pwdfile -n $search_db -a $asmbl_id -j";
    print $logfh "$cmd\n" if $logfh;
    `$cmd`;
}


close $logfh if $logfh;

exit;

sub check_parameters {
    my $options = shift;
    
    ## make sure repository_root was passed and exists exists
    if (! defined $options{database}) {
        print "database not passed!\n";
        pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
    }

    ## make sure pipeline_config was passed and exists exists
    if (! defined $options{pipeline_config}) {
        print "pipeline_config not passed!\n";
        pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
    
    } elsif (! -e $options{pipeline_config}) {
        print "the pipeline_config passed ($options{pipeline_config}) cannot be read or does not exist\n";
        pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
    }
    
    ## make sure the user has a password file in their home directory
    unless ( -e "/home/jorvis/.pwdfile" ) {
        print "\n\ncouldn't find password file /home/jorvis/.pwdfile .  quitting.\n\n";
        exit(1);
    }
    
    ## set any defaults
    $options{log} = 0 if (! defined $options{log});
    $options{database} = lc $options{database};
}

