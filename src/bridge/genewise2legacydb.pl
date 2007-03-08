#!/usr/local/bin/perl

=head1 NAME

genewise2legacydb.pl - uses Brian's geneWiseSearch.dbi to load genewise data
into the legacy database.

=head1 SYNOPSIS

    USAGE: genewise2legacydb.pl 
                --database=sma1
              [  --pipeline_config=/path/to/some/component/pipeline.config ]
              [ --search_db=search_db_name ]
              [ --input_file_list=/path/to/input_file_list ]
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

    Note that is --search_db and --input_file_list are used, this option is not required

B<--search_db,-s>
    The name (not full path) of the search_db, as specified in the path produced
    by move_genewise_results.pl

B<--input_file_list,-i>
    The full path (not just name) of the input file list.  This option must be used for
    version 2 components configured to run off bsml and not off the database.
    For most runs, where genewise has configured its own inputs, this input file is:
    $;OUTPUT_REPOSITORY$;/genewise/$;PIPELINE_ID$;_$;OUTPUT_TOKEN$;/genewise.input_file_listing

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
                          'search_db|s=s',
                          'input_file_list|i=s',
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

my $cfg;
if ($options{pipeline_config}) {
    $cfg = new Config::IniFiles( -file => $options{pipeline_config} );
}

## get the variables from the conf file
my $search_db       = $options{search_db} || $cfg->val( 'parameters', '$;SEARCH_DB$;' ) || $cfg->val( 'parameters genewise_best_loc', '$;SEARCH_DB$;' )  || die "couldn't find SEARCH_DB in config file, MUST be specified using --search_db otherwise!";

## holds the assembly ids to process
my @asmbl_ids;

## did the user specify a list of assembly ids, or just a single one?
if ($options{input_file_list}) {

    my %asmbl_buff;
    open (my $inputlist, $options{input_file_list}) || die "couldn't read input INPUT_FILE_LIST ($options{input_file_list}): $!";
    while (<$inputlist>) {
        if (/\.(\d+)\.\d+\.fsa/) {
            $asmbl_buff{$1}++;
        }
    }
    close $inputlist;
    push (@asmbl_ids, keys %asmbl_buff);

} elsif ( $cfg->val( 'input', '$;ASMBL_LIST_FILE$;' ) || $cfg->val( 'input genewise_best_loc', '$;INPUT_FILE_LIST$;' ) ) {
    my $asmbl_list_file = $cfg->val( 'input', '$;ASMBL_LIST_FILE$;' ) || $cfg->val( 'input genewise_best_loc', '$;ASMBL_LIST_FILE$;' )  || die "couldn't find ASMBL_LIST_FILE in config file";
    
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
    die "ASMBL_LIST_FILE or ASMBL_ID not defined in pipeline.config, MUST be specified via --input_file_list otherwise!";
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
    
    ## make sure database passed
    if (! defined $options{database}) {
        print "database not passed!\n";
        pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
    }

    ## make sure pipeline_config was passed and exists unless we're using the other inputs
    if (! defined $options{pipeline_config} ) {
        if ((! defined $options{search_db} )||( ! defined $options{input_file_list})) {
            print "pipeline_config not passed, and either search_db or input_file_list not given!\n";
            pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
        }
    } elsif (! -e $options{pipeline_config}) {
        print "the pipeline_config passed ($options{pipeline_config}) cannot be read or does not exist\n";
        pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
    }

    ## make sure the input_file_list exists if we are relying on it
    if (defined $options{input_file_list}) {
        unless (-e $options{input_file_list}) {
            print "the input_file_list given does not exist!\n";
            pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
        }
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

