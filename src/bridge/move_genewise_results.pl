#!/usr/local/bin/perl

=head1 NAME

move_genewise_files.pl - copies the files created by the workflow_best_loc component into
the standard legacy output directory structure.

=head1 SYNOPSIS

    USAGE: move_genewise_files.pl 
                --repository_root=/usr/local/annotation/AA1
                --pipeline_config=/path/to/some/component/pipeline.config
              [ --log=/path/to/some/optional/file.log ]

=head1 OPTIONS

B<--repository_root,-r>
    The repository root of some legacy project.  The 'asmbls' directory should be inside this one.

B<--pipeline_config,-p>
    The path to some component's pipeline.config file.  This will be created when running a 
    genewise component in Workflow and should have the following variables within:
    
        $;SEARCH_DBS$;
        $;INPUT_FILE_LIST$;
        $;RAW_OUTPUT_LIST$;

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

This program is used to move the output from the genewise_best_loc component into the
legacy database structure devised by Brian.  Workflow output will contain the .fsa and
.pep files that were used in the analysis, along with the raw genewise alignment file.
These are moved into a directory structure like:

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
                          'repository_root|r=s',
                          'pipeline_config|p=s',
                          'log|l=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

## open the log file if needed:
my $logfh;

if ($options{log}) {
    open($logfh, ">$options{log}") || die "can't open log file: $!";
}

my $cfg = new Config::IniFiles( -file => $options{pipeline_config} );

## get the variables from the conf file
my $search_db       = $cfg->val( 'parameters genewise_best_loc', '$;SEARCH_DB$;' )   || $cfg->val( 'parameters', '$;SEARCH_DB$;' ) || die "couldn't find SEARCH_DB in config file";
my $input_file_list = $cfg->val( 'input genewise_best_loc', '$;INPUT_FILE_LIST$;' )  || $cfg->val( 'input', '$;INPUT_FILE_LIST$;' ) || die "couldn't find INPUT_FILE_LIST in config file";
my $raw_output_list = $cfg->val( 'output genewise_best_loc', '$;RAW_OUTPUT_LIST$;' ) || $cfg->val( 'output', '$;RAW_OUTPUT_LIST$;' ) || die "couldn't find RAW_OUTPUT_LIST in config file";

## for each of the files defined in the input_file_list we want to copy the .fsa and .pep files
##  to the standard legacy output directory structure.
##  the name will be converted from a convention like:
##      aa1.assembly.24834.1.fsa
##      aa1.assembly.24834.1.pep
##  to 
##      24834.1.fsa
##      24834.1.pep
print $logfh "reading input file list: $input_file_list\n" if $logfh;
open(my $ifl_fh, "<$input_file_list") || die "can't read input file list: $!";

my ($root, $asmbl_id, $name, $dir);
while (my $line = <$ifl_fh>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    
    if ( $line =~ m|^(.+/.+\.((\d+)\.\d+))\.fsa$| ) {
        ($root, $name, $asmbl_id) = ($1, $2, $3);
        print $logfh "extracted $asmbl_id and $name from $line\n" if $logfh;
        
        $dir = "$options{repository_root}/asmbls/$asmbl_id/genewise/$search_db";
        
        ## make the directory if it doesn't exist
        if (! -e $dir) {
            print $logfh "executing: mkdir -p $dir\n" if $logfh;
            system("mkdir -p $dir");
        }
        
        print $logfh "executing cp $root.fsa $dir/$name.fsa\n" if $logfh;
        copy("$root.fsa", "$dir/$name.fsa") || die "FAILED: cp $root.fsa $dir/$name.fsa BECAUSE $!";
        print $logfh "executing cp $root.pep $dir/$name.pep\n" if $logfh;
        copy("$root.pep", "$dir/$name.pep") || die "FAILED: cp $root.pep $dir/$name.pep BECAUSE $!";
        
    } else {
        die "unrecognized file path format:\n$line";
    }
}

close $ifl_fh;



## for each of the files defined in the raw_output_list we want to copy it to the standard legacy
##  output directory structure.
##  the name will be converted from a convention like:
##      aa1.assembly.25179.9.genewise_best_loc.raw
##  to 
##      25179.9.genewise
print $logfh "reading raw output list: $raw_output_list\n" if $logfh;
open(my $rol_fh, "<$raw_output_list") || die "can't read raw output list: $!";

while (my $line = <$rol_fh>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    
    if ( $line =~ m|^(.+/.+\.((\d+)\.\d+)\.genewise_best_loc)\.raw$| || $line =~ m|^(.+/.+\.((\d+)\.\d+)\.genewise)\.raw$|) {
        ($root, $name, $asmbl_id) = ($1, $2, $3);
        print $logfh "extracted $asmbl_id and $name from $line\n" if $logfh;
        
        $dir = "$options{repository_root}/asmbls/$asmbl_id/genewise/$search_db";
        
        ## the directory must exist (due to the copies above), so we don't need to check
        ##  for it here.
        
        print $logfh "executing cp $root.raw $dir/$name.genewise\n" if $logfh;
        copy("$root.raw", "$dir/$name.genewise");
        
    } else {
        die "unrecognized file path format:\n$line";
    }    

}

close $rol_fh;



close $logfh if $logfh;

exit;

sub check_parameters {
    my $options = shift;
    
    ## make sure repository_root was passed and exists exists
    if (! defined $options{repository_root}) {
        print "repository_root not passed!\n";
        pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
    
    } elsif (! -e $options{repository_root}) {
        print "the repository_root passed ($options{repository_root}) cannot be read or does not exist\n";
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
    
    ## set any defaults
    $options{log} = 0 if (! defined $options{log});
}

