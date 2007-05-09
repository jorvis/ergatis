#!/usr/local/bin/perl -w

=head1  NAME 

move_hmmpfam_results.pl - moving hmmpfam data into the legacy filesystem schema.

=head1 SYNOPSIS

USAGE: move_hmmpfam_results.pl 
        --input_list=/path/to/hmmpfam.htab.list
        --repository_root=/usr/local/annotation/AA1
      [ --tu_list=/path/to/output.tu.list
        --log=/path/to/some.log
      ]

=head1 OPTIONS

B<--input_list,-i> 
    htab list file from an hmmpfam workflow component run.

B<--repository_root,-r> 
    The project directory, just under which we should find the asmbls directory.

B<--tu_list,-t> 
    This creates a list of the model names processed (legacy convention), which
    can be used by the hmm2_search.condor.dbi script to load the database after
    the files have migrated.

B<--log,-l> 
    optional.  will create a log file with summaries of all actions performed.

B<--help,-h> 
    This help message/documentation.

=head1   DESCRIPTION

The hmmpfam files need to be migrated and renamed according to the following convention:

    $repository_root/asmbls/$asmblid/HMM_searches/CURRENT/$asmblid.m$model.htab
    
for example:

    /usr/local/annotation/BMA1/asmbls/14040/HMM_searches/CURRENT/14040.m00262.htab
    /usr/local/annotation/BMA1/asmbls/14040/HMM_searches/CURRENT/14040.m00263.htab
    /usr/local/annotation/BMA1/asmbls/14040/HMM_searches/CURRENT/14040.m00264.htab
    and so on ...

The contents of the file is also  changed, since column 6 contains a reference to
the model name, which is different between workflow searches and legacy convention.  
Like the file names, this changes the column from a value like:

    bma1.model.15450_00008
    
to 

    15450.m00008


=head1 INPUT

The input is a list of htab files from an hmmpfam workflow component run.

=head1 OUTPUT

This is a file migration script.  There is no other output unless you use the --log option.

=head1 LOADING DATA INTO THE LEGACYDB

This script does not provide this functionality directly, but is required as an 
initial step.  After you've migrated the results, you can load them like this:

    $EGC_SCRIPTS/hmm2_search.condor.dbi -D bma1 -L model_names.list -p $EGC_SCRIPTS/egc_password -x

The -x is important there, as it tells the script just to load data and not actually 
launch any search commands. (Brian's script)

=head1 CONTACT

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
			  'input_list|i=s',
              'repository_root|r=s',
              'tu_list|t=s',
              'log|l=s',
			  'help|h') || pod2usage();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

## play nicely
umask(0000);

## make sure all passed options are peachy
&check_parameters(\%options);

## get the current user
my ($user) = getpwuid($<);

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">$options{log}") || die "can't create log file: $!";
}

## create a TU list if requested
my $tufh;
if (defined $options{tu_list}) {
    open($tufh, ">$options{tu_list}") || die "can't create tu list: $!";
}

## open the input list
open (my $listfh, $options{input_list}) || die "can't read list file\n";

while (my $infile = <$listfh>) {
    chomp $infile;
    next if ($infile =~ /^\s*$/);

    my ($asmblid, $modelnum, $model);

    ## files must end in .hmmpfam.htab
    ## OR be of the form model.(\d+)_(\d+).htab
    if ( ($infile =~ /(\d+)_(\d+)\.hmmpfam.htab/) ||
         ($infile =~ /model\.(\d+)_(\d+)\.htab/)  ) {
        $asmblid = $1;
        $modelnum = $2;
        $model = "$asmblid.m$modelnum";
        
        ## are we creating a TU list?
        if ( $tufh ) {
            print $tufh "$model\n";
        }
        
        _log("processing file $infile");
        
    } else {
        _log("ERROR: file $_ does not match naming convention! quitting.");
        exit(1);
    }
    
    ## make sure there is a directory for this asmbl id
    _check_and_create_dir("$options{repository_root}/asmbls/$asmblid");
    
    ## make sure there is an HMM_searches directory
    _check_and_create_dir("$options{repository_root}/asmbls/$asmblid/HMM_searches");
    
    ## make sure there is a CURRENT directory
    _check_and_create_dir("$options{repository_root}/asmbls/$asmblid/HMM_searches/CURRENT");
    
    ## open the input and output files
    open (my $ifh, "<$infile") || die "can't read input file $infile: $!";
    _log("creating $options{repository_root}/asmbls/$asmblid/HMM_searches/CURRENT/$model.htab");
    open (my $ofh, ">$options{repository_root}/asmbls/$asmblid/HMM_searches/CURRENT/$model.htab") || die "can't create output file: $!";
    
    ## read through the input file, changing the name in column 6 and writing to the output file
    while (my $line = <$ifh>) {
        chomp $line;
        next if ($line =~ /^\s*$/);
        
        my @cols = split("\t", $line);
        if ( $cols[5] ) {
            $cols[5] = $model;
            print $ofh join("\t", @cols), "\n";
                      
        } elsif ($line =~ /No hits above thresholds/i) {
            print $ofh "$line\n";
            
        } else {
            _log("WARNING: unrecognized line: $line");
        }
    }
}

exit(0);

sub _check_and_create_dir {
    my $dir = shift;
    
    if (! -d $dir ) {
        _log("creating $dir");
        mkdir("$dir");
    }    
}

sub _log {
    my $msg = shift;
    
    print $logfh "$msg\n" if $logfh;
}

sub check_parameters {
    
    ## input_list and repository_root are required
    unless ( defined $options{input_list} && $options{repository_root} ) {
        print STDERR "input_list and repository_root options are required\n\n";
        pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
    }
    
    ## make sure input list exists
    unless ( -e $options{input_list} ) {
        print STDERR "\n\ninput_list $options{input_list} does not exist\n\n";
        exit(1);
    }
    
    ## make sure the repository root has an asmbls directory
    unless ( -d "$options{repository_root}/asmbls" ) {
        print STDERR "\n\nthe repository root passed doesn't contain an asmbls directory.\n\n";
        exit(1);
    }
}






