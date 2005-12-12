#!/usr/local/bin/perl -w

=head1  NAME 

move_wublastp_nr_results.pl - load wu-blastp/NR matches into the legacy filesystem schema.

=head1 SYNOPSIS

USAGE: move_wublastp_nr_results.pl 
        --btab_list=/path/to/hmmpfam.btab.list
        --raw_list=/path/to/hmmpfam.raw.list
        --repository_root=/usr/local/annotation/AA1
      [ --log=/path/to/some.log ]

=head1 OPTIONS

B<--btab_list,-b> 
    btab list file from an wu-blastp workflow component run.

B<--raw_list,-w> 
    raw list file from an wu-blastp workflow component run.

B<--repository_root,-r> 
    The project directory, just under which we should find the asmbls directory.

B<--log,-l> 
    optional.  will create a log file with summaries of all actions performed.

B<--help,-h> 
    This help message/documentation.

=head1   DESCRIPTION

Both the btab and raw alignment files from the wu-blastp run need to be migrated and 
renamed according to the following convention:

    $repository_root/asmbls/$asmblid/blastp/$asmblid.m$model.nr.btab
    $repository_root/asmbls/$asmblid/blastp/$asmblid.m$model.nr.gz
    
for example:

    /usr/local/annotation/CNA1/asmbls/179/blastp/179.m00468.nr.btab
    /usr/local/annotation/CNA1/asmbls/179/blastp/179.m00468.nr.gz
    /usr/local/annotation/CNA1/asmbls/179/blastp/179.m00624.nr.btab
    /usr/local/annotation/CNA1/asmbls/179/blastp/179.m00624.nr.gz
    and so on ...

The contents of the file is also  changed, since column 0 of the btab contains a 
reference to the model name, which is different between workflow searches and 
legacy convention.  Like the file names, this changes the column from a value 
like:

    bma1.model.15450_00008
    
to 

    15450.m00008

A simple find/replace is performed on the raw alignment file to make this change.

=head1 INPUT

The input is a list of btab files (--btab_list) and a list of raw alignment files 
(--raw_list) from a wu-blastp workflow component run.  This script will handle the
files properly both if the files withing the raw input list are compressed (using
gzip) or uncompressed.

=head1 OUTPUT

This is a file migration script.  There is no other output unless you use the --log option.

=head1 CONTACT

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
			  'btab_list|b=s',
              'raw_list|w=s',
              'repository_root|r=s',
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

#############################
## first do the btab files ##
#############################

## open the btab list
open (my $listfh, $options{btab_list}) || die "can't read btab list file\n";

while (my $infile = <$listfh>) {
    chomp $infile;
    next if ($infile =~ /^\s*$/);

    my ($asmblid, $modelnum, $model);

    ## files must end in .hmmpfam.htab
    if ( $infile =~ /(\d+)_(\d+)\.wu-blastp.btab/ ) {
        $asmblid = $1;
        $modelnum = $2;
        $model = "$asmblid.m$modelnum";
        
        _log("processing file $infile");
        
    } else {
        _log("ERROR: file $_ does not match naming convention! quitting.");
        exit(1);
    }
    
    ## make sure there is a directory for this asmbl id
    if (! -d "$options{repository_root}/asmbls/$asmblid" ) {
        ## we want this to die, since this should have been made already and
        ## may indicate a problem.
        _log("ERROR: no asmbl directory found for $asmblid, expected $options{repository_root}/asmbls/$asmblid");
        exit(1);
    }
    
    ## make sure there is a blastp directory
    _check_and_create_dir("$options{repository_root}/asmbls/$asmblid/blastp");
    
    ## open the input and output files
    open (my $ifh, "<$infile") || die "can't read input file $infile: $!";
    _log("creating $options{repository_root}/asmbls/$asmblid/blastp/$model.nr.btab");
    open (my $ofh, ">$options{repository_root}/asmbls/$asmblid/blastp/$model.nr.btab") || die "can't create output file: $!";
    
    ## read through the input file, changing the name in column 0 and writing to the output file
    while (my $line = <$ifh>) {
        chomp $line;
        next if ($line =~ /^\s*$/);
        
        my @cols = split("\t", $line);
        if ( $cols[0] =~ /\d+_\d+/) {
            $cols[0] = $model;
            print $ofh join("\t", @cols), "\n";
            
        } else {
            _log("WARNING: unrecognized line: $line");
        }
    }
}

close $listfh;

####################################
## now do the raw alignment files ##
####################################

## open the raw list
open ($listfh, $options{raw_list}) || die "can't read raw list file\n";

while (my $infile = <$listfh>) {
    chomp $infile;
    next if ($infile =~ /^\s*$/);

    my ($asmblid, $modelnum, $model);

    ## files must end in .hmmpfam.htab
    if ( $infile =~ /(\d+)_(\d+)\.wu-blastp.raw/ ) {
        $asmblid = $1;
        $modelnum = $2;
        $model = "$asmblid.m$modelnum";
        
        _log("processing file $infile");
        
    } else {
        _log("ERROR: file $_ does not match naming convention! quitting.");
        exit(1);
    }
    
    ## make sure there is a directory for this asmbl id
    if (! -d "$options{repository_root}/asmbls/$asmblid" ) {
        ## we want this to die, since this should have been made already and
        ## may indicate a problem.
        _log("ERROR: no asmbl directory found for $asmblid, expected $options{repository_root}/asmbls/$asmblid");
        exit(1);
    }
    
    ## make sure there is a blastp directory
    _check_and_create_dir("$options{repository_root}/asmbls/$asmblid/blastp");
    
    ## open the input and output files.  how we open the input file depends on whether
    ## it is compressed.
    my $ifh;
    if ($infile =~ /\.(gz|gzip)$/) {
        open ($ifh, "<:gzip", $infile) || die "can't read input file $infile: $!";
    } else {
        open ($ifh, "<$infile") || die "can't read input file $infile: $!";
    }
    
    
    _log("creating $options{repository_root}/asmbls/$asmblid/blastp/$model.nr.gz");
    open (my $ofh, ">:gzip", "$options{repository_root}/asmbls/$asmblid/blastp/$model.nr.gz") || die "can't create output file: $!";
    
    ## read through the input file, changing any name instances
    while (my $line = <$ifh>) {
        $line =~ s/\S+\.\S+\.${asmblid}_$modelnum/$model/g;
        print $ofh $line;
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
    unless ( defined $options{btab_list} && 
             defined $options{raw_list} &&
             defined $options{repository_root} ) {
        print STDERR "btab_list, raw_list and repository_root options are required\n\n";
        pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
    }
    
    ## make sure input lists exist
    unless ( -e $options{raw_list} ) {
        print STDERR "\n\nraw_list $options{raw_list} does not exist\n\n";
        exit(1);
    }
    
    unless ( -e $options{btab_list} ) {
        print STDERR "\n\nbtab_list $options{btab_list} does not exist\n\n";
        exit(1);
    }
    
    ## make sure the repository root has an asmbls directory
    unless ( -d "$options{repository_root}/asmbls" ) {
        print STDERR "\n\nthe repository root passed doesn't contain an asmbls directory.\n\n";
        exit(1);
    }
}






