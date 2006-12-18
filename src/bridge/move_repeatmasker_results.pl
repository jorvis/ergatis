#!/usr/local/bin/perl -w

=head1  NAME 

move_repeatmasker_results.pl - load repeatmasker data into the legacy filesystem schema.

=head1 SYNOPSIS

USAGE: move_repeatmasker_results.pl 
        --input_list=/path/to/repeatmasker.fasta.list
        --repository_root=/usr/local/annotation/AA1
      [ --log=/path/to/some.log ]

=head1 OPTIONS

B<--input_list,-i> 
    fasta list file from a repeatmasker workflow component run.

B<--repository_root,-r> 
    The project directory, just under which we should find the asmbls directory.

B<--log,-l> 
    optional.  will create a log file with summaries of all actions performed.

B<--help,-h> 
    This help message/documentation.

=head1   DESCRIPTION

The repeatmasker files need to be migrated and renamed according to the following 
convention:

    $repository_root/asmbls/$asmblid/$asmblid.contig.repeatmasker
    
for example:

    /usr/local/annotation/BMA1/asmbls/14040/14040.contig.repeatmasker

The contents of the header line is changed to match the legacy convention like:

    >asmbl_id

=head1 INPUT

The input is a list of repeatmasked fasta files from a repeatmasker workflow 
component run.

=head1 OUTPUT

This is a file migration script.  There is no other output unless you use 
the --log option.

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

## open the input list
open (my $listfh, $options{input_list}) || die "can't read list file\n";

while (my $infile = <$listfh>) {
    chomp $infile;
    
	next if ($infile =~ /^\s*$/);

    my $asmblid;

    if ( $infile =~ /assembly\.(\d+)\./ ) {
        $asmblid = $1;
        _log("processing file $infile");
        
    } else {
        print STDERR "ERROR: file $_ does not match naming convention!  skipping.\n";
        _log("ERROR: file $_ does not match naming convention!  skipping.");
        next;
    }
    
    ## make sure there is a directory for this asmbl id
    if (! -d "$options{repository_root}/asmbls/$asmblid" ) {
        print STDERR "ERROR: no asmbl directory found for $asmblid, expected $options{repository_root}/asmbls/$asmblid\n";
        _log("ERROR: no asmbl directory found for $asmblid, expected $options{repository_root}/asmbls/$asmblid");
        next;
    }

    ## copy the file
    _log("cp $infile $options{repository_root}/asmbls/$asmblid/$asmblid.contig.repeatmasker");
    
    open(my $ifh, "<$infile") || die "can't read input file: $!";
    open(my $ofh, ">$options{repository_root}/asmbls/$asmblid/$asmblid.contig.repeatmasker") || die "can't write output file: $!";
    
    my $header_found = 0;
    
    while (<$ifh>) {
        if ( /\>.*assembly\.(\d+)/ ) {
            print $ofh ">$1\n";
            $header_found = 1;
        } else {
            print $ofh $_;
        }
    }
    
    die "header never found in sequence $infile" unless ($header_found);
}

exit(0);

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






