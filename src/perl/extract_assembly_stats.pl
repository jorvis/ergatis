#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

extract_assembly_stats.pl - Extracts statistics from an assembly file

=head1 SYNOPSIS

USAGE: extract_assembly_stats.pl 
            --input=/path/to/input_file.fa 
            --output=/path/to/output_file.fa
          [ --limit=integer 
            --no50base
	    --help
	    --log
          ]

=head1 OPTIONS

B<--input,-i>
    Required.Input assembly file to be analyzed. Fasta format.
    Files compressed with gzip or bzip2 are automatically decompressed.

B<--output,-o>
    Required. Output file with statistics. Tab-delimited format.

B<--limit,-m>
    Optional. Script ignores contigs smaller than the integer assigned. 

B<--no50base,-n> 
    Optional. Script will use this number as the base genome size for computing the N50 value.

B<--log,-l> 
    Log file

B<--help,-h>
    This help message


=head1  DESCRIPTION

This script extracts statistics from an assembly (fasta) file - contigs or scaffolds.

=head1  INPUT

Input format is fasta. Files compressed with gzip or bzip2 are automatically decompressed.

=head1  OUTPUT

Output is a tab-delimited file with the following information:
File: name of input file
Number: total number of contigs
Total Size: total size of contigs
Min Size: minimum contig size
Max Size: maximum contig size
Average Size: average contig size
N50: size of contig c, such that 50% of the total assembly size is contained in contigs larger than c
Size @ 1Mbp: same as N50 but assuming that genome size is 2Mbp
Number @ 1Mbp: smallest number of contigs that add up to 1Mbp
Size & Number @ 2Mbp, 4Mbp, 10Mbp - same as above but for more values 


=head1  CONTACT

    Kemi Abolude
    kabolude@som.umaryland.edu
    Script adapted from http://www.cbcb.umd.edu/software/metamos/statistics.pl

=cut


use Bio::SeqIO;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Statistics::Descriptive;
use strict;
use Ergatis::Logger;


my %options = ();
my $results = GetOptions (\%options,
			  'input|i=s',
			  'output|o=s',
 			  'limit|m:i',
			  'n50base|n=i',
			  'log|l=s',
			  'debug=s',
			  'help|h') || pod2usage();

my $ifile = $options{'input'};
my $ofile = $options{'output'};
my $limit = $options{'limit'};
my $N50Base = $options{'n50base'};

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();


## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

###code
open(OUT, ">$ofile");

print OUT "File:\t$ifile\n";
    
if ($ifile =~ /\.gz$/){ # gzipped file
    	open(IN, "gzip -dc $ifile |") || $logger->logdie ("Cannot open $ifile\n");
   }elsif ($ifile =~ /\.bz2$/) { #bzipped file
   	open(IN, "bzip2 -dc $ifile |") || $logger->logdie ("Cannot open $ifile\n");
   } else {
   	open(IN, "$ifile") || $logger->logdie ("Cannot open $ifile\n");
   }

my $fr = Bio::SeqIO->new(-fh=>\*IN, -format=> 'fasta');
 
if (! defined $fr){ die ("Cannot parse file\n");}

my @sizes = ();
my $stats = new Statistics::Descriptive::Full;

while (my $seq = $fr->next_seq){
      if (defined $limit && length($seq->seq()) < $limit) { next;}
	push @sizes, length($seq->seq());
	$stats->add_data(length($seq->seq()));
}

close(IN);

print OUT "Number:\t", $stats->count(), "\n"; # count 
print OUT "Total size:\t", $stats->sum(), "\n"; # total size
print OUT "Min size:\t", $stats->min(), "\n"; # min size
print OUT "Max size:\t", $stats->max(), "\n"; # max size
print OUT "Mean size:\t", sprintf("%.2f\n",$stats->mean()); # mean size
print OUT "Median size:\t", $stats->median(), "\n"; # median size
    
# now for the N* statistics
if (! defined $N50Base){
	$N50Base = $stats->sum() / 2;
}

my $total = 0;
my @sizes = sort {$b <=> $a} @sizes;

my $n50 = undef;
my $n1m = undef;
my $n2m = undef;
my $n4m = undef;
my $n10m = undef;
my $s1m = undef;
my $s2m = undef;
my $s4m = undef;
my $s10m = undef;

    for (my $i = 0; $i <= $#sizes; $i++){
	$total += $sizes[$i];
	if ($total >= $N50Base && ! defined $n50){
	    $n50 = $sizes[$i];
	}
	if ($total >= 1000000 && ! defined $n1m){
	    $n1m = $i + 1;
	    $s1m = $sizes[$i];
	}
	if ($total >= 2000000 && ! defined $n2m){
	    $n2m = $i + 1;
	    $s2m = $sizes[$i];
	}
	if ($total >= 4000000 && ! defined $n4m){
	    $n4m = $i + 1;
	    $s4m = $sizes[$i];
	}
	if ($total >= 10000000 && ! defined $n10m){
	    $n10m = $i + 1;
	    $s10m = $sizes[$i];
	}
    }

print OUT  "N50:\t$n50\nSize @ 1Mbp:\t$s1m\nNumber @ 1Mbp:\t$n1m\nSize @ 2Mbp:\t$s2m\nNumber @ 2Mbp:\t$n2m\nSize @ 4Mbp:\t$s4m\nNumber @ 4Mbp:\t$n4m\nSize @ 10Mbp:\t$s10m\nNumber @ 10Mbp:\t$n10m\n\n";

print OUT "File: name of input file
Number: total number of contigs
Total Size: total size of contigs
Min Size: minimum contig size
Max Size: maximum contig size
Average Size: average contig size
N50: size of contig c, such that 50% of the total assembly size is contained in contigs larger than c
Size @ 1Mbp: same as N50 but assuming that genome size is 2Mbp
Number @ 1Mbp: smallest number of contigs that add up to 1Mbp
Size & Number @ 2Mbp, 4Mbp, 10Mbp - same as above but for more values\n"; 

close(OUT);

exit(0);

#Subroutines
sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( input output );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            $logger->logdie ("--$option is a required option");
        }
    }

    ## make sure input file exists
    if (! -e $options{input} ) {
        $logger->logdie("the input file passed ($options{input}) cannot be read or does not exist");
    }
}

sub _pod {
# Used to display the perldoc help.
    pod2usage( {-exitval => 0, -verbose => 2} );
}

 
