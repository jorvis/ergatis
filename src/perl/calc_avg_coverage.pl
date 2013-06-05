#!/usr/local/bin/perl

=head1 NAME

=head1 SYNOPSIS

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use XML::Parser;
use Data::Dumper;
use File::Basename;

my %options;
my $results = GetOptions (\%options,
                          'mpileup=s',
                          'input=s',
                          'help|h'
                          );

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

my $coverages = {};

&read_mpileup();

my $handle;
if($options{input} =~ /.bam/) {
    open($handle, "-|", "/usr/local/bin/samtools view $options{input}") or die "Couldn't open $options{input}\n";
}
else {
    open($handle, "<$options{input}")  or die "Couldn't open $options{input}\n";;
}



while(<$handle>) {
    my @fields = split();

    my $read = $fields[0];
    my $ref = $fields[2];
    my $flag = &parseFlag($fields[1]);

    if($flag->{query_mapped}) {
        my $length = length $fields[9];
        my $start = $fields[3];
        my $stop = $fields[3]+$length;

        my $total = 0;
        for(my $i = $start; $i < $stop;$i++) {
            $total += $coverages->{$ref}->{$i};
        }
        print "$read\t".(sprintf("%.3f",($total/$length)))."\n";
    }

}
sub read_mpileup {

    open IN, "<$options{mpileup}" or die "Unable to open $options{mpileup}\n";
    while(<IN>) {
        chomp;
        my @fields = split();
        $coverages->{$fields[0]}->{$fields[1]} = $fields[3];
    }
    close IN;
}

sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    $str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
    return $str;
}
    
sub parseFlag {
    my $flag = shift;
    my $rawbin = dec2bin($flag);
    my $rev = scalar $rawbin;
    if($rev eq $rawbin) {
        #    print "ERROR $rev $rawbin\n";
    }
    my $bin = sprintf("%011d", $rev); 
    my $final_bin = reverse $bin;
    my $prop = substr($final_bin, 1, 1);
    my $qmap = substr($final_bin, 2, 1);
    my $mmap = substr($final_bin, 3, 1);
    my $qstrand = substr($final_bin, 4, 1);
    my $mstrand = substr($final_bin, 5, 1);

    return {
        'query_mapped' => !$qmap,
        'mate_mapped' => !$mmap
    };
}
