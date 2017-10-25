#!/usr/bin/env perl

use strict;
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use File::Basename;
use Data::Dumper;
use Bio::SeqIO;
use IntervalTree;
use Bio::DB::EUtilities;
use File::Find;

$|++;

my %options = ();
my $results = GetOptions (
    \%options,
    'genbank_file=s',
    'genbank_dir=s',
    'reference=s',
    'bam_file=s',
    'output_dir=s',
    'help|h') || pod2usage();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

# If a genbank file hasn't been passed in then we'll figure out the accession
# looking at the bam header.
if($options{genbank_dir}) {
    &find_gb_file();
}
elsif(!$options{genbank_file}) {
    &get_gb_file();
}


my $tree = IntervalTree->new();
#First read in the genbank file
&read_genbank();
$tree->buildTree();

&read_bam();

sub get_ref_id {
    my @lines = `samtools view -H $options{bam_file}`;

    my $id;
    map {
        if(/^\@SQ\s+SN:gi\|\d+\|\w+\|(\S+)\|/) {
            $id = $1;
        }
    }@lines;
    return $id;
}

sub find_gb_file {
    my $file;

    my $id = &get_ref_id();
    $id =~ s/\..*//;
    print STDERR "Looking for $id.gbk in $options{genbank_dir}\n";
    find(sub { 
        if($File::Find::name =~ /$id.gbk/) { 
            $file = $File::Find::name;
       }},$options{genbank_dir});
    $options{genbank_file} = $file;
}

sub get_gb_file {

    my $id = &get_ref_id();

    my $efetch = Bio::DB::EUtilities->new(
                                          -db => 'nucleotide',
                                          -id => $id,
                                          -rettype => 'gb',
                                          -retmode => 'text',
                                          );
    
    open OUT, ">$options{output_dir}/$id.gbk" or die "Unable to open $options{output_dir}/$id.gbk\n";
    print OUT $efetch->get_Response->content;
    close OUT;

    $options{genbank_file} = "$options{output_dir}/$id.gbk";
}

sub read_bam {
    open(my $handle,"-|", "samtools view $options{bam_file}");

    while(my $line = <$handle>) {
        my @fields = split(/\t/,$line);
        my $flag = &parseFlag($fields[1]);
        if($flag->{query_mapped}) {
            my @overlaps = $tree->searchInterval($fields[3], $fields[3]+$fields[9]);
            if(@overlaps) {
                foreach my $p (@overlaps) {
                    print "$fields[0]\t$p->[2]->{primary}\t$p->[2]->{product}\n";
                }
            }
        }
    }
}

sub read_genbank {

    my $file = "$options{genbank_file}";

    my $stream = Bio::SeqIO->new(-file => $file,
                              -format => 'GenBank');
    
    while(my $seq = $stream->next_seq()) {
        my @feats = $seq->get_all_SeqFeatures();

        foreach my $feat (@feats) {
            my $primary = $feat->primary_tag();
            # Basically we're only looking for things with a product.
            # This includes rRNA,tRNA and CDS I believe
            if($feat->has_tag('product')) {
                my @vals = $feat->get_tag_values('product');
                my $start = $feat->start;
                my $stop = $feat->end;
                my $obj = {
                    'primary' => $primary,
                    'product' => $vals[0],
                    'start' => $feat->start,
                    'end' => $feat->end
                };
                $tree->addInterval($obj,$feat->start,$feat->end);
#                print "$primary\t$vals[0]\t".$feat->start."\t".$feat->end."\n";
            }
#            foreach my $tag ($feat->get_all_tags) {
#                print "$tag\n";

#            }
        }
        
    }
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
