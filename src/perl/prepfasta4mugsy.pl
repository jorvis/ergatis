#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
=head1  NAME 

prepfasta4mugsy.pl - make organism specific fasta files for mugsy

=head1 SYNOPSIS

USAGE:  bsml2fasta.pl 
          --input_list
          --mugsy_map
		  --use_polypeptides
          --output_dir
=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;
use IO::File;
use Digest::MD5 qw(md5 md5_hex md5_base64);

#######
## ubiquitous options parsing and logger creation
my %options = ();
my $results = GetOptions (\%options, 
			  'input_list|i=s',
              'mugsy_map|m=s',
			  'output_dir|o=s',
			  'use_polypeptides|p=i',
                          'checksum_orgs:s',
			  'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

&check_parameters(\%options);

# Setting this to true for now.
$options{checksum_orgs}=1;

# First pull out the sequence ID's for each file
open IN, "<$options{input_list}" or die "Unable to open $options{input_list}\n";
my $headers;
while(<IN>) {
    chomp; 
    push(@$headers, `grep -H \\> $_`);
}
close IN;
my $seq_id_to_file;
map {
    chomp; my($file,$id) = split(/:/,$_); 
    $id =~ s/^\>//;
    $id =~ /^(\S+)/;
    $id = $1;
    $seq_id_to_file->{$id} = $file;
}@$headers;

# Next pull out the sequence ID and organism names from the mugsymap file
open IN2, "<$options{mugsy_map}" or die "Unable to open $options{mugsy_map}\n";
my $org_id_to_seq_ids;

while(<IN2>) {
    chomp;
    my @fields = split(/\t/,$_);
# SAdkins - 1/16/17 - Needed a use-case for polypeptides instead of assembly seq IDs
	my $seq_id_field = ($options{use_polypeptides}) ? $fields[5] : $fields[1];
    if(!$org_id_to_seq_ids->{$fields[7]}) {
        $org_id_to_seq_ids->{$fields[7]} = {};
    }
    $org_id_to_seq_ids->{$fields[7]}->{$seq_id_field} =1;
}
close IN2;

my $id_map = "$options{output_dir}/id_map.txt";
open(FW, "> $id_map") or die "Unable to open file $id_map for writing\n";

foreach my $org_id (keys %$org_id_to_seq_ids) {
    my @files;
    map {
		next if $_ eq "-";	# Found bug where empty polypeptide BSML features get passed as a dash
        if (!defined($seq_id_to_file->{$_})) {
            die "Unable to map $_ to a sequence file\n";
        }
        push(@files,$seq_id_to_file->{$_})
    } keys %{$org_id_to_seq_ids->{$org_id}};
    if(!$options{checksum_orgs}) {
        $org_id =~ s/[\/\.\+\:\;,-]//g;
    }

    my $prev_org_id = $org_id;
    if($options{checksum_orgs}) {
        $org_id = md5_hex($org_id);
    }
    print FW $prev_org_id."\t".$org_id."\n";
# Fix to concatenate FASTA files. The statement below joining the file paths and then running one cat command fails in case the number of arguments to cat exceed ~1500
    my $fsa = "$options{output_dir}/$org_id.fsa";
    if(-e $fsa) {
	`rm $fsa`;
    }
    foreach my $f (@files) {
	my $cat = "cat $f >> $options{output_dir}/$org_id.fsa";
	#print STDERR "$cat\n";
    	`$cat`;
    }
#    my $cat = "cat ".join(" ",@files)." > $options{output_dir}/$org_id.fsa";
}

close(FW);
sub check_parameters {
    my ($options) = @_;
        ## they have to pass some form of input
    unless ($options{input_list}) {
        print STDERR "You must specify input with --input_list";
    }
    unless ($options{mugsy_map}) {
        print STDERR "You must specify a mugsymap with --mugsy_map";
    }
}
