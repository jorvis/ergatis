#!/usr/local/bin/perl
=head1  NAME 

prepfasta4mugsy.pl - make organism specific fasta files for mugsy

=head1 SYNOPSIS

USAGE:  bsml2fasta.pl 
          --input_list
          --mugsy_map
          --output_dir
=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;
use IO::File;

#######
## ubiquitous options parsing and logger creation
my %options = ();
my $results = GetOptions (\%options, 
			  'input_list|i=s',
              'mugsy_map|m=s',
			  'output_dir|o=s',
			  'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

&check_parameters(\%options);

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
    if(!$org_id_to_seq_ids->{$fields[7]}) {
        $org_id_to_seq_ids->{$fields[7]} = {};
    }
    $org_id_to_seq_ids->{$fields[7]}->{$fields[1]} =1;
}
close IN2;

foreach my $org_id (keys %$org_id_to_seq_ids) {
    my @files;
    map {
        if (!defined($seq_id_to_file->{$_})) {
            die "Unable to map $_ to a sequence file\n";
        }
        push(@files,$seq_id_to_file->{$_})
    } keys %{$org_id_to_seq_ids->{$org_id}};
    $org_id =~ s/[\/\.]//g;
    my $cat = "cat ".join(" ",@files)." > $options{output_dir}/$org_id.fsa";
    `$cat`;
}


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
