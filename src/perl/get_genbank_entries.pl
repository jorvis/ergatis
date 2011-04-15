#!/usr/bin/perl

=head1 NAME

get_genbank_entries.pl

=head1 SYNOPSIS

USAGE:  get_genbank_entries.pl (--acc I<accession> | --list I<accession_list> | --infile I<input_file>) [--man]

=head1 OPTIONS

=over 10

=item B<--acc>

Single accession

=item B<--list>

List of accessions

=item B<--infile>

File with list of accessions

=item B<--man>

Display the pod2usage page for this utility

=back

=head1 DESCRIPTION

get_genbank_entries.pl - Retrieves genbank records from genbank.

=cut

use strict;
use Bio::DB::EUtilities;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);;
use Pod::Usage;

my %options=();
my $acc;
my $list;
my $ifile;
my $help;
my $man;

my @ids;

# Process the options list.
GetOptions('acc=s'    => \$acc,
           'list=s'   => \$list,
           'infile=s' => \$ifile,
           'help|h'  => \$help,
           'man'     => \$man
           );
&process_opts();

# Loop through the IDs and pull down the genbank file.
# Write the results to the accession.gbk.
foreach my $id (@ids) {

    my $efetch = Bio::DB::EUtilities->new(
                                          -db => 'nucleotide',
                                          -id => $id,
                                          -rettype => 'gb',
                                          -retmode => 'text',
                                          );
    
    open OUT, ">$id.gbk" or die "Unable to open $id\n";
    print OUT $efetch->get_Response->content;
    close OUT;
    sleep 3; # Required so NCBI does not get mad.
}

sub process_opts {
    &pod2usage() if ($help || !($acc || $list || $ifile));
    &pod2usage(-exitval => 1, -verbose => 2, -output => \*STDOUT) if $man;
    
    if($acc) {
        push(@ids, $acc);
    }
    if($list) {
        push(@ids, split(',', $list));
    }
    if($ifile) {
        open IN, "<$ifile" or die "Unable to open $ifile\n";
        while(<IN>) {
            chomp;
            push(@ids, $_);
        }
    }
}
