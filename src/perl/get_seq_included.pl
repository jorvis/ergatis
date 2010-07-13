#!/usr/bin/perl -w

#get_seqs_noinlist.pl
#input = list and fasta files

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = GetOptions(\%options,
                         'fasta_file|f:s',
                         'fasta_list|fl:s',
                         'seq_id_list|s=s',
                         'output|o=s' ) || die ("Invalid option");

my $list = $options{'seq_id_list'};
my $fasta = $options{'fasta_file'};
my $fasta_list = $options{'fasta_list'};
my $outfile = $options{'output'};                     

push (@fasta_files, $fasta) if (defined($fasta));
## If a list of fasta files are passed in they need to be parsed 
if (defined($fasta_list)) {
    open (FLIST, $fasta_list) or die ("Could not open FASTA list $fasta_list: $!");
    while (my $file = <FLIST>) {
        chomp($file);
        push (@fasta_files, $file);
    }
    close (FLIST);
}

if ( (scalar @fasta_files) < 1 ) { die ("No input files provided."); }

open LIST, "<$list" or die $!;
my %rm;
while (my $line = <LIST>) {
    chomp($line);
    $rm{$line} = 1;
}

open OUT, "> $outfile" or die $!;
local $/ = "\n>";

foreach my $file (@fasta_files) {
    open (IN, $file) or die ("Could not open FASTA file $file: $!");
    while (my $record = <IN>) {
        chomp($record);
        my @rec = split(/\n/, $record);
        my $header = shift @rec;
        my $seq = join("",@rec);
        next W1 unless $seq;
        $header =~ s/^[\s*|>]//g;
        my ($seqname,@desc) = split(/\s+/, $header);
        if ($rm{$seqname}) {
    	    print OUT ">",$record,"\n";
        }
    }

    close (IN);
}
