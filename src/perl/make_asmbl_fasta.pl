#!/usr/local/bin/perl


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BSML::BsmlReader;
use BSML::BsmlParserTwig;
use File::Basename;

my %options = ();
my $results = GetOptions (\%options, 'bsml_dir|b=s', 'asmbl_ids|a=s', 'output|o=s', 'exclude|e=s', 'help|h', 'each_genome|g' );

###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $ASMBL_IDS        = $options{'asmbl_ids'};
my $output_file      = $options{'output'};
my $exclude_asmbl_id = $options{'exclude'};
my $BSML_dir         = $options{'bsml_dir'};
$BSML_dir =~ s/\/+$//;         #remove terminating '/'s
#my $project = $options{'project'};
my $output_dir = dirname($output_file);


if(!defined($ASMBL_IDS) or !$output_dir or !$BSML_dir or exists($options{'help'})) {
    &print_usage();
}

###-------------------------------------------------------###

my $parser = new BSML::BsmlParserTwig;

if(! -d $output_dir) {
    mkdir $output_dir;
    chmod 0777, $output_dir;
}


my @asm_ids;
if($ASMBL_IDS =~ /all/i) {
    my @files = <$BSML_dir/*.bsml>;
    foreach (@files) {
	my $basename = basename($_);
	if($basename =~ /(.+)\.bsml/) {
	    push(@asm_ids, $1);
        }
    }
} else {
    @asm_ids = split(/,/, $ASMBL_IDS);
}

make_assembly_fasta(\@asm_ids);



sub make_assembly_fasta {

    my $assembly_ids = shift;

    open(FILE, ">$output_file") || die "Can't open $output_file due to $!";
    foreach my $asmbl_id (@$assembly_ids) {
	next if($asmbl_id eq $exclude_asmbl_id);
	my $bsml_file = "$BSML_dir/${asmbl_id}.bsml";
	if (-s $bsml_file) {
	    my $reader = BSML::BsmlReader->new();
	    $parser->parse( \$reader, $bsml_file );
	    my $seq_list = $reader->returnAllSequences();
	    foreach my $seq_obj (@$seq_list) {
		if($seq_obj->returnattr('molecule') eq 'dna') {
		    my $seq_id = $seq_obj->returnattr('id');
		    my $raw_seq = $reader->subSequence($seq_id, -1, -1, '0');
		    my $fastaout = &fasta_out($seq_id, $raw_seq);
		    print FILE $fastaout;
                }
            }
	} else {
	    print STDERR "$bsml_file NOT found!!!! Aborting...\n";
        }
    }
    close FILE;  
    chmod 0666, $output_file;
}

sub fasta_out {
#This subroutine takes a sequence name and its sequence and
#outputs a correctly formatted single fasta entry (including newlines).

    my $seq_name = shift;
    my $seq = shift;

    my $fasta=">"."$seq_name"."\n";
    for(my $i=0; $i < length($seq); $i+=60){
	my $seq_fragment = substr($seq, $i, 60);
	$fasta .= "$seq_fragment"."\n";
    }
    return $fasta;

}



sub print_usage {


    print STDERR "SAMPLE USAGE:  make_asmbl_fasta.pl -b bsml_repository -a bsp_3839_assembly -o output_file\n";
    print STDERR "  --bsml_dir          = bsml repository dir\n";
    print STDERR "  --asmbl_ids         = asmbl_id to include in fasta file or 'all'\n";
    print STDERR "  --output            = output file in fasta format\n";
    print STDERR "  --exclude           = asmbl_id to exclude (optional)\n";
    print STDERR "  --help = This help message.\n";
    exit 1;


}
