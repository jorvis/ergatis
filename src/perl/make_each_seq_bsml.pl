#!/usr/local/bin/perl


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BSML::BsmlReader;
use BSML::BsmlParserTwig;
use File::Basename;

my %options = ();
my $results = GetOptions (\%options, 'bsml_dir|b=s', 'output_dir|o=s', 'asmbl_ids|a=s', 'DEBUG', 'asmbl_file=s', 'help|h' );

###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $ASMBL_IDS = $options{'asmbl_ids'};
my $output_dir = $options{'output_dir'};
$output_dir =~ s/\/+$//;       #remove terminating '/'s
my $BSML_dir = $options{'bsml_dir'};
$BSML_dir =~ s/\/+$//;       #remove terminating '/'
my $QUERYPRINT;
my $DEBUG = $options{'DEBUG'} || 0;

my $asmbl_file      = $options{'asmbl_file'};

if(!$BSML_dir or !$output_dir or exists($options{'help'})) {
    &print_usage();
}
if(!$asmbl_file and ! $ASMBL_IDS) {
    print STDERR "Either --asmbl_ids  OR --asmbl_file option is needed\n";
    &print_usage();
}

if($asmbl_file and $ASMBL_IDS) {
    print STDERR " Specify either --asmbl_ids OR --asmbl_file\n"; 
    &print_usage();
}

###-------------------------------------------------------###
my $min_dir = dirname($output_dir);
if(! -d $min_dir) {
    mkdir $min_dir;
    chmod 0777, $min_dir;
}
#create the directory that will hold all the individual peptide files if it doesn't exist
if(! -d $output_dir ) {
    mkdir $output_dir;
}
chmod 0777, $output_dir;

my $result;


my @asm_ids;
if($asmbl_file) {   #asmbl_id will be read from a flat file
    @asm_ids = read_asmbl_file($asmbl_file);
    if(!@asm_ids) {
	print STDERR "No asmbl_ids found in $asmbl_file.  Aborting...\n";
	exit 4;
    }
}else {
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
}




my $parser = new BSML::BsmlParserTwig;



foreach my $asmbl_id (@asm_ids) {
    my $final_output_dir = "$output_dir/".$asmbl_id;
    if(! -d $final_output_dir ) {
	mkdir $final_output_dir;
    } else {
	unlink glob("$final_output_dir/*");
    }
    chmod 0777, $final_output_dir;
    my $bsml_file = "$BSML_dir/${asmbl_id}.bsml";
    if (-s $bsml_file) {
	my $reader = BSML::BsmlReader->new();
	$parser->parse( \$reader, $bsml_file );
	my $extend_seq = $reader->get_all_cds_dna($asmbl_id) ;
	while( my ($gene, $seq) = each %$extend_seq ) {
	    next if(length($seq) < 1); 
	    $gene =~ s/_\d$//;
	    my $protein_seq_id = $reader->cdsIdtoProteinSeqId($gene);
	    my $gene_file = "$final_output_dir/${protein_seq_id}.seq";
	    open(FILE, ">$gene_file") || die "Unable to write to $gene_file due to $!";
	    my $fastaout = &fasta_out($protein_seq_id, $seq);
	    print FILE $fastaout;
	    close FILE;
	    chmod 0777, $gene_file;
	}
    } else {
	print STDERR "$bsml_file NOT found!!!!\n";
    }
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

sub read_asmbl_file {

    my $file = shift;

    my @asmbl_id_list;

    open (IN, "$file")  or die "Unable to read $file due to $!";
    my $line;
    while($line = <IN>) {
	chomp($line);
	next if($line =~ /^\s*$/);
	push(@asmbl_id_list, $line);
    }
    close IN;

    return @asmbl_id_list;

}


sub print_usage {


    print STDERR "SAMPLE USAGE:  make_each_seq_bsml.pl -b bsml_dir -o output_dir -a bsp_3839_assembly\n";
    print STDERR "  --bsml_dir    = dir containing BSML doc\n";
    print STDERR "  --output_dir  = dir to save output to\n";
    print STDERR "  --asmbl_ids  (multiple values can be comma separated)\n";
    print STDERR "               (-a all  grabs all asmbl_ids)\n";    
    print STDERR "  --asmbl_file  = name of the file containing a list of asmbl_ids\n";
    print STDERR "  --help = This help message.\n";
    exit 1;

}






