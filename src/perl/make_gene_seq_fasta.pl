#!/usr/local/bin/perl


=head1  NAME 

make_gene_seq_fasta.pl  - Generates DNA sequences of genes in FASTA format from BSML files

=head1 SYNOPSIS

USAGE:  make_gene_seq_fasta.pl -b bsml_dir -o output_dir -a bsp_3839_assembly -name PNEUMO

=head1 OPTIONS

=over 4

=item *

B<--bsml_dir,-b>   [REQUIRED] Dir containing BSML documents (repository)

=item *

B<--output_dir,-o> [REQUIRED] dir to save output to

=item *

B<--asmbl_ids,-a> Build sequences from this asmbl_id.  Multiple values can be comma separated

=item *

B<--asmbl_file,-i>  name of the file containing a list of asmbl_ids (1 per line)

=item * 

B<--mode, -m>  mode types: 1 (default mode)  = total_fasta; 2 = fasta/genome; 3 = fasta/gene

=item *

B<--name,-n> name of the DNA seq  fasta file. Required in mode 1 ONLY

=item *

B<--consolidate> enable all output files to be saved in $output_dir. Mode 3 ONLY 

=item *

B<--log,-l> Log file

=item *

B<--debug,--D>  Debug level

=item *

B<--help,-h> This help message

=back

=head1   DESCRIPTION

make_gene_seq_fasta.pl is designed to generate multi-fasta files of DNA sequences of genes
given a list of assembly IDs.  The data is obtained through existing
BSML files in some central repository.  This script operates in possible three modes.  
Mode 1, the default mode, generates ONE fasta file containing all the protein sequences
of the user-specified assembly IDs.  In this mode, the name of the fasta file is specified
by --name flag.  By default, the fasta header is in the form ">protein_id".  Mode 2 generates 
one fasta file per assembly_ID.  The file will be named "assembly_id.seq" 
i.e. (bsp_3839_assembly.pep).  Mode 3 generates one fasta file per gene.  The fasta files 
will be stored in the respective asmbl_id directory.  If --consolidate flag is invoked, 
all the fasta files will be stored in a single directory instead of in each assembly subdirectory.
As an alternative to --asmbl_ids flag, you can also supply a file with a list of asmbl_ids
(1 per line).  If you supply both --asmbl_ids and --asmbl_file.  The script will generated
a unique list of asmbl_ids and fetch sequences from the list.

Samples:

1. make 1 fasta file called PNEUMO.seq containing sequences from assembly A and B   
   make_gene_seq_fasta.pl -b bsml_dir -o output_dir -a A,B -name PNEUMO

2. make N fasta files for N assembly IDs stored in a file "list.txt"
   make_gene_seq_fasta.pl -b bsml_dir -o output_dir --asmbl_file list.txt --mode 2

IMPORTANT:

You must specify at least --asmbl_ids OR --asmbl_file flag.
You must specify --name in mode 1 (default mode)

NOTE:  

Calling the script name with NO flags/options or --help will display the syntax requirement.


=cut


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Log::Log4perl qw(get_logger :levels :easy);
use BSML::BsmlReader;
use BSML::BsmlParserTwig;
use File::Basename;
use File::Path;
use Pod::Usage;


umask(0000);

my %options = ();
my $results = GetOptions (\%options, 'bsml_dir|b=s', 'asmbl_ids|a=s', 'consolidate|s', 'name|n=s', 'verbose|v', 'log|l=s', 'debug|D=s',
                                     'output_dir|o=s', 'help|h', 'mode|m=s', 'asmbl_file|i=s', 'man' ) || pod2usage();

###-------------PROCESSING COMMAND LINE OPTIONS-------------###
my $asmbl_ids       = $options{'asmbl_ids'};
my $output_dir      = $options{'output_dir'};
$output_dir =~ s/\/+$//;       #remove terminating '/'s
my $BSML_dir        = $options{'bsml_dir'};
$BSML_dir =~ s/\/+$//;         #remove terminating '/'s
my $asmbl_file      = $options{'asmbl_file'};
my $fasta_name      = $options{'name'};
#mode types: 1 = total_fasta; 2 = fasta/genome; 3 = fasta/gene
my $mode            = $options{'mode'} || 1;    #default always = 1
my $consolidate     = $options{'consolidate'};


my $parser = new BSML::BsmlParserTwig;

&cmd_check();
###-------------------------------------------------------###

my @asm_ids = determine_all_asmbl_id($asmbl_ids, $asmbl_file, $BSML_dir);


if($mode == 1) {
    make_total_seq_fasta(\@asm_ids, $fasta_name);
}elsif($mode == 2) {
    make_seq_fasta_per_genome(\@asm_ids);
}elsif($mode == 3) {
    make_seq_fasta_for_each_gene(\@asm_ids, $consolidate);
} else {
    print STDERR "Unknown Usage Mode. Aborting...\n";
}



sub cmd_check {
#quality check

    if( exists($options{'man'})) {
	pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT});
    }   

    if( exists($options{'help'})) {
	pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT});
    }

    if(!$output_dir or !$BSML_dir) {
	pod2usage({-exitval => 2,  -message => "$0: All the required options are not specified", -verbose => 1, -output => \*STDERR});
    }
    
    if(!$asmbl_file and ! $asmbl_ids) {
	pod2usage(-exitval =>2, -verbose => 1, -message => "$0: Must specify either --asmbl_ids OR --asmbl_file.", -output => \*STDERR);
    }
    
    if( $mode !~ /[123]/ ) {
	pod2usage(-exitval =>2, -verbose => 1, -message => "$0: valid mode =  1 or 2 or 3", -output => \*STDERR);
    }
    
    if($mode =~ /1/ and !defined($fasta_name) ) {
	pod2usage(-exitval =>2, -verbose => 1, -message => "$0: --name required in mode 1", -output => \*STDERR);
    }


    #checking BSML repository directory
    if(! -d $BSML_dir) {
	print STDERR "$BSML_dir NOT found.  Aborting...\n";
	exit 5;
    }
    #check for user defined output directory
    if(! -d $output_dir) {
	mkpath($output_dir) or die "Unable to create $output_dir.  Aborting...\n";
	#chmod 0777, $output_dir;
    }

}

sub  make_seq_fasta_for_each_gene {

    my $assembly_ids = shift;
    my $consolidate_mode = shift;

    if($consolidate_mode) {
	seq_fasta_per_gene_consolidated($assembly_ids);
    } else {
	seq_fasta_per_gene_regular($assembly_ids);
    }

}


sub seq_fasta_per_gene_consolidated {

    my $assembly_ids = shift;

    my $final_output_dir = "$output_dir";
    if(! -d $final_output_dir ) {
	mkpath($final_output_dir); #creates dir with 777 permission
    } else {         #delete old crap if that directory already exists
	opendir(DIR, $final_output_dir) or die "Unable to access $final_output_dir due to $!";
	my @old_files = readdir(DIR);
	@old_files = grep { !(/^\.{1,2}$/) } @old_files;  #get rid of pesky '.' and ' ..'
	@old_files = map { "$final_output_dir/".$_ } @old_files;
	unlink (@old_files);
    }
    chmod 0777, $final_output_dir;

    foreach my $asmbl_id (@$assembly_ids) {
	my $bsml_file = "$BSML_dir/${asmbl_id}.bsml";
	if (-s $bsml_file) {
	    my $reader = BSML::BsmlReader->new();
	    $parser->parse( \$reader, $bsml_file );
	    my $extend_seq = $reader->get_all_cds_dna($asmbl_id) ;
	    while( my ($gene_id, $seq) = each %$extend_seq ) {
		next if(length($seq) < 1); 
		$gene_id =~ s/_\d$//;
		my $protein_seq_id = $reader->cdsIdtoProteinSeqId($gene_id);
		my $gene_file = "$final_output_dir/${protein_seq_id}.seq";
		open(FILE, ">$gene_file") || die "Unable to write to $gene_file due to $!";
		my $fastaout = &fasta_out($protein_seq_id, $seq);
		print FILE $fastaout;
		close FILE;
		chmod 0777, $gene_file;
	    }
	    #$logger->debug("$output_seq_counter sequences for $asmbl_id written to $final_output_dir ");
	} else {
	    print STDERR "$bsml_file NOT found, skipping...\n";
	    #$logger->logwarn("$bsml_file NOT found!  Skipping...");
	}

    }


}


sub seq_fasta_per_gene_regular {

    my $assembly_ids = shift;

    foreach my $asmbl_id (@$assembly_ids) {
	my $final_output_dir = "$output_dir/".$asmbl_id;
	#create the directory that will hold all the individual peptide files if it doesn't exist
	if(! -d $final_output_dir ) {
	    mkpath($final_output_dir); #creates dir with 777 permission
	}else {         #delete old crap if that directory already exists
	    opendir(DIR, $final_output_dir) or die "Unable to access $final_output_dir due to $!";
	    my @old_files = readdir(DIR);
	    @old_files = grep { !(/^\.{1,2}$/) } @old_files;  #get rid of pesky '.' and ' ..'
	    @old_files = map { "$final_output_dir/".$_ } @old_files;
            unlink (@old_files);
	}
	chmod 0777, $final_output_dir;
	my $bsml_file = "$BSML_dir/${asmbl_id}.bsml";
	if (-s $bsml_file) {
	    my $reader = BSML::BsmlReader->new();
	    $parser->parse( \$reader, $bsml_file );
	    my $extend_seq = $reader->get_all_cds_dna($asmbl_id) ;
	    while( my ($gene_id, $seq) = each %$extend_seq ) {
		next if(length($seq) < 1); 
		$gene_id =~ s/_\d$//;
		my $protein_seq_id = $reader->cdsIdtoProteinSeqId($gene_id);
		my $gene_file = "$final_output_dir/${protein_seq_id}.seq";
		open(FILE, ">$gene_file") || die "Unable to write to $gene_file due to $!";
		my $fastaout = &fasta_out($protein_seq_id, $seq);
		print FILE $fastaout;
		close FILE;
		chmod 0777, $gene_file;
	    }
	    #$logger->debug("$output_seq_counter sequences for $asmbl_id written to $final_output_dir ");
	} else {
	    print STDERR "$bsml_file NOT found, skipping...\n";
	    rmtree($final_output_dir);  #delete dir
	    #$logger->logwarn("$bsml_file NOT found!  Skipping...");
	}
    }

}

sub make_seq_fasta_per_genome {
#this mode outputs 1 dna fasta file per assembly ID


    my $assembly_ids = shift;

    #$logger->debug("generating fasta files: 1 assembly per file");
    
    foreach my $asmbl_id (@$assembly_ids) {
	my $seq_file = "$output_dir/$asmbl_id".".seq";
	open(FASTA, ">$seq_file") || die "Unable to write to $seq_file due to $!";
	my $bsml_file = "$BSML_dir/${asmbl_id}.bsml";
	if (-s $bsml_file) {
	    my $reader = BSML::BsmlReader->new();
	    $parser->parse( \$reader, $bsml_file );
	    my $extend_seq = $reader->get_all_cds_dna($asmbl_id) ;
	    while( my ($gene_id, $seq) = each %$extend_seq ) {
		next if(length($seq) < 1); 
		$gene_id =~ s/_\d$//;
		my $protein_seq_id = $reader->cdsIdtoProteinSeqId($gene_id);		
		my $fastaout = &fasta_out($protein_seq_id, $seq);
                print FASTA $fastaout;
            }
	    close FASTA;
	    chmod 0777, $seq_file;
	}else {
	    print STDERR "$bsml_file NOT found, skipping...\n";
	    unlink($seq_file) if (-z $seq_file);  #delete empty file
	    #$logger->logdie("$bsml_file NOT found!  Exiting...");
	}
    }

}



sub make_total_seq_fasta {
#this mode outputs 1 dna fasta file for ALL assembly IDs

    my $assembly_ids = shift;
    my $fasta_name   = shift;

    my $seq_file = "$output_dir/$fasta_name.seq";
    open(FASTA, ">$seq_file") || die "Cant open $seq_file for writing due to $!";
    foreach my $asmbl_id (@$assembly_ids) {
	my $bsml_file = "$BSML_dir/${asmbl_id}.bsml";
	if (-s $bsml_file) {
	    my $reader = BSML::BsmlReader->new();
	    $parser->parse( \$reader, $bsml_file );
	    my $extend_seq = $reader->get_all_cds_dna($asmbl_id) ;
	    while( my ($gene_id, $seq) = each %$extend_seq ) {
		next if(length($seq) < 1); 
		$gene_id =~ s/_\d$//;
		my $protein_seq_id = $reader->cdsIdtoProteinSeqId($gene_id);		
		my $fastaout = &fasta_out($protein_seq_id, $seq);
                print FASTA $fastaout;
            }
	} else {
	    print STDERR "$bsml_file NOT found, skipping...\n";
            #$logger->logdie("$bsml_file NOT found!  Exiting...");
        }
    }
    close FASTA;
    unlink($seq_file) if (-z $seq_file);  #delete empty file
    chmod 0777, $seq_file;

}





sub determine_all_asmbl_id {
#given $asmbl_ids and/or $asmbl_file
#return a list of asmbl_ids
#if both arguments are supplied, will return a list of
#asmbl_ids (duplicates are removed)
    
    my $asmbl_ids  = shift;
    my $asmbl_file = shift;
    my $bsml_dir   = shift;

    my @final_asmbl_ids;
    my (%unique_asmbl_ids, @id_list);

    
    #check asmbl_id passed via command line option
    if(defined($asmbl_ids)) {
	if($asmbl_ids =~ /^all$/) {
	    @id_list = get_all_asmbl_id_via_all($bsml_dir);
	} else{
	    @id_list = split(/,/, $asmbl_ids);
	}
	foreach my $asmbl_id (@id_list) {
	    push(@final_asmbl_ids, $asmbl_id) unless $unique_asmbl_ids{$asmbl_id}++;
	}
    }
       
    #check asmbl_id from a flat file
    if(defined($asmbl_file)) {
	@id_list = read_asmbl_file($asmbl_file);
	foreach my $asmbl_id (@id_list) {
	    push(@final_asmbl_ids, $asmbl_id) unless $unique_asmbl_ids{$asmbl_id}++;
	}
    }

    return @final_asmbl_ids;

}




sub read_asmbl_file {
#parses the flat file containing
#a list of asmbl_ids

    my $file = shift;

    my @asmbl_id_list;

    open (IN, "$file") or die("Unable to read $file due to $!");
    my $line;
    while($line = <IN>) {
	chomp($line);
	next if($line =~ /^\s*$/);
	$line =~ s/^\s+//;
	$line =~ s/\s+$//;
	push(@asmbl_id_list, $line);
    }
    close IN;

    return @asmbl_id_list;

}




sub get_all_asmbl_id_via_all {
#shortcut way of getting all asmbl_ids by user
#specifying --asmbl_id all on command line
#The asmbl_ids are obtained by removing the suffix of 
#the bsml file

    my $bsml_dir = shift;

    my @asm_ids;

    opendir(DIR, $bsml_dir) or die "Unable to access $bsml_dir due to $!";
    while( my $file = readdir(DIR)) {
	next if ($file =~ /^\.{1,2}$/);  #skip  "." ,  ".."
	if($file =~ /(.+)\.bsml$/) {
	    push(@asm_ids, $1);
	}
    }

    return  @asm_ids;

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
