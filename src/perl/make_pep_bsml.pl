#!/usr/local/bin/perl


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BSML::BsmlReader;
use BSML::BsmlParserTwig;
use File::Basename;

my %options = ();
my $results = GetOptions (\%options, 'bsml_dir|b=s', 'asmbl_ids|a=s', 'simple_header|s', 'project|p=s', 'verbose|v',
                                     'output_dir|o=s', 'DEBUG', 'help|h', 'each_file|e', 'each_genome|g' );

###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $ASMBL_IDS       = $options{'asmbl_ids'};
my $simple_header   = $options{'simple_header'};
my $each_file       = $options{'each_file'} || 0;
my $each_genome     = $options{'each_genome'} || 0;
my $project         = $options{'project'};
my $verbose         = $options{'verbose'};
my $output_dir      = $options{'output_dir'};
$output_dir =~ s/\/+$//;       #remove terminating '/'s
my $BSML_dir        = $options{'bsml_dir'};
$BSML_dir =~ s/\/+$//;         #remove terminating '/'s


if(!defined($ASMBL_IDS) or !$output_dir or !$BSML_dir or exists($options{'help'})) {
    &print_usage();
}
if($each_file and $each_genome) {
    print STDERR "--each_genome(-g) and --each_file(-e) CANNOT be invoked at the same time.\n";
    exit 1;
}

if(!$each_file and !$each_genome and !$project) {
    print STDERR "Must specify project name without --each_genome(-g) and --each_file(-e).\n";
    exit 1;
}


###-------------------------------------------------------###
my $min_dir = dirname($output_dir);
if(! -d $min_dir) {
    mkdir $min_dir;
    chmod 0777, $min_dir;
}

if(! -d $output_dir ) {
    mkdir $output_dir;
}
chmod 0777, $output_dir;

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

my $parser = new BSML::BsmlParserTwig;

if($each_genome) {
    make_fasta_for_each_genome(\@asm_ids) if($verbose);
    print "Finished making fasta pep files for EACH genome\n";
}elsif($each_file) {
    make_fasta_for_each_gene(\@asm_ids) if($verbose);
    print "Finished making fasta pep file for each gene\n";
}else {
    make_PNEUMO_pep_for_ALL_genomes(\@asm_ids);
    print "Finished making fasta pep file for all assemblies\n" if($verbose);
}



sub make_fasta_for_each_gene  {

    my $assembly_ids = shift;

    foreach my $asmbl_id (@$assembly_ids) {
	next if($asmbl_id !~ /\d+/);
	my $final_output_dir = "$output_dir/".$asmbl_id;
	#create the directory that will hold all the individual peptide files if it doesn't exist
	if(! -d $final_output_dir ) {
	    mkdir $final_output_dir;
	}else {         #delete old crap if that directory already exists
	    unlink glob("$final_output_dir/*");
	}
	chmod 0777, $final_output_dir;
	my $bsml_file = "$BSML_dir/${asmbl_id}.bsml";
	if (-s $bsml_file) {
	    my $reader = BSML::BsmlReader->new();
	    $parser->parse( \$reader, $bsml_file );
	    my $proteins = $reader->get_all_protein_aa($asmbl_id);
	    while (my ($identifier, $seq) = each %$proteins) {
		$identifier =~ s/_\d$//;
		my $pep_file = "$final_output_dir/$identifier.fsa";
		open(FILE, ">$pep_file") || die "Unable to write to $pep_file due to $!";
		my $prot_name = $identifier;
		my $fastaout = &fasta_out($prot_name, $seq);
		print FILE $fastaout;
		close FILE;
		chmod 0777, $pep_file;
	    }
	} else {
	    print STDERR "The $bsml_file does not exist!\n";
	    exit 5;
	}

    }

}

sub make_fasta_for_each_genome {

    my $assembly_ids = shift;

    foreach my $asmbl_id (@$assembly_ids) {
	my $pep_file = "$output_dir/$asmbl_id".".pep";
	open(FILE, ">$pep_file") || die "Unable to write to $pep_file due to $!";
	my $bsml_file = "$BSML_dir/${asmbl_id}.bsml";
	if (-s $bsml_file) {
	    my $reader = BSML::BsmlReader->new();
	    $parser->parse( \$reader, $bsml_file );
	    print STDERR "Search protein with ID: chado_pneumo_${asmbl_id}\n";
	    my $proteins = $reader->get_all_protein_aa($asmbl_id);
	    while (my ($identifier, $seq) = each %$proteins) {
		$identifier =~ s/_\d$//;	
		my $prot_name = $identifier;
		my $fastaout = &fasta_out($prot_name, $seq);
                print FILE $fastaout;
            }
	    close FILE;
	    chmod 0777, $pep_file;
	}else {
	    print STDERR "$bsml_file NOT found!!!! Aborting...\n";
	    exit 5;
	}
    }

}


sub make_PNEUMO_pep_for_ALL_genomes {

    my $assembly_ids = shift;

    my $pep_file = "$output_dir/$project.pep";
    open(FILE, ">$pep_file") || die "Cant open $pep_file due to $!";
    foreach my $asmbl_id (@$assembly_ids) {
	my $bsml_file = "$BSML_dir/${asmbl_id}.bsml";
	if (-s $bsml_file) {
	    my $reader = BSML::BsmlReader->new();
	    $parser->parse( \$reader, $bsml_file );
	    my $proteins = $reader->get_all_protein_aa($asmbl_id);
	    while (my ($identifier, $seq) = each %$proteins) {
		$identifier =~ s/_\d$//;
		my $prot_name = $simple_header ? "$identifier $asmbl_id" : "${asmbl_id}:$identifier";
		my $fastaout = &fasta_out($prot_name, $seq);
                print FILE $fastaout;
            }
	} else {
	    print STDERR "$bsml_file NOT found!!!! Aborting...\n";
        }
    }
    close FILE;
    chmod 0777, $pep_file;
    #qx(setdb $pep_file);
    #chmod 0777 <$pep_file.*>;
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


    print STDERR "SAMPLE USAGE:  make_pep_bsml.pl -b bsml_dir -o output_dir -a bsp_3839_assembly\n";
    print STDERR "  --bsml_dir    = dir containing BSML doc\n";
    print STDERR "  --output_dir  = dir to save output to\n";
    print STDERR "  --asmbl_ids  (multiple values can be comma separated)\n";
    print STDERR "               (-a all  grabs all asmbl_ids)\n";
    print STDERR "  --each_genome = save fasta file individually for each genome\n";
    print STDERR "  --each_file   = save fasta file individually for each gene\n";
    print STDERR "  --project     = name of the total peptide fasta file\n";
    print STDERR " NOTE* --each_genome and --each_file cannot be invoked concurrently\n";
    print STDERR "  --help = This help message.\n";
    exit 1;

}
