#!/usr/local/bin/perl

=head1  NAME 

gene_pos_bsml.pl  - Generates gene position XML files from BSML files

=head1 SYNOPSIS

USAGE:  gene_pos_bsml.pl -b bsml_dir -a bsp_3839_assembly > out.xml

=head1 OPTIONS

=over 4

=item *

B<--bsml_dir,-b>   [REQUIRED] Dir containing BSML documents (repository)

=item *

B<--asmbl_ids,-a> Build sequences from this asmbl_id.  Multiple values can be comma separated

=item *

B<--asmbl_file,-i>  name of the file containing a list of asmbl_ids (1 per line)

=item * 

B<--min_num_gene>  minimum number of genes an assembly must have (default: 0)

=item *

B<--min_size_gene> minimum length of a gene to be considered  (default: 0)

=item *

B<--log,-l> Log file

=item *

B<--debug,--D>  Debug level

=item *

B<--help,-h> This help message

=back

=head1   DESCRIPTION

gene_pos_bsml.pl is designed to generate XML output file for PEffect program.
The data is obtained through existing BSML files in some central repository.
User can specify the sequences to fetch by --asmbl_ids flag AND/OR --asmbl_file 
flag.  Using the first options, asmbl_ids can be comma separated while the 
latter option requires the asmbl_ids to be stored 1 per line in a file.  
If both options are used, a unique list of asmbl_ids from both sources will
be fetched.  The output is directed to STDOUT.  As an option, the user
can specify the minimum number of genes an asmbl_id must have in order to
be considered via the --min_num_gene flag.  Similarly, the --min_size_gene
denotes the minimum length of a gene must be in order to be considered.
Both options above defaults to zero.

Samples:

1. make an xml file for assembly sequence A and B
   gene_pos_bsml.pl -b bsml_dir -a A,B > output.xml 

IMPORTANT:

You must specify at least --asmbl_ids OR --asmbl_file flag.

NOTE:  

Calling the script name with NO flags/options or --help will display the syntax requirement.


=cut

use strict;
use PEffect::PEffectXML;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BsmlCGCReader;
use BSML::BsmlParserTwig;
use File::Basename;
use Pod::Usage;
use File::Path;


umask(0000);
my %options = ();
my $results = GetOptions (\%options, 'asmbl_ids|a=s', 'bsml_dir|b=s', 'min_num_gene=s', 'man', 
                                     'min_size_gene=s', 'asmbl_file|i=s', 'help|h' ) || pod2usage();
###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $asmbl_ids       = $options{'asmbl_ids'};
my $BSML_dir        = $options{'bsml_dir'};
$BSML_dir =~ s/\/$//;
my $output          = $options{'output'};
my $gene_cutoff     = $options{'min_num_gene'} || '0';
my $gene_size_cutoff = $options{'min_size_gene'} || '0';
my $asmbl_file      = $options{'asmbl_file'};

&cmd_check();
###-------------------------------------------------------###

my $parser = BSML::BsmlParserTwig->new();
my $pexml =  PEffect::PEffectXML->new();


my @asm_ids = determine_all_asmbl_id($asmbl_ids, $asmbl_file, $BSML_dir);

foreach my $asmbl_id (@asm_ids){
    addGenes($asmbl_id, "db_${asmbl_id}", $pexml);
}


#print out the entire xml to STDOUT
my($oref);
$pexml->outputXML(\$oref);
print $oref,"\n";

#----------------------------------------------------------------

sub addGenes {
#This subroutine grabs all the genes' length, orientation and starting coordinate
#for a given asmbl_id and uses PEffectXML module to make xml

    my $asmbl_id = shift;
    my $name = shift;
    my $pexml = shift;

    my $geneID_protID={};
    my $bsml_file = "$BSML_dir/${asmbl_id}.bsml";
    if (-s $bsml_file) {
	my $reader = BsmlCGCReader->new();
	$parser->parse( \$reader, $bsml_file );
	my $rhash = $reader->returnAllIdentifiers();
	$geneID_protID = build_geneID_protID_mapping($rhash);
	my $order = 0;
	my $sorted_genes = get_sorted_gene_position($asmbl_id, $reader);
	#print STDERR "There are ", scalar(@$sorted_genes), " genes for $asmbl_id\n";
	return if(@$sorted_genes < $gene_cutoff);  #filter out assemblies whose gene number is below cutoff
	foreach my $gene (@$sorted_genes) {
	    my $attrref={};
	    my $feat_name = $gene->{'feat_name'};
	    #In prok, only 1 protein id maps to a gene id
	    $feat_name = $geneID_protID->{$feat_name}->[0];
	    $attrref->{'length'} = $gene->{'length'}; 
	    $attrref->{'orient'} = $gene->{'orient'}; 
	    $attrref->{'coord'}  = $gene->{'coord'}; 
	    if($gene->{'length'} > $gene_size_cutoff){
		$pexml->addFeature($feat_name, $attrref, $name, $order);
	    }
	    $order++;
	}
    } else {
	print STDERR "$bsml_file does not exist!!! skipping...\n";
    }
	
}


sub build_geneID_protID_mapping {
#This function builds a mapping between geneID to proteinID. 
#For Prok, it will be 1 to 1, however for Euk, there can be
#more than 1 proteinID to each geneID.
#The returned structure is a hash ref, where key is geneID, value is an array ref of proteinID

    my $rhash = shift;

    my $geneID_protID={};
    foreach my $seqID (keys %$rhash) {
	foreach my $geneID (keys %{ $rhash->{$seqID} }) {
	    foreach my $transcriptID (keys %{ $rhash->{$seqID}->{$geneID} }) {
		push (@{ $geneID_protID->{$geneID} }, $rhash->{$seqID}->{$geneID}->{$transcriptID}->{'proteinId'});
	    }
	}
    }

    return $geneID_protID;

}



sub get_sorted_gene_position {

    my $asmbl_id = shift;
    my $reader   = shift;

    my $gene_pos = $reader->fetch_gene_positions($asmbl_id);
    my $array_ref=[];
    foreach (@$gene_pos) {
	foreach my $gene (keys %$_) {
	    my $end5 = $_->{$gene}->{'startpos'};
	    my $end3 = $_->{$gene}->{'endpos'};
	    my $complement = $_->{$gene}->{'complement'} == 0 ? '-' : '+';
	    my $length = abs($end3 - $end5);
	    my $coord  = $end5 > $end3 ? $end5 : $end3;
	    push( @$array_ref, { 'feat_name' => $gene, 'length' => $length, 'orient' => $complement, 'coord' => $coord } ); 
	}
    }
    my @sorted_ref = sort { $a->{'coord'} <=> $b->{'coord'} || $a->{'feat_name'} cmp $b->{'feat_name'} } @$array_ref;
    
    return (\@sorted_ref);

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
	if($asmbl_ids =~ /^all$/i) {
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


sub cmd_check {
#quality check

    if( exists($options{'man'})) {
	pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT});
    }   

    if( exists($options{'help'})) {
	pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT});
    }

    if(!$BSML_dir) {
	pod2usage({-exitval => 2,  -message => "$0: Must specify --bsml_dir", -verbose => 1, -output => \*STDERR});
    }
    
    if(!$asmbl_file and ! $asmbl_ids) {
	pod2usage(-exitval =>2, -verbose => 1, -message => "$0: Must specify either --asmbl_ids OR --asmbl_file.", -output => \*STDERR);
    }
    
    if( $gene_cutoff !~ /^\d+$/ or $gene_size_cutoff !~ /^\d+$/ ) {
	pod2usage(-exitval =>2, -verbose => 1, -message => "$0: min_num_gene and min_size_gene must be positive integers", -output => \*STDERR);
    }
    
    #checking BSML repository directory
    if(! -d $BSML_dir) {
	print STDERR "$BSML_dir NOT found.  Aborting...\n";
	exit 5;
    }

}






