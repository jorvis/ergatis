#!/usr/local/bin/perl

=head1  NAME 

gene_match_bsml.pl  - Generates pairwise matches for PEffect from BSML files

=head1 SYNOPSIS

USAGE:  gene_match_bsml.pl -d bsml_dir -a bsp_3839_assembly -b gbs_799_assembly -f bsml_btab_dir > out.xml

=head1 OPTIONS

=over 4

=item *

B<--bsml_dir,-d>   [REQUIRED] Dir containing BSML documents (repository)

=item *

B<--bsml_btab_dir,-f> [REQUIRED] Dir containing BSML documents from allvsall

=item *

B<--asmbl_id,-a> [REQUIRED] Build sequences from this asmbl_id.  Multiple values can be comma separated

=item *

B<--asmbl_file,-i>  name of the file containing a list of asmbl_ids (1 per line)

=item * 

B<--match_asmbl_id, -b> [REQUIRED] either a single asmbl_id or 'all'  

=item * 

B<--gene_pos_check> turn on gene_pos file check

=item * 

B<--gene_pos_file> name of the gene position file

=item *

B<--log,-l> Log file

=item *

B<--debug,--D>  Debug level

=item *

B<--help,-h> This help message

=back

=head1   DESCRIPTION

gene_match_bsml.pl is designed to generate pairwise matches in XML for PEffect program.
The data is obtained through existing BSML files in some central repository.
User can specify the sequences to fetch by --asmbl_id flag AND/OR --asmbl_file 
flag.  Using the first options, asmbl_ids can be comma separated while the 
latter option requires the asmbl_ids to be stored 1 per line in a file.  If user supplies
--asmbl_file, ONLY the asmbl_ids in the file will be used.  However, the --asmbl_id flag 
MUST also be invoked with a name for the group of asmbl_ids (i.e. group_1). The reason is
that there must also be a bsml file containing allvsall data with all the asmbl_ids.
The file name should correspond with the --asmbl_id flag (i.e. group_1.allvsall.bsml)  
The output is directed to STDOUT.  As an option, the user can invoke --gene_pos_check.  
This allows the script to examine the gene position xml file and only produce a 
gene_match xml file that contain corresponding information found in the gene_position
xml file.  Since both files are needed to run PEffect, this will ensure that 
both files are good.  In order to use this feature, the user must also specify the name 
of the gene_position file using the --gene_pos_file flag.

Samples:

1. make 1 pairwise xml file between A and B
   gene_match_bsml.pl -d bsml_dir -a A -b B -f bsml_btab_dir > out

2. make 1 xml using file contain asmbl_ids vs 'all'
   *gene_match_bsml.pl -b bsml_dir -i file.txt -b all -a group_A > out
   
   *This assumes there is a file called group_A.allvsall.bsml 

IMPORTANT:

You must specify at least --asmbl_ids OR --asmbl_file flag.
You must specify --asmbl_id if --asmbl_file is invoked.

NOTE:  

Calling the script name with NO flags/options or --help will display the syntax requirement.


=cut


use strict;
use PEffect::PEffectXML;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BsmlCGCReader;
use BSML::BsmlParserTwig;
use Pod::Usage;


my %options = ();
my $results = GetOptions (\%options, 'asmbl_id|a=s', 'verbose|v', 'bsml_btab_dir|f=s', 'gene_pos_file=s', 'man',
                                     'match_asmbl_id|b=s',  'gene_pos_check', 'bsml_dir|d=s', 'asmbl_file|i=s', 'help|h' )|| pod2usage();

###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $ASMBL_IDS     = $options{'asmbl_id'};
my $bsml_btab_dir = $options{'bsml_btab_dir'};
$bsml_btab_dir    =~ s/\/$//;
my $match_asmbl_id = $options{'match_asmbl_id'};
my $verbose       = $options{'verbose'};
my $check_gene_pos = $options{'gene_pos_check'};
my $gene_pos_file  = $options{'gene_pos_file'};
my $BSML_dir = $options{'bsml_dir'};
$BSML_dir =~ s/\/$//;
my $asmbl_file      = $options{'asmbl_file'};

&cmd_check();

###-------------------------------------------------------###
my $valid_asmbl_ids = fetch_valid_asmbl_id($gene_pos_file) if($check_gene_pos);
my ($id_protID, $proteinID_seqID) = build_id_lookups($BSML_dir);

my @asm_ids = determine_all_asmbl_id($ASMBL_IDS, $asmbl_file, $BSML_dir);


my $pexml = new PEffect::PEffectXML();
my $reader;

foreach my $asmbl_id (@asm_ids) {
    print STDERR "processing $asmbl_id\n";
    my @matchfile;
    if($asmbl_file){
	print STDERR "Entering the asmbl_file realm\n";
	my $bsml_btab_file = "$bsml_btab_dir/${ASMBL_IDS}.allvsall.bsml";
	if (!-s $bsml_btab_file) {
	    print STDERR "The $bsml_btab_file does not exist!  Aborting...\n";
	    exit 5;
	} elsif(!$reader) {   #don't reparse same file if already parsed
	    $reader = BsmlCGCReader->new(); 
	    my $parser = BSML::BsmlParserTwig->new();
	    $parser->parse( \$reader, $bsml_btab_file );
	}
    } else {
	print STDERR "Enter the asmbl_id realm\n";
	my $bsml_btab_file = "$bsml_btab_dir/${asmbl_id}.allvsall.bsml";
	if (!-s $bsml_btab_file) {
	    print STDERR "The $bsml_btab_file does not exist!  Skipping...\n";
	    next;
	} else {
	    $reader = BsmlCGCReader->new(); 
	    my $parser = BSML::BsmlParserTwig->new();
	    $parser->parse( \$reader, $bsml_btab_file );
	}
    }
    addMatches($asmbl_id, $match_asmbl_id, $pexml, $reader);
    
}

my($oref);
$pexml->outputXML(\$oref);
print $oref,"\n";




sub addMatches {

    my($asmbl_id, $match_asmbl_id, $pexml, $reader) = @_;

    my $lref;
    if($match_asmbl_id ne 'all') {
	$lref = $reader->fetch_genome_pairwise_matches( $asmbl_id, $match_asmbl_id );
    } else {
	$lref = $reader->fetch_genome_pairwise_matches( $asmbl_id, 'all');
    }

    if($check_gene_pos) {
	if(!exists($valid_asmbl_ids->{$asmbl_id})) {
	    print STDERR "asmbl_id $asmbl_id NOT in gene position xml.  Skipping...\n";
	    return;
        }
	if($match_asmbl_id ne 'all') {
	    if(!exists($valid_asmbl_ids->{$match_asmbl_id})) {
		print STDERR "match_asmbl_id $match_asmbl_id is NOT in gene position xml. Skipping...\n";
		return;
	    }
        }
    }
	    
    
    foreach my $match (@$lref) {
	my $q_feat_name = $match->{'query_gene_name'};
	$q_feat_name = $id_protID->{$q_feat_name} if(exists $id_protID->{$q_feat_name});
        my $m_feat_name = $match->{'match_gene_name'};
	$m_feat_name = $id_protID->{$m_feat_name} if(exists $id_protID->{$m_feat_name});
	my $per_sim     = $match->{'percent_similarity'};
	my $per_id      = $match->{'percent_identity'};
        my $pvalue      = $match->{'pval'};

	#If check_gene_pos option is enabled, the asmbl_id to which each match_gene_name belongs to MUST exist in the gene position xml
	if($check_gene_pos) {              
	    my $seqID = $proteinID_seqID->{$id_protID->{$m_feat_name}};  
	    if(!exists($valid_asmbl_ids->{$seqID})) { #skip if asmbl_id NOT in gene position xml
		#print STDERR "$seqID to which $m_feat_name belongs to is NOT in gene position xml.  Skipping...\n";
                next;
	    }
	}
	$pexml->addAlignment($q_feat_name, $m_feat_name, $per_sim, $per_id, $pvalue);
    }

}


sub build_id_protID_mapping {
#This function builds a mapping between id to proteinID. 
#The returned structure is a hash ref, where key is id, value is proteinID

    my $rhash = shift;
    my $id_protID=shift;
    my $proteinID_seqID=shift;

    foreach my $seqID (keys %$rhash) {
	foreach my $geneID (keys %{ $rhash->{$seqID} }) {
	    foreach my $transcriptID (keys %{ $rhash->{$seqID}->{$geneID} }) {
		my $cdsID = $rhash->{$seqID}->{$geneID}->{$transcriptID}->{'cdsId'};
		my $proteinID = $rhash->{$seqID}->{$geneID}->{$transcriptID}->{'proteinId'};
		my $transcriptID = $rhash->{$seqID}->{$geneID}->{$transcriptID}->{'transcriptId'};
		my $geneID = $rhash->{$seqID}->{$geneID}->{$transcriptID}->{'geneId'};
		$id_protID->{$cdsID} = $proteinID;
		$id_protID->{$transcriptID} = $proteinID;
		$id_protID->{$geneID} = $proteinID;
		$id_protID->{$proteinID} = $proteinID;
		$proteinID_seqID->{$proteinID} = $seqID;
	    }
	}
    }

}

sub build_id_lookups {

    my $BSML_dir = shift;
    my $id_protID = {};
    my $proteinID_seqID = {};

    my @files = <$BSML_dir/*.bsml>;
    my $new_parser = new BSML::BsmlParserTwig;

    foreach my $bsml_doc (@files) {
	if (-s $bsml_doc) {
	    print STDERR "parsing $bsml_doc\n" if($verbose);
	    my $reader = BsmlCGCReader->new();
	    $new_parser->parse( \$reader, $bsml_doc );
	    my $rhash = $reader->returnAllIdentifiers();
	    build_id_protID_mapping($rhash, $id_protID, $proteinID_seqID); 
	} else {
	    print STDERR "Empty $bsml_doc...skipping\n" if($verbose);
        }
    }

    return ($id_protID, $proteinID_seqID);

}





sub determine_all_asmbl_id {
#given $asmbl_ids and/or $asmbl_file
#return a list of asmbl_ids
#if both arguments are supplied, will consider $asmbl_ids 
#ONLY and IF ONLY $asmbl_file is not defined
    
    my $asmbl_ids  = shift;
    my $asmbl_file = shift;
    my $bsml_dir   = shift;

    my @final_asmbl_ids;
    my (%unique_asmbl_ids, @id_list);

    #check asmbl_id from a flat file
    if(defined($asmbl_file)) {
	@id_list = read_asmbl_file($asmbl_file);
	foreach my $asmbl_id (@id_list) {
	    push(@final_asmbl_ids, $asmbl_id) unless $unique_asmbl_ids{$asmbl_id}++;
	}
	return @final_asmbl_ids;
    } else {
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


sub cmd_check {
#quality check

    if( exists($options{'man'})) {
	pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT});
    }   

    if( exists($options{'help'})) {
	pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT});
    }

    if(!$BSML_dir or !$bsml_btab_dir) {
	pod2usage({-exitval => 2,  -message => "$0: Must specify --bsml_dir and --bsml_btab_dir", -verbose => 1, -output => \*STDERR});
    }
    
    if(!$asmbl_file and !$ASMBL_IDS) {
	pod2usage(-exitval =>2, -verbose => 1, -message => "$0: Must specify either --asmbl_ids OR --asmbl_file.", -output => \*STDERR);
    }

    if(!$match_asmbl_id) {
	pod2usage(-exitval =>2, -verbose => 1, -message => "$0: Must specify --match_asmbl_id", -output => \*STDERR);
    }

    if($asmbl_file and !$ASMBL_IDS) {
	pod2usage(-exitval =>2, -verbose => 1, -message => "$0: Must specify --asmbl_id if --asmbl_file is specified", -output => \*STDERR);
    }
    
    if($check_gene_pos and !$gene_pos_file) {
	pod2usage(-exitval =>2, -verbose => 1, -message => "$0: Must specify --gene_pos_file if --gene_pos_check is specified", -output => \*STDERR);
    }


    #checking BSML repository directory
    if(! -d $BSML_dir) {
	print STDERR "$BSML_dir NOT found.  Aborting...\n";
	exit 5;
    }

    #checking BSML btab  directory
    if(! -d $bsml_btab_dir) {
	print STDERR "$bsml_btab_dir NOT found.  Aborting...\n";
	exit 5;
    }


}

sub print_usage {


    print STDERR "SAMPLE USAGE:  gene_match_bsml.pl -a bsp_3839_assembly -b gbs_799_assembly -d bsml_dir -f allvsall.bsml\n";
    print STDERR "  --asmbl_ids  = assembly id\n";
    print STDERR "  --asmbl_file  = name of the file containing a list of asmbl_ids\n";
    print STDERR "  --match_asmbl_ids = assembly id or 'all' \n";
    print STDERR "  --gene_pos_check  = check to see if asmbl is already in gene_pos_file\n";
    print STDERR "  --gene_pos_file  = gene position xml file \n";
    print STDERR "  --bsml_dir = directory containing the BSML documents\n";
    print STDERR "  --bsml_btab_dir(-f) = dir containing bsml encoding btabs\n";
    print STDERR "  --help = This help message.\n";
    exit 1;

}
