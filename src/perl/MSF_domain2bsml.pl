#!/usr/local/bin/perl


=head1  NAME 

MSF_domain2bsml.pl  - convert MSF format alignment files into BSML documents

=head1 SYNOPSIS

USAGE:  MSF_domain2bsml.pl -a ali_dir -o msf.bsml

=head1 OPTIONS

=over 4

=item *

B<--ali_dir,-a>   [REQUIRED] Dir containing MSF format alignment files that ends in *.msf

=item *

B<--output,-o> [REQUIRED] output BSML file

=item *

B<--suffix,-s> suffix of the multiple sequence alignment file. Default(msf)

=item *

B<--program,-p> projects this analyisi is for. Default(clustal)

=item *

B<--help,-h> This help message

=back

=head1   DESCRIPTION

MSF_domain2bsml.pl is designed to convert multiple sequence alignments in MSF format
into BSML documents.  More specifically this script is equippd to parse the multiple 
sequence alignments of domains.  This requires the capture of the start and end 
coordinates of each "subsequence" in the alignment.  As an option, the user can
specify the suffix of the alignment files that the script should parse in the 
directory specified by --ali_dir.  The default suffix is "msf" i.e. family.msf.
In addition, the user can specify what program this script is for.  For example if the
multiple sequence alignment represents COGS, the user can specify that
with the --program flag.  This will be reflected in the analysis object
of the BSML document.

Samples:

1. Convert all MFS formatted files that end in "aln" in ali_dir
   MSF_domain2bsml.pl -a ali_dir -o alignment.bsml -s aln


NOTE:  

Calling the script name with NO flags/options or --help will display the syntax requirement.


=cut


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
BEGIN {
    require '/usr/local/devel/ANNOTATION/cas/lib/site_perl/5.8.5/BSML/BsmlBuilder.pm';
    import BSML::BsmlBuilder;
}
use File::Basename;
use Pod::Usage;


my %options = ();
my $results = GetOptions (\%options, 'ali_dir|a=s', 'output|o=s', 'program|p=s', 
                                     'verbose|v', 'suffix|s=s', 'help|h', 'man') || pod2usage();

###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $output    = $options{'output'};
my $ali_dir   = $options{'ali_dir'};
$ali_dir =~ s/\/+$// if($ali_dir);
my $msf_suffix = $options{'suffix'} || "msf";  #default suffix pattern of alignment files
my $program   = $options{'program'} || "clustal";

&cmd_check();
###-------------------------------------------------------###

my $builder = new BSML::BsmlBuilder;

# add an analysis object to the document. 
my $filename = basename($output);
$builder->createAndAddAnalysis( 'sourcename' => $filename,
				'programversion' => '1.0',
				'program' => $program,
				'bsml_link_url' => 'BsmlTables',
				'bsml_link_relation' => 'MULTIPLE_ALIGNMENTS' );

my $file_found = 0;
my $id=0;
opendir(DIR, $ali_dir) or die "Unable to access $ali_dir due to $!";
while( my $file = readdir(DIR)) {
    next if ($file =~ /^\.{1,2}$/);  #skip  "." ,  ".."
    if($file =~ /^(.+?)(\.msf\.dtab)?\.$msf_suffix$/) {
	$file_found++;
	my $accession = undef;
	my $domain_id = "domain_";
	if($file =~ /^(PF\d{5})\.$msf_suffix$/ || $file =~ /^(TIGR\d{5})\.$msf_suffix$/) { 
	    $accession = $1;
	    $domain_id .= $accession;
	} else {
	    #my $fam_name = $1;
	    $id++;
	    $domain_id .= $id
	}
	my $MSF_alignments = process_MSF_file("$ali_dir/$file");
	if(keys %$MSF_alignments < 1) {  #skip empty msf files
	    #this checks to see if the previous domain_id is from $id++ or $accession
            #and will subtract 1 from $id if previous domain_id is from $id++
	    if(!defined($accession)) { 
		$id--;
            }
	    next;
	}

	my $table = $builder->createAndAddMultipleAlignmentTable('molecule-type' => $MSF_alignments->{'mol_type'});
	$table->addattr( 'id', $domain_id );

	my $summary = $builder->createAndAddAlignmentSummary( 'multipleAlignmentTable' => $table,
							      'seq-type' =>  $MSF_alignments->{'mol_type'},
							      'seq-format' => 'MSF' 
							      );
	my $aln = $builder->createAndAddSequenceAlignment( 'multipleAlignmentTable' => $table );
	my $seqnum=0;
	my $sequences_tag;
	foreach my $seq (keys %{ $MSF_alignments->{'polypeptides'} }) {
	    $seqnum++;
	    my ($pro_acc, $coords) = split(/\//, $seq);
	    $pro_acc =~ s/\|\|$//;                         #remove terminal '||' if exists
	    my $alignment = join ('', @{ $MSF_alignments->{'polypeptides'}->{$seq}->{'alignment'} });
	    my $start  = $MSF_alignments->{'polypeptides'}->{$seq}->{'end5'}; 
	    my $end    = $MSF_alignments->{'polypeptides'}->{$seq}->{'end3'};
	    my $domain_length = $end - $start + 1;
	    #IMPORTANT!!!!
	    #In order to ensure that each seq in a multiple sequence alignment is truly
            #unique, the seq-name and name will be in the form "polypeptide_accession:seqnum"
            #i.e. (ana1.10005.m00234_polypeptide:1). 

	    $builder->createAndAddAlignedSequence( 'alignmentSummary' => $summary,
						   'seqnum' => $seqnum,
						   'length' => $domain_length,
						   'name'   => "$pro_acc:$seqnum", #<--see notes above
                                                   'start'  => $start
                                                 );
	    $builder->createAndAddSequenceData( 'sequenceAlignment' => $aln,
						'seq-name' => "$pro_acc:$seqnum",  #<--see notes above
						'seq-data' => $alignment
                                              );
	    $sequences_tag .= "$seqnum:";
	}
	$aln->addattr( 'sequences', $sequences_tag );

    }

}
if($file_found) {
    $builder->write( $output );
    chmod 0777, $output;
} else {
    print STDERR "No files that ends in \"$msf_suffix\" are found in $ali_dir\n";
}





sub process_MSF_file {

    my $file = shift;


    my $MSF_alignments ={};
    open(MSF, "$file") or die "Unable to open $file due to $!";
    my $line;
    my $msf_type;
    while(defined($line = <MSF>) and $line !~ /^\/\//) {
	if( $line =~ /MSF:\s*([\S]+)\s*Type:\s*([\S]+)\s*Check/) {
	    my $msf_length = $1;
	    return undef if($msf_length == 0);  #abort if align_len = 0

	    if($2 eq 'P') {
		$msf_type = 'polypeptide';
	    }elsif($2 eq 'N') {
		$msf_type = 'nucleotide';
	    }else {
		$msf_type = 'polypeptide';
	    }
	    $MSF_alignments->{'mol_type'} = $msf_type;
	}

	if($line =~ /Name:\s*([\S]+)\s*Len:\s*([\S]+)\s*Check:\s*([\S]+)\s*Weight:\s*([\S]+)/) {
	    my $name    = $1;
	    my $name2   = $name;
	    my $ali_len = $2;
	    my $check   = $3;
	    my $weight  = $4;
	    my ($polypeptide, $coord) = split(/\//, $name);
	    my ($start5, $end3) = split("-", $coord);
	    
	    $MSF_alignments->{'polypeptides'}->{$name}->{'length'} = $ali_len;
	    $MSF_alignments->{'polypeptides'}->{$name}->{'check'}  = $check;
	    $MSF_alignments->{'polypeptides'}->{$name}->{'weight'} = $weight;
	    $MSF_alignments->{'polypeptides'}->{$name}->{'end5'}   = $start5;
	    $MSF_alignments->{'polypeptides'}->{$name}->{'end3'}   = $end3;
	    $MSF_alignments->{'polypeptides'}->{$name}->{'alignment'} = [];
	}
    }

    my $replacements;
    my $spaces;
    while($line = <MSF>) {
	if($line =~ /^([\S]+)/) {
	    my $name = $1;
	    my $name2  = $name;
	    if($name =~ /(\/.+)/) {
		$replacements = $1;
		$spaces = " " x length($replacements);
		$line =~ s/$replacements/$spaces/;
		$name =~ s/$replacements//g;
	    }
	    if(exists($MSF_alignments->{'polypeptides'}->{$name2})) {
		push( @{ $MSF_alignments->{'polypeptides'}->{$name2}->{'alignment'} }, $line );
            } else {
		print STDERR "ERROR, $name is not valid polypeptide name\n";
            }
	}
    }

    return $MSF_alignments;

}
	
	
sub cmd_check {
#quality check

    if( exists($options{'man'})) {
	pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT});
    }   

    if( exists($options{'help'})) {
	pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT});
    }

    if(!$output or !$ali_dir) {
	pod2usage({-exitval => 2,  -message => "$0: All the required options are not specified", -verbose => 1, -output => \*STDERR});
    }
    
    if(! -d $ali_dir) {
	print STDERR "$ali_dir directory NOT found.  Aborting...\n";
	exit 5;
    }

}	    
