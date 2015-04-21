#!/usr/bin/perl

=head1 NAME

download_all.pl - Download all annotation and sequences files related to the databases specified from chado db 

=head1 SYNOPSIS

 USAGE: download_all.pl
       --host=jabba
       --username=user
       --password=password
       --input_file=/path/to/input
       --output_dir=/path/to/output/dir
       [--locus_db=TIGR_moore
       --translation_table=11
       --feature_type=assembly]
       --help


=head1 OPTIONS

B<--host,-h>
    Server where the database is located

B<--username,-u>
    Username for database

B<--password,-p>
    Password for database

B<--input_file,-i>
	Input file

B<--output_dir, -o>
	Output dir location

B<--perl_scripts_dir, -d>
	Optional. The path of the chado2flatfile, chado2aengine_dumper, and other Perl scripts called.  All called scripts must be located in the same directory. (Default = .)

B<--locus_db, -l>
    Optional. The name of the db of the locus identifiers to be used in the output (default = Tigr_moore)

B<--locus_db_version, -L>
	Optional.  The version type of locus db tag the locus_db identifier is (default = public_locus)

B<--translation_table>
    Optional.  Numeric value for translations table (default = 11)

B<--feature_type,-f>
    Optional. Type of feature for which you want to extract sequences, such as assembly or
    polypeptide.  Corresponds to cvterm.name. (default = assembly )

B<--help,-h>
    Print this message

=head1  DESCRIPTION

    Reads information from a chado database and creates
    	a pep file,
    	seq file,
    	coords file,
    	tab delimited annotation file,
	tbl annotation file,
    	gff3,
    	gbk,
    	whole genome
    for all the organism databases specified in the input file.

=head1  INPUT

    The required input are name of host, username, password, output directory and input file.

    Input file has organism/database related information its a tab delimited file with 4 columns species, common name, abbreviation and database
    e.g:
    Escherichia	Escherichia coli 1827-70    gec1827_70
	Escherichia Escherichia coli 2362-75    gec2362
	Escherichia Escherichia coli 2534-86    gec2534_86
	Escherichia Escherichia coli 3030-1 gec3030
	Escherichia Escherichia coli 3431   gec3431
	Escherichia Escherichia coli DEC10A gecDEC10A
	Escherichia Escherichia coli DEC10B gecDEC10B
	Escherichia Escherichia coli DEC10C gecDEC10C
	Escherichia Escherichia coli DEC10D gecDEC10D
	Escherichia Escherichia coli DEC10E gecDEC10E
	Escherichia Escherichia coli DEC10F gecDEC10F
	Escherichia Escherichia coli DEC13D gecDEC13D
	Escherichia Escherichia coli DEC1A  gecDEC1A
	Escherichia Escherichia coli DEC1C  gecDEC1C

=head1 OUTPUT

    All output files are sorted first by molecule located on and then by featureloc.fmax

    pep file:
    Multi fasta file of peptides. Headers have identifiers (feature.uniquename, unless --locus_ids option is used)
    followed by the common name (featureprop.value for 'gene_product_name' cvterm) then followed by the
    molecule uniquenae on which the feature is located.

    ex:
    >hik_1522	glyceraldehyde-3-phosphate dehydrogenase, type I hik.assembly.1
    MAIKIGINGFGRIGRIVFRAAQHRDDIEVVGINDLIDVEYMAYMLKYDSTHGRFDGTVEV
    KDGNLVVNGKTIRVTAERDPANLNWGAIGVDIAVEATGLFLTDETARKHITAGAKKVVLT
    GPSKDATPMFVRGVNFNAYAGQDIVSNASCTTNCLAPLARVVHETFGIKDGLMTTVHATT
    ATQKTVDGPSAKDWRGGRGASQNIIPSSTGAAKAVGKVLPALNGKLTGMAFRVPTPNVSV
    VDLTVNLEKPASYDAIKQAIKDAAEGKTFNGELKGVLGYTEDAVVSTDFNGCALTSVFDA
    DAGIALTDSFVKLVSWYDNETGYSNKVLDLVAHIYNYKG
    >hik_1013	hypothetical protein
    LKKLTALRSVFYYLEFNLRSKNEQGITKNNLKQSIFF
    >hik_593	AMP-binding enzyme family protein
    MNLDLHFVHRIQQQAKTRTNMTALRYKEHGLWRDISWKNFQEQLNQLSRALLAHNIDVQD
    KIAIFAHNMERWTIVDIATLQIRAITVPIYATNTAQQAEFILNHADVKILFVGDQEQYDQ
    TLEIAHHCPKLQKIVAMKSTIQLQQDPLSCTWESFIKTGSNAQQDELTQRLNQKQLSDLF
    TIIYTSGTTGEPKGVMLDYANLAHQLETHDLSLNVTDQDISLSFLPFSHIFERAWAAYIL
    HRGAILCYLEDTNQVRSALTEIRPTLMCAVPRFYEKIYAAVLDKVQKAPKLRQIMFHWAI
    SVGQKYFDLRANNKAIPFLLKKQFALADKLVLSKLRQLLGGRIKMMPCGGAKLEPAIGLF
    FHAIGINIKLGYGMTETTATVSCWHDFQFNPNSIGTLMPKAEVKIGENNEILVRGGMVMK
    GYYKKPEETAQAFTEDGFLKTGDAGEFDEQGNLFITDRIKELMKTSNGKYIAPQYIESKI
    GKDKFIEQIAIIADAKKYVSALIVPCFDSLEEYAKQLNIKYHDRLELLKNSDILKMFEHR
    INAVQKELAHFEQVKKFTLLSQAFSIKLGEITPTLKLRRKVILERYRKQIEAMYHSQEA
    ...

    seq file:
    Same as pep file but with nucleotide sequences

    coords file:
    tab separated file with identifiers (uniquename unless --locus_ids option is used)
    followed by molecule name, end5 and end3.

    ex:
    hik_1522	hik.assembly.1   1	1021
    hik_1013	hik.assembly.1   1029	1143

    ...

    tab delimited annotaiton file:
    cgsp_3829       diguanylate cyclase (GGDEF) domain protein                      GO:0009975,GO:0009966   264,710
    cgsp_3830       AFG1-like ATPase family protein                         157
    cgsp_3831       transposase, Mutator family protein                             154
    cgsp_3832       hypothetical protein                            856
    cgsp_3833       alkaline phosphatase                            703
    cgsp_3834       ATP-dependent Clp protease, ATP-binding subunit ClpX            clpX    GO:0008462,GO:0051082,GO:0006510,GO:0009368,GO:0005524  95,138
    cgsp_3835       hypothetical protein                            856
    cgsp_3836       conserved hypothetical protein                          156
    cgsp_3837       conserved hypothetical protein                          156
    cgsp_3838       HTH-type transcriptional regulator prtR (Pyocin repressor protein)                              129
    cgsp_3839       RNA polymerase sigma-H factor (Sigma-30)                                703
    cgsp_3840       pyrimidine-specific ribonucleoside hydrolase rihB (Cytidine/uridine-specific hydrolase)                         703
    cgsp_3841       tryptophan synthase, alpha subunit      4.2.1.20        trpA    GO:0004834,GO:0000162   70
    cgsp_3842       hypothetical protein                            856
    ...
    ...
    ...


=head1  CONTACT


=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Data::Dumper;
use Time::HiRes;

my %options;

my $results = GetOptions (\%options,
                          'host|h=s',
                          'username|u=s',
                          'password|p=s',
                          'input_file|i=s',
                          'output_dir|o=s',
                          'locus_db|l=s',
						  'perl_scripts_dir|d=s',
						  'locus_db_version|L=s',
                          'translation_table|r=i',
                          'feature_type|f=s',
                          'help'
                          );

&check_options( \%options );

my $cmd = '';
my $start_time = [Time::HiRes::gettimeofday()];

open (FHD, "<", $options{input_file}) or die "Could not open input file $options{input_file}\n";

while (<FHD>){
        chomp $_;
        next if ($_ =~ /^#/);
	next if ($_ =~ /^\s*$/);
	my @info = split(/\t/,$_);

	my $organism = $info[1];
	my $abbr = $info[2];
	my $db = $info[3];
	$organism =~ s/^\s+//;
	$organism =~ s/\s+$//;
	$abbr =~ s/^\s+//;
	$abbr =~ s/\s+$//;
        $db =~ s/^\s+//;
	$db =~ s/\s+$//;
	print "Processing database : $db\n";
	my $out_dir = $options{output_dir}."/".$abbr;
	my $file_name = $out_dir."/".$abbr;
	system("mkdir -p $out_dir");

	print "Processing organism : $organism [$abbr]\n";

	print "\tProcessing and writing pep, seq and coords files.\n";
	#create pep, seq and coord file
	$cmd = "$options{perl_scripts_dir}/chado2flatfile.pl --database=$db --host=$options{host}".
		   " --username=$options{username} --password=$options{password} --database_type=mysql".
		   " --locus_ids=$options{locus_db} --organism='$organism'".
		   " -e $file_name.nucleotide.fsa -f $file_name.polypeptide.fsa -c $file_name.coord.txt";

	print "\t\t".$cmd."\n";
	system($cmd);
	my $diff = Time::HiRes::tv_interval($start_time);
	$start_time = [Time::HiRes::gettimeofday()];
	print "\tCompleted in $diff seconds\n\n";

	print "\tProcessing and writing tbl, gbk and gff files.\n";
	#create annotation tbl file
	$cmd = "$options{perl_scripts_dir}/chado_aengine_dumper.pl --database=$db --server=$options{host}".
		   " --user=$options{username} --password=$options{password} --database_type=mysql".
		   " --output_directory=$out_dir --format='tbl' --locus_db=$options{locus_db}".
		   " --organism='$organism' --translation_table=$options{translation_table}".
		   " --output_filename=$abbr.annotation.tbl --locus_db_version=$options{locus_db_version}";

	print "\t\t".$cmd."\n";
	system($cmd);
	$diff = Time::HiRes::tv_interval($start_time);
	$start_time = [Time::HiRes::gettimeofday()];
	print "\tCompleted in $diff seconds\n";

	#create annotation gbk file
	$cmd = "$options{perl_scripts_dir}/chado_aengine_dumper.pl --database=$db --server=$options{host}".
		   " --user=$options{username} --password=$options{password} --database_type=mysql".
		   " --output_directory=$out_dir --format='gbk' --locus_db=$options{locus_db}".
		   " --organism='$organism' --translation_table=$options{translation_table}".
		   " --output_filename=$abbr.annotation.gbk --locus_db_version=$options{locus_db_version}";

	print "\t\t".$cmd."\n";
	system($cmd);
	$diff = Time::HiRes::tv_interval($start_time);
	$start_time = [Time::HiRes::gettimeofday()];
	print "\tCompleted in $diff seconds\n";

	#create annotation gff file
	$cmd = "$options{perl_scripts_dir}/chado_aengine_dumper.pl --database=$db --server=$options{host}".
		   " --user=$options{username} --password=$options{password} --database_type=mysql".
		   " --output_directory=$out_dir --format='gff' --locus_db=$options{locus_db}".
		   " --organism='$organism' --translation_table=$options{translation_table}".
		   " --output_filename=$abbr.annotation.gff --locus_db_version=$options{locus_db_version}";

	print "\t\t".$cmd."\n";
	system($cmd);
	$diff = Time::HiRes::tv_interval($start_time);
	$start_time = [Time::HiRes::gettimeofday()];
	print "\tCompleted in $diff seconds\n\n";
	
	print "\tProcessing and writing whole genome file.\n";
	#extract sequence file
	$cmd = "$options{perl_scripts_dir}/extract_sequences_by_type.pl --database=$db --host=$options{host}".
		   " --username=$options{username} --password=$options{password} --database_type=mysql".
		   " --feature_type=$options{feature_type}".
		   " --output_file=$file_name.whole_genome.txt --organism='$organism'";

	print "\t\t".$cmd."\n";
	system($cmd);
	$diff = Time::HiRes::tv_interval($start_time);
	$start_time = [Time::HiRes::gettimeofday()];
	print "\tCompleted in $diff seconds\n\n";

	print "\tProcessing and writing annotation tab file.\n";
	#create tab delimieted annotation file
	$cmd = "$options{perl_scripts_dir}/annotation_tab_file.pl --database=$db --host=$options{host}".
		   " --username=$options{username} --password=$options{password} --database_type=mysql".
		   " --output_file=$file_name.annotation.txt --organism='$organism'";

	print "\t\t".$cmd."\n";
	system($cmd);
	$diff = Time::HiRes::tv_interval($start_time);
	print "\tCompleted in $diff seconds\n\n";
}

close FHD;
exit();

#########################################################################
sub check_options {
    my $options = shift;

    ## make sure required arguments were passed
    my @required = qw(host username password input_file output_dir);
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }

	$$options{perl_scripts_dir} = '.' unless defined $$options{perl_scripts_dir};
    $$options{locus_db} = 'TIGR_moore' unless defined $$options{locus_db};
    $$options{locus_db_version} = 'public_locus' unless defined $$options{locus_db_version};
	$$options{feature_type} = 'assembly' unless defined $$options{feature_type};
    $$options{translation_table} = 11 unless defined $$options{translation_table};
}
