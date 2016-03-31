#!/usr/bin/perl

=head1 NAME

run_humann2.pl - runs humann2 and humann2_renorm_table

=head1 SYNOPSIS

USEAGE: run_humann2.pl
       --input_file=/path/to/input/file.txt
       --metaphlan=/path/to/metaphlan/dir
       --metaphlan_options="option1 option2"
       --bowtie2=/path/to/bowtie2/dir
       --diamond=/path/to/diamond/dir
       --nucleotide_database=/path/to/chocophlan/dir
       --protein_database=/path/to/uniref/dir
       --identity_threshold=40
       --evalue=1.0
       --norm=relab
       --force=no
       --humann_dir=/path/to/humann2/dir
       --output_dir=/path/to/output/dir
       --other_opts="other options"
     [ --log
       --help
     ]

=head1 OPTIONS
B<--input_file, -i>
    Required. input file of type {fastq,fastq.gz,fasta,fasta.gz,sam,bam,blastm8,genetable,biom} 

B<--metaphlan, -m>    
    Optional. Directory containing the MetaPhlAn software. [DEFAULT: $PATH]

B<--metaphlan_options, -l>    
    Optional. Options to be provided to the MetaPhlAn software. [DEFAULT: "-t rel_ab"]

B<--bowtie2, -b>    
    Optional. Directory containing the bowtie2 executable. [DEFAULT: $PATH]

B<--diamond, -d>    
    Optional. Directory containing the diamond executable. [DEFAULT: $PATH]

B<--nucleotide_database, -n>    
    Required.  nucleotide database to download

B<--protein_database, -p>
    Required. protein database to download

B<--force, -f>
    Optional. force-download database?

B<--identity_threshold, -t>
    Optional. identity threshold for alignments. Default 40.0

B<--evalue, -e>
    Optional.  the evalue threshold to use with the translated search. Default 1.0.

B<--norm, -r>
    Required.  Normalization scheme: copies per million [cpm], relative abundance [relab]; default=[cpm]. Default relab.

B<--humann_dir, -h>
    Optional. Directory where humann2 is installed.

B<--output_dir, -o>
    Required. Directory for HUMAnN output.

B<--other_opts, -z>
    Optional. Other humann2 options not specified above.

B<--log,-l>
    Log file

B<--help, h>
    This help message

=head1 DESCRIPTION

Runs HUMAnN2 pipeline.

=head1 CONTACT

    Kemi Ifeonu
    kifeonu@som.umaryland.edu

=cut


use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Basename;

my $input_file;
my $bowtie2;
my $diamond;
my $nucleotide_database;
my $protein_database;
my $metaphlan;
my $metaphlan_options;
my $identity_threshold;
my $evalue;
my $norm;
my $force;
my $humann_dir;
my $output_dir;
my $other_opts;
my @dir =();

my %options = ();
my $results = GetOptions (\%options, 
			  'metaphlan|m:s',
                          'metaphlan_options|l:s',
			  'input_file|i=s',
			  'other_opts|z:s',
                          'bowtie2|b:s',
                          'diamond|d:s',
                          'nucleotide_database|n=s',
                          'protein_database|p=s',
			  'force|f:s',
                          'evalue|e:f',
                          'identity_threshold|t:f',
			  'norm|r=s',
			  'humann_dir|d:s',
			  'output_dir|o=s',
                          'log|l=s',
                          'help|h') || pod2usage();

&check_options(\%options);

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">$options{log}") || die "can't create log file: $!";
}

my $str = "";
my $exit_code = 0;

#Download database if it doesn't exist already or if user selects force-download
if ((-d "/mnt/staging/data/chocophlan")&&($force eq "no")) {
	print "\nNucleotide database already exists: /mnt/staging/data/chocophlan\n";
}else{ 
	if (-d "/mnt/staging/data/chocophlan"){
	   system ("rm -r /mnt/staging/data/chocophlan");
	}
	$str = "humann2_databases --download chocophlan $nucleotide_database /mnt/staging/data";
	$exit_code = system($str);
	if($exit_code!=0){
  		print "Command \"$str\" failed with an exit code of $exit_code.\n";
  		exit($exit_code >> 8);
	}else{
		print "\n Chocophlan downloaded.\n";
	}
}
if ((-d "/mnt/staging/data/uniref")&&($force eq "no")) {
        print "\nProtein database already exists: /mnt/staging/data/uniref \n";
}else{          
	if (-d "/mnt/staging/data/uniref"){
           system ("rm -r /mnt/staging/data/uniref");
        }
	$str = "humann2_databases --download uniref $protein_database /mnt/staging/data";
        $exit_code = system ($str);
        if($exit_code!=0){
                print "Command \"$str\" failed with an exit code of $exit_code.\n";
                exit($exit_code >> 8);
        }else{
		print "\n Uniref downloaded.\n";
	}
}

$str = "humann2 --input $input_file --output $output_dir --protein-database=/mnt/staging/data/uniref --nucleotide-database=/mnt/staging/data/chocophlan";
if ($metaphlan ne "") {$str .= " --metaphlan $metaphlan";}
if ($bowtie2 ne "") {$str .= " --bowtie2 $bowtie2";}
if ($diamond ne "") {$str .= " --diamond $diamond";}
if ($metaphlan_options ne "") {$str .= " --metaphlan-options \"$metaphlan_options\"";} 
if ($evalue ne "") {$str .= " --evalue $evalue";}
if ($identity_threshold ne "") {$str .= " --identity-threshold $identity_threshold";}
if ($other_opts ne "") {$str .= " $other_opts";}

#run humann2
print "\n\nRunning...\n$str\n\n";    
$exit_code = system ("$str");
if($exit_code!=0){
	print "Command \"$str\" failed with an exit code of $exit_code.\n";
        exit($exit_code >> 8);
}else{
        print "\n humann2 ran sucessfully.\n";
}


#normalize the abundance output files
my ($sample,$dir,$ext) = fileparse($input_file,'\..*');
$str = "humann2_renorm_table --input ".$output_dir."\/".$sample."_genefamilies.tsv --output ".$output_dir."\/".$sample."_genefamilies_norm.tsv -u ".$norm;
print "\n\nNormalizing output files...\n$str\n\n";
$exit_code = system ("$str");
if($exit_code!=0){
        print "Command \"$str\" failed with an exit code of $exit_code.\n";
        exit($exit_code >> 8);
}else{
        print "\n Genefamilies sucessfully normalized.\n";
}

$str = "humann2_renorm_table --input ".$output_dir."\/".$sample."_pathabundance.tsv --output ".$output_dir."\/".$sample."_pathabundance_norm.tsv -u ".$norm;
print "\n\nNormalizing output files...\n$str\n\n";
$exit_code = system ("$str");
if($exit_code!=0){
        print "Command \"$str\" failed with an exit code of $exit_code.\n";
        exit($exit_code >> 8);
}else{
        print "\n Pathabundance sucessfully normalized.\n";
}



exit(0);


sub _log {
    my $msg = shift;
    print $logfh "$msg\n" if $logfh;
}

sub check_options {
    my ($opts) = @_;

    if( $opts->{'help'} ) {
        pod2usage();
        exit(0);
    }


   if($opts->{'input_file'} && $opts->{'input_file'} ne "") {
        &_die("input_file [$opts->{'input_file'}] does not exist") unless( -e $opts->{'input_file'});
        $input_file = $opts->{'input_file'}; 
   }

   my @reqs = qw( nucleotide_database protein_database output_dir);
    foreach my $req ( @reqs ) {
        die("Option $req is required") unless( exists( $opts->{$req} ) );
    }

    if (exists( $opts-> {'metaphlan_options'})) {
	$metaphlan_options = $opts-> {'metaphlan_options'};
    }else{
	$metaphlan_options = "";
    }	
    
    if (exists( $opts-> {'bowtie2'})) {
	$bowtie2 = $opts-> {'bowtie2'};
    }else{
        $bowtie2 = "";
    }
    
    if (exists( $opts-> {'norm'})) {
	$norm = $opts-> {'norm'};
    }else{
	$norm = "cpm";
    }    

    if (exists( $opts-> {'force'})) {
        $force = $opts-> {'force'};
    }else{
        $force = "no";
    }

    if (exists($opts-> {'diamond'})) {
	$diamond = $opts-> {'diamond'};
    }else{
    	$diamond = "";
    }
    $nucleotide_database = $opts-> {'nucleotide_database'};
    $protein_database = $opts-> {'protein_database'};
    if (exists( $opts-> {'metaphlan'})) {
        $metaphlan = $opts->{'metaphlan'};
    }else{
        $metaphlan = "";
    }

    if (exists( $opts-> {'identity_threshold'})) {
    	$identity_threshold = $opts->{'identity_threshold'};
    }else{
	$identity_threshold = "";
    }    
    
    if (exists( $opts-> {'metaphlan_options'})) {
        $metaphlan_options = $opts->{'metaphlan_options'};
    }else{
        $metaphlan_options = "";
    }


    if (exists( $opts-> {'evalue'})) {
    	$evalue = $opts->{'evalue'};
    }else{
	$evalue = "";
    }
   # $humann_dir = $opts->{'humann_dir'};
    $output_dir = $opts->{'output_dir'};
    
    if (exists( $opts-> {'other_opts'})) {
	$other_opts = $opts->{'other_opts'};
    }else{
	$other_opts = "";
    }


}
