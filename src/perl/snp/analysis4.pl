#! /usr/local/bin/perl -w
#------------------------------------------------------------
#
#
#
#
#
#------------------------------------------------------------
=head1 NAME

euktigr2chado.pl - Migrates Euk legacy datasets to Chado schema

=head1 SYNOPSIS

USAGE:  euktigr2chado.pl -u username -p password -d source_database -t target_database [-a] [-c] [-f]

=head1 OPTIONS

=over 8

=item B<--username,-u>
    
    Database username

=item B<--password,-p>
    
    Database password

=item B<--source_database,-d>
    
    Source database name

=item B<--target_database,-t>
    
    Destination database name

=item B<-a>
    
    Optional  - used to specify a list of assembly ids

=item B<-c>

    Optional -  "clears" deletes all records from the Chado Sequence Module tables

=item B<-f>

    Optional. User can specify whether to produce .out files or to migrate directly into chado database

=item B<--help,-h>

    Print this help


=back

=head1 DESCRIPTION

    euktigr2chado.pl - Migrates Euk legacy datasets to Chado schema

=cut

no strict "refs";
use lib "/home/jravel/lib/";
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Data::Dumper;
use Digest::MD5 qw(md5);
use Log::Log4perl qw(get_logger);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Config::IniFiles;
use Benchmark;

$| = 1;

my ($help, $filename, $contig, $logfile, $coverage, $config_file, $remove_files);

my $result = GetOptions (
			 'logfile|l'         => \$logfile,
			 'help|h'            => \$help,
			 'f=s'               => \$filename,			
			 'c=s'               => \$contig,
			 'cov=s'             => \$coverage,
			 'config_file|C'     => \$config_file,
			 'remove_files|r'    => \$remove_files
			 );

if (!$filename or !$contig){
    &print_usage;
}
if ($help){
    pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT});
}



my $screen_threshold = 'WARN';
if ($verbose){
    $screen_threshold = 'INFO';
}
Log::Log4perl->init(
		    \ qq{
			log4perl.logger                                    = ERROR, A1, Screen
			log4perl.appender.A1                               = Log::Dispatch::File
		        #log4perl.appender.A1.filename                     = euktigr2chado.log
		        log4perl.appender.A1.filename                      = $logfile
		        log4perl.appender.A1.mode                          = write
		        log4perl.appender.A1.Threshold                     = ERROR
		        log4perl.appender.A1.layout                        = Log::Log4perl::Layout::PatternLayout
		        log4perl.appender.A1.layout.ConversionPattern      = %d %p> %F{1}:%L %M - %m%n 
		        log4perl.appender.Screen                           = Log::Dispatch::Screen
		        log4perl.appender.Screen.layout                    = Log::Log4perl::Layout::SimpleLayout
                        #log4perl.appender.Screen.layout.ConversionPattern = %d %p> %F{1}:%L %M - %m%n 
		        log4perl.appender.Screen.Threshold                 = WARN
		        Log::Log4perl::SimpleLayout
		    }
		    );

my $logger = get_logger("snp");
$logger->info("Executing $0");

my $config_hash = &get_config_hash(\$config_file);
if (!defined($config_hash)){
    $logger->logdie("config_hash was not defined");
}

#
# daddy4{} iterates through all files prefixed with "ID" in the <project>_asmbl_seq.dir directory and calls quality_SNP8
# e.g. anthrax_test_asmbl_seq.dir
#
&daddy4(
	filename       => \$filename,
	contig_file    => \$contig,
	config_hash    => $config_hash
	);

#
# quality_SNP8{} 
#


#
# sort_snp_file{} produces output for downstream analysis: post_analysis.pl
# We could choose to omit executing this sort operation
#
&sort_snp_file(
		snp_header_file      => $config_hash->{'snp_header_file'},
		sort_snp_header_file => $config_hash->{'sort_snp_header_file'},
		);

#
# coverage_stat{}: Comments forthcoming
#
&coverage_stat(
	       snp_filename => $config_hash->{'snp_header_file'},
	       out_filename => $coverage_out,
	       coverage     => $contig_hash->{'coverage'}
	       );

#
# indel{}: Comments forthcoming
#
&indel(
       filename     => \$filename,
       contig       => \$contig,
       config_hash  => $config_hash,
       );
#
# indel_info{}:  Comments forthcoming
#
&indel_info(
	    filename     => \$filename,
	    contig       => \$contig,
	    config_hash  => $config_hash,
	    );

#system ("~/src/SNP/Daddy4.pl -F $filename -C $contig");
#system ("sort -k13nr SNP_Header.txt > SNP_sort_Header.txt");
#system ("~/src/SNP/coverage_stat.pl -F SNP_Header.txt -O SNP_cov -C 20");
#system( "~/src/SNP/indel.pl -T $filename -C $contig > indel.txt");
#system( "~/src/SNP/indel_info.pl -T $filename -C $contig");

#
# remove_files{} optionally can delete temporary process files from the system
# uses Perl "unlink" command
#

#if (!$remove_files){
#    &remove_ins_files(
#		      file_list => $config_hash->{'delete_file_list'}
#		      );
#}

print ("End of program $0\n");
print ("Verify $logfile\n");

exit(0);

#--------------------------------------------------------------------------------------------------------------------
#                 END OF MAIN  -- SUBROUTINES FOLLOW
#
#--------------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------
# get_config_hash()
#
#
#----------------------------------------------------------
sub get_config_hash {

    my $logger = get_logger("snp");
    $logger->info("Entered get_config_hash");
    

    my $file = shift;
    if (!defined($file)){
	$logger->logdie("file was not defined");
    }

    my $file_contents = &get_contents($file);
    if (!defined($file_contents)){
	$logger->logdie("file_contents was not defined");
    }

    my %hash;
    foreach my $line (@$file_contents){
	if (!defined($line)){
	    $logger->logdie("line was not defined");
	}
	my ($key, $val) = split([\s]+,$line);
	if (!defined($key)){
	    $logger->logdie("key was not defined for line:$line");
	}
	if (!defined($val)){
	    $logger->logdie("val was not defined for line:$line");
	}
	$hash{$key} = $val;
    }

    return \%hash;

}#end sub get_config_hash()



#----------------------------------------------------------
# coverage_stat()
#
#
#----------------------------------------------------------
sub coverage_stat {

    my logger = get_logger("snp");
    $logger->info("Entered coverage_stat");

    my (%parameter) = @_;
    my $parameter_hash = \%parameter;
    


    my ($filename, $outputfile, $highcoverage);

    #----------------------------------------
    # Extract arguments from parameter hash
    #
    #----------------------------------------
    if (exists $parameter_hash->{'snp_filename'}){
	$filename = $parameter_hash->{'snp_filename'};
    }
    if (exists $parameter_hash->{'out_filename'}){
	$outputfile = $parameter_hash->{'out_filename'};
    }
    if (exists $parameter_hash->{'coverage'}){
	$highcoverage = $parameter_hash->{'coverage'};
    }

    #----------------------------------------
    # Verify whether arguments were defined
    #
    #----------------------------------------
    if (!defined($filename)){
	$logger->logdie("filename was not defined");
    }
    if (!defined($outpufile)){
	$logger->logdie("outputfile was not defined");
    }
    if (!defined($highcoverage)){
	$logger->logdie("highcoverage was not defined");
    }

    my (@filename, $i, $data, @list, @reject, $help);
    my $info = "Coverage count for $filename";

    push(@filename, $info);


    for ($i=0;$i<$highcoverage+1;++$i){

	open(READ, "grep 'COV: $i ' $filename | wc |");   

	my @temp = <READ>;
	close (READ);

	foreach my $line (@temp) {
	    my @temp2 = split(" ", $line);
	    
	    $data = $i."X coverage : $temp2[0]";
	    push(@filename, $data);
	}
    }



    foreach my $line2 (@filename) { 
	print  $line2, "\n";
    }

    my $file = &get_contents(\$filename);
    if (!defined($file)){
	$logger->logdie("file was not defined");
    }


    foreach my $snp (@$file) {

	my @temp3 = split(" ", $snp);

	if ($temp3[12] > 2 and $temp3[14] > 40) {
	    
	    push(@list, $snp);
	} else {
	    push(@reject, $snp);
	}
	
    }
    
    $outputfile = $outputfile . ".quality";

    unless (open(QUALITY, ">$outputfile")) {
	print STDERR "Can't open file!!\n\n";
	exit;
    }
    
    foreach my $sub (@list) { 
	print QUALITY $sub;
    }
    
    close (QUALITY);

}#end

#----------------------------------------------------------
# indel()
#
#
#----------------------------------------------------------
sub indel{

    my logger = get_logger("snp");
    $logger->info("Entered indel");

    my (%parameter) = @_;
    my $parameter_hash = \%parameter;
    

    my ($tag, $genome, $config_hash);
    #----------------------------------------
    # Extract arguments from parameter hash
    #
    #----------------------------------------
    if (exists $parameter_hash->{'filename'}){
	$tag = $parameter_hash->{'filename'};
    }
    if (exists $parameter_hash->{'contig'}){
	$genome = $parameter_hash->{'contig'};
    }
    if (exists $parameter_hash->{'config_hash'}){
	$config_hash = $parameter_hash->{'config_hash'};
    }

    #----------------------------------------
    # Verify whether arguments were defined
    #
    #----------------------------------------
    if (!defined($tag)){
	$logger->logdie("tag was not defined");
    }
    if (!defined($genome)){
	$logger->logdie("genome was not defined");
    }
    if (!defined($config_hash)){
	$logger->logdie("config_hash was not defined");
    }

    my $filename   = $config_hash->{'refins'};
    my $outputfile = $config_hash->{'INDEL_Header'};
    my $filename2  = $config_hash->{'queryins'};

my @filetoparse;
my $id;
my $fileseq;
my @dna;
my $dna;
my @temp3;
my @temp4;
my $directory;
my $fasta;
my $posR;
my $posQ;
my $bpR;
my $bpQ;
my %hash;
my %hash2;
my @list;
my @filetoparse2;
my $posR2;
my $posQ2;
my $bpR2;
my $bpQ2;
my %hash3;
my %hash4;



    $genome = "/home/jravel/src/SNP/" . $genome;
    $directory = $tag . $config_hash->{'asmbl_seq_dir_mask'};

    my $gseq = &get_sequence(\$genome);
    if (!defined($gseq)){
	$logger->logdie("gseq was not defined");
    }


#### READ THE FILE INTO A ARRAY


## 1. READ THE DELETION ON QUERY

### FORMAT: REF    bp      bp     QUERY

########## 67495   t       -       23707
########## 73488   g       -       29700
########## 73489   t       -       29700
########## 73490   a       -       29700
########## 73491   t       -       29700
########## 73492   c       -       29700
########## 73493   t       -       29700
########## 73494   t       -       29700
########## 73495   t       -       29700
########## 29481   t       -       4969


    my $filetoparse = &get_contents(\$filename);
    if (!defined($filetoparse)){
	$logger->logdie("filetoparse was not defined");
    }

## READ THE INSERTION ON 
##EXTRACT INFORMATION FOR DELETION ON

    my $i = 1;
    my $y = 0;

    foreach my $line (@$filetoparse) {
    
	++$y;
	chomp $line;

	## SPLIT EACH FIELD

	my @field = split("\t", $line);


	$posR = $field[0];
	chomp $posR;
	$bpR = $field[1];
	chomp $bpR;
	$bpQ = $field[2];
	chomp $bpQ;
	$posQ = $field[3];
	chomp $posQ;
	

	### CREATE THE HASH
	
	$hash{$y} = $posR;
	
	if ($y == 1) {
	    print $i, "\n";
	    print $line, "\n";
	    push(@{$hash2{$i}}, $line);
	    next;
	    
	}else{
	    
	    if (($hash{$y} - $hash{$y-1}) == 1) {
		print $i, "\n";
		print $line, "\n";
		push(@{$hash2{$i}}, $line);
		
	    }else{
		++$i;
		print $i, "\n";
		print $line, "\n";
		push(@{$hash2{$i}}, $line);
	}
	    
	}
    }
    
### GROUPING THE INDELs
    
    my $z;
    my @temp;
    
    for ($z = 1; $z < $i+1; ++$z) {
	
	if ($z == 1) {
	    @temp3 = split("\t", ${$hash2{$z}}[0]);
	    $id = $temp3[4];
	    $fileseq = "$directory/$id";
	    
	    
	    my $dna = &get_sequence(\$fileseq);
	    if (!defined($dna)){
		$logger->logdie("dna was not defined");
	    }
	    
	}elsif ($z > 1) {
	    
	    @temp3 = split("\t", ${$hash2{$z}}[0]);
	    @temp4 = split("\t", ${$hash2{$z-1}}[0]);
	    
	    if ($temp3[4] eq $temp4[4]) {
		my $dna2 = $dna;
		
	    }elsif ($temp3[4] ne $temp4[4]){
		$id = $temp3[4];
		$fileseq = "$directory/$id";
		
		@dna = get_file_data($fileseq);
		
		$dna =  extract_sequence_from_fasta_data(@dna);
	    }
	}
	
	print "\n\nINDEL-REF-$z\n\n";
	my $w = 1;
	my $count = @{$hash2{$z}}; 
	my $dash = "|" x $count;
	my $dash2 = "-" x $count;
	my $seq2 = quickcut($gseq,$temp3[0]-20, $temp3[0]+$count+20); 
	my $seq3 = quickcut($dna, $temp3[3]-20, $temp3[3]-1);
	my $seq4 = quickcut($dna, $temp3[3], $temp3[3]+20);
	print "$temp3[3] for $count on assembly " . "$temp3[4]\n";
	print "R: " . $seq2, "\n";
	print "Q: " . $seq3 . $dash2 . "$seq4\n";
	push(@list, "INDEL-REF-$z\t$temp3[0]\t$count\t$temp3[3]\t$temp3[4]");
	
    }


    print "\n\n";

## 2. READ THE DELETION ON REF

### FORMAT: REF    bp      bp     QUERY


########## 11603   -       t       4781
########## 18468   -       t       11647
########## 22907   -       t       16087
########## 24644   -       g       17825
########## 58899   -       a       15092
########## 58899   -       t       15093
########## 58899   -       t       15094


    my $filetoparse2 = &get_contents(\$filename2);
    if (!defined($filetoparse2)){
	$logger->logdie("filetoparse2 was not defined");
    }
    
## READ THE INSERTION ON 
##EXTRACT INFORMATION FOR DELETION ON
    
    my $i2 = 1;
    my $y2 = 0;
    foreach my $line2 (@$filetoparse2) {
	
	++$y2;
	chomp $line2;
	
	## SPLIT EACH FIELD
	
	my @field2 = split("\t", $line2);
	
	
	$posR2 = $field2[0];
	chomp $posR2;
	$bpR2 = $field2[1];
	chomp $bpR2;
	$bpQ2 = $field2[2];
	chomp $bpQ2;
	$posQ2 = $field2[3];
	chomp $posQ2;
	
	
	### CREATE THE HASH
	
	$hash3{$y2} = $posQ2;
	
	if ($y2 == 1) {
	    print $i2, "\n";
	    print $line2, "\n";
	    push(@{$hash4{$i2}}, $line2);
	    next;
	    
	}
	else{
	    
	    if (($hash3{$y2} - $hash3{$y2-1}) == 1) {
		print $i2, "\n";
		print $line2, "\n";
		push(@{$hash4{$i2}}, $line2);
		
	    }
	    else{
		++$i2;
		print $i2, "\n";
		print $line2, "\n";
		push(@{$hash4{$i2}}, $line2);
	    }
	    
	}
    }
    
### GROUPING THE INDELs
    
    my $z2;
    my @temp2;
    my @temp5;
    my @temp6;
    my $fileseq2;
    
    for ($z2 = 1; $z2 < $i2+1; ++$z2) {
	
	if ($z2 == 1) {
	    @temp5 = split("\t", ${$hash4{$z2}}[0]);
	    $id = $temp5[4];
	    $fileseq2 = "$directory/$id";
	    
	    
	    @dna = get_file_data($fileseq);
	    
	    $dna = extract_sequence_from_fasta_data(@dna);
	    
	}elsif ($z2 > 1) {
	    
	    @temp5 = split("\t", ${$hash4{$z2}}[0]);
	    @temp6 = split("\t", ${$hash4{$z2-1}}[0]);
	    
	    if ($temp5[4] eq $temp6[4]) {
		
		
	    }elsif ($temp5[4] ne $temp6[4]){
		$id = $temp5[4];
	    $fileseq = "$directory/$id";
		
	    @dna = get_file_data($fileseq);
		
		$dna =  extract_sequence_from_fasta_data(@dna);
	    }
	}
	
	
	
	
	print "\n\nINDEL-QUE-$z2\n\n";
	my $w2 = 1;
	my $count2 = @{$hash4{$z2}}; 
	my $count3 =  @{$hash4{$i2}};
	my $dash3 = "|" x $count2;
	my $dash4 = "-" x $count2;
	
	my $seq5 = quickcut($dna, $temp5[3]-20, $temp5[3]+$count2+20); 
	my $seq6 = quickcut($gseq, $temp5[0]-20, $temp5[0]-1);
	my $seq7 = quickcut($gseq, $temp5[0], $temp5[0]+20);
	print "$temp5[0] for $count2 on assembly " . "$temp5[4]\n";
	print "R: " . $seq6 . $dash4 .  "$seq7\n";
	print "Q: $seq5\n";
	push(@list, "INDEL-QUE-$z2\t$temp5[0]\t$count2\t$temp5[3]\t$temp5[4]");
    }
    
    print "\n\n";
    

##### OUTPUT
    
    unless (open(FILEOUT, ">$outputfile")) {
	print STDERR "Can't open file!!\n\n";
	exit;
    }
    
    foreach my $del (@list) { 
	print FILEOUT $del, "\n";
    }
    
    close (FILEOUT);
    
}#end sub indel()


#----------------------------------------------------------
# indel_info()
#
#
#----------------------------------------------------------
sub indel_info {

    my logger = get_logger("snp");
    $logger->info("Entered indel_info");

    my (%parameter) = @_;
    my $parameter_hash = \%parameter;
    
    my ($tag, $genome, $config_hash);

    #----------------------------------------
    # Extract arguments from parameter hash
    #
    #----------------------------------------
    if (exists $parameter_hash->{'filename'}){
	$tag = $parameter_hash->{'filename'};
    }
    if (exists $parameter_hash->{'contig'}){
	$genome = $parameter_hash->{'contig'};
    }
    if (exists $parameter_hash->{'config_hash'}){
	$config_hash = $parameter_hash->{'config_hash'};
    }

    #----------------------------------------
    # Verify whether arguments were defined
    #
    #----------------------------------------
    if (!defined($tag)){
	$logger->logdie("tag was not defined");
    }
    if (!defined($genome)){
	$logger->logdie("genome was not defined");
    }
    if (!defined($config_hash)){
	$logger->logdie("config_hash was not defined");
    }


    my $filename    = $config_hash->{'refins'};
    my $outputfile  = $config_hash->{'INDEL_Header'};
    my $outputfile2 = $config_hash->{'INDEL_info'};
    my $filename2   = $config_hash->{'queryins'};

    my @filetoparse;
    my $id;
    my $fileseq;
    my @dna;
    my $dna;
    my $dnalen;
    my @temp3;
    my @temp4;
    my $directory;
    my $fasta;
    my $posR;
    my $posQ;
    my $bpR;
    my $bpQ;
    my %hash;
    my %hash2;
    my @list;
    my @filetoparse2;
    my $posR2;
    my $posQ2;
    my $bpR2;
    my $bpQ2;
    my %hash3;
    my %hash4;



    unless(open (INDEL, ">$outputfile2")) {
	print " Cannot open file $outputfile2\n\n";
	exit(0);
    }
    

    $genome = "/home/jravel/src/SNP/" . $genome;
    $directory = $tag . "_asmbl_seq.dir";


    # Extract the sequence data from the genome file
    # and remove blank spaces
    #
    my $gseq = &get_sequence(\$genome);
    if (!defined($gseq)){
	$logger->logdie("gseq was not defined");
    }

#### READ THE FILE INTO A ARRAY


## 1. READ THE DELETION ON QUERY

### FORMAT: REF    bp      bp     QUERY

########## 67495   t       -       23707
########## 73488   g       -       29700
########## 73489   t       -       29700
########## 73490   a       -       29700
########## 73491   t       -       29700
########## 73492   c       -       29700
########## 73493   t       -       29700
########## 73494   t       -       29700
########## 73495   t       -       29700
########## 29481   t       -       4969


    my $filetoparse = &get_contents(\$filename);
    if (!defined($filetoparse)){
	$logger->logdie("filetoparse was not defiend");
    }

## READ THE INSERTION ON 
##EXTRACT INFORMATION FOR DELETION ON

    my $i = 1;
    my $y = 0;
    foreach my $line (@$filetoparse) {
	
	++$y;

	chomp $line;

    ## SPLIT EACH FIELD
	
	my @field = split("\t", $line);


	$posR = $field[0];
	chomp $posR;
	$bpR = $field[1];
	chomp $bpR;
	$bpQ = $field[2];
	chomp $bpQ;
	$posQ = $field[3];
	chomp $posQ;
    

    ### CREATE THE HASH

	$hash{$y} = $posR;
	
	if ($y == 1) {
	    print $i, "\n";
	    print $line, "\n";
	    push(@{$hash2{$i}}, $line);
	    next;
	    
	}
	else{
	    
	    if (($hash{$y} - $hash{$y-1}) == 1) {
		print $i, "\n";
		print $line, "\n";
		push(@{$hash2{$i}}, $line);
		
	    }
	    else{
		++$i;
		print $i, "\n";
		print $line, "\n";
		push(@{$hash2{$i}}, $line);
	    }
	    
	}
    }
    
### GROUPING THE INDELs
    
    my $z;
    my @temp;
    
    for ($z = 1; $z < $i+1; ++$z) {
	
	if ($z == 1) {
	    @temp3 = split("\t", ${$hash2{$z}}[0]);
	    $id = $temp3[4];
	    $fileseq = "$directory/$id";
	    
	    
	    @dna = get_file_data($fileseq);
	    
	    $dna = extract_sequence_from_fasta_data(@dna);
	    
	    $dnalen = length($dna);
	    
	}
	elsif ($z > 1) {
	    
	    @temp3 = split("\t", ${$hash2{$z}}[0]);
	    @temp4 = split("\t", ${$hash2{$z-1}}[0]);
	    
	    if ($temp3[4] eq $temp4[4]) {
		my $dna2 = $dna;
		
	    }
	    elsif ($temp3[4] ne $temp4[4]){
		$id = $temp3[4];
		$fileseq = "$directory/$id";
		
		@dna = get_file_data($fileseq);
		
		$dna =  extract_sequence_from_fasta_data(@dna);
		
		$dnalen = length($dna);
		
	    }
	}
	
	print "\n\nINDEL-REF-$z\n\n";
	my $w = 1;
	my $count = @{$hash2{$z}}; 
	my $dash = "|" x $count;
	my $dash2 = "-" x $count;
	my $end5Q1 = $temp3[3]-20;
	my $end3Q1 = $temp3[3]+20;
	my $end5R1 = $temp3[0]-20;
	my $end3R1 = $temp3[0]+$count+20;
	my $seq2 = quickcut($gseq,$temp3[0]-20, $temp3[0]+$count+20); 
	my $seq3 = quickcut($dna, $temp3[3]-20, $temp3[3]-1);
	my $seq4 = quickcut($dna, $temp3[3], $temp3[3]+20);
	print INDEL "INDEL-REF-$z:$temp3[0] $count bp " . "$temp3[4] ($dnalen bp) INS on REF\n";
	print INDEL "R: $end5R1 - $end3R1 Q: $end5Q1 - $end3Q1\n";
	print INDEL "R: " . $seq2, "\n";
	print INDEL "Q: " . $seq3 . $dash2 . "$seq4\n//\n";
	push(@list, "INDEL-REF-$z\t$temp3[0]\t$count\t$temp3[3]\t$temp3[4]");
	
    }
    
    
    print "\n\n";
    
## 2. READ THE DELETION ON REF
    
### FORMAT: REF    bp      bp     QUERY


########## 11603   -       t       4781
########## 18468   -       t       11647
########## 22907   -       t       16087
########## 24644   -       g       17825
########## 58899   -       a       15092
########## 58899   -       t       15093
########## 58899   -       t       15094


    my $filetoparse2 = &get_contents(\$filename2);
    if (!defined($filetoparse2)){
	$logger->logdie("filetoparse2 was not defined");
    }

## READ THE INSERTION ON 
##EXTRACT INFORMATION FOR DELETION ON

    my $i2 = 1;
    my $y2 = 0;
    foreach my $line2 (@$filetoparse2) {
	
	++$y2;
	
	chomp $line2;
	
    ## SPLIT EACH FIELD
	
	my @field2 = split("\t", $line2);
	
	
	$posR2 = $field2[0];
	chomp $posR2;
	$bpR2 = $field2[1];
	chomp $bpR2;
	$bpQ2 = $field2[2];
	chomp $bpQ2;
	$posQ2 = $field2[3];
	chomp $posQ2;
	
	
	### CREATE THE HASH
	
	$hash3{$y2} = $posQ2;
	
	if ($y2 == 1) {
	    print $i2, "\n";
	    print $line2, "\n";
	    push(@{$hash4{$i2}}, $line2);
	    next;
	    
	}else{
	    
	    if (($hash3{$y2} - $hash3{$y2-1}) == 1) {
		print $i2, "\n";
		print $line2, "\n";
		push(@{$hash4{$i2}}, $line2);
		
	    }else{
		++$i2;
		print $i2, "\n";
		print $line2, "\n";
		push(@{$hash4{$i2}}, $line2);
	    }
	    
	}
    }
    
### GROUPING THE INDELs
    
    my $z2;
    my @temp2;
    my @temp5;
    my @temp6;
    my $fileseq2;
    
    for ($z2 = 1; $z2 < $i2+1; ++$z2) {
	
	if ($z2 == 1) {
	    @temp5 = split("\t", ${$hash4{$z2}}[0]);
	    $id = $temp5[4];
	    $fileseq2 = "$directory/$id";
	    
	
	    @dna = get_file_data($fileseq);
	    
	    $dna = extract_sequence_from_fasta_data(@dna);
	    
	    $dnalen = length($dna);
	    
	}elsif ($z2 > 1) {
	    
	    @temp5 = split("\t", ${$hash4{$z2}}[0]);
	    @temp6 = split("\t", ${$hash4{$z2-1}}[0]);
	    
	    if ($temp5[4] eq $temp6[4]) {
		

	    }elsif ($temp5[4] ne $temp6[4]){
		$id = $temp5[4];
		$fileseq = "$directory/$id";
		
		@dna = get_file_data($fileseq);
		
		$dna =  extract_sequence_from_fasta_data(@dna);
		
		$dnalen = length($dna);
	    }
	}
	
	


	print "\n\nINDEL-QUE-$z2\n\n";
	my $w2 = 1;
	my $count2 = @{$hash4{$z2}}; 
	my $count3 =  @{$hash4{$i2}};
	my $dash3 = "|" x $count2;
	my $dash4 = "-" x $count2;
	my $end5Q = $temp5[3]-20; 
	my $end3Q = $temp5[3]+$count2+20;
	my $end5R = $temp5[0]-20;
	my $end3R = $temp5[0]+20;
	my $seq5 = quickcut($dna, $temp5[3]-20, $temp5[3]+$count2+20); 
	my $seq6 = quickcut($gseq, $temp5[0]-20, $temp5[0]-1);
	my $seq7 = quickcut($gseq, $temp5[0], $temp5[0]+20);
	print INDEL "INDEL-QUE-$z2:$temp5[0] $count2 bp " . "$temp5[4] ($dnalen bp) INS on QUE\n";
	print INDEL "R: $end5R - $end3R Q: $end5Q - $end3Q\n";
	print INDEL "R: " . $seq6 . $dash4 .  "$seq7\n";
	print INDEL "Q: $seq5\n//\n";
	push(@list, "INDEL-QUE-$z2\t$temp5[0]\t$count2\t$temp5[3]\t$temp5[4]");
    }
    
print "\n\n";
    
    
##### OUTPUT
    
    unless (open(FILEOUT, ">$outputfile")) {
	print STDERR "Can't open file!!\n\n";
	exit;
    }
    
    foreach my $del (@list) { 
	print FILEOUT $del, "\n";
    }
    
    close (FILEOUT);
    close (INDEL);
    


}#end sub indel_info()

#----------------------------------------------------------
# daddy4()
#
#
#----------------------------------------------------------
sub daddy4 {

    my logger = get_logger("snp");
    $logger->info("Entered daddy4");

    my (%parameter) = @_;
    my $parameter_hash = \%parameter;
    
    my ($dir_prefix, $contig, $config_hash);


    #----------------------------------------
    # Extract arguments from parameter hash
    #
    #----------------------------------------
    if (exists $parameter_hash->{'dir_prefix'}){
	$dir_prefix = $parameter_hash->{'dir_prefix'};
    }
    if (exists $parameter_hash->{'contig'}){
	$contig = $parameter_hash->{'contig'};
    }
    if (exists $parameter_hash->{'config_hash'}){
	$config_hash = $parameter_hash->{'config_hash'};
    }

    #----------------------------------------
    # Verify whether arguments were defined
    #
    #----------------------------------------
    if (!defined($dir_prefix)){
	$logger->logdie("dir_prefix was not defined");
    }
    if (!defined($contig)){
	$logger->logdie("contig was not defined");
    }
    if (!defined($config_hash)){
	$logger->logdie("config_hash was not defined");
    }

    my ( $quality_file, $align_file, $header_file, $snpout_file, $indelR_file, $indelQ_file, $source);

    #---------------------------------------
    # Extract filenames from config hash
    #
    #---------------------------------------
    if (exists $config_hash->{'quality_file'}){
	$quality_file = $config_hash->{'quality_file'};
    }
    if (exists $config_hash->{'align_file'}){
	$align_file = $config_hash->{'align_file'};
    }
    if (exists $config_hash->{'snp_header'}){
	$header_file = $config_hash->{'snp_header'};
    }
    if (exists $config_hash->{'snpout'}){
	$snpout_file = $config_hash->{'snpout'};
    }
    if (exists $config_hash->{'indelR'}){
	$indelR_file = $config_hash->{'indelR'};
    }
    if (exists $config_hash->{'indelQ'}){
	$indelQ_file = $config_hash->{'indelQ'};
    }
    if (exists $config_hash->{'source'}){
	$source = $config_hash->{'source'};
    }


    my $asmbl_seq_dir = $dir_prefix . ${$config_hash->{'asmbl_seq_dir_mask'}};
    my $mum_align_dir = $dir_prefix . ${$config_hash->{'mum_align_dir_mask'}};
    my $contig_dir    = $dir_prefix . ${$config_hash->{'contig_dir_mask'}};
    my $cov_dir       = $dir_prefix . ${$config_hash->{'cov_dir_mask'}};

    #---------------------------------------------------------
    # Open and read all files from the assembly sequence 
    # directory
    #
    #---------------------------------------------------------
    unless(opendir(DIRECTORY, "$asmbl_seq_dir")) {
	$logger->logdie("Cannot open $directory");
    }
    
    # READ ALL FILES BUT THE . and .. FILES	
	
    @files = grep (!/^\.\.?$/, readdir(DIRECTORY));
    
    closedir(DIRECTORY);

    # GO THROUGH EACH FILENAME AND RUN THE SCRIPT

    foreach my $file (@files) {	
	

	if ($file =~ /^ID/) {
	    chomp $file;
	    print $file, "\n\n";
	    (my $file2 = $file);
	    $file2 =~ s/^ID//;
	    (my $file3 = $file2) =~ s/revcom$//;
	    print "FILE3 $file3\n";

	    #### PIPE THE FILE INTO MUMer
	
	    my $snp_fasta_outfile= "SNP_" . $dir_prefix . ".fasta";

	    my $genome_file = $source . "/" . $contig;
	    my $fasta_file  = $asmbl_seq_dir . "/" . $file;
	    my $info_file   = $mum_align_dir . "/" . $file . ".align";
	    my $contig_file = $contig_dir . "/" . $file3 . ".contig";
	    my $asmbl_id    = $file3;
	    my $t_cov       = $cov_dir . "/" . $file3 . ".tcov";
	    my $out_file    = $snp_fasta_outfile;

	    my $output_ref = quality_SNP8(
					  genome_file => \$genome_file,
					  fasta_file  => \$fasta_file,
					  info_file   => \$info_file,
					  contig_file => \$contig_file,
					  dir_prefix  => \$dir_prefix,
					  asmbl_id    => \$asmbl_id,
					  t_cov       => \$t_cov,
					  out_file    => \$out_file,
					  config_hash => $config_hash
					  );
#	    system("~/src/SNP/quality_SNP8.pl -G ~/src/SNP/$con -F $$asmbl_seq_dir/$file -I $$mum_align_dir/$file.align -C $$contig_dir/$file3.contig -T $$filename -D $file3 -V $$cov_dir/$file3.tcov -O $output");
	
	}#end if
    }#end foreach
}#end sub daddy4{}


#-----------------------------------------------------------
# quality_SNP8()
#
#
#-----------------------------------------------------------
sub quality_SNP8 {
    
    my $logger = get_logger("snp");
    $logger->info("Entered quality_SNP8");

    my %parameter = @_;
    my $parameter_hash = \%parameter;


##### DECLARE VARIABLE

    my ($pos1);
    my ($pos2);
    my ($help);
    my ($fasta); 
    my (@newarray);
    my ($entry2);
    my ($entry3);
    my ($entry4);
    my (@header);
    my ($apos1);
    my ($apos2);
    my ($gpos1);
    my ($gpos2);
    my (@align);
    my (@galign);
    my ($subgseq);
    my ($subaseq);
    my (@quality);
    my ($posM);
    my (@temp);
    my $qual2;
    my $sum;
    my $comment;
    my $len;
    my $coverage;
    my $val;
    my $lval;
    my $lval2;
    my $pval;
    my $elem;
    my $val2;
    my $star;
    my $tcov;
    my @coverage;
    my $newcov;
    my @cov3;
    my $tag;
    my @indelR;
    my @indelQ;
    my @snpseq;


    #--------------------------------------------------------------
    # Extract arguments from parameter hash
    #
    #-------------------------------------------------------------
                                                          # Jacques comments:                                  | Jay's comments: Examples:
    my $genome      = $parameter_hash->{'genome_file'};   # REFERENCE GENOME CHANGE HERE TO ENTER AT PROMPT    | "~src/SNP/GBA6612.1con"
    my $fasta       = $parameter_hash->{'fasta_file'};    # ASSEMBLY FILE IN FASTA FORMAT                      | "anthrax_test_asmbl_seq.dir/ID1090"
    my $readfile    = $parameter_hash->{'info_file'};     # MUMMER .align file                                 | "anthrax_test_mum_align.dir/ID1090.align"
    my $contigfile  = $parameter_hash->{'contig_file'};   # .config FIEL FROM PULLCONTIG WILL NEED TO BE PIPED | "anthrax_test_contig.dir/1090.contig"
    my $tag         = $parameter_hash->{'dir_prefix'};    # e.g. "anthrax_test"                                | "anthrax_test"
    my $asmbl_id    = $parameter_hash->{'asmbl_id'};      # THE ASSEMBLY ID NUMBER FOR CutAsm                  | "1090"
    my $tcov        = $parameter_hash->{'t_cov'};         #                                                    | "anthrax_test_cov.dir/1090.tcov"
    my $outputfile  = $parameter_hash->{'out_file'};      #                                                    | "SNP_anthrax_test.fasta"
    my $config_hash = $parameter_hash->{'config_hash'};   #                                                    | configuration hash

    #--------------------------------------------------------------
    # Extract output file names from configuration hash
    #
    #-------------------------------------------------------------
    my $qualityfile = $config_hash->{'quality_file'};
    my $alignfile   = $config_hash->{'align_file'};
    my $header      = $config_hash->{'header_file'};
    my $snpout      = $config_hash->{'snpout'};
    my $indelR      = $config_hash->{'indelR'};
    my $indelQ      = $config_hash->{'indelQ'};

    #
    # Additional variables
    #
    my $length = $config_hash->{'length'};


    #-------------------------------------------------
    # Verify whether the variables are defined
    #
    #-------------------------------------------------
    if (!defined($fasta)){
	$logger->logdie("fasta was not defined");
    }
    if (!defined($genome)){
	$logger->logdie("genome was not defined");
    }
    if (!defined($readfile)){
	$logger->logdie("readfile was not defined");
    }
    if (!defined($contigfile)){
	$logger->logdie("contigfile was not defined");
    }
    if (!defined($tag)){
	$logger->logdie("tag was not defined");
    }
    if (!defined($asmbl_id)){
	$logger->logdie("asmbl_id was not defined");
    }
    if (!defined($tcov)) {
	$logger->logdie("tcov was not defined");
    }
    if (!defined($outputfile)){
	$logger->logdie("outputfile was not defined");
    }
    if (!defined($qualityfile)){
	$logger->logdie("qualityfile was not defined");
    }
    if (!defined($alignfile)){
	$logger->logdie("alignfile was not defined");
    }
    if (!defined($header)){
	$logger->logdie("header was not defined");
    }
    if (!defined($snpout)){
	$logger->logdie("snpout was not defined");
    }
    if (!defined($indelR)){
	$logger->logdie("indelR was not defined");
    }
    if (!defined($indelQ)){
	$logger->logdie("indelQ was not defined");
    }
    if (!defined($length)){
	$logger->logdie("length was not defined");
    }


    #---------------------------------------------------------------
    # Retrieve the assembly sequence from fasta file
    #---------------------------------------------------------------
    my $assembly_seq = &get_sequence(\$fasta);
    if (!defined($assembly_seq)){
	$logger->logdie("assembly_seq was not defined");
    }

    my $assembly_seq_length = length($assembly_seq);
    print ("$assembly_seq_length\n"); # length of the contig

    #---------------------------------------------------------------
    # Retrieve the genome sequence from fasta file
    #---------------------------------------------------------------
    my $genome_seq = &get_sequence(\$genome);
    if (!defined($genome_seq)){
	$logger->logdie("genome_seq was not defined");
    }

    my $genome_seq_length = length($genome_seq);


    my $tagname = $tag . $config_hash->{'asmbl_seq_dir_mask'}; #### NAME OF THE FOLDER WHERE THE FASTA FILES ARE LOCATED
    
    $fasta =~ s/^$tagname\///;        ### REMOVE THE NAME OF THE FOLDER FROM THE NAME FILE
    
    
    #--------------------------------------------------------------------
    # OPEN THE COVERAGE FILE FOR THIS CONTIG AND STORE IT INTO AN ARRAY
    # THIS ARRAY WILL BE LOOKED UP WHEN THE COVERAGE/QUALITY VALUES WILL 
    # BE FORMATTED
    #
    #--------------------------------------------------------------------
    my $coverage_contents = &get_contents(\$tcov);
    if (!defined($coverage_contents)){
	$logger->logdie("coverage_contents was not defined");
    }

    my ($t, @at);
    
    for ($t=0;$t<@coverage;$t++) {
	
	if ($coverage[$t] =~ /^(\d+)/ && !$at[$1]) {
	    $at[$1] = $t;
	}
    }
    
    #-----------------------------------------------------
    # Retrieve mummer output from readfile
    # with order re-arranged
    #-----------------------------------------------------
    my $change_align = &get_change_align(\$readfile);
    if (!defined($change_align)){
	$logger->logdie("change_align was not defined");
    }

    my (@read);
    my $i=1;
    
    foreach my $array_ref (@$change_align) {
       
	@read = @$array_ref;

	### DOES NOT PROCESS GAPS (INSERTION/DELETION)
	### IF . IN SECOND COLUMN (REF) SKIP
       
	if ($read[1] =~ /-/){
	    my $line = join("\t", @read);
	    $line .= "\t$fasta";
	    push(@indelQ, $line);
	    next;
	}

        ### DOES NOT PROCESS GAPS (INSERTION/DELETION)
	### IF . IN SECOND COLUMN (REF) SKIP

	if ($read[2] =~ /-/){
	    my $line = join("\t", @read);
	    $line .= "\t$fasta";
	    push(@indelR, $line);
	    next;
	}

	### DOES NOT PROCESS AMBIGUOUS BP

	if ($read[2] =~/[mrwsykn]/) {
	    next;
	}

	### LAST LINE SKIP

	if ($read[0] =~ /Total/){
	    next;

	### IF A SNP PROCESS:
    
	}else{
	    
	    ### DETERMINE THE POSITION OF THE SNP COMPARE TO THE ALL ASMBL
	    ### TEST IF ASSEMBLY LESS THAN 500 and CALCULATE THE CUT-OUT ACCORDINGLY
	    
	    if ($read[3] > 500 and $read[3] < ($length2 - 500)){
		$pos1 = $read[3] - 500; ##FOR FASTA FILE 500 bp around the SNP
		$pos2 = $read[3] + 500;
		$comment = "SNP: 501/1001";
	    }elsif ($read[3] < 500 and $read[3] < ($length2 - 500)) {   #### IF THE SNP IS A LESS THAN 500 bp 
		$pos1 = 1;                                              #### FROM THE START OF ASMBL
		$pos2 = $read[3] + 500;
		$comment = "SNP: $read[3]/$pos2";
	    }elsif ($read[3] > 500 and $read[3] > ($length2 - 500)) {   #### IF SNP AT LESS THAN 500 bp
		$pos1 = $read[3] - 500;                                 #### FROM THE END OF ASMBL
		$pos2 = $length2;
		$len = ($length2 - $read[3] + 501);
		$comment = "SNP: 501/$len";
	    }elsif ($read[3] < 500 and $read[3] > ($length2 - 500)) {   ### IF SNP AT LESS THAN 500 bp 
		$pos1 = 1;                                              ### FROM START OR END
		$pos2 = $length2;
		$len = $length2;;
		$comment = "SNP: $read[3]/$len";
	    }elsif ($read[3] < 500 and $length2 < 501) {
		$pos1 = 1;
		$pos2 = $length2;
		$comment = "SNP: $read[3]/$pos2";
	    }
		
	    $apos1 = $read[3] - 20; ## FOR ALIGN FILE 20 bp around the SNP ON ASSEMBLY
	    $apos2 = $read[3] + 20;
	    $gpos1 = $read[0] - 20; ## FOR ALIGN FILE 20 bp around the SNP ON GENOME
	    $gpos2 = $read[0] + 20;

	if ($fasta =~ /revcom/){  ## IF THE FASTA FILE CONTAIN revcom THEN COMPUTE COORDINATE
	    $posM = ($length2 - $read[3] + 1);  ### LENGTH OF ASSEMBLY - SNP pos + 1 (start at 0)

	}else{  ## IF ON FORWARD STRAND JUST KEEP THE POSITION OF THE SNP
	    $posM = $read[3];
        }

##################################################################
####     FORMAT THE OUTPUT
##################################################################
	    

	   
	    print $fasta;
	    print $posM, "\n\n";


	### FORMAT THE FILE TO GET QUALITY VALUES

	    push(@quality, "\n>" . $fasta . "-" . $i . " " .  $read[0] ." ". $read[1] ."-->". $read[2] ." ". $read[3]);
	    push(@snpseq,  ">" . $fasta . "-" . $i . " " .  $read[0] ." ". $read[1] ."-->". $read[2] ." ". $read[3]);

	    ### LOOK UP THE COVERAGE ARRAY AND PULL THE LINE FOR +/- 5 around the SNP

	    @cov3 = @coverage[$at[$posM-5]..$at[$posM+5]];
	    
	    foreach my $covline (@cov3) { ### GO THROUGH EACH LINE
		
		my @cov2 = split(" ", $covline);

		if (!defined $cov2[3]) { $cov2[3] = "NC";}
		
		if (!defined $cov2[4]) { $cov2[4] = "NC";}

		
		if ($cov2[0] == $posM) {  ### TAG THE LINE FOR THE SNP WITH *
		    
		    $pval = $cov2[2];
		    $lval = length($cov2[3]);
		    
		   
		    
		    $newcov = "*" . "$cov2[0] \t $cov2[1] \t $cov2[2] \t $cov2[3]   \t $cov2[4]";
		    push (@quality, $newcov);

		    push (@snpseq, $newcov);

		}else {  ### BASICALLY DO NOTHING HERE BUT ADD THE TAB

		    $newcov = " " . "$cov2[0] \t $cov2[1] \t $cov2[2] \t $cov2[3]   \t $cov2[4]";
		    push (@quality, $newcov);

		}

	    }  
		
	     
	   
	## HEADER

	my $alignhead = ">" . $fasta . "-" . $i." (".$length2."bp) " .  $read[0] ." ". $read[1] ."-->". $read[2] ." ". $read[3]." --- " . $apos1. " to " . $apos2 ." COV: " . $lval . " CB_QVal: ". $pval;;

	push(@galign, $alignhead);

	## GET THE SEQUENCE AROUND THE SNP 20 bp AROUND/ CALL TO QUICKCUT AS A SUBROUTINE IN SUB

	$subgseq = quickcut ($gseq, $gpos1, $gpos2);

	## FORMAT THE ALIGNMENT

	push(@galign, "REF " .  ": " . $subgseq . " " . $gpos2);
	push(@galign, "      |||||||||||||||||||| ||||||||||||||||||||");
	$subaseq = quickcut ($dna, $apos1, $apos2);
	push(@galign, "ASM " .  ": " . $subaseq . " " . $apos2. "\n");
       
	## FORMAT THE FASTA FILE OF SEQUENCE 500 bp AROUND THE SNP

	    ## FORMAT THE HEADER 

	my $entry = ">" . $fasta . "-" . $i." (".$length2."bp) " .  $read[0] ." ". $read[1] ."-->". $read[2] ." ". $read[3]." --- " . $pos1. " to " . $pos2 . " (" . $comment . ") COV: " . $lval . " CB_QVal: ". $pval;
	
	push(@newarray, $entry);

	push(@header, $entry);

            ## GET THE SEQUENCE AROUND THE SNP 500 bp on the ASSEMBLY

	$entry2 = quickcut ($dna, $pos1, $pos2);
	
	   ## FORMAT THE SEQUENCE IN LINES OF 100 bp (CAN BE CHANGED AT $length

	for ( my $pos = 0 ; $pos < length($entry2) ; $pos += $length ) {
        push(@newarray, substr($entry2, $pos, $length));
        }
    

	

	$i++;   ### COUNTER TO GIVE A NUMBER TO THE SNP
    }
    }
	
### PRINT TO FILES

    ### PRINT TO SNP_quality.txt

unless (open(QUALITY, ">>$qualityfile")) {
    print STDERR "Can't open file!!\n\n";
    exit;
}

foreach my $sub (@quality) { 
    print QUALITY $sub, "\n";
}

close (QUALITY);

    ### PRINT TO SNP_align.txt

unless (open(ALIGN, ">>$alignfile")) {
    print STDERR "Can't open file!!\n\n";
    exit;
}

foreach my $sub (@galign) { 
    print ALIGN $sub, "\n";
}

close (ALIGN);

    ### PRINT TO SNP_Header.txt

unless (open(HEAD, ">>$header")) {
    print STDERR "Can't open file!!\n\n";
    exit;
}

foreach my $head (@header) { 
    print HEAD $head, "\n";
}

close (HEAD);

    ### PRINT TO SNP_fasta.txt

unless (open(FILEOUT, ">>$outputfile")) {
    print STDERR "Can't open file!!\n\n";
    exit;
}

foreach my $line2 (@newarray) { 
    print FILEOUT $line2, "\n";
}

close (FILEOUT);

unless (open(SNPOUT, ">>$snpout")) {
    print STDERR "Can't open file!!\n\n";
    exit;
}

foreach my $line5 (@snpseq) { 
    print SNPOUT $line5, "\n";
}

close (FILEOUT);

unless (open(INDELR, ">>$indelR")) {
    print STDERR "Can't open file!!\n\n";
    exit;
}

foreach my $indel (@indelR) { 
    print INDELR $indel, "\n";
}

close (INDELR);

unless (open(INDELQ, ">>$indelQ")) {
    print STDERR "Can't open file!!\n\n";
    exit;
}

foreach my $delin (@indelQ) { 
    print INDELQ $delin, "\n";
}

close (INDELQ);



exit;





}#end sub quality_SNP8()





#----------------------------------------------------------
# remove_ins_files()
#
#
#----------------------------------------------------------
sub remove_ins_files {

    my $logger = get_logger("snp");
    $logger->info("Entered remove_ins_files");
    

    my %parameter = @_;
    my $parameter_hash = \%parameter;
    
    my $delete_list;
    #---------------------------------------------
    # Extract arguments from parameter hash
    #
    #---------------------------------------------
    if (exists $parameter_hash->{'file_list'}){
	$delete_list = $parameter_hash->{'file_list'};
    }

    #---------------------------------------------
    # Verify whether arguments were defined
    #
    #---------------------------------------------
    if (!defined($delete_list)){
	$logger->logdie("delete_list was not defined");
    }


    foreach my $file (@$delete_list){
	if (-e $file){
	    unlink ($file);
	}
    }

}#end sub remove_ins_files()

#----------------------------------------------------------
# sort_snp_file()
#
#
#----------------------------------------------------------
sub sort_snp_file {

    my $logger = get_logger("snp");
    $logger->info("Entered sort_snp_file");


    my %parameter = @_;
    my $parameter_hash = \%parameter;
    my ($snp_in, $snp_out);

    #-------------------------------------------------
    # Extract arguments from parameter hash
    #
    #-------------------------------------------------
    if (exists $parameter_hash->{'snp_header_file'}){
	$snp_in = $parameter_hash->{'snp_header_file'};
    }
    if (exists $parameter_hash->{'sort_snp_header_file'}){
	$snp_out = $parameter_hash->{'sort_snp_header_file'};
    }

    #-------------------------------------------------
    # Verify whether arguments were defined
    #
    #-------------------------------------------------
    if (!defined($snp_in)){
	$logger->logdie("snp_in was not defined");
    }
    if (!defined($snp_out)){
	$logger->logdie("snp_out was not defined");
    }

    &check_file_status(\$snp_in);

    system ("sort -k13nr $snp_in > $snp_out");



}#end sort_snp_file()


#----------------------------------------------------------------------------
# get_sequence()
#
#
#
#----------------------------------------------------------------------------
sub get_sequence {

    my $logger = get_logger("snp");
    $logger->info("Entered get_sequence");

    my $file = shift;
    if (!defined($file)){
	$logger->logdie("file was not defined");
    }

    my $contents_ref = &get_contents($file);
    if (!defined($contents_ref)){
	$logger->logdie("contents_ref was not defined");
    }

    my $cleaned_sequence_ref = &get_clean_sequence($contents_ref);
    if (!defined($cleaned_sequence_ref)){
	$logger->logdie("cleaned_sequence_ref was not defined");
    }

    return $cleaned_sequence_ref;

}#end sub get_sequence()


#----------------------------------------------------------------------------
# get_clean_sequence()
#
#
#
#----------------------------------------------------------------------------
sub get_clean_sequence {

    my $logger = get_logger("snp");
    $logger->info("Entered get_clean_sequence");
    
    my $array_ref = shift;
    if (!defined($array_ref)){
	$logger->logdie("array_ref was not defined");
    }

    # Declare and initialize variables
    my $sequence = '';

    foreach my $line (@$array_ref) {

        # discard blank line
        if ($line =~ /^\s*$/) {
            next;
        # discard comment line
        } elsif($line =~ /^\s*#/) {
            next;
        # discard fasta header line
        } elsif($line =~ /^>/) {
            next;
        # keep line, add to sequence string
        } else {
            $sequence .= $line;
        }
    }

    # remove non-sequence data (in this case, whitespace) from $sequence string
    $sequence =~ s/\s//g;

    return \$sequence;
}




}#end sub get_clean_sequence()



#----------------------------------------------------------------------------
# get_change_align()
#
#
#
#----------------------------------------------------------------------------
sub get_change_align {

    my $logger = get_logger("snp");
    $logger->info("Entered get_change_align");

    my $file = shift;
    if (!defined($file)){
	$logger->logdie("file was not defined");
    }

    my $contents_ref = &get_contents($file);
    if (!defined($contents_ref)){
	$logger->logdie("contents_ref was not defined");
    }

    my @in_array;
    my @out_array;
    my $out_line;
    my $marker = 0;
    my $i =0;

    foreach my $in_line (@$contents_ref) {
    
	if ($in_line =~ /^>/ && $marker == 0) {
	    $marker = 1;
	    next;
	} 
	elsif ( $in_line =~ /Ref/ && $marker == 1) {
	    my @tmp_array;
	    
	    @in_array = split([\s]+, $in_line);

	    $tmp_array[0] = $in_array[2];
	    $tmp_array[1] = $in_array[4];
	    $tmp_array[2] = $in_array[5];
	    $tmp_array[3] = $in_array[3];

	    push (@out_array, \@tmp_array);
	    $i++;
	    next;
	}
	elsif ($in_line =~ /^>/ && $marker == 1) {
	    last;
	}
    }

    push (@out_array, "Total number of errors = $i\n");

    return \@out_array;

}#end sub get_change_align()

#------------------------------------------------------------------------------
# get_contents()
#
#
#
#------------------------------------------------------------------------------
sub get_contents {

    my $logger = get_logger("snp");
    $logger->info("Entered get_contents");

    my $file = shift;
    if (!defined($file)){
	$logger->logdie("file was not defined");
    }

    &check_file_status($file);

    open (IN_FILE, "<$$file") || $logger->logdie("Could not open file: $$file for input");
    my @contents = <IN_FILE>;
    chomp @contents;

    return \@contents;


}#end sub get_contents()
