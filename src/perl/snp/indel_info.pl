#! /usr/local/bin/perl -w

use strict;
use lib "/home/jravel/lib/";
use BeginPerlBioinfo;
use Getopt::Long;
$| = 1;


my $filename = "refins.txt";
my @filetoparse;
my $genome;
my $id;
my $fileseq;
my @dna;
my $dna;
my $dnalen;
my @temp3;
my @temp4;
my $directory;
my $tag;
my $fasta;
my $posR;
my $posQ;
my $bpR;
my $bpQ;
my %hash;
my %hash2;
my @list;
my $outputfile = "INDEL_Header.txt";
my $outputfile2 = "INDEL_info.txt";
my $filename2 = "queryins.txt";
my @filetoparse2;

my $posR2;
my $posQ2;
my $bpR2;
my $bpQ2;
my %hash3;
my %hash4;

my $help;
my $HELPTEXT = "

USAGE -T folder tag -C .1con of reference

Folder Tag; used to create directoryies

.1con: located in the ~/src/SNP/ folder

HELP: -h (This text)
For questions and report bug please Email : jravel\@tigr.org

";


Getopt::Long::config("no_ignore_case");

my $result = GetOptions ("h|help"     => \$help,
			 "T=s"        => \$tag,
			 "C=s"        => \$genome);

if ($result == 0) {
    die ("Command line parsing failed\n")
    }

if ($help) {

    print $HELPTEXT, "\n";
    exit;
}

if (! defined $tag) {
    print STDERR "You must enter a folder tag\n";
    exit(0);
}
if (! defined $genome) {
    print STDERR "You must specify a .1con file\n";
    exit(0);
}


unless(open (INDEL, ">$outputfile2")) {
print " Cannot open file $outputfile2\n\n";
exit(0);
}


$genome = "/home/jravel/src/SNP/" . $genome;

$directory = $tag . "_asmbl_seq.dir";

## FOR THE REFERENCE GENOME

my @genome = get_file_data($genome);



#### Extract the sequence data from the array and remove blank space ####

my $gseq = extract_sequence_from_fasta_data(@genome);




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


@filetoparse = get_file_data($filename);

## READ THE INSERTION ON 


##EXTRACT INFORMATION FOR DELETION ON

my $i = 1;
my $y = 0;
foreach my $line (@filetoparse) {

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

	
	@dna = get_file_data($fileseq);
	
	$dna = extract_sequence_from_fasta_data(@dna);

        $dnalen = length($dna);
	
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


@filetoparse2 = get_file_data($filename2);

## READ THE INSERTION ON 


##EXTRACT INFORMATION FOR DELETION ON

my $i2 = 1;
my $y2 = 0;
foreach my $line2 (@filetoparse2) {

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





    exit;
