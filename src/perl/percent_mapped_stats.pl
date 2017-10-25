#!/usr/bin/env perl 
################################################################################
### This program generates coverage statistics from the sorted by name alignment 
### BAM file across the genic, exonic, intronic, and/or intergenic regions.
################################################################################
###percent_mapped_stats.pl

use strict;
use warnings;
use Cwd;
use Carp;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Spec;

##############################################################################
### Constants
##############################################################################

use constant FALSE => 0;
use constant TRUE  => 1;

use constant VERSION => '0.1.0';
use constant PROGRAM => eval { ($0 =~ m/(\w+\.pl)$/) ? $1 : $0 };

use constant BEDTOOLS_BIN_DIR => '/usr/local/packages/bedtools';
use constant SAMTOOLS_BIN_DIR => '/usr/local/bin';

##############################################################################
### Globals
##############################################################################

my %hCmdLineOption = ();
my $sHelpHeader = "\nThis is ".PROGRAM." version ".VERSION."\n";

my (@aRegions);
my ($sRefFile, $sBamFile, $sAnnotationFile);
my ($sGenicFile, $sGenicBedFile, $sExonicBedFile, $sIntronicBedFile, $sIntergenicBedFile);
my ($sOutDir, $sSizeFile, $sOutFile, $sSortedFile, $sSummaryFile);
my ($Exon_inter,$intergenic_inter,$intron_inter,$genic_inter);
my ($e_file, $in_file, $it_file, $fout, $finalfile, $uniqfile);
my ($sPrefix, $sSampleId, $sRegion, $nTotalMappedReads, $nUniqueMappedReads);
my ($sCmd, $sOpt);
my ($bDebug, $bVerbose);
my $out_flag = 0;

##############################################################################
### Main
##############################################################################


GetOptions( \%hCmdLineOption,
            'reffile|r=s', 'infile|i=s', 'annotation|a=s', 
            'annotationfiletype|t=s', 'feature|f=s', 'attribute|id=s', 'groupby|g=s', 
            'outdir|o=s', 'bedtools_bin_dir|b=s', 'samtools_bin_dir|s=s','org-type|org=s', 
            'verbose|v',
            'debug',
            'help',
            'man') or pod2usage(2);

if ($hCmdLineOption{'help'} || 
	(! defined $hCmdLineOption{'reffile'}) || 
	(! defined $hCmdLineOption{'infile'}) || 
	(! defined $hCmdLineOption{'annotation'}) ||
        (! defined $hCmdLineOption{'org-type'}) || 
	(! defined $hCmdLineOption{'annotationfiletype'})) { 
    pod2usage( -msg => $sHelpHeader, -exitval => 1);
}

if (($hCmdLineOption{'annotationfiletype'} !~ m/^bed$/i) && 
	($hCmdLineOption{'annotationfiletype'} !~ m/^gtf$/i) &&
	($hCmdLineOption{'annotationfiletype'} !~ m/^gff3$/i)) {
    die "\tERROR: Annotation file format needs to be a BED, GTF or GFF3 format file\n";
}

if (($hCmdLineOption{'annotationfiletype'} =~ m/^(gtf|gff3)$/i) && 
	((! defined $hCmdLineOption{'feature'}) || (! defined $hCmdLineOption{'attribute'}))) {
    die "\tERROR: 'feature' and 'attribute' required for GTF and GFF3 annotation file formats\n";
}

if ($hCmdLineOption{'org-type'} ne 'prok' && $hCmdLineOption{'org-type'} ne 'euk'){
    die "\tERROR: Enter organism type: (prokaryote|prok / eukaryote|euk )\n";
}

if ($hCmdLineOption{'org-type'} eq 'euk' && ! defined $hCmdLineOption{'groupby'}) {
    die "\tError: 'group-by' parameter required for eukaryotes.";}

pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};


$bDebug   = (defined $hCmdLineOption{'debug'}) ? TRUE : FALSE;
$bVerbose = (defined $hCmdLineOption{'verbose'}) ? TRUE : FALSE;

$sOutDir = File::Spec->curdir();
if (defined $hCmdLineOption{'outdir'}) {
    $sOutDir = $hCmdLineOption{'outdir'};

    if (! -e $sOutDir) {
        mkdir($hCmdLineOption{'outdir'}) ||
            croak "ERROR! Cannot create output directory\n";
    }
    elsif (! -d $hCmdLineOption{'outdir'}) {
            croak "ERROR! $hCmdLineOption{'outdir'} is not a directory\n";
    }
}
$sOutDir = File::Spec->canonpath($sOutDir);

if (! (defined $hCmdLineOption{'bedtools_bin_dir'}) ) {
	$hCmdLineOption{'bedtools_bin_dir'} = BEDTOOLS_BIN_DIR;
}

if (! (defined $hCmdLineOption{'samtools_bin_dir'}) ) {
	$hCmdLineOption{'samtools_bin_dir'} = SAMTOOLS_BIN_DIR;
}


# process reference file

$hCmdLineOption{'reffile'} = File::Spec->rel2abs($hCmdLineOption{'reffile'});
($_, $_, $sRefFile) = File::Spec->splitpath($hCmdLineOption{'reffile'});

symlink $hCmdLineOption{'reffile'}, "$sOutDir/$sRefFile";

if ( -e "$hCmdLineOption{'reffile'}.fai" ) {
	symlink "$hCmdLineOption{'reffile'}.fai", "$sOutDir/$sRefFile.fai";
}
else {
	($bDebug || $bVerbose) ? 
		print STDERR "\nIndexing $sRefFile .....\n" : ();
	
	$sCmd = $hCmdLineOption{'samtools_bin_dir'}."/samtools faidx".
			" ".$sOutDir."/".$sRefFile;
	
	exec_command($sCmd);
	
	($bDebug || $bVerbose) ? 
		print STDERR "Indexing $sRefFile ..... done\n" : ();
}

# process chromosome size file

$sPrefix = $sRefFile;
$sPrefix =~ s/\.([A-za-z])+$//;

$sSizeFile = File::Spec->rel2abs("$sOutDir/$sPrefix.chromosome.sizes");

($bDebug || $bVerbose) ? 
	print STDERR "\nGenerating $sSizeFile .....\n" : ();

$sCmd = "cut -f1,2 $sOutDir/$sRefFile.fai > $sSizeFile";

exec_command($sCmd);

($bDebug || $bVerbose) ? 
	print STDERR "Generating $sSizeFile ..... done\n" : ();

$hCmdLineOption{'infile'} = File::Spec->rel2abs($hCmdLineOption{'infile'});
($_, $_, $sBamFile) = File::Spec->splitpath($hCmdLineOption{'infile'});

$sSampleId = $sBamFile;
$sSampleId =~ s/\.accepted_hits.sorted_by_name//;
$sSampleId =~ s/\.bam$//;


if ($hCmdLineOption{'org-type'} eq 'euk') {

    ####Exonic Bed File
    $sSortedFile = Process_Annotation(\%hCmdLineOption, $hCmdLineOption{'annotation'}, $hCmdLineOption{'annotationfiletype'}, "exonic", $sOutDir);
   
    ####Genic Bed File
    $sGenicFile = Generate_Gene_BedFile(\%hCmdLineOption, $hCmdLineOption{'annotation'}, $hCmdLineOption{'annotationfiletype'}, "genic", $sOutDir);

    ####Intronic Bed file
    $sIntronicBedFile = Init_OutFileName(\%hCmdLineOption, $sOutDir, $sSortedFile, ".exonic.sorted.bed");
    $sIntronicBedFile .= '.intronic.sorted.bed';
    $sCmd = $hCmdLineOption{'bedtools_bin_dir'}."/bedtools subtract".
	  		    " -s -a ".$sGenicFile.
			    " -b ".$sSortedFile.
			    " > ".$sIntronicBedFile;
    exec_command($sCmd);

    
    if ( -z $sIntronicBedFile) {
    	($bDebug || $bVerbose) ? 
			print STDERR "WARNING! No introns detected in $hCmdLineOption{'annotation'}\n" : ();
        $hCmdLineOption{'org-type'} = "prok";
	$out_flag = 1; 
    }
}
else {    
    $sGenicFile = Process_Annotation(\%hCmdLineOption, $hCmdLineOption{'annotation'}, $hCmdLineOption{'annotationfiletype'}, "genic", $sOutDir);
}

$sGenicBedFile = Init_OutFileName(\%hCmdLineOption, $sOutDir, $sGenicFile, ".bed");
$sGenicBedFile .= '.sanitized.bed';
$sCmd = $hCmdLineOption{'bedtools_bin_dir'}."/bedtools merge".
			" -s -i ".$sGenicFile." > ".$sGenicBedFile;
exec_command($sCmd);

####Intergenic Bed file
$sIntergenicBedFile = Init_OutFileName(\%hCmdLineOption, $sOutDir, $sGenicBedFile, ".genic.sorted.sanitized.bed");
$sIntergenicBedFile .= '.intergenic.sorted.sanitized.bed';
$sCmd = $hCmdLineOption{'bedtools_bin_dir'}."/bedtools complement".
			" -i ".$sGenicBedFile.
			" -g ".$sSizeFile.
			" > ".$sIntergenicBedFile;
		
exec_command($sCmd);



###Bed tools intersect

if ($hCmdLineOption{'org-type'} eq 'euk') {
    $Exon_inter = Init_OutFileName(\%hCmdLineOption, $sOutDir, $sSortedFile, ".bed");
    $Exon_inter .= '.intersect.bed'; 
    bed_intersect(\%hCmdLineOption,$hCmdLineOption{'infile'},$sSortedFile,$Exon_inter); 
    $intron_inter = Init_OutFileName(\%hCmdLineOption, $sOutDir, $sIntronicBedFile, ".bed");
    $intron_inter .= '.intersect.bed';
    bed_intersect(\%hCmdLineOption,$hCmdLineOption{'infile'},$sIntronicBedFile,$intron_inter);  
    $intergenic_inter = Init_OutFileName(\%hCmdLineOption, $sOutDir, $sIntergenicBedFile, ".bed");
    $intergenic_inter .= '.intersect.bed';
    bed_intersect(\%hCmdLineOption,$hCmdLineOption{'infile'},$sIntergenicBedFile,$intergenic_inter);
}
else {
    $genic_inter = Init_OutFileName(\%hCmdLineOption, $sOutDir, $sGenicFile, ".bed");
    $genic_inter .= '.intersect.bed';
    bed_intersect(\%hCmdLineOption,$hCmdLineOption{'infile'},$sGenicFile,$genic_inter);  
    $intergenic_inter = Init_OutFileName(\%hCmdLineOption, $sOutDir, $sIntergenicBedFile, ".bed");
    $intergenic_inter .= '.intersect.bed';
    bed_intersect(\%hCmdLineOption,$hCmdLineOption{'infile'},$sIntergenicBedFile,$intergenic_inter);

}


###Call mapped subroutine to calculate % reads mapping to exon, intron and intergenic region

if ($hCmdLineOption{'org-type'} eq 'euk') {
    
    $finalfile = $sOutDir."/".$sSampleId.".Percent.txt";
    open ($e_file, "<$Exon_inter") or die "Error! Cannot open the file.";
    open ($in_file,"<$intron_inter") or die "Error! Cannot open the file.";
    open ($it_file,"<$intergenic_inter") or die "Error! Cannot open the file.";
    open ($fout,">$finalfile") or die "Error! Cannot open the file.";
    mapped_euk ($e_file, $in_file, $it_file,$fout);

##########For few test runs..Check unique reads obtained using mapped_euk with NH:i:1 option from bam file.
    $uniqfile = $sOutDir."/".$sSampleId.".uniq_reads.txt";
    $sCmd = $hCmdLineOption{'samtools_bin_dir'}."/samtools view -F 260 ".$hCmdLineOption{'infile'}." | grep NH:i:1 | wc -l > ".$uniqfile;
    exec_command($sCmd);

}
else {
    $finalfile= $sOutDir."/".$sSampleId.".Percent.txt";
    open ($e_file, "<$genic_inter") or die "Error! Cannot open the file.";
    open ($it_file,"<$intergenic_inter") or die "Error! Cannot open the file.";
    open ($fout,">$finalfile") or die "Error! Cannot open the file.";
    mapped_prok ($e_file,$it_file,$fout, $out_flag);

}

print "Removing BED files........\n";
$sCmd = "rm ".$sOutDir."/*.bed";
exec_command($sCmd);

exit;


##############################################################################
### Subroutines
##############################################################################

###Percent reads mapping to genic and intergenic region for prokaryotes.
sub mapped_prok {

my  $genic_file = shift;
my  $inter_file = shift;
my  $fout = shift;
my $out_flag = shift;
my  ($g_count,$it_count,$p_genic,$p_inter,$it_line);
my  ($k,$k1,$tag,$read1,$read_g,$read_it);
my  $t_genic = 0;
my  $t_inter = 0;  
my  $uniq_reads = 0;
my  $total_reads = 0;
my  $t_count = 0;
my  %hcount = ();

while (<$genic_file>) {
    chomp($_);
    $read_g = (split (/\t/,$_))[3];
    $read1 = (split (/\//,$read_g))[0];
    $tag = (split (/\//,$read_g))[1];
    if (! (defined $tag)){ $tag = 0;}
    $it_line = <$inter_file>;
    chomp($it_line);
    $read_it = (split (/\t/,$it_line))[3];
    ##Check if reads in genic and intergenic file is same.
    if ($read_g ne $read_it) {
        die "Error: The input files are not correct. Use sorted by name BAM file";
    }
    if (!(exists $hcount{$read1})) {
        foreach $k (keys %hcount) {
            foreach $k1 (keys %{$hcount{$k}}) {
               $total_reads++;                                   # Increment Total Mapped Read Count
               if ($hcount{$k}{$k1}[0] == 1) {$uniq_reads++;}    # Increment Unique Mapped Read Count
               if ($hcount{$k}{$k1}[1] >= 1) {$t_genic++;}       # Increment Genic Mapped Read Count
               elsif ($hcount{$k}{$k1}[2] >= 1) {$t_inter++;}    # Increment Intergenic Mapped Read Count
           }
        }
        %hcount = ();
    }
    
    $hcount{$read1}{$tag} = [0,0,0] if (!(exists $hcount{$read1}{$tag}));

    $g_count = (split (/\t/,$_))[6];             # Read genic map count for read
    $it_count = (split (/\t/,$it_line))[6];      # Read intergenic map count for read
    
    if ($g_count >= 1) {
	$hcount{$read1}{$tag}[1]++;             # Increment Genic Count
    }
    elsif ($it_count >= 1) {
	$hcount{$read1}{$tag}[2]++;             # Increment Intergenic Count
    }
    else {
        print "Unmapped read encountered. Error, if input BAM file was from TopHat. May be encountered in Bowtie output\n";
    }
    
    # Total number of maps for each read
    $hcount{$read1}{$tag}[0] = $hcount{$read1}{$tag}[1] + $hcount{$read1}{$tag}[2];
}

###For the last read in the hash.
foreach $k (keys %hcount) {
    foreach $k1 (keys %{$hcount{$k}}) {
        $total_reads++;                                   # Increment Total Mapped Read Count
        if ($hcount{$k}{$k1}[0] == 1) {$uniq_reads++;}    # Increment Unique Mapped Read Count
        if ($hcount{$k}{$k1}[1] >= 1) {$t_genic++;}       # Increment Genic Mapped Read Count
        elsif ($hcount{$k}{$k1}[2] >= 1) {$t_inter++;}    # Increment Intergenic Mapped Read Count
    }
}
    

close $genic_file;
close $inter_file;


$p_genic = sprintf("%.2f",eval(($t_genic/$total_reads)*100));
$p_inter = sprintf("%.2f",eval(($t_inter/$total_reads)*100));
print $fout "\#Total reads mapped:\t$total_reads\n";
print $fout "\#Uniquely mapped reads:\t$uniq_reads\n";
if ($out_flag == 1) {
    print $fout "\#Exon\tIntron\tIntergenic\n";
    print $fout "$p_genic\t0\t$p_inter\n";
} 
else {
    print $fout "\#Genic\tIntergenic\n";
    print $fout "$p_genic\t$p_inter\n";
}
close $fout;

}

###Percent reads mapping to exon, intron and intergenic region for eukaryotes.
sub mapped_euk {

my  $exon_file = shift;
my  $intron_file = shift;
my  $inter_file = shift;
my  $fout = shift;
my  ($e_count,$in_count,$it_count,$p_exon,$p_intron,$p_inter,$in_line,$it_line);
my  ($k,$k1,$tag,$read1,$read_e,$read_in,$read_it);
my  $t_exon = 0;
my  $t_intron = 0;
my  $t_inter = 0;  
my  $uniq_reads = 0;
my  $total_reads = 0;
my  $t_count = 0;
my  %hcount = ();

while (<$exon_file>) {
    chomp($_);
    $read_e = (split (/\t/,$_))[3];
    $read1 = (split (/\//,$read_e))[0];
    $tag = (split (/\//,$read_e))[1];
    if (! (defined $tag)){ $tag = 0;}
    $in_line = <$intron_file>;
    chomp($in_line);
    $read_in = (split (/\t/,$in_line))[3];
    $it_line = <$inter_file>;
    chomp($it_line);
    $read_it = (split (/\t/,$it_line))[3];
    ##Check if reads in exon,intron and intergenic file is same.
    if ($read_e ne $read_in or $read_e ne $read_it) {
        die "Error: The input files are not correct. Use sorted by name bam file";
    }

    if (!(exists $hcount{$read1})) {
        foreach $k (keys %hcount) {
            foreach $k1 (keys %{$hcount{$k}}) {
               $total_reads++;                                   # Increment Total Mapped Read Count
               if ($hcount{$k}{$k1}[0] == 1) {$uniq_reads++;}    # Increment Unique Mapped Read Count
               if ($hcount{$k}{$k1}[1] >= 1) {$t_exon++;}        # Increment Exonic Mapped Read Count
               elsif ($hcount{$k}{$k1}[2] >= 1) {$t_intron++;}   # Increment Intronic Mapped Read Count
               elsif ($hcount{$k}{$k1}[3] >= 1) {$t_inter++;}    # Increment Intergenic Mapped Read Count
           }
        }
        %hcount = ();
    }
    
    $hcount{$read1}{$tag} = [0,0,0,0] if (!(exists $hcount{$read1}{$tag}));

    $e_count = (split (/\t/,$_))[6];             # Read exonic map count for read
    $in_count = (split (/\t/,$in_line))[6];      # Read intronic map count for read
    $it_count = (split (/\t/,$it_line))[6];      # Read intergenic map count for read
    
    if ($e_count >= 1) {
	$hcount{$read1}{$tag}[1]++;             # Increment Exonic Count
    }
    elsif ($in_count >= 1) {
	$hcount{$read1}{$tag}[2]++;             # Increment Intronic Count
    }
    elsif ($it_count >= 1) {
	$hcount{$read1}{$tag}[3]++;             # Increment Intergenic Count
    }
    else {
        print "Unmapped read encountered. Error, if input BAM file was from TopHat. May be encountered in Bowtie output\n";
    }
    
    # Total number of maps for each read
    $hcount{$read1}{$tag}[0] = $hcount{$read1}{$tag}[1] + $hcount{$read1}{$tag}[2] + $hcount{$read1}{$tag}[3];
}

###For the last read in the hash.
foreach $k (keys %hcount) {
    foreach $k1 (keys %{$hcount{$k}}) {
        $total_reads++;                                   # Increment Total Mapped Read Count
        if ($hcount{$k}{$k1}[0] == 1) {$uniq_reads++;}    # Increment Unique Mapped Read Count
        if ($hcount{$k}{$k1}[1] >= 1) {$t_exon++;}        # Increment Exonic Mapped Read Count
        elsif ($hcount{$k}{$k1}[2] >= 1) {$t_intron++;}   # Increment Intronic Mapped Read Count
        elsif ($hcount{$k}{$k1}[3] >= 1) {$t_inter++;}    # Increment Intergenic Mapped Read Count
    }
}
close $exon_file;
close $intron_file;
close $inter_file;


$p_exon = sprintf("%.2f",eval(($t_exon/$total_reads)*100));
$p_intron = sprintf("%.2f",eval(($t_intron/$total_reads)*100));
$p_inter = sprintf("%.2f",eval(($t_inter/$total_reads)*100));
print $fout "\#Total reads mapped:\t$total_reads\n";
print $fout "\#Uniquely mapped reads:\t$uniq_reads\n";
print $fout "\#Exon\tIntron\tIntergenic\n";
print $fout "$p_exon\t$p_intron\t$p_inter\n";

close $fout;
}

###Bed Intersect
sub bed_intersect {
    my $phCmdLineOption	= shift;
    my $bam_file = shift;
    my $bed_file = shift;
    my $outfile = shift;
    my $sCmd;
    $sCmd = $phCmdLineOption->{'bedtools_bin_dir'}."/bedtools intersect".
			" -abam ".$bam_file.
			" -b ".$bed_file.
			" -bed -c ".
			" > ".$outfile;
   
   exec_command($sCmd); 
}




# exec_command()
#
# Purpose
#   executes shell command
#
# Required Parameters
#   sCmd			= /path/to/output_directory
#   
# Optional Parameters
#   none
#
# Returns
#   nothing
#
# Side Effects
#   none
#
# Assumptions
#   none
#
# Notes
#
sub exec_command {
	my $sCmd = shift;
	
	if ((!(defined $sCmd)) || ($sCmd eq "")) {
		die "\nSubroutine::exec_command : ERROR! Incorrect command!\n";
	}
	
	my $nExitCode;
	
	print STDERR "\n$sCmd\n";
	$nExitCode = system("$sCmd");
	if ($nExitCode != 0) {
		die "\tERROR! Command Failed!\n\t$!\n";
	}
	print STDERR "\n";
	
	return;
}



# Process_Annotation()
#
# Purpose
#   generates a sorted BED file for the user specified annotation file
#
# Required Parameters
#   phCmdLineOption		= pointer to hash containing command line options
#   sAnnotationFile		= path to annotation file
#   sExtension			= file extension for annotation file
#   sRegion				= region tag to be added to the output file
#   sOutDir				= output directory
#
# Optional Parameters
#   none
#
# Returns
#   sBedFile			= sorted BED format file
#
# Side Effects
#   none
#
# Assumptions
#
# Notes
#
sub Process_Annotation {
    my $phCmdLineOption		= shift;
    my $sAnnotationFile		= shift;
    my $sExtension			= shift;
    my $sRegion				= shift;
    my $sOutDir				= shift;

    my $sSubName = (caller(0))[3];

    if (! ((defined $sAnnotationFile) &&
           (defined $sExtension) &&
           (defined $sRegion) &&
           (defined $sOutDir) &&
           (defined $phCmdLineOption))) {
        croak "$sSubName - required parameters missing\n";
    }

    # Local variables
    my $bDebug   = (defined $phCmdLineOption->{'debug'}) ? TRUE : FALSE;
    my $bVerbose = (defined $phCmdLineOption->{'verbose'}) ? TRUE : FALSE;
    
    my ($sBedFile, $sSortedFile);

    # Start

#    ($bDebug) ? print STDERR "In $sSubName\n" : ();

    ($bDebug || $bVerbose) ? print STDERR "\n\tProcessing annotation file $sAnnotationFile ...\n" : ();
	
    if (($sExtension =~ m/^gtf$/i) || 
		($sExtension =~ m/^gff3$/i)) {
		$sBedFile = Process_GFF($phCmdLineOption, $sAnnotationFile, $sRegion,
								$phCmdLineOption->{'feature'}, $phCmdLineOption->{'attribute'}, $sOutDir);
	}
	else {
		$sBedFile = Init_OutFileName($phCmdLineOption, $sOutDir, $sAnnotationFile, ".bed");
    	$sBedFile .= ".$sRegion.bed";
    	
    	$sCmd = "ln -s $sAnnotationFile $sBedFile";
		exec_command($sCmd);
	}
	
    $sSortedFile = Init_OutFileName($phCmdLineOption, $sOutDir, $sBedFile, ".bed");
    $sSortedFile .= '.sorted.bed';
	
	$sCmd = $phCmdLineOption->{'bedtools_bin_dir'}."/bedtools sort".
			" -i ".$sBedFile." > ".$sSortedFile;
	
	exec_command($sCmd);
	
    ($bDebug || $bVerbose) ? print STDERR "\n" : ();

#    ($bDebug) ? print STDERR "Leaving $sSubName\n" : ();
    
    return $sSortedFile;
}

# Generate_Gene_BedFile()
#
# Purpose
#   generates a sorted BED file for the genes from the user specified annotation file
#
# Required Parameters
#   phCmdLineOption		= pointer to hash containing command line options
#   sAnnotationFile		= path to annotation file
#   sExtension			= file extension for annotation file
#   sRegion				= region tag to be added to the output file
#   sOutDir				= output directory
#
# Optional Parameters
#   none
#
# Returns
#   sBedFile			= sorted BED format file
#
# Side Effects
#   none
#
# Assumptions
#
# Notes
#
sub Generate_Gene_BedFile {
    my $phCmdLineOption		= shift;
    my $sAnnotationFile		= shift;
    my $sExtension			= shift;
    my $sRegion				= shift;
    my $sOutDir				= shift;

    my $sSubName = (caller(0))[3];

    if (! ((defined $sAnnotationFile) &&
           (defined $sExtension) &&
           (defined $sRegion) &&
           (defined $sOutDir) &&
           (defined $phCmdLineOption))) {
        croak "$sSubName - required parameters missing\n";
    }

    # Local variables
    my $bDebug   = (defined $phCmdLineOption->{'debug'}) ? TRUE : FALSE;
    my $bVerbose = (defined $phCmdLineOption->{'verbose'}) ? TRUE : FALSE;
    
    my ($sBedFile, $sGroupedFile, $sSortedFile, $sTmpFile);
    my ($sID, $sStrand, $sRefID, $nStart, $nEnd);
    my ($fpBED, $fpTMP);

    # Start

#    ($bDebug) ? print STDERR "In $sSubName\n" : ();

    ($bDebug || $bVerbose) ? print STDERR "\n\tProcessing annotation file $sAnnotationFile ...\n" : ();
	
    if (($sExtension =~ m/^gtf$/i) || 
		($sExtension =~ m/^gff3$/i)) {
		$sBedFile = Process_GFF($phCmdLineOption, $sAnnotationFile, "groupby",
					$phCmdLineOption->{'feature'}, $phCmdLineOption->{'groupby'}, $sOutDir);
	}
	else {
		$sBedFile = Init_OutFileName($phCmdLineOption, $sOutDir, $sAnnotationFile, ".bed");
    	$sBedFile .= ".groupby.bed";
    	
    	$sCmd = "ln -s $sAnnotationFile $sBedFile";
		exec_command($sCmd);
	}
	
	$sTmpFile = Init_OutFileName($phCmdLineOption, $sOutDir, $sBedFile, ".groupby.bed");
    $sTmpFile .= ".sortbygroup.bed";
    
    $sCmd = "sort -k 4,4 -k 6,6 $sBedFile > $sTmpFile";
    exec_command($sCmd);
    
    $sCmd = "mv $sTmpFile $sBedFile";
    exec_command($sCmd);
	
	$sGroupedFile = Init_OutFileName($phCmdLineOption, $sOutDir, $sBedFile, ".groupby.bed");
    $sGroupedFile .= ".$sRegion.bed";
	
	$sCmd = $phCmdLineOption->{'bedtools_bin_dir'}."/bedtools groupby".
			" -i ".$sBedFile.
			" -g 1,4,6".
			" -opCols 2,3".
			" -ops min,max".
			" > ".$sGroupedFile;
	
	exec_command($sCmd);
	
	$sTmpFile = Init_OutFileName($phCmdLineOption, $sOutDir, $sBedFile, ".groupby.bed");
    $sTmpFile .= ".reorder.bed";
    
    open ($fpBED, "$sGroupedFile") or die "\tERROR! Cannot open $sGroupedFile for reading!\n";
    open ($fpTMP, ">$sTmpFile") or die "ERROR! Cannot open $sTmpFile for writing!\n";
    
    while (<$fpBED>) {
    	$_ =~ s/\s+$//;
    	
		($sRefID, $sID, $sStrand, $nStart, $nEnd) = split(/\t/, $_);
		
		print $fpTMP "$sRefID\t$nStart\t$nEnd\t$sID\t.\t$sStrand\n";
    }
    
    close($fpBED);
    close($fpTMP);
    
    $sCmd = "mv $sTmpFile $sGroupedFile";
    exec_command($sCmd);
	
    $sSortedFile = Init_OutFileName($phCmdLineOption, $sOutDir, $sGroupedFile, ".bed");
    $sSortedFile .= '.sorted.bed';
	
	$sCmd = $phCmdLineOption->{'bedtools_bin_dir'}."/bedtools sort".
			" -i ".$sGroupedFile." > ".$sSortedFile;
	
	exec_command($sCmd);
	
    ($bDebug || $bVerbose) ? print STDERR "\n" : ();

#    ($bDebug) ? print STDERR "Leaving $sSubName\n" : ();
    
    return $sSortedFile;
}

# Process_GFF()
#
# Purpose
#   generates a BED file for the user specified GTF or GFF annotation file
#
# Required Parameters
#   phCmdLineOption		= pointer to hash containing command line options
#   sAnnotationFile		= path to annotation file
#   sRegion				= region tag to be added to the output file
#   sFeatureID			= ID to be filtered for from feature (column 3) of the GFF file
#   sAttributeID		= ID to be extracted from the attributes (column 9) of the GFF file
#   sOutDir				= output directory
#
# Optional Parameters
#   none
#
# Returns
#   sBedFile			= BED format file
#
# Side Effects
#   none
#
# Assumptions
#
# Notes
#
sub Process_GFF {
    my $phCmdLineOption		= shift;
    my $sAnnotationFile		= shift;
    my $sRegion				= shift;
    my $sFeatureID			= shift;
    my $sAttributeID		= shift;
    my $sOutDir				= shift;

    my $sSubName = (caller(0))[3];

    if (! ((defined $sAnnotationFile) &&
           (defined $sRegion) &&
           (defined $sFeatureID) &&
           (defined $sAttributeID) &&
           (defined $sOutDir) &&
           (defined $phCmdLineOption))) {
        croak "$sSubName - required parameters missing\n";
    }

    # Local variables
    my $bDebug   = (defined $phCmdLineOption->{'debug'}) ? TRUE : FALSE;
    my $bVerbose = (defined $phCmdLineOption->{'verbose'}) ? TRUE : FALSE;
    my %att=();
    my (@atribs);
    my ($sBedFile);
    my ($sRefID, $sSource, $sFeature, $nStart, $nEnd, $sScore, $sStrand, $sFrame, $sAttributes, $sID,@temp);
    my ($fpGFF, $fpBED);
    my ($nI);

    # Start

#    ($bDebug) ? print STDERR "In $sSubName\n" : ();

    ($bDebug || $bVerbose) ? print STDERR "\n\tProcessing GTF/GFF file $sAnnotationFile ...\n" : ();
    
    # initialize output file spec
    $sBedFile = Init_OutFileName($phCmdLineOption, $sOutDir, $sAnnotationFile, ".".lc($phCmdLineOption->{'annotationfiletype'}));
    $sBedFile .= ".$sRegion.bed";
    
    open($fpGFF, "$sAnnotationFile") or die "\tERROR! Cannot open $sAnnotationFile for reading !!!\n";
    open($fpBED, ">$sBedFile") or die "\tERROR! Cannot open $sBedFile for writing !!!\n";
    
    $nI = 0;
    while (<$fpGFF>) {
    	$_ =~ s/\s+$//;
    	
    	next if ($_ =~ m/^#/);
    	
    	($sRefID, $sSource, $sFeature, $nStart, $nEnd, $sScore, $sStrand, $sFrame, $sAttributes) = split(/\t/, $_);
    	
    	next if ($sFeature ne $sFeatureID);
    	
    	next if ($sStrand !~ m/[-+]/);
    	
    	$sID = "Unknown";
        @atribs = split(/;/,$sAttributes);
        foreach (@atribs) {
	    if ($_=~/\=/) {
                $att{(split (/\=/,$_)) [0]} = (split (/\=/,$_))[1];
            }
            elsif ($_=~/\s+/){
                $_=~s/\s+//g;
                @temp=split (/\"/,$_);
         	$att{$temp[0]} = $temp[1];
            }
            else {die "Unrecognized annotation file";} 
        }
 
        if (exists $att{$sAttributeID}) {
            $sID = $att{$sAttributeID};
        } 
        else {
    		$nI++;
    		$sID .= ".$nI";
    	}

    	print $fpBED $sRefID."\t".($nStart - 1)."\t". $nEnd."\t".$sID;
    	print $fpBED "\t.\t".$sStrand;
    	print $fpBED "\n";
    }
    
    close($fpGFF);
    close($fpBED);

 
    ($bDebug || $bVerbose) ? print STDERR "\n" : ();

#    ($bDebug) ? print STDERR "Leaving $sSubName\n" : ();
    
    return $sBedFile;
}


# Init_OutFileName()
#
# Purpose
#   generates a new filename based on existing filename by removing the
#   extension, and prepending an output directory
#
# Required Parameters
#   phCmdLineOption = pointer to hash containing command line options
#   sOutDir         = output directory
#   sFileName       = existing filename
#   sExtension      = filename extension
#
# Optional Parameters
#   none
#
# Returns
#   sOutFile        = new filename
#
# Side Effects
#   none
#
# Assumptions
#
# Notes
#
sub Init_OutFileName {
    my $phCmdLineOption = shift;
    my $sOutDir         = shift;
    my $sFileName       = shift;
    my $sExtension      = shift;

    my $sSubName = (caller(0))[3];

    if (! ((defined $sFileName) &&
           (defined $sExtension) &&
           (defined $sOutDir) &&
           (defined $phCmdLineOption))) {
        croak "$sSubName - required parameters missing\n";
    }

    # Local variables
    my $bDebug   = (defined $phCmdLineOption->{'debug'}) ? TRUE : FALSE;
    my $bVerbose = (defined $phCmdLineOption->{'verbose'}) ? TRUE : FALSE;
 
    my ($sFile, $sDir);
    my $sOutFile;

    # Start

#    ($bDebug) ? print STDERR "In $sSubName\n" : ();

#    ($bDebug || $bVerbose) ? print STDERR "\n\tInitializing filename...\n" : ();

    ($_, $sDir, $sFile) = File::Spec->splitpath($sFileName);
    $sFile =~ m/^(\S+)$sExtension$/;
	if (defined $1) {
	    $sOutFile = $sOutDir.'/'.$1;
	}
	else {
	    $sOutFile = $sOutDir.'/'.$sFile;
	}
    $sOutFile = File::Spec->canonpath($sOutFile);

#    ($bDebug || $bVerbose) ? print STDERR "\n" : ();

#    ($bDebug) ? print STDERR "Leaving $sSubName\n" : ();
    
    return $sOutFile;
}


##############################################################################
### POD Documentation
##############################################################################

__END__

=head1 NAME

percent_mapped_stats.pl                     -program to generate percentage of reads mapping to exon/intron/intergenic region in euk
					     and genic/intergenic region in prok.

=head1 SYNOPSIS

    percent_mapped_stats.pl                --r <reference_file> --i <bam_file> --a <annotation_file> 
                                           --t <annotation file format(bed|gtf|gff3)> --f <feature> --id <attribute> 
                                           --g <group_by> --o [output_dir] --org [org-type] 
                                           --b [bedtools_bin_dir] --s [samtools_bin_dir] [--v]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

   
    --r <reference_file>                            = path to reference fastA file.

    --i <bam_file>                                  = path to alignment BAM file SORTED BY NAME

    --a <annotation_file>                           = path to annotation file (BED or GTF or GFF3 format file)

    --t <annotation_file_format>                    = annotation file format (bed/gtf/gff3)

    --f <feature>                                   = feature type from column 3 of GTF or GFF3 file

    --id <attribute>                                = attribute id from column 9 of GTF or GFF3 file to be used as region ID

    --g <group_by>                                  = group_by id from column 9 of GTF or GFF3 file to be used to group regions by
                                                      use 'yes' if annotation file format is 'BED'. Required for eukaryotes.

    --org <org-type>                                = Organism type(prok/euk)

    --o [output_dir]                                = output directory. Optional

    --b [bedtools_bin_dir]                          = bedtools binary directory. Optional

    --s [samtools_bin_dir]                          = samtools binary directory. Optional


    --v                                             = generate runtime messages. Optional

=head1 DESCRIPTION

 This script generates percent of reads mapping to Exonic, Intronic and Intergenic region in the genome during reference mapping for either prokaryotes or eukaryotes.

=head1 AUTHOR

 Priti Kumari
 Bioinformatics Software Engineer I
 Institute for Genome Sciences
 University of Maryland
 Baltimore, Maryland 21201

=cut
