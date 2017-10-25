#!/usr/bin/env perl -w
##############################################################################
### This program generates coverage statistics from the alignment BAM file 
### across the genomic, genic, exonic, intronic, and/or intergenic regions
##############################################################################

use strict;
use Cwd;
use Carp;
use Getopt::Long;
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
my ($sGenicBedFile, $sExonicBedFile, $sIntronicBedFile, $sIntergenicBedFile);
my ($sOutDir, $sSizeFile, $sOutFile, $sSortedFile, $sSummaryFile);
my ($sPrefix, $sSampleId, $sRegion, $nTotalMappedReads, $nUniqueMappedReads);
my ($sCmd, $nNumCols);
my ($bDebug, $bVerbose);

##############################################################################
### Main
##############################################################################

GetOptions( \%hCmdLineOption,
            'reffile|r=s', 'infile|i=s', 'annotation|a=s', 'regiontype|rt=s', 
            'annotationfiletype|t=s', 'feature|f=s', 'attribute|id=s', 'groupby|g=s', 
            'outdir|o=s', 'bedtools_bin_dir|b=s', 'samtools_bin_dir|s=s', 
            'verbose|v',
            'debug',
            'help',
            'man') or pod2usage(2);

if ($hCmdLineOption{'help'} || 
	(! defined $hCmdLineOption{'infile'}) || 
	(! defined $hCmdLineOption{'annotation'}) || 
	(! defined $hCmdLineOption{'annotationfiletype'}) || 
	(! defined $hCmdLineOption{'regiontype'})) {
    pod2usage( -msg => $sHelpHeader, -exitval => 1);
}

if (($hCmdLineOption{'annotationfiletype'} !~ m/^bed$/i) && 
	($hCmdLineOption{'annotationfiletype'} !~ m/^gtf$/i) &&
	($hCmdLineOption{'annotationfiletype'} !~ m/^gff3$/i)) {
    die "\tERROR: Annotation file format needs to be a BED, GTF or GFF3 format file\n";
}

if ($hCmdLineOption{'groupby'} eq 'NONE' ){undef $hCmdLineOption{'groupby'};}

if (($hCmdLineOption{'annotationfiletype'} =~ m/^(gtf|gff3)$/i) && 
	((! defined $hCmdLineOption{'feature'}) || (! defined $hCmdLineOption{'attribute'}))) {
    die "\tERROR: 'feature' and 'attribute' required for GTF and GFF3 annotation file formats\n";
}

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

$hCmdLineOption{'infile'} = File::Spec->rel2abs($hCmdLineOption{'infile'});
($_, $_, $sBamFile) = File::Spec->splitpath($hCmdLineOption{'infile'});

$sSampleId = $sBamFile;
$sSampleId =~ s/\.bam$//;

# calculate total number of mapped reads

($bDebug || $bVerbose) ? 
	print STDERR "\nCalculating Total Mapped Reads .....\n" : ();

$sCmd = $hCmdLineOption{'samtools_bin_dir'}."/samtools view".
		" -c -F 4".
		" ".$hCmdLineOption{'infile'};

$nTotalMappedReads = `$sCmd`;
chomp($nTotalMappedReads);

($bDebug || $bVerbose) ? 
	print STDERR "\tTotal Mapped Reads : $nTotalMappedReads\n" : ();

($bDebug || $bVerbose) ? 
	print STDERR "Calculating Total Mapped Reads ..... done\n" : ();

# calculate unique number of mapped reads

($bDebug || $bVerbose) ? 
	print STDERR "\nCalculating Unique Mapped Reads .....\n" : ();

$sCmd = $hCmdLineOption{'samtools_bin_dir'}."/samtools view".
		" -c -F 260".
		" ".$hCmdLineOption{'infile'};

$nUniqueMappedReads = `$sCmd`;
chomp($nUniqueMappedReads);

($bDebug || $bVerbose) ? 
	print STDERR "\tUnique Mapped Reads : $nUniqueMappedReads\n" : ();

($bDebug || $bVerbose) ? 
	print STDERR "Calculating Unique Mapped Reads ..... done\n" : ();

@aRegions = split(/[:,;]/, $hCmdLineOption{'regiontype'});

die "\nERROR! No regions defined for rpkm calculations!\n" if (@aRegions < 1);

@aRegions = sort keys %{{ map { $_ => 1 } @aRegions }};

($bDebug || $bVerbose) ? 
	print STDERR "\nGenerating coverage statistics for the following regions :\n@aRegions\n" : ();

foreach $sRegion (@aRegions) {
	
	if ($sRegion eq "genomic") {
		# process reference file
		
		die "Error! Reference Fasta file not found!\n" if (! (defined $hCmdLineOption{'reffile'}) );
		
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
		if ($hCmdLineOption{'regiontype'} =~ m/genomic/) {
			$sPrefix = $sRefFile;
			$sPrefix =~ s/\.([A-za-z])+$//;
			
			$sSizeFile = File::Spec->rel2abs("$sOutDir/$sPrefix.chromosome.sizes");
			
			($bDebug || $bVerbose) ? 
				print STDERR "\nGenerating $sSizeFile .....\n" : ();
			
			$sCmd = "cut -f1,2 $sOutDir/$sRefFile.fai > $sSizeFile";
			
			exec_command($sCmd);
			
			($bDebug || $bVerbose) ? 
				print STDERR "Generating $sSizeFile ..... done\n" : ();
		}
		
		# genomic coverage statistics
		
		if (! -e "$sOutDir/genomic_coverage") {
			mkdir "$sOutDir/genomic_coverage" or die "\tError : Unable to create directory '$sOutDir/genomic_coverage' .....\n";
		}
		
		($bDebug || $bVerbose) ? 
			print STDERR "\nGenerating genome coverage statistics for $sBamFile .....\n" : ();
		
		$sOutFile = File::Spec->rel2abs("$sOutDir/genomic_coverage/$sSampleId.genomic.coverage.stats.txt");
		$sSummaryFile = File::Spec->rel2abs("$sOutDir/genomic_coverage/$sSampleId.genomic.coverage.summary.txt");
		Genome_Coverage(\%hCmdLineOption, $hCmdLineOption{'infile'}, $sSizeFile, $sOutFile, $sSummaryFile);
		
		($bDebug || $bVerbose) ? 
			print STDERR "Generating genome coverage statistics for $sBamFile ..... done\n" : ();
	}
	
	if ($sRegion eq "genic") {
		# genic coverage statistics
		
		if (defined $hCmdLineOption{'groupby'}) {
			$sGenicBedFile = Generate_Gene_BedFile(\%hCmdLineOption, $hCmdLineOption{'annotation'}, $hCmdLineOption{'annotationfiletype'}, "genic", $sOutDir);
		}
		else {
			$sGenicBedFile = Process_Annotation(\%hCmdLineOption, $hCmdLineOption{'annotation'}, $hCmdLineOption{'annotationfiletype'}, "genic", $sOutDir);
		}
		
		if (! -e "$sOutDir/genic_coverage") {
			mkdir "$sOutDir/genic_coverage" or die "\tError : Unable to create directory '$sOutDir/genic_coverage' .....\n";
		}
		
		($bDebug || $bVerbose) ? 
			print STDERR "\nGenerating gene coverage statistics for $sBamFile .....\n" : ();
		
		$sOutFile = File::Spec->rel2abs("$sOutDir/genic_coverage/$sSampleId.genic.coverage.stats.txt");
		Feature_Coverage(\%hCmdLineOption, $hCmdLineOption{'infile'}, $sGenicBedFile, $sOutFile, $nTotalMappedReads, $nUniqueMappedReads);
		
		($bDebug || $bVerbose) ? 
			print STDERR "Generating gene coverage statistics for $sBamFile ..... done\n" : ();
	}
	
	if ($sRegion eq "exonic") {
		# exonic coverage statistics
		
		$sSortedFile = Process_Annotation(\%hCmdLineOption, $hCmdLineOption{'annotation'}, $hCmdLineOption{'annotationfiletype'}, "exonic", $sOutDir);
		
		$sExonicBedFile = Init_OutFileName(\%hCmdLineOption, $sOutDir, $sSortedFile, ".bed");
	    $sExonicBedFile .= '.sanitized.bed';
		
		$sCmd = $hCmdLineOption{'bedtools_bin_dir'}."/bedtools merge".
				" -i ".$sSortedFile." > ".$sExonicBedFile;
		
		exec_command($sCmd);
		
		if (! -e "$sOutDir/exonic_coverage") {
			mkdir "$sOutDir/exonic_coverage" or die "\tError : Unable to create directory '$sOutDir/exonic_coverage' .....\n";
		}
		
		($bDebug || $bVerbose) ? 
			print STDERR "\nGenerating exon coverage statistics for $sBamFile .....\n" : ();
		
		$sOutFile = File::Spec->rel2abs("$sOutDir/exonic_coverage/$sSampleId.exonic.coverage.stats.txt");
		Feature_Coverage(\%hCmdLineOption, $hCmdLineOption{'infile'}, $sExonicBedFile, $sOutFile, $nTotalMappedReads, $nUniqueMappedReads);
		
		($bDebug || $bVerbose) ? 
			print STDERR "Generating exon coverage statistics for $sBamFile ..... done\n" : ();
	}
}

($bDebug || $bVerbose) ? 
	print "Removing BED files........\n" : ();

unlink glob "$sOutDir/*bed";

($bDebug || $bVerbose) ? 
	print STDERR "\nGenerating coverage statistics for the specified regions ..... done\n" : ();

exit;


##############################################################################
### Subroutines
##############################################################################


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
    
    $sCmd = "sort -k 4,4 $sBedFile > $sTmpFile";
    exec_command($sCmd);
    
    $sCmd = "mv $sTmpFile $sBedFile";
    exec_command($sCmd);
	
	$sGroupedFile = Init_OutFileName($phCmdLineOption, $sOutDir, $sBedFile, ".groupby.bed");
    $sGroupedFile .= ".$sRegion.bed";
	
	$sCmd = $phCmdLineOption->{'bedtools_bin_dir'}."/bedtools groupby".
			" -i ".$sBedFile.
			" -g 1,4".
			" -opCols 2,3,6".
			" -ops min,max,distinct".
			" > ".$sGroupedFile;
	
	exec_command($sCmd);
	
	$sTmpFile = Init_OutFileName($phCmdLineOption, $sOutDir, $sBedFile, ".groupby.bed");
    $sTmpFile .= ".reorder.bed";
    
    open ($fpBED, "$sGroupedFile") or die "\tERROR! Cannot open $sGroupedFile for reading!\n";
    open ($fpTMP, ">$sTmpFile") or die "ERROR! Cannot open $sTmpFile for writing!\n";
    
    while (<$fpBED>) {
    	$_ =~ s/\s+$//;
    	
    	($sRefID, $sID, $nStart, $nEnd, $sStrand) = split(/\t/, $_);
    	
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
    
    my ($sBedFile);
    my ($sRefID, $sSource, $sFeature, $nStart, $nEnd, $sScore, $sStrand, $sFrame, $sAttributes, $sID);
    my (@aAttr, $sAttr, $sKey, $sValue);
    my ($fpGFF, $fpBED);
    my ($nI);
    my ($bFlag);

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
    	
    	$bFlag = 0;
    	$sID = "Unknown";
    	@aAttr = split(/;/, $sAttributes);
    	
    	foreach $sAttr (@aAttr) {
    		$sAttr =~ s/^\s+//;
    		$sAttr =~ s/\s+$//;
    		($sKey, $sValue) = split(/[\s\=]/, $sAttr);
    		
    		if ($sKey eq $sAttributeID) {
    			$sID = $sValue;
    			$sID =~ s/\"//g;
    			$bFlag = 1;
    			last;
    		}
    	}
    	
    	if (!($bFlag)) {
    		$nI++;
    		$sID .= ".$nI";
    	}
    	
    	print $fpBED $sRefID."\t".($nStart - 1)."\t".$nEnd."\t".$sID."\t.\t".$sStrand."\n";
    }
    
    close($fpGFF);
    close($fpBED);
    
    ($bDebug || $bVerbose) ? print STDERR "\n" : ();

#    ($bDebug) ? print STDERR "Leaving $sSubName\n" : ();
    
    return $sBedFile;
}

# Add_RPKMs()
#
# Purpose
#   addition of RPKM values to the coverage statistics file
#
# Required Parameters
#   phCmdLineOption		= pointer to hash containing command line options
#   sStatFile			= path to coverage statistics file
#   nTotalMappedReads	= total number of mapped reads
#   nUniqueMappedReads	= unique number of mapped reads
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
sub Add_RPKMs {
    my $phCmdLineOption		= shift;
    my $sStatFile			= shift;
    my $nTotalMappedReads	= shift;
    my $nUniqueMappedReads	= shift;
    my $sOutDir				= shift;

    my $sSubName = (caller(0))[3];

    if (! ((defined $sStatFile) &&
           (defined $nTotalMappedReads) &&
           (defined $nUniqueMappedReads) &&
           (defined $sOutDir) &&
           (defined $phCmdLineOption))) {
        croak "$sSubName - required parameters missing\n";
    }

    # Local variables
    my $bDebug   = (defined $phCmdLineOption->{'debug'}) ? TRUE : FALSE;
    my $bVerbose = (defined $phCmdLineOption->{'verbose'}) ? TRUE : FALSE;
    
    my ($sTmpFile);
    my (@aBED);
    my ($nRCount, $nLength, $nRPKM);
    my ($fpBED, $fpTMP);
    my ($nC, $nI, $nS, $nT);

    # Start

#    ($bDebug) ? print STDERR "In $sSubName\n" : ();

    ($bDebug || $bVerbose) ? print STDERR "\n\tProcessing statistics file $sStatFile ...\n" : ();
    
    # initialize output file spec
    $sTmpFile = Init_OutFileName($phCmdLineOption, $sOutDir, $sStatFile, ".txt");
    $sTmpFile .= ".rpkm.txt";
    
    open($fpBED, "$sStatFile") or die "\tERROR! Cannot open $sStatFile for reading !!!\n";
    open($fpTMP, ">$sTmpFile") or die "\tERROR! Cannot open $sTmpFile for writing !!!\n";
    
    $nC = $nI = $nS = $nT = 0;
    while (<$fpBED>) {
    	$_ =~ s/\s+$//;
    	
    	next if ($_ =~ m/^#/);
    	
    	@aBED = split(/\t/, $_);
    	
    	$nRCount = $aBED[-4];
    	$nLength = $aBED[2] - $aBED[1];
    	die "#".($nI + 1)." Erroneous feature $aBED[0]:$aBED[1]-$aBED[2] ($nLength <= 0)\n" if ($nLength <= 0);
    	$nRPKM = ($nRCount * (10**9)) / ($nTotalMappedReads * $nLength);
    	$nRPKM = sprintf("%.7f", $nRPKM);
    	push @aBED, $nRPKM;
    	
    	{
    		$" = "\t";
    		print $fpTMP "@aBED\n";
    	}
    	
    	$nI++;
    	$nT++ if ($nRPKM <= 0);
    	$nC += $nRCount;
    	$nS += $nRPKM;
    	
    	if ($nI%200000 == 0) {
    		($bDebug || $bVerbose) ? print STDERR "\r\t\t$nI features processed ..." : ();
    	}
    }
    
    ($bDebug || $bVerbose) ? print STDERR "\r\t\t$nI features processed ..." : ();
    
    print $fpTMP "\n# Number of features : ".$nI." (".eval($nI - $nT)." features completely or partially covered)\n".
    			 "# Total number of mapped reads : ".$nTotalMappedReads." (".$nUniqueMappedReads." unique mapped reads)\n".
    			 "# Number of reads mapped to features : ".$nC."\n".
    			 "# % reads mapped to features : ".sprintf("%.3f", eval(($nC / $nTotalMappedReads) * 100))."\n".
    			 "# Average RPKM across all features : ".sprintf("%.3f", eval($nS / $nI))." (".sprintf("%.3f", eval($nS / ($nI - $nT)))." across covered features)\n";
    
    close($fpBED);
    close($fpTMP);
    
    unlink $sStatFile if ( -e "$sTmpFile" );
    
    ($bDebug || $bVerbose) ? print STDERR "\n" : ();

#    ($bDebug) ? print STDERR "Leaving $sSubName\n" : ();
    
    return;
}

# Genome_Coverage()
#
# Purpose
#   generate genome coverage stats and summary files
#
# Required Parameters
#   phCmdLineOption		= pointer to hash containing command line options
#   sInFile				= path to alignment SAM/BAM file
#   sSizeFile			= path to genome size file
#   sOutFile			= path to output stats file
#   sSummaryFile		= path to output summary file
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
#
# Notes
#
sub Genome_Coverage {
    my $phCmdLineOption		= shift;
    my $sInFile				= shift;
    my $sSizeFile			= shift;
    my $sOutFile			= shift;
    my $sSummaryFile		= shift;

    my $sSubName = (caller(0))[3];

    if (! ((defined $sInFile) &&
           (defined $sSizeFile) &&
           (defined $sOutFile) &&
           (defined $sSummaryFile) &&
           (defined $phCmdLineOption))) {
        croak "$sSubName - required parameters missing\n";
    }

    # Local variables
    my $bDebug   = (defined $phCmdLineOption->{'debug'}) ? TRUE : FALSE;
    my $bVerbose = (defined $phCmdLineOption->{'verbose'}) ? TRUE : FALSE;
    
    my $sCmd;
    
    # Start

#    ($bDebug) ? print STDERR "In $sSubName\n" : ();
	
	$sCmd = $phCmdLineOption->{'bedtools_bin_dir'}."/bedtools genomecov".
			" -split".
			" -ibam ".$sInFile.
			" -g ".$sSizeFile.
			" > ".$sOutFile;
	
	exec_command($sCmd);
	
	$sCmd = "awk '\$2!=0' $sOutFile |".
			" ".$phCmdLineOption->{'bedtools_bin_dir'}."/bedtools groupby".
			" -g 1".
			" -opCols 3,4,5".
			" -ops sum,mean,sum".
			" > ".$sSummaryFile;
	
	exec_command($sCmd);

#    ($bDebug) ? print STDERR "Leaving $sSubName\n" : ();
    
    return;
}

# Feature_Coverage()
#
# Purpose
#   generate feature coverage stats and RPKM files
#
# Required Parameters
#   phCmdLineOption		= pointer to hash containing command line options
#   sInFile				= path to alignment SAM/BAM file
#   sBedFile			= path to feature coords file in BED format
#   sOutFile			= path to output stats file
#   nTotalMappedReads	= total number of mapped reads
#   nUniqueMappedReads	= total number of unique mapped reads
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
#
# Notes
#
sub Feature_Coverage {
    my $phCmdLineOption		= shift;
    my $sInFile				= shift;
    my $sBedFile			= shift;
    my $sOutFile			= shift;
    my $nTotalMappedReads	= shift;
    my $nUniqueMappedReads	= shift;

    my $sSubName = (caller(0))[3];

    if (! ((defined $sInFile) &&
           (defined $sBedFile) &&
           (defined $sOutFile) &&
           (defined $nTotalMappedReads) &&
           (defined $nUniqueMappedReads) &&
           (defined $phCmdLineOption))) {
        croak "$sSubName - required parameters missing\n";
    }

    # Local variables
    my $bDebug   = (defined $phCmdLineOption->{'debug'}) ? TRUE : FALSE;
    my $bVerbose = (defined $phCmdLineOption->{'verbose'}) ? TRUE : FALSE;
    
    my $sOutDir;
    my $sCmd;
    
    # Start

#    ($bDebug) ? print STDERR "In $sSubName\n" : ();
	
	$sCmd = $phCmdLineOption->{'bedtools_bin_dir'}."/bedtools coverage".
			" -split".
			" -abam ".$sInFile.
			" -b ".$sBedFile.
			" > ".$sOutFile;
	
	exec_command($sCmd);
	
	($bDebug || $bVerbose) ? 
		print STDERR "\n\tAdding RPKMs to coverage statistics file .....\n" : ();
	
	($_, $sOutDir, $_) = File::Spec->splitpath($sOutFile);
	
	Add_RPKMs(\%hCmdLineOption, $sOutFile, $nTotalMappedReads, $nUniqueMappedReads, $sOutDir);
	
#    ($bDebug) ? print STDERR "Leaving $sSubName\n" : ();
    
    return;
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

rpkm_coverage_stats.pl - program to generate read counts and coverage statistics across features.

=head1 SYNOPSIS

    rpkm_coverage_stats.pl      --r [reference_file] --i <bam_file> --a <annotation_file> --rt <region_type> 
                                --t <annotation file format (bed|gtf|gff3)> --f <feature> --id <attribute> 
                                --g <group_by> [--o <output_dir>] [--b <bedtools_bin_dir>] [--s <samtools_bin_dir>] 
                                [--v]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

   
    --r <reference_file>                            = path to reference fastA file. Required for 'Genomic Coverage'.

    --i <bam_file>                                  = path to alignment BAM file sorted by position

    --a <annotation_file>                           = path to annotation file (BED or GTF or GFF3 format file)

    --rt <region_type>                              = region to determine coverage for (genomic:genic:exonic) separated by ':' or ',' or ';'

    --t <annotation_file_format>                    = annotation file format (bed/gtf/gff3)

    --f <feature>                                   = feature type from column 3 of GTF or GFF3 file

    --id <attribute>                                = attribute id from column 9 of GTF or GFF3 file to be used as region ID

    --g <group_by>                                  = group_by id from column 9 of GTF or GFF3 file to be used to group regions by
                                                      use 'yes' if annotation file format is 'BED'

    --o <output_dir>                                = output directory. Optional

    --b <bedtools_bin_dir>                          = bedtools binary directory. Optional

    --s <samtools_bin_dir>                          = samtools binary directory. Optional

    --v                                             = generate runtime messages. Optional

=head1 DESCRIPTION

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

PERL5LIB environment variable should be set if Bio::Perl is in non-standard
location

=head1 DEPENDENCIES

Bio::Perl

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

There are no known bugs in this module. Please report problems to Amol Shetty
(ashetty@som.umaryland.edu). Patches are welcome.

=head1 AUTHOR

 Amol Carl Shetty
 Bioinformatics Software Engineer II
 Institute of Genome Sciences
 University of Maryland
 Baltimore, Maryland 21201

=head1 LICENSE AND COPYRIGHT

Copyright (c) 2010 Amol Carl Shetty (<ashetty@som.umaryland.edu>). All rights
reserved.

This program is free software; you can distribute it and/or modify it under
GNU-GPL licenscing terms.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or FITNESS
FOR A PARTICULAR PURPOSE. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=cut
