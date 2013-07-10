#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

################################################################################
### POD Documentation
################################################################################

=head1 NAME

create_euk_rnaseq_pipeline_config.pl - Creates the pipeline.layout and pipeline.config for the
                                       automated eukaryotic rna-seq pipeline

=head1 SYNOPSIS

    create_euk_rnaseq_pipeline_config.pl --s <samples_file> --c <config_file> --r <reference_fasta> 
                                         [--qual <quality_score_format>] [--gtf <annotation_file>] 
                                         [--bowtie_build] [--quality_stats] [--quality_trimming]
	                                 [--split][--alignment] [--bwtidxfile <bowtie_index>] [--visualization] 
                                         [--diff_gene_expr] [--comparison_groups <str>] [--count]  
                                         [--file_type <SAM|BAM>] [--sorted <position|name>] 
                                         [--isoform_analysis] [--include_novel] 
                                         [--diff_isoform_analysis] [--use_ref_gtf]
                                         [--td <template_directory>] [--o <outdir>] [--v] 
                                         [--man] [--help]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --s <samples_file>                = /path/to/samples file with information on all samples to be analyzed.

    --c <config_file>                 = /path/to/config file with parameter information for multiple components.

    --r <reference_fasta>             = /path/to/reference FastA file for all samples. Optional.
    
    --qual <quality_score_format>     = FastQ quality score format (33 or 64). Optional. [33]

    --gtf <annotation_file>           = /path/to/annotation file in GFF or GTF format.
                                        Required with '--diff_gene_expr, --alignment'.

    --annotation_format               = annotation file format (gtf/gff3). Required

    --bowtie_build                    = execute bowtie_build component. Requires '--r'.

    --quality_stats                   = execute fastx_quality_stats component.

    --quality_trimming                = execute fastx_trimming component. Also generates quality statistics.

    --split                           = excute fastq split component for shorter alignement time

      --bwtidxfile <bowtie_index>   = /path/to/bowtie index file for alignment of all samples.
                                        Required for '--alignment' if not specifying '--bowtie_build'.

   --alignment                       = execute tophat alignment component.
                                        Sample file should be in the following format
                                        #Sample_ID<tab>Group_ID<tab>Mate_Pair_1<tab>Mate_Pair_2

        --bwtidxfile <bowtie_index>   = /path/to/bowtie index file for alignment of all samples.
                                        Required for '--alignment' if not specifying '--bowtie_build'.

    --visualization                   = execute bam2bigwig component.
                                        Requires additional information in sample file if not specifying '--alignment'.
                                        Sample file should be in the following format
                                        #Sample_ID<tab>Group_ID<tab>Alignment_BAM_File

        --file_type <SAM|BAM>         = alignment file format (BAM or SAM). [BAM]
                                        Required if not specifying '--alignment' and providing alignment file information in sample file.

        --sorted <position>           = if alignment BAM/SAM file is already sorted by position. [undef]

    --rpkm_analysis                   = execute rpkm coverage analysis component.
                                        Requires additional information in sample file if not specifying '--alignment'.
                                        Sample file should be in the following format
                                        #Sample_ID<tab>Group_ID<tab>Alignment_BAM_File

        --file_type <SAM|BAM>         = alignment file format (BAM or SAM). [BAM]
                                        Required if not specifying '--alignment' and providing alignment file information in sample file.

        --sorted <position>           = if alignment BAM/SAM file is already sorted by position. [undef]

    --diff_gene_expr                  = execute differential gene expression analysis component.
                                        Requires additional information in sample file if not specifying '--alignment'.
                                        Sample file should be in the following format
                                        #Sample_ID<tab>Group_ID<tab>Alignment_File

        --comparison_groups <str>     = string of groups to compare. e.g. "GRP#2vsGRP#1,GRP#3vsGRP#1".
                                        Group Ids SHOULD match Group Names in column 2 of the sample info file

        --count                       = generate count files.
                                        Required if not specifying '--alignment' and providing alignment file information in sample file.

        --file_type <SAM|BAM>         = alignment file format (BAM or SAM).
                                        Required if specifying '--count'. [SAM]

        --sorted <name>               = if alignment SAM file is already sorted by name. [undef]

    --isoform_analysis                = execute isoform identification analysis component.
                                        Requires additional information in sample file if not specifying '--alignment'.
                                        Sample file should be in the following format
                                        #Sample_ID<tab>Group_ID<tab>Alignment_File

        --file_type <SAM|BAM>         = alignment file format (BAM or SAM). [BAM]
                                        Required if not specifying '--alignment' and providing alignment file information in sample file.

        --sorted <position>           = if alignment BAM/SAM file is already sorted by position. [undef]

        --include_novel               = will not use reference annotation to determine isoforms.

    --diff_isoform_analysis           = execute differential isoform expression analysis component.
                                        Requires additional information in sample file if not specifying '--alignment'.
                                        Sample file should be in the following format
                                        #Sample_ID<tab>Group_ID<tab>Alignment_File<tab>Cufflinks_GTF_File

        --comparison_groups <str>     = string of groups to compare. e.g. "GRP#2vsGRP#1,GRP#3vsGRP#1".
                                        Group Ids SHOULD match Group Names in column 2 of the sample info file

        --file_type <SAM|BAM>         = alignment file format (BAM or SAM). [BAM]
                                        Required if not specifying '--alignment' and providing alignment file information in sample file.

        --sorted <position>           = if alignment SAM file is already sorted by position. [undef]

        --use_ref_gtf                 = use reference gtf for with all cuffdiff analysis.

    --td <template_directory>         = /path/to/template directory. Optional. [present working directory]

    --o <output dir>                  = /path/to/output directory. Optional. [present working directory]

    --v                               = generate runtime messages. Optional

=head1  DESCRIPTION
    
    This script will combine a series of components related to the eukaryotic RNA-Seq
    analysis pipeline and create a pipeline.layout and pipeline.config file.
    The config file can then be configured with the correct options and a pipeline
    can be run. 
    
    This script will combine components from the provided templates directory and create
    a new pipeline. The following components will be looked for by this script and a config
    file for each component is expected in the templates directory:
    
    bowtie_build             : generates the bowtie index files for a given reference file.
    fastx_toolkit            : generates the quality stats and the trimmed sequence file(s) for the given sequence file(s).
    tophat                   : generates the TopHat alignment files for single-end or paired-end sequence(s).
    samtools_file_convert    : converts file formats for downstream analysis.
    samtools_alignment_stats : generates the alignment stats from the alignment BAM file(s).
    bam2bigwig               : converts the alignment BAM file(s) to BedGraph and BigWig file(s).
    rpkm_analysis            : generates rpkm coverage analysis utilizing the alignment BAM file(s).
    htseq                    : generates the count files from the alignment SAM file(s) sorted by name.
    deseq                    : generates the differential gene expression analysis results utilizing DESeq software.
    cufflinks                : generates the isoform identification analysis results utilizing the alignment BAM file(s).
    cuffcompare              : generates the isoform comparison analysis results utilizing the isoform GTF file(s).
    cuffdiff                 : generates the differential isoform expression analysis results utilizing the alignment SAM file(s).
 
    There are some restrictions about which component would precede or succeed other components.
    
=head1 AUTHOR

 Amol Carl Shetty
 Bioinformatics Software Engineer II
 Institute of Genome Sciences
 University of Maryland
 Baltimore, Maryland 21201

=cut

################################################################################

use strict;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use XML::Writer;
use Pod::Usage;
use File::Spec;
use FindBin qw($RealBin);

##############################################################################
### Constants
##############################################################################

use constant FALSE => 0;
use constant TRUE  => 1;

use constant BOWTIE_BIN_DIR => '/usr/local/bin';
use constant SAMTOOLS_BIN_DIR => '/usr/local/bin';
use constant TOPHAT_BIN_DIR => '/usr/local/bin';
use constant CUFFLINKS_BIN_DIR => '/usr/local/bin';

use constant VERSION => '2.0.0';
use constant PROGRAM => eval { ($0 =~ m/(\w+\.pl)$/) ? $1 : $0 };

##############################################################################
### Globals
##############################################################################

my @aComponents = ("bowtie_build", "quality_stats", "quality_trimming", "alignment", "split", "visualization", 
				   "rpkm_analysis", "diff_gene_expr", "isoform_analysis", "diff_isoform_analysis");
my %hCmdLineOption = ();
my $sHelpHeader = "\nThis is ".PROGRAM." version ".VERSION."\n";

GetOptions( \%hCmdLineOption,
			'sample_file|s=s', 'config_file|c=s', 'reffile|r=s', 'quality|qual=i', 'gtffile|gtf=s',
			'bowtie_build', 'quality_stats', 'quality_trimming', 'split',  
                        'alignment', 'bwtidxfile=s', 'visualization', 'rpkm_analysis', 'annotation_format=s', 
			'diff_gene_expr', 'comparison_groups=s', 'count', 'file_type=s', 'sorted=s', 
			'isoform_analysis', 'include_novel', 'diff_isoform_analysis', 'use_ref_gtf', 
			'repository_root|rr=s', 'ergatis_ini|ei=s', 
			'outdir|o=s', 'template_dir|td=s',
                        'verbose|v',
                        'debug',
                        'help',
                        'man') or pod2usage(2);

## display documentation
pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};
pod2usage( -msg => $sHelpHeader, -exitval => 1) if $hCmdLineOption{'help'};

## make sure everything passed was peachy
check_parameters(\%hCmdLineOption, \@aComponents);

my (%hConfig, %hParams, %hGroups);
my ($sOutDir, $sTemplateDir);
my ($sPLayout, $sPConfig);
my ($sBwtIndexDir, $sBwtIndexPrefix);
my ($fpPL, $fpPC, $fpLST1, $fpLST2, $fpLST, $fpSMPL, $fpGTF);
my ($sSampleName, $sGroupName, $sRead1File, $sRead2File, @aReadFiles, $sList);
my ($sSamRefFile, $sBamFileList, $sSamFileList, $sBamNameSortList, $sMapStatsList, $sCountsFileList,$Deseq_List, $Cuff_List, $sRpkmFileList);
my ($sGtfFileList, $sGtfFile, $sCuffdiff_SamFileList, $sCuffFileList);
my ($sFeature, $sAttrID);
my (@aComparisons, $sCGrp, $sGrpX, $sGrpY);
my ($sList1File, $sList2File, $sListFile, $sListBamFile, $sListFile1, $sListFile2, $sInfoListFile, $sFile, $sInFile);
my ($sTimeStamp, $sCmd, $sArgs, $nOpt, $nI, $bPE, $bGZ);
my ($nSec, $nMin, $nHour, $nMDay, $nMon, $nYear, $nWDay, $nYDay, $bDST);
my $bDebug   = (defined $hCmdLineOption{'debug'}) ? TRUE : FALSE;
my $bVerbose = (defined $hCmdLineOption{'verbose'}) ? TRUE : FALSE;

################################################################################
### Main
################################################################################

($bDebug || $bVerbose) ? print STDERR "\nInitiating Eukaryotic RNA-Seq Analysis Pipeline .....\n" : undef;

$sOutDir = File::Spec->curdir();
if (defined $hCmdLineOption{'outdir'}) {
    $sOutDir = $hCmdLineOption{'outdir'};

    if (! -e $sOutDir) {
        mkdir($sOutDir) ||
            die "ERROR! Cannot create output directory $sOutDir !!!\n";
    }
    elsif (! -d $sOutDir) {
            die "ERROR! $sOutDir is not a directory !!!\n";
    }
}
$sOutDir = File::Spec->canonpath($sOutDir);
$sOutDir = File::Spec->rel2abs($sOutDir);

$sTemplateDir = $RealBin."/../global_pipeline_templates/Eukaryotic_RNA_Seq_Analysis";
if (defined $hCmdLineOption{'template_dir'}) {
    $sTemplateDir = $hCmdLineOption{'template_dir'};
}
$sTemplateDir = File::Spec->canonpath($sTemplateDir);
$sTemplateDir = File::Spec->rel2abs($sTemplateDir);

($nSec, $nMin, $nHour, $nMDay, $nMon, $nYear, $nWDay, $nYDay, $bDST) = localtime(time);
$sTimeStamp = sprintf("%4d%02d%02d.%02d%02d%02d", ($nYear + 1900), ($nMon + 1), $nMDay, $nHour ,$nMin, $nSec);

$sOutDir .= "/Pipeline_".$sTimeStamp;
mkdir($sOutDir) || die "ERROR! Cannot create output directory $sOutDir !!!\n";

($bDebug || $bVerbose) ? print STDERR "\nProcessing $hCmdLineOption{'config_file'} .....\n" : undef;
read_config(\%hCmdLineOption, \%hConfig);
($bDebug || $bVerbose) ? print STDERR "\nProcessing $hCmdLineOption{'config_file'} ..... done\n" : undef;

($bDebug || $bVerbose) ? print STDERR "\nProcessing $hCmdLineOption{'sample_file'} .....\n" : undef;

if (defined $hCmdLineOption{'split'}) {
	$hCmdLineOption{'alignment'} = 1;
}


if ((defined $hCmdLineOption{'quality_stats'}) || 
	(defined $hCmdLineOption{'quality_trimming'}) || 
	(defined $hCmdLineOption{'alignment'})) {
	
	if (! -e "$sOutDir/raw_reads") {
		mkdir("$sOutDir/raw_reads") ||
			die "ERROR! Cannot create output directory !!!\n";
	}
	elsif (! -d "$sOutDir/raw_reads") {
			die "ERROR! $sOutDir/raw_reads is not a directory !!!\n";
	}
	
	open($fpLST1, ">$sOutDir/first_mates.list") or die "Error! Cannot open $sOutDir/first_mates.list for writing: $!";
	open($fpLST2, ">$sOutDir/second_mates.list") or die "Error! Cannot open $sOutDir/second_mates.list for writing: $!";	
	open($fpSMPL, "<$hCmdLineOption{'sample_file'}") or die "Error! Cannot open $hCmdLineOption{'sample_file'} for reading: $!";
	
	$bPE = $bGZ = FALSE;
	%hGroups = ();
	while (<$fpSMPL>) {
		$_ =~ s/\s+$//;
		next if ($_ =~ /^#/);
		next if ($_ =~ /^$/);
		
		($sSampleName, $sGroupName, $sRead1File, $sRead2File) = split(/\t/, $_);
			 if($bDebug) {
			     print STDERR "$sSampleName\t$sGroupName\t$sRead1File\t$sRead2File\n";
			 }
		die "Error! Missing Reads 1 File !!!\n" if (! (defined $sRead1File));
		
			 $sSampleName =~ s/\s+//g;
			 $sGroupName =~ s/\s+//g;

		@aReadFiles = split(/,/, $sRead1File);
		
		$sList = "";
		for ($nI = 0; $nI < @aReadFiles; $nI++) {
			$sFile = $sSampleName."_1";
			$sFile = $sOutDir."/raw_reads/".$sFile."_".($nI + 1)."_sequence.txt";
			if (is_gzipped(\%hCmdLineOption, $aReadFiles[$nI])) {
				$sFile .= ".gz";
				$bGZ = TRUE;
			}
			$sCmd = "ln -s $aReadFiles[$nI] $sFile";
			($bDebug) ? print STDERR "$sCmd\n" : undef;
			exec_command($sCmd);
			($bDebug) ? print STDERR "\n" : undef;
			
			$sList .= (($sList =~ m/^$/) ? "$sFile" : ",$sFile");
		}
		
		if ((defined $hCmdLineOption{'quality_stats'}) ||
			(defined $hCmdLineOption{'quality_trimming'})) {
			$sList = join("\n", split(/,/, $sList));
		}
		
		print $fpLST1 "$sList\n";
		
		if (defined $sRead2File) {
			@aReadFiles = split(/,/, $sRead2File);
			
			$sList = "";
			for ($nI = 0; $nI < @aReadFiles; $nI++) {
				$sFile = $sSampleName."_2";
				$sFile = $sOutDir."/raw_reads/".$sFile."_".($nI + 1)."_sequence.txt";
				if (is_gzipped(\%hCmdLineOption, $aReadFiles[$nI])) {
					$sFile .= ".gz";
					$bGZ = TRUE;
				}
				$sCmd = "ln -s $aReadFiles[$nI] $sFile";
				($bDebug) ? print STDERR "$sCmd\n" : undef;
				exec_command($sCmd);
				($bDebug) ? print STDERR "\n" : undef;
				
				$sList .= (($sList =~ m/^$/) ? "$sFile" : ",$sFile");
			}
			
			if ((defined $hCmdLineOption{'quality_stats'}) ||
				(defined $hCmdLineOption{'quality_trimming'})) {
				$sList = join("\n", split(/,/, $sList));
			}
			
			print $fpLST2 "$sList\n";
			
			$bPE = TRUE;
		}
		
		if (!(defined $hGroups{$sGroupName})) {
			my @aSampleArray = ();
			$hGroups{$sGroupName} = \@aSampleArray;
		}
		push @{$hGroups{$sGroupName}}, $sSampleName;
	}
	
	close($fpLST1);
	close($fpLST2);
	close($fpSMPL);
	
	$sList1File = "$sOutDir/first_mates.list";
	$sList2File = "$sOutDir/second_mates.list";
}
elsif (defined $hCmdLineOption{'diff_gene_expr'}) {
	die "Error! Comparison groups undefined !!!\n" if (!(defined $hCmdLineOption{'comparison_groups'}));
	
	if (defined $hCmdLineOption{'count'}) {
		if (! -e "$sOutDir/alignment") {
			mkdir("$sOutDir/alignment") ||
				die "ERROR! Cannot create output directory !!!\n";
		}
		elsif (! -d "$sOutDir/alignment") {
				die "ERROR! $sOutDir/alignment is not a directory !!!\n";
		}
		
		if ((defined $hCmdLineOption{'file_type'}) && ($hCmdLineOption{'file_type'} =~ m/BAM/i)) {
			$sBamFileList = $sListFile = "$sOutDir/alignment_bam.list";
		}
		else {
			$sSamFileList = $sListFile = "$sOutDir/alignment_sam.list";
		}
	}
	else {
		if (! -e "$sOutDir/counts") {
			mkdir("$sOutDir/counts") ||
				die "ERROR! Cannot create output directory !!!\n";
		}
		elsif (! -d "$sOutDir/counts") {
				die "ERROR! $sOutDir/counts is not a directory !!!\n";
		}
		
		$sCountsFileList = $sListFile = "$sOutDir/alignment_counts.list";
	}
	
	open($fpLST, ">$sListFile") or die "Error! Cannot open $sListFile for writing: $!";
	open($fpSMPL, "<$hCmdLineOption{'sample_file'}") or die "Error! Cannot open $hCmdLineOption{'sample_file'} for reading: $!";
	
	%hGroups = ();
	while (<$fpSMPL>) {
		$_ =~ s/\s+$//;
		next if ($_ =~ /^#/);
		next if ($_ =~ /^$/);
		
		($sSampleName, $sGroupName, $sInFile) = split(/\t/, $_);
		
		die "Error! Missing Input File !!!\n" if (! (defined $sInFile));
		
		if (defined $hCmdLineOption{'count'}) {
			if ((defined $hCmdLineOption{'file_type'}) && ($hCmdLineOption{'file_type'} =~ m/BAM/i)) {
				$sFile = $sOutDir."/alignment/".$sSampleName.".bam";
			}
			else {
				$sFile = $sOutDir."/alignment/".$sSampleName.".sam";
			}
		}
		else {
			$sFile = $sOutDir."/counts/".$sSampleName.".counts";
		}
		
		$sCmd = "ln -s $sInFile $sFile";
		($bDebug) ? print STDERR "$sCmd\n" : undef;
		exec_command($sCmd);
		($bDebug) ? print STDERR "\n" : undef;
		
		print $fpLST "$sFile\n";
		
		if (!(defined $hGroups{$sGroupName})) {
			my @aSampleArray = ();
			$hGroups{$sGroupName} = \@aSampleArray;
		}
		push @{$hGroups{$sGroupName}}, $sSampleName;
	}
	
	close($fpLST);
	close($fpSMPL);
}
elsif ((defined $hCmdLineOption{'isoform_analysis'}) || (defined $hCmdLineOption{'rpkm_analysis'}) || (defined $hCmdLineOption{'visualization'})) {
	if (! -e "$sOutDir/alignment") {
		mkdir("$sOutDir/alignment") ||
			die "ERROR! Cannot create output directory !!!\n";
	}
	elsif (! -d "$sOutDir/alignment") {
			die "ERROR! $sOutDir/alignment is not a directory !!!\n";
	}
	
	if ((defined $hCmdLineOption{'file_type'}) && ($hCmdLineOption{'file_type'} =~ m/SAM/i)) {
		$sSamFileList = $sListFile = "$sOutDir/alignment_sam.list";
	}
	else {
		$sBamFileList = $sListFile = "$sOutDir/alignment_bam.list";
	}
	
	open($fpLST, ">$sListFile") or die "Error! Cannot open $sListFile for writing: $!";
	open($fpSMPL, "<$hCmdLineOption{'sample_file'}") or die "Error! Cannot open $hCmdLineOption{'sample_file'} for reading: $!";
	
	while (<$fpSMPL>) {
		$_ =~ s/\s+$//;
		next if ($_ =~ /^#/);
		next if ($_ =~ /^$/);
		
		($sSampleName, $sGroupName, $sInFile) = split(/\t/, $_);
		
		die "Error! Missing Input File !!!\n" if (! (defined $sInFile));
		
		if ((defined $hCmdLineOption{'file_type'}) && ($hCmdLineOption{'file_type'} =~ m/SAM/i)) {
			$sFile = $sOutDir."/alignment/".$sSampleName.".sam";
		}
		else {
			$sFile = $sOutDir."/alignment/".$sSampleName.".bam";
		}
	
		$sCmd = "ln -s $sInFile $sFile";
		($bDebug) ? print STDERR "$sCmd\n" : undef;
		exec_command($sCmd);
		($bDebug) ? print STDERR "\n" : undef;
		
		print $fpLST "$sFile\n";
	}
	
	close($fpLST);
	close($fpSMPL);
}
elsif (defined $hCmdLineOption{'diff_isoform_analysis'}) {
	die "Error! Comparison groups undefined !!!\n" if (!(defined $hCmdLineOption{'comparison_groups'}));
	
	if (! -e "$sOutDir/alignment") {
		mkdir("$sOutDir/alignment") ||
			die "ERROR! Cannot create output directory !!!\n";
	}
	elsif (! -d "$sOutDir/alignment") {
			die "ERROR! $sOutDir/alignment is not a directory !!!\n";
	}
	
	if (! -e "$sOutDir/isoform") {
		mkdir("$sOutDir/isoform") ||
			die "ERROR! Cannot create output directory !!!\n";
	}
	elsif (! -d "$sOutDir/isoform") {
			die "ERROR! $sOutDir/isoform is not a directory !!!\n";
	}
	
	if ((defined $hCmdLineOption{'file_type'}) && ($hCmdLineOption{'file_type'} =~ m/SAM/i)) {
		$sSamFileList = $sListFile = "$sOutDir/alignment_sam.list";
	}
	else {
		$sBamFileList = $sListFile = "$sOutDir/alignment_bam.list";
	}
	
	$sGtfFileList = "$sOutDir/cufflinks_gtf.list";
	
	open($fpLST, ">$sListFile") or die "Error! Cannot open $sListFile for writing: $!";
	open($fpGTF, ">$sGtfFileList") or die "Error! Cannot open $sGtfFileList for writing: $!";
	open($fpSMPL, "<$hCmdLineOption{'sample_file'}") or die "Error! Cannot open $hCmdLineOption{'sample_file'} for reading: $!";
	
	while (<$fpSMPL>) {
		$_ =~ s/\s+$//;
		next if ($_ =~ /^#/);
		next if ($_ =~ /^$/);
		
		($sSampleName, $sGroupName, $sInFile, $sGtfFile) = split(/\t/, $_);
		
		die "Error! Missing Input File !!!\n" if ((! (defined $sInFile)) || (! (defined $sGtfFile)));
		
		if ((defined $hCmdLineOption{'file_type'}) && ($hCmdLineOption{'file_type'} =~ m/SAM/i)) {
			$sFile = $sOutDir."/alignment/".$sSampleName.".sam";
		}
		else {
			$sFile = $sOutDir."/alignment/".$sSampleName.".bam";
		}
		
		$sCmd = "ln -s $sInFile $sFile";
		($bDebug) ? print STDERR "$sCmd\n" : undef;
		exec_command($sCmd);
		($bDebug) ? print STDERR "\n" : undef;
		
		print $fpLST "$sFile\n";
		
		$sFile = $sOutDir."/isoform/".$sSampleName.".trancripts.gtf";
		
		$sCmd = "ln -s $sGtfFile $sFile";
		($bDebug) ? print STDERR "$sCmd\n" : undef;
		exec_command($sCmd);
		($bDebug) ? print STDERR "\n" : undef;
		
		print $fpGTF "$sFile\n";
	}
	
	close($fpLST);
	close($fpGTF);
	close($fpSMPL);
}

if (defined $hCmdLineOption{'diff_gene_expr'}) {
	die "Error! Comparison groups undefined !!!\n" if (!(defined $hCmdLineOption{'comparison_groups'}));
	@aComparisons = split(/,/, $hCmdLineOption{'comparison_groups'});
	
	open($fpLST, ">$sOutDir/deseq_sample_info.list") or die "Error! Cannot open $sOutDir/deseq_sample_info.list for writing: $!";
	
	foreach $sCGrp (@aComparisons) {
		($sGrpX, $sGrpY) = split(/vs/, $sCGrp);
		
		open($fpSMPL, ">$sOutDir/$sCGrp\_sample.list") or die "Error! Cannot open $sOutDir/$sCGrp\_sample.list for writing: $!";
		print $fpSMPL "#Sample_Name\tGroup_Name\n";
		
		for ($nI = 0; $nI < @{$hGroups{$sGrpX}}; $nI++) {
			print $fpSMPL $hGroups{$sGrpX}->[$nI]."\t".$sGrpX."\n";
		}
		
		for ($nI = 0; $nI < @{$hGroups{$sGrpY}}; $nI++) {
			print $fpSMPL $hGroups{$sGrpY}->[$nI]."\t".$sGrpY."\n";
		}
		
		close($fpSMPL);
		
		print $fpLST "$sOutDir/$sCGrp\_sample.list\n";
	}
	
	close($fpLST);
}

# Path to pipeline.layout and pipeline.config files
$sPLayout = $sOutDir."/pipeline.".$sTimeStamp.".layout";
$sPConfig = $sOutDir."/pipeline.".$sTimeStamp.".config";

($bDebug || $bVerbose) ? print STDERR "\nGenerating $sPLayout & $sPConfig .....\n" : undef;

# File handle for pipeline.layout and pipeline.config files
open($fpPL, ">$sPLayout") or die "Error! Cannot open $sPLayout for writing: $!";
open($fpPC, ">$sPConfig") or die "Error! Cannot open $sPConfig for writing: $!";

# Since the pipeline.layout is XML, create an XML::Writer
my $oPL = new XML::Writer( 'OUTPUT' => $fpPL, 'DATA_MODE' => 1, 'DATA_INDENT' => 3 );

# Initiate pipeline.layout
init_pipeline_layout($oPL);

if (defined $hCmdLineOption{'bowtie_build'}) {
	die "Error! Reference FastA file undefined !!!\n" if (!(defined $hCmdLineOption{'reffile'}));
	
	if (defined ($hConfig{'bowtie_build'}{'BOWTIE_INDEX_PREFIX'})) {
		$sBwtIndexPrefix = $hConfig{'bowtie_build'}{'BOWTIE_INDEX_PREFIX'}[0];
	}
	else {
		($_, $_, $sBwtIndexPrefix) = File::Spec->splitpath($hCmdLineOption{'reffile'});
		$sBwtIndexPrefix =~ s/.(\w+)$//;
		$hConfig{'bowtie_build'}{'BOWTIE_INDEX_PREFIX'}[0] = $sBwtIndexPrefix;
		$hConfig{'bowtie_build'}{'BOWTIE_INDEX_PREFIX'}[1] = "bowtie index prefix";
	}
	
	###	Add Bowtie Build Component & Parameters ###
	init_component($oPL, "serial");
		include_component_layout($oPL, $sTemplateDir, "bowtie_build", "reference");
	complete_component($oPL);
	
	%hParams = ();
	$hParams{'INPUT_FILE'} = ["$hCmdLineOption{'reffile'}", "path to reference FastA file"];
	config2params(\%hParams, \%hConfig, 'bowtie_build');
	add_config_section($fpPC, "bowtie_build", "reference");
	add_config_parameters($fpPC, \%hParams);
	
	$sBwtIndexDir = '$;REPOSITORY_ROOT$;/output_repository/bowtie_build/$;PIPELINEID$;_reference/i1/g1';
}

if (defined $hCmdLineOption{'quality_stats'}) {
	###	Add FastX Quality Stats Component ###
	init_component($oPL, "serial");
		init_component($oPL, "parallel");
			include_component_layout($oPL, $sTemplateDir, "fastx_quality_stats", "read1");
			include_component_layout($oPL, $sTemplateDir, "fastx_quality_stats", "read2") if ($bPE);
		complete_component($oPL);
	complete_component($oPL);
	
	###	Add FastX Quality Stats Parameters for First Mates ###
	%hParams = ();
	$hParams{'INPUT_FILE_LIST'} = ["$sList1File", "path to list of input FastQ/A sequence files"];
	$hParams{'QUALITY_STRING'} = ["$hCmdLineOption{'quality'}", "quality string type for FastQ files (33 or 64)"];
	
	$sArgs = "";
	$sArgs = $hConfig{'fastx_quality_stats'}{'OTHER_ARGS'}[0] if (defined $hConfig{'fastx_quality_stats'}{'OTHER_ARGS'});
	$sArgs .= " --z" if (($bGZ) && ($sArgs !~ m/--z/));
	$hConfig{'fastx_quality_stats'}{'OTHER_ARGS'}[0] = $sArgs;
	$hConfig{'fastx_quality_stats'}{'OTHER_ARGS'}[1] = "Additional Arguments" if (! defined $hConfig{'fastx_quality_stats'}{'OTHER_ARGS'});
	config2params(\%hParams, \%hConfig, 'fastx_quality_stats');
	add_config_section($fpPC, "fastx_quality_stats", "read1");
	add_config_parameters($fpPC, \%hParams);
	
	if ($bPE) {
		###	Add FastX Quality Stats Parameters for Second Mates ###
		%hParams = ();
		$hParams{'INPUT_FILE_LIST'} = ["$sList2File", "path to list of input FastQ/A sequence files"];
		$hParams{'QUALITY_STRING'} = ["$hCmdLineOption{'quality'}", "quality string type for FastQ files (33 or 64)"];
		
		$sArgs = "";
		$sArgs = $hConfig{'fastx_quality_stats'}{'OTHER_ARGS'}[0] if (defined $hConfig{'fastx_quality_stats'}{'OTHER_ARGS'});
		$sArgs .= " --z" if (($bGZ) && ($sArgs !~ m/--z/));
		$hConfig{'fastx_quality_stats'}{'OTHER_ARGS'}[0] = $sArgs;
		$hConfig{'fastx_quality_stats'}{'OTHER_ARGS'}[1] = "Additional Arguments" if (! defined $hConfig{'fastx_quality_stats'}{'OTHER_ARGS'});
		config2params(\%hParams, \%hConfig, 'fastx_quality_stats');
		add_config_section($fpPC, "fastx_quality_stats", "read2");
		add_config_parameters($fpPC, \%hParams);
	}
}

if (defined $hCmdLineOption{'quality_trimming'}) {
	###	Add FastX Quality Trimming Component ###
	init_component($oPL, "serial");
		init_component($oPL, "parallel");
			include_component_layout($oPL, $sTemplateDir, "fastx_trimming", "read1");
			include_component_layout($oPL, $sTemplateDir, "fastx_trimming", "read2") if ($bPE);;
		complete_component($oPL);
	complete_component($oPL);
	
	###	Add FastX Quality Stats Parameters for First Mates ###
	%hParams = ();
	$hParams{'INPUT_FILE_LIST'} = ["$sList1File", "path to list of input FastQ/A sequence files"];
	$hParams{'QUALITY_STRING'} = ["$hCmdLineOption{'quality'}", "quality string type for FastQ files (33 or 64)"];
	
	$sArgs = "";
	$sArgs = $hConfig{'fastx_trimming'}{'OTHER_ARGS'}[0] if (defined $hConfig{'fastx_trimming'}{'OTHER_ARGS'});
	$sArgs .= " --qs" if ((! defined $hConfig{'fastx_trimming'}{'LAST_BASE'}) && ($sArgs !~ m/--qs/));
	$sArgs .= " --z" if (($bGZ) && ($sArgs !~ m/--z/));
	$hConfig{'fastx_trimming'}{'OTHER_ARGS'}[0] = $sArgs;
	$hConfig{'fastx_trimming'}{'OTHER_ARGS'}[1] = "Additional Arguments" if (! defined $hConfig{'fastx_trimming'}{'OTHER_ARGS'});
	config2params(\%hParams, \%hConfig, 'fastx_trimming');
	add_config_section($fpPC, "fastx_trimming", "read1");
	add_config_parameters($fpPC, \%hParams);
	
	$sList1File = '$;REPOSITORY_ROOT$;/output_repository/fastx_trimming/$;PIPELINEID$;_read1/fastx_trimming.trimmed.sequence.list';
	
	if ($bPE) {
		###	Add FastX Quality Stats Parameters for Second Mates ###
		%hParams = ();
		$hParams{'INPUT_FILE_LIST'} = ["$sList2File", "path to list of input FastQ/A sequence files"];
		$hParams{'QUALITY_STRING'} = ["$hCmdLineOption{'quality'}", "quality string type for FastQ files (33 or 64)"];
		
		$sArgs = "";
		$sArgs = $hConfig{'fastx_trimming'}{'OTHER_ARGS'}[0] if (defined $hConfig{'fastx_trimming'}{'OTHER_ARGS'});
		$sArgs .= " --qs" if ((! defined $hConfig{'fastx_trimming'}{'LAST_BASE'}) && ($sArgs !~ m/--qs/));
		$sArgs .= " --z" if (($bGZ) && ($sArgs !~ m/--z/));
		$hConfig{'fastx_trimming'}{'OTHER_ARGS'}[0] = $sArgs;
		$hConfig{'fastx_trimming'}{'OTHER_ARGS'}[1] = "Additional Arguments" if (! defined $hConfig{'fastx_trimming'}{'OTHER_ARGS'});
		config2params(\%hParams, \%hConfig, 'fastx_trimming');
		add_config_section($fpPC, "fastx_trimming", "read2");
		add_config_parameters($fpPC, \%hParams);
		
		$sList2File = '$;REPOSITORY_ROOT$;/output_repository/fastx_trimming/$;PIPELINEID$;_read2/fastx_trimming.trimmed.sequence.list';
	}
}

if (defined $hCmdLineOption{'reffile'}) {
	###	Add Samtools Reference Index Component & Parameters ###
	init_component($oPL, "serial");
		include_component_layout($oPL, $sTemplateDir, "samtools_reference_index", "reference");
	complete_component($oPL);
	
	%hParams = ();
	$hParams{'INPUT_FILE'} = ["$hCmdLineOption{'reffile'}", "path to reference FastA file"];
	add_config_section($fpPC, "samtools_reference_index", "reference");
	add_config_parameters($fpPC, \%hParams);
	
	($_, $_, $sSamRefFile) = File::Spec->splitpath($hCmdLineOption{'reffile'});
	$sSamRefFile = '$;REPOSITORY_ROOT$;/output_repository/samtools_reference_index/$;PIPELINEID$;_reference/i1/g1/'.$sSamRefFile;
}

if (defined $hCmdLineOption{'alignment'} || defined $hCmdLineOption{'split'}) {


	if (defined $hCmdLineOption{'split'}) {
		init_component($oPL, "parallel");

		init_component($oPL, "serial");
		if ($bPE) {
			###	Add Create Pairwise List Component & Parameters ###
			include_component_layout($oPL, $sTemplateDir, "create_paired_list_file", "list");
			
			%hParams = ();
			$hParams{'LIST_FILE_1'} = ["$sList1File", "path to list file of input file 1s"];
			$hParams{'LIST_FILE_2'} = ["$sList2File", "path to list file of input file 2s"];
			$hParams{'SAMPLE_INFO'} = ["$hCmdLineOption{'sample_file'}", "path to sample info file with information on all samples to be analyzed"] if (defined $hCmdLineOption{'quality_trimming'});
			add_config_section($fpPC, "create_paired_list_file", "list");
			add_config_parameters($fpPC, \%hParams);
			
			$sListFile = '$;REPOSITORY_ROOT$;/output_repository/create_paired_list_file/$;PIPELINEID$;_list/paired_input_file.list';
		}
		else {
			$sListFile = $sList1File;
		}
		
		###	Add FastQC Stats  Components ###
		include_component_layout($oPL, $sTemplateDir, "fastqc_stats", "fastqc");
		complete_component($oPL);
		
		###	Add FastQC Stats Parameters ###
		%hParams = ();
		$hParams{'INPUT_FILE_LIST'} = ["$sListFile", "path to list file consisting of tab separated first mate and second mate sequence files"];
		add_config_section($fpPC, "fastqc_stats", "fastqc");
		add_config_parameters($fpPC, \%hParams);
		
  
		init_component($oPL, "serial");
		if ($bPE) {

			init_component($oPL, "parallel");
			include_component_layout($oPL, $sTemplateDir, "split_fastq", "first_mate");
			include_component_layout($oPL, $sTemplateDir, "split_fastq", "second_mate");
			complete_component($oPL);

			#Add parameters to split_fastq_first_mate
			%hParams = ();
			$hParams{'INPUT_FILE_LIST'} = ["$sList1File", "path to list file of input file 1s"];
                        config2params(\%hParams, \%hConfig, 'split_fastq');
			add_config_section($fpPC, "split_fastq", "first_mate");
			add_config_parameters($fpPC, \%hParams);
			$sListFile1 = '$;REPOSITORY_ROOT$;/output_repository/split_fastq/$;PIPELINEID$;_first_mate/split_fastq.fastq.list';
			$sInfoListFile = '$;REPOSITORY_ROOT$;/output_repository/split_fastq/$;PIPELINEID$;_first_mate/split_fastq.fastq.info.list';			

                       #Add parameters to split_fastq_second_mate
			%hParams = ();
			$hParams{'INPUT_FILE_LIST'} = ["$sList2File", "path to list file of input file 2s"];
                        config2params(\%hParams, \%hConfig, 'split_fastq');
			add_config_section($fpPC, "split_fastq", "second_mate");
			add_config_parameters($fpPC, \%hParams);
			$sListFile2 = '$;REPOSITORY_ROOT$;/output_repository/split_fastq/$;PIPELINEID$;_second_mate/split_fastq.fastq.list';

			###	Add Create Pairwise List Component & Parameters ###
			include_component_layout($oPL, $sTemplateDir, "create_paired_list_file", "split");
			
			%hParams = ();
			$hParams{'LIST_FILE_1'} = ["$sListFile1", "path to split_list file of 1s"];
			$hParams{'LIST_FILE_2'} = ["$sListFile2", "path to split_list file of 2s"];
			add_config_section($fpPC, "create_paired_list_file", "split");
			add_config_parameters($fpPC, \%hParams);
			$sListFile = '$;REPOSITORY_ROOT$;/output_repository/create_paired_list_file/$;PIPELINEID$;_split/paired_input_file.list';
		}
		else {

			include_component_layout($oPL, $sTemplateDir, "split_fastq", "first_mate");

			#Add parameters to split_fastq_first_mate
			%hParams = ();
			$hParams{'INPUT_FILE_LIST'} = ["$sList1File", "path to list file of input file 1s"];
			$hParams{'SEQ_NUMBER'} = ["20000000", "path to list file of input file 1s"];
			add_config_section($fpPC, "split_fastq", "first_mate");
			add_config_parameters($fpPC, \%hParams);
			$sListFile1 = '$;REPOSITORY_ROOT$;/output_repository/split_fastq/$;PIPELINEID$;_first_mate/split_fastq.fastq.list';
			$sInfoListFile = '$;REPOSITORY_ROOT$;/output_repository/split_fastq/$;PIPELINEID$;_first_mate/split_fastq.fastq.info.list';			

			$sListFile = $sListFile1;
		}
		

		if ((!(defined $sBwtIndexDir)) || (!(defined $sBwtIndexPrefix))) {
			die "Error! Bowtie Index undefined !!!\n" if (!(defined $hCmdLineOption{'bwtidxfile'}));
			
			($_, $sBwtIndexDir, $sBwtIndexPrefix) = File::Spec->splitpath($hCmdLineOption{'bwtidxfile'});
		}
	
	
		###	Tophat Components ###
		include_component_layout($oPL, $sTemplateDir, "tophat", "alignment");
		
		###	Add TopHat Parameters ###
		%hParams = ();
		$hParams{'INPUT_FILE_LIST'} = ["$sListFile", "path to list file consisting of tab separated first mate and second mate sequence files"];
		$hParams{'BOWTIE_INDEX_DIR'} = ["$sBwtIndexDir", "path to bowtie package binary directory"];
		$hParams{'BOWTIE_INDEX_PREFIX'} = ["$sBwtIndexPrefix", "bowtie index prefix"];
	
		$sArgs = "";
		$sArgs = $hConfig{'tophat'}{'OTHER_ARGS'}[0] if (defined $hConfig{'tophat'}{'OTHER_ARGS'});
		$sArgs .= " --solexa1.3-quals" if (($hCmdLineOption{'quality'} == 64) && ($sArgs !~ m/--solexa1.3-quals/));
		$hConfig{'tophat'}{'OTHER_ARGS'}[0] = $sArgs;
		$hConfig{'tophat'}{'OTHER_ARGS'}[1] = "Additional Arguments" if (! defined $hConfig{'tophat'}{'OTHER_ARGS'});
		config2params(\%hParams, \%hConfig, 'tophat');
		add_config_section($fpPC, "tophat", "alignment");
		add_config_parameters($fpPC, \%hParams);
		
		$sListBamFile = '$;REPOSITORY_ROOT$;/output_repository/tophat/$;PIPELINEID$;_alignment/tophat.bam.list';


		###	BAM merge Components ###
		include_component_layout($oPL, $sTemplateDir, "bam_merge", "list");
		
		###	Add BAM Merge Parameters ###
		%hParams = ();
		$hParams{'INPUT_FILE_LIST'} = ["$sInfoListFile", "path to list of filea of fastq split information"];
		$hParams{'INPUT_BAM_LIST'} = ["$sListBamFile", "path to list of bam file"];

		add_config_section($fpPC, "bam_merge", "list");
		add_config_parameters($fpPC, \%hParams);
		
		$sListFile = '$;REPOSITORY_ROOT$;/output_repository/bam_merge/$;PIPELINEID$;_list/bam_merge.file.list';

		complete_component($oPL);
		
		complete_component($oPL);

	}
	else {
		init_component($oPL, "serial");
		if ($bPE) {
			###	Add Create Pairwise List Component & Parameters ###
			include_component_layout($oPL, $sTemplateDir, "create_paired_list_file", "list");
			
			%hParams = ();
			$hParams{'LIST_FILE_1'} = ["$sList1File", "path to list file of input file 1s"];
			$hParams{'LIST_FILE_2'} = ["$sList2File", "path to list file of input file 2s"];
			$hParams{'SAMPLE_INFO'} = ["$hCmdLineOption{'sample_file'}", "path to sample info file with information on all samples to be analyzed"] if (defined $hCmdLineOption{'quality_trimming'});
			add_config_section($fpPC, "create_paired_list_file", "list");
			add_config_parameters($fpPC, \%hParams);
			
			$sListFile = '$;REPOSITORY_ROOT$;/output_repository/create_paired_list_file/$;PIPELINEID$;_list/paired_input_file.list';
		}
		else {
			$sListFile = $sList1File;
		}
		
		if ((!(defined $sBwtIndexDir)) || (!(defined $sBwtIndexPrefix))) {
			die "Error! Bowtie Index undefined !!!\n" if (!(defined $hCmdLineOption{'bwtidxfile'}));
			
			($_, $sBwtIndexDir, $sBwtIndexPrefix) = File::Spec->splitpath($hCmdLineOption{'bwtidxfile'});
		}
		
		###	Add FastQC Stats and Tophat Components ###
		init_component($oPL, "parallel");
		include_component_layout($oPL, $sTemplateDir, "fastqc_stats", "fastqc");
		include_component_layout($oPL, $sTemplateDir, "tophat", "alignment");
		complete_component($oPL);
		
		###	Add FastQC Stats Parameters ###
		%hParams = ();
		$hParams{'INPUT_FILE_LIST'} = ["$sListFile", "path to list file consisting of tab separated first mate and second mate sequence files"];
		add_config_section($fpPC, "fastqc_stats", "fastqc");
		add_config_parameters($fpPC, \%hParams);
		
		###	Add TopHat Parameters ###
		%hParams = ();
		$hParams{'INPUT_FILE_LIST'} = ["$sListFile", "path to list file consisting of tab separated first mate and second mate sequence files"];
		$hParams{'BOWTIE_INDEX_DIR'} = ["$sBwtIndexDir", "path to bowtie package binary directory"];
		$hParams{'BOWTIE_INDEX_PREFIX'} = ["$sBwtIndexPrefix", "bowtie index prefix"];
	
		$sArgs = "";
		$sArgs = $hConfig{'tophat'}{'OTHER_ARGS'}[0] if (defined $hConfig{'tophat'}{'OTHER_ARGS'});
		$sArgs .= " --solexa1.3-quals" if (($hCmdLineOption{'quality'} == 64) && ($sArgs !~ m/--solexa1.3-quals/));
		$hConfig{'tophat'}{'OTHER_ARGS'}[0] = $sArgs;
		$hConfig{'tophat'}{'OTHER_ARGS'}[1] = "Additional Arguments" if (! defined $hConfig{'tophat'}{'OTHER_ARGS'});
		config2params(\%hParams, \%hConfig, 'tophat');
		add_config_section($fpPC, "tophat", "alignment");
		add_config_parameters($fpPC, \%hParams);
		
		$sListFile = '$;REPOSITORY_ROOT$;/output_repository/tophat/$;PIPELINEID$;_alignment/tophat.bam.list';
		
		complete_component($oPL);
		
	}
	
	init_component($oPL, "serial");
	
	###	Add Samtools File Conversion Component & Parameters ###
	init_component($oPL, "parallel");
	include_component_layout($oPL, $sTemplateDir, "samtools_file_convert", "sorted_position");
	include_component_layout($oPL, $sTemplateDir, "samtools_file_convert", "sorted_name");
	complete_component($oPL);
	
	%hParams = ();
	$hParams{'INPUT_FILE_LIST'} = ["$sListFile", "path to list of alignment files"];
	$hParams{'INPUT_FILE_FORMAT'} = ["BAM", "input alignment file format (BAM or SAM)"];
	$nOpt = "12";
	$hParams{'OPTIONS'} = ["$nOpt", "string of options for file conversion (eg : 123). 1 - BAM to sorted BAM, 2 - sorted BAM to indexed BAM, 3 - BAM to SAM, and 4 - SAM to BAM"];
	add_config_section($fpPC, "samtools_file_convert", "sorted_position");
	add_config_parameters($fpPC, \%hParams);
	
	$sBamFileList = '$;REPOSITORY_ROOT$;/output_repository/samtools_file_convert/$;PIPELINEID$;_sorted_position/samtools_file_convert.sorted_by_position_bam.list';
	
	%hParams = ();
	$hParams{'INPUT_FILE_LIST'} = ["$sListFile", "path to list of alignment files"];
	$hParams{'INPUT_FILE_FORMAT'} = ["BAM", "input alignment file format (BAM or SAM)"];
	$hParams{'SAMTOOLS_SORT_PARAMETERS'} = ["-n", "samtools sort parameters"];
	$nOpt = "13";
	$hParams{'OPTIONS'} = ["$nOpt", "string of options for file conversion (eg : 123). 1 - BAM to sorted BAM, 2 - sorted BAM to indexed BAM, 3 - BAM to SAM, and 4 - SAM to BAM"];
	add_config_section($fpPC, "samtools_file_convert", "sorted_name");
	add_config_parameters($fpPC, \%hParams);
	
	$sSamFileList = '$;REPOSITORY_ROOT$;/output_repository/samtools_file_convert/$;PIPELINEID$;_sorted_name/samtools_file_convert.sorted_by_name_sam.list';
	$sBamNameSortList = '$;REPOSITORY_ROOT$;/output_repository/samtools_file_convert/$;PIPELINEID$;_sorted_name/samtools_file_convert.sorted_by_name_bam.list';
	
	include_component_layout($oPL, $sTemplateDir, "samtools_alignment_stats", "alignment_stats");
	
	%hParams = ();
	$hParams{'INPUT_FILE_LIST'} = ["$sBamFileList", "path to list of alignment BAM files"];
	add_config_section($fpPC, "samtools_alignment_stats", "alignment_stats");
	add_config_parameters($fpPC, \%hParams);
	
	$sMapStatsList = '$;REPOSITORY_ROOT$;/output_repository/samtools_alignment_stats/$;PIPELINEID$;_alignment_stats/samtools_alignment_stats.mapstats.list';
	
	###	Add TopHat Stats Component & Parameters below ###

	if  (!(defined $hCmdLineOption{'split'})){
		include_component_layout($oPL, $sTemplateDir, "align_tophat_stats", "tophat_stats");
	
		%hParams = ();
		$hParams{'INPUT_FILE'} = ["$sListFile", "path to list of alignment BAM files"];
		add_config_section($fpPC, "align_tophat_stats", "tophat_stats");
		add_config_parameters($fpPC, \%hParams);
	}
	else{
		include_component_layout($oPL, $sTemplateDir, "align_tophat_split_stats", "tophat_stats");
	
		%hParams = ();
		$hParams{'INPUT_FILE'} = ["$sListFile", "path to list of merged BAM files"];
		$hParams{'SPLITBAM_LIST'} = ["$sListBamFile", "path to list of BAM files from alignment of split fastq file"];
		add_config_section($fpPC, "align_tophat_split_stats", "tophat_stats");
		add_config_parameters($fpPC, \%hParams);

	}
	
	if ((defined $hCmdLineOption{'gtffile'}) && (defined $hCmdLineOption{'annotation_format'}) ) {
		###	Add Percent Mapped Component & Parameters below ###
		include_component_layout($oPL, $sTemplateDir, "percent_mapped_stats", "percent_mapped");
		
		%hParams = ();
		$hParams{'INPUT_FILE_LIST'} = ["$sBamNameSortList", "path to list of alignment BAM files"];
		$hParams{'REFERENCE_FASTA'} = ["$sSamRefFile", "path to reference FastA file"];
		$hParams{'ANNOTATION_FILE'} = ["$hCmdLineOption{'gtffile'}", "path to annotation file (BED or GTF or GFF3 format file)"];
		$hParams{'ANNO_FORMAT'} = ["$hCmdLineOption{'annotation_format'}", "annotation file format (bed/gtf/gff3)"];
		$hParams{'ORG_TYPE'} = ["euk", "Organism type (prok/euk)"];
		
		config2params(\%hParams, \%hConfig, 'percent_mapped_stats');
		add_config_section($fpPC, "percent_mapped_stats", "percent_mapped");
		add_config_parameters($fpPC, \%hParams);
	}
	
	complete_component($oPL);
}

if ( (defined $hCmdLineOption{'diff_gene_expr'}) || (defined $hCmdLineOption{'visualization'}) || (defined $hCmdLineOption{'rpkm_analysis'}) ||
	 ((defined $hCmdLineOption{'isoform_analysis'}) || (defined $hCmdLineOption{'diff_isoform_analysis'})) ) {
	init_component($oPL, "parallel");
}

if (defined $hCmdLineOption{'visualization'}) {
	init_component($oPL, "serial");
	
	die "Error! Reference file undefined !!!\n" if (!(defined $hCmdLineOption{'reffile'}));
	
	###	Add Visualization Component & Parameters below ###
	init_component($oPL, "serial");
		if (! defined $hCmdLineOption{'alignment'}) {
			if (!((defined $hCmdLineOption{'sorted'}) && ($hCmdLineOption{'sorted'} =~ m/position/))) {
				include_component_layout($oPL, $sTemplateDir, "samtools_file_convert", "sorted_position");
			}
		}
		include_component_layout($oPL, $sTemplateDir, "bam2bigwig", "visualization");
	complete_component($oPL);
	
	if (! defined $hCmdLineOption{'alignment'}) {
		if (!((defined $hCmdLineOption{'sorted'}) && ($hCmdLineOption{'sorted'} =~ m/position/))) {
			%hParams = ();
			if ((defined $hCmdLineOption{'file_type'}) && ($hCmdLineOption{'file_type'} =~ m/SAM/i)) {
				$hParams{'INPUT_FILE_LIST'} = ["$sSamFileList", "path to list of alignment files"];
				$hParams{'INPUT_FILE_FORMAT'} = ["SAM", "input alignment file format (BAM or SAM)"];
				$nOpt = "412";
				$hParams{'OPTIONS'} = ["$nOpt", "string of options for file conversion (eg : 123). 1 - BAM to sorted BAM, 2 - sorted BAM to indexed BAM, 3 - BAM to SAM, and 4 - SAM to BAM"];
			}
			else {
				$hParams{'INPUT_FILE_LIST'} = ["$sBamFileList", "path to list of alignment files"];
				$hParams{'INPUT_FILE_FORMAT'} = ["BAM", "input alignment file format (BAM or SAM)"];
				$nOpt = "12";
				$hParams{'OPTIONS'} = ["$nOpt", "string of options for file conversion (eg : 123). 1 - BAM to sorted BAM, 2 - sorted BAM to indexed BAM, 3 - BAM to SAM, and 4 - SAM to BAM"];
			}
			add_config_section($fpPC, "samtools_file_convert", "sorted_position");
			add_config_parameters($fpPC, \%hParams);
			
			$sBamFileList = '$;REPOSITORY_ROOT$;/output_repository/samtools_file_convert/$;PIPELINEID$;_sorted_position/samtools_file_convert.sorted_by_position_bam.list';
		}
	}
	
	%hParams = ();
	$hParams{'INPUT_FILE_LIST'} = ["$sBamFileList", "path to list of alignment BAM files"];
	$hParams{'REFERENCE_FASTA'} = ["$sSamRefFile", "path to reference FastA file"];
	
	config2params(\%hParams, \%hConfig, 'bam2bigwig');
	add_config_section($fpPC, "bam2bigwig", "visualization");
	add_config_parameters($fpPC, \%hParams);
	
	complete_component($oPL);
}

if (defined $hCmdLineOption{'rpkm_analysis'}) {
	init_component($oPL, "serial");
	
	die "Error! Reference file undefined !!!\n" if (!(defined $hCmdLineOption{'reffile'}));
	die "Error! Annotation file undefined !!!\n" if (!(defined $hCmdLineOption{'gtffile'}));
	
	###	Add RPKM Analysis Component & Parameters below ###
	init_component($oPL, "serial");
		if (! defined $hCmdLineOption{'alignment'}) {
			if (!((defined $hCmdLineOption{'sorted'}) && ($hCmdLineOption{'sorted'} =~ m/position/))) {
				include_component_layout($oPL, $sTemplateDir, "samtools_file_convert", "sorted_position");
			}
		}
		include_component_layout($oPL, $sTemplateDir, "rpkm_coverage_stats", "rpkm_cvg");
	        include_component_layout($oPL, $sTemplateDir, "wrapper_align", "wrap");
  	        include_component_layout($oPL, $sTemplateDir, "expression_plots", "rpkm"); 
	complete_component($oPL);
	
	if (! defined $hCmdLineOption{'alignment'}) {
		if (!((defined $hCmdLineOption{'sorted'}) && ($hCmdLineOption{'sorted'} =~ m/position/))) {
			%hParams = ();
			if ((defined $hCmdLineOption{'file_type'}) && ($hCmdLineOption{'file_type'} =~ m/SAM/i)) {
				$hParams{'INPUT_FILE_LIST'} = ["$sSamFileList", "path to list of alignment files"];
				$hParams{'INPUT_FILE_FORMAT'} = ["SAM", "input alignment file format (BAM or SAM)"];
				$nOpt = "412";
				$hParams{'OPTIONS'} = ["$nOpt", "string of options for file conversion (eg : 123). 1 - BAM to sorted BAM, 2 - sorted BAM to indexed BAM, 3 - BAM to SAM, and 4 - SAM to BAM"];
			}
			else {
				$hParams{'INPUT_FILE_LIST'} = ["$sBamFileList", "path to list of alignment files"];
				$hParams{'INPUT_FILE_FORMAT'} = ["BAM", "input alignment file format (BAM or SAM)"];
				$nOpt = "12";
				$hParams{'OPTIONS'} = ["$nOpt", "string of options for file conversion (eg : 123). 1 - BAM to sorted BAM, 2 - sorted BAM to indexed BAM, 3 - BAM to SAM, and 4 - SAM to BAM"];
			}
			add_config_section($fpPC, "samtools_file_convert", "sorted_position");
			add_config_parameters($fpPC, \%hParams);
			
			$sBamFileList = '$;REPOSITORY_ROOT$;/output_repository/samtools_file_convert/$;PIPELINEID$;_sorted_position/samtools_file_convert.sorted_by_position_bam.list';
		}
	}
	
	%hParams = ();
	$hParams{'INPUT_FILE_LIST'} = ["$sBamFileList", "path to list of alignment BAM files"];
	$hParams{'REFERENCE_FASTA'} = ["$sSamRefFile", "path to reference FastA file"];
	$hParams{'ANNOTATION_FILE'} = ["$hCmdLineOption{'gtffile'}", "path to annotation file (BED or GTF or GFF3 format file)"];
	$hParams{'ANNOTATION_FILE_TYPE'} = ["$hCmdLineOption{'annotation_format'}", "annotation file format (bed/gtf/gff3)"];
	$hParams{'REGION_TYPE'} = ["genic", "region to determine coverage for (genomic:genic:exonic) separated by ':' or ',' or ';'"];
	
	config2params(\%hParams, \%hConfig, 'rpkm_coverage_stats');
	add_config_section($fpPC, "rpkm_coverage_stats", "rpkm_cvg");
	add_config_parameters($fpPC, \%hParams);
	
	$sRpkmFileList = '$;REPOSITORY_ROOT$;/output_repository/rpkm_coverage_stats/$;PIPELINEID$;_rpkm_cvg/rpkm_coverage_stats.rpkm.stats.list';
	%hParams = ();
	$hParams{'PIPELINE_ID'} = ["\$;PIPELINEID\$;", "ergatis pipeline id"];
	$hParams{'OUTPUT_REPOSITORY'} = ["\$;REPOSITORY_ROOT\$;/output_repository", "pipeline output repository"];
	add_config_section($fpPC, "wrapper_align", "wrap");
	add_config_parameters($fpPC, \%hParams);
	
	%hParams = ();
	
	$hParams{'INPUT_FILE'} = ["$sRpkmFileList", "path to list of rpkm coverage file"];
	add_config_section($fpPC, "expression_plots", "rpkm");
	add_config_parameters($fpPC, \%hParams);

	complete_component($oPL);
}
elsif (defined $hCmdLineOption{'alignment'}) {
	init_component($oPL, "serial");
		include_component_layout($oPL, $sTemplateDir, "wrapper_align", "wrap");
	complete_component($oPL);
	
	%hParams = ();
	$hParams{'PIPELINE_ID'} = ["\$;PIPELINEID\$;", "ergatis pipeline id"];
	$hParams{'OUTPUT_REPOSITORY'} = ["\$;REPOSITORY_ROOT\$;/output_repository", "pipeline output repository"];
	add_config_section($fpPC, "wrapper_align", "wrap");
	add_config_parameters($fpPC, \%hParams);
}

if (defined $hCmdLineOption{'diff_gene_expr'}) {
	die "Error! Comparison groups undefined !!!\n" if (!(defined $hCmdLineOption{'comparison_groups'}));
	
	init_component($oPL, "serial");
	
	if ((defined $hCmdLineOption{'count'}) || (defined $hCmdLineOption{'alignment'})) {
		die "Error! Annotation file undefined !!!\n" if (!(defined $hCmdLineOption{'gtffile'}));
		
		###	Add HTSeq Component & Parameters below ###
		init_component($oPL, "serial");
			if (! defined $hCmdLineOption{'alignment'}) {
				if (!((defined $hCmdLineOption{'sorted'}) && ($hCmdLineOption{'sorted'} =~ m/name/))) {
					include_component_layout($oPL, $sTemplateDir, "samtools_file_convert", "sorted_name") ;
				}
			}
			include_component_layout($oPL, $sTemplateDir, "htseq", "exon_counts");
		complete_component($oPL);
		
		if (! defined $hCmdLineOption{'alignment'}) {
			if (!((defined $hCmdLineOption{'sorted'}) && ($hCmdLineOption{'sorted'} =~ m/name/))) {
				%hParams = ();
				if ((defined $hCmdLineOption{'file_type'}) && ($hCmdLineOption{'file_type'} =~ m/BAM/i)) {
					$hParams{'INPUT_FILE_LIST'} = ["$sBamFileList", "path to list of alignment files"];
					$hParams{'INPUT_FILE_FORMAT'} = ["BAM", "input alignment file format (BAM or SAM)"];
					$nOpt = "13";
					$hParams{'OPTIONS'} = ["$nOpt", "string of options for file conversion (eg : 123). 1 - BAM to sorted BAM, 2 - sorted BAM to indexed BAM, 3 - BAM to SAM, and 4 - SAM to BAM"];
				}
				else {
					$hParams{'INPUT_FILE_LIST'} = ["$sSamFileList", "path to list of alignment files"];
					$hParams{'INPUT_FILE_FORMAT'} = ["SAM", "input alignment file format (BAM or SAM)"];
					$nOpt = "413";
					$hParams{'OPTIONS'} = ["$nOpt", "string of options for file conversion (eg : 123). 1 - BAM to sorted BAM, 2 - sorted BAM to indexed BAM, 3 - BAM to SAM, and 4 - SAM to BAM"];
				}
				$hParams{'SAMTOOLS_SORT_PARAMETERS'} = ["-n", "samtools sort parameters"];
				add_config_section($fpPC, "samtools_file_convert", "sorted_name");
				add_config_parameters($fpPC, \%hParams);
				
				$sSamFileList = '$;REPOSITORY_ROOT$;/output_repository/samtools_file_convert/$;PIPELINEID$;_sorted_name/samtools_file_convert.sorted_by_name_sam.list';
			}
		}
		
		%hParams = ();
		$hParams{'INPUT_FILE_LIST'} = ["$sSamFileList", "path to list of sorted-by-name alignment SAM files"];
		$hParams{'ANNOTATION_FILE'} = ["$hCmdLineOption{'gtffile'}", "path to annotation file in GFF or GTF format"];
		config2params(\%hParams, \%hConfig, 'htseq');
		add_config_section($fpPC, "htseq", "exon_counts");
		add_config_parameters($fpPC, \%hParams);
		
		$sCountsFileList = '$;REPOSITORY_ROOT$;/output_repository/htseq/$;PIPELINEID$;_exon_counts/htseq.counts.list';
	}
	
	###	Add DESeq Component & Parameters below ###
	init_component($oPL, "serial");
	        init_component($oPL,"parallel");
		      include_component_layout($oPL, $sTemplateDir, "deseq", "differential_expression");
	              include_component_layout($oPL, $sTemplateDir, "edgeR", "edgeR_diff_expression");
	
                      ##Add Deseq parameters
	              %hParams = ();
                      $hParams{'INPUT_FILE_LIST'} = ["$sOutDir/deseq_sample_info.list", "path to list of tab-delimited sample information files"];
	              $hParams{'LIST_FILE'} = ["$sCountsFileList", "path to list file of HTSeq alignment count files"];
	              config2params(\%hParams, \%hConfig, 'deseq');
	              add_config_section($fpPC, "deseq", "differential_expression");
	              add_config_parameters($fpPC, \%hParams);
                      ##Add EdgeR parameters
                      %hParams = ();
	              $hParams{'INPUT_FILE_LIST'} = ["$sOutDir/deseq_sample_info.list", "path to list of tab-delimited sample information files"];
	              $hParams{'LIST_FILE'} = ["$sCountsFileList", "path to list file of HTSeq alignment count files"];
	              config2params(\%hParams, \%hConfig, 'edgeR');
	              add_config_section($fpPC, "edgeR", "edgeR_diff_expression");
	              add_config_parameters($fpPC, \%hParams);
	        complete_component($oPL);

                include_component_layout($oPL, $sTemplateDir, "filter_deseq", "filter_de");
	        include_component_layout($oPL, $sTemplateDir, "expression_plots", "deseq");
                ##Add Deseq Filter component.
	        $Deseq_List = '$;REPOSITORY_ROOT$;/output_repository/deseq/$;PIPELINEID$;_differential_expression/deseq.table.list';
	        %hParams = ();
	        $hParams{'INPUT_FILE'} = [$Deseq_List,"path to output list file of deseq"];
	        config2params(\%hParams, \%hConfig, 'filter_deseq');
	        add_config_section($fpPC, "filter_deseq", "filter_de");
	        add_config_parameters($fpPC, \%hParams);
	
	        %hParams = ();

	        $hParams{'INPUT_FILE'} = [$Deseq_List,"path to output list file of deseq"];
	        config2params(\%hParams, \%hConfig, 'expression_plots');
	        add_config_section($fpPC, "expression_plots", "deseq");
	        add_config_parameters($fpPC, \%hParams);

        complete_component($oPL);

	complete_component($oPL);    

}


    

if ((defined $hCmdLineOption{'isoform_analysis'}) || (defined $hCmdLineOption{'diff_isoform_analysis'})) {
	init_component($oPL, "serial");
	
	if (defined $hCmdLineOption{'isoform_analysis'}) {
		###	Add Cufflinks Component & Parameters below ###
		init_component($oPL, "serial");
			if (! defined $hCmdLineOption{'alignment'}) {
				if (!((defined $hCmdLineOption{'sorted'}) && ($hCmdLineOption{'sorted'} =~ m/position/))) {
					include_component_layout($oPL, $sTemplateDir, "samtools_file_convert", "sorted_position");
				}
			}
			include_component_layout($oPL, $sTemplateDir, "cufflinks", "isoform");
		complete_component($oPL);
		
		if (! defined $hCmdLineOption{'alignment'}) {
			if (!((defined $hCmdLineOption{'sorted'}) && ($hCmdLineOption{'sorted'} =~ m/position/))) {
				%hParams = ();
				if ((defined $hCmdLineOption{'file_type'}) && ($hCmdLineOption{'file_type'} =~ m/SAM/i)) {
					$hParams{'INPUT_FILE_LIST'} = ["$sSamFileList", "path to list of alignment files"];
					$hParams{'INPUT_FILE_FORMAT'} = ["SAM", "input alignment file format (BAM or SAM)"];
					$nOpt = "412";
					$hParams{'OPTIONS'} = ["$nOpt", "string of options for file conversion (eg : 123). 1 - BAM to sorted BAM, 2 - sorted BAM to indexed BAM, 3 - BAM to SAM, and 4 - SAM to BAM"];
				}
				else {
					$hParams{'INPUT_FILE_LIST'} = ["$sBamFileList", "path to list of alignment files"];
					$hParams{'INPUT_FILE_FORMAT'} = ["BAM", "input alignment file format (BAM or SAM)"];
					$nOpt = "12";
					$hParams{'OPTIONS'} = ["$nOpt", "string of options for file conversion (eg : 123). 1 - BAM to sorted BAM, 2 - sorted BAM to indexed BAM, 3 - BAM to SAM, and 4 - SAM to BAM"];
				}
				add_config_section($fpPC, "samtools_file_convert", "sorted_position");
				add_config_parameters($fpPC, \%hParams);
				
				$sBamFileList = '$;REPOSITORY_ROOT$;/output_repository/samtools_file_convert/$;PIPELINEID$;_sorted_position/samtools_file_convert.sorted_by_position_bam.list';
			}
		}
		
		%hParams = ();
		$hParams{'INPUT_FILE_LIST'} = ["$sBamFileList", "path to list of sorted-by-position alignment SAM or BAM files"];
		$hParams{'ANNOTATION_FILE'} = ["$hCmdLineOption{'gtffile'}", "path to annotation file in GFF or GTF format"] if ((! defined $hCmdLineOption{'include_novel'}) && (defined $hCmdLineOption{'gtffile'}));
		config2params(\%hParams, \%hConfig, 'cufflinks');
		add_config_section($fpPC, "cufflinks", "isoform");
		add_config_parameters($fpPC, \%hParams);
		
		$sGtfFileList = '$;REPOSITORY_ROOT$;/output_repository/cufflinks/$;PIPELINEID$;_isoform/cufflinks.transcripts.gtf.list';
	}
	
	if (defined $hCmdLineOption{'diff_isoform_analysis'}) {
		die "Error! Comparison groups undefined !!!\n" if (!(defined $hCmdLineOption{'comparison_groups'}));
		
		init_component($oPL, "serial");
		
		###	Add CuffCompare Component & Parameters below ###
		if ((! defined $hCmdLineOption{'alignment'}) && (! defined $hCmdLineOption{'isoform_analysis'})) {
			if (!((defined $hCmdLineOption{'sorted'}) && ($hCmdLineOption{'sorted'} =~ m/position/))) {
				include_component_layout($oPL, $sTemplateDir, "samtools_file_convert", "sorted_position");
			}
			else {
				if ((defined $hCmdLineOption{'file_type'}) && ($hCmdLineOption{'file_type'} =~ m/SAM/i)) {
					include_component_layout($oPL, $sTemplateDir, "samtools_file_convert", "sorted_position");
				}
			}
		}
		include_component_layout($oPL, $sTemplateDir, "create_cuffsuite_files", "cuffcompare");
		include_component_layout($oPL, $sTemplateDir, "cuffcompare", "comparison");
		
		if ((! defined $hCmdLineOption{'alignment'}) && (! defined $hCmdLineOption{'isoform_analysis'})) {
			if (!((defined $hCmdLineOption{'sorted'}) && ($hCmdLineOption{'sorted'} =~ m/position/))) {
				%hParams = ();
				if ((defined $hCmdLineOption{'file_type'}) && ($hCmdLineOption{'file_type'} =~ m/SAM/i)) {
					$hParams{'INPUT_FILE_LIST'} = ["$sSamFileList", "path to list of alignment files"];
					$hParams{'INPUT_FILE_FORMAT'} = ["SAM", "input alignment file format (BAM or SAM)"];
					$nOpt = "412";
					$hParams{'OPTIONS'} = ["$nOpt", "string of options for file conversion (eg : 123). 1 - BAM to sorted BAM, 2 - sorted BAM to indexed BAM, 3 - BAM to SAM, and 4 - SAM to BAM"];
				}
				else {
					$hParams{'INPUT_FILE_LIST'} = ["$sBamFileList", "path to list of alignment files"];
					$hParams{'INPUT_FILE_FORMAT'} = ["BAM", "input alignment file format (BAM or SAM)"];
					$nOpt = "12";
					$hParams{'OPTIONS'} = ["$nOpt", "string of options for file conversion (eg : 123). 1 - BAM to sorted BAM, 2 - sorted BAM to indexed BAM, 3 - BAM to SAM, and 4 - SAM to BAM"];
				}
				add_config_section($fpPC, "samtools_file_convert", "sorted_position");
				add_config_parameters($fpPC, \%hParams);
				
				$sBamFileList = '$;REPOSITORY_ROOT$;/output_repository/samtools_file_convert/$;PIPELINEID$;_sorted_position/samtools_file_convert.sorted_by_position_bam.list';
			}
			else {
				if ((defined $hCmdLineOption{'file_type'}) && ($hCmdLineOption{'file_type'} =~ m/SAM/i)) {
					%hParams = ();
					$hParams{'INPUT_FILE_LIST'} = ["$sSamFileList", "path to list of alignment files"];
					$hParams{'INPUT_FILE_FORMAT'} = ["SAM", "input alignment file format (BAM or SAM)"];
					$nOpt = "4";
					$hParams{'OPTIONS'} = ["$nOpt", "string of options for file conversion (eg : 123). 1 - BAM to sorted BAM, 2 - sorted BAM to indexed BAM, 3 - BAM to SAM, and 4 - SAM to BAM"];
					add_config_section($fpPC, "samtools_file_convert", "sorted_position");
					add_config_parameters($fpPC, \%hParams);
					
					$sBamFileList = '$;REPOSITORY_ROOT$;/output_repository/samtools_file_convert/$;PIPELINEID$;_sorted_position/samtools_file_convert.sorted_by_position_bam.list';
				}
			}
		}
		
		%hParams = ();
		$hParams{'SAMPLE_INFO'} = ["$hCmdLineOption{'sample_file'}", "path to sample info file with information on all samples to be analyzed"];
		$hParams{'CUFF_PROG'} = ["Cuffcompare", "Cuffsuite program (Cuffcompare or Cuffdiff) to create files for"];
		$hParams{'GRP_COMP'} = ["$hCmdLineOption{'comparison_groups'}", "string of groups to compare. e.g. \"GRP#2vsGRP#1,GRP#3vsGRP#1\""];
		$hParams{'GTF_LISTFILE'} = ["$sGtfFileList", "path to list file with information on all cufflinks or cuffcompare GTF files"];
		add_config_section($fpPC, "create_cuffsuite_files", "cuffcompare");
		add_config_parameters($fpPC, \%hParams);
		
		$sCuffFileList = '$;REPOSITORY_ROOT$;/output_repository/create_cuffsuite_files/$;PIPELINEID$;_cuffcompare/cuffsuite_input_file.list';
		
		%hParams = ();
		$hParams{'INPUT_FILE_LIST'} = ["$sCuffFileList", "path to list of Cufflinks GTF files or list of sets of Cufflinks GTF files"];
		$hParams{'ANNOTATION_FILE'} = ["$hCmdLineOption{'gtffile'}", "path to annotation file in GFF or GTF format"] if (defined $hCmdLineOption{'gtffile'});
		config2params(\%hParams, \%hConfig, 'cuffcompare');
		add_config_section($fpPC, "cuffcompare", "comparison");
		add_config_parameters($fpPC, \%hParams);
		
		$sGtfFileList = '$;REPOSITORY_ROOT$;/output_repository/cuffcompare/$;PIPELINEID$;_comparison/cuffcompare.combined.gtf.list';
		
		###	Add CuffDiff Component & Parameters below ###
		include_component_layout($oPL, $sTemplateDir, "create_cuffsuite_files", "cuffdiff");
		include_component_layout($oPL, $sTemplateDir, "cuffdiff", "differential_expression");
		include_component_layout($oPL, $sTemplateDir, "cuffdiff_filter", "filter_cuff");

		%hParams = ();
		$hParams{'SAMPLE_INFO'} = ["$hCmdLineOption{'sample_file'}", "path to sample info file with information on all samples to be analyzed"];
		$hParams{'CUFF_PROG'} = ["Cuffdiff", "Cuffsuite program (Cuffcompare or Cuffdiff) to create files for"];
		$hParams{'GRP_COMP'} = ["$hCmdLineOption{'comparison_groups'}", "string of groups to compare. e.g. \"GRP#2vsGRP#1,GRP#3vsGRP#1\""];
		if (defined $hCmdLineOption{'use_ref_gtf'}) {
			die "Error! Annotation file undefined !!!\n" if (!(defined $hCmdLineOption{'gtffile'}));
			$hParams{'GTFFILE'} = ["$hCmdLineOption{'gtffile'}", "path to single reference GTF file to be used with all cuffdiff comparisons"];
		}
		else {
			$hParams{'GTF_LISTFILE'} = ["$sGtfFileList", "path to list file with information on all cufflinks or cuffcompare GTF files"];
		}
		$hParams{'SAM_LISTFILE'} = ["$sBamFileList", "path to list info file with information on all alignment SAM files sorted by position"];
		add_config_section($fpPC, "create_cuffsuite_files", "cuffdiff");
		add_config_parameters($fpPC, \%hParams);
		
		$sCuffFileList = '$;REPOSITORY_ROOT$;/output_repository/create_cuffsuite_files/$;PIPELINEID$;_cuffdiff/cuffsuite_input_file.list';
		
		%hParams = ();
		$hParams{'INPUT_FILE_LIST'} = ["$sCuffFileList", "path to list of tab-delimited lists of aligmment SAM file(s) and respective annotation GTF files"];
		config2params(\%hParams, \%hConfig, 'cuffdiff');
		add_config_section($fpPC, "cuffdiff", "differential_expression");
		add_config_parameters($fpPC, \%hParams);

		$Cuff_List = '$;REPOSITORY_ROOT$;/output_repository/cuffdiff/$;PIPELINEID$;_differential_expression/cuffdiff.isoform.diff.list';
		%hParams = ();
		$hParams{'INPUT_FILE'} = [$Cuff_List,"path to output list file of cuffdiff"];
		config2params(\%hParams, \%hConfig, 'cuffdiff_filter');
		add_config_section($fpPC, "cuffdiff_filter", "filter_cuff");
		add_config_parameters($fpPC, \%hParams);
		
		complete_component($oPL);

		init_component($oPL, "parallel");
		include_component_layout($oPL, $sTemplateDir, "expression_plots", "cuffdiff");
		include_component_layout($oPL, $sTemplateDir, "cummerbund", "cummer_cuff");
		%hParams = ();
                
	        $hParams{'INPUT_FILE'} = [$Cuff_List,"path to output list file of cuffdiff"];
	        config2params(\%hParams, \%hConfig, 'expression_plots');
	        add_config_section($fpPC, "expression_plots", "cuffdiff");
	        add_config_parameters($fpPC, \%hParams);

		%hParams = ();
                
	        $hParams{'INPUT_FILE'} = [$Cuff_List,"path to output list file of cuffdiff"];
	        config2params(\%hParams, \%hConfig, 'cummerbund');
	        add_config_section($fpPC, "cummerbund", "cummer_cuff");
	        add_config_parameters($fpPC, \%hParams);

		
		complete_component($oPL);

	}
	
	complete_component($oPL);
}

if ( (defined $hCmdLineOption{'diff_gene_expr'}) || (defined $hCmdLineOption{'visualization'}) || (defined $hCmdLineOption{'rpkm_analysis'}) ||
	 ((defined $hCmdLineOption{'isoform_analysis'}) || (defined $hCmdLineOption{'diff_isoform_analysis'})) ) {
	complete_component($oPL);
}

# Complete pipeline.layout
complete_pipeline_layout($oPL);

# End XML::Writer object
$oPL->end();

close($fpPL);
close($fpPC);

($bDebug || $bVerbose) ? print STDERR "\nGenerating $sPLayout & $sPConfig ..... done\n" : undef;

($bDebug || $bVerbose) ? print STDERR "\nProcessing $hCmdLineOption{'sample_file'} ..... done\n" : undef;

if (defined $hCmdLineOption{'repository_root'}) {
	($bDebug || $bVerbose) ? print STDERR "\nInitiation of pipeline on ergatis .....\n" : undef;
	
	$sCmd = $RealBin."/run_rnaseq_pipeline.pl".
			" --layout ".$sPLayout.
			" --config ".$sPConfig.
			" --repository_root ".$hCmdLineOption{'repository_root'};
	
	$sCmd .= " --ergatis_config ".$hCmdLineOption{'ergatis_ini'};
	
	exec_command($sCmd);
	
	($bDebug || $bVerbose) ? print STDERR "\nInitiation of pipeline on ergatis ..... done\n" : undef;
}

################################################################################
### Subroutines
################################################################################

sub check_parameters {
    my $phOptions		= shift;
    my $paComponents	= shift;
    
    ## make sure input parameters are provided
    if ((!(defined $phOptions->{'sample_file'})) || 
    	(!(defined $phOptions->{'config_file'}))) {
    	print STDERR "\nError! Both Sample File and Config File are required!\n";
		pod2usage( -msg => $sHelpHeader, -exitval => 1);
	}
	
	$phOptions->{'reffile'} = File::Spec->rel2abs($phOptions->{'reffile'}) if (defined $phOptions->{'reffile'});
	$phOptions->{'gtffile'} = File::Spec->rel2abs($phOptions->{'gtffile'}) if (defined $phOptions->{'gtffile'});
	$phOptions->{'bwtidxfile'} = File::Spec->rel2abs($phOptions->{'bwtidxfile'}) if (defined $phOptions->{'bwtidxfile'});
	
	my ($nI, $bFlag);
	
	$bFlag = 0;
	
	for ($nI = 0; $nI < @{$paComponents}; $nI++) {
		$bFlag = 1 if (defined $phOptions->{$paComponents->[$nI]});
	}
	
	if (!($bFlag)) {
		print STDERR "\nAnalysis Component Options :\n";
		for ($nI = 0; $nI < @{$paComponents}; $nI++) {
			print STDERR "\t--$paComponents->[$nI]\n";
		}
		die "\nSelect one or more of the above options and specify the required inputs ...\n";
	}
	
	## handle some defaults
    $phOptions->{'quality'} = 33 if (! (defined $phOptions->{'quality'}) );
    
    return;
}

sub is_gzipped {
    my $phOptions	= shift;
    my $sFile		= shift;
    
    ## make sure config file and config hash are provided
    if ((!(defined $phOptions)) || 
    	(!(defined $sFile))) {
    	die "Error! In subroutine is_gzipped: Incomplete parameter list !!!\n";
	}
	
	my ($sDir, $sLinkFile, $bFlag);
	my ($sCmd, $sMsg);
	
	$bFlag = 0;
	
	if ( -l $sFile ) {
		$sLinkFile = readlink $sFile;
		$sLinkFile =~ s/\s+$//;
		
		($_, $sDir, $_) = File::Spec->splitpath($sFile);
		$sLinkFile = File::Spec->rel2abs($sLinkFile, $sDir);
		$bFlag = is_gzipped($phOptions, $sLinkFile);
	}
	else {
		$sCmd = "file $sFile";
		$sMsg = `$sCmd`;
		if ($sMsg =~ m/gzip/i) {
			$bFlag = 1;
		}
	}
    
    return $bFlag;
}

sub exec_command {
	my $sCmd = shift;
	
	if ((!(defined $sCmd)) || ($sCmd eq "")) {
		die "\nSubroutine::exec_command : ERROR! Incorrect command!\n";
	}
	
	my $nExitCode;
	
	#print STDERR "$sCmd\n";
	$nExitCode = system("$sCmd");
	if ($nExitCode != 0) {
		die "\tERROR! Command Failed!\n\t$!\n";
	}
	print STDERR "\n";
}


sub read_config {
    my $phOptions	= shift;
    my $phConfig	= shift;
    
    ## make sure config file and config hash are provided
    if ((!(defined $phOptions->{'config_file'})) || 
    	(!(defined $phConfig))) {
    	die "Error! In subroutine read_config: Incomplete parameter list !!!\n";
	}
	
	my ($sConfigFile);
	my ($sComponent, $sParam, $sValue, $sDesc);
	my ($fpCFG);
	
	$sConfigFile = $RealBin."/Eukaryotic.RNASeq.pipeline.config";
	if (defined $phOptions->{'config_file'}) {
		$sConfigFile = $phOptions->{'config_file'};
	}
	open($fpCFG, "<$sConfigFile") or die "Error! Cannot open $sConfigFile for reading: $!";
	
	$sComponent = $sParam = $sValue = $sDesc = "";
	while (<$fpCFG>) {
		$_ =~ s/\s+$//;
		next if ($_ =~ /^#/);
		next if ($_ =~ /^$/);
		
		if ($_ =~ m/^\[(\S+)\]$/) {
			$sComponent = $1;
			next;
		}
		elsif ($_ =~ m/^;;\s*(.*)/) {
			$sDesc .= "$1.";
			next;
		}
		elsif ($_ =~ m/\$;(\S+)\$;\s*=\s*(.*)/) {
			$sParam = $1;
			$sValue = $2;
			
			if ((defined $sValue) && ($sValue !~ m/^\s*$/)) {
				$phConfig->{$sComponent}{$sParam} = ["$sValue", "$sDesc"];
			}
			
			$sParam = $sValue = $sDesc = "";
			next;
		}
	}
	
	close($fpCFG);
	    
    return;
}

sub config2params {
    my $phParams	= shift;
    my $phConfig	= shift;
    my $sComponent	= shift;
    
    ## make sure config file and config hash are provided
    if ((!(defined $phParams)) || 
    	(!(defined $phConfig)) || 
    	(!(defined $sComponent))) {
    	die "Error! In subroutine config2params: Incomplete parameter list !!!\n";
	}
	
	my ($sParam, $sValue, $sDesc);
	
	foreach $sParam (keys %{$phConfig->{$sComponent}}) {
		$sValue = $phConfig->{$sComponent}{$sParam}[0];
		$sDesc  = $phConfig->{$sComponent}{$sParam}[1];
		
		$phParams->{"$sParam"} = ["$sValue", "$sDesc"];
	}
    
    return;
}

sub init_pipeline_layout {
    my $oLayout		= shift;
    
    ## make sure XML::Writer object provided
    if (!(defined $oLayout)) {
		die "Error! In subroutine init_pipeline_layout: Incomplete parameter list !!!\n";
	}
    
    $oLayout->startTag("commandSetRoot", "xmlns:xsi" => "http://www.w3.org/2001/XMLSchema-instance", 
    				   "xsi:schemaLocation" => "commandSet.xsd", "type" => "instance");
    
    $oLayout->startTag("commandSet", "type" => "serial" );
    
    $oLayout->dataElement("state", "incomplete");
    
    $oLayout->dataElement("name", "start pipeline:");
    
    return;
}

sub complete_pipeline_layout {
	my $oLayout		= shift;
    
    ## make sure XML::Writer object provided
    if (!(defined $oLayout)) {
		die "Error! In subroutine complete_pipeline_layout: Incomplete parameter list !!!\n";
	}
    
    $oLayout->endTag("commandSet");
    
    $oLayout->endTag("commandSetRoot");
    
    return;
}

sub init_component {
	my $oLayout		= shift;
	my $sType		= shift;
	
	## make sure XML::Writer object and command set type provided
    if ((!(defined $oLayout)) || 
    	(!(defined $sType)) ) {
		die "Error! In subroutine init_component: Incomplete parameter list !!!\n";
	}
	
	$oLayout->startTag("commandSet", "type" => "$sType" );
	
	$oLayout->dataElement("state", "incomplete");
    
    return;
}

sub include_component_layout {
	my $oLayout		= shift;
	my $sTDir		= shift;
	my $sComponent	= shift;
	my $sToken		= shift;
	
	## make sure XML::Writer object and component name provided
    if ((!(defined $oLayout)) || 
    	(!(defined $sComponent)) ) {
		die "Error! In subroutine add_component: Incomplete parameter list !!!\n";
	}
	
	$sToken = "default" if (!(defined $sToken));
	
	my $sCLayout = "$sTDir/$sComponent.$sToken.pipeline.layout";
	
	die "Error! Cannot find $sCLayout !!!\n"
		if (! -e "$sCLayout");
	
	$oLayout->emptyTag("INCLUDE", 'file' => "$sCLayout");
    
    return;
}

sub complete_component {
	my $oLayout		= shift;
    
    ## make sure XML::Writer object provided
    if (!(defined $oLayout)) {
		die "Error! In subroutine complete_component: Incomplete parameter list !!!\n";
	}
	
	$oLayout->endTag("commandSet");
    
    return;
}

sub add_config_section {
	my $fpConfig		= shift;
	my $sComponent		= shift;
	my $sToken			= shift;
	
	## make sure XML::Writer object and component name provided
    if ((!(defined $fpConfig)) || 
    	(!(defined $sComponent)) ) {
		die "Error! In subroutine add_config_section: Incomplete parameter list !!!\n";
	}
	
	$sToken = "default" if (!(defined $sToken));
	
	print $fpConfig "[".$sComponent." ".$sToken."]\n";
    
    return;
}

sub add_config_parameters {
	my $fpConfig		= shift;
	my $phParams		= shift;
	
	## make sure XML::Writer object and component name provided
    if ((!(defined $fpConfig)) || 
    	(!(defined $phParams)) ) {
		die "Error! In subroutine add_config_parameters: Incomplete parameter list !!!\n";
	}
	
	my ($sKey, $sValue, $sDesc);
	
	foreach $sKey (sort keys %{$phParams}) {
		$sValue = $phParams->{$sKey}[0];
		$sDesc = $phParams->{$sKey}[1];
		
		print $fpConfig ";; $sDesc\n".
						"\$;".$sKey."\$;".
						" = ".$sValue."\n";
	}
	
	print $fpConfig "\n";
    
    return;
}

################################################################################
