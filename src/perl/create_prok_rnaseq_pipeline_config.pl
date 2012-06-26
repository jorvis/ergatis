#!/usr/bin/perl

################################################################################
### POD Documentation
################################################################################

=head1 NAME

create_prok_rnaseq_pipeline_config.pl - Creates the pipeline.layout and pipeline.config for the
                                        automated prokaryotic rna-seq pipeline

=head1 SYNOPSIS

    create_prok_rnaseq_pipeline_config.pl --s <samples_file> --c <config_file> [--r <reference_fasta>] 
                                          [--qual <quality_score_format>] [--gtf <annotation_file>] 
                                          [--bowtie_build] [--quality_stats] [--quality_trimming] 
                                          [--alignment] [--bwtidxfile <bowtie_index>] [--visualization] 
                                          [--diff_gene_expr] [--comparison_groups <str>] [--count]  
                                          [--file_type <SAM|BAM>] [--sorted <position|name>] 
                                          [--td <template_directory>] [--o <outdir>] [--v] 
                                          [--man] [--help]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --s <samples_file>               = /path/to/samples file with information on all samples to be analyzed.

    --c <config_file>                = /path/to/config file with parameter information for multiple components.

    --r <reference_fasta>            = /path/to/reference FastA file for all samples. Optional.
    
    --qual <quality_score_format>    = FastQ quality score format (33 or 64). Optional. [33]

    --gtf <annotation_file>          = /path/to/annotation file in GFF or GTF format.
                                       Required with '--diff_gene_expr'.

    --bowtie_build                   = execute bowtie_build component. Requires '--r'.

    --quality_stats                  = execute fastx_quality_stats component.

    --quality_trimming               = execute fastx_trimming component. Also generates quality statistics.

    --alignment                      = execute bowtie alignment component.
                                       Sample file should be in the following format
                                       #Sample_ID<tab>Group_ID<tab>Mate_Pair_1<tab>Mate_Pair_2

        --bwtidxfile <bowtie_index>  = /path/to/bowtie index file for alignment of all samples.
                                       Required for '--alignment' if not specifying '--bowtie_build'.

    --visualization                   = execute bam2bigwig component.
                                        Requires additional information in sample file if not specifying '--alignment'.
                                        Sample file should be in the following format
                                        #Sample_ID<tab>Group_ID<tab>Alignment_BAM_File

        --file_type <SAM|BAM>         = alignment file format (BAM or SAM). [BAM]
                                        Required if not specifying '--alignment' and providing alignment file information in sample file.

        --sorted <position>           = if alignment BAM/SAM file is already sorted by position. [undef]

    --diff_gene_expr                 = execute differential gene expression analysis component.
                                       Requires additional information in sample file if not specifying '--alignment'.
                                       Sample file should be in the following format
                                       #Sample_ID<tab>Group_ID<tab>Alignment_File

        --comparison_groups <str>    = string of groups to compare. e.g. "GRP#2vsGRP#1,GRP#3vsGRP#1".
                                       Group Ids SHOULD match Group Names in column 2 of the sample info file

        --count                      = generate count files.
                                       Required if not specifying '--alignment' and 
                                       providing alignment SAM file information in sample file.

        --file_type <SAM|BAM>        = alignment file format (BAM or SAM).
                                       Required if specifying '--count'. [SAM]

        --sorted <name>              = if alignment SAM file is already sorted by name. [undef]

    --td <template_directory>        = /path/to/template directory. Optional. [present working directory]

    --o <output dir>                 = /path/to/output directory. Optional. [present working directory]

    --v                              = generate runtime messages. Optional

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
    bowtie                   : generates the Bowtie alignment files for single-end or paired-end sequence(s).
    samtools_file_convert    : converts file formats for downstream analysis.
    samtools_alignment_stats : generates the alignment stats from the alignment BAM file(s).
    bam2bigwig               : converts the alignment BAM file(s) to BedGraph and BigWig file(s).
    htseq                    : generates the count files from the alignment SAM file(s) sorted by name.
    deseq                    : generates the differential gen expression analysis results utilizing DESeq software.
 
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

use constant VERSION => '1.0.0';
use constant PROGRAM => eval { ($0 =~ m/(\w+\.pl)$/) ? $1 : $0 };

##############################################################################
### Globals
##############################################################################

my @aComponents = ("bowtie_build", "quality_stats", "quality_trimming", "alignment", "visualization",
				   "diff_gene_expr");
my %hCmdLineOption = ();
my $sHelpHeader = "\nThis is ".PROGRAM." version ".VERSION."\n";

GetOptions( \%hCmdLineOption,
			'sample_file|s=s', 'config_file|c=s', 'reffile|r=s', 'quality|qual=i', 'gtffile|gtf=s',
			'bowtie_build', 'quality_stats', 'quality_trimming', 'alignment', 'bwtidxfile=s', 
			'visualization', 
			'diff_gene_expr', 'comparison_groups=s', 'count', 'file_type=s', 'sorted=s', 
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
my ($fpPL, $fpPC, $fpLST1, $fpLST2, $fpLST, $fpSMPL);
my ($sSampleName, $sGroupName, $sRead1File, $sRead2File, @aReadFiles, $sList);
my ($sSamRefFile, $sBamFileList, $sSamFileList, $sMapStatsList, $sCountsFileList);
my ($sFeature, $sAttrID);
my (@aComparisons, $sCGrp, $sGrpX, $sGrpY);
my ($sList1File, $sList2File, $sListFile, $sFile, $sInFile);
my ($sTimeStamp, $sCmd, $sArgs, $nOpt, $nI, $bPE, $bGZ);
my ($nSec, $nMin, $nHour, $nMDay, $nMon, $nYear, $nWDay, $nYDay, $bDST);
my $bDebug   = (defined $hCmdLineOption{'debug'}) ? TRUE : FALSE;
my $bVerbose = (defined $hCmdLineOption{'verbose'}) ? TRUE : FALSE;

################################################################################
### Main
################################################################################

($bDebug || $bVerbose) ? print STDERR "\nInitiating Prokaryotic RNA-Seq Analysis Pipeline .....\n" : undef;
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

$sTemplateDir = $RealBin."/Prokaryotic_RNA_Seq_Analysis";
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
if ((defined $hCmdLineOption{'alignment'})) {
	
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
		
		die "Error! Missing Reads 1 File !!!\n" if (! (defined $sRead1File));
		
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
elsif (defined $hCmdLineOption{'visualization'}) {
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
	init_component($oPL, "serial");
		init_component($oPL, "parallel");
			include_component_layout($oPL, $sTemplateDir, "fastx_quality_stats", "read1");
			include_component_layout($oPL, $sTemplateDir, "fastx_quality_stats", "read2") if ($bPE);
		complete_component($oPL);
	complete_component($oPL);
	
	%hParams = ();
	$hParams{'INPUT_FILE_LIST'} = ["$sList1File", "path to list of input FastQ/A sequence files"];
	$hParams{'QUALITY_STRING'} = ["$hCmdLineOption{'quality'}", "quality string type for FastQ files (33 or 64)"];
	
	$sArgs = "";
	$sArgs = $hConfig{'fastx_quality_stats'}{'OTHER_ARGS'}[0] if (defined $hConfig{'fastx_quality_stats'}{'OTHER_ARGS'});
	$sArgs .= " --z" if (($bGZ) && ($sArgs !~ m/--z/));
	$hConfig{'fastx_quality_stats'}{'OTHER_ARGS'}[0] = $sArgs;
	config2params(\%hParams, \%hConfig, 'fastx_quality_stats');
	$hConfig{'fastx_quality_stats'}{'OTHER_ARGS'}[1] = "Additional Arguments" if (! defined $hConfig{'fastx_quality_stats'}{'OTHER_ARGS'});
	add_config_section($fpPC, "fastx_quality_stats", "read1");
	add_config_parameters($fpPC, \%hParams);
	
	if ($bPE) {
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
	init_component($oPL, "serial");
		init_component($oPL, "parallel");
			include_component_layout($oPL, $sTemplateDir, "fastx_trimming", "read1");
			include_component_layout($oPL, $sTemplateDir, "fastx_trimming", "read2") if ($bPE);;
		complete_component($oPL);
	complete_component($oPL);
	
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

if (defined $hCmdLineOption{'alignment'}) {
	init_component($oPL, "serial");
	
	if ($bPE) {
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
	
	include_component_layout($oPL, $sTemplateDir, "bowtie", "alignment");
	
	%hParams = ();
	$hParams{'INPUT_FILE_LIST'} = ["$sListFile", "path to list file consisting of tab separated first mate and second mate sequence files"];
	$hParams{'BOWTIE_INDEX_DIR'} = ["$sBwtIndexDir", "path to bowtie package binary directory"];
	$hParams{'BOWTIE_INDEX_PREFIX'} = ["$sBwtIndexPrefix", "bowtie index prefix"];
	
	$sArgs = "";
	$sArgs = $hConfig{'bowtie'}{'OTHER_PARAMETERS'}[0] if (defined $hConfig{'bowtie'}{'OTHER_PARAMETERS'});
	$sArgs .= " --sam" if (($sArgs !~ m/--sam/));
	$sArgs .= " --solexa1.3-quals" if (($hCmdLineOption{'quality'} == 64) && ($sArgs !~ m/--solexa1.3-quals/));
	$hConfig{'bowtie'}{'OTHER_PARAMETERS'}[0] = $sArgs;
	$hConfig{'bowtie'}{'OTHER_PARAMETERS'}[1] = "Additional Parameters" if (! defined $hConfig{'bowtie'}{'OTHER_PARAMETERS'});
	
	$sArgs  = "";
	$sArgs = $hConfig{'bowtie'}{'OTHER_ARGS'}[0] if (defined $hConfig{'bowtie'}{'OTHER_ARGS'});
	$sArgs .= " --gzip" if (($bGZ) && ($sArgs !~ m/--gzip/));
	$hConfig{'bowtie'}{'OTHER_ARGS'}[0] = $sArgs;
	$hConfig{'bowtie'}{'OTHER_ARGS'}[1] = "Additional Arguments" if (! defined $hConfig{'bowtie'}{'OTHER_ARGS'});
	config2params(\%hParams, \%hConfig, 'bowtie');
	add_config_section($fpPC, "bowtie", "alignment");
	add_config_parameters($fpPC, \%hParams);
	
	$sListFile = '$;REPOSITORY_ROOT$;/output_repository/bowtie/$;PIPELINEID$;_alignment/bowtie.sam.list';
	
	complete_component($oPL);
	
	init_component($oPL, "serial");
		
	init_component($oPL, "parallel");
		include_component_layout($oPL, $sTemplateDir, "samtools_file_convert", "sorted_position");
		include_component_layout($oPL, $sTemplateDir, "samtools_file_convert", "sorted_name");
	complete_component($oPL);
	
	%hParams = ();
	$hParams{'INPUT_FILE_LIST'} = ["$sListFile", "path to list of alignment files"];
	$hParams{'INPUT_FILE_FORMAT'} = ["SAM", "input alignment file format (BAM or SAM)"];
	$nOpt = "412";
	$hParams{'OPTIONS'} = ["$nOpt", "string of options for file conversion (eg : 123). 1 - BAM to sorted BAM, 2 - sorted BAM to indexed BAM, 3 - BAM to SAM, and 4 - SAM to BAM"];
	add_config_section($fpPC, "samtools_file_convert", "sorted_position");
	add_config_parameters($fpPC, \%hParams);
	
	$sBamFileList = '$;REPOSITORY_ROOT$;/output_repository/samtools_file_convert/$;PIPELINEID$;_sorted_position/samtools_file_convert.sorted_by_position_bam.list';
	
	%hParams = ();
	$hParams{'INPUT_FILE_LIST'} = ["$sListFile", "path to list of alignment files"];
	$hParams{'INPUT_FILE_FORMAT'} = ["SAM", "input alignment file format (BAM or SAM)"];
	$hParams{'SAMTOOLS_SORT_PARAMETERS'} = ["-n", "samtools sort parameters"];
	$nOpt = "413";
	$hParams{'OPTIONS'} = ["$nOpt", "string of options for file conversion (eg : 123). 1 - BAM to sorted BAM, 2 - sorted BAM to indexed BAM, 3 - BAM to SAM, and 4 - SAM to BAM"];
	add_config_section($fpPC, "samtools_file_convert", "sorted_name");
	add_config_parameters($fpPC, \%hParams);
	
	$sSamFileList = '$;REPOSITORY_ROOT$;/output_repository/samtools_file_convert/$;PIPELINEID$;_sorted_name/samtools_file_convert.sorted_by_name_sam.list';
	
	include_component_layout($oPL, $sTemplateDir, "samtools_alignment_stats", "alignment_stats");
	
	%hParams = ();
	$hParams{'INPUT_FILE_LIST'} = ["$sBamFileList", "path to list of alignment BAM files"];
	add_config_section($fpPC, "samtools_alignment_stats", "alignment_stats");
	add_config_parameters($fpPC, \%hParams);
	
	$sMapStatsList = '$;REPOSITORY_ROOT$;/output_repository/samtools_alignment_stats/$;PIPELINEID$;_alignment_stats/samtools_alignment_stats.mapstats.list';
	
	complete_component($oPL);
}

if ( (defined $hCmdLineOption{'diff_gene_expr'}) || (defined $hCmdLineOption{'visualization'}) ) {
	init_component($oPL, "parallel");
}

if (defined $hCmdLineOption{'visualization'}) {
	init_component($oPL, "serial");
	
	die "Error! Annotation file undefined !!!\n" if (!(defined $hCmdLineOption{'reffile'}));
	
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

if (defined $hCmdLineOption{'diff_gene_expr'}) {
	die "Error! Comparison groups undefined !!!\n" if (!(defined $hCmdLineOption{'comparison_groups'}));
	
	init_component($oPL, "serial");
	
	if ((defined $hCmdLineOption{'count'}) || (defined $hCmdLineOption{'alignment'})) {
		die "Error! Annotation file undefined !!!\n" if (!(defined $hCmdLineOption{'gtffile'}));
		
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
	
	init_component($oPL, "serial");
		include_component_layout($oPL, $sTemplateDir, "deseq", "differential_expression");
	complete_component($oPL);
	
	%hParams = ();
	$hParams{'INPUT_FILE_LIST'} = ["$sOutDir/deseq_sample_info.list", "path to list of tab-delimited sample information files"];
	$hParams{'LIST_FILE'} = ["$sCountsFileList", "path to list file of HTSeq alignment count files"];
	config2params(\%hParams, \%hConfig, 'deseq');
	add_config_section($fpPC, "deseq", "differential_expression");
	add_config_parameters($fpPC, \%hParams);
	
	complete_component($oPL);
}

if ( (defined $hCmdLineOption{'diff_gene_expr'}) || (defined $hCmdLineOption{'visualization'}) ) {
	complete_component($oPL, "parallel");
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
		print STDERR "\nSelect one or more of the above options and specify the required inputs ...\n";
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
	
	print STDERR "$sCmd\n";
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
	
	$sConfigFile = $RealBin."/Prokaryotic.RNASeq.pipeline.config";
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
