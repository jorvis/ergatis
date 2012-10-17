
#!/usr/bin/perl 

################################################################################
### POD Documentation
################################################################################

=head1 NAME

wrapper.pl          -script to combines various statistics output and creates a PDF.

=head1 SYNOPSIS

    wrapper.pl      --p pipeline_id --r out_rep 
                    [--o outdir] [--v] [--help]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS


    --r <out rep>         = Path to the output repository.

    --c <pipeline id>     = ID of the pipeline whose summary is required 

    --o <output dir>      = /path/to/output directory. Optional. [present working directory]

    --v                   = generate runtime messages. Optional
 
    --help                = Help document.

=head1 DESCRIPTION

The script generates summary of the alignment statistics.

=head1 AUTHOR

 Priti Kumari
 Bioinformatics Software Engineer 
 Institute for Genome Sciences
 University of Maryland
 Baltimore, Maryland 21201

=cut

################################################################################


###############################################################################
###Modules
###############################################################################
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
use File::Spec;
use lib "/local/projects-dr/ifx_core/package-tcreasy/lib/perl5";
use PDF::Create;

##############################################################################
###Constants
##############################################################################

use constant FALSE => 0;
use constant TRUE  => 1;

use constant PROGRAM => eval { ($0 =~ m/(\w+\.pl)$/) ? $1 : $0 };

###############################################################################
##Globals
###############################################################################

my %hCmdLineOption =();
my $sHelpHeader = "\nThis is ".PROGRAM."\n";
my ($sOutDir,$prefix,$path, $f1,$read,$a,$val);
my ($filehandle,$key,$fout,$file,$flag,$ffile);
my (%sample,%final);
my @arr;

################################################################################
### Main
################################################################################

GetOptions( \%hCmdLineOption,'outrep|r=s','pipeline|p=s',
            'outdir|o=s',
	    'verbose|v','help','man') or pod2usage(2);

## display documentation
pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};
pod2usage( -msg => $sHelpHeader, -exitval => 1) if $hCmdLineOption{'help'};

## Check for parameters.
check_parameters(\%hCmdLineOption);



my $bDebug   = (defined $hCmdLineOption{'debug'}) ? TRUE : FALSE;
my $bVerbose = (defined $hCmdLineOption{'verbose'}) ? TRUE : FALSE;


$sOutDir = File::Spec->curdir();
if (defined $hCmdLineOption{'outdir'}) {
    $sOutDir = $hCmdLineOption{'outdir'};

    if (! -e $sOutDir) {
        mkdir($hCmdLineOption{'outdir'}) ||
            die "ERROR! Cannot create output directory\n";
    }
    elsif (! -d $hCmdLineOption{'outdir'}) {
            die "ERROR! $hCmdLineOption{'outdir'} is not a directory\n";
    }
}

$sOutDir = File::Spec->canonpath($sOutDir);
$hCmdLineOption{'outrep'} = File::Spec->canonpath($hCmdLineOption{'outrep'});

mkdir ("$sOutDir/Quality");

$path =  $hCmdLineOption{'outrep'}."/fastqc_stats/".$hCmdLineOption{'pipeline'}."_fastqc/fastqc_stats.base.png.list";
if( -e $path ) {
    open($filehandle,"<$path") or die "Cannot open png file list";
    while (<$filehandle>) {
	chomp ($_);
	($f1, $path, $prefix) = File::Spec->splitpath($_);
	@arr = split(/\./,$prefix);
	$prefix = $arr[0];
	$prefix =~ s/_1_sequence//;
        @arr = split(/_/,$prefix);
	$prefix = '';
	for ($a = 0;$a < (scalar @arr-1);$a++){
	    if ($a!=0){$prefix.='_';}
	    $prefix.=$arr[$a];
	}
	$f1 = scalar @arr-1;
	$read = $arr[$f1];
	$final{$hCmdLineOption{'pipeline'}}{$prefix}{'base_png'}{$read} = $_;
    }
}
$path =  $hCmdLineOption{'outrep'}."/fastqc_stats/".$hCmdLineOption{'pipeline'}."_fastqc/fastqc_stats.length.png.list";
if( -e $path ) {
    open($filehandle,"<$path") or die "Cannot open seq length png file list";
    while (<$filehandle>) {
	chomp ($_);
	($f1, $path, $prefix) = File::Spec->splitpath($_);
	@arr = split(/\./,$prefix);
	$prefix = $arr[0];
	$prefix =~ s/_1_sequence//;
        @arr = split(/_/,$prefix);
	$prefix = '';
	for ($a = 0;$a < (scalar @arr-1);$a++){
	    if ($a!=0){$prefix.='_';}
	    $prefix.=$arr[$a];
	}
	$f1 = scalar @arr-1;
	$read = $arr[$f1];
	$final{$hCmdLineOption{'pipeline'}}{$prefix}{'length_png'}{$read} = $_;
    }
}

$path =  $hCmdLineOption{'outrep'}."/percent_mapped_stats/".$hCmdLineOption{'pipeline'}."_percent_mapped/percent_mapped_stats.txt.list";
if( -e $path ) {
    open($filehandle,"<$path") or die "Cannot open percent_mapped_stats.list file";
    while (<$filehandle>) {
	chomp($_);
	($f1, $path, $prefix) = File::Spec->splitpath($_);
	if (!exists $sample{$hCmdLineOption{'pipeline'}}) {
	    $sample{$hCmdLineOption{'pipeline'}}{'percent'} = [$_];
        }
	else {
	    push (@{$sample{$hCmdLineOption{'pipeline'}}{'percent'}} , $_);
	}
    }
    close $filehandle;
}


$path = $hCmdLineOption{'outrep'}."/align_tophat_stats/".$hCmdLineOption{'pipeline'}."_tophat_stats/align_tophat_stats.txt.list";
if( -e $path ) { 
    open($filehandle,"<$path") or die "Cannot open Tophat_stats.list file";
    $flag = 't';
    while (<$filehandle>) {
	chomp($_);
	($f1, $path, $prefix) = File::Spec->splitpath($_);
	if (exists $sample{$hCmdLineOption{'pipeline'}}) {
	    $sample{$hCmdLineOption{'pipeline'}}{'tophat'} = $_;
	}
    }
    close $filehandle;
}

$path =  $hCmdLineOption{'outrep'}."/align_bowtie_stats/".$hCmdLineOption{'pipeline'}."_bowtie_stats/align_bowtie_stats.txt.list";
if( -e $path ) {
    open($filehandle,"<$path") or die "Cannot open Bowtie_stats.list file";
    $flag = 'b';
    while (<$filehandle>) {
	chomp($_);
	($f1, $path, $prefix) = File::Spec->splitpath($_);
	if (exists $sample{$hCmdLineOption{'pipeline'}}) {
	    $sample{$hCmdLineOption{'pipeline'}}{'bowtie'} = $_;
	}
    }
    close $filehandle;
}	

$path =  $hCmdLineOption{'outrep'}."/rpkm_coverage_stats/".$hCmdLineOption{'pipeline'}."_rpkm_cvg/rpkm_coverage_stats.rpkm.stats.list";
if (-e $path) {
    open($filehandle,"<$path") or die "Cannot open Rpkm file list.";
    while (<$filehandle>){
	chomp($_);
	($f1, $path, $prefix) = File::Spec->splitpath($_);
	if (!exists $sample{$hCmdLineOption{'pipeline'}}) {
	    $sample{$hCmdLineOption{'pipeline'}}{'rpkm'} = [$_];
        }
	else {
	    push (@{$sample{$hCmdLineOption{'pipeline'}}{'rpkm'}} , $_);
	}
    }
    close $filehandle;
}	


foreach $key (keys %sample) {
    if (exists $sample{$key}{'percent'}) {
	foreach $a (@{$sample{$key}{'percent'}}) {
	    open($file,"<$a");
	    ($f1, $path, $prefix) = File::Spec->splitpath($a);
	    @arr = split(/\./,$prefix);
	    $prefix = $arr[0];
	    $prefix =~ s/_1//;
	    while (<$file>) {
		chomp($_);
		if (! exists $final{$key}{$prefix}{'percent'}){
		    $final{$key}{$prefix}{'percent'} = [$_];
		}
		else {
		    push (@{$final{$key}{$prefix}{'percent'}}, $_);
		}
	    }
	    close $file;
	}
	
    }
 if (exists $sample{$key}{'rpkm'}) {
	foreach $a (@{$sample{$key}{'rpkm'}}) {
	    open($file,"<$a");
	    ($f1, $path, $prefix) = File::Spec->splitpath($a);
	    @arr = split(/\./,$prefix);
	    $prefix = $arr[0];
	    $prefix =~ s/_1//;
	    while (<$file>) {
		chomp($_);
		if ($_ =~ /^\#/) {
		    if (! exists $final{$key}{$prefix}{'rpkm'}){
			$final{$key}{$prefix}{'rpkm'} = [$_];
		    }
		    else {
			push (@{$final{$key}{$prefix}{'rpkm'}}, $_);
		    }
		}
	    }
	    close $file;
	}
	
    }


    if (exists $sample{$key}{'bowtie'}) {
	open($file,"<$sample{$key}{'bowtie'}");
	while (<$file>) {
	    if($_=~/^\#Sample/){
		@arr = split(/\t/,$_);
		$prefix = $arr[1];
		$prefix =~ s/_1//;
		<$file>;
		$val = <$file>;
		chomp ($prefix);
		chomp ($val);
		if (exists $final{$key}{$prefix}) {
		    $final{$key}{$prefix}{'bowtie'} = $val;
		}
	    }
	}
	close $file;
    }
		
    if (exists $sample{$key}{'tophat'}) {
	open($file,"<$sample{$key}{'tophat'}");
	while (<$file>) {
	    if($_=~/^\#Sample/){
		@arr = split(/\t/,$_);
		$prefix = $arr[1];
		$prefix =~ s/_1//;
		<$file>;
		$val = <$file>;
		chomp($prefix);
		chomp($val);
		if (exists $final{$key}{$prefix}) {
		    $final{$key}{$prefix}{'tophat'} = $val;
		}
	    }
	}
	close $file;
    }    
}

open ($fout,">$sOutDir/$hCmdLineOption{pipeline}.txt") or die "Error Cannot open output file.";
open ($ffile,">$sOutDir/$hCmdLineOption{pipeline}.pngs.list") or die "Error Cannot open output file.";

foreach (keys %final) {
    print $fout "\#Pipeline id\t$_\n";
    if ($flag eq 'b') {print $fout "\#Sample_ID\tTotal.Reads\tTotal.Mapped.Reads\tPercent.Mapped.Reads\tPercent.Properly.Paired\tUniquely.Mapped.Reads\tPercent.Genic\tPercent.Intergenic\tTotal.Features\tFeatures.With.Coverage\tAvg.RPKM\n";}
    if ($flag eq 't') {print $fout "\#Sample_ID\tTotal.Reads\tTotal.Mapped.Reads\tPercent.Mapped.Reads\tPercent.Properly.Paired\tUniquely.Mapped.Reads\tPercent.Exonic\tPercent.Intronic\tPercent.Intergenic\tTotal.Features\tFeatures.With.Coverage\tAvg.RPKM\n";}
    foreach $f1 (keys %{$final{$_}}) {
	if (exists $final{$_}{$f1}{'bowtie'}){	    
	    print $fout "$f1\t$final{$_}{$f1}{'bowtie'}\t";
	}
	if (exists $final{$_}{$f1}{'tophat'}){
	    print $fout "$f1\t$final{$_}{$f1}{'tophat'}\t";
	}
	if (exists $final{$_}{$f1}{'percent'}){
	    $val = (split(/:/,$final{$_}{$f1}{'percent'}[1]))[1];
	    $val =~ s/\s//g;
	    print $fout "$val\t";
	    @arr = split(/\t/,$final{$_}{$f1}{'percent'}[3]);
	    foreach $a (@arr) {
		print $fout "$a\t";
	    }
	}
	if (exists $final{$_}{$f1}{'rpkm'}){
	    $val = (split(/[\:\(]/,$final{$_}{$f1}{'rpkm'}[0]))[1];
	    $val =~ s/\s//g;
	    print $fout "$val\t";
	    $val = (split(/[\:]/,$final{$_}{$f1}{'rpkm'}[2]))[1];
	    $val =~ s/\s//g;
	    print $fout "$val\t";
	    $val = (split(/[\:\(]/,$final{$_}{$f1}{'rpkm'}[4]))[1];
	    $val =~ s/\s//g;
	    print $fout "$val\t";
	    
	}
	print $fout "\n";
        if (exists $final{$_}{$f1}{'base_png'}){
	    if (-e $final{$_}{$f1}{'base_png'}{1}){
		print $ffile "$f1\tRead1\tBase_Qualtiy\t$sOutDir/Quality/$f1.R1.base_png\n";
		symlink  "$final{$_}{$f1}{'base_png'}{1}","$sOutDir/Quality/$f1.R1.base_png";
	    }
	    if (-e $final{$_}{$f1}{'base_png'}{2}){
		print $ffile "$f1\tRead2\tBase_Quality\t$sOutDir/Quality/$f1.R2.base_png\n";
		symlink  "$final{$_}{$f1}{'base_png'}{2}","$sOutDir/Quality/$f1.R2.base_png";
	    }
	    
	}
	if (exists $final{$_}{$f1}{'length_png'}){
	    if (-e $final{$_}{$f1}{'length_png'}{1}){
		print $ffile "$f1\tRead1\tRead_Length\t$sOutDir/Quality/$f1.R1.length_png\n";
		symlink  "$final{$_}{$f1}{'length_png'}{1}","$sOutDir/Quality/$f1.R1.length_png";
	    }
	    if (-e $final{$_}{$f1}{'base_png'}{2}){
		print $ffile "$f1\tRead2\tRead_Length\t$sOutDir/Quality/$f1.R2.length_png\n";
		symlink  "$final{$_}{$f1}{'length_png'}{2}","$sOutDir/Quality/$f1.R2.length_png";
	    }	
    	}
    }
} 
    

close $fout;
close $ffile;



######################
##Subroutine
#####################
sub check_parameters {
     my $phOptions = shift;    
     ## make sure input files provided
     if (! (defined $phOptions->{'outrep'}) ) {
                 pod2usage( -msg => $sHelpHeader, -exitval => 1);
     }

     if (! (defined $phOptions->{'pipeline'}) ) {
                pod2usage( -msg => $sHelpHeader, -exitval => 1);
     }
}



