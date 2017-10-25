#!/usr/bin/env perl -w

eval 'exec /usr/bin/env perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell
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

use constant SAMTOOLS_BIN_DIR => '/usr/local/bin';

##############################################################################
### Globals
##############################################################################

my %hCmdLineOption = ();
my $sHelpHeader = "\nThis is ".PROGRAM." version ".VERSION."\n";

my (@inBamFile);
my ( $BamInfoIN,$BamListIN ); 
my ($sBamInfoFile, $bamFileBase, $bamFileNum, $sPrefix, $splitInfo, $sBamFileList);
my ($sOutDir, $outBamFile);
my ($sCmd, $nNumCols);
my ($bDebug, $bVerbose);

##############################################################################
### Main
##############################################################################

GetOptions( \%hCmdLineOption,
            'inBamFileList|l=s', 'sBamInfoFile|i=s', 
			'outdir|o=s', 'samtools_bin_dir|s=s', 
            'verbose|v',
            'debug',
            'help',
            'man') or pod2usage(2);

if ($hCmdLineOption{'help'} || 
    (! defined $hCmdLineOption{'inBamFileList'})|| (! defined $hCmdLineOption{'sBamInfoFile'}) ){
    pod2usage( -msg => $sHelpHeader, -exitval => 1);
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

#Get split information
$sBamInfoFile = $hCmdLineOption{'sBamInfoFile'};
open($BamInfoIN, "<$sBamInfoFile") or die "\t Cannot open $sBamInfoFile for reading\n";
$splitInfo = <$BamInfoIN>;
#print "$splitInfo\n";
chomp($splitInfo);
( $bamFileBase, $bamFileNum ) = split ("\t", $splitInfo);
$bamFileBase =~ s/.*\///;
close $BamInfoIN;

#Get bamfile to be merged
$sBamFileList = $hCmdLineOption{'inBamFileList'};

open($BamListIN, "<$sBamFileList") or die "\t Cannot open $sBamFileList for reading\n";
while(<$BamListIN>){
	chomp;

	if ( $_ =~  /$bamFileBase/){
		push @inBamFile, $_;
	}
}
close $BamListIN;

#Check split bam file number
if ($bamFileNum != ($#inBamFile + 1)){
	print STDERR "\nSplit File number doesn't match. Exit .....\n" ;

}


#Get merged bam file name as regular alignment bamfile  
$sPrefix = $bamFileBase;
$sPrefix =~ s/.1_1_dequence.*//;
$sPrefix =~ s/_R1_.*//;
$outBamFile = $sOutDir."/".$sPrefix.".accepted_hits.bam";

# start merging process
($bDebug || $bVerbose) ? 
	print STDERR "\nMerging bam files .....\n" : ();

if (! (defined $hCmdLineOption{'samtools_bin_dir'}) ) {
	$hCmdLineOption{'samtools_bin_dir'} = SAMTOOLS_BIN_DIR;
}

$sCmd = $hCmdLineOption{'samtools_bin_dir'}."/samtools merge".
		" -f".
		" $outBamFile";

for (my $i = 0; $i <= $#inBamFile; $i++) {
	$sCmd .= " $inBamFile[$i]";
}

system($sCmd);
			
($bDebug || $bVerbose) ? 
	print STDERR "Generating file: $sPrefix.accepted_hits.bam  ..... done\n" : ();


##############################################################################
### POD Documentation
##############################################################################

__END__

=head1 NAME

	bam_merge.pl - program to merge bam file from alignment of splitted fastq file 

=head1 SYNOPSIS

    bam_merge.pl      --l [bam_file_list] --i <fastq_file_split_info> 
                      [--o <output_dir>] [--b <bedtools_bin_dir>][--s <samtools_bin_dir>] 
                      [--v]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

   
    --l <bam_file_list>                             = path to splitted BAM file list

    --i < fastq file split infomation>              = path to fastq file split information

    --o <output_dir>                                = output directory. Optional

    --s <samtools_bin_dir>                          = samtools binary directory. Optional

    --v                                             = generate runtime messages. Optional

=head1 DESCRIPTION

=head1 DIAGNOSTICS


=head1 BUGS AND LIMITATIONS

There are no known bugs in this module. Please report problems to Amol Shetty
(ashetty@som.umaryland.edu). Patches are welcome.

=head1 AUTHOR

 Weizhong Chang
 Bioinformatics Software Engineer
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
