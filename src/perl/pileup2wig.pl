#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

##############################################################################
### This program converts a pileup file into wig format file
##############################################################################

use strict;
use warnings;
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

use constant VERSION => '1.0.0';
use constant PROGRAM => eval { ($0 =~ m/(\w+\.pl)$/) ? $1 : $0 };

##############################################################################
### Globals
##############################################################################

my %hCmdLineOption = ();
my $sHelpHeader = "\nThis is ".PROGRAM." version ".VERSION."\n";

my ($sPileupFile, $sWigFile, $sOutDir);
my (@aTmpArray);
my ($sChr, $nPos, $nCvg, $sPrevChr, $nPrevPos);
my ($fpPUP, $fpWIG);
my ($bDebug, $bVerbose);

##############################################################################
### Main
##############################################################################

GetOptions( \%hCmdLineOption,
            'pileup|p=s', 'outdir|o=s',
            'verbose|v',
            'debug',
            'help',
            'man') or pod2usage(2);

pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};

pod2usage( -msg => $sHelpHeader, -exitval => 1) if $hCmdLineOption{'help'};

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

$bDebug   = (defined $hCmdLineOption{'debug'}) ? TRUE : FALSE;
$bVerbose = (defined $hCmdLineOption{'verbose'}) ? TRUE : FALSE;

if (! (defined $hCmdLineOption{'pileup'}) ) {
    pod2usage( -msg => $sHelpHeader, -exitval => 1);
}

$sPileupFile = $hCmdLineOption{'pileup'};

($_, $_, $sWigFile) = File::Spec->splitpath($sPileupFile);
$sWigFile =~ s/.txt$//;
$sWigFile = $sOutDir."/".$sWigFile.".wig";

($bDebug || $bVerbose) ?
	print STDERR "Converting $sPileupFile to $sWigFile .....\n" : ();

open ($fpPUP, "<$sPileupFile") or die "\tERROR : Cannot open $sPileupFile for reading .....\n";

open ($fpWIG, ">$sWigFile") or die "\tERROR : Cannot open $sWigFile for writing .....\n";

$sPrevChr = "";
$nPrevPos = 0;

while (<$fpPUP>) {
	$_ =~ s/\s+$//;
	
	@aTmpArray = split(/\t/, $_);
	$sChr = $aTmpArray[0];
	$nPos = $aTmpArray[1];
	$nCvg = $aTmpArray[$#aTmpArray - 2];
	
	if ($sChr ne $sPrevChr) { 
		print $fpWIG "fixedStep chrom=$sChr start=$nPos step=1\n$nCvg\n";
	}
	else {
		if ($nPos == ($nPrevPos + 1)) {
			print $fpWIG "$nCvg\n";
		}
		else {
			print $fpWIG "fixedStep chrom=$sChr start=$nPos step=1\n$nCvg\n";
		}
	}
	
	$sPrevChr=$sChr;
	$nPrevPos=$nPos;
	
	($bDebug || $bVerbose) ?
		print STDERR "\r\tProcessing $sChr:$nPos ....." : ();
}

close($fpWIG);
close($fpPUP);

($bDebug || $bVerbose) ?
	print STDERR "\n\nDone! Stop looking at the log file!\n\n" : ();

##############################################################################
### POD Documentation
##############################################################################

__END__

=head1 NAME

pileup2wig.pl - program to convert a pileup file to a wig format file.

=head1 SYNOPSIS

    pileup2wig.pl --p <pileup file> [--o <outdir>] [--v]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --p <pileup file>      = /path/to/pileup_file.

    --o <outdir>           = /path/to/output directory. Optional

    --v                    = generate runtime messages. Optional

=head1 DESCRIPTION

This program converts a pileup file into wig format file

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

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

Copyright (c) 2011 Amol Carl Shetty (<ashetty@som.umaryland.edu>). All rights
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



