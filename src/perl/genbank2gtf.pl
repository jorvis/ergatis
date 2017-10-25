#!/usr/bin/env perl -w
##############################################################################
### This program converts a genbank file to a GTF format file
##############################################################################

use strict;
use Cwd;
use Carp;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Spec;
use Bio::SeqIO;

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

my (@aFiles, @aInputFiles);
my ($sOutDir, $sFile);
my ($fpIN);
my %hUniq;

##############################################################################
### Main
##############################################################################

GetOptions( \%hCmdLineOption,
            'genbank|g=s', 'listfile|l=s', 'outdir|o=s', 'source|s=s',
            'fasta|f', 'verbose|v',
            'debug',
            'help',
            'man') or pod2usage(2);

pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};

if ($hCmdLineOption{'help'} ||
   (!((defined $hCmdLineOption{'genbank'}) || (defined $hCmdLineOption{'listfile'})))) {
    pod2usage( -msg => $sHelpHeader, -exitval => 1);
}


@aInputFiles = ();
if (defined $hCmdLineOption{'genbank'}) {
	if (-d $hCmdLineOption{'genbank'}) {
    	@aInputFiles = glob($hCmdLineOption{'genbank'}."/*.gb*");
	}
	else {
    	@aInputFiles = glob($hCmdLineOption{'genbank'});
	}
}

if (defined $hCmdLineOption{'listfile'}) {
	open ($fpIN, "<$hCmdLineOption{'listfile'}") || 
		die "\tCannot open $hCmdLineOption{'listfile'} for reading .....\n";
	
	while (<$fpIN>) {
    	$_ =~ s/\s+$//;
    	next if ($_ =~ /^\#/);
    	next if ($_ =~ /^$/);
    	push @aInputFiles, $_;
    }
    
    close($fpIN);
}

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

# process the input files
foreach $sFile (@aInputFiles) {

    GBK2GTF(\%hCmdLineOption, $sFile, $sOutDir);
    
}

exit;


##############################################################################
### Subroutines
##############################################################################


# GBK2GTF()
#
# Purpose
#   Converts genbank file to GTF format file
#
# Required Parameters
#   phCmdLineOption = pointer to hash of command line options
#   sInFile         = /path/to/bed_file
#   sOutDir         = directory to store results
#
# Optional Parameters
#   None
#
# Returns
#   None
#
# Side Effects
#   None
#
# Assumptions
#	sInFile is a genbank file to create a Bio::SeqIO object
#
# Notes
#
sub GBK2GTF {
    my $phCmdLineOption = shift;
    my $sInFile         = shift;
    my $sOutDir         = shift;
    
    my $sSubName = (caller(0))[3];

    # check for required parameters
    if (! ((defined $phCmdLineOption) &&
           (defined $sInFile) &&
           (defined $sOutDir))) {
        croak "$sSubName - required parameters missing\n";
    }

    # Local variables
	    
    my $bDebug   = (defined $phCmdLineOption->{'debug'}) ? TRUE : FALSE;
    my $bVerbose = (defined $phCmdLineOption->{'verbose'}) ? TRUE : FALSE;
    my $sSource  = (defined $phCmdLineOption->{'source'}) ? 
    					$phCmdLineOption->{'source'} : "GenBank";
    
    my ($oFile, $oSeq, $oFeature, $oLocation, $oFasta);
    my ($sFastaFile, $sGtfFile);
    my (@aTmpArray, @aAttributes);
    my ($sAccID, $nStart, $nEnd, $sStrand);
    my ($sGID, $sTID, $sCID, $sEID, $sID, $sGeneName, $sTranscriptName);
    my ($nI, $nM, $nC, $nX);
    my ($fpGTF);
    my ($bFlag);
    
    ($bDebug || $bVerbose) ? print STDERR "Parsing $sInFile .....\n" : ();
    
    die "\tCannot open $sInFile .....\n" if (! -e $sInFile);
    
    # Bio::SeqIO object creation
    $oFile = Bio::SeqIO->new(-file => $sInFile,
                             -format => 'genbank');
    
    while ( $oSeq = $oFile->next_seq() ) {
    	$sAccID = $oSeq->id;
    	
    	if (defined $hCmdLineOption{'fasta'}) {
    		$sFastaFile = $sOutDir."/".$sAccID.".gb.fa";
    		$oFasta = Bio::SeqIO->new(-file => ">$sFastaFile",
                                	  -format => 'fasta');
            
            $oFasta->write_seq($oSeq);
    	}
    	
    	$sGtfFile = $sOutDir."/".$sAccID.".gb.gtf";
    	open ($fpGTF, ">$sGtfFile") || 
			die "\tCannot open $sGtfFile for writing .....\n";
		
    	$nI = $bFlag = 0;
    	foreach $oFeature ($oSeq->get_SeqFeatures) {
    		next if ( $oFeature->primary_tag eq 'source' );
    		if ( $oFeature->primary_tag eq 'gene' ) {
    			$bFlag = $nM = $nC = 0;
    			if ($oFeature->has_tag('gene')) {
    				@aTmpArray = $oFeature->get_tag_values("gene");
    				$sGeneName = $aTmpArray[0];
    				$sGeneName =~ s/\s+.*//;
    			}
    			else {
    				$nI++;
    				$sGeneName = sprintf("orf%8d", $nI);
    			}
    			$sGID = $sGeneName;
    			
    			@aAttributes = ("ID", $sGID, "gene_name", $sGeneName);
    			
    			Print_Annotation(\%hCmdLineOption, $sAccID, $sSource, $oFeature->primary_tag, $oFeature->location->start,
    							 $oFeature->location->end, $oFeature->location->strand, \@aAttributes, $fpGTF);
    		}
    		elsif ( $oFeature->primary_tag eq 'mRNA' ) {
    			$nM++;
    			if ($oFeature->has_tag('gene')) {
    				@aTmpArray = $oFeature->get_tag_values("gene");
    				$sGeneName = $aTmpArray[0];
    				$sGeneName =~ s/\s+.*//;
    			}
    			else {
    				$sGeneName = $sGID;
    			}
    			$sTID = $sGeneName."-T".$nM;
    			
    			if ($oFeature->has_tag('transcript_id')) {
    				@aTmpArray = $oFeature->get_tag_values("transcript_id");
    				$sTranscriptName = $aTmpArray[0];
    				$sTranscriptName =~ s/\s+.*//;
    			}
    			else {
    				$sTranscriptName = $sTID;
    			}
    			
    			@aAttributes = ("ID", $sTID, "gene_name", $sGeneName, "transcript_name", $sTranscriptName);
    			
    			if ( $oFeature->location->isa('Bio::Location::SplitLocationI') ) {
    				$bFlag = 1;
    				
    				$nStart = $nEnd = 0;
    				for $oLocation ( $oFeature->location->sub_Location ) {
    					$nStart = $oLocation->start if ($nStart == 0);
    					$nEnd = $oLocation->end;
    				}
    				
    				Print_Annotation(\%hCmdLineOption, $sAccID, $sSource, $oFeature->primary_tag, $nStart,
    							 	 $nEnd, $oFeature->location->strand, \@aAttributes, $fpGTF);
    				
    				$nX = 0;
         			for $oLocation ( $oFeature->location->sub_Location ) {
         				$nX++;
         				$sEID = $sTID."-E".$nX;
         				@aAttributes = ("ID", $sEID, "gene_name", $sGeneName, "transcript_name", $sTranscriptName);
         				
         				Print_Annotation(\%hCmdLineOption, $sAccID, $sSource, "exon", $oLocation->start,
    							 $oLocation->end, $oFeature->location->strand, \@aAttributes, $fpGTF);
					}
				}
				else {
					Print_Annotation(\%hCmdLineOption, $sAccID, $sSource, $oFeature->primary_tag, $oFeature->location->start,
    							 	 $oFeature->location->end, $oFeature->location->strand, \@aAttributes, $fpGTF);
				}
			}
    		elsif ( ($oFeature->primary_tag eq 'CDS') && ($bFlag == 0) ) {
    			$nC++;
    			if ($oFeature->has_tag('gene')) {
    				@aTmpArray = $oFeature->get_tag_values("gene");
    				$sGeneName = $aTmpArray[0];
    				$sGeneName =~ s/\s+.*//;
    			}
    			else {
    				$sGeneName = $sGID;
    			}
    			$sTID = $sGeneName."-T".$nM;
    			$sCID = $sGeneName."-P".$nC;
    			
    			if ($oFeature->has_tag('transcript_id')) {
    				@aTmpArray = $oFeature->get_tag_values("transcript_id");
    				$sTranscriptName = $aTmpArray[0];
    				$sTranscriptName =~ s/\s+.*//;
    			}
    			else {
    				$sTranscriptName = $sTID;
    			}
    			
    			@aAttributes = ("ID", $sCID, "gene_name", $sGeneName, "transcript_name", $sTranscriptName);
    			
    			if ( $oFeature->location->isa('Bio::Location::SplitLocationI') ) {
    				$nStart = $nEnd = 0;
    				for $oLocation ( $oFeature->location->sub_Location ) {
    					$nStart = $oLocation->start if ($nStart == 0);
    					$nEnd = $oLocation->end;
    				}
    				
    				Print_Annotation(\%hCmdLineOption, $sAccID, $sSource, $oFeature->primary_tag, $nStart,
    							 	 $nEnd, $oFeature->location->strand, \@aAttributes, $fpGTF);
    				
    				$nX = 0;
         			for $oLocation ( $oFeature->location->sub_Location ) {
         				$nX++;
         				$sEID = $sTID."-E".$nX;
         				@aAttributes = ("ID", $sEID, "gene_name", $sGeneName, "transcript_name", $sTranscriptName);
         				
         				Print_Annotation(\%hCmdLineOption, $sAccID, $sSource, "exon", $oLocation->start,
    							 $oLocation->end, $oFeature->location->strand, \@aAttributes, $fpGTF);
					}
				}
				else {
					Print_Annotation(\%hCmdLineOption, $sAccID, $sSource, $oFeature->primary_tag, $oFeature->location->start,
    							 	 $oFeature->location->end, $oFeature->location->strand, \@aAttributes, $fpGTF);
				}
			}
			else {
				if ($oFeature->has_tag('gene')) {
    				@aTmpArray = $oFeature->get_tag_values("gene");
    				$sGeneName = $aTmpArray[0];
    				$sGeneName =~ s/\s+.*//;
    			}
    			else {
    				$nI++;
    				$sGeneName = sprintf("orf%8d", $nI);
    			}
    			$sID = $sGeneName."-X";
    			
    			@aAttributes = ("ID", $sID, "gene_name", $sGeneName);
    			
    			Print_Annotation(\%hCmdLineOption, $sAccID, $sSource, $oFeature->primary_tag, $oFeature->location->start,
    							 $oFeature->location->end, $oFeature->location->strand, \@aAttributes, $fpGTF);
			}
    	}
    	
    	close($fpGTF);
    }
    
    ($bDebug || $bVerbose) ? print STDERR "Parsing $sInFile ..... Done\n" : ();
    
    ($bDebug || $bVerbose) ? print STDERR "\n" : ();
    
	($bDebug) ? print STDERR "Leaving $sSubName\n" : ();

    return;
}

# Print_Annotation()
#
# Purpose
#   prints annotation to output file or putput screen
#
# Required Parameters
#   phCmdLineOption = pointer to hash containing command line options
#   sAccID			= Accession id
#	sSource			= source
#	sFeature		= feature tag
#	nStart			= feature start position
#	nEnd			= feature end position
#   sStrand			= feature strand
#   paAttribute     = pointer to array of attributes
#	fpOUT			= pointer to output file
#
# Optional Parameters
#   none
#
# Returns
#   none
#
# Side Effects
#   writes to external file or screen
#
# Assumptions
#
# Notes
#
sub Print_Annotation {
	my $phCmdLineOption	= shift;
	my $sAccID			= shift;
	my $sSource			= shift;
	my $sFeature		= shift;
	my $nStart			= shift;
	my $nEnd			= shift;
	my $sStrand			= shift;
	my $paAttribute		= shift;
	my $fpOUT			= shift;

    my $sSubName = (caller(0))[3];

    if (! ((defined $sAccID) &&
           (defined $sSource) &&
           (defined $sFeature) &&
           (defined $nStart) &&
           (defined $nEnd) &&
           (defined $sStrand) &&
           (defined $paAttribute) &&
           (defined $phCmdLineOption))) {
        croak "$sSubName - required parameters missing\n";
    }

    # Local variables
    my $bDebug   = (defined $phCmdLineOption->{'debug'}) ? TRUE : FALSE;
    my $bVerbose = (defined $phCmdLineOption->{'verbose'}) ? TRUE : FALSE;
    
    my ($nI);
    
    # Start

    ($bDebug) ? print STDERR "In $sSubName\n" : ();

    if (defined $fpOUT) {
    	print $fpOUT $sAccID."\t".
    				 $sSource."\t".
    				 $sFeature."\t".
    				 $nStart."\t".
    				 $nEnd."\t".
    				 ".\t";
    	
    	if ($sStrand == 1) { print $fpOUT "+\t"; }
    	elsif ($sStrand == -1) { print $fpOUT "-\t"; }
    	else { print $fpOUT ".\t"; }
    	
    	print $fpOUT ".\t";
    	
    	for ($nI = 0; $nI < @{$paAttribute}; $nI=$nI+2) {
    		print $fpOUT $$paAttribute[$nI]." "."\"".$$paAttribute[$nI+1]."\"; ";
    	}
    	
    	print $fpOUT "\n";
    }
    else {
    	print STDERR $sAccID."\t".
    				 $sSource."\t".
    				 $sFeature."\t".
    				 $nStart."\t".
    				 $nEnd."\t".
    				 ".\t";
    	
    	if ($sStrand == 1) { print STDERR "+\t"; }
    	elsif ($sStrand == -1) { print STDERR "-\t"; }
    	else { print STDERR ".\t"; }
    	
    	print STDERR ".\t";
    	
    	for ($nI = 0; $nI < @{$paAttribute}; $nI=$nI+2) {
    		print STDERR $$paAttribute[$nI]." "."\"".$$paAttribute[$nI+1]."\"; ";
    	}
    	
    	print STDERR "\n";
    }

    ($bDebug) ? print STDERR "Leaving $sSubName\n" : ();
    
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

    ($bDebug) ? print STDERR "In $sSubName\n" : ();

    ($_, $sDir, $sFile) = File::Spec->splitpath($sFileName);
    if ($sOutDir ne File::Spec->curdir()) {
        $sDir = '';
    }
	$sFile =~ m/^(\S+)($sExtension)$/;
	if (defined $1) {
	    $sOutFile = $sOutDir.'/'.$sDir.'/'.$1;
	}
	else {
	    $sOutFile = $sOutDir.'/'.$sDir.'/'.$sFile;
	}
    $sOutFile = File::Spec->canonpath($sOutFile);

    ($bDebug) ? print STDERR "Leaving $sSubName\n" : ();
    
    return $sOutFile
}


##############################################################################
### POD Documentation
##############################################################################

__END__

=head1 NAME

genbank2gtf.pl - program to generate a GTF format file from a GenBank file.

=head1 SYNOPSIS

    genbank2gtf.pl --g <genbank> --l <listfile> [--o <outdir>] [--s <source>] [--f] [--v]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --g <genbank>   = /path/to/genbank_file. If a directory, it will
                      be searched for files with .gb extension.
                      Wildcards can be used if enclosed in quotes,
                      e.g. --i "NLGN*.bed".

    --l <listfile>  = List file of genbank files to be converted to GTF files.

    --s <source>    = Source of genbank file (e.g. NCBI, EMBL, etc). [Genbank]

    --o <outdir>    = /path/to/output directory. Optional.

    --f             = generate FastA sequence file from Genbank file. Optional.

    --v             = generate runtime messages. Optional

=head1 DESCRIPTION

The program parses a genbank file and generates a GTF file from the genbank annotation section.

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

Copyright (c) 2013 Amol Carl Shetty (<ashetty@som.umaryland.edu>). All rights
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
