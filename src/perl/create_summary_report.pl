#!/usr/local/bin/perl -w

#################################################################################
#										#
# Name	      : create_summary_report.pl					#
# Version     : 1.0								#
# Project     : CloVR Microbe Pipeline						#
# Description : Script to generate a summary report for the Microbe pipeline	#
# Author      : Sonia Agrawal							#
# Date        : February 11, 2013						#
#										#
#################################################################################

use strict;
# Extended processing of command line options
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
# Prints usage message from embedded pod documentation
use Pod::Usage;
# Processes huge XML(BSML) documents in tree mode
use XML::Twig;
use Math::Round;

###########
# GLOBALS #
###########
my %cmdLineArgs = ();
# Log file handle;
my $logfh;
my ($ERROR, $WARN, $DEBUG) = (1, 2, 3);
my ($file, $contig, $twig);
my @files = ();
my %contigs = ();
my @lengths = ();
my @sorted_len = ();
my $median = 0;
my $tcds = 0;
my $trrna = 0;
my $ttrna = 0;
my $tlength = 0;

################
# MAIN PROGRAM #
################
GetOptions(\%cmdLineArgs,
	   'bsml_file|f=s',
	   'bsml_list|b=s',
	   'bsml_dir|d=s',
	   'bsml_ext|e=s',
	   'output_file|o=s',
	   'log|l=s', 
	   'help|h'
	  ) or pod2usage();

pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ) if ($cmdLineArgs{'help'});

&checkCmdLineArgs();

open(SR, "> $cmdLineArgs{'output_file'}") or printLogMsg($ERROR, "ERROR! : Could not open file $cmdLineArgs{'output_file'} for writing. Reason: $!");
foreach $file (@files) {
	$contigs{$file}{'CDS'} = 0;
	$contigs{$file}{'tRNA'} = 0;
	$contigs{$file}{'rRNA'} = 0;
	$twig = XML::Twig->new();
	$twig->setTwigRoots({'Feature-group' => sub { processFeature($file, @_); }, 'Sequence' => sub { processSeq($file, @_); }});
	$twig->parsefile($file);
	$twig->purge();
	
}
$twig = XML::Twig->new();
$twig->setTwigRoots({'Organism' => sub { processGenome($files[$#files], @_);}});
$twig->parsefile($files[$#files]);
$twig->purge();

print SR "************************* SUMMARY REPORT FOR GENOME $contigs{$files[$#files]}{'organism'} *************************\n\n"; 
print SR "Contig Id\tSize(in Kbp) > 5Kbp\tNumber of CDS\n";
print SR "------------------------------------------------------\n";
foreach $contig (keys %contigs) {
	if(($contigs{$contig}{'length'}/1000) > 5) {
		print SR $contigs{$contig}{'id'}."\t".nearest(0.01,($contigs{$contig}{'length'}/1000))."\t".$contigs{$contig}{'CDS'}."\n";
	}
	$tcds += $contigs{$contig}{'CDS'};
	$trrna += $contigs{$contig}{'rRNA'};
	$ttrna += $contigs{$contig}{'tRNA'};
	$tlength += $contigs{$contig}{'length'};
	push(@lengths,$contigs{$contig}{'length'});
}
@sorted_len = sort {$a <=> $b} @lengths;
#Odd
if(@sorted_len%2) {
	$median = $sorted_len[int(@sorted_len/2)];		
} else {
	$median = ($sorted_len[(@sorted_len/2)] + $sorted_len[(@sorted_len/2) - 1])/2;
}
print SR "\nTotal number of contigs/scaffolds : ".@files."\n";
print SR "Mean contig size : ".nearest(0.01,(($tlength/@files)/1000))." Kbp\n";
print SR "Median contig size : ".nearest(0.01,($median/1000))." Kbp\n";
print SR "Total genome size : ".nearest(0.01,$tlength/1000000)." Mbp\n";
print SR "Total number of CDS : $tcds\n";
print SR "Total number of tRNAs : $ttrna\n";
print SR "Total number of rRNAs : $trrna\n";

close(SR);
###############
# SUBROUTINES #
###############

# Description   : Used to extract genome species and genus from the BSML file. 
# Parameters    : NA
# Returns       : NA
# Modifications :
sub processGenome {
	my ($filename, $twig, $elt) = @_;
	$contigs{$filename}{'organism'} = $elt->att('genus')." ".$elt->att('species');
}

####################################################################################################################################################

# Description   : Used to process an assembly sequence, extract its length and id from the BSML file. 
# Parameters    : NA
# Returns       : NA
# Modifications :
sub processSeq {
	my ($filename, $twig, $elt) = @_;
	return if $elt->att('class') ne "assembly";
	$contigs{$filename}{'id'} = $elt->att('id');
	$contigs{$filename}{'length'} = $elt->att('length');
}

####################################################################################################################################################

# Description   : Used to find number of CDS, rRNA and tRNA in a contig
# Parameters    : NA
# Returns       : NA
# Modifications :
sub processFeature {
	my ($filename, $twig, $elt) = @_;
	my ($feat_member, $feat_type);
	foreach $feat_member ($elt->children('Feature-group-member')) {
		$feat_type = $feat_member->att('feature-type');
		if($feat_type eq "CDS") {
			$contigs{$filename}{'CDS'} += 1;
		}
		if($feat_type eq "rRNA") {
			$contigs{$filename}{'rRNA'}++;
		}
		if($feat_type eq "tRNA") {
			$contigs{$filename}{'tRNA'}++;
		}
	}
}

####################################################################################################################################################

# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing. 
# Parameters    : NA
# Returns       : NA
# Modifications :

sub checkCmdLineArgs {
	my ($option, $line);
	if(exists($cmdLineArgs{'log'})) {
		open($logfh, "> $cmdLineArgs{'log'}") or die "Could not open $cmdLineArgs{'log'} file for writing.Reason : $!\n"
	}
	my @required = qw(output_file);
        foreach $option(@required) {
                if(!defined($cmdLineArgs{$option})) {
                        printLogMsg($ERROR,"ERROR! : Required option $option not passed");
                }   
        }
	if($cmdLineArgs{'bsml_list'}) {
		open(BL, "< $cmdLineArgs{'bsml_list'}") or printLogMsg($ERROR, "ERROR! : Could not open file $cmdLineArgs{'bsml_list'} for reading. Reason: $!");
		while($line = <BL>) {
			chomp($line);
			push(@files, $line);
		}
	} elsif($cmdLineArgs{'bsml_file'}) {
		push(@files, $cmdLineArgs{'bsml_file'});
	} elsif($cmdLineArgs{'bsml_dir'}) {
		if($cmdLineArgs{'bsml_ext'}) {
			opendir(BD, "$cmdLineArgs{'bsml_dir'}") or printLogMsg($ERROR, "ERROR! : Could not open directory $cmdLineArgs{'bsml_dir'} for reading. Reason: $!");
			@files = grep { -f "$cmdLineArgs{'bsml_dir'}/$_" && /\.$cmdLineArgs{'bsml_ext'}$/ } readdir BD;	
			closedir(BD);
			@files = map { $cmdLineArgs{'bsml_dir'}."/".$_ } @files;
		} else {
			printLogMsg($ERROR,"ERROR! : Required option bsml_ext not passed when passing a bsml_dir");
		}
	} else {
		printLogMsg($ERROR,"ERROR! : One of the required options bsml_file or bsml_list or bsml_dir not passed");
	}
}

####################################################################################################################################################

# Description   : Used to handle logging of messages(errors and warnings) during the execution of the script
# Parameters    : level = can be ERROR, WARNING or INFO
#		  msg   = msg to be printed in the log file or to STDERR
# Returns       : NA
# Modifications : 

sub printLogMsg {
	my ($level, $msg) = @_;
	if( $level <= $DEBUG ) {
		print STDERR "$msg\n";
		print $logfh "$msg\n" if(defined($logfh));
		die "" if($level == $ERROR);
	}	
}

__END__

#####################
# POD DOCUMENTATION #
#####################

=head1 NAME

create_summary_report.pl - Script to generate a summary report for the Microbe pipeline

=head1 SYNOPSIS

USAGE : perl create_summary_report -bsml_file <input bsml file> -output_file <output file> [ -bsml_list <input bsml list> -bsml_dir <input bsml dir> -bsml_ext <bsml extension> -log <log file> -help ]

	parameters in [] are optional

=head1 OPTIONS

	-bsml_file	=	Input BSML file containing the annotation pipeline summary
	
	-output_file	=	Output summary report file

	-bsml_list	=	List of BSML files containing annotation pipeline summary for each molecule

	-bsml_dir	=	Directory of BSML files containing annotation pipeline summary for each molecule

	-bsml_ext	=	File extension to be used in case bsml_dir is provided as an argument. Default: bsml

	-log		=	Log file

=head1 DESCRIPTION

The script is used to generate a summary report for the clovr microbial pipeline.

=head1 INPUT

Pipeline summary BSML file

=head1 OUTPUT

Summary report file conating following information:
1. Contig Ids	Size(Kbp)	Number of CDS in the contig - This list only includes contigs greater than 5Kbp
2. Total number of contigs/scaffolds in the genome.
3. Mean contig size in Kbp
4. Median contig size in Kbp
5. Total genome size in Mbp
6. Total numbmer of CDS predicted in the genome
7. Total number of tRNAs predicted in the genome
8. Total number of rRNAs predicted in the genome 

=head1 AUTHOR

	Sonia Agrawal
	Bioinformatics Software Engineer II
	Institute for Genome Sciences
	School of Medicine
	University of Maryland
	Baltimore, MD - 21201
	sagrawal@som.umaryland.edu

==cut
