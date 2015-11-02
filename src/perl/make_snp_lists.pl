#!/usr/local/bin/perl -w

#####################################################################################
#
# Name	      :	make_snp_lists.pl
# Version     :	1.0
# Project     :	SNP Verification Pipeline
# Description : 
# Author      : Sonia Agrawal
# Date        :	November 14, 2011
#
#####################################################################################

use strict;
# Extended processing of command line options
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
# Prints usage message from embedded pod documentation
use Pod::Usage;
# Include packages here 
use File::Basename;
use Bio::SeqIO;

#############
# CONSTANTS #
#############


###########
# GLOBALS #
###########
my %cmdLineArgs = ();
my ($ref_name, $coords_file);
# Log file handle;
my $logfh;
my ($ERROR, $WARN, $DEBUG) = (1, 2, 3);

################
# MAIN PROGRAM #
################
GetOptions(\%cmdLineArgs,
	   'ref_genbank|r=s',
	   'snp_positions|s=s',
	   'blast_list|b=s',
	   'output_dir|o=s',
	   'log|l=s', 
	   'help|h'
	  ) or pod2usage();

pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ) if ($cmdLineArgs{'help'});

&checkCmdLineArgs();

if (-e $cmdLineArgs{'ref_genbank'}) {
	my ($gbk_base,$gbk_dir,$gbk_ext) = fileparse($cmdLineArgs{'ref_genbank'},qr/\.[^.]*/);
	$ref_name = $gbk_base;
} else {
	printLogMsg($ERROR, "File $cmdLineArgs{'ref_genbank'} does not exist.");
}
$coords_file = &parse_ref_genbank($cmdLineArgs{'ref_genbank'}, $cmdLineArgs{'output_dir'});
open(FR, "< $cmdLineArgs{'blast_list'}") or printLogMsg($ERROR, "Could not open file $cmdLineArgs{'blast_list'} for reading.\nReason : $!");
my $blast_sublist = $cmdLineArgs{'output_dir'}."/blast_sublist_".$ref_name.".txt";
open(FW, "> $blast_sublist") or printLogMsg($ERROR, "Could not open file $blast_sublist for writing.\nReason : $!");
while(<FR>) {
	my $line = $_;
	chomp($line);
	next if ($line =~ /^\s*$/);
	next if ($line =~ /^#/);
	if (-e $line) {
		my ($file_base,$file_dir,$file_ext) = fileparse($line,qr/\.[^.]*/);
		$file_base =~ /_refmol_(.+)/;
		my $gbk_id = $1;
		if ($gbk_id eq $ref_name) {
			print FW $line."\n";
		}
	}
}
close(FW);
close(FR);

my $snp_verify_list = $cmdLineArgs{'output_dir'}."/snp_verify_ip_".$ref_name.".order";
open(FOUT, "> $snp_verify_list") or printLogMsg($ERROR, "Could not open file $snp_verify_list for writing.\nReason : $!");
if(-e $blast_sublist) {
	print FOUT $coords_file."\n";
	print FOUT $cmdLineArgs{'snp_positions'}."\n";
	print FOUT $blast_sublist."\n";
	close(FOUT);
} else {
	printLogMsg($ERROR, "File $blast_sublist was not created.");
}

###############
# SUBROUTINES #
###############

# Description   : 
# Parameters    : Reference genome's genBank file path
#               : Output directory to store reference genome annotation file
# Returns       : NA
# Modifications :

sub parse_ref_genbank {
        my ($genbank, $result_dir) = @_; 
        my ($ref_base,$ref_dir,$ref_ext) = fileparse($genbank,qr/\.[^.]*/);
        my $annot_file = "$result_dir/".$ref_base.".coords";
        open(FTAB, "> $annot_file") or printLogMsg($ERROR, "Could not open $annot_file for writing.\nReason : $!");
        my $seqio_obj = Bio::SeqIO->new(-file => "$genbank", -format => "GenBank" );
        while (my $seqobj = $seqio_obj->next_seq()) {
                for my $feat_object ($seqobj->get_SeqFeatures) {
                        if ($feat_object->primary_tag eq "CDS") {
                                for my $locus ($feat_object->get_tag_values("locus_tag")) {
#                                        $locus =~ s/\_//;
                                        print FTAB "$locus\t";
                                }
                                if ($feat_object->location->strand == -1) {
                                        print FTAB $feat_object->location->end."\t".$feat_object->location->start."\t";
                                } else {
                                       print FTAB $feat_object->location->start."\t".$feat_object->location->end."\t";
                                }
				if($feat_object->has_tag("product")) {
                                	for my $product ($feat_object->get_tag_values("product")) {
                                        	print FTAB "$product\n";
                                	}
				} else {
					print FTAB "NA\n";
				}
                        }
                }
        }
        close(FTAB);
	return($annot_file);
}

####################################################################################################################################################

# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing. 
# Parameters    : NA
# Returns       : NA
# Modifications :

sub checkCmdLineArgs {
	my @required = qw(ref_genbank snp_positions blast_list output_dir);
	foreach my $option(@required) {
		if(!defined($cmdLineArgs{$option})) {
			printLogMsg($ERROR, "ERROR! : Required option $option not passed");
		}
	}
	if(exists($cmdLineArgs{'log'})) {
		open($logfh, "> $cmdLineArgs{'log'}") or printLogMsg($ERROR, "Could not open $cmdLineArgs{'log'} file for writing.\nReason : $!");
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
		die "$msg\n" if($level == $ERROR);
	}	
}

__END__

#####################
# POD DOCUMENTATION #
#####################

=head1 NAME

# Name of the script and a 1 line desc

=head1 SYNOPSIS

# USAGE : 

	parameters in [] are optional

=head1 OPTIONS



=head1 DESCRIPTION



=head1 INPUT



=head1 OUTPUT



=head1 AUTHOR

	Sonia Agrawal
	Bioinformatics Software Engineer II
	Institute for Genome Sciences
	School of Medicine
	University of Maryland
	Baltimore, MD - 21201
	sagrawal@som.umaryland.edu

==cut
