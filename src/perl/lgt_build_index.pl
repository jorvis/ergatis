#!/usr/local/bin/perl -w

#########################################################################################
#											#
# Name	      : lgt_build_index.pl							#
# Version     : 1.0									#
# Project     : LGT Seek Pipeline							#
# Description : Script to build a BWA index for the reference sequence			#
# Author      : Sonia Agrawal								#
# Date        : September 2, 2015							#
#											#
#########################################################################################

use strict;

# Extended processing of command line options
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

# Prints usage message from embedded pod documentation
use Pod::Usage;

# Include packages here
use File::Basename;

#############
# CONSTANTS #
#############

my $EXT = ".fsa";

###########
# GLOBALS #
###########
my %hCmdLineArgs = ();

# Log file handle;
my $fhLog;
my ( $sFileBase, $sFileDir, $sFileExt, $sCmd, $nExitCode );
my ( $ERROR, $WARN, $DEBUG ) = ( 1, 2, 3 );

################
# MAIN PROGRAM #
################
GetOptions( \%hCmdLineArgs, 'reference|r=s', 'algo|a=s', 'bwa_path|b=s',
    'output_dir|o=s', 'log|l=s', 'help|h' )
  or pod2usage();

pod2usage( { -exitval => 0, -verbose => 2, -output => \*STDERR } )
  if ( $hCmdLineArgs{'help'} );

checkCmdLineArgs( \%hCmdLineArgs );


my @files;

if ( ( -e $hCmdLineArgs{'reference'} ) && ( -r $hCmdLineArgs{'reference'} ) ) {
	if ($hCmdLineArgs{'reference'} =~ /\.list$/){
		# If option is a list of refs, copy all refs as well as the list
		open "FH", $hCmdLineArgs{'reference'} || printLogMsg( $ERROR, "ERROR : Can't open list file for reading: $!");
		push @files, <FH>;
		( $sFileBase, $sFileDir, $sFileExt ) =
			  fileparse( $hCmdLineArgs{'reference'}, qr/\.[^.]*/ );
		$sCmd = "ln -sf " . $hCmdLineArgs{'reference'} . " " . $hCmdLineArgs{'output_dir'} . "/" . $sFileBase . $EXT;
	} else {
		# Otherwise just copy the ref
		push @files, $hCmdLineArgs{'reference'};
	}

	foreach my $ref (@files){
		chomp $ref;
		( $sFileBase, $sFileDir, $sFileExt ) =
  			fileparse( $ref, qr/\.[^.]*/ );
		$sCmd =
	        "ln -sf "
	      . $ref . " "
	      . $hCmdLineArgs{'output_dir'} . "/"
	      . $sFileBase
	      . $EXT;
	    printLogMsg( $DEBUG,
	        "INFO : Creating symlink to reference file $ref in output directory $hCmdLineArgs{'output_dir'}.\nINFO : Command : $sCmd"
	    );
	    $nExitCode = system($sCmd);
	    if ( $nExitCode == 0 ) {
	        printLogMsg( $DEBUG,
	            "INFO : Symlink to reference file $ref created in output directory $hCmdLineArgs{'output_dir'}"
	        );
	    } else {
	        printLogMsg( $ERROR,
	            "ERROR : Symlink to reference file $ref could not be created"
	        );
	    }
	    if ( !exists( $hCmdLineArgs{'algo'} ) ) {
	        $sCmd =
	            $hCmdLineArgs{'bwa_path'}
	          . " index -p "
	          . $hCmdLineArgs{'output_dir'} . "/"
	          . $sFileBase
	          . $EXT . " "
	          . $hCmdLineArgs{'output_dir'} . "/"
	          . $sFileBase
	          . $EXT;
	    } else {
	        $sCmd =
	            $hCmdLineArgs{'bwa_path'}
	          . " index -p "
	          . $hCmdLineArgs{'output_dir'} . "/"
	          . $sFileBase
	          . $EXT . " -a "
	          . $hCmdLineArgs{'algo'} . " "
	          . $hCmdLineArgs{'output_dir'} . "/"
	          . $sFileBase
	          . $EXT;
	    }
	}
} else {
    printLogMsg( $ERROR,
        "ERROR : Reference file $hCmdLineArgs{'reference'} does not exist or is not readable for creating index."
    );
}

printLogMsg( $DEBUG,
    "INFO : Starting to index passed reference $hCmdLineArgs{'reference'} using BWA.\nINFO : Command : $sCmd"
);
$nExitCode = system($sCmd);
if ( $nExitCode == 0 ) {
    printLogMsg( $DEBUG,
        "INFO : Reference file $hCmdLineArgs{'reference'} indexing completed in directory $hCmdLineArgs{'output_dir'}"
    );
} else {
    printLogMsg( $ERROR,
        "ERROR : Reference file $hCmdLineArgs{'reference'} indexing failed" );
}

###############
# SUBROUTINES #
###############

# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing.
# Parameters    : NA
# Returns       : NA
# Modifications :

sub checkCmdLineArgs {
    my ($phCmdLineArgs) = @_;
    my $sOption;
    my @aRequired = ();
    if ( exists( $phCmdLineArgs->{'log'} ) ) {
        open( $fhLog, "> $phCmdLineArgs->{'log'}" )
          or die
          "Could not open $phCmdLineArgs->{'log'} file for writing.Reason : $!\n";
    }
    @aRequired = qw(reference output_dir bwa_path);
    foreach my $sOption (@aRequired) {
        if ( !defined( $phCmdLineArgs->{$sOption} ) ) {
            printLogMsg( $ERROR,
                "ERROR! : Required option $sOption not passed" );
        }
    }
}

####################################################################################################################################################

# Description   : Used to handle logging of messages(errors and warnings) during the execution of the script
# Parameters    : level = can be ERROR, WARNING or INFO
#		  msg   = msg to be printed in the log file or to STDERR
# Returns       : NA
# Modifications :

sub printLogMsg {
    my ( $nLevel, $sMsg ) = @_;
    if ( $nLevel <= $DEBUG ) {
        print STDERR "$sMsg\n";
        print $fhLog "$sMsg\n" if ( defined($fhLog) );
        die "" if ( $nLevel == $ERROR );
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

	<-r|--reference>	:	Path to reference fasta file to be indexed

	<-b|--bwa_path>		:	Path to BWA executable

	<-a|--algo>			:	Name of Algorithm to be used
							Use 'is' for genome database sizes projected to be < 2Gb
							Use 'btwsw' for larger genome databases
	
	<-o|--output_dir>	:	Path to output directory

=head1 DESCRIPTION



=head1 INPUT

A fasta file of the reference sequence to be indexed is required

=head1 OUTPUT

An indexed fasta sequence (and supporting files) are created in the output directory specified

=head1 AUTHOR

	Sonia Agrawal
	Bioinformatics Software Engineer II
	Institute for Genome Sciences
	School of Medicine
	University of Maryland
	Baltimore, MD - 21201
	sagrawal@som.umaryland.edu

==cut
