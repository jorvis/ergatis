#!/usr/local/bin/perl -w

=head1 NAME

 	replace_defline_with_locusID.pl - script that replaces the FastA header with the locus ID of that protein or nucleotide.

=head1 SYNOPSIS

# USAGE : 

	perl replace_defline_with_locusID.pl -f file.fsa -b file2.bsml -o output.fsa

=head1 OPTIONS

	--bsml_list (-b)	-	BSML list file for BSML files with locus IDs (typically output from the pipeline_summary Ergatis component)

	--fasta (-f)	-	The fasta file whose headers need to be replaced

	--output (-o)	-	The path of the output file

	--log (-l)	-	writes to a log file

	--help (-h)	-	displays this documentation

=head1 DESCRIPTION

	This script will map the identifier within the header of the FastA file to the ID attribute of the Feature element within the BSML file.  Using the start and end positions located within the Interval-loc element, the "identifier" attribute (locus_id) from within the "Cross-reference" element of the "gene" class will be used to replace the original FastA header.

	The purpose of this is to use the locus identifier as the ID within the blastable database that will be accessed via Manatee

=head1 INPUT

	A multi-fasta file
	A BSML list for BSML files that includes a Feature element with class "gene".  Withiin this class must be a "Cross-reference" element with "identifier" attribute

=head1 OUTPUT

	A multi-fasta file

=head1 CONTACT

	Shaun Adkins
	sadkins@som.umaryland.edu

=cut

use strict;
use warnings;
# Extended processing of command line options
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
# Prints usage message from embedded pod documentation
use Pod::Usage;
# Include packages here 
use XML::Twig;
use Data::Dumper;
use File::OpenFile qw(open_file);
$|++;

#############
# CONSTANTS #
#############
sub map_locus_to_uniquename();
sub replace_fasta_defline();

my %locus_id;
my %feature_index;

my $list;
my $output;
my $fasta_in;
my @input_files;

###########
# GLOBALS #
###########
my %options;
# Log file handle;
my $logfh;
my ($ERROR, $WARN, $DEBUG) = (1, 2, 3);

################
# MAIN PROGRAM #
################
GetOptions(\%options,
	   'bsml_list|b=s',
	   'fasta|f=s',
	   'output|o=s',
	   'log|l=s', 
	   'help|h'
	  ) or pod2usage();

pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ) if ($options{'help'});

&check_parameters();

my $in = open_file( $list, 'in' );	#parse the bsml list and store each bsml into an array
chomp( my @tmp = <$in> );
close($in);
push( @input_files, @tmp );

map_locus_to_uniquename();
replace_fasta_defline();

###############
# SUBROUTINES #
###############

sub map_locus_to_uniquename() {

# use XML::Twig to parse out entries from the Feature-table element from the BSML file
    my $twig = new XML::Twig( 'twig_roots' => {
        'Feature-table' => \&get_sequence_features
    });
    foreach my $bsml (@input_files) {	#iterate through the bsml files and parse
        my $fh= open_file( $bsml, "in" );
        $twig->parse( $fh );
    }

}

sub get_sequence_features {
    my ($twig, $elem) = @_;

    my @features = $elem->children('Feature');	#Access each Feature element within Feature-table
    my $total_features = scalar(@features);
    printLogMsg(1, "There are no features present in this BSML file\n") if( $total_features == 0 );
    
    foreach my $feature ( @features ) {
        my $interval_loc = &get_interval_loc( $feature );
        my @cross_references = &get_cross_references( $feature );

        $feature_index{ $feature->att('id') } = {
            'interval_loc' => $interval_loc,
            'cross-references' => \@cross_references
        };

#  Making a hash of locus identifiers with start and end coords
	my $id = $feature_index{ $feature->att('id') };
	if (scalar(@{$id->{'cross-references'}}) != 0) {
	    $locus_id{ $id->{'cross-references'}->[0]->{'identifier'} }{'start'} = $id->{'interval_loc'}->{'startpos'};
	    $locus_id{ $id->{'cross-references'}->[0]->{'identifier'} }{'end'} = $id->{'interval_loc'}->{'endpos'};
    	}
    }
}

sub get_cross_references {
    my ($elem) = @_;
    my @retval;
# map Cross-reference element attributes to the retval array and return
    map { my $atts = $_->atts; delete( $atts->{'id'} ); push(@retval, $atts); }  $elem->children('Cross-reference');	 
    return @retval;	# This contains the locus_id as Cross-reference->[0]->{'identifier'}
}

sub get_interval_loc {
    my ($elem) = @_;
    my $retval = {};
    my $int_loc = $elem->first_child('Interval-loc');
    $retval = $int_loc->atts() if( $int_loc );
    return $retval;
}

sub replace_fasta_defline() {
    open FA, "<$fasta_in" or die "Can't open $fasta_in for reading: $!\n";
    open OUT, ">$output" or die "Can't open $output for writing: $!\n";

    while (<FA>) {
	my $line = $_;
	chomp $line;
	if ($line =~ />(\w+\.\w+\.\d+\.\d+)/) {	#the unique name within the header
	    my $id = $1;
#compare the start and end coords of the ID with those mapped in the locus_id hash
	    foreach my $l_id (keys %locus_id) {
		if ($locus_id{$l_id}->{'start'} == $feature_index{$id}->{'interval_loc'}->{'startpos'}
				&& $locus_id{$l_id}->{'end'} == $feature_index{$id}->{'interval_loc'}->{'endpos'} ) {
		    $line =~ s/>$id.*/>$l_id/;	#replace the unique name with the locus identifier
									#encounted some cases were the input unique name was printed twice on the FastA header
		    last;
		}  
	    }
	}
#we print all header and sequence lines to the new file at the end of each line being read.
	print OUT $line, "\n";
    }

    close FA;
    close OUT;
}



# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing. 
# Parameters    : NA
# Returns       : NA
# Modifications :

sub check_parameters{

	if(exists($options{'log'})) {
		open($logfh, "> $options{'log'}") or die "Could not open $options{'log'} file for writing: $!\n"
	}
		
	if ( exists($options{'bsml_list'}) && exists($options{'fasta'}) && exists($options{'output'}) ) {
		$list = $options{'bsml_list'};
		$output =  $options{'output'};
		$fasta_in = $options{'fasta'};
	} else {
		printLogMsg(1, "Options --bsml_list, --fasta, and --output are required\n");
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

