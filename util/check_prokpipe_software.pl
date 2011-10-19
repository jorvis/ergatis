#!/usr/bin/perl

=head1 NAME

check_pipeline_template_configs.pl - Verifies software used by the Ergatis 
	prokaryotic pipeline is installed in the specified path.

=head1 SYNOPSIS

USAGE: create_evidence_file_mapping.pl
    --software_config=/path/to/software.config/
  [ --help ]

=head1 OPTIONS

B<--software_config,-s>

    REQUIRED Path to Ergatis software.config file

B<--help,-h>

    Print this message

=head1  DESCRIPTION

    This program checks the Ergatis software.config file to verify if the software
    file path associated with a given prokaryotic pipeline component exists and is
    installed.  Note this is run before a pipeline is instantiated.
 
=head1  INPUT
    
    Just the pipeline layout file and the Ergatis software.config file

=head1 OUTPUT

    No output files..just says "Everything is all OK" if all software filepaths are correct

=head1  CONTACT

    Shaun Adkins
    sadkins@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

sub check_for_software($);
sub get_options();
sub _pod($$);

my $software_config;
my @components;
get_options; 

@components = qw(ncbi-blastx ncbi-blastp glimmer3 hmmpfam3 lipop 
		create_pseudomolecules ber prodigal ps_scan RNAmmer signalp 
		tmhmm translate_sequence tRNAscan-SE formatdb p_func 
		parse_evidence pipeline_summary train_for_glimmer3_iteration);

#print $components[0], "\n";

check_for_software ($software_config);
print "Everything came out OK\n";
exit 0;

#Uses component array to check software config file and verify paths to software
#	the component uses are still valid
sub check_for_software($) {
    my $software_config = shift;
    my $sc_line;
    my $sc_line2;
    
    foreach my $component (@components) {
	print "$component\n";
	open SH, "<$software_config" or die "Can't open $software_config: $!\n";
	while (<SH>) {
	    chomp;
	    $sc_line = $_;
	    #print "$sc_line\n";
	    if ($sc_line =~ /^\[component $component\]/) {	
		while (! (($sc_line2 = <SH>) =~ /^\s*$/)) {	
		#read until file has a blank line (which separates components)
		    chomp $sc_line2;
		    #print "$sc_line2\n";
		    if ($sc_line2 =~ /=(\/.+)$/) {
			print "PATH: $1\n";
			die ("Software $1 for component $_ cannot be found: $!\n")
				unless (-e $1);
		    }
		}
		next;
	    }
	}
	close SH   
    }

#check for a few extra files that are listed under 'common inputs'
    open SH, "<$software_config" or die "Can't open $software_config: $!\n";
    while(<SH>) {
	chomp;
	if ($_ =~ /\$;(HMM_ALL|HMM_LIB_DB|DB_UNIPROT_100|UNIREF100_LOOKUP_FILE)\$;=(\/.+)$/) {
	    #print "$_\n";
	    #print "$1\t$2\n";
	    die ("Input $2 cannot be found: $!\n")
		unless (-e $2);	
	}
    }
    close SH;
}

sub get_options() {
    my %options = ();
    my $results = GetOptions (\%options, 
                              'software_config|s=s',
                              'help|h');

    if( $options{'help'} ){
        pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
    }

    if( $options{'software_config'} ) {
        $software_config = $options{'software_config'};
        die("File does not exist: $software_config")
            unless( -e $software_config && -f $software_config );
    } else {
        _pod("The software.config filepath is required.",1);
    }
   
}

sub _pod($$) {
    my ($message, $exitval) = @_;
    $message = " " unless( $message );
    $exitval = 0 unless( $exitval );
    die;
}
