#!/usr/bin/perl

=head1 NAME

check_pipeline_template_configs.pl - Verifies software used by the Ergatis pipeline
	is installed in the specified path.

=head1 SYNOPSIS

USAGE: create_evidence_file_mapping.pl
    --pipeline_layout=/path/to/pipeline.layout
    --software_config=/path/to/software.config/
  [ --help ]

=head1 OPTIONS

B<--pipeline_layout,-p>

    REQUIRED. Path to pipeline.layout file.

B<--software_config,-s>

    REQUIRED Path to Ergatis software.config file

B<--help,-h>

    Print this message

=head1  DESCRIPTION

    This program checks the program layout file to determine what components are
    required for this pipeline.  Next for each component the program finds each
    software path or data path associated and verifies that the paths exist or
    dies if they aren't.
 
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
sub collect_unique_components();
sub get_options();
sub _pod($$);

my $software_config;
my $pipeline_layout;
my @components;
get_options; 


collect_unique_components;
check_for_software ($software_config);
print "Everything came out OK\n";
exit 0;

#Runs through the pipeline layout file and adds unique components to an array
sub collect_unique_components() {
    my $pl_line;
    my %seen = ();
    open PL, "<$pipeline_layout" or die "Can't open $pipeline_layout: $!\n";
    while(<PL>) {
	chomp;
	$pl_line = $_;
	#print "$pl_line\n";
	while ($pl_line =~ /<name>([^\s\.]+)\.[^\s\.]+<\/name>/g) {
	    #print "BEF: $1\n";
	    push (@components, $1) unless $seen{$1}++;	#add element to array if it isn't in array already	
	}
    }
    close PL;
}

#Uses component array to check software config file and verify paths to software
#	the component uses are still valid
sub check_for_software($) {
    my $software_config = shift;
    my $sc_line;
    my $sc_line2;
    
    foreach my $component (@components) {
	#print "AFTER: $component\n";
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
			#print "PATH: $1\n";
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
			      'pipeline_layout|p=s',
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

    if ($options{'pipeline_layout'}) {
	$pipeline_layout = $options{'pipeline_layout'};
        die("File does not exist: $pipeline_layout")
            unless( -e $pipeline_layout && -f $pipeline_layout );
    } else {
	_pod("The pipeline.layout filepath for the current pipeline is required.",1);
    }
   
}

sub _pod($$) {
    my ($message, $exitval) = @_;
    $message = " " unless( $message );
    $exitval = 0 unless( $exitval );
    die;
}
