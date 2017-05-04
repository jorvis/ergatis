#!/usr/bin/env perl

=head1 NAME

orthomcl_adjust_fasta.pl - Wrapper for the orthomclAdjustFasta script

=head1 SYNOPSIS

 USAGE: orthomcl_adjust_fasta.pl
       --input_file=/path/to/some/input.file
	   --input_list=/path/to/some/input.list
       --output_dir=/path/to/output/dir
     [ --log=/path/to/file.log
	   --mapping_file=/genome/mapping/file.txt
	   --unique_id_field=1
	   --orthomcl_bin=/path/to/orthoMCL/bin/
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>

B<--input_file,-I>

B<--output_dir,-o>

B<--mapping_file,-m>
	Tab-delimited file with two fields: 1) basename of fasta file 2) Unique ID prefix to append to adjusted FASTA files
	If not provided, each ID prefix will be "ID" + an incremented number

B<--unique_id_field, -u>
	Field position (first position is 1) that has the unique ID for each protein.  Fields are separated by whitespace or '|'
	Default is 1

B<--orthomcl_bin -b>
	Path to orthoMCL bin directory.  Default is /usr/local/packages/orthomcl/bin

B<--log,-l>
    Logfile.

B<--debug,-d>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

 DESCRIPTION
 
=head1  INPUT

    Describe the input

=head1 OUTPUT

    Describe the output

=head1  CONTACT

    Shaun Adkins
    sadkins@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
my @files;

my $ORTHOMCL_BIN = "/usr/local/packages/orthomcl/bin";
my $ID_FIELD = 1;
my $PREFIX="ID";
my $prefix_count=0;
####################################################

my %options;
my %file2prefix;

# Allow program to run as module for unit testing if necessary (change extension to .pm)
main() unless caller();

sub main {
    my $results = GetOptions (\%options,
                         "input_file|i=s",
						 "input_list|I=s",
                         "output_dir|o=s",
						 "orthomcl_bin|b:$ORTHOMCL_BIN",
						 "mapping_file|m=s",
						 "unique_id_field|u:1",
                         "log|l=s",
                         "debug|d=s",
                         "help|h"
                          );

    &check_options(\%options);
	parse_mapping_file($options{'mapping_file'});
	# Create the working directory for the orthoMCL exe to operate in
	my $output_dir = $options{'output_dir'} . '/compliantFasta/';
	mkdir($output_dir); chdir($output_dir);
	run_orthomcl_cmd(\@files, \%file2prefix);
}

# Parse mapping file into File/Prefix value pairs
# NOTE - I wrote that the basename of the file should be the first field entry, but technically anything that can be regexed in the list of files 
sub parse_mapping_file {
	my $map_file = shift;

	open MAP, $map_file or die("Cannot open mapping file $map_file for reading: $!");
	while (<MAP>){
		chomp;
		my @fields = split('\t');
		# File => Prefix
		$file2prefix{$fields[0]} = $fields[1];
	}
	close MAP;
}

sub run_orthomcl_cmd {
	my ($file_list, $f2p) = @_;

	foreach my $file (@$file_list) {
		my $count_flag = 0;
		my $prefix = $PREFIX . $prefix_count;
		foreach my $f (keys %$f2p) {
			if ($file =~ /$f/) {
				$prefix = $f2p->{$f};
				$count_flag = 1;
				last;
			}
		}
		$prefix_count++ if ($count_flag);
		
		my $cmd = $options{'orthomcl_bin'}."/orthomclAdjustFasta $prefix $file " . $options{'unique_id_field'};
		&_log($DEBUG, "CMD - $cmd");
		system($cmd);
	}
}

sub check_options {
    my $opts = shift;
    if( $opts->{'help'} ) {
        &_pod;
    }

    if( $opts->{'log'} ) {
        open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
    }

    $debug = $opts->{'debug'} if( $opts->{'debug'} );

    foreach my $req ( qw(output_dir) ) {
        &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
    }

	if (! ($opts->{'input_file'} || $opts->{'input_list'}) ) {
		&_log($ERROR, "Either option 'input_file' or 'input_list' is required");
	}

	if ($opts->{'input_file'} && $opts->{'input_list'} ) {
		&_log($ERROR, "Please pass only one option of 'input_file' or 'input_list'... not both");
	}

	if( $opts->{'input_file'} ) {
		push @files, $opts->{'input_file'};
	}

	if( $opts->{'input_list'} ) {
    	open FH, "<".$opts->{'input_list'} or die("Error in opening the file,". $opts->{'input_list'} ." :$!\n");
		while( my $file = <FH> ) {
			chomp $file;
			if( -e $file ) { push @files, $file; }
			else { &_log($WARN, "WARNING: $file - No such file exists\n"); }
		}
		close(FH);
	}

}

sub _log {
    my ($level, $msg) = @_;
    if( $level <= $debug ) {
        print STDERR "$msg\n";
    }
    print $logfh "$msg\n" if( defined( $logfh ) );
    exit(1) if( $level == $ERROR );
}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
