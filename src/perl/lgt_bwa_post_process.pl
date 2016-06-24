#!/usr/bin/env perl
=head1 NAME

lgt_bwa_post_process.pl - Description

=head1 SYNOPSIS

 USAGE: lgt_bwa_post_process.pl
       --input_file=/path/to/some/input.file
       --output=/path/to/transterm.file
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--donor_file,-d>
	Filtered BAM file mapping to a donor reference genome

B<--donor_file_list, -D>
	List of BAM files mapping to a donor reference genome.  Filename must have a shared identifier with the recipient genome

B<--recipient_file, -r>
	Filtered BAM file mapping to a recipient or host reference genome

B<--recipient_file_list, -R>
	List of BAM files mapping to a recipient reference genome.  Filename must have a shared identifier with the donor genome

B<--output_dir,-o>

B<--prefix,-p>
	Prefix name to give output files

B<--samtools_path,-s>
    Path to samtools executable

B<--log,-l>
    Logfile.

B<--debug,-d>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

 	This is a wrapper script that will call LGT::LGTSeek->bwaPostProcess to handle the post processing.

=head1  INPUT

    Either:
    a) Two filtered BAM files.  One mapping to a donor reference genome and the other to the recipient
    or:
    b) Two lists of filtered BAM files.  One list has files mapped to donor reference genomes and the other has them mapped to the recipient

=head1 OUTPUT

    Three new BAM files

    1) Microbiome file - Both ends of the paired reads mapped to the donor genome but both ends did not map to the recipient

    2) LGT_donor file
    3) LGT_recipient file
    - Potential lateral gene transfer (LGT) is defined as pairs of reads where exactly one read maps to the donor (other is unmapped), and the other maps to the recipient (first is unmapped).  Thus we output alignments for both the donor and recipient.

=head1  CONTACT

    Shaun Adkins
    sadkins@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Basename;
use LGT::LGTSeek;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
####################################################
my $prefix;
my @donor_files;
my @recipient_files;
my $donor_only = 0;
my $host_only = 0;
my %LGT_groups;
my $output_dir;
my %options;
my $count_id = 0;

my $results = GetOptions (\%options,
                     "donor_file|d=s",
                     "donor_file_list|D=s",
                     "recipient_file|r=s",
                     "recipient_file_list|R=s",
                     "samtools_path|s=s",
                     "output_dir|o=s",
					 "prefix|p=s",
                     "log|l=s",
                     "debug|d=s",
                     "help|h"
                      );

&check_options(\%options);

$prefix = (defined $options{'prefix'}) ? $options{'prefix'} : "post_process";

# Need to take a different approach depending on if we have a host file or not.
if ($donor_only || $host_only) {
	# First, lets give our input files a generic name to avoid being redundant with code
	my @files = @donor_files if $donor_only;
	@files = @recipient_files if $host_only;

	foreach my $file (@files){
		# Assign a number to the output file to make it unique
		my $curr_prefix = $prefix;
		$curr_prefix = $prefix .  "_" . $count_id++ if (scalar @files > 1);

		# Instatitate an LGTSeek object
		my $lgt_obj = LGT::LGTSeek->new( {
				'samtools_bin' => $options{'samtools_path'},
				'output_dir' => $options{'output_dir'},
				'verbose' => 1
			} );

		# Call bwaPostProcess using the right config, depending on donor-only or host-only
		my $pp_data = $lgt_obj->bwaPostProcess( {
				'donor_bam'	=> $file,
				'output_prefix' => $curr_prefix
			} ) if $donor_only;
		$pp_data = $lgt_obj->bwaPostProcess( {
				'host_bam'	=> $file,
				'output_prefix' => $curr_prefix
			} ) if $host_only;

		# Create the list of counts
		my $counts_file = $output_dir . "/" . $curr_prefix . ".counts";
		my (@header, @vals);
	    map {
	        push( @header, $_ );
	        my $foo = $pp_data->{counts}->{$_} ? $pp_data->{counts}->{$_} : 0;
	        push( @vals, $foo );
	    } ( 'total', 'paired', 'single', 'none' );
		LGT::LGTSeek->print_tab( $counts_file, \@header, \@vals );
	}

} else {
	# Before doing post-processing, ensure donor and host file counts are equal, and match them up
	my $matching_files = find_matching_files(\@donor_files, \@recipient_files);
	# Group reads from donor/recipient files for each mapping type
	foreach my $r (keys %$matching_files) {
		my $curr_prefix = $prefix;
		$curr_prefix = $prefix . "_" . $count_id++ if (scalar keys %$matching_files > 1);
	
		my $lgt_obj = LGT::LGTSeek->new( {
				'samtools_bin' => $options{'samtools_path'},
				'output_dir' => $options{'output_dir'},
				'verbose' => 1
			} );
		
		my $d = $matching_files->{$r};
		# Returns a hash of group assignemt counts and the files, but not really necessary to bring out here.
		my $pp_data = $lgt_obj->bwaPostProcess( {
				'donor_bam'	=> $d,
				'host_bam'	=> $r,
				'output_prefix' => $curr_prefix
			} );
		
		# Take list of counts and write them to file
		my $counts_file = $output_dir . "/" . $curr_prefix . ".counts";
		my (@header, @vals);
	    map {
	        push( @header, $_ );
	        my $foo = $pp_data->{counts}->{$_} ? $pp_data->{counts}->{$_} : 0;
	        push( @vals, $foo );
	    } ( 'total', 'host', 'no_map', 'all_map', 'single_map', 'integration_site_host', 'integration_site_donor', 'microbiome', 'lgt' );
		LGT::LGTSeek->print_tab( $counts_file, \@header, \@vals );
	};
}

exit(0);

# Finds the matching donor and recipient files and returns in a hash
sub find_matching_files {
    my ($d_arr, $r_arr) = @_;
    my %m_files;

	# The donor BAMs were created from recipient BAMs, so the names are longer.
	# So we need to grep for the recipient basename in the donor file

    if (scalar @$d_arr == 1 && scalar @$r_arr == 1) {
        $m_files{$r_arr->[0]} = $d_arr->[0];
        return \%m_files;
    }

    foreach my $r_file (@$r_arr){
        # TODO: May need to refine this extension pattern
        my ($r_base, $r_dir, $r_ext) = fileparse($r_file, qr/\.\w*\.?bam/);
        my @grepped = grep {$_ =~ /$r_base/} @$d_arr;
        &_log($ERROR, "Found more than 1 potential donor BAM file match for the recipient BAM file $r_file.  File basenames need to be 1-to-1 matching : $!") if scalar @grepped > 1;
        &_log($ERROR, "Found no match in list of donor BAM files to the recipient BAM file $r_file.  Check the file basenames and ensure match is present : $!") if scalar @grepped < 1;
        $m_files{$r_file} = $grepped[0];

    }
    return \%m_files;
}

# Process read options
sub check_options {
   my $opts = shift;
   if( $opts->{'help'} ) {
       &_pod;
   }

   if( $opts->{'log'} ) {
       open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
   }

   $debug = $opts->{'debug'} if( $opts->{'debug'} );

   foreach my $req ( qw(samtools_path output_dir) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }

   # Gather files together
   if($opts->{donor_file}) {
       push @donor_files, $opts->{donor_file};
   } elsif($opts->{donor_file_list}) {
       @donor_files = `cat $opts->{donor_file_list}`;
   }

    if($opts->{recipient_file}) {
       push @recipient_files, $opts->{recipient_file};
    } elsif($opts->{recipient_file_list}) {
       @recipient_files = `cat $opts->{recipient_file_list}`;
    }
    if (scalar @recipient_files == 0) {
		print STDOUT "Assuming this is a donor-file only run.\n";
		$donor_only = 1;
	}
    if (scalar @donor_files == 0) {
		print STDOUT "Assuming this is a host-file only run.\n";
		$host_only = 1;
	}
   	if ( !($host_only || $donor_only) && scalar @donor_files != scalar @recipient_files) {
       	&_log($ERROR, "ERROR : Number of donor files to recipient files is not equal.  Please check both lists. $!");
	}
   $output_dir = $opts->{'output_dir'};
}

sub _log {
   my ($level, $msg) = @_;
   if( $level <= $debug ) {
      print STDOUT "$msg\n";
   }
   print $logfh "$msg\n" if( defined( $logfh ) );
   exit(1) if( $level == $ERROR );
}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
