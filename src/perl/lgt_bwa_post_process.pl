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

 DESCRIPTION

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
    - Potential lateral gene transfer (LGT) is defined as pairs of reads where exactly one read maps to the donor (other is unmapped),
        and the other maps to the recipient (first is unmapped).  Thus we output alignments for both the donor and recipient.
=head1  CONTACT

    Shaun Adkins
    sadkins@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Basename;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
####################################################
my $prefix;
my @donor_files;
my @recipient_files;
my %LGT_groups;
my $output_dir;
my %options;

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

my $matching_files = find_matching_files(\@donor_files, \@recipient_files);

# Group reads from donor/recipient files for each mapping type
foreach my $r (keys %$matching_files) {
    my $d = $matching_files->{$r};
    classify($d, "donor");
    classify($r, "recipient");

    # Grab headers of BAM files
    my $d_head = `$options{'samtools_path'} view -H $d`;
    my $r_head = `$options{'samtools_path'} view -H $r`;
    # TODO:  Append additional information to headers if necessary (ID step, PG program VN version
	
	$prefix = $options{'prefix'} if (defined $options{'prefix'});

    # Open some filehandles to be used, and print headers
    open (my $lgtd, "| $options{'samtools_path'} view -S -b -o $output_dir/$prefix.lgt_donor.bam -") || &_log($ERROR, "Unable to open outfile");
    print $lgtd $d_head . "\n";
    open (my $lgtr, "| $options{'samtools_path'} view -S -b -o $output_dir/$prefix.lgt_recipient.bam -") || &_log($ERROR, "Unable to open outfile");
    print $lgtr $r_head . "\n";
    open (my $mb_donor, "| $options{'samtools_path'} view -S -b -o $output_dir/$prefix.microbiome.bam -") || &_log($ERROR, "Unable to open outfile");
    print $mb_donor $d_head . "\n";

    # print new BAM files based on classification groups
    foreach my $id (keys %LGT_groups){
        if ($LGT_groups{$id}{'donor'}{'group'} = 'MM' && $LGT_groups{$id}{'recipient'}{'group'} = 'UU') {
            print $mb_donor $LGT_groups{$id}{'donor'}{'line1'} . "\n";
            print $mb_donor $LGT_groups{$id}{'donor'}{'line2'} . "\n";
        }
        if (($LGT_groups{$id}{'donor'}{'group'} = 'UM' && $LGT_groups{$id}{'recipient'}{'group'} = 'MU') ||
            ($LGT_groups{$id}{'donor'}{'group'} = 'MU' && $LGT_groups{$id}{'recipient'}{'group'} = 'UM')) {
            print $lgtd $LGT_groups{$id}{'donor'}{'line1'} . "\n";
            print $lgtd $LGT_groups{$id}{'donor'}{'line2'} . "\n";
            print $lgtr $LGT_groups{$id}{'recipient'}{'line1'} . "\n";
            print $lgtr $LGT_groups{$id}{'recipient'}{'line2'} . "\n";
        }
    }

    %LGT_groups = ();
    close $lgtd;
    close $lgtr;
    close $mb_donor;
}


sub classify {
    my $file = shift;
    my $reference = shift;
    # Open the sorted BAM file for reading
    open(my $bam_fh, "$options{'samtools_path'} view " . $file." |") or &_log($ERROR, "ERROR : main :: Could not open BAM file $file for reading.\nReason : $!");
    my $read_1;
    my $first_line_read = 0;
    while(my $line = <$bam_fh>) {
    	chomp($line);
    	# Keeping track if current line is first or second mate pair
    	if($first_line_read == 0) {
    		$read_1 = $line;
    		$first_line_read = 1;
    		next;
    	}
    	$first_line_read = 0;
    	my $read_2 = $line;
    	my ($query_1, $bit_flag_1, $cigar_1) = (split /\t/, $read_1)[ 0, 1, 5];
    	my ($query_2, $bit_flag_2, $cigar_2) = (split /\t/, $read_2)[ 0, 1, 5];

    	my $stat_r1 = parse_flag($bit_flag_1);
    	my $stat_r2 = parse_flag($bit_flag_2);

    	# Want to keep UU, MU, and UM reads.  Keep MM if specified (reference genome is donor for example)
        if(! $stat_r1->{'qunmapped'} && ! $stat_r2->{'qunmapped'}) {
    		$LGT_groups{$query_1}{$reference}{'group'} = "MM"; # Both reads have the same name
    	}
        if ($stat_r1->{'qunmapped'} || $stat_r2->{'qunmapped'}) {
    	    if (! $stat_r1->{'qunmapped'} && $stat_r2->{'qunmapped'}) {
                 $LGT_groups{$query_1}{$reference}{'group'} = "MU";
             }
    	    if ($stat_r1->{'qunmapped'} && ! $stat_r2->{'qunmapped'})
            {
                $LGT_groups{$query_1}{$reference}{'group'} = "UM";
            }
    	    if ($stat_r1->{'qunmapped'} && $stat_r2->{'qunmapped'})
            {
                $LGT_groups{$query_1}{$reference}{'group'} = "UU";
            }
    	}
        # Store read lines into hash so we don't have to iterate to retrieve again later
        $LGT_groups{$query_1}{$reference}{'line1'} = $read_1;
        $LGT_groups{$query_1}{$reference}{'line2'} = $read_2;
    }
    close $bam_fh;
}

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
        $prefix = $r_base if (! $options{'prefix'});
        my @grepped = grep {$_ =~ /$r_base/} @$d_arr;
        &_log($ERROR, "Found more than 1 potential donor BAM file match for the recipient BAM file $r_file.  File basenames need to be 1-to-1 matching : $!") if scalar @grepped > 1;
        &_log($ERROR, "Found no match in list of donor BAM files to the recipient BAM file $r_file.  Check the file basenames and ensure match is present : $!") if scalar @grepped < 1;
        $m_files{$r_file} = $grepped[0];

    }
    return \%m_files;
}

# Parse the bitwise flag from SAM output
sub parse_flag {
	my $flag = shift;
	my ($left_bin, $bin);
	$left_bin = unpack("B32", pack("N", $flag));
	$left_bin =~ s/^0+(?=\d)//;    # otherwise you'll get leading zeros
	$bin = sprintf("%012d", $left_bin);
	my $right_bin = reverse($left_bin);
    my $bit_hash = {
			'paired' => substr($right_bin, 0, 1),
        	'propermap' => substr($right_bin, 1, 1),
        	'qunmapped' => substr($right_bin, 2, 1),
        	'munmapped' => substr($right_bin, 3, 1),
        	'qrev' => substr($right_bin, 4, 1),
        	'mrev' => substr($right_bin, 5, 1),
        	'firstpair' => substr($right_bin, 6, 1),
        	'secondpair' => substr($right_bin, 7, 1),
        	'scndryalign' => substr($right_bin, 8, 1),
        	'failqual' => substr($right_bin, 9, 1),
        	'pcrdup' => substr($right_bin, 10, 1),
        	'supplealign' => substr($right_bin, 11, 1)
	};
	return $bit_hash;
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

   if (scalar @donor_files != scalar @recipient_files) {
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
