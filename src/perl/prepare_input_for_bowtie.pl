#!/usr/bin/env perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";


=head1 NAME

prepare_input_for_bwotie.pl - Collates all input files together and creates a
                              single input file list passed to the bowtie 
                              component

=head1 SYNOPSIS

USAGE: ./prepare_input_for_bowtie.pl --paired_end_reads=/name/of/paired/end/tags
                                     --non_paired_end_reads=/name/of/non/paired/end/tags
                                     --output_file=/path/to/desired/outupt/file
                                     
=head1 OPTIONS

B<--paired_end_reads, -p>
    A list of one or many (comma-separated) tags that point to paired-end sequence data.
    
B<--non_paired_end_reads, -n>
    A set of file lists containing non-paired end data; one file per line.
    
B<--output_file, -o>
    The desired output file collating all reads into one list.
    
B<--help, -h>
    Print perldocs for this script.
    
=head1 DESCRIPTION

This file exists to prepare input data that will be fed into the bowtie component 
when running under the CLoVR VM appliance. Takes a list of tags for both 
paired-end data and non-paired-end data and collates them into one large file list
to be passed to the bowtie component. 

Special attention should be paid to how paired-end data is handled. A paired-end
data set should have both files (pairs) tagged as one tag. A new tag should be created
for any subsequent paired-end data sets.

 It is not advised to run this script outside of the CloVR environment. Al                                                                           

 =head1 INPUT

 A comma-separated list of tag names for both paired-end data and non-paired-end data.

 =head1 OUTPUT

 One file list combining both paired-end data and non-paired end data:

/tmp/sequence.fastq
/tmp/sequence1.fastq,/tmp/sequence2.fastq

=head1 CONTACT

    Cesar Arze
    carze@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Basename;
use Ergatis::Logger;
umask(0000);
my $logger;

my %options = &_parse_options();
my $paired_end_tag_list = $options{'paired_end_reads'};
my $non_paired_tag_list = $options{'non_paired_end_reads'};
my $output_file = $options{'output_file'};

my @merged_files = merge_file_lists($paired_end_tag_list, $non_paired_tag_list);
print_new_file_list($output_file, @merged_files);

#########################################################################
#                                                                       #
#                           SUBROUTINES                                 #
#                                                                       #
#########################################################################

#--------------------------------------------------------------------- 
# Creates a new file list containing all files found in the two
# sets of tag lists passed in.
#---------------------------------------------------------------------
sub print_new_file_list {
    my ($out_file, @seq_files) = @_;
    open ("SEQOUT", "> $out_file") or $logger->logdie("Could not write to $out_file: $!");

    ## TODO: Make use of xargs here so we can avoid an overflow issue if too many files
    ## have been tagged.
    foreach my $file (@seq_files) {
        chomp($file);
        # In order for this file not to fail in CloVR we need to replace the comma 
        # with a space (rsync can handle a space sep'd string but not comma-sep'd)
        $file =~ s/,/ /; 
        print SEQOUT $file . "\n";
    }

    close (SEQOUT);
}
    
#--------------------------------------------------------------------- 
# Iterates over both list of tags and extracts all the files necessary.
# Handles paired-end data if present.  
#---------------------------------------------------------------------
sub merge_file_lists {
    my ($pair_tags, $nonpair_tags) = @_;
    my @files = ();

    # Our paired-end data needs to be handled to make sure that all
    # pairs are on the same line delimited by a comma:
    #
    # sequence_1_1.fastq,sequence_1_2.fastq
    # sequence_2_1.fastq,sequence_2_2.fastq
    push(@files, convert_tags_list_to_file_list($pair_tags, 1)) if (defined($pair_tags));
    push(@files, convert_tags_list_to_file_list($nonpair_tags, 0)) if (defined($nonpair_tags));

    _verify_files(@files);
    return @files;
}

sub convert_tags_list_to_file_list {
    my ($tag_str, $is_paired) = @_;
    my @pair_list = ();
    my @tags = ();

    if ($tag_str =~ /,/) {
        @tags = split(/,/, $tag_str);
    } else {
        push(@tags, $tag_str);
    }
                                
    foreach my $tag (@tags) {
        my $cmd = "vp-describe-dataset --tag-name=$tag | grep FILE | cut -f 2";
        my $files = _run_system_cmd($cmd);

        if ($is_paired) {
            $files =~ s/\n/,/;
            push(@pair_list, $files);
        } else {
           push(@pair_list, split(/\n/, $files));
        }
    }

    return @pair_list;
}

#--------------------------------------------------------------------- 
# Parse command-line arguments                                       
#---------------------------------------------------------------------
sub _parse_options {
    my %opts = ();

    GetOptions(\%opts,
                'paired_end_reads|p=s',
                'non_paired_end_reads|n=s',
                'output_file|o=s',
                'help' ) || pod2usage();

    if ($opts{'help'}) {
        pod2usage ( { -exitval => 0, -verbose => 2, -output => \*STDERR } );
    }
    
    my $logfile = Ergatis::Logger::get_default_logfilename();
    my $debug = 4;
    $logger = new Ergatis::Logger( 'LOG_FILE'   =>  $logfile,
                                   'LOG_LEVEL'  =>  $debug );
    $logger = Ergatis::Logger::get_logger();
    

    unless ((defined $opts{'paired_end_reads'}) || (defined $opts{'non_paired_end_reads'})) {
        $logger->logdie("Please provide one or more CloVR tags (comma-separated) containing fastq data file");
    }
    defined ($opts{'output_file'}) || $logger->logdie("Please specify an output file");

    return %opts
}

#--------------------------------------------------------------------- 
# Run system command verifying if the command completed succesfully
# and returns output written to STDOUT
#---------------------------------------------------------------------
sub _run_system_cmd {
    my $cmd = shift;

    my $output = `$cmd`;
    my $ret_val = $?;

    if ($ret_val >> 8 != 0) {
        $logger->logdie("Could not execute command $cmd: $!");
    }

    return $output;
}

#--------------------------------------------------------------------- 
# Verifies that all files in the passed in array exist and are 
# readable and non-zero
#---------------------------------------------------------------------
sub _verify_files {
    my @files = @_;

    $logger->logdie("No files were found in provided tags") if (scalar @files == 0);
    
    foreach my $file (@files) {
        if ($file =~ /,/) {
            my @mates = split(/,/, $file);
            _verify_files(@mates);
            next;
        } 

        chomp($file);
        unless (-e $file) {
            $logger->logdie("File $file does not exist: $!");
        }

        unless (-r $file) {
            $logger->logdie("File $file is not readable: $!");
        }

        unless (-s $file) {
            $logger->logdie("File $file is zero-size: $!");
        }
    }
}
