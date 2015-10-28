#!/usr/bin/env perl
package run_prinseq;    # Change this to reflect name of script.  If you want to use as module, make sure to save as .pm instead of .pl

=head1 NAME

run_prinseq.pl - Run Prinseq to filter out duplicate and low-complexity sequences

=head1 SYNOPSIS

 USAGE: perl_template.pl
       --input_file=/path/to/some/input.file
       --input_base=/path/to/some/base
       --output_dir=/path/to/output/dir
	   --prinseq_path=/path/to/prinseq.pl
	   --samtools_path=/path/to/samtools/exec
	   --picard_path=/path/to/picard/dir
     [ --rm_duplicates
	   --rm_low_complexity
       --lc_method='dust'
       --lc_threshold=7
       --tmp_dir/path/to/scratch
	   --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>

B<--input_base,-I>
    Path to files that share a given basename.  Useful for grabbing paired fastq sequences

B<--output_dir,-o>

B<--tmp_dir,-t>

B<--prinseq_path,-p>
	The path to the prinseq-lite Perl script

B<--samtools_path,-s>
	The path to the SAMtools executable

B<--picard_path,-P>
	The path to the Picard JAR file

B<--rm_duplicates,-d>
	If flag is present, will remove duplications of sequences

B<--rm_low_complexity,-c>
	If flag is present, will remove sequences deemed to have low complexity

B<--lc_method,-m>
    Method used to filter out low complexity seqs.  Default is 'dust'

B<--lc_threshold,-t>
    Threshold for filtering low complexity seqs.  Default is 7.

B<--log,-l>
    Logfile.

B<--debug>
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
my $DEFAULT_LC_THRESHOLD = 7;
my $DEFAULT_LC_METHOD = "dust";
my $JAVA_PATH = "/usr/bin/java";
####################################################

my %options;

my $lc_method;
my $lc_threshold;

my $samtools;
my $prinseq;
my $picard;

my $input_file;
my $input_base;
my $out_dir;
my $tmp_dir;

# Allow program to run as module for unit testing if necessary
main() unless caller();

sub main {
	my $results = GetOptions (\%options,
                         "input_file|i=s",
                         "input_base|I=s",
                         "output_dir|o=s",
                         "tmp_dir|t=s",
						 "prinseq_path|p=s",
						 "samtools_path|s=s",
						 "picard_path|p=s",
						 "rm_duplicates|d",
						 "rm_low_complexity|c",
                         "lc_method|m=s",
                         "lc_threshold|t=i",
                         "log|l=s",
                         "debug=s",
                         "help|h"
                          );

    &check_options(\%options);

    $lc_method = (defined $options{'lc_method'}) ? $options{'lc_method'} : $DEFAULT_LC_METHOD;
    $lc_threshold = (defined $options{'lc_threshold'}) ? $options{'lc_threshold'} : $DEFAULT_LC_THRESHOLD;
    $samtools = $options{'samtools_bin'};
    $picard = $options{'picard_path'};
    $prinseq = $options{'prinseq_path'};

    $input_file = $options{'input_file'};
    $input_base = $options{'input_base'};
    $out_dir = $options{'output_dir'}
    $tmp_dir = (defined $options{'tmp_dir'}) ? $options{'tmp_dir'} : $out_dir;

    my $dups_outfile = rm_duplicates($input_file, $samtools, $picard, $tmp_dir) if ($options{'rm_duplicates'};
    my ($lc_1_outfile, $lc_2_outfile) = remove_low_complexity($input_base, $prinseq, $lc_method, $lc_threshold, $tmp_dir) if ($options{'rm_low_complexity');
    my $bad_ids_outfile = merge_bad_ids($dups_outfile, $lc_1_outfile, $lc_2_outfile);
    $filtered = filter_bam_by_ids(
        {   input_bam      => $original_bam,
            output_dir     => $output_dir,
            header_comment => "\@CO\tID:PrinSeq-filtered\tPG:LGTseek\tVN:$VERSION",
            bad_list       => $bad_ids_outfile
        }
    );
}

# Use concatenated fastq files in order to remove duplicates
sub rm_duplicates {
    my ($bam_file, $samtools, $picard, $tmp_dir)
    my $sorted_bam_file;

    #TODO: Currently accepts single BAM.  Needs to accept combined fastq from sam2fasta component

    my ( $name, $path, $suff ) = fileparse( $bam_file, qr/.bam/ );
    &_log($DEBUG,"========= Deduplication Filtering =========");
    # Sort BAM by read names if the header is valid.
     if ( &_run_cmd("$samtools view -H $bam_file | head -n 1") !~ /queryname/ ) {
         &_run_cmd("$samtools sort -n $bam_file $tmp_dir/$name");
         $sorted_bam_file = "$tmp_dir/$name\.bam";
     }
     # Use Picard-tools to mark dups, samtools to view new file, and filter dups into separate file.
     &_run_cmd("$JAVA_PATH -jar $picard FixMateInformation I=$sorted_bam_file TMP_DIR=$tmp_dir SO=coordinate ASSUME_SORTED=1 VALIDATION_STRINGENCY=SILENT");
     &_run_cmd("$JAVA_PATH -jar $picard MarkDuplicates I=$sorted_bam_file TMP_DIR=$tmp_dir OUTPUT=$tmp_dir/$name\_dedup.bam METRICS_FILE=$tmp_dir/$name\_dedup-metrics.txt REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT");
     open( my $BAM, "$samtools view $tmp_dir/$name\_dedup.bam |" ) or &_log($ERROR, "ERROR : Unable to open the deduped bam: $tmp_dir/$name\_dedup.bam");
     my $dups_outfile = $tmp_dir . "/$name\_derep_bad_ids.out"
     open( my $BADS, "> $dups_outfile" ) or &_log($ERROR, "ERROR : Unable to open output file for dedup bad ids: $tmp_dir/$name\_derep_bad_ids.out");
     while ( my $bam_line = &_read_bam_line($BAM) ) {
         if ( $bam_line->{flag}->{pcrdup} ) { print $BADS "$bam_line->{id}\n"; }
     }
     close $BAM;
     close $BADS;

     return $dups_outfile;
}

# Use single-read fastq files for low_complexity filtering
sub rm_low_complexity {
    ($input_base, $prinseq, $lc_method, $lc_threshold, $tmp_dir);
    my ( $name, $path, $suff ) = fileparse( $input_base );  # $suff will be blank
    &_log($DEBUG,"========= Low Complexity Filter ========");
    # Run prinseq for low complexity filtering
    &_run_cmd("perl $prinseq --fastq=$tmp_dir/$name\_1.fastq --out_good null --out_bad=$tmp_dir/$name\_lc_1_bad -lc_method $lc_method -lc_threshold $lc_threshold");

    my $lc_1_outfile = $tmp_dir . "/$name\_lc_1_bad_ids.out";
    if ( -e "$tmp_dir/$name\_lc_1_bad.fastq" ) {
        # Pull out bad ids
        &_log($DEBUG,"========= Pull Low-Cmplx-1 Bad ID's ========");
        &_run_cmd("perl -e 'while(<>){s/\\@//;s/\\/\\d//;print;<>;<>;<>;}' $tmp_dir/$name\_lc_1_bad.fastq > $tmp_dir/$name\_lc_1_bad_ids.out");
    }
    else {
        &_log($DEBUG,"Didn't find any low complexity sequences in read 1");
        &_run_cmd("touch $lc_1_outfile");
    }

    # Run prinseq for low complexity filtering
    &_run_cmd("perl $prinseq_bin --fastq=$tmp_dir/$name\_2.fastq --out_good null --out_bad=$tmp_dir/$name\_lc_2_bad -lc_method $lc_method -lc_threshold $lc_threshold");

    my $lc_2_outfile = $tmp_dir . "/$name\_lc_1_bad_ids.out";
    if ( -e "$tmp_dir/$name\_lc_2_bad.fastq" ) {
        # Pull out bad ids
        &_log($DEBUG,"========= Pull Low-Cmplx-2 Bad ID's ========");
        &_run_cmd("perl -e 'while(<>){s/\\@//;s/\\/\\d//;print;<>;<>;<>;}' $tmp_dir/$name\_lc_2_bad.fastq > $tmp_dir/$name\_lc_2_bad_ids.out");
    }
    else {
        &_log($DEBUG,"Didn't find any low complexity sequences in read 2");
        &_run_cmd("touch $lc_2_outfile");
    }

    return ($lc_1_outfile, $lc_2_outfile);
}

# Merge IDs obtained from filtering dups and low_complexity
sub merge_bad_ids {
    my ($dups_outfile, $lc_1_outfile, $lc_2_outfile) = @_;
    $cmd = "cat";
    if ( $dedup == 1 )        { $cmd = $cmd . " $dups_outfile"; }
    if ( $rm_low_cmplx == 1 ) { $cmd = $cmd . " $lc_1_outfile $lc_2_outfile"; }
    my $bad_ids_outfile = $tmp_dir . "/$name\_prinseq-bad_ids.out";
    $cmd = $cmd . " | sort -u > $bad_ids_outfile";
    &_run_cmd($cmd);
    return $bad_ids_outfile;
}

# Lastly, filter out identified duplicates and low-complexity seqs
# Input Hash keys:
#   input_bam      : the input bam File
#   output_dir     : the output directory
#   header_comment : Additional header information for the SAM/BAM header
#   bad_list       : File of ids to be filtered out due to duplication or low_complexity

sub filter_bam_by_ids {
    my ( $config ) = @_;
     &_log($DEBUG,"======== &filter_bam_by_ids: Start ========");
     if ( !$config->{input_bam} ) {
         &_log($ERROR,"ERROR : Must pass &filter_bam_by_ids an input_bam => <BAM_TO_FILTER>");
     }
     # Check for a non_empty file
     if ( ! -s $config->{input_bam} ) {
         &_log($WARN, "WARNING : Cannot filter ids from an empty input bam: $config->{input_bam})";
         return { count => 0, file => $config->{input_bam} };
     }
     ## Setup hash of ids
     my $bad_ids  = {};
     my %found_ids;
     if ( !$config->{bad_list} ) {
         &_log($ERROR, "ERROR : Must pass &filter_bam_by_ids a file with a list of reads to filter on. Use bad_list => /path/to/list");
     } elsif ( $config->{bad_list} )  {
         $bad_ids  = &_read_ids( $config->{bad_list} );
     }
     ## Setup input and output bams
     my $input = $config->{input_bam};
     my ( $fn, $path, $suf ) = fileparse( $input, ".bam" );
     my $out_dir = defined $config->{output_dir}    ? $config->{output_dir}    : $path;
     my $prefix  = $fn;
     my $suffix  = "filtered";
     my $out     = "$out_dir/$prefix\_$suffix\.bam";
     my $header  = &_run_cmd("samtools view -H $input");
     open( my $in, "-|", "samtools view $input" )
         or &_log($ERROR, "ERROR : SUB-filter_bam_by_ids : Can't open input bam: $input because: $!");
     open( my $fh, "| samtools view -S - -bo $out" )
         or &_log($ERROR"ERROR : SUB-filter_bam_by_ids : Can't open output bam: $out because: $!");
     print $fh "$header";

     if ( defined $config->{header_comment} ) {
         my @split_header = split( /\n/, $header );
         my $last_pg      = $split_header[-1];
         my $last_pg_id   = ( split /\t/, $last_pg )[1];
         $last_pg_id =~ s/ID\:/PP\:/;
         chomp( $config->{header_comment} );
         print $fh "$config->{header_comment}\t$last_pg_id\n";
     }
     # Read through BAM file and print out those alignments not in the bad IDs list
     while (<$in>) {
         chomp;
         my @fields = split(/\t/);
         if ( $config->{bad_list}  && !$bad_ids->{ $fields[0] } ) {
             print $fh "$_\n";
             $found_ids{ $fields[0] }++;
         }
     }
     close $in;
     close $fh;
     my $count = 0;
     $count++ foreach ( keys %found_ids );
     &_log($DEBUG, "======== &filter_bam_by_ids: Finished ========");
     return {
         count => $count,
         bam   => $out
     };
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

   foreach my $req ( qw(output_dir prinseq_path samtools_path picard_path) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }

   if (! definted $opts->{'input_file'} && ! defined ($opts->{'input_base'} ) {
       &_log($ERROR, "Either an input file or input basename must be supplied");
   }
}

# Private subroutine to read the BAM file line by line
sub _read_bam_line {
    my $fh = shift;
    my $line = <$fh>;
    if ( !$line ) { last; }
    chomp($line);
    my ( $id, $flag, $cigar, $sequence ) = ( split /\t/, $line )[ 0, 1, 5, 9 ];
    $id =~ s/(.+)\/\d+/$1/;
    my $converted_flag = &_parseFlag($flag);
    return {
        id       => $id,
        flag     => $converted_flag,
        sequence => $sequence,
        cigar    => $cigar,
        line     => $line,
    };
}

# Private subroutine to make hash of filtered out (bad) ids
sub _read_ids {
    my $list = shift;
    if ( !$list ) { &_log($ERROR, "ERROR : SUB-_read_ids : Must pass a list"); }
    if ( ! -s $list ) {
        &_log($WARNING, "WARNING :  SUB-_read_ids : Input: $list is empty");
    }
    my %hash;
    open IN, "<$list" or &_log($ERROR, "ERROR : SUB-_read_ids unable to open: $list because: $!");
    while (<IN>) {
        chomp;
        $hash{$_} = 1;
    }
    close IN;
    return \%hash;
}

sub _parseFlag {
    my $int    = shift;
    my $rawbin = unpack("B32", pack("N", $int));
    $rawbin =~ s/^0+(?=\d)//;    # otherwise you'll get leading zeros
    my $bin = sprintf( "%012d", $rawbin );
    my $final_bin = reverse $bin;
    return {
        'paired'        => substr( $final_bin, 0,  1 ),
        'proper'        => substr( $final_bin, 1,  1 ),
        'qunmapped'     => substr( $final_bin, 2,  1 ),
        'munmapped'     => substr( $final_bin, 3,  1 ),
        'qrev'          => substr( $final_bin, 4,  1 ),
        'mrev'          => substr( $final_bin, 5,  1 ),
        'first'         => substr( $final_bin, 6,  1 ),
        'last'          => substr( $final_bin, 7,  1 ),
        'secondary'     => substr( $final_bin, 8,  1 ),
        'failqual'      => substr( $final_bin, 9,  1 ),
        'pcrdup'        => substr( $final_bin, 10, 1 ),
        'supplementary' => substr( $final_bin, 11, 1 ),
    };
}

sub _run_cmd {
    my ( $self, $cmd ) = @_;

    &_log($DEBUG,"CMD: $cmd");
    my $res = `$cmd`;
    if ($?) {
        &_log($ERROR,"FAIL_CMD: $cmd died with message: $res");
    }
    return $res;
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
