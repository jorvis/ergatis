#!/usr/bin/perl

=head1  NAME 

lgt_bwa.pl

=head1 SYNOPSIS


      
=head1 OPTIONS

=over 8

This help message

=back

=head1   DESCRIPTION


=head1 INPUT



=head1 OUTPUT


=head1 CONTACT

    David Riley
    driley@som.umaryland.edu

=cut
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
              'input_dir|i=s',
              'input_base|ib=s',
              'ref_file|r=s',
              'ref_file_list|rl=s',
              'output_dir|o=s',
              'mismatch|mm=s',
              'max_gaps|mg=s',
              'max_gap_extentions|mge=s',
              'open_gap_penalty|og=s',
              'extend_gap_penalty|eg=s',
              'threads|t=s',
              'use_bwasw',
              'num_aligns|na=s',
              'bwa_path|b=s',
              'cleanup_sai',
              'help|h');

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

my $ref_files = [];
my $options_string = '';
## make sure all passed options are peachy
&check_parameters(\%options);

# Check if the data are paired or not
my $PAIRED = 0;
$PAIRED=1 if(-e "$options{output_dir}/$options{input_base}_1.fastq");

if($options{use_bwasw}) {
    die "bwasw is not implemented yet\n";
}
else {
    &run_bwa();
}


sub run_bwa {

    foreach my $ref (@$ref_files) {
        chomp $ref;
        $ref =~ /.*\/([^\/]+)\.[^\.]+$/;
        my $refname = $1;
        # In here if we're paired
        if($PAIRED) {

            my $in1 = "$options{input_dir}/$options{input_base}_1.fastq";
            my $out1 = "$options{output_dir}/$refname\_$options{input_base}_aln_sa1.sai";
            my $in2 = "$options{input_dir}/$options{input_base}_2.fastq";
            my $out2 = "$options{output_dir}/$refname\_$options{input_base}_aln_sa2.sai";
            
            # Run the first one through aln
            my $cmd = "$options{bwa_path} aln $options_string $ref $in1 > $out1";
            print "Running: $cmd\n";
            system($cmd) == 0 or die "Unable to run $cmd\n";
            
            # Run the second one through aln
            $cmd = "$options{bwa_path} aln $options_string $ref $in2 > $out2";
            print "Running: $cmd\n";
            system($cmd) == 0 or die "Unable to run $cmd\n";

            # Run sampe
            $cmd = "$options{bwa_path} sampe -n $options{num_aligns} $ref \"$out1\" \"$out2\" \"$in1\" \"$in2\" > $options{output_dir}/$refname\_$options{input_base}.sam";
            print "Running: $cmd\n";
            system($cmd) == 0 or die "Unable to run $cmd\n";


            if($options{cleanup_sai}) {
                $cmd = "rm -f $out1 $out2";
                print "Running: $cmd\n";
                system($cmd) == 0 or die "Unable to run $cmd\n";
            }
        }
        
        # In here if we aren't paired
        else {

            my $in = "$options{input_dir}/$options{input_base}.fastq";
            my $out = "$options{output_dir}/$refname\_$options{input_base}_aln_sa.sai";
            my $cmd = "$options{bwa_path} aln $options_string $ref $in > $out";

            print "Running: $cmd\n";
            system($cmd) == 0 or die "Unable to run $cmd\n";

            $cmd = "$options{bwa_path} samse -n $options{num_aligns} $ref \"$out\" \"$in\" > $options{output_dir}/$refname\_$options{input_base}.sam";
            print "Running: $cmd\n";
            system($cmd) == 0 or die "Unable to run $cmd\n";

            if($options{cleanup_sai}) {
                $cmd = "rm -f $out";
                print "Running: $cmd\n";
                system($cmd) == 0 or die "Unable to run $cmd\n";
            }
        }
    }

}
sub check_parameters {
    
    ## make sure input file exists
    if (! -e $options{'input_dir'}) { print STDERR "Input invalid\n"; pod2usage( {-exitval=>0, -verbose => 0, -output => \*STDOUT})};

    if($options{ref_file}) {
        my @files = split(/,/,$options{ref_file});
        $ref_files = \@files;
    }
    elsif($options{ref_file_list}) {
        my @lines = `cat $options{ref_file_list}`;
        $ref_files = \@lines;
    }
    else {
        die "No reference file specified in ref_file or ref_file_list"
    }

    my $options_to_param_keys = {
        'mismatch' => '-M',
        'max_gaps' => '-o',
        'max_gap_extensions' => '-e',
        'open_gap_penalty' => '-O',
        'extend_gap_penalty' => '-E',
        'threads' => '-t'
    };

    my $opts = [];
    foreach my $key (keys %$options_to_param_keys) {
        if($options{$key}) {
            push(@$opts, "$options_to_param_keys->{$key} $options{$key}");
        }
    }
    $options_string = join(" ",@$opts);
    
    return 1;
}
