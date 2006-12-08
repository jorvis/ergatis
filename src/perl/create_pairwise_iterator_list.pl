#!/usr/local/bin/perl

$| = 1;

use warnings;
use strict;

use File::Basename;
use File::Path;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my ($output, $input_file_list, $input_file, $input_directory, $input_directory_extension, $help);

    GetOptions('output_iter_list|o=s'             => \$output,
               'input_file_list|i=s'    => \$input_file_list,
               'input_file|f=s'         => \$input_file,
               'input_directory|d=s'    => \$input_directory,
               'input_directory_extension|x=s', => \$input_directory_extension,
               'help|man|h|m|?'         => \$help) || die "\n\nProblems processing the options\n\n";

if ($input_file_list && !-e $input_file_list) {
    die "Specified input list does not exist";
}

if ($input_file || $input_directory) {
    die "Create pairwise iterator list currently doesn't support input_file or input_directory";
}

unless ($output) {
    die "Output file name must be provided with --output flag";
}

open (my $outfh, ">$output") || die "Couldn't open file '$output' for writing: $!";

open (my $list, $input_file_list) || die "Couldn't open input list '$input_file_list' for reading: $!";

my @input_sequences;

my $dieflag;
while (<$list>) {
    chomp;
    if (!-e $_) {
        print STDERR "Specified input sequence '$_' does not exist\n";
        $dieflag = 1;
    } else {
        push (@input_sequences, $_);
    }
}
close $list;

if ($dieflag) {
    die "Some input files are missing";
}

print $outfh "\$;PAIR_COUNT\$;\t\$;I_FILE_BASE1\$;\t\$;I_FILE_NAME1\$;\t\$;I_FILE_PATH1\$;\t\$;I_FILE_EXT1\$;\t\$;I_DIR1\$;\t\$;I_FILE_BASE2\$;\t\$;I_FILE_NAME2\$;\t\$;I_FILE_PATH2\$;\t\$;I_FILE_EXT2\$;\t\$;I_DIR2\$;\n";
my $counter = 0;
while (scalar @input_sequences > 1) {
    my ($query, @subjects) = @input_sequences;
    foreach my $subject (@subjects) {
        $counter++;
        my $i_file_name1 = basename($query);
        my $i_file_name2 = basename($subject);
        $i_file_name1 =~ /(.*)\.([^\.]+)$/;
        my $i_file_base1 = $1;
        my $i_file_ext1 = $2;
        $i_file_name2 =~ /(.*)\.([^\.]+)$/;
        my $i_file_base2 = $1;
        my $i_file_ext2 = $2;
        my $i_dir1 = dirname($i_file_name1);
        my $i_dir2 = dirname($i_file_name2);
        
        print $outfh "pair_$counter\t$i_file_base1\t$i_file_name1\t$query\t$i_file_ext1\t$i_dir1\t$i_file_base2\t$i_file_name2\t$subject\t$i_file_ext2\t$i_dir2\n";
    }
    shift @input_sequences;
}
