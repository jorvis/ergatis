#!/usr/bin/perl

$| = 1;

use warnings;
use strict;

use File::Basename;
use File::Path;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my ($output, $ptt_file_list, $nuc_file_list, $prot_file_list, $help);

    GetOptions('output_iter_list|o=s' => \$output,
               'ptt_file_list|p=s'    => \$ptt_file_list,
               'nuc_file_list|n=s'    => \$nuc_file_list,
               'prot_file_list|a=s'   => \$prot_file_list,
               'help|man|h|m|?'       => \$help) || die "\n\nProblems processing the options\n\n";

if (!-e $ptt_file_list) {
    die "Specified input list does not exist";
}
if (!-e $nuc_file_list) {
    die "Specified input list does not exist";
}
if (!-e $prot_file_list) {
    die "Specified input list does not exist";
}
unless ($output) {
    die "Output file name must be provided with --output flag";
}

my $seq_id_to_files = {};

open (my $outfh, ">$output") || die "Couldn't open file '$output' for writing: $!";

&process_list_file($ptt_file_list,'ptt');
&process_list_file($nuc_file_list,'nuc');
&process_list_file($prot_file_list,'prot');

print $outfh "\$;GROUP_COUNT\$;\t\$;PTT_FILE_BASE\$;\t\$;PTT_FILE_NAME\$;\t\$;PTT_FILE_PATH\$;\t\$;NUC_FILE_BASE\$;\t\$;NUC_FILE_NAME\$;\t\$;NUC_FILE_PATH\$;\t\$;PROT_FILE_BASE\$;\t\$;PROT_FILE_NAME\$;\t\$;PROT_FILE_PATH\$;\n";
my $count =0;
my $types = ['ptt','nuc','prot'];
foreach my $seq_id (keys %$seq_id_to_files) {
    $count++;
    my $files = $seq_id_to_files->{$seq_id};
    my @fields;
    foreach my $type (@$types) {
        if(!(scalar keys %{$files->{$type}} ==3)) {
            die "Problem mapping $seq_id\n";
        }
        push(@fields, ($files->{$type}->{'base'},$files->{$type}->{'name'},$files->{$type}->{'path'}));
    }
    print $outfh join("\t",($count,@fields))."\n";
}

sub process_list_file {
    my $file = shift;
    my $type = shift;
    open (my $fh, $file) || die "Couldn't open input list '$file' for reading: $!";

    foreach my $line (<$fh>) {
        chomp $line;
        my $seq_id = basename($line);
        $seq_id =~ s/\.[^\.]+$//; # HACK here since we don't know the extensions.
        $seq_id_to_files->{$seq_id}->{$type} = {'path'=> $line,
                                                'base' => basename($line),
                                                'name' => $seq_id};
    }
    close $fh;
}



#my $counter = 0;
#while (scalar @input_sequences > 1) {
#    my ($query, @subjects) = @input_sequences;
#    foreach my $subject (@subjects) {
#        $counter++;
#        my $i_file_name1 = basename($query);
#        my $i_file_name2 = basename($subject);
#        $i_file_name1 =~ /(.*)\.([^\.]+)$/;
#        my $i_file_base1 = $1;
#        my $i_file_ext1 = $2;
#        $i_file_name2 =~ /(.*)\.([^\.]+)$/;
#        my $i_file_base2 = $1;
#        my $i_file_ext2 = $2;
#        my $i_dir1 = dirname($i_file_name1);
#        my $i_dir2 = dirname($i_file_name2);
#        
#        print $outfh "pair_$counter\t$i_file_base1\t$i_file_name1\t$query\t$i_file_ext1\t$i_dir1\t$i_file_base2\t$i_file_name2\t$subject\t$i_file_ext2\t$i_dir2\n";
#    }
#    shift @input_sequences;
#}
