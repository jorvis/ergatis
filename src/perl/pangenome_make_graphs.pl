#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
=head1  NAME

pangenome_do_R.pl - Runs a specified R file on the specified input file.

=head1 SYNOPSIS

  USAGE: pangenome_do_R.pl
    --input=/path/to/pangenome.table.txt
    --r_script=/path/to/some_R_script.R
    --output_path=/path/to/output/
    --title='Genus species'
    [ --log=/path/to/some/log ]

=head1 OPTIONS

B<--input_table,-i>
    The pangenome.table.txt file that comes from pangenome_do_analysis.pl or the pangenome.output
    file that comes from pangenome_make_pangenome.pl

B<--r_script,-r>
    The R script to run on the input pangenome_table

B<--output_path,-o>
    Path to which output files will be written.

B<--title,-g>
    Path to which output files will be written.

B<--help,-h>
    This help message/documentation.

=head1   DESCRIPTION

    The pangenome analysis script creates an array of BLAST results data which is then
    processed to create pangenome data.

=head1 INPUT

    The input should be a list of files containing serialized BLAST results array data.

=head1 OUTPUT

    There is no output unless you use the --log option.

=cut

use Pod::Usage;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = GetOptions( \%options,
                          'input_table|i=s',
                          'r_script|r=s',
			  'r_exec_path|p=s',
                          'output_path|o=s',
                          'title|g=s',
                          'help|h') || pod2usage();



pod2usage if $options{'help'};


my $title = $options{'title'};
my $pangenome_table = $options{'input_table'};
my $r_script = $options{'r_script'};
my $output_path = $options{'output_path'};
my $R_EXEC_PATH = $options{'r_exec_path'};

my $r_filename = $r_script;
$r_filename =~  s/^.*\/([^\/]*)/$1/;

open (IN, "$r_script") || die "Couldn't open R script '$r_script': $!";

my $input_r = "$output_path/$r_filename"."in";

my $ps_file = $r_filename;
$ps_file =~ s/\.R/\.ps/;
my $eps_file= $ps_file;
$eps_file =~ s/\.ps/\.eps/;

open (OUT, ">$input_r");
while (<IN>) {
    s/###TITLE###/$title/;
    s/###input_file###/$pangenome_table/;
    s/###output_path###/$output_path\//;
    print OUT;
}
close OUT;
close IN;

if(system("$R_EXEC_PATH CMD BATCH $input_r $output_path/".$r_filename."out")==0 ) {
    my @r_files = `ls $output_path/*.ps`;
    foreach my $file(@r_files) {
        chomp $file;
        my $eps_file = $file;
        $eps_file =~ s/\.ps/\.eps/;
        print STDERR "converting $file to eps\n";
        system("ps2epsi $file $eps_file" )==0 or warn "Unable to run ps2epsi command $!\n";
    }
}
else {
    warn "Unable to run the R command $!\n";
}

print STDERR "done.\n";
exit();
