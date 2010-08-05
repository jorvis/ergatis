#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config pass_through);

my $tophat_home="/usr/local/stow/tophat-1.0.12";
my $bowtie_home="/opt/opt-packages/bowtie-0.12.0";
my $ergatis_base = "/opt/ergatis";

my ($bowtie_index_dir, $extension, $directory);

&GetOptions("extension=s"         => \$extension,
            "directory=s"         => \$directory,
            "bowtie_index_dir=s" => \$bowtie_index_dir,
            );

if ( exists $ENV{"SCHEMA_DOCS_DIR"} ) {
    $ENV{"SCHEMA_DOCS_DIR"} = "";
}
if ( exists $ENV{"WORKFLOW_WRAPPERS_DIR"} ) {
    $ENV{"WORKFLOW_WRAPPERS_DIR"} = "$ergatis_base/bin";
}
if ( exists $ENV{"WORKFLOW_DOCS_DIR"} ) {
    $ENV{"WORKFLOW_DOCS_DIR"} = "";
}

fix_environment();
my $csv = build_csv_filelist($directory, $extension);

my $tophat_bin = "$tophat_home/bin/tophat";
my @tophat_args = (@ARGV, $bowtie_index_dir, $csv);
#print $tophat_bin . " " . join(" ", @tophat_args) . "\n";

# Invoke the Tophat program and pass it the arguments that it requires for
# operation.
exec ( $tophat_bin , @tophat_args);


# THis subroutine mimics the fixups evident in other Ergatis scripts and
# wrappers. Most importantly, it adjusts PATH environment variable to include
# the paths to both bowtie and tophat. WHen Tophat runs, it needs the bowtie
# binaries to be available in the PATH for invocation.

sub fix_environment {
    umask 0000;

    delete $ENV{"PERL5LIB"};
    delete $ENV{"LD_LIBRARY_PATH"};

    $ENV{"LANG"} = "C";
    $ENV{"LC_ALL"} = "C";
    $ENV{"PERL_MOD_DIR"} = "$ergatis_base/lib/5.8.8";
    $ENV{"PERL5LIB"} = "$ergatis_base/lib/perl5/";
    $ENV{"BOWTIE_HOME"} = $bowtie_home;
    $ENV{"PATH"} = $ENV{"PATH"} . ":$bowtie_home:$tophat_home/bin";
}

sub build_csv_filelist {
    my ($directory, $extension) = @_;
    my @files;
   opendir(DIR, $directory) or die "Cannot open directory $directory: $!";
   my @files = grep {-f "$directory/$_" } readdir(DIR);
   closedir DIR;

   @files = grep { /\.$extension$/ } @files;
   @files = map { $directory . "/" . $_ } @files;
   my $csv = join(",", @files);
   return $csv;
}
