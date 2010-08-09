#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config pass_through);
use FindBin qw($Bin);
use File::Basename;
use Cwd qw(abs_path);

my $cufflinks_home = "/usr/local/packages/cufflinks";
#my $ergatis_base = "/usr/local/projects/ergatis/package-devel";
my $ergatis_base = abs_path(dirname($Bin));

my $sam_list;

&GetOptions("sam-list=s"         => \$sam_list,
            );

fix_environment();

my $sam_file = get_sam_file($sam_list);

my $cufflinks_bin = "$cufflinks_home/cufflinks";
my @cufflinks_args = (@ARGV, $sam_file);

print "$cufflinks_bin " . @cufflinks_args . "\n";

# Invoke the Tophat program and pass it the arguments that it requires for
# operation.
exec ( $cufflinks_bin , @cufflinks_args);

sub get_sam_file {
    my $list_file = shift;

    if (! (-f $list_file && -r $list_file)) {
        warn "The sam list file, $list_file, does not appear to exist or is unreadable.\n";
        exit 2;
    }

    open (FILE, "<", $list_file) or die "Unable to open the list file $list_file: $!";
    my @lines = (<FILE>);
    close FILE or die "Unable to close filehandle to $list_file: $!";

    my $sam_file;
    if (scalar @lines) {
        $sam_file = $lines[0];
        chomp($sam_file);
    }

    if ( ! -f $sam_file) {
        warn "File $sam_file does not appear to exist.";
        exit 2;
    }
    return $sam_file;
}

# THis subroutine mimics the fixups evident in other Ergatis scripts and
# wrappers. Most importantly, it adjusts PATH environment variable to include
# the paths to both bowtie and tophat. WHen Tophat runs, it needs the bowtie
# binaries to be available in the PATH for invocation.

sub fix_environment {
    umask 0000;

    delete $ENV{"PERL5LIB"};
    delete $ENV{"LD_LIBRARY_PATH"};

    if ( exists $ENV{"SCHEMA_DOCS_DIR"} ) {
        $ENV{"SCHEMA_DOCS_DIR"} = "";
    }
    if ( exists $ENV{"WORKFLOW_WRAPPERS_DIR"} ) {
        $ENV{"WORKFLOW_WRAPPERS_DIR"} = "$ergatis_base/bin";
    }
    if ( exists $ENV{"WORKFLOW_DOCS_DIR"} ) {
        $ENV{"WORKFLOW_DOCS_DIR"} = "";
    }

    $ENV{"LANG"} = "C";
    $ENV{"LC_ALL"} = "C";
    $ENV{"PERL_MOD_DIR"} = "$ergatis_base/lib/5.8.8";
    $ENV{"PERL5LIB"} = "$ergatis_base/lib/perl5/";
}
