#!/usr/local/bin/perl -w

=head1  NAME 

load_multi_aat_data.pl - load all AAT components within a single pipeline.

=head1 SYNOPSIS

USAGE: load_multi_aat_data.pl 
        --pipeline_id=453
        --repository_root=/usr/local/annotation/AA1
        --database=aa1

=head1 OPTIONS

B<--input_list,-i> 
    htab list file from an hmmpfam workflow component run.

B<--repository_root,-r> 
    The project directory, just under which we should find the output_repository directory.

B<--database,-d> 
    The name of the Sybase database to load.

B<--help,-h> 
    This help message/documentation.

=head1   DESCRIPTION

    This script will scan the output directories for any aat output from the given
    pipeline.  Please note that it will not currently load data for coordinate
    adjusted output, which we have been trying to avoid creating.

=head1 CONTACT

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use File::Basename;
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
			  'pipeline_id|p=s',
              'repository_root|r=s',
              'database|d=s',
			  'help|h') || pod2usage();

my $loader = '/usr/local/devel/ANNOTATION/ard/current/bin/bsml2legacydb';

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

## play nicely
umask(0000);

## make sure all passed options are peachy
&check_parameters(\%options);

my @output_lists;

## check for any aat_aa runs under this pipeline:
for ( glob "$options{repository_root}/output_repository/aat_aa/$options{pipeline_id}*" ) {
    print STDERR "adding $_/aat_aa.bsml.list to load list\n";
    push @output_lists, "$_/aat_aa.bsml.list";
}

## check for any aat_na runs under this pipeline:
for ( glob "$options{repository_root}/output_repository/aat_na/$options{pipeline_id}*" ) {
    print STDERR "adding $_/aat_na.bsml.list to load list\n";
    push @output_lists, "$_/aat_na.bsml.list";
}

foreach my $list ( @output_lists ) {
    my $output_dir = dirname($list);

    my $cmd = "$loader --input_list=$list --database=$options{database} --debug=4 --log=$output_dir/loader.log";
    print STDERR "currently executing: $cmd\n\n";
    `$cmd`;
}



exit(0);


sub check_parameters {
    
    ## pipeline_id and repository_root are required
    unless ( defined $options{pipeline_id} && $options{repository_root} ) {
        print STDERR "--pipeline_id and --repository_root options are required\n\n";
        pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
    }
    
    ## database is required.  force it to lower case.
    if ( defined $options{database} ) {
        $options{database} = lc $options{database};
    } else {
        print STDERR "\n\n--database is a required option\n\n";
        exit(1);
    }
    
    ## make sure the repository root has an output_repository directory
    unless ( -d "$options{repository_root}/output_repository" ) {
        print STDERR "\n\nthe repository root passed doesn't contain an output_repository directory.\n\n";
        exit(1);
    }
    
    ## make sure the loader exists
    unless ( -e $loader ) {
        print STDERR "\n\nloader $loader does not exist.  please specify another within the code\n\n";
        exit(1);
    }
}






