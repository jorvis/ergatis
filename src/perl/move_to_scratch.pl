#!/usr/bin/perl

=head1 move_to_scratch.pl - Script designed to move unwanted files to the scratch repository

=head1 DESCRIPTION

	This script will take files created by various components, specified by file extension and component, and move them into the user's specified scratch repository.  This is intended to shrink space from Ergatis prokaryotic pipeline runs by not keeping files that wouldn't be useful for error diagnosis nor used by components later in the program.

=head1 INPUT

	--input_list - A tab-delimited file listing the 'component.token' name and the file extension to move

create_pseudomolecules.default	.multi.fasta
create_pseudomolecules.default	.fasta.orig
glimmer3.iter1	.detail
hmmpfam3.pre	.htab

	--scratch_dir - The scratch directory you are using (including the move_to_scratch folder)

	--output_dir - The output repository directory (including the move_to_scratch/pipe##_default/ directories)

=head1 SYNOPSIS

perl ./move_to_scratch.pl -i /path/to/input/list -s /local/scratch/username/move_to_scratch/1234_default -o /repository/path/output_repository/move_to_scratch/1234_default

=head1 OUTPUT

	A Stderr file that lists if there are problems with the moving process

=head1 CONTACT

Shaun Adkins
sadkins@som.umaryland.edu

=cut

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Ergatis::Logger;


my %options = ();
my $results = GetOptions (\%options,
		'output_dir|o=s',
		'input_list|i=s',
		'scratch_dir|s=s',
		'log|l=s',
		'debug|b=s',
                'help') || pod2usage();



## Display documentation
if( $options{'help'} ){
        pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## Getting the log file
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

## Verifying parameters
&check_parameters(\%options);

my $pipeline;

my @paths = parse_list($options{'input_list'}, $options{'output_dir'});
move_files(\@paths, $options{'scratch_dir'});

exit(0);

# Break down list of components/extensions and derive paths of files to save to an array
sub parse_list {
    my $list = shift;
    my $output_d = shift;

    my @path;
    my @files;
    my $repository;

    if ($output_d =~ /(.+)\/output_repository\/move_to_scratch\/(\d+)_(\w+)/) {
	$repository = $1;
	$pipeline = $2;
    }

    open LIST, $list or $logger->logdie("Cannot open the specified input list file $list");
    while (<LIST>) {
	chomp;
	my @parse = split;
	my ($component, $token) = split(/\./, $parse[0]);
	my $path = $repository . '/output_repository/' . "$component/$pipeline" . "_$token/";

	push (@files, `find $path -name '*$parse[1]' -type f`);
    }
    close LIST;

    return @files;
}


# Move files from array to scratch and (optional) create symlinks to those in the the original path
sub move_files {
    my $files = shift;
    my $scratch = shift;

    foreach my $f (@{$files}) {
    	chomp $f;
	my $base = `basename $f`;

	system("mv $f $scratch/$base");

#	system("ln -sf $scratch/$base $f");
    }

    #make symlink to original location (ln -s)
    #of course in 2 weeks the scratch contents are deleted making the symlinks useless
}

sub check_parameters {
    my $options = shift;
    unless ($options{output_dir} && $options{input_list} && $options{scratch_dir}) {
	$logger->logdie("All the manadatory parameters should be passed");
    }

}
