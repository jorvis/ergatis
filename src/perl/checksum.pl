#! /local/perl/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME

checksum.pl

=head1 SYNOPSIS

Checksum.pl is a utility to generate md5 checksums for a set of input files. Input is either from a file containing a list of filenames (one per line) or STDIN (one file name per line). Output consists of a file containing each file name and md5 checksum, one per line. 

USAGE: checksum.pl -f checksum_filelist.txt -o CHECKSUM.file

=head1 OPTIONS
=over 4

B<--filelist, -f> [OPTIONAL]   input file containing a list of filenames, one per line

B<--outfile, -o> [REQUIRED]  output file for md5 checksums

B<--help, -h> [OPTIONAL] program help

=cut

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

my %options = ();

my $results = GetOptions( \%options, 'filelist|f=s', 'outfile|o=s', 'help|h' );

if( $options{'help'} )
{
    pod2usage({-exitval=>0, -verbose => 2, -output => \*STDOUT});
}

if( !($options{'outfile'}) )
{
    print STDERR "ERROR: User must specify an output file...\n";
    exit(1);
}

open( OUTFILE, ">$options{'outfile'}" ) or die "Could not open $options{'outfile'}\n";

if( $options{'filelist'} )
{
    open( INFILE, $options{'filelist'} ) or die "Could not open $options{'filelist'}\n";

    while( my $line = <INFILE> )
    {
	chomp $line;

	my $checksum = `/usr/local/bin/md5 $line`;
	
	print OUTFILE "$line\t$checksum";
    }

    close( INFILE );
}
else
{
    while( my $line = <STDIN> )
    {
	chomp $line;

	my $checksum = `/usr/local/bin/md5 $line`;
	
	print OUTFILE "$line\t$checksum";
    }
}

close( OUTFILE );

