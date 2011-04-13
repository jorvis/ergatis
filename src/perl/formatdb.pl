#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use warnings;
use Getopt::Long;

my $options = {};

GetOptions( $options, 'input_file_list|l:s', 'input_files|f:s', 'formatdb_path|p=s', 'protein=s', 'database_name|n=s' );

my $files = "";

if( $$options{'input_files'} ) {
	$files .= $$options{'input_files'}." ";
} 

if( $$options{'input_file_list'} ) {
	$files .= get_files( $$options{'input_file_list'} );
}

my $command = "$$options{'formatdb_path'} -p $$options{'protein'} -oT -i '$files' -n $$options{'database_name'}";

print $command,"\n";
system( $command ) == 0 or die "Error in executing the command, $command, $!\n";

exit $?;


sub get_files {
	my ($list) = @_;
	open( FH, "<$list" ) or die "Error in opening the file, $list, $!\n";
	my @files = ();
	while( my $file = <FH> ) {
		chomp $file;
		if( -e $file ) { push @files, $file; }
		else { print STDERR "$file : No such file exists\n"; }
	}
	return join($", @files);
}
