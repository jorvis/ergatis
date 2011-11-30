#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Basename;

my $options = {};

GetOptions( $options, 'input_files|f=s', 'formatdb_path|p=s', 'protein=s', 'database_name|n=s', 'output_dir|o=s', 'other_args|a=s' );

my $files = "";
my $database_name;
my $formatdb_args = "";

if( $$options{'input_files'} ) {
	$files = $$options{'input_files'};
} 

#if(exists($$options{'input_file_list'})) {
#	$files .= get_files( $$options{'input_file_list'} );
#}

if(exists($$options{'other_args'})) {
	$formatdb_args = $$options{'other_args'};
}

if (exists($$options{'database_name'})) {
	$database_name = $$options{'output_dir'}."/".$$options{'database_name'};
} else {
	my ($file_base,$file_dir,$file_ext) = fileparse($$options{'input_files'},qr/\.[^.]*/);
	$database_name = $$options{'output_dir'}."/".$file_base;
}
my $log_file = $database_name."_formatdb.log";
my $command = "$$options{'formatdb_path'} -p $$options{'protein'} -o T -i '$files' -n $database_name -l $log_file $formatdb_args";

print STDERR $command,"\n";
system( $command ) == 0 or die "Error in executing the command, $command, $!\n";

exit $?;


#sub get_files {
#	my ($list) = @_;
#	open( FH, "<$list" ) or die "Error in opening the file, $list, $!\n";
#	my @files = ();
#	while( my $file = <FH> ) {
#		chomp $file;
#		if( -e $file ) { push @files, $file; }
#		else { print STDERR "$file : No such file exists\n"; }
#	}
#	return join($", @files);
#}

sub check_parameters {
	my @required = qw(input_files formatdb_path protein output_dir);
	foreach my $option(@required) {
		if(!defined($$options{$option})) {
			die "ERROR! : Required option $option not passed\n";
		}
	}
}
