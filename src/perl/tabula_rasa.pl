#!/usr/local/packages/perl-5.8.5/bin/perl
                                                                               
eval 'exec /usr/local/packages/perl-5.8.5/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %options = ();
my $results = GetOptions( \%options,
			  'input_command=s',
			  'iter_file_path=s',
			  'iter_file_base=s',
			  'iter_file_ext=s',
			  'iter_file_name=s',
			  'iter_dir=s',
			  'output_directory=s'
			  ) || die('Unsupported option');

my $cmd = $options{input_command} || die '$;INPUT_COMMAND$; unspecified';

# for reasons unclear to me, key needs to be
# VAR not $;VAR$; (and handle $;'s in regexp)
my %replacements = (
		    'ITER_FILE_PATH' => $options{iter_file_path},
		    'ITER_FILE_BASE' => $options{iter_file_base},
		    'ITER_FILE_EXT'  => $options{iter_file_ext},
		    'ITER_FILE_NAME' => $options{iter_file_name},
		    'ITER_DIR'       => $options{iter_dir},
		    'OUTPUT_DIRECTORY'     => $options{output_directory}
		    );

# print "Original cmd: $cmd\n";
foreach (keys %replacements) {
#    print "$_ => $replacements{$_}\n";
    $cmd =~ s/\$;$_\$;/$replacements{$_}/g;
}
# print "Revised cmd:  $cmd\n";

die "Unresolved variables in command ($cmd)" if ($cmd =~ /\$;/);

exec $cmd;
