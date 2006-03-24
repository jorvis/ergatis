#!/usr/bin/perl

use strict;
#use warnings;
use File::Basename;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options=();
GetOptions(\%options,
        "ini_template|i=s",
		"ergatis_install_root|e=s",
		"man|m",
        "help|h",
);

## We'll make sure that the ergatis install path ends with a /
unless ($options{'ergatis_install_root'} =~ /\/$/) {
	$options{'ergatis_install_root'} .= "/";
}

my $ini_file = basename($options{'ini_template'});

unless(-e $options{'ini_template'} && $ini_file =~ /^([^\.]+)/) {
			die $options{'ini_template'}." is not a valid ini template file.\n";
}
my $component = $1;
my $component_ini = $component."conf.ini";

unless(-e $options{'ergatis_install_root'}."docs/$component_ini") {
	die $options{'ergatis_install_root'}." is not a valid ergatis install directory.\n";
}

my $component_ini_path = $options{'ergatis_install_root'}."docs/".$component_ini;

my @template_variables;
my @template_variables_lines;
my %template_variables;
my $template_revision;
my @component_conf_variables;
my @component_conf_variables_lines;
my %component_conf_variables;
my $component_conf_revision;
my $count;
my $err_flag = 0;

print STDERR "# Processing $ini_file / $component_ini\n";
open (IN, $options{'ini_template'}) || die "Couldn't open $ini_file for reading";

print STDERR "# Parsing $ini_file\n";
$count = 0;
while (<IN>) {
	$count++;
	chomp;
	if (/(\$;[A-Z0-9_]+\$[^;])/ || /(\$;[A-Z0-9_]+[^\$];)/ || /([^\$];[A-Z0-9_]+\$;)/ || /(\$[^;]?[A-Z0-9_]+\$;)/ || /(\$[A-Z0-9_]+\$)/ || /(;\$[A-Z0-9_]+;\$)/) {
		print STDERR "## Malformed variable $1 at line $count\n";
		$err_flag = 1;
	}
	#while (/\$;([^;]+)\$;/g) {
	while (/\$;([^;]+)\$;\s*=/g) {
		$template_variables{$1} = 1;
		push(@template_variables, $1);
		push(@template_variables_lines, $count." | ".$_);
		if ($template_variables[$#template_variables] eq 'REVISION') {
			/Revision:([^\$]+)/i;
			$template_revision = $1;
			$template_revision =~ s/\s+//g;
		}
	}
}
close IN;

open (IN, $component_ini_path) || die "Couldn't open $component_ini for reading";

print STDERR "# Parsing $component_ini\n";
$count = 0;
while (<IN>) {
	$count++;
	chomp;
	if (/(\$;[A-Z0-9_]+\$[^;])/ || /(\$;[A-Z0-9_]+[^\$];)/ || /([^\$];[A-Z0-9_]+\$;)/ || /(\$[^;]?[A-Z0-9_]+\$;)/ || /(\$[A-Z0-9_]+\$)/ || /(;\$[A-Z0-9_]+;\$)/) {
		print STDERR "## Malformed variable $1 at line $count\n";
		$err_flag = 1;
	}
	#while (/\$;([^;]+)\$;/g) {
	while (/\$;([^;]+)\$;\s*=/g) {
		$component_conf_variables{$1} = 1;
		push(@component_conf_variables, $1);
		push(@component_conf_variables_lines, $count." | ".$_);
		if ($component_conf_variables[$#component_conf_variables] eq 'REVISION') {
			/Revision:([^\$]+)/i;
			$component_conf_revision = $1;
			$component_conf_revision =~ s/\s+//g;
		}
	}
}
close IN;

unless ($template_revision == $component_conf_revision) {
	print STDERR "## Revision numbers of template and component conf do not match!\n";
	if ($template_revision eq '') {
		$template_revision = 'N/A';
	}
	if ($component_conf_revision eq '') {
		$component_conf_revision = 'N/A';
	}
	print STDERR "##   Template revision: \t$template_revision\n";
	print STDERR "##   Component revision:\t$component_conf_revision\n";
	$err_flag = 1;
}

my $head_flag_a = 0;

foreach my $k(keys(%template_variables)) {
	unless (defined($component_conf_variables{$k})) {
		$err_flag = 2;
		unless ($head_flag_a) {
			print STDERR "##  Variables unique to $ini_file\n";
			$head_flag_a = 1;
		}
		print STDERR "    [\$;$k\$;]\n";
		my $var_count = scalar(@template_variables);
		for (my $i = 0; $i < $var_count; $i++) {
			if ($template_variables[$i] eq $k) {
				print STDERR "    ".$template_variables_lines[$i]."\n";
			}
		}
	}
}

my $head_flag_b = 0;
foreach my $k(keys(%component_conf_variables)) {
	unless (defined($template_variables{$k})) {
		$err_flag = 2;
		unless ($head_flag_b) {
			print STDERR "##  Variables unique to $component_ini\n";
			$head_flag_b = 1;
		}
		print STDERR "    [\$;$k\$;]\n";
		my $var_count = scalar(@component_conf_variables);
		for (my $i = 0; $i < $var_count; $i++) {
			if ($component_conf_variables[$i] eq $k) {
				print STDERR "    ".$component_conf_variables_lines[$i]."\n";
			}
		}
	}
}
if ($err_flag > 1) {
	print STDERR "# Variable mismatch in $ini_file / $component_ini\n";
}

#if ($err_flag == 2) {
#	system('sed -i -r ');
#}

exit($err_flag);
