#!/usr/bin/env perl

=head1 NAME

B<replace_software_config_values.pl> Replace software config values in an ergatis software.config with values from a template

=head1 SYNOPSIS

 USAGE: replace_software_config_values.pl
        --template_config=/path/to/template/software.config
        --new_config=/path/to/unconfigured/software.config
      [ --log=/path/to/file.log
        --debug=3
        --help ]

=head1 OPTIONS

B<--template_config,-t>
    A configured software config

B<--new_config,-n>
    An unconfigured software config.

B<--log,-l>
    Logfile.

B<--debug,-d>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

    Will read in values from a template software config and will write the values of matching variables 
    to the new file. Variables which are not found within the new software.config or variables which are 
    not present in the template will be ignored

=head1  CONTACT

    Kevin Galens
    kgalens@gmail.com

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Config::IniFiles;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
####################################################

my %options;
my $results = GetOptions (\%options,
                         "template_config|t=s",
                         "new_config|n=s",
                         "log|l=s",
                         "debug|d=s",
                         "help|h"
                          );

&check_options(\%options);

&replace_software_config_values( $options{'template_config'}, $options{'new_config'} );

sub replace_software_config_values {
  my ($template_file, $new_config_file) = @_;

  print "Replacing $new_config_file with values from $template_file\n";
    
  my $new_config = new Config::IniFiles( -file => $new_config_file ) or
	die("Failed to open new software config: $!");

  my $template_config = new Config::IniFiles( -file => $template_file ) or
	die("Failed to open template software config: $!");

  if (defined $template_config) {
	for my $section ( $new_config->Sections() ) {
	  if ( $template_config->SectionExists($section) ) {		
		for my $parameter ( $new_config->Parameters($section) ) {
		  if ( defined $template_config->val($section, $parameter) ) {
			$new_config->setval( $section, $parameter, $template_config->val($section, $parameter) );
		  }
		}
	  }
	}
  }
    
  $new_config->RewriteConfig;
}

sub check_options {

   my $opts = shift;

   if( $opts->{'help'} ) {
       &_pod;
   }

   if( $opts->{'log'} ) {
       open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
   } else {
	 $logfh = \*STDOUT;
   }

   foreach my $req ( qw(template_config new_config) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }

   $debug = $opts->{'debug'} if( $opts->{'debug'} );
}

sub _log {
   my ($level, $msg) = @_;
   if( $level <= $debug ) {
	 print $logfh "$msg\n" if( defined( $logfh ) );
   }
   exit(1) if( $level == $ERROR );
}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
