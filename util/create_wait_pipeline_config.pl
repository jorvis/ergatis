#!/usr/bin/env perl
package main;    # Change this to reflect name of script.  If you want to use as module, make sure to save as .pm instead of .pl

=head1 NAME

create_wait_pipeline_config.pl - Will create a pipeline.layout and pipeline.config for the 'wait test' pipeline

=head1 SYNOPSIS

 USAGE: create_wait_pipeline_config.pl
       --input_file=/path/to/some/input.file
       --output=/path/to/transterm.file
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--template_directory,-t>
	Path of the template configuration and layout files used to create the pipeline config and layout files.

B<--output_directory, -o>
	Directory to write the pipeline.config and pipeline.layout files

B<--log,-l>
    Logfile.

B<--debug,-d>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

 DESCRIPTION

=head1  INPUT

    Describe the input

=head1 OUTPUT

    Describe the output

=head1  CONTACT

    Shaun Adkins
    sadkins@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Basename;
use XML::Writer;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;

my $outdir = ".";
my $template_directory = "/local/devel/sadkins/repo/workflow/project_saved_templates";
####################################################

my %options;
my $pipelines = {
		'wait' => 'wait_test',
};


# Allow program to run as module for unit testing if necessary
main() unless caller();
exit(0);

sub main {
	my $results = GetOptions (\%options,
						  "template_directory|t=s",
						  "output_directory|o=s",
						  "log|l=s",
						  "debug=i",
						  "help"
					);

    &check_options(\%options);

	my $pipeline_layout;
	my $pipeline_config;
	# The file that will be written to
	$pipeline_layout = $outdir."/wait.layout";
	$pipeline_config = $outdir."/wait.config";

	# File handles for files to be written
	open( my $plfh, "> $pipeline_layout") or &_log($ERROR, "Could not open $pipeline_layout for writing: $!");
	# Since the pipeline.layout is XML, create an XML::Writer
	my $layout_writer = new XML::Writer( 'OUTPUT' => $plfh, 'DATA_MODE' => 1, 'DATA_INDENT' => 3 );
	# Write the pipeline.layout file
	&write_pipeline_layout( $layout_writer, sub {
		my ($writer) = @_;
		&write_include($writer, $pipelines->{'wait'});
	});
	# end the writer
	$layout_writer->end();

	my %config;
	# Write the pipeline config file
	&add_config( \%config, $pipelines->{'wait'} );
	# open config file for writing
	open( my $pcfh, "> $pipeline_config") or &_log($ERROR, "Could not open $pipeline_config for writing: $!");
	# Write the config
	&write_config( \%config, $pcfh );

	# close the file handles
	close($plfh);
	close($pcfh);

	my $mode = 0755;
	chmod $mode, $pipeline_config;
	chmod $mode, $pipeline_layout;

	print "Wrote $pipeline_layout and $pipeline_config for 'wait' pipeline\n";
}

### UTILITY SUBROUTINES ###

sub write_config {
    my ($config, $fh) = @_;

    # Make sure this section is first
    &write_section( 'global', $config->{'global'}, $fh );
    delete( $config->{'global'} );

    foreach my $section ( keys %{$config} ) {
        &write_section( $section, $config->{$section}, $fh );
    }
}

sub write_section {
    my ($section, $config, $fh) = @_;
    print $fh "[$section]\n";

    foreach my $k ( sort keys %{$config} ) {
      print $fh "$k=$config->{$k}\n";
    }
    print $fh "\n";
}

sub add_config {
    my ($config, $subpipeline, $config_name) = @_;
	print $template_directory, "\t", $subpipeline, "\n";
    my $pc = "$template_directory/$subpipeline/pipeline.config";
    $pc = "$template_directory/$subpipeline/$config_name" if( $config_name );
    open(IN, "< $pc") or &_log($ERROR, "Could not open $pc for reading: $!");

    my $section;
    while(my $line = <IN> ) {
        chomp( $line );
        next if( $line =~ /^\s*$/ || $line =~ /^\;/ );

        if( $line =~ /^\[(.*)\]/ ) {
            $section = $1;
        } elsif( $line =~ /(\$\;.*\$\;)\s*=\s*(.*)/ ) {
            &_log($ERROR, "Did not find section before line $line") unless( $section );
            $config->{$section} = {} unless( exists( $config->{$section} ) );
            $config->{$section}->{$1} = $2;
        }

    }

    close(IN);
}

sub write_parallel_commandSet {
    my ($writer, $block) = @_;
    $writer->startTag("commandSet", 'type' => 'parallel');
    $writer->dataElement("state","incomplete");
    $block->($writer);
    $writer->endTag("commandSet");
}

sub write_pipeline_layout {
    my ($writer, $body) = @_;
    $writer->startTag("commandSetRoot",
                      "xmlns:xsi" => "http://www.w3.org/2001/XMLSchema-instance",
                      "xsi:schemaLocation" => "commandSet.xsd",
                      "type" => "instance" );

    $writer->startTag("commandSet",
                      "type" => "serial" );

    $writer->dataElement("state", "incomplete");
    $writer->dataElement("name", "start pipeline:");

    $body->($writer);

    $writer->endTag("commandSet");
    $writer->endTag("commandSetRoot");
}

sub write_include {
    my ($writer, $sub_pipeline, $pipeline_layout) = @_;
    $pipeline_layout = "pipeline.layout" unless( $pipeline_layout );
    my $sublayout = $template_directory."/$sub_pipeline/$pipeline_layout";
    &_log($ERROR, "Could not find sub pipeline layout $sublayout\n") unless( -e $sublayout );
    $writer->emptyTag("INCLUDE", 'file' => "$sublayout");
}

sub check_options {
   	my $opts = shift;
   	if( $opts->{'help'} ) {
       &_pod;
   	}

   	if( $opts->{'log'} ) {
       open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
   	}

   	$debug = $opts->{'debug'} if( $opts->{'debug'} );

   	$outdir = $opts->{'output_directory'} if( $opts->{'output_directory'} );
   	$template_directory = $opts->{'template_directory'} if( $opts->{'template_directory'} );

}

sub _log {
   my ($level, $msg) = @_;
   if( $level <= $debug ) {
      print STDOUT "$msg\n";
   }
   print $logfh "$msg\n" if( defined( $logfh ) );
   exit(1) if( $level == $ERROR );
}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
