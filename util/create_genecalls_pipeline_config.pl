#!/usr/bin/env perl

=head1 NAME

create_genecalls_pipeline_config.pl - Will create a pipeline.layout and pipeline.config for the Prokaryotic Annotation Genecalls v2 pipeline

=head1 SYNOPSIS

 USAGE: create_genecalls_pipeline_config.pl
     [ --templates_directory|-t=/path/to/ergatis/global_pipeline_templates
       --output_directory=/path/to/directory
       --load|-l=1
       --log|-L=/path/to/file.log
       --debug|-D=3
       --help
     ]

=head1 OPTIONS

B<--templates_directory,-t>
    The directory get the templates from

B<--output_directory,-o>
    The directory of where to write the output pipeline.layout and pipeline.config file.
    Defaults to current directory.

B<--load,-l>
    Default = 0
    To load the database into a chado instance

B<--log,-L>
    Logfile.

B<--debug,-D>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION
    
    This script will combine a series of sub-pipelines related to the prokaryotic
    annotation engine genecalls pipeline and create a pipeline.layout and pipeline.config file.
    
    This script will include sub-pipelines from the provided templates_directory and create
    a new pipeline. The following pipelines will be looked for by this script and a directory for
    each is expected in the templates_directory:
    
    Prokaryotic_Annotation_Pipeline_Genecalls_v2
    Prokaryotic_Annotation_Pipeline_Chado_Loading
    
=head1  INPUT
    No input files are required. Just the options and an output directory if
    necessary.

=head1 OUTPUT


=head1  CONTACT

    Shaun Adkins
    sadkins@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use XML::Writer;
use Pod::Usage;
use Data::Dumper;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
my $outdir = ".";
my $template_directory = "/local/projects/ergatis/package-nightly/global_pipeline_templates";
my %included_subpipelines = ();
####################################################

## Maps sub-pipeline names
my $pipelines = {
				 'genecalls' => 'Prokaryotic_Annotation_Pipeline_Genecalls_v2',
				 'load' => 'Prokaryotic_Annotation_Pipeline_Chado_Loading',
				};

my %options;
my $results = GetOptions (\%options,
                          "load|l=s",
						  "template_directory|t=s",
						  "output_directory|o=s",
						  "log|L=s",
						  "debug|d=s",
						  "help|h"
						 );

&check_options(\%options);

# The file that will be written to
my $pipeline_layout = $outdir."/pipeline.$$.layout";
my $pipeline_config = $outdir."/pipeline.$$.config";

# File handles for files to be written
open( my $plfh, "> $pipeline_layout") or &_log($ERROR, "Could not open $pipeline_layout for writing: $!");

# Since the pipeline.layout is XML, create an XML::Writer
my $layout_writer = new XML::Writer( 'OUTPUT' => $plfh, 'DATA_MODE' => 1, 'DATA_INDENT' => 3 );

# Write the pipeline.layout file
&write_pipeline_layout( $layout_writer, sub {
    my ($writer) = @_;
    &write_include($writer, $pipelines->{'genecalls'});
    &write_include($writer, $pipelines->{'load'}) if $included_subpipelines{'load'};
});

# end the writer
$layout_writer->end();

my %config;

&add_config( \%config, $pipelines->{'genecalls'} );
&add_config( \%config, $pipelines->{'load'} ) if $included_subpipelines{'load'};

# open config file for writing
open( my $pcfh, "> $pipeline_config") or &_log($ERROR, "Could not open $pipeline_config for writing: $!");

# Write the config
&write_config( \%config, $pcfh );
# close the file handles
close($plfh);
close($pcfh);


print "Wrote $pipeline_layout and $pipeline_config\n";
exit 0;

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
   
   &_pod if( $opts->{'help'} );
   open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)") if( $opts->{'log'} );
   
    $outdir = $opts->{'output_directory'} if( $opts->{'output_directory'} );
    $template_directory = $opts->{'template_directory'} if( $opts->{'template_directory'} );
    $included_subpipelines{'load'} = 1 if( $opts->{'load'} );
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
