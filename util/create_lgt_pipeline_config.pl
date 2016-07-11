#!/usr/bin/env perl
package main;    # Change this to reflect name of script.  If you want to use as module, make sure to save as .pm instead of .pl

=head1 NAME

create_lgt_pipeline_config.pl - Will create a pipeline.layout and pipeline.config for selected sub-pipelines of
    the automated lateral genome transfer pipeline

=head1 SYNOPSIS

 USAGE: create_lgt_pipeline_config.pl
       --input_file=/path/to/some/input.file
       --output=/path/to/transterm.file
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>

B<--output_file,-o>

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
my $template_directory = "/local/projects/ergatis/package-latest/global_pipeline_templates";
my %included_subpipelines = ();
my $donor_only = 0;
my $host_only = 0;
####################################################

my %options;
my $pipelines = {
	    'sra' => 'LGT_Seek_Pipeline_SRA',
		'indexing' => 'LGT_Seek_Pipeline_BWA_Index',
		'lgtseek' => 'LGT_Seek_Pipeline'
};


# Allow program to run as module for unit testing if necessary
main() unless caller();
exit(0);

sub main {
	my $results = GetOptions (\%options,
						  #"lgtseek|l=i",
						  "sra_id|s=s",
						  "host_reference|h=s",
						  "donor_reference|d=s",
						  "refseq_reference|r=s",
                          "template_directory|t=s",
                          "output_directory|o=s",
						  "no_pipeline_id|p",
                          "log|L=s",
                          "debug=i",
                          "help"
                          );

    &check_options(\%options);

	my $pipeline_layout;
	my $pipeline_config;
	# The file that will be written to
	if ($options{no_pipeline_id}){
		$pipeline_layout = $outdir."/lgt.layout";
		$pipeline_config = $outdir."/lgt.config";
	} else {
		$pipeline_layout = $outdir."/pipeline.$$.layout";
		$pipeline_config = $outdir."/pipeline.$$.config";
	}
	# File handles for files to be written
	open( my $plfh, "> $pipeline_layout") or &_log($ERROR, "Could not open $pipeline_layout for writing: $!");

	# Since the pipeline.layout is XML, create an XML::Writer
	my $layout_writer = new XML::Writer( 'OUTPUT' => $plfh, 'DATA_MODE' => 1, 'DATA_INDENT' => 3 );

	# Write the pipeline.layout file
	&write_pipeline_layout( $layout_writer, sub {
		my ($writer) = @_;
   		&write_include($writer, $pipelines->{'sra'}) if( $included_subpipelines{'sra'} );
		# Use the right layout file if this run is donor-only, or both donor/host alignment
		if ($donor_only) {
   			&write_include($writer, $pipelines->{'indexing'}, "pipeline.donor_only.layout") if( $included_subpipelines{'indexing'} );
   			&write_include($writer, $pipelines->{'lgtseek'}, "pipeline.donor_only.layout") if( $included_subpipelines{'lgtseek'} );
		} elsif ($host_only) {
   			&write_include($writer, $pipelines->{'indexing'}, "pipeline.host_only.layout") if( $included_subpipelines{'indexing'} );
   			&write_include($writer, $pipelines->{'lgtseek'}, "pipeline.host_only.layout") if( $included_subpipelines{'lgtseek'} );
		} else {
   			&write_include($writer, $pipelines->{'indexing'}) if( $included_subpipelines{'indexing'} );
   			&write_include($writer, $pipelines->{'lgtseek'}) if( $included_subpipelines{'lgtseek'} );
		}
	});

	# end the writer
	$layout_writer->end();

	my %config;

	# Write the pipeline config file
	foreach my $sp ( keys %included_subpipelines ) {
		if ($sp eq 'sra') {
			&add_config( \%config, $pipelines->{ $sp } );
		} else {
			if ($donor_only) {
			    &add_config( \%config, $pipelines->{ $sp }, "pipeline.donor_only.config");
			} elsif ($host_only) {
			    &add_config( \%config, $pipelines->{ $sp }, "pipeline.host_only.config");
			} else {
			    &add_config( \%config, $pipelines->{ $sp } );
			}
		}
	}

	# Now do some global config modification
	if ($options{refseq_reference} =~ /list$/) {
		$config{"global"}->{'$;REFSEQ_LIST$;'} = $options{refseq_reference};
		# Each individual genome is small enough to use 'is' instead of 'btwsw'
		$config{"lgt_build_bwa_index refseq"}->{'$;ALGORITHM$;'} = "is";
	} else {
		$config{"global"}->{'$;REFSEQ_REFERENCE$;'} = $options{refseq_reference};
	}

	$config{"global"}->{'$;SRA_RUN_ID$;'} = $options{sra_id};
	# Default use case (good donor and good host), we just want two specific list files.  For donor and host-only cases, we want other specific list files
	$config{"lgt_bwa_post_process default"}->{'$;SKIP_WF_COMMAND$;'} = 'create single map BAM file list,create no map BAM file list';

	unless ($donor_only) {
	# Only add host-relevant info to config if we are aligning to a host
		if ($options{host_reference} =~ /list$/) {
			$config{"global"}->{'$;HOST_LIST$;'} = $options{host_reference};
		} else {
			$config{"global"}->{'$;HOST_REFERENCE$;'} = $options{host_reference};
		}
	} else {
		# In donor-only alignment cases, we do not keep the 'MM' matches, so no microbiome run
		$config{"lgt_bwa donor"}->{'$;QUERY_FILE$;'} = '$;REPOSITORY_ROOT$;/output_repository/sra2fastq/$;PIPELINEID$;_default/sra2fastq.list';
		$config{"lgt_bwa_post_process default"}->{'$;RECIPIENT_FILE_LIST$;'} = '';
		$config{"lgt_bwa_post_process default"}->{'$;SKIP_WF_COMMAND$;'} = 'create LGT BAM file list,create microbiome BAM file list,create no map BAM file list';
		$config{"filter_dups_lc_seqs lgt"}->{'$;INPUT_FILE_LIST$;'} = '$;REPOSITORY_ROOT$;/output_repository/lgt_bwa_post_process/$;PIPELINEID$;_default/lgt_bwa_post_process.single_map.bam.list';
	}

	unless ($host_only) {
	# Only add donor-relevant info to config if we are aligning to a donor
		if ($options{donor_reference} =~/list$/) {
			$config{"global"}->{'$;DONOR_LIST$;'} = $options{donor_reference};
		} else {
			$config{"global"}->{'$;DONOR_REFERENCE$;'} = $options{donor_reference};
		}
	} else {
		# In a host-only run, we do not keep the 'MM' matches, so no microbiome run
		$config{"lgt_bwa recipient"}->{'$;QUERY_FILE$;'} = '$;REPOSITORY_ROOT$;/output_repository/sra2fastq/$;PIPELINEID$;_default/sra2fastq.list';
		$config{"lgt_bwa_post_process default"}->{'$;DONOR_FILE_LIST$;'} = '';
		$config{"lgt_bwa_post_process default"}->{'$;SKIP_WF_COMMAND$;'} = 'create LGT BAM file list,create microbiome BAM file list';
		$config{"filter_dups_lc_seqs lgt"}->{'$;INPUT_FILE_LIST$;'} = '$;REPOSITORY_ROOT$;/output_repository/lgt_bwa_post_process/$;PIPELINEID$;_default/lgt_bwa_post_process.single_map.bam.list';
		$config{"filter_dups_lc_seqs mb"}->{'$;INPUT_FILE_LIST$;'} = '$;REPOSITORY_ROOT$;/output_repository/lgt_bwa_post_process/$;PIPELINEID$;_default/lgt_bwa_post_process.no_map.bam.list';
	}

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

	print "Wrote $pipeline_layout and $pipeline_config for LGT pipeline\n";
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

   	foreach my $req ( qw(sra_id refseq_reference) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   	}

   	$outdir = $opts->{'output_directory'} if( $opts->{'output_directory'} );
   	$template_directory = $opts->{'template_directory'} if( $opts->{'template_directory'} );
   	$included_subpipelines{lgtseek} = 1;	# LGTSeek is required... duh!
   	# Currently these will be default until I make the pipeline more sophisticated in accepting fastq files in place of the SRA template and/or previously created BWA indexes in place of that template
   	$included_subpipelines{sra} = 1;
   	$included_subpipelines{indexing} = 1;

	# If donor reference is not present, then we have a host-only run
	$host_only = 1 unless ($opts->{'donor_reference'});

	# If host reference is not present, then we have a donor-only run
	$donor_only = 1 unless ($opts->{'host_reference'});

	&_log($ERROR, "Cannot specify both 'donor_only' and 'host_only' options.  Choose either, or none") if ($donor_only && $host_only);

   print STDOUT "Perform alignments to the donor reference only.\n" if ($donor_only);
   print STDOUT "Perform alignments to the host reference only.\n" if ($host_only);

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
