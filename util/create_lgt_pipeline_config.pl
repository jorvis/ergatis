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

B<--sra_id,-s>
	Valid ID from the Sequence Read Archive

B<--bam_input,-B>
	Valid path to a BAM input file.  Either this or the SRA ID must be provided

B<--donor_reference,-d>
	Path to the donor reference fasta file, or list file (ends in .list).  If not provided, the script assumes this is a host-only LGTSeek run.  If the reference has already been indexed by BWA, the index files must be in the same directory as the reference(s).

B<--host_reference,-h>
	Path to the donor reference fasta file, or list file (ends in .list).  If not provided, the script assumes this is a donor-only LGTSeek run.If the reference has already been indexed by BWA, the index files must be in the same directory as the reference(s).

B<--refseq_reference,-r>
	Path to the RefSeq reference fasta file, or list file (ends in .list).  If the reference has already been indexed by BWA, the index files must be in the same directory as the reference(s).

B<--build_indexes,-B>
	If the flag is enabled, will build indexes using BWA in the pipeline.  If you are using pre-build indexes it is important they are compatible with the version of BWA running in the pipeline (0.7.12 for internal Ergatis, 0.7.15 for Docker Ergatis)

B<--template_directory,-t>
	Path of the template configuration and layout files used to create the pipeline config and layout files.

B<--output_directory, -o>
	Directory to write the pipeline.config and pipeline.layout files

B<--no_pipeline_ids, -p>
	If the flag is enabled, do not add process IDs to the pipeline config and layout file names and just use whatever lgt.config and lgt.layout.

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
						  "sra_id|s=s",
						  "host_reference|h=s",
						  "donor_reference|d=s",
						  "refseq_reference|r=s",
						  "build_indexes|B",
						  "bam_input|b",
						  "template_directory|t=s",
						  "output_directory|o=s",
						  "data_directory|O=s",
						  "no_pipeline_id|p",
						  "log|l=s",
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

	# If the starting point is BAM input, then use that.
	if ($options{bam_input}) {
		$config{"lgt_bwa recipient"}->{'$;QUERY_FILE$;'} = $options{bam_input};
		$config{"lgt_bwa recipient"}->{'$;PAIRED;'} = 1;
	} else {
		$config{"global"}->{'$;SRA_RUN_ID$;'} = $options{sra_id};
	}
	# Default use case (good donor and good host), we just want two specific list files.  For donor and host-only cases, we want other specific list files
	$config{"lgt_bwa_post_process default"}->{'$;SKIP_WF_COMMAND$;'} = 'create single-map BAM file list,create no-map BAM file list';

	if ($donor_only) {
		# In donor-only alignment cases, we do not keep the 'MM' matches, so no microbiome run
		if ($options{bam_input}) {
			$config{"lgt_bwa donor"}->{'$;QUERY_FILE$;'} = $options{bam_input};
			$config{"lgt_bwa donor"}->{'$;PAIRED;'} = 1;
		} else {
			$config{"lgt_bwa donor"}->{'$;QUERY_FILE$;'} = '$;REPOSITORY_ROOT$;/output_repository/sra2fastq/$;PIPELINEID$;_default/sra2fastq.list';
		}
		$config{"lgt_bwa_post_process default"}->{'$;RECIPIENT_FILE_LIST$;'} = '';
		$config{"lgt_bwa_post_process default"}->{'$;SKIP_WF_COMMAND$;'} = 'create LGT host BAM file list,create microbiome BAM file list,create no-map BAM file list';
		$config{"filter_dups_lc_seqs lgt"}->{'$;INPUT_FILE_LIST$;'} = '$;REPOSITORY_ROOT$;/output_repository/lgt_bwa_post_process/$;PIPELINEID$;_default/lgt_bwa_post_process.single_map.bam.list';
	} else {
		# Only add host-relevant info to config if we are aligning to a host
		if ($options{host_reference} =~ /list$/) {
			$config{"global"}->{'$;HOST_LIST$;'} = $options{host_reference};
		} else {
			$config{"global"}->{'$;HOST_REFERENCE$;'} = $options{host_reference};
		}
		# The mpileup component needs the host reference to serve as a reference here too
		$config{'lgt_mpileup lgt'}->{'$;FASTA_REFERENCE$;'} = $options{host_reference};
	}

	if ($host_only) {
		# In a host-only run, we do not keep the 'MM' matches, so no microbiome run

		if ($options{bam_input}) {
			$config{"lgt_bwa recipient"}->{'$;QUERY_FILE$;'} = $options{bam_input};
			$config{"lgt_bwa recipient"}->{'$;PAIRED;'} = 1;			
		} else {
			$config{"lgt_bwa recipient"}->{'$;QUERY_FILE$;'} = '$;REPOSITORY_ROOT$;/output_repository/sra2fastq/$;PIPELINEID$;_default/sra2fastq.list';
		}	# I think this if/else block is not necessary.
		$config{"lgt_bwa_post_process default"}->{'$;DONOR_FILE_LIST$;'} = '';
		$config{"lgt_bwa_post_process default"}->{'$;SKIP_WF_COMMAND$;'} = 'create LGT host BAM file list,create microbiome BAM file list';
		$config{"filter_dups_lc_seqs lgt"}->{'$;INPUT_FILE_LIST$;'} = '$;REPOSITORY_ROOT$;/output_repository/lgt_bwa_post_process/$;PIPELINEID$;_default/lgt_bwa_post_process.single_map.bam.list';
		$config{"filter_dups_lc_seqs mb"}->{'$;INPUT_FILE_LIST$;'} = '$;REPOSITORY_ROOT$;/output_repository/lgt_bwa_post_process/$;PIPELINEID$;_default/lgt_bwa_post_process.no_map.bam.list';
	} else {
		# Only add donor-relevant info to config if we are aligning to a donor
		if ($options{donor_reference} =~/list$/) {
			$config{"global"}->{'$;DONOR_LIST$;'} = $options{donor_reference};
		} else {
			$config{"global"}->{'$;DONOR_REFERENCE$;'} = $options{donor_reference};
		}
	}

	# If we are indexing references in the pipeline, we need to change some config inputs
	if ($included_subpipelines{'indexing'}) {

		# Change the Refseq reference for lgt_bwa
		$config{'lgt_bwa lgt'}->{'$;INPUT_FILE$;'} = '';
		$config{'lgt_bwa lgt'}->{'$;INPUT_FILE_LIST$;'} = '$;REPOSITORY_ROOT$;/output_repository/lgt_build_bwa_index/$;PIPELINEID$;_refseq/lgt_build_bwa_index.fsa.list';
		unless ($donor_only) {
			$config{'lgt_bwa mb'}->{'$;INPUT_FILE$;'} = '';
			$config{'lgt_bwa mb'}->{'$;INPUT_FILE_LIST$;'} = '$;REPOSITORY_ROOT$;/output_repository/lgt_build_bwa_index/$;PIPELINEID$;_refseq/lgt_build_bwa_index.fsa.list';

			# Change the host reference for lgt_bwa
			$config{'lgt_bwa recipient'}->{'$;INPUT_FILE$;'} = '';
			$config{'lgt_bwa recipient'}->{'$;INPUT_FILE_LIST$;'} = '$;REPOSITORY_ROOT$;/output_repository/lgt_build_bwa_index/$;PIPELINEID$;_recipient/lgt_build_bwa_index.fsa.list';

			$config{'lgt_bwa validation'}->{'$;INPUT_FILE$;'} = '';
			$config{'lgt_bwa validation'}->{'$;INPUT_FILE_LIST$;'} = '$;REPOSITORY_ROOT$;/output_repository/lgt_build_bwa_index/$;PIPELINEID$;_recipient/lgt_build_bwa_index.fsa.list';
		}

		unless ($host_only) {
			# Change the donor reference for lgt_bwa
			$config{'lgt_bwa donor'}->{'$;INPUT_FILE$;'} = '';
			$config{'lgt_bwa donor'}->{'$;INPUT_FILE_LIST$;'} = '$;REPOSITORY_ROOT$;/output_repository/lgt_build_bwa_index/$;PIPELINEID$;_donor/lgt_build_bwa_index.fsa.list';
		}
	} else {
		# If not building indexes, delete reference to lgt_build_bwa_index config
		delete $config{"lgt_build_bwa_index refseq"};
	}

# If we are passing a directory to store important output files, then change a few parameters
if ($options{data_directory}){
	if ($included_subpipelines{'sra'}){
		$config{'global'}->{'$;DATA_DIR$;'} = $options{data_directory};
	}
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

	&_log($ERROR, "ERROR - Cannot specify both an SRA ID and a BAM input file.  Choose one.") if ($opts->{sra_id} && $opts->{bam_id});
	&_log($ERROR, "ERROR - Must specify either an SRA ID or a BAM input file.") if (!($opts->{sra_id} || $opts->{bam_id}));

	if ($opts->{'sra_id'}) {
   		$included_subpipelines{sra} = 1;
	}
	if ($opts->{'bam_input'}) {
		$included_subpipelines{sra} = 0;
	}

	# If donor reference is not present, then we have a host-only run
	$host_only = 1 unless ($opts->{'donor_reference'});

	# If host reference is not present, then we have a donor-only run
	$donor_only = 1 unless ($opts->{'host_reference'});

	# If we need to build BWA reference indexes, then set option
	$included_subpipelines{indexing} = 1 if ( $opts->{'build_indexes'} );

	&_log($ERROR, "ERROR - Need either a host_reference, a donor_reference or both provided") if ($donor_only && $host_only);

   print STDOUT "Perform alignments to the donor reference only.\n" if ($donor_only);
   print STDOUT "Perform alignments to the host reference only.\n" if ($host_only);
   print STDOUT "Perform BWA reference indexing in pipeline.\n" if ($included_subpipelines{indexing});
   print STDOUT "Starting point is BAM input.\n" if $opts->{bam_input};
   print STDOUT "Starting point is SRA ID. \n" if $opts->{sra_id};

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
