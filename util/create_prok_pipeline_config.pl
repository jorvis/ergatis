#!/usr/bin/env perl

=head1 NAME

create_prok_pipeline_config.pl - Will create a pipeline.layout and pipeline.config for selected sub-pipelines of 
    the automated prokaryotic annotation pipeline

=head1 SYNOPSIS

 USAGE: create_prok_pipeline_config.pl
       --pseudomolecule|-p
       --rna_predictions|-r
       --gene_prediction|-g
       --annotation|-a
       --load|-l
     [ --templates_directory|-t=/path/to/ergatis/global_pipeline_templates
       --output_directory=/path/to/directory
       --log|-L=/path/to/file.log
       --debug|-D=3
       --help
     ]

=head1 OPTIONS

B<--pseudomolecule,-p>
    Default = 0. For annotation on a pseudomolecule

B<--rna_predictions,-r>
    Default = 1. To include RNA predictions
    
B<--gene_prediction,-g>
    Possible options: [glimmer, prodigal, none].
    Default: glimmer
    If option none used, --pseudomolecule option is ignored 
    
B<--annotation,-a>
    Default = 1
    To include annotation of gene structures.
    If this is not chosen, --load option will be ignored.

B<--load,-l>
    Default = 0
    To load the database into a chado instance
    
B<--templates_directory,-t>
    The directory get the templates from

B<--output_directory,-o>
    The directory of where to write the output pipeline.layout and pipeline.config file.
    Defaults to current directory.

B<--log,-L>
    Logfile.

B<--debug,-D>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION
    
    This script will combine a series of sub-pipelines related to the prokaryotic
    annotation engine pipeline and create a pipeline.layout and pipeline.config file.
    The config file can then be configured with the correct options and a pipeline
    can be run. 
    
    This script will include sub-pipelines from the provided templates_directory and create
    a new pipeline. The following pipelines will be looked for by this script and a directory for
    each is expected in the templates_directory:
    
    Prokaryotic_Annotation_Pipeline_Pseudomolecule_Generation
    Prokaryotic_Annotation_Pipeline_Input_Processing
    Prokaryotic_Annotation_Pipeline_RNA_Predictions
    Prokaryotic_Annotation_Pipeline_Glimmer_Predictions
    Prokaryotic_Annotation_Pipeline_Prodigal_Predictions
    Prokaryotic_Annotation_Pipeline_Gene_Annotations
    Prokaryotic_Annotation_Pipeline_Chado_Loading

    If you run the script without any parameters, it will create a config file which contains the following
    sub-pipelines:

    Prokaryotic_Annotation_Pipeline_Input_Processing
    Prokaryotic_Annotation_Pipeline_RNA_Predictions
    Prokaryotic_Annotation_Pipeline_Glimmer_Predictions
    Prokaryotic_Annotation_Pipeline_Gene_Annotations
 
    There are some restrictions about which parts can be included with others. For example,
    to include pseudomolecule generation you also have to include gene predictions. If gene predictions 
    is set to 'none', then the --pseudomolecule option will be ignored. If annotation is not
    included, the --load option will be ignored.
    
=head1  INPUT
    No input files are required. Just the options and an output directory if
    necessary.

=head1 OUTPUT


=head1  CONTACT

    Kevin Galens
    kgalens@gmail.com

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
my $valid_gene_prediction_algorithms = {
    'glimmer' => 1,
    'prodigal'  => 1
};
my $gene_prediction = "glimmer";
my $template_directory = "/local/projects/ergatis/package-v1r30/global_pipeline_templates";
my %included_subpipelines = ();
####################################################

## Maps sub-pipeline names
my $pipelines = {
    'pseudomolecule' => 'Prokaryotic_Annotation_Pipeline_Pseudomolecule_Generation',
    'input_processing' => 'Prokaryotic_Annotation_Pipeline_Input_Processing',
    'rna_prediction'  => 'Prokaryotic_Annotation_Pipeline_RNA_Predictions',
    'glimmer' => 'Prokaryotic_Annotation_Pipeline_Glimmer_Predictions',
    'prodigal'  => 'Prokaryotic_Annotation_Pipeline_Prodigal_Predictions',
    'annotation'  => 'Prokaryotic_Annotation_Pipeline_Gene_Annotations',
    'load' => 'Prokaryotic_Annotation_Pipeline_Chado_Loading'
};

my %options;
my $results = GetOptions (\%options,
                            "pseudomolecule|p=s",
                            "rna_prediction|r=s",
                            "gene_prediction|g=s",
                            "annotation|a=s",
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
   &write_include($writer, $pipelines->{'pseudomolecule'}) if( $included_subpipelines{'pseudomolecule'} );
   &write_include($writer, $pipelines->{'input_processing'}) if( $included_subpipelines{'input_processing'} );
   
   if( $included_subpipelines{'rna_prediction'} || $included_subpipelines{'gene_prediction'} ne 'none' ) {
       &write_parallel_commandSet( $writer, sub {
           my ($writer) = @_;
           &write_include($writer, $pipelines->{'rna_prediction'}) if( $included_subpipelines{'rna_prediction'} );
           &write_include($writer, $pipelines->{$included_subpipelines{'gene_prediction'}}) 
               unless( $included_subpipelines{'gene_prediction'} eq 'none' );
       });
   }
   
   &write_include($writer, $pipelines->{'annotation'}) if( $included_subpipelines{'annotation'} );
   if( $included_subpipelines{'load'} ) {
       if( $included_subpipelines{'pseudomolecule'} ) {
           &write_include($writer, $pipelines->{'load'}, 'pipeline_pmarks.layout');
       } else {
           &write_include($writer, $pipelines->{'load'});
       }
   }
});

# end the writer
$layout_writer->end();

my %config;

# Write the pipeline config file
foreach my $sp ( keys %included_subpipelines ) {
    if( $sp eq 'gene_prediction' && $included_subpipelines{$sp} ne 'none' ) {
        print "Adding config for $sp\n";
        &add_config( \%config, $pipelines->{$included_subpipelines{$sp}} );
    } elsif( $sp eq 'load' && $included_subpipelines{$sp} && $included_subpipelines{'pseudomolecule'} ) {
        print "Adding config for $sp\n";
        &add_config( \%config, $pipelines->{ $sp }, "pipeline_pmarks.config" );
    } elsif( $sp ne 'gene_prediction' && $included_subpipelines{$sp} ) {
        print "Adding config for $sp\n";
        &add_config( \%config, $pipelines->{ $sp } );
    }
}

# we can be a little bit smarter here and configure some of the connections for the user.
# If entry point is gene prediction:
#    - if annotation
#       - set input of annotation to output of gene prediction
if( $included_subpipelines{'pseudomolecule'} ) {
     # delete the $;INPUT_PROCESSING_FASTA_LIST$;
    delete( $config{'global'}->{'$;INPUT_FSA_LIST$;'} );
    delete( $config{'global'}->{'$;INPUT_FSA_FILE$;'} );
    
    
    foreach my $section ( keys %config ) {
        foreach my $key ( keys %{$config{$section}} ) {
            if( $config{$section}->{$key} eq '$;INPUT_FSA_LIST$;' ) {
                $config{$section}->{$key} = '$;REPOSITORY_ROOT$;/output_repository/create_pseudomolecules'.
                    '/$;PIPELINEID$;_default/create_pseudomolecules.list';
            }
            if( $config{$section}->{$key} eq '$;INPUT_FSA_FILE$;' ) {
                $config{$section}->{$key} = '';
            }
        }
    }   
}


if( $included_subpipelines{'gene_prediction'} ne 'none' && $included_subpipelines{'annotation'} ) {
    # remove INPUT_BSML_LIST from global section
    delete( $config{'global'}->{'$;INPUT_BSML_LIST$;'} );
    

    # anywhere else you find the $;INPUT_BSML_LIST$; variable, replace it with the output of 
    # promote_gene_prediction.promote_prediction output list
    foreach my $section ( keys %config ) {
        foreach my $key ( keys %{$config{$section}} ) {
            if( $config{$section}->{$key} eq '$;INPUT_BSML_LIST$;' ) {
                $config{$section}->{$key} = '$;REPOSITORY_ROOT$;/output_repository/promote_gene_prediction'.
                    '/$;PIPELINEID$;_promote_prediction/promote_gene_prediction.bsml.list';
            }
        }
    }
}

# open config file for writing
open( my $pcfh, "> $pipeline_config") or &_log($ERROR, "Could not open $pipeline_config for writing: $!");

# Write the config
&write_config( \%config, $pcfh );

# close the file handles
close($plfh);
close($pcfh);


print "Wrote $pipeline_layout and $pipeline_config\n";

sub write_config {
    my ($config, $fh) = @_;
    
    # Make sure this section is first
    &write_section( 'global', $config{'global'}, $fh );
    delete( $config{'global'} );

    foreach my $section ( keys %{$config} ) {
        &write_section( $section, $config->{$section}, $fh );
    }
}

sub write_section {
    my ($section, $config, $fh) = @_;
    print $fh "[$section]\n";
    while( my ($k,$v) = each( %{$config} ) ) {
        print $fh "$k=$v\n";
    }
    print $fh "\n";
}

sub add_config {
    my ($config, $subpipeline, $config_name) = @_;
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
   
   # Make sure the value for the gene_prediction option is valid
   if( $opts->{'gene_prediction'} ) {
       &_log($ERROR, 
             "Value for gene_prediction option must of one of [none|".(keys %{$valid_gene_prediction_algorithms}).join("|")."]")
           unless( exists( $valid_gene_prediction_algorithms->{$opts->{'gene_prediction'}} ) 
                   || lc($opts->{'gene_prediction'}) eq 'none' );
       $gene_prediction = lc($opts->{'gene_prediction'});
   } elsif( exists( $opts->{'gene_prediction'} ) && $opts->{'gene_prediction'} == 0 ) {
       $gene_prediction = 'none';
   }
   
   $outdir = $opts->{'output_directory'} if( $opts->{'output_directory'} );
   $template_directory = $opts->{'template_directory'} if( $opts->{'template_directory'} );

   # This is where some logic goes about what parts actually get included into the pipeline.
   $included_subpipelines{'pseudomolecule'} = 1 if( $opts->{'pseudomolecule'} && $gene_prediction ne 'none' );
   $included_subpipelines{'rna_prediction'} = 1 unless( exists( $opts->{'rna_prediction'} ) && $opts->{'rna_prediction'} == 0 );
   $included_subpipelines{'gene_prediction'} = $gene_prediction;
   $included_subpipelines{'annotation'} = 1; # Annotation is required 
   $included_subpipelines{'load'} = 1 if( $opts->{'load'} );
      
   $included_subpipelines{'input_processing'} = 1 if( $included_subpipelines{'rna_prediction'} || 
                                                      $included_subpipelines{'gene_prediction'} ne 'none' );


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
