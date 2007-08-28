package Prokaryotic::Pipeline::Statistics;

=head1 NAME

    Statistics.pm - A module designed to gather statistics for a prokaryotic pipeline run

=head1 Description

    Defines accessor methods to interesting statistics of a prokaryotic pipeline to monitor its progress.
    Including gene coverage, gene count and evidence gathering.

=head2 Constructor and initialization

    my $prok = new Prokaryotic::Pipeline::Statistics( {
        'pipeline_id' => 9234,
        'repository_root' => '/usr/local/annotation/MOORE'
    });

=head2 Class and object methods

=over 4

=cut

use strict;
use warnings;
use XML::Twig;
use Data::Dumper;
use File::OpenFile qw( open_file );

sub new {
    my $class = shift;
    my $self = {};

    #Bless it.
    bless $self, $class;

    #Initialize the data type.
    $self->_init($_[0]);

    #And Return
    return $self;
}


=item $prok->get_gene_coverage()

B<Description:> returns length of the genome covered by genes

B<Parameters:> $force = non zero value will force a reparse of bsml

B<Returns:> Length of genome covered by genes

=cut 
sub get_gene_coverage {
    my ($self, $force) = @_;
    my $retval;
    $retval = $self->{'gene_coverage'} if( $self->{'gene_coverage'} && !$force );

    unless( $retval ) {
        
        #Find the list file
        my $list_file = $self->_find_gene_bsml_list();
        die( "Cannot find a gene bsml list file to calculate gene coverage")
            unless( $list_file );

        my $lfh = open_file( $list_file, 'in' );
        chomp( my @files = <$lfh> );
        close($lfh);

        foreach my $bsml_file ( @files ) {
            $retval += $self->_get_non_overlapping_total_gene_length( $bsml_file );
        }

    }

    die("Did not find any gene coverage") unless($retval);

    return $retval;
}

sub get_gene_count {
    my ($self, $force) = @_;
    my $retval = 0;
    $retval = $self->{'gene_count'} if( $self->{'gene_count'} && !$force );

    unless( $retval ) {
        
        my $list_file = $self->_find_gene_bsml_list();
        last unless( $list_file );

        my $lfh = open_file( $list_file, 'in' );
        chomp( my @files = <$lfh> );

        foreach my $file ( @files ) {
            my $bfh = open_file( $file, "in" );

            my $seq_id = $1 if( $file =~ m|/(([^/\.]+\.){3})| );
            $seq_id =~ s/\.$//;
            my $asmbl_count = 0;

            while( my $line = <$bfh> ) {
                if( $line =~ /<Feature\s.*class=\"gene\"/ ) {
                    $retval++;
                    $asmbl_count++;
                }
            }
            close( $bfh );

            $self->{'gene_count_assemblies'}->{$seq_id} = $asmbl_count;
        }
        close( $lfh );

    }

    die("Did not find any genes for organism $self->{'genome_name'}")
        unless( $retval );

    $self->{'gene_count'} = $retval;
    return $retval;
}

sub get_genome_length {
    my ($self, $force) = @_;
    my $length = 0;

    unless( $force ) {
        $length = $self->{'genome_length'};
        return $length if( $length );
    }

    #Check in the input to split multifasta.  Get the config file
    #for the run.
    my $final_config = $self->{'repository_root'}."/workflow/runtime/split_multifasta/".
        $self->{'pipeline_id'}."_default/split_multifasta.default.final.config";

    unless( -e $final_config ) {
        warn("Could not find configuration file for split_multifasta run.  Expected: ".
             $final_config);
        return 0;
    }
    
    my $cfh = open_file( $final_config, 'in' );
    
    my $output_list;
    while( my $line = <$cfh> ) {
        chomp $line;
        
        if( $line =~ /\$\;OUTPUT_LIST\$\;=(.*)$/ ) {
            $output_list = $1;
            last;
        }

    }
    close( $cfh );

    die( "Could not find output list for split_multifasta.default" ) 
        unless( $output_list && -e $output_list);

    my $ffh = open_file( $output_list, 'in' );

    while( my $line = <$ffh> ) {
        chomp $line;
        $length += $self->_get_length_from_fasta( $line );
    }


    die( "Could not find genome length for organism $self->{'genome_name'}") unless( $length );
    return $length;
}

#This was in the stub I made, but I don't know what I meant by it.  
#Maybe I will figure it out along the way.
#
#This is why we comment first, code second.
sub get_genome_source {

}

#Ugh.  This is HUGE. HUUUUGE. Maybe it should be it's own module.
sub get_annotation_statistics {

}

sub get_ber_results_count {
    my ($self) = @_;
    my $pre_overlap_count = $self->_count_evidence_files( 'ber', 'pre_overlap_analysis',
                                                          'btab' );
    my $post_overlap_count = $self->_count_evidence_files( 'ber', 'post_overlap_analysis',
                                                           'btab' );
    my $total_ber_results = $pre_overlap_count + $post_overlap_count;
    return $total_ber_results;    
}

sub get_hmm_results_count {
    my ($self) = @_;
    my $pre_overlap_count = $self->_count_evidence_files( 'hmmpfam', 'pre_overlap_analysis',
                                                          'htab' );
    my $post_overlap_count = $self->_count_evidence_files( 'hmmpfam', 'post_overlap_analysis',
                                                           'htab' );
    my $total_hmmpfam_results = $pre_overlap_count + $post_overlap_count;
    return $total_hmmpfam_results;    
}

sub get_genome_name {
    my ($self) = @_;
    my $organism = "";
    my $abbreviation = "";

    if( $self->{'genome_name'} ) {
        return $self->{'genome_name'};
    }

    #This can be found in the pipeline summary config file among other places
    #Since the pipeline summary config file will determine what gets put in the 
    #output, I'll report that one.
    my $ps_config = $self->{'repository_root'}."/workflow/runtime/pipeline_summary/".
        $self->{'pipeline_id'}."_default/pipeline_summary.default.user.config";

    unless( -e $ps_config ) {
        warn("pipeline_summary config file ( $ps_config ) does not exist.  It should");
        return 0;
    }

    my $pfh = open_file( $ps_config, 'in' );
    while( my $line = <$pfh> ) {
        chomp $line;
        $organism = $1 if( $line =~ /\$\;ORGANISM\$\;=(.*)$/ );
        $abbreviation = $1 if( $line =~ /\$\;LOCUS_PREFIX\$\;=(.*)$/ );
    }

    die("Could not find the key \$\;ORGANISM\$\; in pipeline_summary config file ".
         "[$ps_config]") unless( $organism );
    $self->{'genome_name'} = $organism;

    
    die("Could not find the key \$\;LOCUS_PREFIX\$\; in pipeline summary config file ".
        "[$ps_config] and it should be there.") unless( $abbreviation );

    $self->{'abbreviation'} = $abbreviation;
    return $organism;
    
}

sub get_genome_abbreviation {
    my ($self) = @_;

    if( $self->{'abbreviation'} ) {
        return $self->{'abbreviation'};
    }

    $self->get_genome_name;
    
    die("Could not find genome abbreviation") unless( $self->{'abbreviation'} );

    return $self->{'abbreviation'};
    
}

sub get_asmbl_info {
    my ($self) = @_;
    my $retval;

    if( $self->{'asmbl_info'} ) {
        return $self->{'asmbl_info'};
    }

    #Find the config file for the split_multifasta component
    my $split_config = $self->{'repository_root'}."/workflow/runtime/split_multifasta/".
        $self->{'pipeline_id'}."_default/split_multifasta.default.final.config";

    warn("Could not find the split_multifasta.default [$split_config] config file to retrieve asmbl info")
        unless( -e $split_config );
    
    #Find the input_fsa file
    my $cfh = open_file( $split_config, 'in' );

    my $fsa_file;
    while( my $line = <$cfh> ) {
        chomp $line;
        if( $line =~ /\$\;INPUT_FILE\$\;=(.*)$/ ) {
            $fsa_file = $1;
            last;
        }
    }
    close $cfh;

    die("Could not find fsa file in split_multifasta config file") 
        unless( $fsa_file );

    my $cmd = "residues $fsa_file";
    
    open( CMD, "$cmd |") or die("Could not run $cmd ($!)");
    my @residues = <CMD>;
    close( CMD );

    foreach my $residue_line ( @residues ) {
        my ( $id, $length ) = split( /\s+/, $residue_line );
        $retval->{$id} = $length;
    }

    return $retval;

}

sub get_rna_count {
    my ($self) = @_;
    my ($tRNA_retval) = (0);

    if( $self->{'tRNA_count'} ) {
        return $self->{'tRNA_count'};
    }

    #This will count tRNA-scan output.
    my $tRNA_lists = $self->is_component_finished( 'tRNAscan-SE', 'find_tRNA' );
    unless( $tRNA_lists ) {
        warn("tRNA analysis has not finished yet");
        $tRNA_retval = 0;
    } else {
        my $list_file = $tRNA_lists->{'bsml'};

        my $lfh = open_file( $list_file, 'in' );

        while( my $bsml_file = <$lfh> ) {
            chomp( $bsml_file );

            my $bfh = open_file( $bsml_file, 'in' );
            my @lines = grep( { /Feature\s.*class=\"tRNA\"/ } <$bfh> );
            
            $tRNA_retval += scalar(@lines);

        }
        
    }

    die("Found no tRNAs.") unless( $tRNA_retval );

    $self->{'tRNA_count'} = $tRNA_retval;
    return $tRNA_retval;
            
}

sub is_component_finished {
    my ($self, $component, $output_token ) = @_;
    my $retval;

    #Get the directory
    my $component_out_directory = $self->{'repository_root'}."/output_repository/$component/".
        $self->{'pipeline_id'}."_$output_token";

    #Make sure that it's a valid directory
    unless( -d $component_out_directory ) {
        die("Directory [$component_out_directory] does not exist");
    }

    #Open the dir and find list files
    opendir( DIR, $component_out_directory );
    my @list_files = grep { /.*list/ } readdir(DIR);
    closedir( DIR );

    if( @list_files > 0 ) {

        foreach my $list_file ( @list_files ) {
            $retval->{$1} = 
                $component_out_directory."/".$list_file if( $list_file =~ /([^\.]+)\.list(\.gz)?$/ );
        }
        
    }
    
    return $retval;
}


############ PRIVATE ######################
sub _init {
    my ($self, $args) = @_;

    $self->{'pipeline_id'} = $args->{'pipeline_id'};
    $self->{'repository_root'} = $args->{'repository_root'};

    #This defines the components that generate output which describe
    #genes in bsml
    $self->{'_gene_bsml_components'} = 
        ['intergenic_analysis', 'promote_gene_prediction','auto_gene_curation', 'glimmer3'];  
    
}
sub _count_evidence_files {
    my ($self, $component, $output_token, $file_type ) = @_;
    $file_type = 'bsml' unless( $file_type );
    my $retval = 0;
    my $total = 0;

    if( $self->{"${component}_${output_token}_${file_type}_nonzero"} ) {
        return $self->{"${component}_${output_token}_${file_type}_nonzero"};
    }

    #Find the finished computes
    #If the list is present
    my $analysis = $self->{'repository_root'}."/output_repository/$component/".
        $self->{'pipeline_id'}."_$output_token/$component.$file_type.list";

    unless( $analysis ) {
        warn("Could not find analysis list (analysis probably not finished");
        $self->{"${component}_${output_token}_${file_type}_nonzero"} = 0;
        $self->{"${component}_${output_token}_${file_type}_total"} = 0;
        return 0;
    }

    my $afh = open_file( $analysis, 'in' );

    while( my $file = <$afh> ) {
        chomp $file;
        $total++;
        $retval++ unless( -z $file );
    }

    $self->{"${component}_${output_token}_${file_type}_nonzero"} = $retval;
    $self->{"${component}_${output_token}_${file_type}_total"}   = $total;

    return $retval;
}

sub _find_gene_bsml_list {
    my ($self) = @_;

    #This will return the latest list of bsml files that contain
    #gene information ( from a completed component ).
    my @gene_descriptors = @{$self->{'_gene_bsml_components'}};
    
    #Open pipeline.xml
    my $pipeline_xml =  $self->{'repository_root'}."/workflow/runtime/pipeline/".
        $self->{'pipeline_id'}."/pipeline.xml";

    unless( -e $pipeline_xml ) {
        warn( "Could not find pipeline.xml: Repository root = ".$self->{'repository_root'}.
              " Pipeline id = ".$self->{'pipeline_id'}.". Expected $pipeline_xml");
        return 0;
    }

    my $fh = open_file( $pipeline_xml, 'in' );
    my $last_gene_descriptor;
    my $last_output_token;

    while ( my $line = <$fh> ) {
        local $" = "|";
        if( $line =~ /\<name\>(@gene_descriptors)\.(\w+)\<\/name\>/ ) {
            my $curr_component = $1;
            my $curr_out_token = $2;

            #Now find the state element
            while( my $find_state_elem = <$fh> ) {
                if( $find_state_elem =~ /<state>([^<]+)<\/state>/ ) {
                    if( $1 eq 'complete' ) {
                        $last_gene_descriptor = $curr_component;
                        $last_output_token = $curr_out_token;
                        last;
                    }
                }
            }

        }
    } 


    #Find the output list file (SHOULD GET THIS FROM THE CONFIG FILE)
    my $list_file = $self->{'repository_root'}."/output_repository/$last_gene_descriptor/".
        $self->{'pipeline_id'}."_$last_output_token/$last_gene_descriptor.bsml.list";

    unless( -e $list_file ) {
        #Auto gene curation has a slightly different style of naming the list file
        $list_file .= ".all" if( -e $list_file.".all" );
    }

    unless( -e $list_file ) {
        #If we still don't have a list file, warn and return nothing
        warn("Could not find list file for latest gene describing component".
             " Expected $list_file");
        $list_file = 0;
    }

    return $list_file;
        
 
}
sub _get_length_from_fasta {
    my ($self, $fsa) = @_;
    my $total_length = 0;
    
    unless( -e $fsa ) {
        warn("Alleged fsa file [$fsa] does not exist");
        return 0;
    }

    my $ffh = open_file( $fsa, 'in' );
    
    my ($seq_id, $seq) = ("","");
    while( my $line = <$ffh> ) {
        chomp $line;
        if( $line =~ /^>([^\s]+)/ ) {
            $seq_id = $1;

            if( length( $seq ) ) {
                $self->{$seq_id} = length( $seq );
                $total_length += length( $seq );
                $seq = "";
            }
        } else {
            unless( $seq_id ) {
                warn("Could not parse sequence id from fasta [$fsa]");
                return 0;
            }
            
            $seq .= $line;
        }
    }

    #make sure we get the last sequence
    $self->{$seq_id} = length( $seq );
    $total_length += length( $seq );

    return $total_length;
    
    
}

sub _get_non_overlapping_total_gene_length {
    my ($self, $bsml_file) = @_;
    my @coords = ();

    my $fh = open_file( $bsml_file, 'in' );

    my $flag=0;
    my ($start, $stop);
    while( <$fh> ) {
        $flag = 1 if( /Feature.*class=\"gene\"/ );
        if( $flag && /((startpos=\"(\d+).*endpos=\"(\d+).*)|(endpos=\"(\d+).*startpos=\"(\d+).*))/) {
            ($start, $stop) = ($3) ? ($3,$4) : ($6, $7 );
            
            die("Didn't get start and stop") unless( $start && $stop && $start != 0 && $stop != 0 );

            my @tmp = sort { $a <=> $b } ( $start, $stop );
            push( @coords, \@tmp );
            $flag = 0;
        }
    }

    close($fh);
    
    my @sorted = sort { $a->[0] <=> $b->[0] } @coords;
    my $length_total = 0;

    my ($first, $second);
    for( my $i = 1; $i < @sorted; $i++ ) {
        ($first, $second) = ($sorted[$i-1], $sorted[$i]) unless( $first && $second );
        my @all = sort { $a <=> $b } ( @{$first}, @{$second} );
        my $overlap = ( ($first->[1] - $first->[0]) + ($second->[1] - $second->[0]) - ( $all[3] - $all[0] ) );

        #This happens when the genes overlap
        if( $overlap > 0 ) {
            $first = [$all[0],$all[3]];
        } else {
            $length_total += ($first->[1] - $first->[0]);
            $first = $second;
        }

        $second = $sorted[$i+1];

    }

    $length_total += $first->[1] - $first->[0];
    return $length_total;
}
1;
