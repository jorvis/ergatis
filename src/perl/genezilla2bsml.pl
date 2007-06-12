#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

genezilla2bsml.pl - convert genezilla GFF output to BSML

=head1 SYNOPSIS

USAGE: genezilla2bsml.pl 
        --input_file=/path/to/genezilla.output.file.raw
        --output=/path/to/output.bsml
        --project=aa1 
        --fasta_file=/path/to/somefile.fsa 
        --id_repository=/path/to/repository
        --sourcename=sourcename
      

=head1 OPTIONS

B<--input_file,-f> 
    Input file from an genezilla scan.  -i, --input_list, will take
    in a list of input files, all of which will be stored in a single
    output bsml.

B<--output,-o> 
    Output BSML file

B<--project,-p> 
    Project ID.  Used in creating feature ids. 

B<--fasta_file,-a>
    Needed tp create a Seq-data-import element referencing this path.

B<--id_repository, -r>
    path to --project's id_repository

B<--sourcename,-s>
    Value to be used for the sourcename, required for the analysis section.  Most
    often should be the output_directory for the run.

B<--log,-l> 
    Log file

B<--help,-h> 
    This help message

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

=head1   DESCRIPTION

This script is used to convert the output from an genezilla search into BSML.

=head1 INPUT

You define the input file using the --input option.  This file does not need any
special file extension.  The regular output of genezilla looks like this:

    aa1.assembly.20285	genezilla	poly-A-signal	2654	2660	.	-	.	transgrp=1;
    aa1.assembly.20285	genezilla	final-exon	2748	2774	19.9	-	0	transgrp=1;
    aa1.assembly.20285	genezilla	initial-exon	3431	3958	19.05	-	0	transgrp=1;
    aa1.assembly.20285	genezilla	initial-exon	14359	14389	3.318	+	0	transgrp=2;
    aa1.assembly.20285	genezilla	final-exon	14445	14602	266	+	1	transgrp=2;
    aa1.assembly.20285	genezilla	poly-A-signal	15873	15879	.	-	.	transgrp=3;
    aa1.assembly.20285	genezilla	single-exon	15904	16293	38.33	-	0	transgrp=3;

which has the general format of:

    <seqname> <source> <feature> <start> <end> <score> <strand> <frame> transgrp=N;

This script currently assumes that the input GFF file only contains the output from
the analysis of a single FASTA sequence.  

=head1 OUTPUT

Base positions from the input file are renumbered so that positions start at zero.

=head1 CONTACT

    Jason Inman
    jinman@tigr.org

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

use Ergatis::Logger;
use Ergatis::IdGenerator;
use BSML::GenePredictionBsml;
use Chado::Gene;

my @inputFiles;
my $project;
my $output;
my $idMaker;
my $bsml;
my $data;
my $inputFsa;
my $sourcename;
my $debug;
my $length;

my %options = ();
my $results = GetOptions (\%options, 
			    'input_list|i=s',
                'input_file|f=s',
                'output|o=s',
                'project|p=s',
                'id_repository|r=s',
                'fasta_input|a=s',
                'sourcename|s=s',
                'log|l=s',
                'command_id=s',       ## passed by workflow
                'logconf=s',          ## passed by workflow (not used)
                'debug=s',
			    'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

## make sure all passed options are peachy
&check_parameters(\%options);

foreach my $file (@inputFiles) {

    $data = &parseGenezillaData($file);
    $bsml = &generateBsml($data);

}

$bsml->writeBsml($output);

exit(0);

# Store all the info we need in the data hash, 
sub parseGenezillaData {

    my $inFile = shift;

    ## open the input file for parsing
    open (my $ifh,"< $inFile") || $logger->logdie("can't open input file for reading");

    my $source_seq_name = '';
    my $current_group_name = '';
    my $current_transcript_id = '';
    my $previous_group_name = '';
    my $first_in_group = 1;
    my $last_gene_found = 0;
    my $comp_val;
    my $poly_end = '';
    my $min_exon_low = '';
    my $max_exon_high = '';
    my $genes;
    my %tmp;
    my @group_members = ();

    ## go through the file

    while (<$ifh>) {

        ## skip comment lines
        next if (/^##/);        

        ## set flag for having passed out of gene data:
        $last_gene_found++ if /^#/;

        my @cols = ();

        ## grab/modify a few values for later use
        unless ($last_gene_found) {

            chomp;

            @cols = split(/\t/);

            ## store the source sequence name
            $source_seq_name = $cols[0];

            ## count in interbase:
            $cols[3]--;

            ## store complement flag (1 is true, this is on reverse)
            $comp_val = ($cols[6] eq '+') ? 0 : 1;

            ## Get the group id.
            if ($cols[8] =~ /transgrp=(\d+)\;/) {
                $current_group_name = $1;
            } else {
                $logger->logdie("unrecognized format in attributes column: $cols[8]");
            }

            if ($cols[2] =~ /poly-A-signal/) {
                $poly_end = ($comp_val) ? $cols[3] : $cols[4];
            }

        }

        ## if column 9 (group, $cols[8]) is defined and is different than the last one we need to
        ## add the gene to our list, after also adding the polypeptide and transcript.
        ## Also, if we hit a # then we should add the gene, too, as it indicates we've
        ## reached the end of the gene calls.
        if ( ($last_gene_found) || (!$first_in_group && $current_group_name ne $previous_group_name) ) {

            ## Take care of the polypeptide
            my $type = 'polypeptide';
            my $typeid = $idMaker->next_id( 'type' => $type, 'project' => $project );
            %tmp = ('typeid'=> $typeid,
                    'low'   => $min_exon_low,
                    'high'  => $max_exon_high,
                    'comp'  => $comp_val,
                    'type'  => $type);

            push @group_members, {%tmp};

            ## Take care of the transcript.  Make sure we account for the possibility
            ## of a poly-A-signal.
            $type = 'transcript';
            my $low = 0;
            my $high = 0;

            if ($poly_end ne '') {
                $low = ($comp_val) ? $poly_end : $min_exon_low;
                $high = ($comp_val) ? $max_exon_high : $poly_end;
            } else {
                $low = $min_exon_low;
                $high = $max_exon_high;
            }
            $typeid = $idMaker->next_id( 'type' => $type, 'project' => $project);
            %tmp = ('typeid' => $typeid,
                    'low'    => $low,
                    'high'   => $high,
                    'comp'   => $comp_val,
                    'type'   => $type);

            push @group_members, {%tmp};

            ## Create the gene object here
            my $tmpGene = new Chado::Gene( $idMaker->next_id( 'type'  => 'gene',
                                                       'project' => $project ),
                                    $low, $high, $comp_val, $source_seq_name );

            ## Add the features to the gene object here:
            foreach my $feat (@group_members) {

                $tmpGene->addFeature($$feat{'typeid'},$$feat{'low'},$$feat{'high'},
                                       $$feat{'comp'},$$feat{'type'});

            }
 
            ## Form the group here.                             
            my $count = $tmpGene->addToGroup($tmpGene->getId, { 'all' => 1 });
            $logger->logdie("Nothing added to group") unless ($count);

            push (@{$genes}, $tmpGene);

            last if $last_gene_found;

            # reset flag for to generate new Chado::Gene object with next line
            $first_in_group = 1;

        }

        ## remember this group name
        $previous_group_name = $current_group_name;

        if ($first_in_group) {
            # we'll need to make a new Chado::Gene for each new transgrp set.

            $first_in_group = 0;
            $min_exon_low = '';
            $max_exon_high = '';
            $comp_val = '';
            $poly_end = '';
            @group_members = ();
        }

        # Handle exons:
        if ($cols[2] =~ /-exon$/) {
            foreach my $type( qw(exon CDS) ) {
                my $typeid = $idMaker->next_id( 'type' => $type,
                                                'project' => $project);
                %tmp = ('typeid' => $typeid,
                        'low'    => $cols[3],
                        'high'   => $cols[4],
                        'comp'   => $comp_val,
                        'type'   => $type);
                push @group_members, {%tmp};

            }

            # Adjust the boundaries for the transcript and possibly, polypeptide
            $min_exon_low = ($min_exon_low eq '') ? $cols[3] :
                            (($min_exon_low < $cols[3] ) ? $min_exon_low : $cols[3] );
            $max_exon_high = ($max_exon_high eq '') ? $cols[4] :
                             (($max_exon_high > $cols[4]) ? $max_exon_high : $cols[4]);

        }
       
        # Handle poly-A-signal sequences 
        if ($cols[2] =~ /poly-A-signal/) {
                my $type = 'polyA_signal_sequence';
                my $typeid = $idMaker->next_id( 'type' => $type,
                                                'project' => $project);
                %tmp = ('typeid' => $typeid, 
                        'low'    => $cols[3], 
                        'high'   => $cols[4], 
                        'comp'   => $comp_val, 
                        'type'   => $type);
                push @group_members, {%tmp};
        }

    }

    return $genes;

}

sub generateBsml {
    my $data = shift;

    #Create the document
    my $doc = new BSML::GenePredictionBsml( 'genezilla', $sourcename);

    foreach my $gene(@{$data}) {
        $doc->addGene($gene);
    }

    my $seqId;
    open(IN, "< $inputFsa") or $logger->logdie("Unable to open $inputFsa");
    while(<IN>) {
        #assume it's a single fasta file
        if(/^>([^\s+]+)/) {
            $seqId = $1;
            last;
        }
    }
    close(IN);

    my $addedTo = $doc->setFasta('', $inputFsa);
    $logger->logdie("$seqId was not a sequence associated with the gene") unless($addedTo);

    return $doc;

}


sub check_parameters {

    my $options = shift;

    my $error = "";

    if($options{'input_list'}) {
        $error .= "Option input_list ($options{'input_list'}) does not exist\n"
            unless(-e $options{'input_list'});
        open(IN, "< $options{'input_list'}") || $logger->logdie("Unable to open $options{'input_list'} ($!)");
        @inputFiles = <IN>;
        close(IN);
    }
    if($options{'input_file'}) {
        $error .= "Option input_file ($options{'input_file'}) does not exist\n" unless(-e $options{'input_file'});
        push(@inputFiles, $options{'input_file'});
    }
    unless($options{'input_list'} || $options{'input_file'}) {
        $error .= "Either option input_list or input_file is required\n";
    }

    unless($options{'project'}) {
        $error .= $options{'project'};
    } else {
        $project = $options{'project'};
    }

    unless($options{'output'}) {
        $error .= "Option output is required\n";
    } else {
        $output = $options{'output'};
    }

    unless($options{'id_repository'}) {
        $error .= "Option id_repository is required.  Please see Ergatis::IdGenerator ".
            "for details.\n";
    } else {
        $idMaker = new Ergatis::IdGenerator( 'id_repository' => $options{'id_repository'} );
        $idMaker->set_pool_size( 'exon'        => 20,
                                 'transcript'  => 20,
                                 'gene'        => 20,
                                 'polypeptide' => 20,
                                 'CDS'         => 20 );

    }

    unless($options{'fasta_input'}) {
        $error .= "Option fasta_input is required\n";
    } else {
        $error .= "$options{'fasta_input'} (fasta_input) does not exist\n"
            unless(-e $options{'fasta_input'});
        $inputFsa = $options{'fasta_input'};

        #Find the length of the sequence
        open(IN, "<$options{'fasta_input'}") or
            $logger->logdie("Can't open $options->{'fasta_input'}");
        my $seq = "";
        my $curId;
        while(<IN>) {
            chomp;
            if(/^>(\S+)/) {
                if(defined($curId)) {
                    $length->{$curId} = length($seq);
                    $seq = "";
                }
                $curId = $1;
            } else {
                $seq.= $_;
            }

        }

        $length->{$curId} = length($seq);
        close(IN);
    }

    ## get sourcedir for the analysis section:
    if ($options{'sourcename'}) {
        $sourcename = $options{'sourcename'};
    } else {
        $error .= "--sourcename is a required option.\n";
    }

    if($options{'debug'}) {
        $debug = $options{'debug'};
    }

    unless($error eq "") {
        $logger->logdie($error);
    }

}
