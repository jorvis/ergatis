#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

twinscan2bsml.pl - convert twinscan GTF output to BSML

=head1 SYNOPSIS

USAGE: twinscan2bsml.pl 
        --input_file=/path/to/twinscan.raw.file 
        --output=/path/to/output.bsml
        --project=aa1 
        --fasta_input=/path/to/somefile.fsa 
        --id_repository=/path/to/repository
        --sourcename=sourcename
        --programversion='current'

=head1 OPTIONS

B<--input_file,-f> 
    Input file file from a twinscan run.  -i, --input_list, will take in a list
    of input files, all of which will be stored in a single output bsml.

B<--output,-o> 
    Output BSML file (will be created, must not exist)

B<--project,-p> 
    Project ID.  Used in creating feature ids. 

B<--fasta_input,-a>
    Needed to create a Seq-data-import element referencing this path.

B<--id_repository,-r>
    path to --project's id_repository

B<--programversion,-v>
    Version string to be used as value for the analysis attribute 'programversion'

B<--sourcename,-s>
    Sourcename string to be used as value for the analysis attribute 'sourcename'
    Due to a silly hack in analysis2bsml.pl, you might need to put /dummy/dir at
    the end of it.

B<--log,-l> 
    Log file

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used to convert the output from a twinscan search into BSML.

=head1 INPUT

You define the input file using the --input_file option.  This file does not need any
special file extension.  The regular output of twinscan, less about a dozen comment
lines of header, looks like this:

eda1.assembly.1.1.0.seq iscan   CDS 62 293      6764    +       1       gene_id "eda1.assembly.1.1.0.seq.001"; transcript_id "eda1.assembly.1.1.0.seq.001.1";
eda1.assembly.1.1.0.seq iscan   stop_codon      294     296     .       +       0       gene_id "eda1.assembly.1.1.0.seq.001"; transcript_id "eda1.assembly.1.1.0.seq.001.1";
eda1.assembly.1.1.0.seq iscan   start_codon     460     462     .       +       0       gene_id "eda1.assembly.1.1.0.seq.002"; transcript_id "eda1.assembly.1.1.0.seq.002.1";
eda1.assembly.1.1.0.seq iscan   CDS 460 502     903     +       0       gene_id "eda1.assembly.1.1.0.seq.002"; transcript_id "eda1.assembly.1.1.0.seq.002.1";
eda1.assembly.1.1.0.seq iscan   CDS 513 565     427     +       2       gene_id "eda1.assembly.1.1.0.seq.002"; transcript_id "eda1.assembly.1.1.0.seq.002.1";
eda1.assembly.1.1.0.seq iscan   stop_codon      566     568     .       +       0       gene_id "eda1.assembly.1.1.0.seq.002"; transcript_id "eda1.assembly.1.1.0.seq.002.1";
eda1.assembly.1.1.0.seq iscan   start_codon     742     744     .       +       0       gene_id "eda1.assembly.1.1.0.seq.003"; transcript_id "eda1.assembly.1.1.0.seq.003.1";
eda1.assembly.1.1.0.seq iscan   CDS 742 802     411     +       0       gene_id "eda1.assembly.1.1.0.seq.003"; transcript_id "eda1.assembly.1.1.0.seq.003.1";
eda1.assembly.1.1.0.seq iscan   CDS 826 926     2058    +       2       gene_id "eda1.assembly.1.1.0.seq.003"; transcript_id "eda1.assembly.1.1.0.seq.003.1";
eda1.assembly.1.1.0.seq iscan   stop_codon      927     929     .       +       0       gene_id "eda1.assembly.1.1.0.seq.003"; transcript_id "eda1.assembly.1.1.0.seq.003.1";
eda1.assembly.1.1.0.seq iscan   stop_codon      1016    1018    .       -       0       gene_id "eda1.assembly.1.1.0.seq.004"; transcript_id "eda1.assembly.1.1.0.seq.004.1";
eda1.assembly.1.1.0.seq iscan   CDS 1019        1146    890     -       2       gene_id "eda1.assembly.1.1.0.seq.004"; transcript_id "eda1.assembly.1.1.0.seq.004.1";



which has the general format of:

    <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes]

where attributes contains:
    "gene_id", referring to the source input sequence followed by a unique gene id.
    "transcript_id", referring to the gene_id followed by a transcript id, which will
        always be '.1' for twinscan, until multiple predictions per gene (that is,
        alternative splices) are implemented.

This script currently assumes that the input GTF file only contains the output from
the analysis of a single FASTA sequence.

Also it should be noted that all transcripts will have, at a minimum, 2 rows: a CDS
and a stop_codon.  If there are at least 2 CDS predictions, an additional start_codon
is added to the results.  

=head1 OUTPUT


Base positions from the input file are renumbered so that positions start at zero and
reflect interbase numbering.  Also, the stop coordinates of each terminal CDS is
extended 3bp to include the stop codon.

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

### Some globals
my @inputFiles;
my $project;
my $output;
my $sourcename;
my $idMaker;
my $bsml;
my $data;
my $inputFsa;
my $debug;
my $length;
my $programversion;

my %options = ();
my $results = GetOptions (\%options, 
			    'input_list|i=s',
                'input_file|f=s',
                'output|o=s',
                'project|p=s',
                'id_repository|r=s',
                'fasta_input|a=s',
                'sourcename|s=s',
                'programversion|v=s',
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

    $data = &parseTwinscanData($file);
    $bsml = &generateBsml($data);

}
$bsml->writeBsml($output);

exit(0);


sub parseTwinscanData {

    my $inFile = shift;

    ## open the input file for parsing
    open (my $ifh, "< $inFile") || $logger->logdie("Can't open input file for reading");

    my $source_seq_name = '';
    my $current_group_name = '';
    my $current_transcript_id = '';
    my $previous_group_name = '';
    my $first_feat_found = 0;
    my $last_feat_found = 0;
    my $comp_val;
    my $min_exon_low;
    my $max_exon_high;
    my $genes;
    my %tmp;
    my @group_members = ();

    ## go through the data now

    while (<$ifh>) {

        ## skip comment lines
        next if (/^#/);

        ## set flag for having found last gene data;
        $last_feat_found++ if eof($ifh);

        ## grab/modify a few values for later use
        chomp;

        my @cols = split(/\t/);

        ## Get the gene id, ignore the transcript id for now.  It's simple the 
        ## gene_id . '.1' now, anyway.
        if ($cols[8] =~ /gene_id "(\S+);/) {
            $current_group_name = $1;
        } else {
            $logger->logdie("unrecognized format in attributes column: $cols[8]");
        }

        ## If this is a new gene id AND isn't the first one found, add the previous
        ## gene to out data hash.
        if ($current_group_name ne $previous_group_name) {

            if ($first_feat_found) {
                &add_group(\@group_members, $min_exon_low, $max_exon_high,
                             $comp_val, \$genes, $source_seq_name);
            } else {
                $first_feat_found++;
            }

            $min_exon_low = '';
            $max_exon_high = '';            
            $comp_val = '';
            @group_members = ();

            # make this group the one to look for:
            $previous_group_name = $current_group_name;

        }

        ## store the source sequence name
        $source_seq_name = $cols[0];

        ## count in interbase:
        $cols[3]--;

        ##store complement flas (1 is true, meaining feature lies on reverse strand)
        $comp_val = ($cols[6] eq '+') ? 0 : 1;

        ## build our temporary object:
        # skip start and stop codon entries, as they are pointless for us.
        unless ($cols[2] =~ /_codon/) {

                %tmp = ('low'    => $cols[3],
                        'high'   => $cols[4],
                        'comp'   => $comp_val );

                push @group_members, {%tmp};

        }

        # Adjust the min/max coords for the gene  (this will capture the extra three
        # bases from the stop codon, too, so all we'll have to modify later is the
        # terminal exon's last position.
        $min_exon_low = ($min_exon_low eq '') ? $cols[3] :
                        (($min_exon_low < $cols[3] ) ? $min_exon_low : $cols[3] );
        $max_exon_high = ($max_exon_high eq '') ? $cols[4] :
                         (($max_exon_high > $cols[4]) ? $max_exon_high : $cols[4] );

        ## If this is the last feat in the file, we should add this gene.
        if ($last_feat_found) {
            &add_group(\@group_members, $min_exon_low, $max_exon_high,
                             $comp_val, \$genes, $source_seq_name);
        }

    }

    return $genes;

}

sub add_group {

    my $arr_ref = shift;
    my @grp_mems = @{$arr_ref};

    my ($min, $max, $comp, $genes, $source_seq_name) = @_;

    ## "Add" the stop codon by increasing the terminal exon 3 bases
    if ($comp) {
        ${$grp_mems[0]}{'low'} -= 3;
    } else {
        ${$grp_mems[-1]}{'high'} += 3;
    }

    ## Create a new gene object:
    my $tmpGene = new Chado::Gene ( $idMaker->next_id( 'type' => 'gene', 'project' => $project ),
                                    $min, $max, $comp, $source_seq_name );

    ## Add the polypeptide/transcript:
    foreach my $type( qw( transcript polypeptide ) ) {

        $tmpGene->addFeature( $idMaker->next_id( 'type' => $type, 'project' => $project ),
                                    $min, $max, $comp, $type );

    }

    ## Now add the exons
    foreach my $feat (@grp_mems) {
        foreach my $type( qw( exon CDS ) ) {
        
            $tmpGene->addFeature( $idMaker->next_id( 'type' => $type, 'project' => $project),
                                $$feat{'low'}, $$feat{'high'}, $$feat{'comp'}, $type );

        }

    }

    ## Form the group now:
    my $count = $tmpGene->addToGroup($tmpGene->getId, {'all' => 1});
    $logger->logdie("Nothing added to group") unless ($count);

    push @{$$genes}, $tmpGene;

}

sub generateBsml {
    my $data = shift;
 
    #Create the document
    my $doc = new BSML::GenePredictionBsml( 'twinscan', $sourcename, $programversion);

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

    # Check for input file(s)
    if($options{'input_list'}) {
        $error .= "Option input_list ($options{'input_list'}) does not exist\n"
            unless(-e $options{'input_list'});
        open(IN, "< $options{'input_list'}") || &_die("Unable to open $options{'input_list'} ($!)");
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

    # make sure we got a project
    unless ($options{'project'}) {
        $error .= "Option project is required.\n";
    } else {
        $project = $options{'project'};
    }

    # Check for output file
    unless($options{'output'}) {
        $error .= "Option output is required.\n";
    } else {
        $output = $options{'output'};
    }

    # Set up the id generator
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
   
    # Check for input fasta sequence:
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

    ## get programversion for the analysis section
    if ($options{'programversion'}) {
        $programversion = $options{'programversion'};
    } else {
        $programversion = 'current';
    }

    if($options{'debug'}) {
        $debug = $options{'debug'};
    }

    unless($error eq "") {
        $logger->logdie($error);
    }

}

