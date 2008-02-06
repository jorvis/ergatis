#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

genscan2bsml.pl - convert genscan output to BSML

=head1 SYNOPSIS

USAGE: genscan2bsml.pl 
        --input_file=/path/to/genscan.raw.file 
        --output=/path/to/output.bsml
        --project=aa1 
        --fasta_input=/path/to/somefile.fsa 
        --id_repository=/path/to/repository
        --sourcename=sourcename
        --programversion='current'

=head1 OPTIONS

B<--input_file,-f> 
    Input file file from a genscan run.  -i, --input_list, will take in a list
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

This script is used to convert the output from a genscan search into BSML.

=head1 INPUT

You define the input file using the --input_file option.  This file does not need any
special file extension.  The regular prediction output of genscan looks like this:
(Excluding header lines)

Gn.Ex Type S .Begin ...End .Len Fr Ph I/Ac Do/T CodRg P.... Tscr..
----- ---- - ------ ------ ---- -- -- ---- ---- ----- ----- ------

 1.01 Term +   5651   6289  639  1  0   81   35   400 0.481  34.31
 1.02 PlyA +   6539   6544    6                              -2.62

 2.00 Prom +   6570   6609   40                              -6.39
 2.01 Init +   8942   9175  234  1  0   66   57   714 0.999  67.88
 2.02 Intr +   9209  11608 2400  1  0    9    6  2378 0.397 225.56
 2.03 Intr +  11642  12440  799  1  1   50  -60   947 0.007  74.06
 2.04 Intr +  12604  12679   76  2  1   43   87    13 0.006  -2.13
 2.05 Term +  13098  13374  277  0  1   51   46   232 0.275  14.16
 2.06 PlyA +  13643  13648    6                               0.77

which has the general format shown below.  (Taken from
 http://genome.imim.es/courses/Bioinformatics2003_genefinding/results/genscan.html)

Gn.Ex : gene number, exon number (for reference)
Type  : Init = Initial exon (ATG to 5' splice site)
        Intr = Internal exon (3' splice site to 5' splice site)
        Term = Terminal exon (3' splice site to stop codon)
        Sngl = Single-exon gene (ATG to stop)
        Prom = Promoter (TATA box / initation site)
        PlyA = poly-A signal (consensus: AATAAA)
S     : DNA strand (+ = input strand; - = opposite strand)
Begin : beginning of exon or signal (numbered on input strand)
End   : end point of exon or signal (numbered on input strand)
Len   : length of exon or signal (bp)
Fr    : reading frame (a forward strand codon ending at x has frame x mod 3)
Ph    : net phase of exon (exon length modulo 3)
I/Ac  : initiation signal or 3' splice site score (tenth bit units)
Do/T  : 5' splice site or termination signal score (tenth bit units)
CodRg : coding region score (tenth bit units)
P     : probability of exon (sum over all parses containing exon)
Tscr  : exon score (depends on length, I/Ac, Do/T and CodRg scores)

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

    $data = &parseGenscanData($file);
    $bsml = &generateBsml($data);

}
$bsml->writeBsml($output);

exit(0);


sub parseGenscanData {

    my $inFile = shift;

    ## open the input file for parsing
    open (my $ifh, "< $inFile") || $logger->logdie("Can't open input file for reading");

    my $source_seq_name = '';
    my $group_name = '';
    my $comp_val;
    my $min_exon_low = '';
    my $max_exon_high = '';
    my $genes;
    my %tmp;
    my @group_members = ();
print "getting into the file\n";
    while (<$ifh>) {
        print $_;
        if ($_ =~ /^Sequence (\S+) :/) {
            $source_seq_name = $1;
            last;
        }
    }

    ## skip the rest of the headers
    while (<$ifh>) {
        last if ($_ =~ /^[- ]+\n/); 
    }
    <$ifh>;

print "data points here we come!\n";
    ## go through the data now
    while (<$ifh>) {
        ## We're done if we encounter this line:
        last if ($_ =~ /Predicted peptide sequence/);

        chomp;
        ## remove a leading space from the early lines.  This is sometimes
        ## inserted to keep formatting 'pretty'.  Honestly, who's reading the
        ## raw files?  Come ON!
        $_ =~ s/^ +//;

        ## add the group if we come across a blank line
        if ($_ eq '') {
            &add_group(\@group_members, $min_exon_low, $max_exon_high,
                        $comp_val, \$genes, $source_seq_name);
            @group_members = ();
            $min_exon_low = '';
            $max_exon_high = '';
            next;
        }

        ## See description of the columns above.
        my @cols = split(/\s+/);

        ##store complement flag (1 is true, meaining feature lies on reverse strand)
        $comp_val = ($cols[2] eq '+') ? 0 : 1;

        # set our low/high for this feature
        my $low  = ($comp_val) ? $cols[4] : $cols[3];
        my $high = ($comp_val) ? $cols[3] : $cols[4];

        ## count in interbase:
        $low--;

        ## build our temporary object for this element
        %tmp = ('type'   => $cols[1],
                'low'    => $low,
                'high'   => $high,
                'comp'   => $comp_val );
        push @group_members, {%tmp};

        # Adjust the min/max coords for the gene and polypeptide
        unless ($cols[1] =~ /Prom|PlyA/) {

            # set up initial low/highs or compare against existing low/highs
            $min_exon_low = ($min_exon_low eq '') ? $low :
                            (($min_exon_low < $low ) ? $min_exon_low : $low );
            $max_exon_high = ($max_exon_high eq '') ? $high :
                             (($max_exon_high > $high) ? $max_exon_high : $high );

        }


    }

    return $genes;

}

sub add_group {
print "adding gfrp\n";
    my ($grp_mems, $min, $max, $comp, $genes, $source_seq_name) = @_;

    ## Create a new gene object:
    my $tmpGene = new Chado::Gene ( $idMaker->next_id( 'type' => 'gene', 'project' => $project ),
                                    $min, $max, $comp, $source_seq_name );

    ## Add the polypeptide/transcript:
    foreach my $type( qw( transcript polypeptide ) ) {

        $tmpGene->addFeature( $idMaker->next_id( 'type' => $type, 'project' => $project ),
                                    $min, $max, $comp, $type );

    }

    ## Now add the exons and other features
    foreach my $feat (@$grp_mems) {

        my $mem_type = $feat->{'type'};

        if ($mem_type =~ /Init|Intr|Term|Sngl/) {

            foreach my $type ( qw( exon CDS ) ) {
        
                $tmpGene->addFeature( $idMaker->next_id( 'type' => $type, 'project' => $project),
                                      $$feat{'low'}, $$feat{'high'}, $$feat{'comp'}, $type );

            }

        } elsif ($mem_type eq 'Prom') {

            $tmpGene->addFeature( $idMaker->next_id( 'type' => 'promoter', 'project' => $project),
                                  $$feat{'low'}, $$feat{'high'}, $$feat{'comp'}, 'promoter' );

        } elsif ($mem_type eq 'PlyA') {

            $tmpGene->addFeature( $idMaker->next_id( 'type' => 'polyA_signal_sequence', 'project' => $project),
                                  $$feat{'low'}, $$feat{'high'}, $$feat{'comp'}, 'polyA_signal_sequence' );

        } else {
            $logger->logdie("unrecognized feature type: $mem_type\n");
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
    my $doc = new BSML::GenePredictionBsml( 'genscan', $sourcename, $programversion);

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
                                 'CDS'         => 20,
                                 'promoter'    => 20,
                                 'polyA_signal_sequence' => 20,
                                 );

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

