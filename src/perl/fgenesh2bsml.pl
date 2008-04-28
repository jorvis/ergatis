#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

fgenesh2bsml.pl - convert fgenesh output to BSML

=head1 SYNOPSIS

USAGE: fgenesh2bsml.pl 
        --input_file=/path/to/fgenesh.raw.file 
        --output=/path/to/output.bsml
        --project=aa1 
        --fasta_input=/path/to/somefile.fsa 
        --id_repository=/path/to/repository
        --sourcename=sourcename
        --programversion='current'

=head1 OPTIONS

B<--input_file,-f> 
    Input file file from a fgenesh run.  -i, --input_list, will take in a list
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

This script is used to convert the output from a fgenesh search into BSML.

=head1 INPUT

You define the input file using the --input_file option.  This file does not need any
special file extension.  The regular prediction output of fgenesh looks like this:
(Excluding header lines)

=over

      Gn S   Type   Start       End   Score        ORF           Len
      -- -   ----   -----       ---   -----        ---           ---
       1 +   TSS    19447             -7.15  
       1 +   CDSf   19541 -   19632   16.12   19541 -   19630     90
       1 +   CDSi   19755 -   19977   14.12   19756 -   19977    222
       1 +   CDSl   20833 -   20961    2.99   20833 -   20961    129
       1 +   PolA   21055              1.05  

       2 +   TSS    34437             -7.15  
       2 +   CDSf   34531 -   34622   15.25   34531 -   34620     90
       2 +   CDSi   34745 -   34967   20.74   34746 -   34967    222
       2 +   CDSl   35854 -   35982    5.59   35854 -   35982    129
       2 +   PolA   36043              1.05  

       3 +   TSS    39373             -7.15  
       3 +   CDSf   39467 -   39558   15.25   39467 -   39556     90
       3 +   CDSi   39681 -   39903   20.74   39682 -   39903    222
       3 +   CDSl   40770 -   40898    5.74   40770 -   40898    129
       3 +   PolA   40959              1.05  

       4 +   TSS    44415             -8.75  
       4 +   CDSf   45995 -   46151   16.01   45995 -   46150    156
       4 +   CDSl   46997 -   47100    2.71   46999 -   47100    102
       4 +   PolA   47243              1.05  

       5 +   TSS    54703             -4.45  
       5 +   CDSf   54790 -   54881   13.41   54790 -   54879     90
       5 +   CDSi   55010 -   55232   14.20   55011 -   55232    222
       5 +   CDSl   56131 -   56259    3.87   56131 -   56259    129
       5 +   PolA   56365              1.05  

       6 +   TSS    62100             -6.65  
       6 +   CDSf   62187 -   62278   13.59   62187 -   62276     90
       6 +   CDSi   62409 -   62631   19.50   62410 -   62631    222
       6 +   CDSl   63482 -   63610   10.23   63482 -   63610    129
       6 +   PolA   63718              1.05  

       7 +   TSS    68088             -9.45  
       7 +   CDSo   68183 -   68428   14.87   68183 -   68428    246
       7 +   PolA   68509              1.05  

    Predicted protein(s):
    >ID  1   3 exon (s)  19541  -  20961    147 aa, chain +
    MVHFTAEEKAAVTSLWSKMNVEEAGGEALGRLLVVYPWTQRFFDSFGNLSSPSAILGNPK
    VKAHGKKVLTSFGDAIKNMDNLKPAFAKLSELHCDKLHVDPENFKLLGNVMVIILATHFG
    KEFTPEVQAAWQKLVSAVAIALAHKYH
    >ID  2   3 exon (s)  34531  -  35982    147 aa, chain +
    MGHFTEEDKATITSLWGKVNVEDAGGETLGRLLVVYPWTQRFFDSFGNLSSASAIMGNPK
    VKAHGKKVLTSLGDAIKHLDDLKGTFAQLSELHCDKLHVDPENFKLLGNVLVTVLAIHFG
    KEFTPEVQASWQKMVTGVASALSSRYH
    >ID  3   3 exon (s)  39467  -  40898    147 aa, chain +
    MGHFTEEDKATITSLWGKVNVEDAGGETLGRLLVVYPWTQRFFDSFGNLSSASAIMGNPK
    VKAHGKKVLTSLGDAIKHLDDLKGTFAQLSELHCDKLHVDPENFKLLGNVLVTVLAIHFG
    KEFTPEVQASWQKMVTAVASALSSRYH
    >ID  4   2 exon (s)  45995  -  47100     86 aa, chain +
    MGNPKVKAHGKKVLISFGKAVMLTDDLKGTFATLSDLHCNKLHVDPENFLVSTLRQRDID
    CFGNPLQRGFYPTDTGFLAVTNKCCG
    >ID  5   3 exon (s)  54790  -  56259    147 aa, chain +
    MVHLTPEEKTAVNALWGKVNVDAVGGEALGRLLVVYPWTQRFFESFGDLSSPDAVMGNPK
    VKAHGKKVLGAFSDGLAHLDNLKGTFSQLSELHCDKLHVDPENFRLLGNVLVCVLARNFG
    KEFTPQMQAAYQKVVAGVANALAHKYH
    >ID  6   3 exon (s)  62187  -  63610    147 aa, chain +
    MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK
    VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFG
    KEFTPPVQAAYQKVVAGVANALAHKYH
    >ID  7   1 exon (s)  68183  -  68428     81 aa, chain +
    MEQSWAENDFDELREEGFRRSNYSKLKEEVRTNGKEVKNFEKKLDEWITRITNAQKSLKD
    LMELKTKAGELRDKYTSLSNR

=back

The prediction output has the general form of:
(modified from http://linux1.softberry.com/berry.phtml?topic=fgenesh&group=help&subgroup=gfind
 to make sense)

    Gn    - predicted gene number, starting from start of sequence;
    S     - DNA strand (+ for direct or - for complementary);
    Type  - type of coding sequence:
          CDSf - First (Starting with Start codon),
          CDSi - internal (internal exon),
          CDSl - last coding segment, ending with stop codon),
          TSS  - Position of transcription start (TATA-box position and score);
          PolA - Position of Polyadenylation site
    Start - Start position of the Feature;
    End   - End Position of the Feature;
    Score - score for this Feature;
    ORF   - start/end positions where the first complete codon starts and the last codon ends.

=head1 OUTPUT

Base positions from the input file are renumbered so that positions start at zero and
reflect interbase numbering.  

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

    $data = &parseFgeneshData($file);
    $bsml = &generateBsml($data);

}
$bsml->writeBsml($output);

exit(0);


sub parseFgeneshData {

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

    while (<$ifh>) {
        if ($_ =~ /Seq name: (\S+)/) {
            $source_seq_name = $1;
            last;
        }
    }

    ## skip the rest of the headers
    while (<$ifh>) {
        last if ($_ =~ /^   G Str   Feature   Start        End    Score           ORF           Len/); 
    }
    <$ifh>;

    ## go through the data now
    while (<$ifh>) {

        ## We're done if we encounter this line:
        last if ($_ =~ /^Predicted protein/);

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
        ## Because we don't get start and stop at every line we have to be a little tricky
        my ($strand, $type, $low, $high);
        if (/TSS|PolA/) {

            ($strand, $type, $low) = (split(/\s+/))[1,2,3];

            ## have to handle the $high since it's not inlcuded (results are base-indexed)
            ## We'll handle conversion to gap-indexing (interbase) in a little bit.
            $high = $low;

        } elsif (/CDS/) {

            ($strand, $type, $low, $high) = (split(/\s+/))[1,3,4,6];

        }

        ##store complement flag (1 is true, meaining feature lies on reverse strand)
        $comp_val = ($strand eq '+') ? 0 : 1;

        ## count in interbase (Note that fgenesh reports all coords as l_end - r_end
        $low--;

        ## build our temporary object for this element
        %tmp = ('type'   => $type,
                'low'    => $low,
                'high'   => $high,
                'comp'   => $comp_val );
        push @group_members, {%tmp};

        # Adjust the min/max coords for the gene and polypeptide
        unless ($type =~ /TSS|PolA/) {

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

        ## All kinds of CDS (CDSf, CDSi, CDSl, CDSo) get handled the same way.
        if ($mem_type =~ /CDS/) {

            foreach my $type ( qw( exon CDS ) ) {
        
                $tmpGene->addFeature( $idMaker->next_id( 'type' => $type, 'project' => $project),
                                      $$feat{'low'}, $$feat{'high'}, $$feat{'comp'}, $type );

            }

        } elsif ($mem_type eq 'TSS') {

            $tmpGene->addFeature( $idMaker->next_id( 
                                        'type' => 'transcription_start_site',
                                        'project' => $project
                                        ),
                                  $$feat{'low'},
                                  $$feat{'high'},
                                  $$feat{'comp'},
                                  'transcription_start_site' );

        } elsif ($mem_type eq 'PolA') {

            $tmpGene->addFeature( $idMaker->next_id(
                                         'type' => 'polyA_signal_sequence', 
                                         'project' => $project),
                                  $$feat{'low'},
                                  $$feat{'high'},
                                  $$feat{'comp'},
                                  'polyA_signal_sequence' );

        } else {
            ## we need to at least consider each feature type in the gene object.
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
    my $doc = new BSML::GenePredictionBsml( 'fgenesh', $sourcename, $programversion);

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

