#!/usr/local/bin/perl

eval 'exec /local/packages/perl-5.8.8/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
use lib (@INC,$ENV{"PERL_MOD_DIR"});
#use lib (@INC,"/usr/local/devel/ANNOTATION/ard/chado-v1r12b1-jinman/lib/5.8.8/");
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

genewise2bsml.pl

=head1 SYNOPSIS

USAGE: genewise2bsml.pl
             --input_file|-i
             --input_seq|-s
             --output|-o
             --project|-p
             --id_repository|r
             --help|-h
             --log|-l

=head1 OPTIONS

B<--input_file,-i>
    Input is a genewise raw file.

B<--input_seq,-s>
    File used as input to genewise and containing the genomic fasta sequence.

B<--output,-o>
    Output BSML file name.

B<--project,-p>
    Project/database name used to create feature identifiers.

B<--id_repository,-r>
    Path to the id repository for the named project.

B<--help,-h>
    This help documentation

B<--log,-l>
    Path to intended log file.

=head1 DESCRIPTION

This script parses GFF output from genewise and writes it out 
as BSML suitable for import into CHADO/legacy DBs.  Specifically,
this 2bsml script is meant to accompany the genewise executable in
the genewise_best_loc pipeline.  It expects the input fasta to contain
the coordinates mapping it back to the larger genome molecule, as such is the
way genewise_best_loc works.  It splits the genomic data into regions around previously
described aat search hits, then uses those proteins and the respective regions
they hit on as the input for genewise.  Because of this, we must stitch
the hits back into the overall assemblies by adjusting the coords
to compensate for the truncation of the sequence by the script
'prepare_for_genwise_best_loc.pl'.  


=head1 INPUT

Current assumption is that input files are syntactically correct with respect 
to the GFF spec, (at least as far as that produced by genewise with -gff goes) 
and that they are encoding predicted gene features.

=head1 OUTPUT

Output will be BSML encoding the geneid predicted gene models.  By the time
the data is in bsml, it is 'standalone' and no more tricks need be played
with it.  It can be loaded by a standard loader like bsml2legacydb or bsml2chado.

=head1 CONTACT

    Jason Inman 
    jinman@jcvi.org

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Data::Dumper;
use Chado::Gene;
use BSML::GenePredictionBsml;
use BSML::BsmlBuilder;
use Ergatis::Logger;
use Ergatis::IdGenerator;

my $input;
my $project;
my $output;
my $idMaker;
my $bsml;
my $data;
my $input_seq;
my $seq_id;

my %options = ();
my $results = GetOptions ( \%options,
                            'input_file|i=s',
                            'output|o=s',
                            'project|p=s',
                            'id_repository|r=s',
                            'input_seq|s=s', 
                            'log|l=s',
                            'debug=s',
                            'help|h',
                          ) || pod2usage();

&pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if ($options{'help'});

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

&check_params(\%options);

$data = &parseGenewiseData($input);
$bsml = &generateBsml($data);

$bsml->writeBsml($output);

exit(0);

#### Subroutines ###

sub parseGenewiseData {

    my $file  = shift;

    ## Get orientation:
    my ($strand, $lend, $chunk_len) = &get_info_from_src;

    open(IN, "<$file") or $logger->logdie("Unable to open raw genewise output: $file: $!");

    my $genes;

    my $section = 0;
    my $match_id = '';
    my $tmp;

    while (<IN>) {
    
        # Skip ahead until we reach the section with the evidence lines
        if (/^\/\//) {
            $section++;
            next;
        }
        next if ($section < 2);

        # Now, the fun begins.
        chomp;
        my @fields = split(/\t/,$_);

        # Get the id and score
        my $type_id = $fields[2];
        my $score = $fields[5];

        # Get coords settled.  Start by getting the start and end.
        my $start  = $strand ? $fields[4] : $fields[3];
        my $end    = $strand ? $fields[3] : $fields[4];

        # If we're reversed, switch the distance between the front and back
        # of the source sequence to make adjusting for the offset possible.
        if ($strand) {
            $start = $chunk_len - $start + 1;
            $end   = $chunk_len - $end   + 1;
        }

        ## Next, compensate for the offset from the beginning of the
        #  actual assembly sequence
        $start = $lend + $start - 1;
        $end   = $lend + $end   - 1;

        # Finally, convert to base 0.
        $start--;

        # If this is a new 'match', we'll start a new 'gene'.        
        if ($match_id ne $fields[8]) {

            # If this isn't our first time through here, add what is in $tmp to
            # the array in $genes
            if ($match_id) {

                my $count = $tmp->addToGroup($tmp->getId, { 'all' => 1 });
                $logger->logdie("Nothing was added to group") unless($count);

                push(@{$genes}, $tmp);
                undef $tmp;

            }

            $match_id = $fields[8];

            $logger->logdie("Unexpected start of new prediction at line $.") 
                unless $type_id eq 'match';
            $tmp = new Chado::Gene( $idMaker->next_id( 'type' => 'gene',
                                                       'project' => $project ),
                                    $start, $end, $strand, $seq_id, $score );

            # Add transcript and polypeptide for the match
            foreach my $type ('transcript', 'polypeptide') {
                my $featid = $idMaker->next_id ('type' => $type, 'project' => $project );
                $tmp->addFeature($featid, $start, $end, $strand, $type);
            }

        } elsif ($type_id eq 'cds') {

            # if this is a 'cds' we'll add a 'cds' and 'exon'
            foreach my $type ('exon', 'cds') {
                my $featid = $idMaker->next_id ('type' => $type, 'project' => $project );
                $tmp->addFeature($featid, $start, $end, $strand, $type);
            }

        } elsif ($type_id eq 'intron') {
            next;
        } else {
            $logger->warn("Found unrecognized feature type. Skipping: $type_id at line $.");
        }

    }
    close IN;

    # Catch the last gene...
    if ($match_id) {

        my $count = $tmp->addToGroup($tmp->getId, { 'all' => 1 });
        $logger->logdie("Nothing was added to group") unless($count);

        push(@{$genes}, $tmp);

    }

    # and send them on the way!
    return $genes;

}

sub generateBsml {

    my $data = shift;

    #Create the document
    my $doc = new BSML::GenePredictionBsml( 'genewise', $input_seq );

    foreach my $gene(@{$data}) {
        $doc->addGene($gene);
    }

    my $addedTo = $doc->setFasta($seq_id, $input_seq);
    $logger->logdie("$seq_id was not a sequence associated with the gene") unless($addedTo);

    return $doc;

}

sub check_params {
# Make sure we have what we expect to get

    my $error = '';

    if ($options{'input_file'}) {
        $input = $options{'input_file'};
    } else {
        $error .= "input genewise raw was not defined with --input\n";
    }

    if ($options{'input_seq'}) {

        $input_seq = $options{'input_seq'};
        open(IN, "< $input_seq") or $logger->logdie("Unable to open $input_seq");
        while(<IN>) {
            #assume it's a single fasta file
            if(/^>([^\s+]+)/) {
                $seq_id = $1;
                last;
            }
        }
        close(IN); 

    } else {
        $error .= "input sequence was not defined with --input_seq\n";
    }

    if ($options{'output'}) {
        $output = $options{'output'};
    } else {
        $error .= "output BSML filename was not defined with --output\n";
    }

    if ($options{'project'}) {
        $project = $options{'project'};
    } else {
        $error .= "You must specify a project name with --project\n";
    }

    if ($options{'id_repository'}) {
        $idMaker = new Ergatis::IdGenerator('id_repository' => $options{'id_repository'});
        $idMaker->set_pool_size( 'exon'        => 20,
                                 'transcript'  => 20,
                                 'gene'        => 20,
                                 'polypeptide' => 20,
                                 'CDS'         => 20 );
    } else {
        $error .= "You must specify an id repository with --id_repository\n";
    }

    if ($error) {
        $logger->logdie($error);
    }

}

sub get_info_from_src {
# pull the coords of the genome source and determine if these are on the forward
# or reverse strands.

    my $ori;    

    # open the input sequence
    open (my $src,"<$input_seq") ||
            $logger->logdie("Unable to open genomic source: $input_seq: $!");

    # yank out the defline
    my $defline;
    while (<$src>) {
        if (/^>(\S+)/) {
            $defline = $1;
            last;
        }
        
    }
    close $src;

    # We need to have a defline.
    die "No defline found in input sequence: $input_seq\n";
    
    # pull out the coords and determine the orientation, also determine the
    # 'chunk_len' for adjusting coords on the reverse
    my ($lend, $rend, $chunk_len);
    if ($defline =~ /\.?\d+\.(\d+)\.(\d+)$/) {

        my ($end5, $end3) = ($1, $2);

        die "Can't get coords from defline: $defline\n" unless ($end5 && $end3);
        $ori = ($end5 < $end3) ? 0 : 1;

        ($lend, $rend) = sort {$a <=> $b} ($end5, $end3);
        $chunk_len = $rend - $lend + 1; 

    }

    return $ori, $lend, $chunk_len;

}
