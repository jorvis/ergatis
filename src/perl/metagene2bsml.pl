#!/usr/local/bin/perl

use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME

metagene2bsml.pl  - converts raw metagene output into BSML

=head1 SYNOPSIS

USAGE:  metagene2bsml.pl -i raw_metagene_output.txt -f /path/to/fasta_input.fsa -o /path/to/outfile.bsml [-l /path/to/logfile.log] --project project_name --id_repository /path/to/id_repository

=head1 OPTIONS

B<--input_file,-i>
    The input raw output from metagene.

B<--fasta_input,-f>
    The full path to the fasta sequence file that was input to metagene.
    
B<--output,-o>
    The full path to the BSML output file to be created. 

B<--debug,-d>
    Debug level.  Use a large number to turn on verbose debugging.

B<--log,-l>
    Log file

B<--help,-h>
    This help message

=head1   DESCRIPTION

Create BSML document from metagene output. Will extract nucleotide sequences for the
predicted ORFs.

=head1  CONTACT
        Brett Whitty
        bwhitty@tigr.org

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Ergatis::Logger;
use Ergatis::IdGenerator;
use XML::Twig;
use Data::Dumper;
use BSML::BsmlBuilder;
use File::Basename;
use BSML::GenePredictionBsml;

my %options = ();
my $results = GetOptions (\%options,
                          'input_file|i=s',
                          'fasta_input|f=s',
                          'output|o=s',
                          'id_repository=s',
                          'analysis_id=s',
                          'project=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();
                     
if ($options{'help'}) {
    pod2usage(verbose => 2);
}

if (! $options{'analysis_id'}) {
    $options{'analysis_id'} = 'metagene_analysis';
}

if (!$options{'id_repository'}) {
    pod2usage("must provided --id_repository");
}
if (!$options{'project'}) {
    pod2usage("project name must be provided with --project");
}
my $id_gen = Ergatis::IdGenerator->new('id_repository' => $options{'id_repository'});

unless ($options{'input_file'}) {
    pod2usage("raw metagene output must be provided with --input_file");
}

unless ($options{'output'}) {
    pod2usage("must specify an output bsml file with --output");
}
my $output_dir = dirname($options{'output'});
if ($output_dir eq '') {
    $output_dir = '.';
}
my $output_bsml_file = $options{'output'};
my $output_fsa_file = $output_dir."/".basename($options{'output'}, ".bsml").".fsa";

my $in_fh = open_file_read($options{'input_file'});

my $orfs = {};

my $doc = new BSML::BsmlBuilder();

my $orf_count = 0;
while (<$in_fh>) {
    chomp;
    if (/^# (.*)/) {
        my $defline = $1;
        $defline =~ /^(\S+)/;
        my $id = $1;
        
        ## read % gc line
        my $temp = <$in_fh>;
        chomp $temp;
        $temp =~ /^# gc = (.*)/ || die "Failed parsing gc line";
        my $gc = $1;
        
        ## read model line
        $temp = <$in_fh>;
        $temp =~ /^# (.*)/ || die "Failed parsing model prediction";
        my $model = $1;
        
        ## read predicted orfs
        while (<$in_fh>) {
            chomp;
            if ($_ eq '') {
                last;
            } else {
                my (
                    $startpos, 
                    $endpos, 
                    $strand, 
                    $frame, 
                    $score, 
                    $state,
                   ) = split("\t", $_);
                
                ## set flags for partial ORFs
                my ($five_prime_partial, $three_prime_partial);
                if ($state eq 'complete') {
                    ($five_prime_partial, $three_prime_partial) = (0, 0);
                } elsif ($state eq 'partial (lack both ends)') {
                    ($five_prime_partial, $three_prime_partial) = (1, 1);
                } elsif ($state eq 'partial (lack 5\'-end)') {
                    ($five_prime_partial, $three_prime_partial) = (1, 0);
                } elsif ($state eq 'partial (lack 3\'-end)') {
                    ($five_prime_partial, $three_prime_partial) = (0, 1);
                }
               
                push (@{$orfs->{$id}},  {    
                     #    'defline'   =>  $defline, ## not going to store the defline to save memory
                                                    ## grab it from the fasta file later
                    'gc'          =>  $gc,
                    'model'       =>  $model,
                    'startpos'    =>  $startpos,
                    'endpos'      =>  $endpos,
                    'complement'  =>  ($strand eq '+') ? 0 : 1,
                    'frame'       =>  $frame,
                    'score'       =>  $score,
                    '5_partial'   =>  $five_prime_partial,
                    '3_partial'   =>  $three_prime_partial,
                               });
                $orf_count++;                  
#                print "$gc\t$model\t$state\n";
            }
        }
    }
}

$id_gen->set_pool_size('ORF' => $orf_count, 'CDS' => $orf_count);

my $seq_refs = {};
my $parent_seq = {};

my $out_fsa_fh = open_file_write($output_fsa_file);

foreach my $seq_id (keys(%{$orfs})) {
    my $seq = $doc->createAndAddExtendedSequenceN(  
                                                'id' => $seq_id,
                                                'identifier' => $seq_id,
                                                'molecule' => 'na',
                                                'class' => 'assembly',
                                                 );
    $seq_refs->{$seq_id} = $seq;                                                 
    $doc->createAndAddSeqDataImportN(  
                                        'seq'           => $seq,
                                        'format'        => 'bsml',
                                        'source'        => $options{'fasta_input'},
                                        'id'            => '',
                                        'identifier'    => $seq_id,
                                    );
    $seq->addBsmlLink('analysis', '#'.$options{'analysis_id'}, 'input_of');

    my $ft = $doc->createAndAddFeatureTable($seq);

    foreach my $orf (@{$orfs->{$seq_id}}) {
                
        my $orf_id = $id_gen->next_id(
                        'project' => $options{'project'}, 
                        'type'    => 'ORF'
                                     );

        $orf->{'id'} = $orf_id;

        $parent_seq->{$orf_id} = $seq_id;        
                                     
        my $feat = $doc->createAndAddFeature(
                                                $ft, 
                                                $orf_id, 
                                                $orf_id, 
                                                'ORF'
                                            );
        $doc->createAndAddIntervalLoc(
                                        $feat,
                                        $orf->{'startpos'}, # start
                                        $orf->{'endpos'}, # end
                                        $orf->{'complement'}, # complement
                                        $orf->{'5_partial'}, # startopen
                                        $orf->{'3_partial'}, # endopen
                                                );

        ## This should be added to the interval-loc instead of the feature
        $doc->createAndAddBsmlAttribute(
                                            $feat, 
                                            'phase', 
                                            $orf->{'frame'}
                                       );
        $doc->createAndAddBsmlAttribute(
                                            $feat, 
                                            'score', 
                                            $orf->{'score'}
                                       );
        $feat->addBsmlLink('analysis', '#'.$options{'analysis_id'}, 'computed_by');
        $feat->addBsmlLink('sequence', "#${orf_id}_seq");
        
        my $orf_seq = $doc->createAndAddExtendedSequenceN(  
                                    'id' => $orf_id."_seq",
                                    'molecule' => 'na',
                                    'class' => 'ORF',
                                                         );
        $seq_refs->{$orf_id} = $orf_seq;                                                         
        $doc->createAndAddSeqDataImportN(   
                                    'seq'           => $orf_seq,
                                    'format'        => 'fasta',
                                    'source'        => $output_fsa_file,
                                    'id'            => '',
                                    'identifier'    => $orf_id,
                                        );
                                         
    }
} 

## add seq-data imports to sequence elements
my $fasta_fh = open_file_read($options{'fasta_input'});

my $sequence = '';
my $sequence_defline = '';
my $sequence_id = '';
while (<$fasta_fh>) {
    chomp;
    if (/^#/) { next;}
    if (/^>((\S+).*)/) {
        my $current_id = $2;
        my $current_defline = $1;
       
        if ($sequence ne '') {
            add_sequence_to_bsml($sequence_id, $sequence_defline, $out_fsa_fh);
        }
        $sequence_defline = $current_defline;
        $sequence_id = $current_id;
        $sequence = '';
    } else {
        $sequence .= $_;
    }
}
add_sequence_to_bsml($sequence_id, $sequence_defline, $out_fsa_fh);

$doc->write($output_bsml_file, '', 1);

## adds input sequences (if not already added), deflines, and predicted orf sequences
sub add_sequence_to_bsml {
    my ($sequence_id, $sequence_defline, $out_fsa_fh) = @_;
    
        my $seq_ref;
    if ($seq_refs->{$sequence_id}) {
        $seq_ref = $seq_refs->{$sequence_id};
    } else {
        $seq_ref = $doc->createAndAddExtendedSequenceN(  
                                    'id' => $sequence_id,
                                    'identifier' => $sequence_id,
                                    'molecule' => 'na',
                                    'class' => 'assembly',
                                                      );
        $doc->createAndAddSeqDataImportN(  
                                    'seq'           => $seq_ref,
                                    'format'        => 'fasta',
                                    'source'        => $options{'fasta_input'},
                                    'id'            => '',
                                    'identifier'    => $sequence_id,
                                        );
        $seq_ref->addBsmlLink('analysis', '#'.$options{'analysis_id'}, 'input_of');
    }
    $doc->createAndAddBsmlAttribute(
                                    $seq_ref, 
                                    'defline', 
                                    $sequence_defline,
                                   );
    ## write ORF sequences                                       
    foreach my $orf(@{$orfs->{$sequence_id}}) {
        extract_and_write_sequence($out_fsa_fh, \$sequence, $orf);
    }
}

## extracts an ORF from its parent sequence and writes it to a fasta file
sub extract_and_write_sequence {
    my ($out_fsa_fh, $seq_ref, $orf) = @_;

    my $seq = substr(${$seq_ref}, $orf->{'startpos'}, $orf->{'endpos'} - $orf->{'startpos'});
  
    ## reverse complement orf if on complement 
    if ($orf->{'complement'} == 1) {
        $seq = reverse_complement_dna($seq); 
    }
  
    ## remove untranslated bases from the 5' end 
    $seq = substr($seq, $orf->{'frame'});
    
    my $outfile = $output_dir."/".$orf->{'id'}.".fsa";
    
    write_seq_to_fasta($out_fsa_fh, $orf->{'id'}, $seq); 
}

## writes sequence to a fasta file
sub write_seq_to_fasta {
    my ($out_fsa_fh, $header, $sequence) = @_;
    
    $sequence =~ s/\W+//g;
    $sequence =~ s/(.{1,60})/$1\n/g;
        
    print $out_fsa_fh ">$header\n$sequence";
}

## opens a filehandle for writing
sub open_file_write {
    my ($filename, $gzip_mode) = @_;

    my $fh;

    if ($gzip_mode) {
        open ($fh, ">:gzip", $filename.'.gz') || die "Couldn't open '$filename' for writin
g: $!";
    } else {
        open ($fh, ">$filename") || die "Couldn't open '$filename' for writing: $!";
    }

    return $fh
}

## opens a filehandle for reading
sub open_file_read {
    my ($filename) = @_;

    my $fh;

    if (! -e $filename && -e $filename.'.gz') {
        $filename .= '.gz';
    }
    if ($filename =~ /\.gz$|\.gzip$/) {
        open ($fh, "<:gzip", $filename) || die "Couldn't open '$filename' for reading: $!"
;
    } else {
        open ($fh, "<$filename") || die "Couldn't open '$filename' for reading: $!";
    }

    return $fh
}

## reverse complement orfs
sub reverse_complement_dna {
        my ($r_seq) = @_;
        $r_seq =~ tr/AaCcGgTtMmRrWwSsYyKkVvHhDdBb/TtGgCcAaKkYyWwSsRrMmBbDdHhVv/;
        $r_seq = reverse($r_seq);
        return $r_seq;
}
