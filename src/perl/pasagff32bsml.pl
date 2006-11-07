#!/usr/local/bin/perl

=head1 NAME

pasagff32bsml.pl

=head1 SYNOPSIS

USAGE: pasagff32bsml.pl
           --input|-i
           --output|-o
           --project|-p
           --help|-h

=head1 OPTIONS

B<--input,-i>
    Input gff3 format file created by PASA.

B<--output,-o>
    Output BSML file name.

B<--project,-p>
    Project/database name used to create feature identifiers.

B<--help,-h>
    This help documentation

=head1 DESCRIPTION

This script parses GFF output from PASA and writes it out 
as BSML suitable for import into CHADO DBs.

=head1 INPUT

Input should be GFF3 format output from PASA.

=head1 OUTPUT

Output will be BSML encoding the PASA output.

=head1 CONTACT

    Brett Whitty
    bwhitty@tigr.org

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use BSML::BsmlBuilder;
#use Ergatis::IdGenerator;

my ($input, $help, $man, $output, $project, $id_repository, $analysis, $mapping);

GetOptions (
             'help|h'              => \$help,
             'man|m'               => \$man,
             'output|o=s'          => \$output,
             'input|i=s'           => \$input,
#             'project|p=s'         => \$project,
#             'id_repository=s'     => \$id_repository,
             'analysis=s'          => \$analysis,
             'mapping:s'           => \$mapping,
           ) || pod2usage();

&pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if ($man);
&pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT}) if ($help);


if (!$input){
    pod2usage("input gff_file was not defined with --input");
}

if (!$output){
    pod2usage("output BSML filename was not defined with --output");
}

if (!$project) {
        pod2usage("You must specify a project name with --project");
}

if (!$id_repository) {
    pod2usage("You must specify an id repository with --id_repository");
}

if ($mapping) {
    unless (-e $mapping) {
        pod2usage("Specified mapping file does not exist");
    }
}

$project =~ tr/A-Z/a-z/;

$analysis = 'pasa_'.$analysis;

my $idcreator = Ergatis::IdGenerator->new('id_repository' => $id_repository);

my $global_id_counter=0;
my $nodes = {};
my @root_nodes = ();

my %id_map;
if ($mapping) {
    open (IN, $mapping) || die "couldn't open mapping file for reading";
    while (<IN>) {
        chomp;
        my @t = split ("\t", $_);
        $id_map{$t[0]} = $t[1];
    }
    close IN;
}
open (IN, $input) || die "couldn't open input file '$input'";

while (<IN>) {
    chomp;
    
    if (/^#/) {
        #process header line
    } elsif (/^>/ || /^#.*fasta.*/i) {
        last;
    } else {
        my $record = parse_record($_);
        
        ## correct numbering for interbase (see parse_record)
        $record->{'_end'}++;
        
        ## if it has no parents it's a root node    
        if (!$record->{'Parent'}) {
            push(@root_nodes, $record);
        }
    
        ## store records in a set of arrays indexed by ID
        if ($record->{'ID'}) {
            push(@{$nodes->{$record->{'ID'}}->{'records'}},$record);
        }
    }
    
}
close IN;

my @chains;
foreach my $node_id(keys(%{$nodes})) {
    my $chain;
    my $refstart;
    my $refend;
    foreach my $node_hash(@{$nodes->{$node_id}->{'records'}}) {

        if (!defined($refstart)) {
            $refstart   = $node_hash->{'_start'};
            $refend     = $node_hash->{'_end'};
        }
        if ($refstart > $node_hash->{'_start'}) {
            $refstart = $node_hash->{'_start'};
        }
        if ($refend < $node_hash->{'_end'}) {
            $refend = $node_hash->{'_end'};
        }
        
        my @comp = split(" ", $node_hash->{'Target'}->[0]);
        $comp[1]--; # fix numbering for interbase
        
        if ($mapping) {
            if (!defined($id_map{$node_hash->{'_seqid'}})) {
                die "mapping for refseq id '$node_hash->{_seqid}' was not defined";
            }
            $chain->{'refseq'}          = $id_map{$node_hash->{'_seqid'}};
        } else {
            $chain->{'refseq'}          = $node_hash->{'_seqid'};
        }
        $chain->{'compseq'}         = $comp[0];
        $chain->{'compdatabase'}    = $node_hash->{'_source'};
        $chain->{'compidentifier'}  = $comp[0];

        my $score = $node_hash->{'_score'} eq '.' ? '' : $node_hash->{'_score'};
        
        push (  @{$chain->{'percent_identity'}},
                $score
             );
        
        push (  @{$chain->{'refpos'}},
                $node_hash->{'_start'}
             );
        push (  @{$chain->{'runlength'}},
                $node_hash->{'_end'} - $node_hash->{'_start'}
             );
        push (  @{$chain->{'refcomplement'}},
                $node_hash->{'_strand'}
             );
        
        push (  @{$chain->{'comppos'}},
                $comp[1]
             );
        push (  @{$chain->{'comprunlength'}},
                $comp[2] - $comp[1]
             );
        my $compcomplement = $comp[3] eq '+' ? 0 : 1; 
        push (  @{$chain->{'compcomplement'}},
                $compcomplement
             );
    }
    $chain->{'refstart'}    = $refstart;
    $chain->{'refend'}      = $refend;
    $chain->{'reflength'}   = $refend - $refstart;

    push (@chains, $chain);
}

my $doc = new BSML::BsmlBuilder();
foreach my $chain(@chains) {
    if( !( $doc->returnBsmlSequenceByIDR($chain->{'refseq'}) )){
        my $seq = $doc->createAndAddSequence(
                                        $chain->{'refseq'},  # query sequence
                                        $chain->{'refseq'},  # query sequence
                                        '',                  # length 
                                        'dna',               # mol-type
                                        '',                  # class 
                                            );
        $seq->addBsmlLink('analysis', '#'.$analysis.'_analysis', 'input_of');
    }
    if( !( $doc->returnBsmlSequenceByIDR($chain->{'compseq'}) )){
        my $seq = $doc->createAndAddSequence(
                                        $chain->{'compseq'},    #subject id
                                        $chain->{'compseq'},    #subject id
                                        '',                     #length
                                        'dna',                  #mol-type
                                        '',                     #class
                                                );
        $seq->addBsmlLink('analysis', '#'.$analysis.'_analysis', 'input_of');   
        $doc->createAndAddCrossReference(
                        'parent'          => $seq,
                        'database'        => $chain->{'compdatabase'},
                        'identifier'      => $chain->{'compidentifier'},
#                        'identifier-type' => '',         
                                        );

    }
    my $aln = $doc->createAndAddSequencePairAlignment( 
                                        'refseq'    => $chain->{'refseq'},
                                        'refstart'  => $chain->{'refstart'},
                                        'refend'    => $chain->{'refend'},
                                        'reflength' => $chain->{'reflength'},
                                        'compseq'   => $chain->{'compseq'},
                                        'class'     => 'match',
                                        'method'    => 'PASA',
                                                     );
    $aln->addBsmlLink('analysis', '#'.$analysis.'_analysis', 'computed_by');   
    
    my $spr_count = scalar(@{$chain->{'refpos'}});
    
    foreach (my $i=0; $i < $spr_count; $i++) {
        my $s = $doc->createAndAddSequencePairRun(
                'alignment_pair'    => $aln,
                'refpos'            => $chain->{'refpos'}->[$i],
                'runlength'         => $chain->{'runlength'}->[$i],
                'comprunlength'     => $chain->{'comprunlength'}->[$i],
                'comppos'           => $chain->{'comppos'}->[$i],
                'refcomplement'     => $chain->{'refcomplement'}->[$i],
                'compcomplement'    => $chain->{'compcomplement'}->[$i],
                                                 );
        $doc->createAndAddBsmlAttribute(
                                            $s, 
                                            'class', 
                                            'match_part'
                                       );
        $doc->createAndAddBsmlAttribute(
                                            $s, 
                                            'percent_identity', 
                                            $chain->{'percent_identity'}->[$i]
                                       );
    }
 
}

# add the analysis element
$doc->createAndAddAnalysis(
                            id => $analysis.'_analysis',
                            sourcename => $output,
                          );

$doc->write($output);

exit();

sub parse_record {
    my ($line) = @_;

    my $record = {};
    my $attrib_hash = {};
    my $record_id = '';
    my $parent_id = '';
    
    my @cols = split("\t");

    ## adjust both positions so that we are numbering from zero
    $cols[3]--;
    $cols[4]--;

    ## change the + and - symbols in strand column to 0 and 1, respectively
    if ($cols[6] eq '+') {
        $cols[6] = 0;
    } elsif ($cols[6] eq '-') {
        $cols[6] = 1;
    } else {
        die("unknown value ($cols[6]) in strand column.  expected + or -.");
    }

    $record->{'_seqid'}     = $cols[0];
    $record->{'_source'}    = $cols[1];
    $record->{'_type'}      = $cols[2];
    $record->{'_start'}     = $cols[3];
    $record->{'_end'}       = $cols[4];
    $record->{'_score'}     = $cols[5];
    $record->{'_strand'}    = $cols[6];
    $record->{'_phase'}     = $cols[7];
    
    my @attribs = split(";", $cols[8]);
    foreach my $attrib(@attribs) {
        my ($type, $val) = split("=", $attrib);
        my @vals = split(",", $val);
        $record->{$type}=\@vals;
    }

    if (!defined($record->{'ID'})) {
        #print STDERR "no ID defined for record\n";
    } else {
        $record->{'ID'} = $record->{'ID'}->[0];
    }

    if (!defined($record->{'Parent'})) {
        #print STDERR "no Parent attribute defined for record\n";
    } else {
        foreach my $parent(@{$record->{'Parent'}}) {
            ## store record's hash reference as a child of Parent
            if (!defined($nodes->{$parent}->{'children'})) {
                $nodes->{$parent}->{'children'}->{$record->{'_type'}} = [];
            }
            push (@{$nodes->{$parent}->{'children'}->{$record->{'_type'}}},$record);
        }
    }

    return $record;
}
