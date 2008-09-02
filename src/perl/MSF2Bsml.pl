#!/usr/bin/perl

BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

MSF2Bsml.pl - Convert a multiple sequence alignment file (such as clustalw) to BSML.

=head1 SYNOPSIS

USAGE: MSF2Bsml.pl
            --msffile=/path/to/somefile.clw
            --analysisconf=/path/to/some.config
            --output=/path/to/output.bsml
          [ --log=/path/to/some.log
            --dnd_file/path/to/some.dnd
            --debug=4
          ]

=head1 OPTIONS

B<--msffile,-f>
    Input multiple sequence alignment file, such as the .clw file created by clustalw.

B<--output,-o>
    The full path to the BSML file that will be created.

B<--analysis_conf,-a>
    Optional. The BSML file created contains an Analysis element which eventually record all 
    the parameters in the analysis.  Usually this should point to the Ergatis pipeline.config 
    file for the component/analysis preceeding the clustalw step (such as jaccard or cogs).
    If omitted, the id attribute within the Analysis element will default to 'clustalw_analysis'.

B<--dnd_file,-d>
    Optional.  The full path to a dnd (tree) file corresponding to the given alignment.
    If passed, this will store the newick tree as a string in BSML.

B<--fastafile>
    Optional.  The full path to the multifasta input file containing the sequences that
    were aligned to generate the MSF file. Will be used to populate Seq-data-import in
    Sequence elements if provided.
    
B<--debug> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

Converts an MSF file, created by programs such as clustalw and muscle, into
BSML with the option of including a corresponding newick tree.

=head1  INPUT

An MSF file is the main required input, and has a format like this:

    PileUp



       MSF:  219  Type: P    Check:  8452   .. 

     Name: ntrp01_2_NTORF0805_polypeptide oo  Len:  219  Check:  6977  Weight:  25.1
     Name: ntrt01_1_NTORF0807_polypeptide oo  Len:  219  Check:  4764  Weight:  25.2
     Name: ntrc01_2_NTORF1325_polypeptide oo  Len:  219  Check:  5582  Weight:  24.1
     Name: ntru01_1_NTORF1351_polypeptide oo  Len:  219  Check:  5364  Weight:  25.4
     Name: got_3002_ORFB02050_polypeptide oo  Len:  219  Check:  5765  Weight:  45.4

    //



    ntrp01_2_NTORF0805_polypeptide      MFEKYIMYLK NLIFFQFIVY FFFISLTILI IKNFQQEYSK SILDKQVAQE 
    ntrt01_1_NTORF0807_polypeptide      MFEKYIMYLK NLIFFQFIIY FFFIFLTIWI IKNFQQEYSK SILDKQAAQE 
    ntrc01_2_NTORF1325_polypeptide      MFEKYIMYLK NLIFFQFIIY FFFISLTIWI IKIFQQEYSK SILDKQVLQE 
    ntru01_1_NTORF1351_polypeptide      MFEKYIMYLK NLIFFQFIIY FFFISLTIWI IKNFQQEYSK SILDKQVSQE 
    got_3002_ORFB02050_polypeptide      MINQRILYLK NLICFKIVLY GIIILLCHIY GSYLYAHF.. EILNTKIQKN 


    ntrp01_2_NTORF0805_polypeptide      NLTEEVLKLY SVINSKEEIL ES..YKKYVA LSVPKNSVSC LNYQELIPRI 
    ntrt01_1_NTORF0807_polypeptide      NLTEEVLKLY SVINSKEEIL ES..YKKYVA LSVP.NSVRR FNYQELIPRI 
    ntrc01_2_NTORF1325_polypeptide      NLTEEVLKLY SVINSKEEIL ES..YKKYVD LSVP.SSVSY LNYQELIPKI 
    ntru01_1_NTORF1351_polypeptide      NLTEEVLKLY SVINSKEEIL ES..YKKYVD LSVP.SNVSC LNYQELIPKI 
    got_3002_ORFB02050_polypeptide      KKHLEILQN. .KYNSAIELL NNSALQKKIK PILQQKLSQP YDRNKLIEHY 

It can be created from components such as cogs, clustalw and muscle.  Optional input includes
a DND file containing the distance information gleaned from the alignment above.  Its format
is like:

    (
    got_3002_ORFB02050_polypeptide:0.77227,
    (
    ntrc01_2_NTORF1325_polypeptide:0.01364,
    ntru01_1_NTORF1351_polypeptide:0.01953)
    :0.02066,
    (
    ntrp01_2_NTORF0805_polypeptide:0.01970,
    ntrt01_1_NTORF0807_polypeptide:0.03244)
    :0.02309);

If passed with the --dnd_file option, this string will be included in the output BSML.

=cut


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Config::IniFiles;
use Ergatis::Logger;
use BSML::BsmlBuilder;


my %options = ();
my $results = GetOptions (\%options, 
              'msffile|f=s',
              'fastafile=s',
              'dnd_file|d=s',
              'output|o=s',
              'analysis_conf|a=s',
              'log|l=s',
              'debug=s',
              'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Ergatis::Logger::get_logger();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDERR} );
}

&check_parameters(\%options);

$logger->info("Instantiating the BSML builder object");
my $builder = new BSML::BsmlBuilder;
$logger->logdie("builder was not defined") if (!defined($builder));

my $analysis_name = &get_analysis_name($options{analysis_conf});
my $analysis_id = $analysis_name . '_analysis';
my $MSF_alignments = process_MSF_file("$options{msffile}");

if ($MSF_alignments->{'mol_type'} eq 'polypeptide') {
    $MSF_alignments->{'mol_type'} = 'protein';
} else {
    $MSF_alignments->{'mol_type'} = 'nucleotide';
}

my $seq_class = {
                        'protein'       => 'polypeptide',
                        'nucleotide'    => 'nucleotide',
                   };
my $seq_molecule = {
                        'protein'       => 'aa',
                        'nucleotide'    => 'na',
                   };

if (keys %$MSF_alignments > 1) {   #skip empty msf files
    my $table = $builder->createAndAddMultipleAlignmentTable('molecule-type' => $MSF_alignments->{'mol_type'},
                                                             'class' => 'match' );
    
    if ( defined $options{dnd_file} ) {
        open(my $dnd_fh, "<$options{dnd_file}") || die "can't read DND file: $!\n";
        undef $/;
        my $newick_tree = <$dnd_fh>;
           $newick_tree =~ s/\s*//g;
        $/ = "\n";
        
        $table->addBsmlAttr( 'newick_tree', $newick_tree );
    }
    
    $table->addattr('class', 'match');
    $logger->logdie("table was not defined") if (!defined($table));
    
    my $summary = $builder->createAndAddAlignmentSummary( 
                              'multipleAlignmentTable' => $table,
                              'seq-type'               => $MSF_alignments->{'mol_type'},
                              'seq-format'             => 'msf'
                              );
    $logger->logdie("summary was not defined") if (!defined($summary));
    
    
    my $aln = $builder->createAndAddSequenceAlignment( 'multipleAlignmentTable' => $table );
    $logger->logdie("aln was not defined") if (!defined($aln));
    $table->addBsmlLink('analysis', '#'."$analysis_id", 'computed_by');
    my $seqnum=0;
    my $sequences_tag;
    
    foreach my $seq (keys %{ $MSF_alignments->{'polypeptide'} }) {
        $logger->logdie("seq was not defined") if (!defined($seq));

        $seqnum++;

        my $alignment = join ('', @{ $MSF_alignments->{'polypeptide'}->{$seq}->{'alignment'} });
        $logger->logdie("alignment was not defined") if (!defined($alignment));


        my $align_length = $MSF_alignments->{'polypeptide'}->{$seq}->{'length'} if ((exists $MSF_alignments->{'polypeptide'}->{$seq}->{'length'}) and (defined($MSF_alignments->{'polypeptide'}->{$seq}->{'length'})));

        $logger->logdie("align_length was not defined") if (!defined($align_length));

        ## Add sequence stub
        if(!($builder->returnBsmlSequenceByIDR($seq))){
            my $seq_stub = $builder->createAndAddSequence( 
                                $seq,           #id
                                $seq,           #title
                                '',             #length
                                $seq_molecule->{$MSF_alignments->{'mol_type'}}, #molecule
                                $seq_class->{$MSF_alignments->{'mol_type'}},                  #class
                                                         );
            if ($options{'fastafile'}) {
                $builder->createAndAddSeqDataImport(
                                    $seq_stub, 
                                    'fasta', 
                                    $options{'fastafile'}, 
                                    '', 
                                    $seq
                                                   );
            }
            $seq_stub->addBsmlLink('analysis', '#' . $analysis_id , 'input_of');
        }

        ####### THIS BLOCK OF CODE SEEMS TO BE OBSOLETE DUE TO IDGENERATOR
        ##
        ##
        
        #IMPORTANT!!!!
        #In order to ensure that each seq in a multiple sequence alignment is truly
        #unique, the seq-name and name will be in the form "polypeptide_accession:seqnum"
        #i.e. (ana1.10005.m00234_polypeptide:1). 

        # don't apply this so-called repair to non-underscore-delimited (i.e., new) ids
        if ($seq !~ /\.polypeptide\./) {
            if (length($seq) == 30 && $seq !~ /_polypeptide$/){
                $logger->warn("polypeptide identifier '$seq' was truncated by clustalw, repairing now");
                $seq =~ s/_[^_]+$/_polypeptide/;
            }
            
            if ($alignment !~ /_polypeptide\s/){
                $logger->warn("polypeptide identifier in the alignment of '$seq' was truncated by clustalw, repairing now");
                $alignment =~ s/_[^_\s]+\s/_polypeptide /g;
            }
        }

        ##
        ##
        ###### END BLOCK
        
        $builder->createAndAddAlignedSequence(
                              'alignmentSummary' => $summary,
                              'seqnum'           => $seqnum,
                              'length'           => $align_length,
                              'name'             => "$seq:$seqnum"
                              );

        $builder->createAndAddSequenceData(
                           'sequenceAlignment' => $aln,
                           'seq-name'          => "$seq:$seqnum",
                           'seq-data'          => $alignment
                           ); 
        $sequences_tag .= "$seqnum:";
    }
    $aln->addattr( 'sequences', $sequences_tag );
}

my $algorithm = 'unknown';
my $program = 'unknown';

if ($analysis_name) {
    $algorithm = $analysis_name;
    $program = $analysis_name;

} elsif ( $options{msffile} =~ /\.clw/ ) {
    $algorithm = 'clustalw';
    $program = 'clustalw';
}

## add the analysis element
$builder->createAndAddAnalysis(
    id => $analysis_id,
    sourcename => $options{'output'},
    algorithm => $algorithm,
    program => $program
);

$builder->write( $options{'output'} );

sub check_parameters{
    my ($options) = @_;
    
    ## make sure required arguments were passed
    my @required = qw( msffile output );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
    
    ## handle some defaults
    #$options{optional_argument2}   = 'foo'  unless ($options{optional_argument2});
}

#-------------------------------------------------------------
# process_MSF_file()
#
#-------------------------------------------------------------
sub process_MSF_file {

    $logger->debug("Entered process_MSF_file") if $logger->is_debug;

    my $file = shift;
    $logger->logdie("file was not defined") if (!defined($file));

    my $MSF_alignments ={};
    open(MSF, "$file") or die "Unable to open $file due to $!";
    my $line;
    my $msf_type;
    while(defined($line = <MSF>) and $line !~ /^\/\//) {
    if( $line =~ /MSF:\s*([\S]+)\s*Type:\s*([\S]+)\s*Check/) {
        my $msf_length = $1;
        return undef if($msf_length == 0);   #abort if align_len = 0

        if($2 eq 'P') {
        $msf_type = 'polypeptide';
        }elsif($2 eq 'N') {
        $msf_type = 'nucleotide';
        }else {
        $msf_type = 'polypeptide';
        }
        $MSF_alignments->{'mol_type'} = $msf_type;
    }
    #if($line =~ /Name:\s+([\S]+)\s+[o]{2}\s+Len:\s+([\S]+)\s+Check:\s+([\S]+)\s+Weight:\s+([\S]+)/) {

    if($line =~ /Name:\s*([\S]+)\s*[o]*\s*Len:\s*([\S]+)\s*Check:\s*([\S]+)\s*Weight:\s*([\S]+)/) {
        my $name    = $1;
        my $ali_len = $2;
        my $check   = $3;
        my $weight  = $4;
        
        $MSF_alignments->{'polypeptide'}->{$name}->{'length'} = $ali_len;
        $MSF_alignments->{'polypeptide'}->{$name}->{'check'}  = $check;
        $MSF_alignments->{'polypeptide'}->{$name}->{'weight'} = $weight;
        $MSF_alignments->{'polypeptide'}->{$name}->{'alignment'} = [];
    }
    }

    my $replacements;
    my $spaces;
    while($line = <MSF>) {
    if($line =~ /^([\S]+)/) {
        my $name = $1;
        if(exists($MSF_alignments->{'polypeptide'}->{$name})) {
        push( @{ $MSF_alignments->{'polypeptide'}->{$name}->{'alignment'} }, $line );
            } else {
        print STDERR "ERROR, $name is not valid polypeptide name for $file\n";
        exit;
            }
    }
    }

    return $MSF_alignments;

}#end sub process_MSF_file()


sub get_analysis_name{
    my($conf) = @_;
    my $analysis_name = "clustalw_analysis";
    
    if(-e $conf){
        my $cfg = new Config::IniFiles( -file => $conf);
        $analysis_name = $cfg->val( 'component', '$;COMPONENT_NAME$;' ) || $analysis_name;
    }
    return $analysis_name;
}
