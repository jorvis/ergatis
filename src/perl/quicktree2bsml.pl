#!/local/packages/perl-5.8.8/bin/perl

BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

qiucktree2bsml.pl - Convert output from quicktree to BSML.

=head1 SYNOPSIS

USAGE: quicktree2bsml.pl
            --quicktree_output=/path/to/output/of/quicktree
            --quicktree_input=/path/to/input/of/quicktree
            --output=/path/to/output.bsml
          [ --log=/path/to/some.log
            --debug=4
          ]

=head1 OPTIONS

B<--quicktree_input,-i>
    Input multiple sequence alignment file, used as quicktree input.  Alternatively,
    this may be the distance matrix used as input.  If the input is a MSA file, the
    alignment is stored in the bsml.  If only given a distance matrix, we can at the
    moment only store the sequence ids, with no seq-data-import and no alignment info.

B<--quicktree_output,-j>
    Quicktree output in Newick (New Hampshire) format.

B<--output,-o>
    The full path to the BSML file that will be created.

B<--debug> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

Converts a Stockholm formatted Multiple Sequence alignment (or distance matrix)
and the Newick formatted output of quicktree from that alignment into BSML.


=head1  INPUT

Input consists of two files:

The first file is the quicktree input, either a Stockholm format alignment file 
which has a format like this:

# STOCKHOLM 1.0
#=GF ID    GP120
#=GF AC    PF00516
#=GF DE    Envelope glycoprotein GP120
#=GF GA    -267.6 -267.6
#=GF NC    -267.7 -267.7
#=GF TC    -267.6 -267.6

ENV_SIVM1/24-528  QYVTVFYGVPAWRNATIPLFCATKNR.......DTWGTTQCLPDNDDYSE
ENV_HV2NZ/24-502  QFVTVFYGIPAWRNASIPLFCATKNR.......DTWGTIQCLPDNDDYQE
ENV_HV2G1/23-502  QYVTVFYGVPVWRNASIPLFCATKNR.......DTWGTIQCKPDNDDYQE
ENV_HV2D1/24-501  QYVTVFYGIPAWRNASIPLFCATKNR.......DTWGTIQCLPDNDDYQE

ENV_SIVM1/24-528  LALN.VTESFDAWE..NTVTEQAIEDVWQLFETSIKPCVKLSPLCITMRC
ENV_HV2NZ/24-502  ITLN.VTEAFDAWN..NTVTEQAVEDVWNLFETSIKPCVKLTPLCVAMNC
ENV_HV2G1/23-502  ITLN.VTEAFDAWD..NTVTEQAVEDVWSLFETSIKPCVKLTPLCVAMSC
ENV_HV2D1/24-501  ITLN.VTEAFDAWD..NTVTEQAIEDVWRLFETSIKPCVKLTPLCVAMNC

ENV_SIVM1/24-528  NKSETDKWGLTKSSTTTASTTTTTTAKSV.......ETRDIVNETSPCVV
ENV_HV2NZ/24-502  TR...........N.....MTTWTGRTDT.......QNITIINDTS.HAR
ENV_HV2G1/23-502  N........STTNN.....TTT...TGST.......TGMSEINETS.PSY
ENV_HV2D1/24-501  NITSGTTATPSP........................PNITIIDENSTCIG

ENV_SIVM1/24-528  TSRNKR
ENV_HV2NZ/24-502  Q.RHTR
ENV_HV2G1/23-502  V.RNKR
ENV_HV2D1/24-501  V.RNKR

... or a distance matrix, something like this:

ENV_SIVM1/   0.00000   0.27848   0.24948   0.28059
ENV_HV2NZ/   0.27848   0.00000   0.19533   0.19313
ENV_HV2G1/   0.24948   0.19533   0.00000   0.16094
ENV_HV2D1/   0.28059   0.19313   0.16094   0.00000

The second file contains the quicktree output, in New Hampshire (Newick) format.
Newick trees look something like:

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

=head1 CONTACT

    Jason Inman
    jinman@tigr.org

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Config::IniFiles;
use Ergatis::Logger;
use BSML::BsmlBuilder;


my $q_input;
my $q_output;
my $output;
my $class = '';
my $input_type;
my %options= ();

my $results = GetOptions (\%options, 
                'quicktree_input|i=s',
                'quicktree_output|j=s',
                'output|o=s',
                'class|c=s',
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

# Prepare the input data
my $analysis_id = 'quicktree_analysis';
my $seq_data = process_quicktree_input_file($q_input);

if ( $seq_data->{'mol_type'} eq 'polypeptide') {
    $seq_data->{'mol_type'} = 'protein';
} else {
    $seq_data->{'mol_type'} = 'nucleotide';
}

my $seq_class = {
                        'protein'       => 'polypeptide',
                        'nucleotide'    => 'nucleotide',
                   };
my $seq_molecule = {
                        'protein'       => 'aa',
                        'nucleotide'    => 'na',
                   };

if (keys %$seq_data > 1) {   #skip empty msf files

    my $table = $builder->createAndAddMultipleAlignmentTable('molecule-type' => $seq_data->{'mol_type'},
                                                             'class' => 'match' );

    # take care of the tree
    open(my $tree_fh, "< $q_output") || die "can't read quicktree output: $!\n";
    undef $/;
    my $newick_tree = <$tree_fh>;
       $newick_tree =~ s/\n//g;
    $/ = "\n";
    
    $table->addBsmlAttr( 'newick_tree', $newick_tree );

    
    $table->addattr('class', 'match');
    $logger->logdie("table was not defined") if (!defined($table));
    
    my $summary = $builder->createAndAddAlignmentSummary( 
                              'multipleAlignmentTable' => $table,
                              'seq-type'               => $seq_data->{'mol_type'},
                              'seq-format'             => 'msf'
                              );
    $logger->logdie("summary was not defined") if (!defined($summary));
    
    
    my $aln = $builder->createAndAddSequenceAlignment( 'multipleAlignmentTable' => $table );

    $logger->logdie("aln was not defined") if (!defined($aln));
    $table->addBsmlLink('analysis', '#'."$analysis_id", 'computed_by');
    my $seqnum=0;
    my $sequences_tag;

    
    foreach my $seq (keys %{ $seq_data->{'polypeptide'} }) {
   
        $logger->logdie("seq was not defined") if (!defined($seq));

        $seqnum++;

        my $alignment;
        my $aln_length = 0;

        if ($input_type eq 'stockholm') {
            $alignment = join ('', @{ $seq_data->{'polypeptide'}->{$seq}->{'alignment'} });
            $logger->logdie("alignment was not defined") if (!defined($alignment));
            $_ = $alignment;
            $aln_length = length($alignment) - tr/\n//;
        }

        ## Add sequence stub
        if (!($builder->returnBsmlSequenceByIDR($seq))) {
            my $seq_stub = $builder->createAndAddSequence( 
                                $seq,           #id
                                $seq,           #title
                                '',             #length
                                $seq_molecule->{$seq_data->{'mol_type'}}, #molecule
                                $seq_class->{$seq_data->{'mol_type'}},                  #class
                                                         );

            $seq_stub->addBsmlLink('analysis', '#' . $analysis_id , 'input_of');

        }

        $builder->createAndAddAlignedSequence(
                                'alignmentSummary'  => $summary,
                                'seqnum'            => $seqnum,
                                'name'              => "$seq:$seqnum",
                                'length'            => $aln_length,
                              );

        my @cAASD_list = ('sequenceAlignment',$aln,'seq-name',"$seq:$seqnum");
        push (@cAASD_list,('seq-data',$alignment)) if ($input_type eq 'stockholm');
        $builder->createAndAddSequenceData( @cAASD_list );

        $sequences_tag .= "$seqnum:";

    }
    $aln->addattr( 'sequences', $sequences_tag );
}

my $algorithm = 'quicktree';
my $program   = 'quicktree';

## add the analysis element
$builder->createAndAddAnalysis(
    id => $analysis_id,
    sourcename => $output,
    algorithm => $algorithm,
    program => $program
);

$builder->write( $output );

$logger->debug("FINISHED");

exit(0);


sub check_parameters{

    my $error = '';

    if (exists $options{'quicktree_input'}) {
        $q_input = $options{'quicktree_input'};
    } else {
        $error .= "quicktree_input must be given.\n";
    }

    if (exists $options{'quicktree_output'}) {
        $q_output = $options{'quicktree_output'};
    } else {
        $error .= "quicktree_output must be given.\n";
    }

    if (exists $options{'output'}) {
        $output = $options{'output'};
    } else {
        $error .= "output must be given.\n";
    }

    if (exists $options{'class'}) {
        $class = $options{'class'};
        unless ($class eq 'polypeptide' || $class eq 'assembly') {
            $error .= "unexpected vallue for option --class.\n";
        }
    } else {
        $error .= "class must be given.\n";
    }

   $logger->logdie($error) if $error;
    
}

#-------------------------------------------------------------
# process_quicktree_input_file()
#
#-------------------------------------------------------------
sub process_quicktree_input_file {

    $logger->debug("Entered process_quicktree_input_file") if $logger->is_debug;

    my $file = shift;
    $logger->logdie("file was not defined") if (!defined($file));

    my $seq_data; 
    open(MSF, "$file") or die "Unable to open $file due to $!";

    my $first_line = <MSF>;
  
    if ($first_line =~ /^\t\d+/) {
        $input_type = 'dsm';
    } elsif ($first_line =~ /# STOCKHOLM/) {
        $input_type = 'stockholm';
    } else {
        $logger->logdie("Unrecognized quicktree input! Expecting distance matrix or stockholm MSA");
    }

    $seq_data->{'mol_type'} = $class;

    if ($input_type eq 'dsm') {
        &process_dsm(*MSF,$seq_data);
    } else {
        &process_stockholm(*MSF,$seq_data);
    }

    close MSF;

    return $seq_data;

}#end sub process_quicktree_input_file()

sub process_dsm {
# collect information from the dsm.

    my $handle = shift;
    my $seq_data = shift;


    while (<$handle>) {
        if ($_ =~ /^(\S+)\s/) {
            undef $seq_data->{$class}->{$1};
        } else {
            $logger->logdie("unexpected line format in $q_input");
        }
    }

    return $seq_data;

}

sub process_stockholm {
# collect information from the msa in stockholm format

    my $handle = shift;
    my $seq_data = shift;

    while (<$handle>) {

        # skip auxilliary information, if it's here
        last if (/\/\//);
        next if (/#(.*)/);
        next if (/^$/);

        if($_ =~ /(\S+)\s+(\S+\n)/) {
 
            my $name      = $1;
            my $alignment = $2;
        
            push @{$seq_data->{$class}->{$name}->{'alignment'}}, $2;

        }
 
    }

    return $seq_data;

}
