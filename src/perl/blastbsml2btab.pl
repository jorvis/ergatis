#!/usr/bin/perl

=head1  NAME 

blastbsml2btab.pl - Converts blast bsml into btab

=head1 SYNOPSIS

USAGE: hmmpfam2bsml.pl 
        --input=/path/to/somefile.blast.bsml
        --output=/path/to/somefile.blast.btab
     [  --log=/path/to/some.log
        --debug=4 
        --help
      ]

=head1 OPTIONS

B<--input,-i> 
    BSML file from blast run.

B<--output,-o> 
    Output btab file

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--help,-h> 
    This help message

=head1   DESCRIPTION

    This script is used to convert blast bsml into the btab format (described below)

=head1 INPUT

    The input file to this script will be a blast bsml file.  Generally, in the format 
    shown here:

<?xml version="1.0"?>

<Bsml>
  <Definitions>
    <Sequences>
      <Sequence length="57" class="polypeptide" title="stpe.polypeptide.254786074.1" id="stpe.polypeptide.254786074.1" molecule="aa">
        <Attribute name="defline" content="stpe.polypeptide.254786074.1"></Attribute>
        <Seq-data-import source="/usr/local/projects/aengine/output_repository/translate_sequence/254740466_translate/i1/g1/stpe.polypeptide.254786074.1.fsa" identifier="stpe.polypeptide.254786074.1" format="fasta" id="Bsml0"></Seq-data-import>
        <Link rel="analysis" href="#wu-blastp_analysis" role="input_of"></Link>
      </Sequence>
      <Sequence length="57" class="polypeptide" title="hypothetical protein taxon:279808 {Staphylococcus haemolyticus JCSC1435;} (exp=0; wgp=1; cg=1; closed=1; pub=1; rf_status=provisional;)" id="RF_YP_252946.1_70726032_NC_007168" molecule="aa">
        <Attribute name="primary_label" content="RF|YP_252946.1|70726032|NC_007168"></Attribute>
        <Seq-data-import source="AllGroup.niaa" identifier="RF|YP_252946.1|70726032|NC_007168" format="fasta" id="Bsml1"></Seq-data-import>
      </Sequence>
    </Sequences>
    <Tables id="BsmlTables">
    <Seq-pair-alignment refend="57" compseq="RF_YP_252946.1_70726032_NC_007168" compxref="AllGroup.niaa:RF|YP_252946.1|70726032|NC_0071
68" refseq="stpe.polypeptide.254786074.1" refstart="0" reflength="57" class="match" refxref="/usr/local/projects/aengine/output_repositor
y/translate_sequence/254740466_translate/i1/g1/stpe.polypeptide.254786074.1.fsa:stpe.polypeptide.254786074.1" method="BLASTP">
        <Attribute name="percent_coverage_compseq" content="96.5"></Attribute>
        <Attribute name="percent_coverage_refseq" content="96.5"></Attribute>
        <Attribute name="percent_identity" content="66.1"></Attribute>
        <Attribute name="percent_similarity" content="82.1"></Attribute>
        <Seq-pair-run refcomplement="0" runprob="1.9e-13" comppos="0" refpos="0" runlength="55" compcomplement="0" comprunlength="55" run
score="70.2">
          <Attribute name="class" content="match_part"></Attribute>
          <Attribute name="p_value" content="1.9e-13"></Attribute>
          <Attribute name="percent_coverage_compseq" content="96.5"></Attribute>
          <Attribute name="percent_coverage_refseq" content="96.5"></Attribute>
          <Attribute name="percent_identity" content="66.1"></Attribute>
          <Attribute name="percent_similarity" content="82.1"></Attribute>
        </Seq-pair-run>
        <Link rel="analysis" href="#wu-blastp_analysis" role="computed_by"></Link>
      </Seq-pair-alignment>
      </Tables>
    </Definitions>
</Bsml>

=head1 OUTPUT

    The output btab file is a tab delimited file with the columns as described below:

    0: query_id
    1: ??
    2: query_length
    3: algorithm
    4: match (subject) db
    5: match id
    6: query start
    7: query stop
    8: match start
    9: match stop
    10: percent identity
    11: percent similarity
    12: score (raw)
    13: bit score
    14: ??
    15: match name
    16: blast frame (??)
    17: query strand
    18: subject protein length
    19: e-value
    20: p-value

=head1 CONTACT

    Kevin Galens
    kgalens@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Ergatis::Logger;
use XML::Twig;
use File::OpenFile qw(open_file);
use Data::Dumper;

########## GLOBALS ##########
my $input;             #input file
my $output;            #output file
my $data;              #holds parsed info
#############################

my %options = ();
my $results = GetOptions (\%options, 
              'input|i=s',
              'output|o=s',
              'log|l=s',
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

my $in = open_file( $input, "in" );
my $twig = new XML::Twig( 'twig_handlers' => {
    'Sequence' => \&sequence_handler,
    'Seq-pair-alignment' => \&spa_handler,
    'Analysis' => \&analysis_handler,
    } );

$twig->parse($in);
close($in);

&print_btab( $output );
print "$output\n";

sub sequence_handler {
    my ($twig, $el) = @_;
    
    my $db;
    my $title;
    my $real_id;

    #get the id
    my ($att) = $el->find_nodes('Attribute[@name="primary_label"]');
    if( $att ) {

        #this means it's a subject sequence and we should get the db name
        my $sdi = $el->first_child('Seq-data-import');
        $logger->logdie("Could not parse Seq-data-import from sequence ".$el->att('id')) unless( $sdi );
        $db = $sdi->att('source');
        $title = $el->att('title');
        $real_id = $att->att('content');

    }

    my $id = $el->att('id');
    my $len = $el->att('length');

    $data->{'sequences'}->{$id}->{'length'} = $len;
    $data->{'sequences'}->{$id}->{'database'} = $db if( $db );
    $data->{'sequences'}->{$id}->{'title'} = $title if( $title );
    $data->{'sequences'}->{$id}->{'real_id'} = $real_id if( $real_id );
}

sub spa_handler {
    my ($twig, $el) = @_;
    my $spa = {};

    my $query = $el->att('refseq');
    my $match = $el->att('compseq');

    my @sprs = $el->children( 'Seq-pair-run' );
    foreach my $spr ( @sprs ) {
        $spa->{'query'} = $query;
        $spa->{'match'} = $match;
        $spa->{'query_start'} = $spr->att('refpos') + 1;
        $spa->{'query_stop'}  = $spr->att('refpos') + $spr->att('runlength');
        $spa->{'match_start'} = $spr->att('comppos') + 1;
        $spa->{'match_stop'} =  $spr->att('comppos') + $spr->att('comprunlength');
        
        my ($per_ident_att) = $spr->find_nodes('Attribute[@name="percent_identity"]');
        $logger->logdie("Could not parse Attribute conting percent identity from alignment ".
                        "$query and $match") unless( $per_ident_att );
        $spa->{'percent_identity'} = $per_ident_att->att('content');

        my ($per_sim_att) = $spr->find_nodes('Attribute[@name="percent_similarity"]');
        $logger->logdie("Could not parse Attribute conting percent similarity from alignment ".
                        "$query and $match") unless( $per_sim_att );
        $spa->{'percent_similarity'} = $per_sim_att->att('content');


        my ($p_val_att) = $spr->find_nodes('Attribute[@name="p_value"]');
        $logger->logdie("Could not parse Attribute conting p-value from alignment ".
                        "$query and $match") unless( $p_val_att );
        $spa->{'p_value'} = $p_val_att->att('content');
        
        $spa->{'bit_score'} = $spr->att('runscore');
        $spa->{'e_value'} = $spr->att('runprob');

        push( @{$data->{'spas'}}, $spa );
    }

}

sub analysis_handler {
    my ($twig, $el) = @_;
    my ($algorithm_att) = $el->find_nodes('Attribute[@name="algorithm"]');
    $data->{'algorithm'} = $algorithm_att->att('content');
}

sub print_btab {
    my ($outfile) = @_;

    my $out = open_file( $outfile, 'out' );
    
    my @line;
    foreach my $spa ( @{$data->{'spas'}} ) {
        my $query = $spa->{'query'};
        my $match = $spa->{'match'};

        $line[0] = $query;
        $line[1] = "";
        $line[2] = $data->{'sequences'}->{$query}->{'length'};
        $line[3] = $data->{'algorithm'};
        $line[4] = $data->{'sequences'}->{$match}->{'database'};
        $line[5] = $data->{'sequences'}->{$match}->{'real_id'};
        $line[6] = $spa->{'query_start'};
        $line[7] = $spa->{'query_stop'};
        $line[8] = $spa->{'match_start'};
        $line[9] = $spa->{'match_stop'};
        $line[10] = $spa->{'percent_identity'};
        $line[11] = $spa->{'percent_similarity'};
        $line[12] = "";
        $line[13] = $spa->{'bit_score'};
        $line[14] = "";
        $line[15] = $data->{'sequences'}->{$match}->{'title'};
        $line[16] = 0;
        $line[17] = "";
        $line[18] = $data->{'sequences'}->{$match}->{'length'};
        $line[19] = $spa->{'e_value'};
        $line[20] = $spa->{'p_value'};

        print $out join("\t",@line);
        print $out "\n";
    }
}

sub check_parameters {
    my $opts = shift;

    if( $opts->{'input'} ) {
        $input = $opts->{'input'};
    } else {
        $logger->logdie("Option --input is required");
    }

    if( $opts->{'output'} ) {
        $output = $opts->{'output'};
    } else {
        $logger->logdie("Option --output is required");
    }
}
