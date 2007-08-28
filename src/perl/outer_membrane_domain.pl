#!/usr/local/bin/perl

use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

outer_membrane_domain.pl - looks for carboxy-terminal motif correlated with sorting to outer membrane

=head1 SYNOPSIS

USAGE: 
    outer_membrane_domain.pl 
        --input_fsa=/path/to/aa.fsa
        --output_bsml=/path/to/output.bsml
        --id_repository=/path/to/valid_id_repository
        --sourcename=/path/to
      [ --project=AFU
        --log
        --debug
        --help ]
      

=head1 OPTIONS

B<--input_fsa,-i>
    REQUIRED. Can be a single or multi fasta file.

B<--output_bsml,-o>
    REQUIRED. The output bsml file.

B<--id_repository,-r>
    REQUIRED. Id repository for use by Ergatis::IdGenerator.pm
    See Ergatis::IdGenerator.pm for more details.

B<--project,-p>
    OPTIONAL. The project (used for id generation). If left blank
    script will attempt to parse from existing ids.  If can't, program
    will error.

B<--log,-l>
    Logfile.

B<--debug,-d>
    Higher number = more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

    his program looks for a motif described in (insert reference here) that correlates with 
    proteins found in the outer membrane. The motif: the carboxy-terminal residue must be 
    aromatic (eg phenylalanine, tryptophan, or tyrosine) and at least two more of the 10 
    carboxy-terminal residues must also be aromatic.

=head1  INPUT

    Single or multi fasta file.

=head1 OUTPUT

    Bsml.  I'll [hopefully] fill this in later when that gets figured out.

=head1  CONTACT

    Kevin Galens
    kgalens@jcvi.org

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use BSML::BsmlBuilder;
use Ergatis::IdGenerator;
use Ergatis::Logger;

####### GLOBALS AND CONSTANTS ###########
my $input;                              #Holds the input file
my $project;                            #The project (ex aa1)
my $output;                             #Output file
my $idMaker;                            #The Ergatis::IdGenerator
my $bsml;                               #BSML::BsmlBuilder object object.
my $data;                               #Holds predicted data
my $debug;                              #The debug variable
my $sourcename;                         #For the analysis section of the bsml
my $type = 'outer_membrane_domain';     #The type of feature this will predict
#########################################

my %options = ();
my $results = GetOptions (\%options, 
                          'input_fsa|i=s',
                          'output_bsml|o=s',
                          'project|p=s',
                          'id_repository|r=s',
                          'sourcename|s=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || &_pod;

#Setup the logger
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

# Check the options.
&check_parameters(\%options);

#Gather all the sequences, find the length and store the last ten aa's
open( FSA, "< $input ") or $logger->logdie("Unable to open $input");

my ($seq, $id);
while( my $line = <FSA> ) {
    if( $line =~ /^>(\S+)/ ) {
        $data->{$id} = { 'length' => length( $seq ),
                         'end'    => substr( $seq, length($seq)-10 ) 
                         } if($seq);
        $id = $1;
        $seq = "";
        
    } else {
        chomp $line;
        $seq .= $line;
    }
}
$data->{$id} = { 'length' => length( $seq ),
                 'end'    => substr( $seq, length($seq)-10 ) 
                 } if( $seq );

#Generate the bsml
$bsml = &generateBsml( $data );

$bsml->write($output);
print "Wrote $output\n" if( $logger->is_debug );


######################## SUB ROUTINES #######################################
sub generateBsml {
    my $data = shift;
    my $doc = new BSML::BsmlBuilder();

    #Each of the sequences is going to be a poly peptide
    foreach my $pep ( keys %{$data} ) {
        
        my $seq  = $doc->createAndAddSequence( $pep, $pep, $data->{$pep}->{'length'}, 'aa');
        my $sdi  = $doc->createAndAddSeqDataImport( $seq, 'fasta', $input, '', $pep );
        my $link = $doc->createAndAddLink( $seq, 'analysis', '#outer_membrane_domain_analysis',
                                           'input_of' );

        #Find the project if it was not specified
        if( !$project || $project eq 'parse' ) {
            $project = $1 if( $pep =~ /^([^\.])+/ );
        }

        #############################################################
        #                      SEARCH FOR MOTIF                     #
        #############################################################
        
        if( $data->{$pep}->{'end'} =~ /([FWY].*){3,}/ && 
            $data->{$pep}->{'end'} =~ /[FWY]$/ ) {
            my $feature_table = $doc->createAndAddFeatureTable( $seq );
            $logger->logdie("project has not been determined.  Something is wrong.")
                unless( $project );
            my $feat_id = $idMaker->next_id( 'project' => $project,
                                             'type' => $type );

            my $feature = $doc->createAndAddFeature($feature_table, $feat_id, $feat_id, $type);
            my $int_loc = $doc->createAndAddIntervalLoc( $feature, $data->{$pep}->{'length'} - 10,
                                                         $data->{$pep}->{'length'}, 0 );
            my $link = $doc->createAndAddLink( $feature, 'analysis', '#outer_membrane_domain_analysis',
                                               'computed_by' );
        }
        
    }

    #Add the analysis
    $doc->createAndAddAnalysis( 'id' => 'outer_membrane_domain_analysis', 
                                'sourcename' => $sourcename,
                                'algorithm'  => 'outer_membrane_domain',
                                'programversion' => 'current'  );
    
    return $doc;
    
}

sub check_parameters {
    my $opts = shift;

    &_pod if($opts->{'help'});

    if($opts->{'input_fsa'}) {
        $logger->logdie("Option input_file ($options{'input_file'}) does not exist") 
            unless( -e $opts->{'input_fsa'});
        $input = $opts->{'input_fsa'};
    } else {
        $logger->logdie("Option input_file is required");
    }

    unless($opts->{'output_bsml'}) {
        $logger->logdie("Option output is required");
    } else {
        $output = $opts->{'output_bsml'};
    }

    unless($opts->{'project'}) {
        $logger->logdie("Option project is required");
    } else {
        $project = $opts->{'project'};
    }

    unless($opts->{'id_repository'}) {
        $logger->logdie("Option id_repository is required.  Please see Ergatis::IdGenerator ".
            "for details.");
    } else {
        $idMaker = new Ergatis::IdGenerator( 'id_repository' => $opts->{'id_repository'} );
        $idMaker->set_pool_size( 'outer_membrane_domain' => 20);                                 
    }

    unless( $opts->{'sourcename'} ) {
        $logger->logdie("Option sourcename is required");
    } else {
        $sourcename = $opts->{'sourcename'};
    }
    
    if($opts->{'debug'}) {
        $debug = $opts->{'debug'};
    }
    
}

sub _pod {   
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
