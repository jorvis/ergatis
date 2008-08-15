#!/usr/bin/perl

=head1 NAME

assign_ec_from_priam_rps.pl - assigns ec numbers to gene predictions based on rpsblast against priam
    profile database.

=head1 SYNOPSIS

 USAGE: assign_ec_from_rps.pl
       --input_file=/path/to/some/rpsblast.bsml
       --output=/path/to/ec_numbers.bsml
       --evalue_cutoff=1e-10
     [ --log=/path/to/file.log
       --debug=4
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>
    BSML output generated from blastbtab2bsml.pl (run on rps blast output).  Will expect
    priam profiles as subject sequences

B<--output,-o>
    Path to BSML output file.  

B<--evalue_cutoff,-e>
    Will not consider any rps hits with an evalue higher than this [default: 1e-10]

B<--log,-l>
    Logfile.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

 DESCRIPTION

=head1  CONTACT

    Kevin Galens
    kgalens@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use MLDBM 'DB_File';
use Pod::Usage;
use Ergatis::Logger;
use XML::Twig;
use File::OpenFile qw(open_file);
use BSML::BsmlBuilder;
use Data::Dumper;

## globals #####
my $logger;
my $evalue_cutoff = 1e-10;
my $input_query_fsa;
my $input_file;
my $output;
my $ec_numbers = {};
my %seq_data_import_info;
###############

&check_options();

my $rps_hits = &parse_rps_bsml_file( $input_file );

 SEQUENCE:
    foreach my $query_seq ( keys %{$rps_hits} ) {
        $ec_numbers->{$query_seq} = {};
        
        my $hits = $rps_hits->{$query_seq};
        
        ## skip if we didnt' have any hits
        next SEQUENCE if( @{$hits} == 0 );

        my @sources;

        ## sort by e_value 
        my @sorted_hits = sort { $a->{'e_value'} <=> $b->{'e_value'} } @{$hits};

        my $top_hit = shift( @sorted_hits );
        my $ec_num = $top_hit->{'priam_id'};
        push( @sources, $ec_num );
        $ec_num =~ s/^(\d+)p//;

        $ec_numbers->{$query_seq}->{'ec_number'} = $ec_num;
        $ec_numbers->{$query_seq}->{'source'} = join(", ", @sources);
    }

## now that we know the ec number, print the bsml
&print_bsml( $output, $ec_numbers );


sub print_bsml {
    my ( $output, $ec_numbers ) = @_;
    my $doc = new BSML::BsmlBuilder;

    my $analysis_name = "assign_ec_from_rps_analysis";

    foreach my $query_seq ( keys %{$ec_numbers} ) {

        # assuming the input sequences to rpsblast were proteins
        my $seq = $doc->createAndAddSequence( $query_seq, $query_seq, '', 'aa', 'polypeptide' );
        my $link = $doc->createAndAddLink( $seq, 'analysis', '#'.$analysis_name, 'input_of' );

        my ($format, $source, $identifier) = &get_seq_data_import_info( $query_seq );
        
        my $sdi = $doc->createAndAddSeqDataImport( $seq, $format, $source, '', $identifier );

        ## if we've found an ec number for the query sequence, add an Attribute-list
        if( exists( $ec_numbers->{$query_seq}->{'ec_number'} ) ) {

            my $att_list = [ { 'name' => 'EC',
                               'content' => $ec_numbers->{$query_seq}->{'ec_number'},
                           },
                             { 'name' => 'IEA',
                               'content' => $ec_numbers->{$query_seq}->{'source'},
                           }];

            $seq->{'BsmlAttributeList'} = [$att_list];
            
        }        
    }

    my $analysis = $doc->createAndAddAnalysis( 'id' => $analysis_name,
                                               'program' => 'assign_ec_from_rps',
                                               'programversion' => 'current',
                                               'algorithm' => 'assign_ec_from_rps' );

    $doc->write($output);
    print "wrote: $output\n";
}

sub parse_rps_bsml_file {
    my $input_file = shift;

    my $hits = {};

    my $twig = new XML::Twig( 'twig_handlers' => {
        'Seq-pair-alignment' => sub { &handle_seq_pair_alignment(@_, $hits) },
        'Sequence' => sub { &handle_sequence(@_, $hits) },
    });

    my $in = open_file( $input_file, 'in' );
    $twig->parse($in);
    close($in);
    
    return $hits;             
}

sub handle_seq_pair_alignment {
    my ($twig, $spa_elem, $hits) = @_;

    ## grab the sequence ids
    my $priam_id = $spa_elem->att('compseq');
    my $query_id = $spa_elem->att('refseq');

    ## remove _ from begining of priam id (if there).
    ## this is added so that it is a valid bsml id (can't start with a number).
    $priam_id =~ s/^\_//;

    foreach my $spr_elem ( $spa_elem->children('Seq-pair-run') ) {
        
        my $spr;

        ## scoring info
        $spr->{'e_value'} = $spr_elem->att('runprob');
        $spr->{'score'}   = $spr_elem->att('runscore');
        
        $spr->{'percent_coverage_refseq'} = 
            ($spr_elem->find_nodes('Attribute[@name="percent_coverage_refseq"]'))[0]->att('content');
        $spr->{'percent_coverage_compseq'} = 
            ($spr_elem->find_nodes('Attribute[@name="percent_coverage_compseq"]'))[0]->att('content');
        $spr->{'percent_identity'} = 
            ($spr_elem->find_nodes('Attribute[@name="percent_identity"]'))[0]->att('content');

        $spr->{'priam_id'} = $priam_id;
        $spr->{'query_id'} = $query_id;

        push( @{$hits->{$spr->{'query_id'}}}, $spr);

    }        
}

sub handle_sequence {
    my ($twig, $seq_elem, $hits) = @_;
    
    ## grab and store the seq-data import information where available
    my $seq_id = $seq_elem->att('id');
    my $sdi_elem = $seq_elem->first_child('Seq-data-import');

    if( $sdi_elem ) {
        &store_seq_data_import_info( $seq_id,
                                     $sdi_elem->att('format'),
                                     $sdi_elem->att('source'),
                                     $sdi_elem->att('identifier') );

        $hits->{$seq_id} = [] unless( exists( $hits->{$seq_id} ) );
    }
    
}

sub store_seq_data_import_info {
    my ($id, $format, $source, $identifier) = @_;

    $seq_data_import_info{$id} = {
        'format' => $format,
        'source' => $source,
        'identifier' => $identifier,
    };
}


sub get_seq_data_import_info {
    my ($id) = @_;
    my @retval;

    if( exists( $seq_data_import_info{$id} ) ) {
        my $info = $seq_data_import_info{$id};
        @retval = ( $info->{'format'},
                    $info->{'source'},
                    $info->{'identifier'}, );
    }

    return @retval;

}

sub check_options {
    
    my %options;
    my $results = GetOptions (\%options,
                              'input_file|i=s',
                              'output|o=s',
                              'evalue_cutoff|e=s',
                              'help|h',
                              'debug|d=s',
                              'log|l=s');

    #setup the logger
    my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
    $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
    $logger = $logger->get_logger();

    if( $options{'help'} ) {
        &_pod;
    }

    my @required_opts = qw( input_file output );

    #required options check
    foreach my $req ( @required_opts ) {
        unless( $options{$req} ) {
            print STDERR "Option $req is required\n";
            &_pod;
        }
    }

    if( $options{'evalue_cutoff'} ) {
        $evalue_cutoff = $options{'evalue_cutoff'};
    }

    $input_file = $options{'input_file'};
    $output = $options{'output'};
}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => *STDERR} );
}
