#!/usr/local/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

combine_ncRNA_alignments.pl - creates a bsml document from mulitple wu-blastn results ( against various ncRNA databases ).

=head1 SYNOPSIS

USAGE:  combine_ncRNA_alignments.pl --bsml_input input.bsml --other_bsml_lists lists/one.list,lists/two.list 
    --output_file output.bsml --id_repository valid_id_repository/ --project ORG

=head1 OPTIONS

B<--bsml_input,-b>
    Input wu-blastn bsml file (not a list). 

B<--other_bsml_lists,-l>
    Comma separated list of wu-blastn lists.

B<--output_file,-o>
    The name of the output file

B<--id_repository,-r>
    Valid id_repository.  See Ergatis::Id_Repository.

B<--project,-p>
    [Optional].  Used in id generation.  If left blank (or keyword parse is used) will parse
    the name from the sequence id.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1   DESCRIPTION

    

=head1  CONTACT

    Kevin Galens
    kgalens@jcvi.org

=cut

use lib (@INC, "/usr/local/devel/ANNOTATION/ard/current/lib/5.8.8/");

use strict;
use warnings;
use IntervalTree;
use Getopt::Long;
use XML::Twig;
use BSML::BsmlBuilder;
use Ergatis::Logger;
use Ergatis::IdGenerator;
use Data::Dumper;

########################################### SETUP #
my %options;
GetOptions( \%options,
            'bsml_input|b=s',
            'other_bsml_lists|l=s',
            'output_file|o=s',
            'hit_count_cutoff|c=s',
            'id_repository|r=s',
            'project|p=s',
            'log|l=s',
            'debug|d=s',
            'help|h=s');

########################################### GLOBALS / CONSTANTS #
my $logger;
my @input_files;
my @other_files;
my $output_file;
my $rna_type = "";
my $final_predictions = {};
my $hit_count_cutoff = 10;
my $sequences = {};
my %DB_HANDLER = (
    'RNA_16S.db'   => \&rRNA_16s,
    'RNA_23S.db'   => \&rRNA_23s,
    'RNA_5S.db'    => \&rRNA_5s,
    'RNA_SRP.db'   => \&SRP_RNA,
    'RNA_rnpB.db'  => \&rnpB_RNA,
    'RNA_tmRNA.db' => \&tmRNA 
    );
my $id_generator;
my $project;
my $direction = {};        #Holds the directions of the ncRNA genes.
################################################################
    

##Logging Setup
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
$logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

&check_parameters( \%options );

########################################## MAIN #
foreach my $in_file( @input_files ) {
    
    my $rna = {};
    my $twig = new XML::Twig( 'twig_handlers' => {
        'Sequence' => sub { 
            return unless( $_[1]->att('id') =~ /assembly/ );
            $sequences->{$_[1]->att('id')}->{'source'} = $_[1]->first_child('Seq-data-import')->att('source');
            $sequences->{$_[1]->att('id')}->{'identifier'} = $_[1]->first_child('Seq-data-import')->att('identifier');
            $sequences->{$_[1]->att('id')}->{'length'} = $_[1]->att('length') if( $_[1]->att('length') );
        },
                                                   
        'Seq-pair-alignment' => sub { &handle_spa( @_, $rna ) }
    });

    $twig->parsefile( $in_file );

    #For each sequence
    foreach my $seq_id ( keys %{$rna} ) {
        $rna->{$seq_id}->{'tree'}->buildTree();

        my $ids = &get_best_predictions( $rna->{$seq_id}->{'tree'},
                                         $rna->{$seq_id}->{'matches'});

        for ( @{$ids} ) {
            push( @{$rna->{$seq_id}->{'matches'}->{$_}}, $_ );
            push( @{$final_predictions->{$seq_id}->{$rna_type}}, 
                  $rna->{$seq_id}->{'matches'}->{$_});
        }
    }
}

&print_bsml( $final_predictions, $output_file );

########################################## SUB ROUTINES #
sub check_parameters {
    my $opts = shift;

    my $base;
    unless( $opts->{'bsml_input'} ) {
        $logger->logdie("option bsml_input is required");
    } else {
        push( @input_files, $opts->{'bsml_input'});
        $base = $1 if( $opts->{'bsml_input'} =~ m|.*/(([^/\.]+\.){3})[^/]+$|);
        $base =~ s/\.$//;
        $project = $1 if( $base =~ /^([^\.]+)/ );
    }

    if( $opts->{'other_bsml_lists'} ) {
        my @lists = split( /[,\s]+/,  $opts->{'other_bsml_lists'});
        foreach my $list ( split(/[,\s]+/, $opts->{'other_bsml_lists'}) ) {
            open(IN, "< $list") or
                $logger->logdie("Can't open $list ($!)");
            chomp( @other_files = <IN> );
            close(IN);
            
            for ( @other_files ) {
                push( @input_files, $_ ) if( m|/$base\.| );
            }
        }
    }

    unless( $opts->{'output_file'} ) {
        $logger->logdie("option output_file is required");
    } else {
        $output_file = $opts->{'output_file'};
    }

    unless( $opts->{'id_repository'} ) {
        $logger->logdie("Option id_repository is required.  See Ergatis::IdGenerator.pm ".
                        "for more details");
    }
    $id_generator = new Ergatis::IdGenerator( 'id_repository' => $opts->{'id_repository'} );

    if( $opts->{'project'} && $opts->{'project'} ne 'parse' ) {
        $project = $opts->{'project'};
    }

    $hit_count_cutoff = $opts->{'hit_count_cutoff'} if( $opts->{'hit_count_cutoff'} );
    

}
sub get_best_predictions {
    my ($tree, $matches) = @_;

    #contains all the ids of matches already contained within an overlap.
    #so don't search overlaps for these sequences again.
    my %found;
    my @predictions;

    while ( 1 ) {

        #Retrieve a match to search with.
        #Just cycle through array and pick first match that's hasn't been found yet
        my $search;
        while( my ($compseq, $match) = each( %{$matches} ) ) {
            next if( exists( $found{$compseq} ) );
            $search = $match;
        }
        last unless( $search );

        my $difference = 0;
        my @overlaps;
        my @first_overlaps = $tree->searchInterval( $search->[0], $search->[1] );

        my ($tmp_left, $tmp_right);
        my ($left_acc, $right_acc);
        do {
            my $count = scalar( @first_overlaps );

            my ($left, $right);
            foreach my $overlap ( @first_overlaps ) {
                $left_acc = $overlap->[2] if( !$left || $overlap->[0] < $left );
                $left = $overlap->[0] if( !$left || $overlap->[0] < $left );
                $right_acc = $overlap->[2]  if( !$right || $overlap->[1] > $right );
                $right = $overlap->[1] if( !$right || $overlap->[1] > $right );
            }
            ($tmp_left, $tmp_right) = ($left,$right);

            @overlaps = $tree->searchInterval( $left, $right );
            
            #Filter the overlaps based on what we've seen before:
            my @tmp;
            foreach my $ol ( @overlaps ) {
                next if( exists( $found{$ol->[2]} ) );
                push( @tmp, $ol );
            }
            @overlaps = @tmp;

            my $new_count = scalar( @overlaps );
            $difference = $new_count - $count;

            @first_overlaps = @overlaps if( $difference );

        } while( $difference != 0 );


        my %tmp_hash = map( ($_->[2],$matches->{$_->[2]}->[2]), @overlaps );
        %found = ( %tmp_hash, %found );

        #find the best score
        my @tmp;
        while( my ($k,$v) = each(%tmp_hash) ) {push(@tmp,[$k,$v])}

        my @sorted = sort { $b->[1] <=> $a->[1] } @tmp;
        my $best_match = shift(@sorted);

        push( @predictions, $best_match->[0]  ) 
            if( scalar( @overlaps ) > $hit_count_cutoff );

    }

    return \@predictions;
    
}

sub handle_spa {
    my ($twig, $spa_elem, $rna) = @_;

    my $seq_id = $spa_elem->att('refseq');
    
    #If there isn't a tree yet, make one
    unless( exists( $rna->{$seq_id}->{'tree'} ) ) {
        $rna->{$seq_id}->{'tree'} = new IntervalTree();
    }

    $rna_type = $1 if( $spa_elem->att('compxref') =~ /^([^:]+):/);

    my @sprs = $spa_elem->children('Seq-pair-run');
    my $compseq = $spa_elem->att('compseq');
    my $id_num = 1;
    foreach my $spr ( @sprs ) {
        my $start = $spr->att('refpos');
        my $stop = $start + $spr->att('runlength');
        my $score = $spr->att('runscore');
        $direction->{$compseq."_$id_num"} = $spr->att('refcomplement');
        $rna->{$seq_id}->{'tree'}->addInterval( $compseq."_$id_num", $start, $stop );
        $rna->{$seq_id}->{'matches'}->{$compseq."_$id_num"} = [$start, $stop, $score];
        $id_num++;
    }
    
    
}

sub print_bsml {
    my ($predictions, $out_file) = @_;

    my $doc = new BSML::BsmlBuilder;

    foreach my $seq ( keys %{$predictions} ) {
        my $seq_elem = $doc->createAndAddSequence( $seq, $seq, $sequences->{$seq}->{'length'}, 'na', 'assembly' );
        my $seq_data_import = $doc->createAndAddSeqDataImport( $seq_elem, 'fasta', $sequences->{$seq}->{'source'}, '', $sequences->{$seq}->{'identifier'} );
        my $link = $doc->createAndAddLink( $seq_elem, 'analysis', '#combine_ncRNA_alignments_analysis', 'computed_by' );
        my $feature_table = $doc->createAndAddFeatureTable( $seq_elem );

        foreach my $db ( keys %{$predictions->{$seq}} ) {
            foreach my $feature ( @{$predictions->{$seq}->{$db}} ) {
                my $hit_id = $feature->[3];
                $hit_id =~ s/(.*)_\d+/$1/;
                my $tmp_feature = [$feature->[0], $feature->[1], $direction->{$feature->[3]}, $feature->[3]];
                my $transcript_type = $db;
                $transcript_type =~ s/\.db$//;
                &print_feature_to_bsml( $tmp_feature, $seq_elem, 'ncRNA', $transcript_type, $doc, $feature_table );
            }
        }
    }

    my $sourcename = $1 if( $output_file =~ m|(.*/\d{4}_[^/]+)| );
    $sourcename = $1 if( $output_file =~ m|(.*)/[^/]+$| );
    $logger->logdie("Can't parse sourcename for analysis element [$output_file]") unless( $sourcename );
    my $analysis = $doc->createAndAddAnalysis('id'         => "combine_ncRNA_alignments",
                                              'sourcename' => $sourcename,
                                              'algorithm'  => 'combine_ncRNA_alignments',
                                              'program'    => 'combine_ncRNA_alignments',
                                              'version'    => '1.0' );
    $doc->write( $out_file );
    print "Wrote to $out_file\n";

}

sub print_feature_to_bsml {
    my ($feature, $seq, $class, $gpn, $doc, $feature_table) = @_;
    
    #Generate the ids
    my %ids;
    ($ids{'gene'}, $ids{'exon'}, $ids{$class}) = 
        ( $id_generator->next_id( 'project' => $project, 'type' => 'gene' ),
          $id_generator->next_id( 'project' => $project, 'type' => 'exon' ),
          $id_generator->next_id( 'project' => $project, 'type' => $class ));
                                    
    
    my $fg = $doc->createAndAddFeatureGroup( $seq, '', $ids{'gene'} );

    foreach my $type ( ( $class, 'exon', 'gene' ) ) {
        my $bsml_feature = $doc->createAndAddFeatureWithLoc( $feature_table, $ids{$type}, $ids{$type}, $type, '', '',
                                          $feature->[0], $feature->[1], $feature->[2] );
        $doc->createAndAddLink( $bsml_feature, 'analysis', '#combine_ncRNA_alignments_analysis', 'computed_by' );

        if( $type eq $class ) {
            $doc->createAndAddBsmlAttribute( $bsml_feature, 'gene_product_name', $gpn );
            $doc->createAndAddBsmlAttribute( $bsml_feature, 'gene_product_name_source', $feature->[3] );
        }

        $doc->createAndAddFeatureGroupMember( $fg, $ids{$type}, $type );
    }
}

################################################### DB_HANDLER FUNCTIONS ##
sub default_handler {
    print "In the default handler\n";
}
sub rRNA_16s {
    my ($doc, $rnas, $seq, $ft) = @_;

    print Dumper( $rnas );
    print $direction->{$rnas->[3]}."\n";
}
sub rRNA_23s {
    my ($doc, $rnas, $seq, $ft) = @_;
    
    print Dumper( $rnas );
    print $direction->{$rnas->[3]}."\n";
}
sub rRNA_5s {
    my ($doc, $rnas, $seq, $ft) = @_;

    print Dumper( $rnas );
    print $direction->{$rnas->[3]}."\n";
}       
sub tmRNA {
    my ($doc, $rnas, $seq, $ft) = @_;

    print Dumper( $rnas );
    print $direction->{$rnas->[3]}."\n";
}
sub SRP_RNA {
    my ($doc, $rnas, $seq, $ft) = @_;

    print Dumper( $rnas );
    print $direction->{$rnas->[3]}."\n";
}
sub rnpB_RNA {
    my ($doc, $rnas, $seq, $ft) = @_;

    
}
                

        
          
