#!/usr/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

#This is where the perl doc goes.

################ MODULES #####################
use strict;
use warnings;
use XML::Twig;
use BSML::BsmlParserTwig;
use BSML::BsmlBuilder;
use IntervalTree;
use Ergatis::Logger;
use Ergatis::IdGenerator;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Data::Dumper;
#############################################

####################### GLOBALS ##############################
my $bsml_file;                      #Input bsml file
my $data = {};                      #Holds data about new orfs
my %sequence;                       #Holds sequence information ( fasta, length, etc )
my $minimum_orf_length = 300;       #Configurable param ( default: 100 )
my $id_generator;                   #Creates unique ids
my $output;                         #Output file name
my $project = 'parse';              #Used in id generation
my $sourcename;

#This should be read in somehow so it can be changed.
#But for now it's hardcoded until I find the files.
my $translation_table = {
    'start' => { 
        'forward' => { 'TTG' => 1, 'CTG' => 1, 'ATT' => 1,
                       'ATC' => 1, 'ATA' => 1, 'ATG' => 1,
                       'GTG' => 1 },
        'reverse' => { 'CAA' => 1, 'CAG' => 1, 'AAT' => 1,
                       'GAT' => 1, 'TAT' => 1, 'CAT' => 1,
                       'CAC' => 1 } },
    'stop' => {
        'forward' => { 'TAG' => 1, 'TAA' => 1, 'TGA' => 1 },
        'reverse' => { 'CTA' => 1, 'TTA' => 1, 'TCA' => 1 } },
};
###############################################################

############################ OPTIONS ##########################
my %options = ();
my $results = GetOptions (\%options, 
                          'input_bsml|i=s',
                          'output|o=s',
                          'project|p=s',
                          'id_repository|r=s',
                          'minimum_orf_length|m=s',
                          'sourcename|s=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || &_pod;
###############################################################

############# OPTIONS CHECKING AND LOG SETUP ##################
#Setup the logger
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

# Check the options.
&check_parameters(\%options);
#################################################################

my $intergenic_regions = &get_intergenic_regions($bsml_file );

my @final_orfs;
my $longest_orf = 0;
my $count_orfs = 0;
foreach my $seq ( keys %{$intergenic_regions} ) {

    foreach my $region ( @{$intergenic_regions->{$seq}} ) {
        my $subseq = substr( $sequence{$seq}->{'sequence'}, $region->[0], $region->[1] - $region->[0] );

        my $orfs = &find_orfs( $subseq, $region->[0] );

        #Filter out small orfs
        foreach my $orf ( @{$orfs} ) {
            $longest_orf = $orf->[1] - $orf->[0] if( $orf->[1] - $orf->[0] > $longest_orf );
            push( @final_orfs, $orf ) if( $orf->[1] - $orf->[0] >= $minimum_orf_length );
        }
        
    }

    $data->{$seq} = \@final_orfs;
    $count_orfs += scalar( @final_orfs );

    print "$count_orfs\n";
}

#Grab the correct number of ids
$id_generator->set_pool_size( 'exon' => $count_orfs,
                              'transcript' => $count_orfs,
                              'gene' => $count_orfs,
                              'CDS' => $count_orfs,
                              'polypeptide' => $count_orfs );

&print_bsml( $data, $bsml_file, $output );

######################### SUB ROUTINES #########################
sub print_bsml {
    my ($data, $bsml, $out) = @_;

    my $parser = new BSML::BsmlParserTwig;
    my $doc = new BSML::BsmlBuilder;

    $parser->parse( \$doc, $bsml );

    foreach my $seq ( keys %{$data} ) {
        
        #Try and retrieve the sequence.  If we can't, something went wrong, so die.
        my $bsml_seq = $doc->returnBsmlSequenceByIDR( $seq );
        $logger->logdie("Could not find sequence $seq in bsml file $bsml") unless( $bsml_seq );
        my $feature_table = $bsml_seq->returnBsmlFeatureTableR(0);

        foreach my $coord ( @{$data->{$seq}} ) {

            my %id_holder;

            foreach my $type ( qw( exon transcript gene CDS polypeptide ) ) {
                my $tmp_id = $id_generator->next_id( 'project' => $project,
                                                     'type' => $type );
                my $feat = $doc->createAndAddFeature( $feature_table, $tmp_id, '', $type );
                my ($start,$stop,$comp) = ( $coord->[0] < $coord->[1] ) ? ($coord->[0],$coord->[1],0) :
                    ($coord->[1],$coord->[0],1);

                my $int_loc = $doc->createAndAddIntervalLoc( $feat, $start, $stop, $comp );
                my $link = $doc->createAndAddLink( $feat, 'analysis', '#intergenic_analysis_analysis', 'computed_by' );

                $id_holder{$type} = $tmp_id;
                
            }

            my $fg = $doc->createAndAddFeatureGroup( $bsml_seq, '', $id_holder{'gene'} );

            foreach my $type ( qw( exon transcript gene CDS polypeptide ) ) {
                my $fgm = $doc->createAndAddFeatureGroupMember( $fg, $id_holder{$type}, $type );
            }
        }        
    }
    
    my @analysis_ids;
    foreach my $ba ( @{$doc->{'BsmlAnalyses'}} ) {
        my $aId = $ba->{'attr'}->{'id'};
        print "Found $aId\n";
        push( @analysis_ids, $ba->{'attr'}->{'id'});
    }
    $doc->{'BsmlAnalyses'} = [];

    my $analysis = $doc->createAndAddAnalysis( 'id' => 'intergenic_analysis_analysis',
                                               'algorithm' => 'intergenic_analysis',
                                               'version' => 'current',
                                               'program' => 'intergenic_analysis');

    $doc->createAndAddBsmlAttribute( $analysis, 'sourcename', $sourcename ) if($sourcename);
    

    $doc->write($out);

    #All the old links were pointing to the old analysis element.
    #So change them.
    local( *IN, *OUT );
    open( IN, "< $out") or $logger->logdie("Could not open $out ($!)");
    open( OUT, "> $out.tmp") or $logger->logdie("Could not open $out.tmp for writing ($!)");
    
    local $" = "|";
    while( my $line = <IN> ) {
        $line =~ s/(@analysis_ids)/intergenic_analysis_analysis/g;
        print OUT $line;
    }

    close IN;
    close OUT;

}

sub check_parameters {
    my $opts = shift;

    &_pod if( $opts->{'help'} );

    #Input bsml is required.
    if( $opts->{'input_bsml'} ) {
        $logger->logdie("bsml file: $opts->{'input_bsml'} does not exist") 
            unless( -e $opts->{'input_bsml'} );
        $bsml_file = $opts->{'input_bsml'};
    } else {
        $logger->logdie("option --input_bsml is required");
    }

    #output is required
    if( $opts->{'output'} ) {
        $output = $opts->{'output'};
    } else {
        $logger->logdie("option --ouput is required");
    }

    #project is optional
    if( $opts->{'project'} ) {
        $project = $opts->{'project'};
    }
    $project = &parse_project( $bsml_file ) if($project eq 'parse');
    $logger->logdie("Could not parse a value for project") unless( $project );

    #id repository is required
    if( $opts->{'id_repository'} ) {
        $id_generator = new Ergatis::IdGenerator( 'id_repository' => $opts->{'id_repository'} );
    } else {
        $logger->logdie("option --id_repository is required");
    }

    #minimum orf length is optional
    if( $opts->{'minimum_orf_length'} ) {
        $minimum_orf_length = $opts->{'minimum_orf_length'};
    }

    #sourcename is optional
    if( $opts->{'sourcename'} ) {
        $sourcename = $opts->{'sourcename'};
    }

}

sub parse_project {
    my $bsml = shift;
    my $proj;


    my $twig = new XML::Twig( 'twig_handlers' => {
        'Sequence' => sub {
            my $id = $_[1]->att('id');
            $logger->logdie("Could not parse id from sequence in $bsml.")
                unless( $id );
            ($proj) = split( /\./, $id );
        }
    } );

    $twig->parsefile( $bsml );

    return $proj;
}


#Takes a bsml file name and will return hash of gene regions ( non-overlapping set of intervals )
# $return_hash = { 'Sequence_id' => [ [start,stop,strand], [start,stop,strand] ... ] };
sub get_intergenic_regions {
    my ($bsml_file) = @_;
    my $retval;

    my $fh = &open_file( $bsml_file, 'in' );

    my $i_tree = new IntervalTree();  #Used to check overlaps
    my $gene_set = {};                #Holds the gene coordinates
    my $sequence_length = "";

    my $twig = new XML::Twig( 'twig_handlers' => {
        'Sequence' => sub { &_handle_sequence( @_, \$sequence_length, $i_tree, $gene_set ) },
    } );

    $twig->parse( $fh );
    close( $fh );

    $i_tree->buildTree;

    #Merge predictions
    foreach my $seq ( keys %{$gene_set} ) {
        my $merged = &merge_predictions( $i_tree, $gene_set->{$seq} );
        my @sorted = sort { $a->[0] <=> $b->[1] } @{$merged};


        if( @sorted == 0 ) {
            
            die(" Can't find sequence length for $seq" ) unless( $sequence_length );

            push( @{$retval->{$seq}}, [0,$sequence_length] );

        } else {
            
            push( @{$retval->{$seq}}, [ 0, $sorted[0]->[0] ] ) unless( $sorted[0]->[0] <= 0 );

            for( my $i = 0; $i < (scalar(@sorted) - 1); $i++ ) {
                push( @{$retval->{$seq}}, [$sorted[$i]->[1], $sorted[$i+1]->[0]] );
            }

            push( @{$retval->{$seq}}, [ $sorted[-1]->[1], $sequence_length ] )
                  unless( $sorted[-1]->[1] >= $sequence_length);

        }        
        
    }

    return $retval;
    
}

#Handler for the twig defined in sub get_intergenic_regions
sub _handle_features {
    my ( $twig, $feat, $i_tree, $gene_set ) = @_;
    my $feat_id = $feat->att('id');
    
    my $int_loc = $feat->first_child('Interval-loc');
    die("Unable to parse Interval-loc element from $feat_id") 
        unless( $int_loc );

    my ($start, $stop, $comp) = ($int_loc->att('startpos'),
                                 $int_loc->att('endpos'),
                                 $int_loc->att('complement') );

    die("startpos greater than endpos: $feat_id") if( $start > $stop );

    $i_tree->addInterval( $feat_id, $start, $stop );
    $gene_set->{$feat_id} = [$start, $stop, $comp];
   
    
}

#Get the sequence length (handler method for twig defined in get_intergenic_regions );
sub _handle_sequence {
    my( $twig, $seq_elem, $seq_len, $i_tree, $gene_set ) = @_;

    my $seq_id = $seq_elem->att('id');
    my $class = $seq_elem->att('class');

    #We don't want CDS or polypeptide sequences
    return if( $class eq 'CDS' || $class eq 'polypeptide' );

    my $len = $seq_elem->att('length');

    #If the attribute isn't defined get the length from the seq-data-import fasta
    unless( $len ) {
        my $sdi = $seq_elem->first_child('Seq-data-import');
        die("Could not parse Seq-data-import element from Sequence $seq_id")
            unless( $sdi );
        my $source = $sdi->att('source');
        die("Could not parse source from Seq-data-import (child of Sequence $seq_id)")
            unless( $source );
        my $ident = $sdi->att('identifier');
        die("Could not parse identifier from Seq-data-import (child of Sequence $seq_id)")
            unless( $ident );

        $sequence{$seq_id}->{'fasta'} = $source;
        ($len, $sequence{$seq_id}->{'sequence'}) = &get_fsa( $source, $ident );
    }

    die("Can't parse the length from $seq_id") unless( $len );
    $$seq_len =  $len;

    #Get all the gene features
    my @gene_features = $seq_elem->find_nodes('//Feature[@class="gene"]');

    my $features = {};

    foreach my $gene ( @gene_features ) {
        &_handle_features( $twig, $gene, $i_tree, $features );
    }

    $gene_set->{$seq_id} = $features;
}

#get_length_from_fsa
sub get_fsa {
    my ($source, $ident) = @_;
    
    my $fh = &open_file($source, 'in' );
    
    my $flag = 0;
    my $seq = "";
    while( my $line = <$fh> ) {
        chomp($line);

        if( $line =~ /^>(.*)/ ) {
            my $defline = $1;
            if( !$flag && $defline =~ /($ident)/ ) {
                $flag = 1;
            } elsif( $flag ) {
                last;
            }
            
        } elsif( $flag ) {
            $seq .= $line;
        }

    }

    my $retval = length($seq);
    
    close( $fh );
    return ( $retval, $seq );
}


#Flatten the gene predictions (combine those genes that overlap)
sub merge_predictions {
    my ($tree, $genes) = @_;

    #contains all the ids of matches already contained within an overlap.
    #so don't search overlaps for these sequences again.
    my %found;
    my @predictions;

    while ( 1 ) {

        #Retrieve a match to search with.
        #Just cycle through the array and pick first match that hasn't been found yet
        my $search;
        my ($gene_id, $gene);
        while( ($gene_id, $gene) = each( %{$genes} ) ) {
            next if( exists( $found{$gene_id} ) );
            $search = $gene;
        }
        last unless( $search );

        my $difference;
        my $prev_count = 1;
        do {

            my @overlaps = $tree->searchInterval( $search->[0], $search->[1] );

            
            my ($left, $right);
            foreach my $overlap ( @overlaps ) {
                $left = $overlap->[0] if( !$left || $overlap->[0] < $left );
                $right = $overlap->[1] if( !$right || $overlap->[1] > $right );
                $found{$overlap->[2]} = 1;
            }
            
            my $cur_count = scalar( @overlaps );
            $difference = $cur_count - $prev_count;
            
            $prev_count = $cur_count;
            
            $search = [$left, $right];
            
        } while( $difference != 0 );
        
        push( @predictions, $search );

    }

    return \@predictions;
    
}

#Pass in a sequence and you get back a multi-level hash of 
#start and stop codon positions, using $translation_table.
#Optional argument offset.  Will add this offset to the final
#coordinates returned.
sub find_orfs {
    my ($sequence, $offset) = @_;
    $offset = 0 unless( $offset );
    
    my $orfs;

    my @starts;

    my @orfs;

    for( my $i = 0; $i < length( $sequence ); $i++ ) {

        my $codon = substr( $sequence, $i, 3 );

        if( defined( $starts[ $i % 3 ] ) ) {
            if( exists( $translation_table->{'stop'}->{'forward'}->{$codon} ) ) {
                push( @orfs, [ $starts[ $i % 3 ] + $offset, $i + 3 +$offset ] );
                $starts[ $i % 3 ] = undef;
            }
        } else {
            if( exists( $translation_table->{'start'}->{'forward'}->{$codon} ) ) {
                $starts[ $i % 3 ] = $i;
            }
        }
    }

    @starts = undef;

    for( my $i = length( $sequence ); $i >= 0; $i-- ) {

        my $codon = substr( $sequence, $i, 3 );
        
        if( defined( $starts[ $i % 3 ] ) ) {
            if( exists( $translation_table->{'stop'}->{'reverse'}->{$codon} ) ) {
                push( @orfs, [ $starts[ $i % 3 ] + 3 + $offset, $i + $offset ] );
                $starts[ $i % 3 ] = undef;
            }
        } else {
            if( exists( $translation_table->{'start'}->{'reverse'}->{$codon} ) ) {
                $starts[ $i % 3 ] = $i;
            }
        }
        

    }
    return \@orfs;

}

#Wrapper for open in case the file is zipped.
sub open_file {
    my $file = shift;
    my $direction = shift;
    my $fh;
    
    if( $direction eq 'out' ) {
        if( $file =~ /\.gz$/ ) {
            open( $fh, ">:gzip", $file ) or die("can't open $file ($!)");
        } else {
            open( $fh, "> $file" ) or die("Can't open $file ($!)");
        }
    } elsif( $direction eq 'in' ) {
        
        if( -e $file ) {
            
            if( $file =~ /\.gz$/ ) {
                open( $fh, "<:gzip", $file ) or die("can't open $file ($!)");
            } else {
                open( $fh, "< $file") or die("can't open $file ($!)");
            } 
        } elsif( -e $file.".gz" ) {
            my $tmp = $file.".gz";
            open( $fh, "<:gzip", $tmp ) or die("Can't open $tmp ($!)");
        } else {
            die("Could not find $file or a gz version");
        }

    } else {
        die("Please specifiy a direction.  'in' or 'out'");
    }

    return $fh;
}

sub _pod {   
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
