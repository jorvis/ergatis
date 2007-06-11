#!/usr/local/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

use strict;
use warnings;
use XML::Twig;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Ergatis::Logger;
use Data::Dumper;

my $trans_poly_lookup = {};
my @input_files;
my $output_file;
my %options;

GetOptions(\%options, 
           'input_bsml|i=s',
           'output_file|o=s',
           'log|l=s',
           'debug=s',
           'help|h');

#Setup the logger
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                 'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();
&check_parameters( \%options );

#################################################################################

my $info = {};
foreach my $infile ( @input_files ) {
    
    my $twig = new XML::Twig( 'twig_handlers' => {
        'Sequence[@class="assembly"]' => sub {
            my $seq_elem = $_[1];
            my $id = $seq_elem->att('id');
            my $length;
            $length = $seq_elem->att('length') if( $seq_elem->att('length') );
            $length = &get_length_from_seq_elem( $seq_elem ) unless( $length );
            $logger->logdie("Unable to get length of molecule $id") unless( $length );
            $info->{$id}->{'length'} = $length;

            $info->{$id}->{'polypeptides'} = &get_polypeptides_from_seq_elem( $seq_elem );

            print "Found ".scalar( @{$info->{$id}->{'polypeptides'}} )." polypeptides for sequence $id in $infile.\n";
        } } );

    $twig->parsefile( $infile );
}

&generate_output_structure( $info );

#####################################################################################
sub generate_output_structure {
    my $info = shift;

    open( OUT, "> $output_file") or 
        $logger->logdie("Unable to open $output_file for writing ($!)");

    foreach my $seq_id ( keys %{$info} ) {
        
        foreach my $polypeptide ( @{$info->{$seq_id}->{'polypeptides'}} ) {
            print OUT "$seq_id\t".$info->{$seq_id}->{'length'}."\t".$trans_poly_lookup->{$polypeptide->{'feat_name'}};
            foreach my $term ( qw( end5 end3 com_name ) ) {
                print OUT "\t".$polypeptide->{$term};
            }
            print OUT "\n";
        }
    }
    
    close(OUT);
    print "Wrote to $output_file\n";

}
sub get_polypeptides_from_seq_elem {
    my $seq_elem = shift;
    my $polypeptides = [];

    my $seq_id = $seq_elem->att('id');

    my @feature_groups = $seq_elem->find_nodes('//Feature-group');
    print "Found ".scalar( @feature_groups )." feature-groups\n";
    
    foreach my $fg ( @feature_groups ) {
        my $group_set = $fg->att('group-set');
        my @feature_group_members = $fg->children( 'Feature-group-member' );
        my $tmp = {};
        foreach my $fgm ( @feature_group_members ) {
            next unless( $fgm->att('feature-type') =~ /(transcript|polypeptide)/);
            $tmp->{ $fgm->att('feature-type') } = $fgm->att('featref');
            
        }
        
        if( exists( $tmp->{'transcript'} ) && exists( $tmp->{'polypeptide'} ) ) {
            $trans_poly_lookup->{ $tmp->{'transcript'} } = $tmp->{'polypeptide'};
        }
    }

    my @features = $seq_elem->find_nodes('//Feature[@class="transcript"]');

    foreach my $feature ( @features ) {
        my $id = $feature->att('id');
        $logger->logdie("Couldn't find transcript id $id in lookup hash") unless( exists( $trans_poly_lookup->{$id} ));
        my $poly_id = $trans_poly_lookup->{$id};
        
        #get coords
        my $int_loc = $feature->first_child( 'Interval-loc' );
        my ($lend, $rend, $complement) = 
            ( $int_loc->att('startpos'), $int_loc->att('endpos'), $int_loc->att('complement') );

        my ($end5, $end3) = ($complement) ? ($lend, $rend) : ($rend, $lend);
        $logger->logdie("end5 and end3 not present $id") unless( $end5 && $end3 );

        my @gene_product_name = $feature->find_nodes( 'Attribute[@name="gene_product_name"]' );
        @gene_product_name = $feature->find_nodes( 'Attribute[@name="com_name"]' ) unless( @gene_product_name > 0 );

        my $com_name = "hypothetical protein";
        $com_name = $gene_product_name[0]->att('content') if( $gene_product_name[0] );

        push( @{$polypeptides}, {
            'feat_name' => $id,
            'end5'      => $end5,
            'end3'      => $end3,
            'com_name'  => $com_name } );
    }

    return $polypeptides;

}

sub get_length_from_seq_elem {
    my $seq_elem = shift;
    my $length;

    my $id = $seq_elem->att('id');
    my $sdi_elem = $seq_elem->first_child( 'Seq-data-import' );
    my ($identifier, $source) = ($sdi_elem->att('identifier'), $sdi_elem->att('source'));
    $identifier =~ s/^>//;

    $logger->logdie("Could not parse correct information from seq data import source for sequence $id")
        unless( $source );

    open( IN, "< $source") or $logger->logdie("Unable to open $source ($!)");

    
    my $seq = "";
    if( $identifier ) {
        my $flag = 0;
        while( my $line = <IN> ) {
            chomp $line;
            if( $line =~ /^>$identifier/ && !$flag) {
                $flag = 1;
            } elsif( $flag ) {
                $seq .= $line;
            } elsif( $flag && $line =~ /^>/ ) {
                last;
            }
        }
    } else {
        my $flag = 0;
        while( <IN> ) {
            chomp;
            if( /^>/ && !$flag ) {
                $flag = 1;
            } elsif( $flag ) {
                $seq .= $_;
            } elsif( $flag && /^>/ ) {
                last;
            }
        }
    }
    close(IN);

    $length = length( $seq );
    print "length of $id is $length\n";
    return $length;
    
}

sub check_parameters {
    my $opts = shift;

    if( $opts->{'input_bsml'} ) {
        open(IN, "< $opts->{'input_bsml'}") or 
            $logger->logdie("Unable to open $opts->{'input_bsml'} ($!)");

        chomp( @input_files = <IN> );
        @input_files = ( $opts->{'input_bsml'} ) if( $input_files[0] =~ /^\s*</ );
        close( IN );
    } else {
        $logger->logdie("Option input_bsml is required");
    }

    if( $opts->{'output_file'} ) {
        $output_file = $opts->{'output_file'};
    } else {
        $logger->logdie("Option output_file is required ($!)");
    }
    
}
