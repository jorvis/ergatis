#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use File::OpenFile qw(open_file);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use PFunc::AnnotationCollection;
use PFunc::Annotation;

##### GLOBALS ######
my @inputs;
my $output;

&check_options();

# open the output file
my $ofh = open_file( $output, 'out' );

# create the collection
my $an_col = new PFunc::AnnotationCollection('output' => $ofh);

my ($current_feat_id, $next_feat_id);

# cycle through each input file
foreach my $input ( @inputs ) {

    my $in = open_file( $input, 'in' );

    while( my $line = <$in> ) {
        chomp( $line );
        next if( $line =~ /^\s*$/ );

        # create the annotation object
        my $annotation = PFunc::Annotation::to_object( $line );
        
        $next_feat_id = $annotation->get_feature_id;

        $current_feat_id = $next_feat_id unless( defined( $current_feat_id ) );

        if( $next_feat_id ne $current_feat_id ) {
            # write current data to file
            $an_col->flush();
            $current_feat_id = $next_feat_id;
        }

        #And add the line to the collection as an annotation object
        $an_col->add_annotation( PFunc::Annotation::to_object($line) );
    }
    
    close($in);
}

if( $an_col->annotations_in_buffer ) {
    $an_col->flush();
}
print "Printed ".$an_col->total_annotations_printed()." annotations in the collection\n";

close( $ofh );


sub check_options {
    my %options;
    my $results = GetOptions (\%options,
                              'input_files|i=s',
                              'input_lists|n=s',
                              'output|o=s',
                              'help|h',
                              );

    if( $options{'input_files'} ) {
        my @input_files = split( /[,\s]+/, $options{'input_files'} );
        push(@inputs, @input_files);
    }

    if($options{'input_lists'}) {
        my @input_lists = split( /[,\s]+/, $options{'input_lists'} );
        foreach my $list( @input_lists ) {
            my $in = open_file( $list, 'in' );
            chomp( my @tmp = <$in> );
            push( @inputs, @tmp );
        }
    }

    die("Option --input_files or --input_lists is required. I need input.")
        unless( @inputs > 0 );

    if( $options{'output'} ) {
        $output = $options{'output'};
    } else {
        die("Option --output is required");
    }
}

