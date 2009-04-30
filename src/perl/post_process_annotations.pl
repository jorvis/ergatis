#!/usr/bin/perl

=head1 NAME

post_process_annotations.pl - Contains functionality to refine final annotations.

=head1 SYNOPSIS

 USAGE: post_process_annotations.pl
       --input_file=/path/to/file.tab
       --output=/path/to/output.tab
     [ --log=/path/to/file.log
       --debug=4
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>
    The input tab file to parse.  Should be output from assign_annotations.pl. In general, a file containing
    Annotation.pm objects ->to_string representations

B<--output,-o>
    Output tab file.

B<--log,-l>
    Logfile.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

    Will process various annotations assigned to the features including:

    1. Name clean up as described in CommonNameAnnotation.pm
    2. TIGR role assignment by common name keywords
         DESCRIPTION or link to module doing the work
 
=head1  INPUT

    A file based on the Annotation->to_string implementation.  See Annotation.pm for details.  

=head1 OUTPUT

    Same format file, with the annotations changed where described in description

=head1  CONTACT

    Kevin Galens
    kgalens@som.umaryland.edu

=cut


use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::OpenFile qw(open_file);
use PFunc::Annotation;
use PFunc::CommonNameAnnotation qw(clean_common_name);
use PFunc::TIGRRolesAnnotation qw(assign_tigr_roles_by_keyword);
use Data::Dumper;

##### Globals #####
my $input_file;
my $output;
###################

&check_options();

my $in = open_file( $input_file, 'in' );
my $out = open_file( $output, 'out' );

while( my $line = <$in> ) {
    chomp( $line );

    my $annotation = PFunc::Annotation::to_object( $line );

    &post_process( $annotation );
}

close($in);
close($out);

sub post_process {
    my ($annotation) = @_;

    ## check for ec numbers in the common name and add them as annotations
    &check_for_ec_numbers( $annotation );

    ## format ec numbers
    &format_ec_numbers( $annotation );

    ## clean up the common name with the CommonNameAnnotation module
    my ($gene_product_name, $gpn_source, $gpn_source_type) = $annotation->get_gene_product_name;
    $annotation->set_gene_product_name( clean_common_name( $gene_product_name->[0] ),
                                        $gpn_source, $gpn_source_type );

    ## can we assign any TIGR roles based on common name keyword? It helps if this is after the
    ## clean up function
    my ($tigr_roles, $tr_source, $tr_source_type) = $annotation->get_TIGR_Role;
    my $new_tigr_roles;
    if( $new_tigr_roles = assign_tigr_roles_by_keyword( $gene_product_name->[0], $tigr_roles ) ) {
        $annotation->set_TIGR_Role( assign_tigr_roles_by_keyword( $gene_product_name->[0], $tigr_roles ),
                                    'by_keyword', 'by_keyword' );
    }

    print $out $annotation->to_string."\n";
}

sub format_ec_numbers {
    my ($annotation) = @_;
    my ($ecs, $source, $type) = $annotation->get_EC;
    my @new_ecs;
    my $change = 0;

    foreach my $ec ( @{$ecs} ) {
        if( $ec =~ /(([\d\-]+\.){3}[\d\-]+)/ ) {
            my $tmp_ec = $1;
            if( $ec ne $tmp_ec ) {
                $change = 1;
                $ec = $tmp_ec;
            }
            push( @new_ecs, $ec );
        }
    }

    if( $change ) {
        $annotation->set_EC(\@new_ecs, $source, $type);
    }
}

sub check_for_ec_numbers {
    my ($annotation) = @_;

    my ($gene_product_name_array, $source, $type) = $annotation->get_gene_product_name( 'gene_product_name' );
    my $gene_product_name = $gene_product_name_array->[0];

    if( !defined( $gene_product_name ) ) {
        print Dumper( $annotation );
        die("Don't have a gene_product_name");
        
    }

    my $feature_id = $annotation->get_feature_id;
    
    while( $gene_product_name =~ /ec\s+([\d\.\-]+)/gi ) {

        my $ec = $1;
        if( $annotation->has_annotation( 'EC' ) ) {
            $annotation->add_EC( $ec );
        } else {
            $annotation->set_EC( $ec, $source, $type );
        }
    }
}

sub check_options {
    my %options;
    my $results = GetOptions (\%options,
                              'input_file|i=s',
                              'output|o=s',
                              'log|l=s',
                              'debug|d=s',
                              'help|h',
                              );

    if( $options{'help'} ) {
        &_pod;
    }

    if( $options{'input_file'} ) {
        $input_file = $options{'input_file'};
    } else {
        die("Option --input_file is required");
    }

    if( $options{'output'} ) {
        $output = $options{'output'};
    } else {
        die("Option --output is required");
    }

}
