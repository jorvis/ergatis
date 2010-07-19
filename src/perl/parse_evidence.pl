#!/usr/bin/perl

=head1 NAME

parse_evidence.pl - Description

=head1 SYNOPSIS

 USAGE: parse_evidence.pl
       --input_file=/path/to/file.txt
       --input_list=/path/to/evidence.list
       --output=/path/to/output.tab
       --evidence_type=Pfunc
       --database_path=/path/to/db_dir
       --other_options="option1=value1 option2=value2..."
     [ --log=/path/to/file.log
       --debug=4
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>
    The input files to parse.  Can be a comma separated list of input files.

B<--input_list,-n>
    input list.  Can be a comma separated list of evidence lists. 

B<--output,-o>
    Output tab file

B<--evidence_type,-e>
    The type of evidence being parsed.

B<--database_path,-b>
    Path which holds various databases which the parsers may require.  Not required, but an individual parser
    will fail if it requires this path specifically.  For requirements, please see specific parser files in
    EvidenceParser namespace.

B<--other_options,-t>
    Allows the user to send extra options to the parser.  For specific details on options accepted
    and what they do, see documentation for the individual parser used. Should be a space separated list
    of key value pairs, quoted.

B<--log,-l>
    Logfile.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

    Will parse evidence files and create an annotation tab file for use in combining multiple 
    evidence sets into a final annotation.
 
=head1  INPUT
    The input files should be specific to the format defined by the individual parser.  See EvidenceParser.pm
    for details.

=head1 OUTPUT
    The output is a tab file based on the Annotation->to_string implementation.  See Annotation.pm for details.

=head1  CONTACT

    Kevin Galens
    kgalens@gmail.com

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::OpenFile qw(open_file);
use PFunc::EvidenceParser::ValidEvidenceTypes;

our %valid_evidence_types;
my @input_files;
my @bsml_files = ();
my $annotation_class;
my $append = 0;
my $parser;
my $ev_type;
my $output;

&check_options();

$parser->evidence( @input_files );
$parser->annotate_on( $annotation_class ) if( defined( $annotation_class ) );
$parser->bsml( @bsml_files ) if( @bsml_files );
$parser->parse();
$parser->to_file($output, $append);

sub check_options {
    
    my %options;
    my $results = GetOptions (\%options,
                              'input_file|i=s',
                              'input_list|n=s',
                              'output|o=s',
                              'evidence_type|e=s',
                              'bsml_feature_lookup_list|bl=s',
                              'bsml_feature_lookup_file|bf=s',
                              'append_mode|a=s',
                              'annotate_on|c=s',
                              'database_path|d=s',
                              'other_options|t=s',
                              'help|h',
                              );
   
    if( $options{'help'} ) {
        &_pod;
    }

    my @reqs = qw( output evidence_type );
    foreach my $req ( @reqs ) {
        die("Option $req is required") unless( $options{$req} );
    }

    $output = $options{'output'};

    $ev_type = $options{'evidence_type'};
    die("Evidence type $ev_type is not supported. See EvidenceParser::ValidEvidenceTypes.pm for ".
        "supported evidence types") 
        unless( exists( $valid_evidence_types{$ev_type} ) );

    ## grab all the input
    if( $options{'input_file'} ) {
        my @files = split(",", $options{'input_file'} );
        push( @input_files, @files );
    }
    if( $options{'input_list'} ) {
        my @lists = split(",", $options{'input_list'} );
        foreach my $list ( @lists ) {
            my $in = open_file( $list, 'in' );
            chomp( my @tmp = <$in> );
            close($in);
            push( @input_files, @tmp );
        }
    }

    if( @input_files == 0 ) {
        die("Please specify input using either the --input_file or ".
            "--input_list options");
    }

    if( $options{'bsml_feature_lookup_list'} ) {
        my @lists = split(",", $options{'bsml_feature_lookup_list'} );
        foreach my $list ( @lists ) {
            my $in = open_file( $list, 'in' );
            chomp( my @tmp = <$in> );
            close($in);
            push( @bsml_files, @tmp );
        }
    }
    if( $options{'bsml_feature_lookup_file'} ) {        
        my @files = split(",", $options{'bsml_feature_lookup_file'} );
        push( @bsml_files, @files );
    }
    
    if( $options{'annotate_on'} ) {
          die("Either --bsml_feature_lookup_file or --bsml_feature_lookup_list options are required ".
              "when using the --annotate_on option" ) if( @bsml_files == 0 );
          $annotation_class = $options{'annotate_on'};
    }

    my $db_path;
    if( $options{'database_path'} ) {
        $db_path = $options{'database_path'};
    }

    my @args = ();
    push(@args, "database_path => \"$db_path\"") if $db_path;

    ## parse any other user provided options and send them to the
    ## parser object
    if( $options{'other_options'} ) {
        my @pairs = split(/[,\s]+/, $options{'other_options'} );
        foreach my $p (@pairs) {
            my ($key, $value) = split(/[\s=]+/, $p );
            push(@args, "$key => \"$value\"");
        }
    }

    my $arg_string = join(", ", @args);
    my $command = "\$parser = new $valid_evidence_types{ $ev_type }( $arg_string )";
    print "$command\n";
    
    ## create the correct parser type
    eval("require $valid_evidence_types{$ev_type}") or 
        die("Could not require $valid_evidence_types{$ev_type}. Perhaps ".
            "EvidenceParser::ValidEvidenceTypes.pm is configured incorrectly [$@]");
    eval("\$parser = new $valid_evidence_types{ $ev_type }( $arg_string )") or 
        die("Could not create parser object ($valid_evidence_types{ $ev_type }). [$@]");
    

}
sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2} );
}

