#!/usr/bin/perl

=head1 NAME

lipop2bsml.pl - Takes raw lipop output (short) and prints to bsml

=head1 SYNOPSIS

 USAGE: lipop2bsml.pl
       --input_file=/path/to/some/lipop.raw
       --input_fasta=/path/to/polypeptide.fsa
       --id_repository=/path/to/valid_id_repo/
       --project=db1
       --output=/path/to/lipop.bsml
     [ --analysis_id=lipoP_analysis
       --sourcename=/path/to/all_output_dir
       --log=/path/to/file.log
       --debug=4
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>
    Raw output from a LipoP run (using the -short option)

B<--input_fasta,-f>
    The input fasta file used to generate the lipop output

B<--id_repository,-d>
    A valid id repository (used for Ergatis::IdGenerator)

B<--project,-p>
    The project (database name).  Used in id generation

B<--output,-o>
    Path to the output bsml file

B<--analysis_id,-a>
    Optional. Default = lipoP_analysis

B<--sourcename,-s>
    Optional. Path to where the output for the analysis was located (raw output). 
    Default is parse the directory from the input_file (and traverse 2 dirs up)

B<--log,-l>
    Logfile.


B<--help,-h>
    Print this message

=head1  DESCRIPTION

    Example input:

    # gnmM04c.polypeptide.451723460.1 SpII score=34.6562 margin=19.0125 cleavage=19-20 Pos+2=S

    Example output:
    
    <Sequence id="gnmM04c.polypeptide.451723460.1" title="gnmM04c.polypeptide.451723460.1" ...>
      <Feature id="gnmM04c.lipoprotein_signal_peptide.3458954893.1" ...>
         <Interval-loc startpos="0" endpos="19" complement="0">
      </Feature>
      <Feature id="gnmM04c.cleavage_site.54890549.1" ...>
         <Site-loc sitepos="19" complement="0"></Site-loc>
      </Feature>
    </Sequence>
 
=head1  INPUT
    Describe the input

=head1 OUTPUT
    Describe the output

=head1  CONTACT

    Kevin Galens
    kgalens@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use BSML::BsmlBuilder;
use Ergatis::IdGenerator;
use File::OpenFile qw(open_file);
use Pod::Usage;

my $options = &check_options();
my $idgen;
my $analysis_id = "lipoP_analysis";
my $sourcename;

#grab the input
my $in = open_file( $options->{'input_file'}, 'in' );
chomp( my @lines = <$in> );
close($in);

#there should only be one line
if( @lines != 1 ) {
    die("There should only be one line in the LipoP raw output when the -short option is used");
}

my (undef, $id, $type, $score, $margin, $cleavage, $pos) = split( /\s+/, shift @lines );

my $doc = new BSML::BsmlBuilder;

my $seq = $doc->createAndAddSequence( $id, $id, '', 'aa', 'polypeptide' );

#validate the input fasta file to make sure it contains the sequence we are looking for
&validate_input_fasta( $options->{'input_fasta'}, $id );

#create the seq data import and analysis link
my $sdi = $doc->createAndAddSeqDataImport( $seq, 'fasta', $options->{'input_fasta'}, '', $id );
my $s_link = $doc->createAndAddLink( $seq, 'analysis', '#'.$analysis_id, 'input_of' );

#we are only encoding the cleavage site predictions (type either SpI or SpII)
if( $type eq 'SpI' || $type eq 'SpII' ) {
    
    #add a feature table
    my $ft = $doc->createAndAddFeatureTable( $seq );

    ## SpI is a signal peptide 
    my $signal_peptide_type;
    if( $type eq 'SpI' ) {
        $signal_peptide_type = 'signal_peptide';
    ## SpII is a lipoprotein signal peptide
    } elsif( $type eq 'SpII' ) {
        $signal_peptide_type = 'lipoprotein_signal_peptide';
    }

    #create the signal peptide feature, interval loc and analysis link
    my $sp_id = $idgen->next_id( 'type' => $signal_peptide_type, 'project' => $options->{'project'} );
    my $sp_feat = $doc->createAndAddFeature( $ft, $sp_id, $sp_id, $signal_peptide_type );
    my $endpos = $1 if( $cleavage =~ /(\d+)-(\d+)/ );
    my $sp_il = $doc->createAndAddIntervalLoc( $sp_feat, 0, $endpos, 0 );
    my $sp_link = $doc->createAndAddLink( $sp_feat, 'analysis', "#".$analysis_id, 'computed_by' );
    

    #create the cleavage site feature, site loc, and analysis link
    my $cs_id = $idgen->next_id( 'type' => 'cleavage_site', 'project' => $options->{'project'} );
    my $cs_feat = $doc->createAndAddFeature( $ft, $cs_id, $cs_id, 'cleavage_site' );
    my $cs_sl = $doc->createAndAddSiteLoc( $cs_feat, $endpos, 0 );
    my $cs_link = $doc->createAndAddLink( $cs_feat, 'analysis', '#'.$analysis_id, 'computed_by' );

}


#add the analysis
my $an = $doc->createAndAddAnalysis( 'id' => $analysis_id,
                                     'algorithm' => 'lipoP',
                                     'sourcename' => $sourcename, 
                                     'program' => 'lipoP',
                                     'programversion' => 'current');

#write the bsml file
$doc->write( $options->{'output'} );
print "Finished writing $options->{'output'}\n";

sub validate_input_fasta {
    my ($file, $id) = @_;
    my $in = open_file( $file, 'in' );
    my $found = 0;
    while( <$in> ) {
        if( /^>(\S+)/ ) {
            if( $1 eq $id ) {
                $found = 1;
                last;
            }
        }
    }
    close($in);
    
    die("Could not find sequence for $id in $file") unless( $found );
}

sub check_options {
    my %options;
    my $results = GetOptions (\%options,
                              'input_file|i=s',
                              'input_fasta|f=s',
                              'id_repository|d=s',
                              'project|p=s',
                              'output|o=s',
                              'analysis_id|a=s',
                              'sourcename|s=s',
                              'log|l=s',
                              'help|h',
                              );

   if( $options{'help'} ) {
       &_pod;
   }

    my @reqs = qw( input_file input_fasta id_repository project output );
    foreach my $req ( @reqs ) {
        die("Option $req is required") unless( $options{$req} );
    }

    # create IdGenerator.  Not setting pooling because there will only be a max of one id
    # of each type generated
    $idgen = new Ergatis::IdGenerator( 'id_repository' => $options{'id_repository'} );

    #if analysis_id was specified, replace the default
    $analysis_id = $options{'analysis_id'} if( $options{'analysis_id'} );

    #parse the sourcename if it was not specified
    if( $options{'sourcename'} ) {
        $sourcename = $options{'sourcename'};
    } else {
        $sourcename = $1 if( $options{'input_file'} =~ m|(.*)/[^/]+/[^/]+/.*| );
        die("Could not parse sourcename out of input_file [$options{'input_file'}]")
            unless( $sourcename );
    }

    return \%options;
}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
