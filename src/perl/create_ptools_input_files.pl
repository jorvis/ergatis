#!/usr/bin/perl

=head1 NAME

create_ptools_input_files.pl 

=head1 SYNOPSIS

  USAGE: prepare_for_auto_annotate.pl  
    --gbk_dir=/path/to/dir/with/molecule.gbk
    --organism_name='Genus species strain'
    --database=db123
    --taxon_id=12334
    --output_dir=/path/to/output/dir
    [ --help ]

=head1 OPTIONS

B<--gbk_dir,-g>
    Path to a directory containing the organism's gbk files

B<--organism_name,-n>
    Name of the organism. Must contain Genus and species, strain is optional

B<--database,-d>
    Name of the database

B<--taxon_id,-t>
    NCBI taxon id.  If unknown, enter closest related.

B<--output_dir,-o>
    Directory to print the two output files:
    organism-params.dat
    genetic-elements.dat

B<--log,-l>
    Logfile name

B<--help,-h>
    Print this message

=head1  DESCRIPTION

    Will print 2 files based on input information:
    
    organism-params.dat:
    ID      DB123
    STORAGE FILE
    NAME    Genus species strain
    ABBREV-NAME     G. species 
    DOMAIN  TAX-12334
    RANK    |strain|
    CODON-TABLE     11
    MITO-CODON-TABLE        0
    DBNAME  Db123Cyc
    NCBI-TAXON-ID   12334

    genetic-elements.dat
    ID      CHROM-1
    NAME    db123.assembly.1
    TYPE    :CHRSM
    CIRCULAR?       N
    ANNOT-FILE      db123.assembly.1.gbk
    SEQ-FILE        db123.assembly.1.fsa
    //

=head1  INPUT

    Input GBK files should be named according to molecule and contain only one molecule per file.
    The fasta file used for input to pathway-tools should be named the same as the gbk files found.

=head1 OUTPUT

    See above

=head1  CONTACT

    Kevin Galens
    kgalens@tigr.org

=cut

##################### Modules #############################################
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Ergatis::Logger;
use Pod::Usage;
use Data::Dumper;
############################################################################
$|++;


########################### GLOBALS AND CONSTANTS ##################################
my @gbk_file_names;
my $db;
my $taxon_id;
my $org_name;
my $output_dir;
####################################################################################

###################### OPTION GATHERING AND LOG SETUP ##############################
my %options = ();
my $results = GetOptions (\%options, 
                          'gbk_dir|g=s',
                          'organism_name|n=s',
                          'database|d=s',
                          'taxon_id|t=s',
                          'output_dir|o=s',
                          'log|l=s',
                          'help|h') || &_pod;

#Setup the logger
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

#Check the options
&checkParameters(\%options);
###################################################################################

############################## MAIN ###############################################
my $org_param_file = $output_dir."/organism-params.dat";
&make_organism_params_file( $org_param_file, $db, $org_name, $taxon_id );

my $gen_ele_file = $output_dir."/genetic-elements.dat";
&make_genetic_elements_file( $gen_ele_file, \@gbk_file_names );
###################################################################################


############################## SUB ROUTINES #######################################

#Name:  make_organism_params_file
#Desc:  Creates the organism-params.dat file required for ptools batch input
#Args:  $file: file name to write to
#         $db: db name
#         $org_name: Organism name
#         $taxon_id: NCBI Taxon ID
#Rets:  True on success
sub make_organism_params_file {
    my ($file, $db, $org_name, $taxon_id) = @_;

    my @org_words = split(/\s+/, $org_name);
    my $genus = shift( @org_words );
    my $first_initial = $1 if( $genus =~ /(\w)/ );
    my $org_abbrev = uc($first_initial).". $org_words[0]";
    
    open( OUT, ">$file" ) or die("Can't open $file for writing ($!)");
    print OUT "ID\t".uc($db)."\n";
    print OUT "STORAGE\tFILE\n";
    print OUT "NAME\t$org_name\n";
    print OUT "ABBREV-NAME\t$org_abbrev\n";
    print OUT "DOMAIN\tTAX-$taxon_id\n";
    print OUT "RANK\t|strain|\n";
    print OUT "CODON-TABLE\t11\n";
    print OUT "MITO-CODON-TABLE\t0\n";
    print OUT "DBNAME\t".ucfirst(lc($db))."Cyc\n";
    print OUT "NCBI-TAXON-ID\t$taxon_id\n";
    close(OUT);
    
    return 1;
}

#Name:  make_genetic_elements_file
#Desc:  Creates the genetic-elements.dat file
#Args:  $file: file to be written to
#         $gbk_names: array ref of gbk_names
#Rets:  True on success
sub make_genetic_elements_file {
    my ($file, $gbk_names) = @_;

    open(OUT, "> $file") or die("Can't open $file for writing ($!)");
    for( my $i = 0; $i < @{$gbk_names}; $i++ ) {
        my $basename = $1 if( $gbk_names->[$i] =~ m|/?([^/]+)\.gbk| );

        print OUT "ID\tCHROM-".($i+1)."\n";
        print OUT "NAME\t$basename\n";
        print OUT "TYPE\t:CHRSM\n";
        print OUT "CIRCULAR?\tN\n";
        print OUT "ANNOT-FILE\t$basename.gbk\n";
        print OUT "SEQ-FILE\t$basename.fsa\n";
        print OUT "//\n";

    }
    close(OUT);

    return 1;
}

#Name: _get_gbk_file_names
#Desc: Will grab the file names for all gbk files in a directory
#Args: Directory name
#Rets: array of file names
sub _get_gbk_file_names {
    my ($dir) = @_;
    opendir( DIR, $dir ) or die("Could not open $dir");
    my @filenames = grep { /.*gbk$/ } readdir( DIR );
    closedir( DIR );
    return @filenames;
}   

#Name:    checkParameters
#Desc:    Will check the input options and die if something isn't good.
#Params:  $opts: hashref, contains the options (from getopts).
#Returns: Nothing
sub checkParameters {
    my $opts = shift;
    
    foreach my $required_parameter( qw(gbk_dir database taxon_id organism_name output_dir) ) {
        die("Option $required_parameter required") if( !$opts->{$required_parameter} );
    }

    @gbk_file_names = &_get_gbk_file_names( $opts->{'gbk_dir'} );
    $db = $opts->{'database'};
    $taxon_id = $opts->{'taxon_id'};
    $org_name = $opts->{'organism_name'};
    $output_dir = $opts->{'output_dir'};
}
##### EOF #####
