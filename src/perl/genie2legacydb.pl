#!/usr/local/bin/perl

=head1  NAME 

genie2legacydb.pl - load genie BSML file into the legacy database schema

=head1 SYNOPSIS

USAGE:  genie2legacydb.pl
            --bsml_file=/path/to/somemapfile.bsml
            --list_file=/path/to/somefile.list
            --database=aa1
          [
            --debug=4
            --log=/path/to/somefile.log
          ]

=head1 OPTIONS

B<--bsml_file,-b> 
    Path to a single BSML genie result file.

B<--list_file,-i> 
    Path to a list of genie BSML output files.

B<--database,-d> 
    Sybase database name for the project in which these genie predictions
    will be loaded.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--help,-h> 
    This help message

=head1 DESCRIPTION

This script is used to load gene predictions from the genie program, in BSML format,
into one of the legacy databases.  Such BSML files are most likely generated from
the genie workflow component.

=head1 INPUT

The input is a BSML result file.  See the genie component documentation for the
elements that should be within this document.  

The name of each assembly is pulled from the BSML file names processed, which should
be in the format like:

    aa1.assembly.16296.0.genie.bsml
    
where 16296 is the asmbl id.  The only part of the name searched is ".assembly.(?)".

=head1 OUTPUT

This file generates no output on the file system unless a log file is specified.  It
does perform inserts into the project database passed.

Running this script causes a lot of unnecessary output on STDOUT.  This is done by the
EGC libraries used, so don't blame me.

=head1 CONTACT

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
BEGIN {
    require '/usr/local/devel/ANNOTATION/ard/chado-v1r5b1/lib/site_perl/5.8.5/Workflow/Logger.pm';
    import Workflow::Logger;
}
use XML::Twig;
use DBI;
use lib '/usr/local/devel/ANNOTATION/Euk_modules/bin';
use Gene_obj;
use Annot_prediction_loader;

#######
## ubiquitous options parsing and logger creation
my %options = ();
my $results = GetOptions (\%options, 
                            'bsml_file|b=s',
                            'list_file|i=s',
                            'database|d=s',
                            'debug=s',
                            'log|l=s',
                            'help|h') || pod2usage();

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = Workflow::Logger::get_logger();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

&check_parameters(\%options);

## create the database connection
my $dbh = DBI->connect("dbi:Sybase:server=SYBTIGR; packetSize=8092", 'egc', 'egcpwd', {RaiseError => 1, PrintError => 1});
my $result = $dbh->do("use $options{database}");

## gather the files we are going to process
my @files;
if (defined $options{bsml_file}) {
    $logger->debug("adding $options{bsml_file} to input list") if ($logger->is_debug);
    push @files, $options{bsml_file};
} elsif (defined $options{list_file}) {
    open(my $listfh, "<$options{list_file}") || $logger->logdie("can't read list file $options{list_file} : $!");

    while (<$listfh>) {
        chomp;
        next if (/^\s*$/);
        $logger->debug("adding $_ to input list") if ($logger->is_debug);
        push @files, $_;
    }
}

## because the coordinates are given within the Features, we need to store
##  them all so that the data is available when we start processing the 
##  Feature-groups.
my %coords;
my @genes;

## process each file
my $file;
for (@files) {
    $file = $_;
    $logger->debug("processing $file") if ($logger->is_debug);
    
    ## get the assembly id from the filename, such as:
    ##  aa1.assembly.(15588).0.genie.bsml
    my $asmbl_id;
    if ($file =~ /\.assembly\.(\d+)/) {
        $asmbl_id = $1;
    } else {
        $logger->logdie("couldn't get asmbl_id from file $file");
    }
    
    ## create the twig
    my $twig = XML::Twig->new(
                                twig_roots  => { 'Feature' => \&processFeature,
                                                 'Feature-group' => \&processFeatureGroup,
                                               }
                             );
    $twig->parsefile($file);
    
    ## load all the genes found in this file
    my $loader = new Annot_prediction_loader($dbh, $asmbl_id, 'genie');
    $loader->load_predictions(@genes) if (scalar(@genes));
    
    undef %coords;
    undef @genes;
}



exit(0);

sub check_parameters {
    my ($options) = @_;

    ## bsml_file or list_file must have been passed
    unless ( $options{bsml_file} || $options{list_file} ) {
        $logger->logdie("bsml_file or list_file must be passed");
    }

    ## check the map dir, if passed
    if ( $options{bsml_file} && ! -e $options{bsml_file} ) {
        $logger->logdie("bsml_file was passed but does not exist");
    }

    ## check map file
    if ( $options{list_file} && ! -e $options{list_file} ) {
        $logger->logdie("list_file either not passed or does not exist");
    }

    ## make sure a database was passed
    unless ($options{database}) {
        $logger->logdie("database must be passed!");
    }

    if(0){
        pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});    
    }
}


sub processFeature { 
    my ($twig, $feat) = @_;
    my $id = $feat->{att}->{id} || $logger->logdie("feature found with no id!");
    
    ## skip this if it isn't in the class 'exon'
    if ($feat->{att}->{class} eq 'exon') {
        ## the positions are given within the Interval-loc
        my $il = $feat->first_child('Interval-loc') || $logger->logdie("feature $id has no Interval-loc!");
        
        ## start and end pos both need to be incremented, because we number from
        #   0 in BSML but 1 in the legacydb
        my ($startpos, $endpos);
        if (defined $il->{att}->{startpos}) {
            $startpos = $il->{att}->{startpos};
            $startpos++;
        } else {
            $logger->logdie("feature $id in file $file has no startpos");
        }
            
        if (defined $il->{att}->{endpos}) {
            $endpos = $il->{att}->{endpos};
            $endpos++;
        } else {
            $logger->logdie("feature $id has no endpos");
        }
        
        ## store the coordinates for later lookup
        if ($il->{att}->{complement} == 0) {
            $coords{$id} = { end5 => $startpos, end3 => $endpos};
            
        } elsif ($il->{att}->{complement} == 1) {
            $coords{$id} = { end5 => $endpos, end3 => $startpos};
        
        } else {
            $logger->logdie("feature $id has an unknown complement value (" . $il->{att}->{complement} . ")!");
        }
        
    } else {
        $logger->debug("skipping feature $id because class != 'exon'") if ($logger->is_debug);
    }
}


sub processFeatureGroup {
    my ($twig, $fg) = @_;
    
    ## group-set and id of the Feature-group are ignored
    ## get each of the Feature-group-members in this Feature-group
    my @members = $fg->children('Feature-group-member');
    
    my %exon_coords;
    
    for my $member (@members) {
        ## skip this if it isn't an exon
        next unless ($member->{att}->{'feature-type'} eq 'exon');
        my $id = $member->{att}->{featref} || $logger->logdie("Feature-group-member in file $file has no featref!");
        
        $exon_coords{$coords{$id}->{end5}} = $coords{$id}->{end3};
    }

    ## now load this exon group into the database.  brian's method here
    
    # create a gene object
    my $gene = new Gene_obj();
    $gene->populate_gene_obj(\%exon_coords, \%exon_coords);
    
    push @genes, $gene;
}









