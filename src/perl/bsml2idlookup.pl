#!/usr/local/bin/perl

=head1  NAME 


=head1 SYNOPSIS

USAGE:  

=head1 OPTIONS

=head1 DESCRIPTION

=head1 INPUT

=head1 OUTPUT

=head1 CONTACT

=cut

use strict;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Workflow::Logger;
use XML::Twig;

use BSML::BsmlReader;
use BSML::BsmlParserTwig;
use BSML::BsmlParserSerialSearch;
use MLDBM "DB_File";

#######
## ubiquitous options parsing and logger creation
my %options = ();
#######
## ubiquitous options parsing and logger creation
my %options = ();
my $results = GetOptions (\%options, 
                            'bsml_input|i=s',
                            'bsml_list|s=s',
                            'bsml_file|b=s',        ## deprecated
                            'bsml_dir|d=s',         ## deprecated
                            'output_file|o=s',
                            'debug=s',
                            'log|l=s',
                            'output|o=s',
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

#######
## gather the list of files we're going to processes.  the user can pass file or
#  directory names using --bsml_input, which can be a single file, a list of 
#  filenames, a directory of files, or a list of directories of files (lists are
#  comma-separated.  it can even be a mixed list of file and dir names.  fancy!
#  a file containing a list of file/dir names should be passed with --bsml_list.  
my @files;
my @list_elements = ();

## did the user pass a list?  if so, get each file/dir names out of it, check it, and add it
if ( $options{bsml_list} ) {
    for my $list ( split(/,/, $options{bsml_list}) ) {
        $list =~ s/\s//g;
        
        if (-f $list && -s $list) {
            open (my $fh, $list) || $logger->logdie("Can't open file $list");
            for ( <$fh> ) {
                chomp;
                if (-e $_ && -s $_) {
                    push @list_elements, $_;
                } else {
                    $logger->warn("Error reading $_ from list $list") if ($logger->is_warn);
                }
            }
        } else {
            $logger->warn("Error reading $list") if ($logger->is_warn);
        }
    }
}

## loop through each input thing passed
for my $thing ( split(/,/, ($options{bsml_input} || '') ),
                split(/,/, ($options{bsml_file}  || '') ),  ## backwards compatibility
                split(/,/, ($options{bsml_dir}   || '') ),  ## backwards compatibility
                @list_elements) {
    $thing =~ s/\s//g;
    next unless ($thing);
    
    ## is this a directory?
    if (-d $thing) {
        opendir(my $dh, $thing) || $logger->warn("Unable to access $thing due to $!");

        for my $file (readdir $dh) {
            &add_file("$thing/$file");
        }

    ## else it's probably a file, make sure it exists and has size and add it, else warn
    } elsif (! &add_file($thing) ) {
        $logger->warn("Error reading $thing") if ($logger->is_warn);
    }
}

## die if no files found
if ( scalar @files == 0 ) {
    $logger->logdie("No files found");
}


#######
## parse out sequences from each file
my $mfh;
my $title;

my %lookup;
	
tie %lookup, 'MLDBM', $options{'output'} or $logger->logdie("Can't tie $options{'output'}");

for my $file ( @files ) {
    my $seqParser = new BSML::BsmlParserSerialSearch(
						     SequenceCallBack =>sub 
						     {
							 my $seqRef = shift;
							 my $seq_id = $seqRef->returnattr('id') || $seqRef->returnattr('identifier');
							 #support for deprecated link of ASSEMBLY attribute
							 $lookup{$seqRef->returnattr('id')} = $seqRef->returnBsmlAttr('ASSEMBLY');
							 if(defined $seqRef && defined $seqRef->returnBsmlFeatureTableListR->[0]){
							     my $features = $seqRef->returnBsmlFeatureTableListR->[0]->returnBsmlFeatureListR();
							     foreach my $feat (@$features){
								 $lookup{$feat->returnattr('id')} = $seq_id;
							     }
							 }
						     });
    $logger->debug("Parsing file $file") if ($logger->debug);
    $seqParser->parse( $file );
}

exit;	    

sub add_file {
    ## adds a file to the list of those to process, checking to make sure it
    #  exists and is populated.
    my $file = shift;
    
    ## only do .bsml files
    return 0 unless ( $file =~ /\.bsml$/ );
    
    if (-e $file && -s $file) {
        $logger->debug("Adding file $file for processing") if ($logger->is_debug);
        push @files, $file;
    } else {
        $logger->warn("Error reading file $file") if ($logger->is_warn);
        return 0;
    }
    
    return 1;
}

sub check_parameters {
    my ($options) = @_;
    
    ## they have to pass some form of input
    unless ($options{bsml_input} || $options{bsml_list} ||
            $options{bsml_file}  || $options{bsml_dir}) {
        $logger->logdie("You must specify input with --bsml_input --bsml_list --bsml_file or --bsml_dir");
    }

    ## output is required
    unless ( $options{output} ) {
        $logger->logdie("You must specify an output directory or file with --output");
    }

    if(0){
        pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});    
    }
}
