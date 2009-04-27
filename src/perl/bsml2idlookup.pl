#!/usr/bin/perl

#Creates a tied lookup of the form
#$lookup{$featureid} = $sequenceid
#from a BSML file like
#<Sequence id=$sequenceid>
# ....
#  <Feature id=$featureid>
#

use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use XML::Parser;
use Ergatis::Logger;
use DB_File;

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
                            'output_format=s', #tie[default],tab
                            'debug=s',
                            'log|l=s',
                            'output|o=s',
			    'fasta_list|f=s',
                            'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = Ergatis::Logger::get_logger();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

&check_parameters(\%options);

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

#######
## parse out sequences from each file

my %protein_fastas = ();
my %lookup;

if($options{output_format} eq 'tie'){
    $DB_BTREE->{'cachesize'} = 100000000;	
    unlink $options{'output'};
    my $dbtie = tie %lookup, 'DB_File', $options{'output'},  O_RDWR|O_CREAT, 0660, $DB_BTREE 
	or $logger->logdie("Can't tie $options{'output'}");
}
elsif($options{output_format} eq 'tab'){
    open OUTFILE,">$options{'output'}" or $logger->logdie("Can't open output file $options{'output'}");
}
else{
    $logger->logdie("Invalid output format $options{output_format}. Valid formats are 'tie' and 'tab'");
}

my $currid;
my $currclass;
my $is_tiedoutput = ($options{output_format} eq 'tie') ? 1 : 0;

my $funcs = {'Feature'=>
		 sub {
		     my ($expat,$elt,%params) = @_;
		     if($is_tiedoutput){
			 $lookup{$params{'id'}} = $currid;
		     }
		     else{
			 print OUTFILE "$params{'id'} $currid\n";
		     }
		 },
	     'Sequence'=>
		 sub {
		     my ($expat,$elt,%params) = @_;
		     $currid = $params{'id'};
		     $currclass = $params{'class'};
		 },
	     'Seq-data-import'=>
		 sub {
		     my ($expat,$elt,%params) = @_;
		     if($currclass eq 'polypeptide'){
			 ++$protein_fastas{$params{'source'}} if length $params{'source'};
		     }
		 }
	 };

my $parsedfiles=0;
for my $file ( @files ) {
    my $x = new XML::Parser(Handlers => 
			    {
				Start =>
				    sub {
					#$_[1] is the name of the element
					if(exists $funcs->{$_[1]}){
					    $funcs->{$_[1]}(@_);
					}
				    }
				}
			    );
    if (!(-e $file) && -e "$file.gz") {
	$file .= ".gz";
    }
    if(-e $file){
	print STDERR "Parsing $file\n";
    	my $ifh;
	if ($file =~ /\.(gz|gzip)$/) {
	    open ($ifh, "<:gzip", $file) || $logger->logdie("can't read input file $file: $!");
	} else {
	    open ($ifh, "<$file") || $logger->logdie("can't read input file $file: $!");
	}
	$x->parse( $ifh );
	close $ifh;
	$parsedfiles++;
    }
    else{
        $logger->logdie("Can't read jaccard bsml file $file");
    }
}

if($parsedfiles==0){
    $logger->logdie("No files found");
}

if ($options{fasta_list}){
    open FASTA_OUT, ">>$options{fasta_list}" or
	die "Error writing fasta list $options{fasta_list}: $!";
    foreach my $fasta (keys %protein_fastas) {
	print FASTA_OUT "$fasta\n";
    }
}

if($options{output_format} eq 'tab'){
    close OUTFILE;
}

sub add_file {
    ## adds a file to the list of those to process, checking to make sure it
    #  exists and is populated.
    my $file = shift;
    
    ## only do .bsml files
    if (!(-e $file) && -e "$file.gz") {
	$file .= ".gz";
    }
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

    if(!exists $options{output_format}){
	$options{output_format} = 'tie';
    }
    else{
	$logger->logdie("Invalid output format $options{output_format}. Valid formats are 'tie' and 'tab'") if($options{output_format} ne 'tie' && $options{output_format} ne 'tab');
    }

    if(0){
        pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});    
    }
}
