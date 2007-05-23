#!/usr/local/bin/perl

use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use XML::Parser;
use Ergatis::Logger;

#Example invocation
#/usr/local/devel/ANNOTATION/ard/current/lib/ bsml2featurerelationships.pl --bsml_file /usr/local/annotation/PATHEMA/output_repository/legacy2bsml/2663_b_anthracis/gba_6615_assembly.prok.bsml --output=/tmp/test --output_order=gene,polypeptide,cds --add_assembly 
#output
#gba_6615_ORF00003_gene  gba_6615_ORF00003_CDS   gba_6615_ORF00003_polypeptide   gba_6615_assembly 
#gba_6615_ORF00004_gene  gba_6615_ORF00004_CDS   gba_6615_ORF00004_polypeptide   gba_6615_assembly
#gba_6615_ORF00005_gene  gba_6615_ORF00005_CDS   gba_6615_ORF00005_polypeptide   gba_6615_assembly
#...

#######
## ubiquitous options parsing and logger creation
my %options = ();

#
#Output format specifies class types
#For lookups, first class is key and subsequent classes are comma delimeted
#For tab output, format specifies the tab order
#The class ALL
my $results = GetOptions (\%options, 
                            'bsml_input|i=s',
                            'bsml_list|s=s',
                            'bsml_file|b=s',        ## deprecated
                            'bsml_dir|d=s',         ## deprecated
                            'output_order|t=s',     #eg. "polypeptide,cds". default alphabetical
                            'add_assembly|a=s',
                            'debug=s',
                            'log|l=s',
                            'output|o=s',
                            'fasta_list=s',
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

my $output_order;
if($options{'output_order'}){
    my $num=0;
    my @types = split(/,/,$options{'output_order'});
    foreach my $t (@types){
	$output_order->{lc($t)} = $num++;
    }
}

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

my %lookup;

my $currid;
my $currclass;
my $currgroup;
my @groups;

my $funcs = {'Feature'=>
		 sub {
		     my ($expat,$elt,%params) = @_;
		     $lookup{$params{'id'}} = $currid if($options{add_assembly});
		 },
	     'Sequence'=>
		 sub {
		     my ($expat,$elt,%params) = @_;
		     $currid = $params{'id'};
		     $currclass = $params{'class'};
		 },
	     'Feature-group'=>
		 sub {
		     my ($expat,$elt,%params) = @_;
		     if($currgroup){
			 push @groups,$currgroup;
		     }
		     $currgroup = {};
		 },
	     'Feature-group-member'=>
		 sub {
		     my ($expat,$elt,%params) = @_;
		     if(exists $options{output_order}){
			 if(exists $output_order->{lc($params{'feature-type'})}){
			     $currgroup->{$params{'feature-type'}}->{$params{'featref'}}++;
			 }
		     }
		     else{
			 $currgroup->{$params{'feature-type'}}->{$params{'featref'}}++;
		     }
		 }
	 };


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
    $x->parsefile( $file );
}

push( @groups, $currgroup);

print "Opened the output file for writing\n";
open OUTFILE,">$options{'output'}" or $logger->logdie("Can't open output file $options{'output'}");

print "There are ".@groups." groups\n";
foreach my $group (@groups){
    my @outline;
    my @types = keys %$group;
    foreach my $type (
		      #sort based on sort order
		      sort {
			  if(exists $options{'output_order'}){
			      $output_order->{$a} <=> $output_order->{$b};
			  }
			  else{
			      $a cmp $b;
			  }
		      }
		      @types) {
	#if multiple features for a type, concatenate with a ','
	push @outline,join(',',keys %{$group->{$type}});
    }
    #push assembly id if available
    my @vals = keys %{$group->{$types[0]}};
    if(exists $lookup{$vals[0]} && $options{add_assembly}){
	push @outline,$lookup{$vals[0]};
    }
    print "Printing to file\n";
    print OUTFILE join("\t",@outline),"\n";
}
close OUTFILE;

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
