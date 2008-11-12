#!/usr/local/bin/perl
=head1  NAME 

create_chado2bsml_iterator_list.pl will take a new-line separated list of assembly identifiers and create an ergatis workflow iterator list file

=head1 SYNOPSIS

USAGE:  create_chado2bsml_iterator_list.pl --control_file=inputfile --output_file=outputfile

=head1 OPTIONS

=item *

B<--control_file> Input file containing new-line separated list of assembly identifiers

=item *

B<--output_file> Output iterator file that will be processed by the chado2bsml workflow component

=item *

B<--debug,-d> Debug level.  Use a large number to turn on verbose debugging. 

=item *

B<--logfile> Ergatis Workflow::Logger logfile

=item *

B<--help,-h> This help message

=head1   DESCRIPTION

=cut

#use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);


BEGIN {
    require '/usr/local/devel/ANNOTATION/cas/lib/site_perl/5.8.5/Workflow/Logger.pm';
    import Workflow::Logger;
}
use XML::Twig;
use File::Basename;
use Annotation::Util2;

umask(0000);

my %options = ();


my $results = GetOptions (\%options, 
                          'control_file=s', 
                          'bsmlfilelist=s', 
			  'output_file=s',
			  'logfile=s',
                          'debug=s', 
                          'help|h' ) || pod2usage();

my $logfile;

if (!defined($options{logfile})){

    $logfile = '/tmp/' . File::Basename::basename($0) . '.log';

} else {

    $logfile = $options{logfile};
}

my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Workflow::Logger::get_logger();

my $outfile;

if (defined($options{'output_file'})){

    $outfile = $options{output_file};

} else {

    $outfile = '/tmp/' . File::Basename::basename($0) . '.out';
}


my $asmblCtr=0;
my $asmblLookup={};
my $currentAsmblId;

if (defined($options{bsmlfilelist})){

    my $bsmlfilelist = $options{bsmlfilelist};

    if (!Annotation::Util2::checkInputFileStatus($bsmlfilelist)){
	$logger->logdie("Detected some problem with the bsml list file '$bsmlfilelist'");
    }

    my $contents = Annotation::Util2::getFileContentsArrayRef($bsmlfilelist);
    if (!defined($contents)){
	$logger->logdie("Could not retrieve contents of bsml list file '$bsmlfilelist'");
    }

    my $bsmlFileCtr=0;
    foreach my $bsmlfile (@{$contents}){
	parseBSMLFile($bsmlfile);
	$bsmlFileCtr++;


	if (! exists $asmblLookup->{$currentAsmblId}){
	    $asmblCtr++;
	    $asmblLookup->{$currentAsmblId} = $bsmlfile;
	} else {
	    $logger->logdie("Found duplicate assembly id '$currentAsmblId' ".
			    "while processing BSML file '$bsmlfile'.  The ".
			    "previous file that had the same id was ".
			    "'$asmblLookup->{$currentAsmblId}'");
	}		
	
    }

    print "Parsed '$bsmlFileCtr' BSML files\n";

} elsif (defined($options{control_file})){
    
    my $controlfile = $options{control_file};

    if (! Annotation::Util2::checkInputFileStatus($controlfile)){
	$logger->logdie("Detected some problem with control file '$controlfile'");
    }

    my $contents = Annotation::Util2::getFileContentsArrayRef($controlfile);
    if (!defined($contents)){
	$logger->logdie("Could not retrieve contents of control file '$controlfile'");
    }

    my $lineCtr=0;

    foreach my $line (@{$contents}){

	chomp $line;
	$lineCtr++;

	if ($line =~ /^\s*$/){
	    next; ## skip blank lines
	}

	if ($line =~ /^\#+/){
	    next;
	}

	$line =~ s/^\s+//; ## remove leading white spaces
	$line =~ s/\s+$//; ## remove trailing white spaces

	if (! exists $asmblLookup->{$line}){
	    $asmblCtr++;
	    $asmblLookup->{$line}++;
	} else {
	    $logger->warn("Found duplicate asmbl '$line' at line '$lineCtr' of control_file '$options{'control_file'}'");
	}		
    }

    print "Parsed '$lineCtr' lines in control file '$controlfile'\n";

} else {
    $logger->logdie("You must specify either --bsmlfilelist or --control_file");
}


if ($asmblCtr>0){

    open (OUTFILE, ">$outfile") || $logger->logdie("Could not open output_file '$outfile' in write mode: $!");

    print OUTFILE  '$;ASMBL$;' . "\n";

    foreach my $asmblId (sort keys %{$asmblLookup} ){
	print OUTFILE "$asmblId\n";
    }
    
    close OUTFILE;

} else {
    $logger->logdie("Did not find any assembly identifiers in control_file '$options{'control_file'}'");
}

print "$0 execution completed\n";
print "The output file is '$outfile'\n";
print "The log file is '$logfile'\n";
exit(0);
						     
#---------------------------------------------------------------------------------------------------------
#
#                           END OF MAIN  --  SUBROUTINES FOLLOW
#
#---------------------------------------------------------------------------------------------------------
 

sub parseBSMLFile {

    my ($bsmlfile) = @_;

    ## Function: parse input BSML file for sequence and gene information
    my $ifh;

    if (stat "$bsmlfile.gz"){
	$bsmlfile .= ".gz";
    }
    
    if ($bsmlfile =~ /\.(gz|gzip)$/) {
	open ($ifh, "<:gzip", $bsmlfile) || $logger->logdie("Could not open BSML file '$bsmlfile' ".
							    "in read mode:$!");
    }

    my $twig = new XML::Twig(   TwigRoots => { 'Sequence' => 1},
				TwigHandlers => { 'Sequence' => \&sequenceCallBack }
				);
    
    $twig->parsefile($bsmlfile);    # build the twig
}

sub sequenceCallBack {

    my ($twig, $bsmlSequence) = @_;

    my $class = $bsmlSequence->{'att'}->{'class'};
    if (!defined($class)){
	$logger->logdie("class was not defined");
    }
    if ($class ne 'assembly'){
	last;
    }

    if (! exists $bsmlSequence->{'att'}->{'id'}){
	$logger->logdie("id attribute does not exist for BSML <Sequence>");
    }

    $currentAsmblId =  $bsmlSequence->{'att'}->{'id'};
}
