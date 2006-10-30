#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

dummy.pl - do nothing

=head1 SYNOPSIS

USAGE:  dummy.pl --debug debug_level --log log_file

=head1 OPTIONS

=item *

B<--debug,-d> Debug level.  Use a large number to turn on verbose debugging. 

=item *

B<--log,-l> Log file

=item *

B<--help,-h> This help message

=head1   DESCRIPTION

=cut


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use PEffect::PEffectXML;

use File::Basename;
use File::Path;
use Pod::Usage;
BEGIN {
use Workflow::Logger;
use BSML::BsmlReader;
use BSML::BsmlParserSerialSearch;
}
use DB_File;

my %options = ();
my $results = GetOptions (\%options, 
			  'bsml_file|b=s', 
			  'bsml_list|s=s', 
			  'bsml_dir|d=s',
			  'asmbl_lookup|a=s',
			  'output|o=s', 
			  'scorefield|c=s',
			  'rankfield|r=s',
			  'ranktype|t=s',
			  'log|l=s',
			  'debug=s',
			  'help|h') || pod2usage();


my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Workflow::Logger::get_logger();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDERR} );
}

$options{'scorefield'} = 'runscore' if($options{'scorefield'} eq "");
$options{'rankfield'} = 'runprob' if($options{'rankfield'} eq "");
$options{'rankfield'} = 0 if($options{'ranktype'} eq "");

&check_parameters(\%options);

my @files;
if ($options{'bsml_file'} ne ""){
    if((-s $options{'bsml_file'})) {
	push @files,$options{'bsml_file'};
    }
    else{
	$logger->warn("Error reading file $options{'bsml_file'}");
    }
}
if ($options{'bsml_list'}){
    open FILE, "$options{'bsml_list'}";
    while (my $filename=<FILE>){
	chomp $filename;
	if(-s $filename) {
	    push @files,$filename;
	}
        else{
	    $logger->warn("Error reading file $options{'bsml_file'}");
	}
    }
}

if ($options{'bsml_dir'} && (-r $options{'bsml_dir'})) {
    opendir(DIR, $options{'bsml_dir'}) or $logger->warn("Unable to access $options{'bsml_dir'} due to $!");
    while( my $filename = readdir(DIR)) {
	if($filename =~ /(.+)\.bsml$/) {
	    if(-s "$options{'bsml_dir'}/$filename"){
		$logger->debug("Adding file $options{'bsml_dir'}/$filename to list of files to process") if($logger->is_debug());
		push (@files, "$options{'bsml_dir'}/$filename");
	    }
	    else{
		$logger->warn("Error reading file $options{'bsml_dir'}/$filename");
	    }
        }
    }
}
else{
    $logger->warn("Error reading directory $options{'bsml_dir'}");
}


my $pexml =  PEffect::PEffectXML->new();

my $asmbllookup = build_asmbl_lookup($options{'asmbl_lookup'});

foreach my $file (@files){
    my $parser = new BSML::BsmlParserSerialSearch( MultipleAlignmentCallBack => 
						   sub{
						       my $aln = shift;
						       multipleAlignmentHandler($aln,$pexml);
						   },
						   AlignmentCallBack => 
						   sub{
						       my $aln = shift;
						       alignmentHandler($aln,$asmbllookup,$pexml);
						   });

    $logger->debug("Parsing entire file $file") if($logger->debug());

    $parser->parse( $file );

}

my($oref);
$pexml->outputXML(\$oref);
open FILE, "+>$options{'output'}" or $logger->logdie("Can't open file $options{'output'}");
print FILE $oref,"\n";
close FILE;

sub check_parameters{
    my ($options) = @_;
    
    if(0){
	pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});    
    }
}

sub build_asmbl_lookup{
    my $reader = shift;
    my %lookup;
    $logger->debug("Reading lookup $reader") if($logger->is_debug());
    tie %lookup, 'DB_File', $reader, undef, undef, $DB_BTREE or $logger->logdie("Can't tie $reader");
    $logger->debug("Found ".scalar(keys %lookup)." keys") if($logger->is_debug());
    return \%lookup;
}

sub multipleAlignmentHandler
{
    my $aln = shift;
    my $pexml = shift;
    my $bsml_reader = new BSML::BsmlReader();

    my $maln_ref = $bsml_reader->readMultipleAlignmentTable($aln);

    foreach my $alnSum ( @{ $maln_ref->{'AlignmentSummaries'} } )
    {
	my $seqCount = 0;
	foreach my $seq1 ( @{ $alnSum->{'AlignedSequences'} } )
	{
	    foreach my $seq2 ( @{ $alnSum->{'AlignedSequences'} } ){
		my $name1 = $seq1->{'name'};
		my $name2 = $seq2->{'name'};
		if($name1 ne $name2){
		    $name1 =~ s/:\d+$//;
		    $name2 =~ s/:\d+$//;
		    my $score = 100;
		    $logger->debug("Adding alignment between $name2 and $name2 with score $score") if($logger->is_debug());
		    $pexml->addAlignment($name2, $name1, $score);
		}
	    }
	}
    }
}

sub alignmentHandler
{
    my $aln = shift;
    my $asmbllookup = shift;
    my $pexml = shift;

    my $refseq = $aln->returnattr( 'refseq' );
    my $compseq = $aln->returnattr( 'compseq' );

    if($asmbllookup->{$refseq} ne $asmbllookup->{$compseq}){
	$logger->debug("Alignment between $refseq : $asmbllookup->{$refseq} and $compseq : $asmbllookup->{$compseq}") if($logger->is_debug());

	my $seqpairruns = $aln->returnBsmlSeqPairRunListR;

	my @sortedseqpairruns;
	if($seqpairruns->[0]->{'BsmlAttr'}->{$options{'rankfield'}}->[0]){
	    @sortedseqpairruns = sort { $a->{'BsmlAttr'}->{$options{'rankfield'}}->[0] <=> $b->{'BsmlAttr'}->{$options{'rankfield'}}->[0]  } @$seqpairruns;
	}
	elsif($seqpairruns->[0]->returnattr($options{'rankfield'})){
	    @sortedseqpairruns = sort { $a->returnattr($options{'rankfield'}) <=> $b->returnattr($options{'rankfield'})  } @$seqpairruns;
	}
	else{
	    $logger->logdie("Bad rankfield $options{'rankfield'}");
	    die;
	}
	    
	my $bestrun;
	if($options{'ranktype'} == 1){
	    $bestrun = @sortedseqpairruns[0];
	}
	else{
	    $bestrun = @sortedseqpairruns[$#sortedseqpairruns];
	}
	

	my $score;
	if($bestrun->{'BsmlAttr'}->{$options{'scorefield'}}->[0]){
	   $score = $bestrun->{'BsmlAttr'}->{$options{'scorefield'}}->[0];
	}
	elsif($bestrun->returnattr($options{'scorefield'})){
	   $score = $bestrun->returnattr($options{'scorefield'});
	}
	else{
	    $logger->logdie("Bad scorefield $options{'scorefield'}");
	    die;
	}   
	 
	$logger->debug("Best run between $refseq $compseq of ".scalar(@sortedseqpairruns). " with rankfield: $options{'rankfield'} scorefield:$options{'scorefield'} bestscore: $score") if($logger->is_debug());

	$logger->debug("Adding alignment between $refseq $compseq with score $score") if($logger->is_debug());
	$pexml->addAlignment($refseq, $compseq, $score);
    }
}
