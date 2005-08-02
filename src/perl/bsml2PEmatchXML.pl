#!/usr/local/bin/perl

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

=back

=head1   DESCRIPTION

=cut


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use PEffect::PEffectXML;
use BSML::BsmlReader;
use BSML::BsmlParserTwig;

use File::Basename;
use File::Path;
use Pod::Usage;
use Workflow::Logger;
use MLDBM "DB_File";

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

foreach my $file (@files){
    my $parser = new BSML::BsmlParserTwig();
    my $reader = new BSML::BsmlReader();
    $logger->debug("Parsing entire file $file") if($logger->debug());
    $parser->parse( \$reader, $file );

    my $asmbllookup = build_asmbl_lookup($options{'asmbl_lookup'});

    foreach my $maln (@{$reader->returnMultipleAlignmentTables()}){
	my $mtable = $reader->readMultipleAlignmentTable($maln);
	
	my $sums = $mtable->{'AlignmentSummaries'};
	foreach my $sum (@{$sums}){
	    foreach my $seq1 (@{$sum->{'AlignedSequences'}}){
		foreach my $seq2 (@{$sum->{'AlignedSequences'}}){
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

    foreach my $aln (@{$reader->returnAllSeqPairAlignmentsListR()}){
	my $spaln = $reader->readSeqPairAlignment($aln);
	$logger->debug("Alingment between $spaln->{'refseq'} : $asmbllookup->{$spaln->{'refseq'}} and $spaln->{'compseq'} : $asmbllookup->{$spaln->{'compseq'}}") if($logger->is_debug());
	if($asmbllookup->{$spaln->{'refseq'}} ne $asmbllookup->{$spaln->{'compseq'}}){
	    my $seqpairruns = $spaln->{'seqPairRuns'};
	    my @sortedseqpairruns = sort { $a->{$options{'rankfield'}} <=> $b->{$options{'rankfield'}}  } @$seqpairruns;
	    my $bestrun;
	    if($options{'ranktype'} == 1){
		$bestrun = @sortedseqpairruns[0];
	    }
	    else{
		$bestrun = @sortedseqpairruns[$#sortedseqpairruns];
	    }
	    $logger->debug("Best run between $spaln->{'refseq'} $spaln->{'compseq'} of ".scalar(@sortedseqpairruns). " with rankfield: $bestrun->{$options{'rankfield'}} scorefield:$bestrun->{$options{'scorefield'}}") if($logger->is_debug());
	    my $score = $bestrun->{$options{'scorefield'}};
	    $logger->debug("Adding alignment between $spaln->{'refseq'} $spaln->{'compseq'} with score $score") if($logger->is_debug());
	    $pexml->addAlignment($spaln->{'refseq'}, $spaln->{'compseq'}, $score);
	}
    }
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
    eval {
	tie %lookup, 'MLDBM', $reader;
    };
    $logger->debug("Found ".scalar(keys %lookup)." keys") if($logger->is_debug());
    print STDERR "Found ".scalar(keys %lookup)." keys\n";
    return \%lookup;
}
