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
use BSML::BsmlParserSerialSearch;

use File::Basename;
use File::Path;
use Pod::Usage;
use Workflow::Logger;

my %options = ();
my $results = GetOptions (\%options, 
			  'bsml_file|b=s', 
			  'bsml_list|s=s', 
			  'bsml_dir|d=s', 
			  'output|o=s', 
			  'gene_count_cutoff|g=s', 
			  'gene_length_cutoff|c=s', 
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

$options{'gene_count_cutoff'} = 2 if($options{'gene_count_cutoff'} eq "");

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
    my $seqParser = new BSML::BsmlParserSerialSearch(
						     SequenceCallBack =>sub 
						     {
							 my $seqRef = shift; 
							 my $order = 0; 
							 if($seqRef->returnattr('molecule') eq "dna"){
							     my $sorted_genes = get_sorted_gene_position($seqRef);
							     if(scalar(@$sorted_genes) >= $options{'gene_count_cutoff'}){
								 foreach my $gene (@$sorted_genes) {
								     my $attrref={};
								     my $feat_name = $gene->{'feat_name'};
								     $attrref->{'length'} = $gene->{'length'}; 
								     $attrref->{'orient'} = $gene->{'orient'}; 
								     $attrref->{'coord'}  = $gene->{'coord'}; 
								     if($gene->{'length'} >= $options{'gene_length_cutoff'}){
									 $pexml->addFeature($feat_name, $attrref, $seqRef->returnattr('id'), $order);
								     }
								     else{
									 $logger->warn("$feat_name length $gene->{'length'} is less than $options{'gene_length_cutoff'}!!! skipping...");
								     }	
								     $order++;
								 }
							     }
							 }
						     });
    $seqParser->parse($file);
}

my($oref);
$pexml->outputXML(\$oref);
open FILE, "+>$options{'output'}" or $logger->logdie("Can't open file $options{'output'}");
print FILE $oref,"\n";
close FILE;

sub get_sorted_gene_position {

    my $seq = shift;
    my $reader = new BSML::BsmlReader();
    my $feats = $reader->readFeatures($seq);

    $logger->debug("Reading sequence $seq->{'attr'}->{'id'}") if($logger->is_debug());

    my $gene_pos = [];

    foreach my $feat (@$feats){
	if(lc($feat->{'class'}) eq "cds"){
	    $logger->debug("Reading cds $feat->{'id'}") if($logger->is_debug());
	    my $loc = $feat->{'locations'}->[0];
	    my $coord = $loc->{'startpos'} > $loc->{'endpos'} ? $loc->{'startpos'} : $loc->{'endpos'};
	    $logger->debug("Coord start:$loc->{'startpos'} end:$loc->{'endpos'}") if($logger->is_debug());
	    my $length = abs($loc->{'endpos'} - $loc->{'startpos'});
	    my $complement = $loc->{'complement'} == 1 ? '+' : '-'; 
	    $logger->debug("Storing position:$coord length:$length complement:$complement for $feat->{'id'}") if($logger->is_debug());
	    my $proteinId;
	    foreach my $link (@{$feat->{'bsmllinks'}}){
		# a link of type 'SEQ' links to the protein sequence object
		if( $link->{'rel'} eq 'SEQ' ){
		    $proteinId = $link->{'href'};
		    $proteinId =~ s/#//;
		    $logger->debug("Storing protein id $proteinId associated with cds $feat->{'id'}") if($logger->is_debug());
		    if($proteinId eq ""){
			$logger->logdie("Can't find associted protein with cds $feat->{'id'}");
		    }
		}
	    }
	    push @$gene_pos, { 'feat_name' => $proteinId, 'length' => $length, 'orient' => $complement, 'coord' => $coord } if($proteinId); 
	}
	elsif(lc($feat->{'class'}) eq "protein"){
	    $logger->debug("Reading protein $feat->{'id'}") if($logger->is_debug());
	    my $loc = $feat->{'locations'}->[0];
	    my $coord = $loc->{'startpos'} > $loc->{'endpos'} ? $loc->{'startpos'} : $loc->{'endpos'};
	    $logger->debug("Coord start:$loc->{'startpos'} end:$loc->{'endpos'}") if($logger->is_debug());
	    my $length = abs($loc->{'endpos'} - $loc->{'startpos'});
	    my $complement = $loc->{'complement'} == 1 ? '+' : '-'; 
	    $logger->debug("Storing position:$coord length:$length complement:$complement for $feat->{'id'}") if($logger->is_debug());
	    my $proteinId = $feat->returnattr('id');
	    push @$gene_pos, { 'feat_name' => $proteinId, 'length' => $length, 'orient' => $complement, 'coord' => $coord } if($proteinId); 
	}
	    
    }

    my @sorted_ref = sort { $a->{'coord'} <=> $b->{'coord'}} @$gene_pos;
    
    return (\@sorted_ref);

}

sub check_parameters{
    my ($options) = @_;
    
    if(0){
	pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});    
    }
}
