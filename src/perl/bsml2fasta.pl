#!/usr/local/bin/perl

=head1  NAME 

bsml2fasta.pl  - convert BSML features to fasta format

=head1 SYNOPSIS

USAGE:  bsml2fasta.pl -d bsml_dir -o output.fsa -f bsml_file -t protein

=head1 OPTIONS

=over 4

=item *

B<--bsml_dir,-d> [required] Directory containing BSML documents.  The
script will look for files with a .bsml suffix.

=item *

B<--bsml_file,-b> [required] BSML file

=item *

B<--output,-o> [required] Output fasta file or file prefix.  If
--format option is set to 'single', this option will provide a prefix
for the output file.  If --format option is et to 'multi', this option
will provide the output filename.

=item *

B<--type,-t> [required] Type of sequence to parse from the BSML
file(s).  Either 'assembly','protein','cds'.

=item *

B<--format,-f> [optional:default=multi] File format to produce. Either
'single' fasta file or 'multi' fasta file.  For 'single' fasta files,
the --output parameter will provide a prefix for the filename.  The
first word of the header line will be used as the rest of the filename
and end in a .fsa extension.

=item *

B<--debug,-d> Debug level.  Use a large number to turn on verbose debugging. 

=item *

B<--log,-l> Log file

=item *

B<--help,-h> This help message

=back

=head1   DESCRIPTION

Samples:

1. convert directory of bsml files to a multi-fasta file

    bsml2fasta -d bsml_dir -o allseqs.pep -t protein -f 'multi'

NOTE:  

Calling the script name with NO flags/options or --help will display the syntax requirement.

=cut


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
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
			  'format|f=s', 
			  'type|t=s', 
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

&check_parameters(\%options);

my @files;
if ($options{'bsml_file'} ne ""){
    my @inputfiles = split(/,/,$options{'bsml_file'});
    foreach my $file (@inputfiles){
	if((-s $file)) {
	    $logger->debug("Adding file $file for processing") if($logger->is_debug());
	    push @files,$file;
	}
	else{
	    $logger->warn("Error reading file $file");
	}
    }
}
if ($options{'bsml_list'}){
    open FILE, "$options{'bsml_list'}" or $logger->logdie("Can't open file $options{'bsml_list'}");
    while (my $filename=<FILE>){
	chomp $filename;
	if(-s $filename) {
	    $logger->debug("Adding file $filename for processing") if($logger->is_debug());
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
		$logger->debug("Adding file $options{'bsml_dir'}/$filename for processing") if($logger->is_debug());
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



foreach my $file (@files){

    
    if($options{'type'} eq "cds"){
	my $parser = new BSML::BsmlParserTwig();
	my $reader = new BSML::BsmlReader();
	$logger->debug("Parsing entire file $file") if($logger->debug());
	$parser->parse( \$reader, $file );
	&write_cds_fasta($reader,$options{'output'},$options{'format'});
    }
    else{
	my $lookup;
	$logger->debug("Parsing file $file with serial parse") if($logger->is_debug());
	my $reader = new BSML::BsmlReader();
	if($options{'type'} eq "protein"){
	    my $seqParser = new BSML::BsmlParserSerialSearch( ReadFeatureTables => 0,
							      SequenceCallBack =>sub 
							      {
								  my $seqRef = shift;
								  my $id = $seqRef->returnattr( 'id' );
								  my $type = $seqRef->returnattr( 'molecule' );
								  if($type eq "aa"){
								      $logger->debug("Retrieving sequence $id from reader")  if($logger->is_debug());
								      my $seqdat = $reader->subSequence($seqRef, -1, -1, '0');
								      $logger->debug("Storing sequence $id in lookup")  if($logger->is_debug());
								      $lookup->{$type}->{$id} = $seqdat;
								  }
							      }); 
	    $logger->debug("Beginning parse")  if($logger->is_debug());
	    $seqParser->parse( $file );
	    $logger->debug("Ending parse")  if($logger->is_debug());
	    &write_protein_fasta($lookup,$options{'output'},$options{'format'});
	}
	elsif($options{'type'} eq "assembly"){
	    my $seqParser = new BSML::BsmlParserSerialSearch( ReadFeatureTables => 0,
							      SequenceCallBack =>sub 
							  {
							      my $seqRef = shift;
							      my $id = $seqRef->returnattr( 'id' );
							      my $type = $seqRef->returnattr( 'molecule' );
							      if($type eq "dna"){
								  $logger->debug("Retrieving sequence $id from reader")  if($logger->is_debug());
								  my $seqdat = $reader->subSequence($seqRef, -1, -1, '0');
								  $logger->debug("Storing sequence $id in lookup")  if($logger->is_debug());
								  $lookup->{$type}->{$id} = $seqdat;
							      }
							  }); 
	    	
	    $logger->debug("Beginning parse")  if($logger->is_debug());
	    $seqParser->parse( $file );
	    $logger->debug("Ending parse")  if($logger->is_debug());
	    &write_assembly_fasta($lookup,$options{'output'},$options{'format'});
	}
	else{
	    $logger->logdie("Unsupported type $options{'type'}. Try protein, assembly, or cds.");
	}
    }
}

sub write_protein_fasta{
    my($lookup,$output,$format) = @_;
    if($format eq "multi"){
	my $multifile = "$output";
	open(MULTI_FILE, ">>$output") || $logger->logdie("Unable to write to file $output due to $!");
	$logger->debug("Writing output to multi-fasta file $output") if($logger->is_debug());
    }
    
    my $seqs = $lookup->{'aa'};
    if(scalar (keys %$seqs)){
	foreach my $seq_id (keys %$seqs) {
	    $logger->debug("Processing $seq_id") if($logger->is_debug());
	    my $raw_seq = $seqs->{$seq_id};
	    if($format eq "single"){
		my $outputfile = "$output${seq_id}.fsa";
		open(SINGLE_FILE, ">$outputfile") || $logger->logdie("Unable to write to file $outputfile due to $!");
		$logger->debug("Writing output to single-fasta file $outputfile") if($logger->is_debug());
		&fasta_out($seq_id,$raw_seq,\*SINGLE_FILE);
		close SINGLE_FILE;
	    }
	    elsif($format eq "multi"){
		&fasta_out($seq_id,$raw_seq,\*MULTI_FILE);
	    }	
	}
    }
    else{
	if($format eq "multi"){
	    print MULTI_FILE ">\n";
	}
    }
    if($format eq "multi"){
	close MULTI_FILE;
    }
}

sub write_assembly_fasta{
    my($lookup,$output,$format) = @_;
    if($format eq "multi"){
	my $multifile = "$output";
	open(MULTI_FILE, ">>$output") || $logger->logdie("Unable to write to file $output due to $!");
	$logger->debug("Writing output to multi-fasta file $output") if($logger->is_debug());
    }   
    my $seqs = $lookup->{'dna'};
    foreach my $seq_id (keys %$seqs) {
	$logger->debug("Processing $seq_id") if($logger->is_debug());
	my $raw_seq = $seqs->{$seq_id};
	if($format eq "single"){
	    my $outputfile = "$output${seq_id}.fsa";
	    open(SINGLE_FILE, ">$outputfile") || $logger->logdir("Unable to write to file $outputfile due to $!");
	    $logger->debug("Writing output to single-fasta file $outputfile") if($logger->is_debug());
	    &fasta_out($seq_id,$raw_seq,\*SINGLE_FILE);
	    close SINGLE_FILE;
	}
	elsif($format eq "multi"){
	    &fasta_out($seq_id,$raw_seq,\*MULTI_FILE);
	}	
    }
    if($format eq "multi"){
	close MULTI_FILE;
    }
}

		
sub write_cds_fasta{
    my($reader,$output,$format) = @_;
    
    if($format eq "multi"){
	my $multifile = "$output";
	open(MULTI_FILE, ">>$output") || $logger->logdie("Unable to write to file $output due to $!");
	$logger->debug("Writing output to multi-fasta file $output") if($logger->is_debug());
    }
    
    my $seq_list = $reader->returnAllSequences();
    foreach my $seq_obj (@$seq_list) {
	if($seq_obj->returnattr('molecule') eq 'dna') {
	    my $cds_seq = $reader->get_all_cds_dna($seq_obj->returnattr('id')) ;
	    while( my ($gene, $raw_seq) = each %$cds_seq ) {
#why don't we use the cdsid?
		my $protein_seq_id = $reader->cdsIdtoProteinSeqId($gene);
		$logger->debug("Processing $protein_seq_id") if($logger->is_debug());
		if($format eq "single"){
		    my $gene_file = "$output${protein_seq_id}.fsa";
		    open(SINGLE_FILE, ">$gene_file") || $logger->logdie("Unable to write to $gene_file due to $!");
		    $logger->debug("Writing output to single-fasta file $gene_file") if($logger->is_debug());
		    &fasta_out($protein_seq_id,$raw_seq,\*SINGLE_FILE);
		    close SINGLE_FILE;
		}
		elsif($format eq "multi"){
		    &fasta_out($protein_seq_id,$raw_seq,\*MULTI_FILE);
		}
	    }
	}
    }
    if($format eq "multi"){
	close MULTI_FILE;
    }

}
sub fasta_out {
#This subroutine takes a sequence name and its sequence and
#outputs a correctly formatted single fasta entry.

    my $seq_name = shift;
    my $seq = shift;
    my $fh = shift;

    if(length($seq) >0){
	$logger->debug("Writing sequence header:$seq_name length:".length($seq)) if($logger->is_debug());
	print $fh ">$seq_name\n";
	for(my $i=0; $i < length($seq); $i+=60){
	    my $seq_fragment = substr($seq, $i, 60);
	    print $fh $seq_fragment,"\n";
	}
    }
    else{
	$logger->warn("Skipping sequence $seq because length is ".length($seq));
    }

    return 1;
}

sub check_parameters{
    my ($options) = @_;
    
    if(!($options{'format'})){
	$options{'format'} = "multi";
    }
    $logger->debug("Fasta file format is $options{'format'}");

    if($options{'format'} ne "single" && $options{'format'} ne "multi"){
	$logger->fatal("Problem with --format $options{'format'}");
	pod2usage({-exitval => 2,  -message => "$0: --format must be 'single' or 'multi'", -verbose => 1, -output => \*STDERR});    
    }
}

