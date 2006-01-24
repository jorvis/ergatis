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

=head1   DESCRIPTION

=cut


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Config::IniFiles;

BEGIN {
use Workflow::Logger;
use BSML::BsmlBuilder;
}

my %options = ();
my $results = GetOptions (\%options, 
			  'msffile|f=s',
			  'output|o=s',
			  'log|l=s',
			  'debug=s',
			  'conf=s',
			  'analysis_conf=s',
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

$logger->info("Instantiating the BSML builder object");
my $builder = new BSML::BsmlBuilder;
$logger->logdie("builder was not defined") if (!defined($builder));
my $analysis_name = &get_analysis_name($options{'analysis_conf'});
my $MSF_alignments = process_MSF_file("$options{'msffile'}");
if($MSF_alignments->{'mol_type'} eq 'polypeptide'){
    $MSF_alignments->{'mol_type'} = 'protein';
}
else{
    $MSF_alignments->{'mol_type'} = 'nucleotide';
}
if(keys %$MSF_alignments > 1){   #skip empty msf files
    my $table = $builder->createAndAddMultipleAlignmentTable('molecule-type' => $MSF_alignments->{'mol_type'},
							     );
    $table->addattr('class', 'match');
    $logger->logdie("table was not defined") if (!defined($table));
    
    my $summary = $builder->createAndAddAlignmentSummary( 
							  'multipleAlignmentTable' => $table,
							  'seq-type'               => $MSF_alignments->{'mol_type'},
							  'seq-format'             => 'msf'
							  );
    $logger->logdie("summary was not defined") if (!defined($summary));
    
    
    my $aln = $builder->createAndAddSequenceAlignment( 'multipleAlignmentTable' => $table );
    $logger->logdie("aln was not defined") if (!defined($aln));
    $table->addBsmlLink('analysis', '#'."$analysis_name", 'computed_by');
    my $seqnum=0;
    my $sequences_tag;
    
    foreach my $seq (keys %{ $MSF_alignments->{'polypeptide'} }) {
	$logger->logdie("seq was not defined") if (!defined($seq));
	
	$seqnum++;
	
	my $alignment = join ('', @{ $MSF_alignments->{'polypeptide'}->{$seq}->{'alignment'} });
	$logger->logdie("alignment was not defined") if (!defined($alignment));
	
	
	my $align_length = $MSF_alignments->{'polypeptide'}->{$seq}->{'length'} if ((exists $MSF_alignments->{'polypeptide'}->{$seq}->{'length'}) and (defined($MSF_alignments->{'polypeptide'}->{$seq}->{'length'})));
	
	$logger->logdie("align_length was not defined") if (!defined($align_length));
	
	#IMPORTANT!!!!
	#In order to ensure that each seq in a multiple sequence alignment is truly
	#unique, the seq-name and name will be in the form "polypeptide_accession:seqnum"
	#i.e. (ana1.10005.m00234_polypeptide:1). 
	

	#
	# bugzilla case 1979
	# sundaram@tigr.org
	# 2005.07.12
	# Polypeptide identifiers in clustalw output are being truncated.
	#
	if (length($seq) == 30 && $seq !~ /_polypeptide$/){
	    $logger->warn("polypeptide identifier '$seq' was truncated by clustalw, repairing now");
	    $seq =~ s/_[^_]+$/_polypeptide/;
	}
	if ($alignment !~ /_polypeptide\s/){
	    $logger->warn("polypeptide identifier in the alignment of '$seq' was truncated by clustalw, repairing now");
	    $alignment =~ s/_[^_\s]+\s/_polypeptide /g;
	}

	$builder->createAndAddAlignedSequence(
					      'alignmentSummary' => $summary,
					      'seqnum'           => $seqnum,
					      'length'           => $align_length,
					      'name'             => "$seq:$seqnum"
					      );
	
	$builder->createAndAddSequenceData(
					   'sequenceAlignment' => $aln,
					   'seq-name'          => "$seq:$seqnum",
					   'seq-data'          => $alignment
					   ); 
	$sequences_tag .= "$seqnum:";
    }
    $aln->addattr( 'sequences', $sequences_tag );
}

## add the analysis element
$builder->createAndAddAnalysis(
			   id => $analysis_name,
			   sourcename => $options{'output'},
			   );

$builder->write( $options{'output'} );

sub check_parameters{
    my ($options) = @_;
    
    if(0){
	pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});    
    }
}

#-------------------------------------------------------------
# process_MSF_file()
#
#-------------------------------------------------------------
sub process_MSF_file {

    $logger->debug("Entered process_MSF_file") if $logger->is_debug;

    my $file = shift;
    $logger->logdie("file was not defined") if (!defined($file));

    my $MSF_alignments ={};
    open(MSF, "$file") or die "Unable to open $file due to $!";
    my $line;
    my $msf_type;
    while(defined($line = <MSF>) and $line !~ /^\/\//) {
	if( $line =~ /MSF:\s*([\S]+)\s*Type:\s*([\S]+)\s*Check/) {
	    my $msf_length = $1;
	    return undef if($msf_length == 0);   #abort if align_len = 0

	    if($2 eq 'P') {
		$msf_type = 'polypeptide';
	    }elsif($2 eq 'N') {
		$msf_type = 'nucleotide';
	    }else {
		$msf_type = 'polypeptide';
	    }
	    $MSF_alignments->{'mol_type'} = $msf_type;
	}
	#if($line =~ /Name:\s+([\S]+)\s+[o]{2}\s+Len:\s+([\S]+)\s+Check:\s+([\S]+)\s+Weight:\s+([\S]+)/) {

	if($line =~ /Name:\s*([\S]+)\s*[o]{2}\s*Len:\s*([\S]+)\s*Check:\s*([\S]+)\s*Weight:\s*([\S]+)/) {
	    my $name    = $1;
	    my $ali_len = $2;
	    my $check   = $3;
	    my $weight  = $4;
	    
	    $MSF_alignments->{'polypeptide'}->{$name}->{'length'} = $ali_len;
	    $MSF_alignments->{'polypeptide'}->{$name}->{'check'}  = $check;
	    $MSF_alignments->{'polypeptide'}->{$name}->{'weight'} = $weight;
	    $MSF_alignments->{'polypeptide'}->{$name}->{'alignment'} = [];
	}
    }

    my $replacements;
    my $spaces;
    while($line = <MSF>) {
	if($line =~ /^([\S]+)/) {
	    my $name = $1;
	    if(exists($MSF_alignments->{'polypeptide'}->{$name})) {
		push( @{ $MSF_alignments->{'polypeptide'}->{$name}->{'alignment'} }, $line );
            } else {
		print STDERR "ERROR, $name is not valid polypeptide name for $file\n";
		exit;
            }
	}
    }

    return $MSF_alignments;

}#end sub process_MSF_file()


sub get_analysis_name{
    my($conf) = @_;
    my $analysis_name = "clustalw_analysis";
    if(-e $conf){
	my $cfg = new Config::IniFiles( -file => $conf);
	my @sections = $cfg->Sections();
	my $name;
	for my $section (@sections) {
	    $name = $cfg->val($section,'$;NAME$;');
	    if($name ne ""){
		last;
	    }
	}
	if($name ne ""){
	    $analysis_name = $name."_analysis";
	}
    }
    return $analysis_name;
}
