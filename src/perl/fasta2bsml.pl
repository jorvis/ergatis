#!/local/perl/bin/perl

=head1 NAME

fasta2bsml.pl

=head1 SYNOPSIS

USAGE:  

=head1 OPTIONS

=over 8

=item 
    
=back

=head1 DESCRIPTION

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use BSML::BsmlBuilder;
use Workflow::Logger;

my %options = ();
my $results = GetOptions (\%options,
			  'fasta_file|f=s',
			  'fasta_list|l=s',
			  'fasta_dir|d=s',
			  'output|o=s', 
			  'genus|g=s',
			  'species|s=s',
			  'strain|S=s',
			  'source_database|u=s', 
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


my $doc = new BSML::BsmlBuilder(); 

my $organismelt = $doc->createAndAddOrganism( 
					   'genome'  => $doc->createAndAddGenome(),
					   'genus'   => $options{'genus'},
					   'species' => $options{'species'}
					);

my $strainelt = $doc->createAndAddStrain( 
					  'organism'        => $organismelt,
					  'name'            => $options{'strain'}, 
					  'database'        => $options{'source_database'},
					  'source_database' => $options{'source_database'}
				    );				   

my @files;
if ($options{'fasta_file'} ne ""){
    my @inputfiles = split(/,/,$options{'fasta_file'});
    foreach my $file (@inputfiles){
	$file =~ s/\s//g;
	if($file ne ""){
	    if((-s $file)) {
		$logger->debug("Adding file $file for processing") if($logger->is_debug());
		push @files,$file;
	    }
	    else{
		$logger->warn("Error reading file $file");
	    }
	}
    }
}
if ($options{'fasta_list'} ne ""){
    my @filelists = split(/,/,$options{'fasta_list'});
    foreach my $filelist (@filelists){
	$filelist =~ s/\s//g;
	if($filelist ne ""){
	    open FILE, "$filelist" or $logger->logdie("Can't open file $filelist");
	    while (my $filename=<FILE>){
		chomp $filename;
		if(-s $filename) {
		    $logger->debug("Adding file $filename for processing") if($logger->is_debug());
		    push @files,$filename;
		}
		else{
		    $logger->warn("Error reading file $filename");
		}
	    }
	}
    }
}

if ($options{'fasta_dir'} ne "") {
    my @fastadirs = split(/,/,$options{'fasta_dir'});
    foreach my $dir (@fastadirs){
	$dir =~ s/\s//g;
	if($dir ne ""){
	    if(-r $dir){
		opendir(DIR, $dir) or $logger->warn("Unable to access $dir due to $!");
		while( my $filename = readdir(DIR)) {
		    if($filename =~ /(.+)\.fasta$/ || $filename =~ /(.+)\.pep$/ || $filename =~ /(.+)\.fsa$/) {
			if(-s "$dir/$filename"){
			    $logger->debug("Adding file $dir/$filename for processing") if($logger->is_debug());
			    push (@files, "$dir/$filename");
			}
			else{
			    $logger->warn("Error reading file $dir/$filename");
			}
		    }
		}
	    }
	    else{
		$logger->warn("Error reading directory $dir");
	    }
	}
    }
}

if(scalar(@files)==0){
    $logger->logdie("No files found");
}




foreach my $fasta_file (@files){

    my ($uid, $seq);
    
    
    open(F_IN, $fasta_file) || die ("no can find: $fasta_file");
    my $first_time = 1;
    while(<F_IN>) {
	s/\n//;
	if(/^>/) {
	    s/>//;
	    s/\s.*//;
	    $uid = $_;
	    
	    $uid = "_" . $uid if ($uid =~ /^[0-9]/);
	    
	    if ($first_time != 1) {
		add_stuff($doc, $uid, $seq, $fasta_file);
	    }
	    
	    $first_time = 0;
	    my $seq;
	}
	else {
	    $seq .= $_;
	}
    }
}

$doc->write("$options{'output'}");

if(! -e "$options{'output'}"){
    die ("File not created $options{'output'}");
}

sub add_stuff {
    my($doc, $uid, $seq, $fasta_file) = @_;
    my($seq_length);

    $seq_length = length($seq);

    my $asmseq = $doc->createAndAddExtendedSequenceN( 'id' => $uid, 
						      'title' => '', 
						      'length' => $seq_length, 
						      'molecule' => 'aa', 
						      'locus' => '', 
						      'dbsource' => '', 
						      'icAcckey' => '', 
						      'strand' => '');

    $doc->createAndAddSeqDataImport($asmseq, 
				    'fasta', 
				    $fasta_file, 
				    $uid);
}

sub check_parameters{
    my ($options) = @_;
    
}

