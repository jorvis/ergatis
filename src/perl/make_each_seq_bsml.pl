#!/usr/local/bin/perl



=head1  NAME 

make_each_seq_bsml.pl  - Generates DNA sequences in FASTA format

=head1 SYNOPSIS

USAGE:  make_each_seq_bsml.pl -b bsml_dir -o output_dir -a bsp_3839_assembly

=head1 OPTIONS

=over 4

=item

B<--bsml_dir,-b>   [REQUIRED] Dir containing BSML documents (repository)

=item *

B<--output_dir,-o> [REQUIRED] dir to save output to

=item *

B<--asmbl_ids,-a> Build sequences from this asmbl_id.  Multiple values can be comma separated

=item *

B<--asmbl_file>  name of the file containing a list of asmbl_ids

=item * 

B<--project,-p> save all sequences in this subdirectory (default in each asmbl_id's dir)

=item *

B<--log,-l> Log file

=item *

B<--debug,--D>  Debug level

=item *

B<--help,-h> This help message

=back

=head1   DESCRIPTION

make_each_seq_bsml.pl is designed to generate multi-fasta files of DNA sequences
for each gene given a list of assembly IDs.  The data is obtained through existing
BSML files in some central repository.  The output files are 1 gene per fasta file.
If -p flag is not used, the output files will be stored in the asmbl_id sub_dir under
the output_dir.  If -p is specified, all the fasta files will be stored in the directory
specified by -p which sits in the directory specified by -o.  You can pass a single
asmbl_id  to -a flag or a list of asmbl_ids separated by comma.  In addition, you can use
the --asmbl_file flag to specify a external file containing a list of asmbl_ids (1 per line)
to genrate sequences from.

IMPORTANT:

You must specify EITHER --asmbl_ids OR --asmbl_file flag. Not BOTH.

NOTE:  

Calling the script name with NO flags/options or --help will display the syntax requirement.


=cut


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Log::Log4perl qw(get_logger :levels :easy);
use BSML::BsmlReader;
use BSML::BsmlParserTwig;
use File::Basename;
use Pod::Usage;


my %options = ();
my $results = GetOptions (\%options, 'bsml_dir|b=s', 'output_dir|o=s', 'project|p=s', 'log|l=s', 'debug|D=s',
                                     'asmbl_ids|a=s', 'asmbl_file=s', 'help|h' ) || pod2usage();

###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $ASMBL_IDS  = $options{'asmbl_ids'};
my $output_dir = $options{'output_dir'};
$output_dir =~ s/\/+$//;       #remove terminating '/'s
my $BSML_dir   = $options{'bsml_dir'};
$BSML_dir =~ s/\/+$//;       #remove terminating '/'
my $debug      = $options{'debug'};
my $log        = $options{'log'};
my $project    = $options{'project'};
my $asmbl_file = $options{'asmbl_file'};

if( exists($options{'help'})) {
    pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT});
}


if(!$BSML_dir or !$output_dir) {
    pod2usage({-exitval => 1,  -message => "$0: All the required options are not specified", -verbose => 1, -output => \*STDERR});

}
if(!$asmbl_file and ! $ASMBL_IDS) {
    pod2usage(-verbose => 1, -message => "$0: Must specify either --asmbl_ids OR --asmbl_file.", -output => \*STDERR)
}

if($asmbl_file and $ASMBL_IDS) {
    pod2usage(-verbose => 1, -message => "$0: Must specify either --asmbl_ids OR --asmbl_file. NOT BOTH.", -output => \*STDERR) 
}

Log::Log4perl-> Log::Log4perl::init_and_watch($ENV{LOG4PERL_CONF}) if($ENV{LOG4PERL_CONF});
my $logger = get_logger('papyrus::pe');
$logger->level($INFO);
$logger->more_logging($debug) if($debug);

# Define a file appender or a screen appender
if($log){
    my $file_appender = Log::Log4perl::Appender->new(
						     "Log::Dispatch::File",
						     mode => "append",
						     filename  => $log);
    
    my $layout = 
	Log::Log4perl::Layout::PatternLayout->new(
						  "%d %p> %F{1}:%L %M - %m%n");
    $file_appender->layout($layout);
    $logger->add_appender($file_appender);
}else{
    my $screen_appender = Log::Log4perl::Appender->new(
						       "Log::Dispatch::Screen");	
    
    $logger->add_appender($screen_appender);
}
###-------------------------------------------------------###
my $min_dir = dirname($output_dir);
if(! -d $min_dir) {
    mkdir $min_dir;
    chmod 0777, $min_dir;
}
#create the directory that will hold all the individual peptide files if it doesn't exist
if(! -d $output_dir ) {
    mkdir $output_dir;
    $logger->debug("$output_dir not exists, creating $output_dir");
}
chmod 0777, $output_dir;

my $result;


my @asm_ids;
if($asmbl_file) {   #asmbl_id will be read from a flat file
    $logger->debug("Getting asmbl_ids from $asmbl_file");
    @asm_ids = read_asmbl_file($asmbl_file);
    if(!@asm_ids) {
	$logger->logdie("No asmbl_ids found in $asmbl_file.  Aborting");
    }
}else {
    $logger->debug("Getting asmbl_ids from -a flag");
    if($ASMBL_IDS =~ /all/i) {
	my @files = <$BSML_dir/*.bsml>;
	foreach (@files) {
	    my $basename = basename($_);
	    if($basename =~ /(.+)\.bsml/) {
		push(@asm_ids, $1);
	    }
	}
    } else {
	@asm_ids = split(/,/, $ASMBL_IDS);
    }
}
my $parser = new BSML::BsmlParserTwig;

if($project) {
    $logger->debug("Outputting all files into $output_dir/$project");
    consolidated_output(\@asm_ids);
}
else {
    $logger->debug("Outputting all files into asmbl_id subdir under $output_dir");
    regular_output(\@asm_ids);
} 

sub consolidated_output {

    my $assembly_ids = shift;
    my $final_output_dir = "$output_dir/".$project;
    if(! -d $final_output_dir ) {
	mkdir $final_output_dir;
    } else {
	unlink glob("$final_output_dir/*");
    }
    chmod 0777, $final_output_dir;
    foreach my $asmbl_id (@$assembly_ids) {
	my $bsml_file = "$BSML_dir/${asmbl_id}.bsml";
	if (-s $bsml_file) {
	    my $reader = BSML::BsmlReader->new();
	    $parser->parse( \$reader, $bsml_file );
	    my $extend_seq = $reader->get_all_cds_dna($asmbl_id) ;
	    my $output_seq_counter=0;
	    while( my ($gene, $seq) = each %$extend_seq ) {
		next if(length($seq) < 1); 
		$output_seq_counter++;
		$gene =~ s/_\d$//;
		my $protein_seq_id = $reader->cdsIdtoProteinSeqId($gene);
		my $gene_file = "$final_output_dir/${protein_seq_id}.seq";
		open(FILE, ">$gene_file") || die "Unable to write to $gene_file due to $!";
		my $fastaout = &fasta_out($protein_seq_id, $seq);
		print FILE $fastaout;
		close FILE;
		chmod 0777, $gene_file;
	    }
	    $logger->debug("$output_seq_counter sequences for $asmbl_id written to $final_output_dir ");
	} else {
	    $logger->logwarn("$bsml_file NOT found!  Skipping...");
	}

    }

}

sub regular_output {

    my $assembly_ids = shift;

    foreach my $asmbl_id (@$assembly_ids) {
	my $final_output_dir = "$output_dir/".$asmbl_id;
	if(! -d $final_output_dir ) {
	    mkdir $final_output_dir;
	} else {
	    unlink glob("$final_output_dir/*");
	}
	chmod 0777, $final_output_dir;
	my $bsml_file = "$BSML_dir/${asmbl_id}.bsml";
	if (-s $bsml_file) {
	    my $reader = BSML::BsmlReader->new();
	    $parser->parse( \$reader, $bsml_file );
	    my $extend_seq = $reader->get_all_cds_dna($asmbl_id) ;
	    my $output_seq_counter=0;
	    while( my ($gene, $seq) = each %$extend_seq ) {
		next if(length($seq) < 1); 
		$output_seq_counter++;
		$gene =~ s/_\d$//;
		my $protein_seq_id = $reader->cdsIdtoProteinSeqId($gene);
		my $gene_file = "$final_output_dir/${protein_seq_id}.seq";
		open(FILE, ">$gene_file") || die "Unable to write to $gene_file due to $!";
		my $fastaout = &fasta_out($protein_seq_id, $seq);
		print FILE $fastaout;
		close FILE;
		chmod 0777, $gene_file;
	    }
	    $logger->debug("$output_seq_counter sequences for $asmbl_id written to $final_output_dir ");
	} else {
	    $logger->logwarn("$bsml_file NOT found!  Skipping...");
	}
    }

}


sub fasta_out {
#This subroutine takes a sequence name and its sequence and
#outputs a correctly formatted single fasta entry (including newlines).

    my $seq_name = shift;
    my $seq = shift;

    my $fasta=">"."$seq_name"."\n";
    for(my $i=0; $i < length($seq); $i+=60){
	my $seq_fragment = substr($seq, $i, 60);
	$fasta .= "$seq_fragment"."\n";
    }
    return $fasta;

}

sub read_asmbl_file {

    my $file = shift;

    my @asmbl_id_list;

    open (IN, "$file")  or $logger->logdie("Unable to open and read $file!");
    my $line;
    while($line = <IN>) {
	chomp($line);
	next if($line =~ /^\s*$/);
	push(@asmbl_id_list, $line);
    }
    close IN;

    return @asmbl_id_list;

}







