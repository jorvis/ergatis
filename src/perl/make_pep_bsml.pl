#!/usr/local/bin/perl

=head1  NAME 

make_pep_bsml.pl  - Generates protein sequences in FASTA format

=head1 SYNOPSIS

USAGE:  make_pep_bsml.pl -b bsml_dir -o output_dir -a bsp_3839_assembly -g

=head1 OPTIONS

=over 4

=item *

B<--bsml_dir,-b>   [REQUIRED] Dir containing BSML documents (repository)

=item *

B<--output_dir,-o> [REQUIRED] dir to save output to

=item *

B<--asmbl_ids,-a> Build sequences from this asmbl_id.  Multiple values can be comma separated

=item *

B<--asmbl_file>  name of the file containing a list of asmbl_ids

=item * 

B<--project,-p> name of the TOTAL peptide fasta file. Required if -g -e options ARE NOT invoked

=item *

B<--each_genome,-g> save fasta file one asmbl_id per fasta file

=item *

B<--each_file,-e> save fasta file one sequence per fasta file

=item *

B<--log,-l> Log file

=item *

B<--debug,--D>  Debug level

=item *

B<--help,-h> This help message

=back

=head1   DESCRIPTION

make_pep_bsml.pl is designed to generate multi-fasta files of protein sequences
for each gene given a list of assembly IDs.  The data is obtained through existing
BSML files in some central repository.  There are three modes:  -g outputs all the
sequences for 1 assembly into a fasta file and stores that file in the arguments
of --output_dir.  -e outputs each sequence into a fasta 
file and all the fasta files into the a dir called with the asmbl_id name which sits
in the --output_dir argument. Without -e or -g, a -p flag must be used to specify the
name of the pep file.  This file will be stored in the argument of --output_dir.  The
last mode is ideal for generating sequences from multiple assembly ids and outputting
them into 1 file (e.g. PNEUMOl.pep)

IMPORTANT:

You must specify EITHER --asmbl_ids OR --asmbl_file flag. Not BOTH.
You must specify --project if -g and -e ARE NOT invoked

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
my $results = GetOptions (\%options, 'bsml_dir|b=s', 'asmbl_ids|a=s', 'simple_header|s', 'project|p=s', 'verbose|v', 'log|l=s', 'debug|D=s',
                                     'output_dir|o=s', 'help|h', 'each_file|e', 'each_genome|g', 'asmbl_file=s' ) || pod2usage();

###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $ASMBL_IDS       = $options{'asmbl_ids'};
my $simple_header   = $options{'simple_header'};
my $each_file       = $options{'each_file'} || 0;
my $each_genome     = $options{'each_genome'} || 0;
my $project         = $options{'project'};
my $verbose         = $options{'verbose'};
my $output_dir      = $options{'output_dir'};
$output_dir =~ s/\/+$//;       #remove terminating '/'s
my $BSML_dir        = $options{'bsml_dir'};
$BSML_dir =~ s/\/+$//;         #remove terminating '/'s
my $debug           = $options{'debug'};
my $log             = $options{'log'};

my $asmbl_file      = $options{'asmbl_file'};

if( exists($options{'help'})) {
    pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT});
}

if(!$output_dir or !$BSML_dir) {
    pod2usage({-exitval => 1,  -message => "$0: All the required options are not specified", -verbose => 1, -output => \*STDERR});
}
if(!$asmbl_file and ! $ASMBL_IDS) {
    pod2usage(-verbose => 1, -message => "$0: Must specify either --asmbl_ids OR --asmbl_file.", -output => \*STDERR);
}

if($asmbl_file and $ASMBL_IDS) {
    pod2usage(-verbose => 1, -message => "$0: Must specify either --asmbl_ids OR --asmbl_file. NOT BOTH.", -output => \*STDERR);
    
}

if($each_file and $each_genome) {
    pod2usage(-verbose => 1, -message => "$0: --each_genome(-g) and --each_file(-e) CANNOT be invoked at the same time", -output => \*STDERR) ;
}

if(!$each_file and !$each_genome and !$project) {
    pod2usage(-verbose => 1, -message => "$0: Must specify project name without --each_genome(-g) and --each_file(-e)",  -output => \*STDERR) ;
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
    
    my $layout = Log::Log4perl::Layout::PatternLayout->new("%d %p> %F{1}:%L %M - %m%n");
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

if(! -d $output_dir ) {
    mkdir $output_dir;
}
chmod 0777, $output_dir;

my @asm_ids;
if($asmbl_file) {   #asmbl_id will be read from a flat file
    $logger->debug("Getting asmbl_ids from $asmbl_file");
    @asm_ids = read_asmbl_file($asmbl_file);
    if(!@asm_ids) {
	$logger->logdie("No asmbl_ids found in $asmbl_file.  Aborting");
    }
} else {   #asmbl_id will be read from --asmbl_ids flag
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
$logger->debug("asmbl_id(s) to do: @asm_ids");

my $parser = new BSML::BsmlParserTwig;

if($each_genome) {
    make_fasta_for_each_genome(\@asm_ids);
    print "Finished making fasta pep files for EACH genome\n" if($verbose);
}elsif($each_file) {
    make_fasta_for_each_gene(\@asm_ids);
    print "Finished making fasta pep file for each gene\n" if($verbose);
}else {
    make_pep_for_ALL_genomes(\@asm_ids);
    print "Finished making fasta pep file for all assemblies\n" if($verbose);
}



sub make_fasta_for_each_gene  {

    my $assembly_ids = shift;

    $logger->debug("generating fasta files: 1 sequence per file");

    foreach my $asmbl_id (@$assembly_ids) {
	next if($asmbl_id !~ /\d+/);
	my $final_output_dir = "$output_dir/".$asmbl_id;
	#create the directory that will hold all the individual peptide files if it doesn't exist
	if(! -d $final_output_dir ) {
	    mkdir $final_output_dir;
	}else {         #delete old crap if that directory already exists
	    unlink glob("$final_output_dir/*");
	}
	chmod 0777, $final_output_dir;
	my $bsml_file = "$BSML_dir/${asmbl_id}.bsml";
	if (-s $bsml_file) {
	    my $reader = BSML::BsmlReader->new();
	    $parser->parse( \$reader, $bsml_file );
	    my $proteins = $reader->get_all_protein_aa($asmbl_id);
	    while (my ($identifier, $seq) = each %$proteins) {
		$identifier =~ s/_\d$//;
		my $pep_file = "$final_output_dir/$identifier.fsa";
		open(FILE, ">$pep_file") || die "Unable to write to $pep_file due to $!";
		my $prot_name = $identifier;
		my $fastaout = &fasta_out($prot_name, $seq);
		print FILE $fastaout;
		close FILE;
		chmod 0777, $pep_file;
	    }
	} else {
	    $logger->logdie("$bsml_file NOT found!  Exiting...");
	}

    }

}

sub make_fasta_for_each_genome {

    my $assembly_ids = shift;

    $logger->debug("generating fasta files: 1 assembly per file");
    
    foreach my $asmbl_id (@$assembly_ids) {
	my $pep_file = "$output_dir/$asmbl_id".".pep";
	open(FILE, ">$pep_file") || die "Unable to write to $pep_file due to $!";
	my $bsml_file = "$BSML_dir/${asmbl_id}.bsml";
	if (-s $bsml_file) {
	    my $reader = BSML::BsmlReader->new();
	    $parser->parse( \$reader, $bsml_file );
	    my $proteins = $reader->get_all_protein_aa($asmbl_id);
	    while (my ($identifier, $seq) = each %$proteins) {
		$identifier =~ s/_\d$//;	
		my $prot_name = $identifier;
		my $fastaout = &fasta_out($prot_name, $seq);
                print FILE $fastaout;
            }
	    close FILE;
	    chmod 0777, $pep_file;
	}else {
	    $logger->logdie("$bsml_file NOT found!  Exiting...");
	}
    }

}


sub make_pep_for_ALL_genomes {

    my $assembly_ids = shift;

    $logger->debug("generating fasta files: multiple assemblies per file");
    if($simple_header) {
	$logger->debug("using simple header");
    } else {
	$logger->debug("using complex header");
    }
	  

    my $pep_file = "$output_dir/$project.pep";
    open(FILE, ">$pep_file") || die "Cant open $pep_file due to $!";
    foreach my $asmbl_id (@$assembly_ids) {
	my $bsml_file = "$BSML_dir/${asmbl_id}.bsml";
	if (-s $bsml_file) {
	    my $reader = BSML::BsmlReader->new();
	    $parser->parse( \$reader, $bsml_file );
	    my $proteins = $reader->get_all_protein_aa($asmbl_id);
	    while (my ($identifier, $seq) = each %$proteins) {
		$identifier =~ s/_\d$//;
		my $prot_name = $simple_header ? "$identifier $asmbl_id" : "${asmbl_id}:$identifier";
		my $fastaout = &fasta_out($prot_name, $seq);
                print FILE $fastaout;
            }
	} else {
	    $logger->logdie("$bsml_file NOT found!  Exiting...");
        }
    }
    close FILE;
    chmod 0777, $pep_file;
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

    open (IN, "$file")  or $logger->logdie("Unable to read $file due to $!");
    my $line;
    while($line = <IN>) {
	chomp($line);
	next if($line =~ /^\s*$/);
	$line =~ s/^\s+//;
	$line =~ s/\s+$//;
	push(@asmbl_id_list, $line);
    }
    close IN;

    return @asmbl_id_list;

}



sub print_usage {


    print STDERR "SAMPLE USAGE:  make_pep_bsml.pl -b bsml_dir -o output_dir -a bsp_3839_assembly\n";
    print STDERR "  --bsml_dir    = dir containing BSML doc\n";
    print STDERR "  --output_dir  = dir to save output to\n";
    print STDERR "  --asmbl_ids  (multiple values can be comma separated)\n";
    print STDERR "               (-a all  grabs all asmbl_ids)\n";
    print STDERR "  --each_genome = save fasta file individually for each genome\n";
    print STDERR "  --each_file   = save fasta file individually for each gene\n";
    print STDERR "  --asmbl_file  = name of the file containing a list of asmbl_ids\n";
    print STDERR "  --project     = name of the total peptide fasta file\n";
    print STDERR " NOTE* --each_genome and --each_file cannot be invoked concurrently\n";
    print STDERR "       --Specify Either --asmbl_ids OR --asmbl_file \n";
    print STDERR "  --help = This help message.\n";
    exit 1;

}
