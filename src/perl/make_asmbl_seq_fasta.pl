#!/usr/local/bin/perl

=head1  NAME 

make_asmbl_seq_fasta.pl  - Generates assembly sequences in FASTA format from BSML files

=head1 SYNOPSIS

USAGE:  make_asmbl_seq_fasta.pl -b bsml_dir -o outputfile -a bsp_3839_assembly

=head1 OPTIONS

=over 4

=item *

B<--bsml_dir,-b>   [REQUIRED] Dir containing BSML documents (repository)

=item *

B<--output,-o>   [REQUIRED] name of the output file

=item *

B<--asmbl_ids,-a> Build sequences from this asmbl_id.  Multiple values can be comma separated

=item *

B<--asmbl_file,-f>  name of the file containing a list of asmbl_ids (1 per line)

=item * 

B<--exclude_asmbl_ids,-e>  Exclude the asmbl_id(s).   Multiple values can be comma separated

=item *

B<--exclude_file>  name of the file containing a list of asmbl_ids (1 per line) to exclude

=item *

B<--log,-l> Log file

=item *

B<--debug,--D>  Debug level

=item *

B<--help,-h> This help message

=back

=head1   DESCRIPTION

make_asmbl_seq_fasta.pl is designed to generate multi-fasta files of assembly sequences
given a list of assembly IDs.  The data is obtained through existing
BSML files in some central repository.  User can specify the sequences to fetch by
--asmbl_ids flag AND/OR --asmbl_file flag.  Using the first options, asmbl_ids can be
comma separated while the latter option requires the asmbl_ids to be stored 1 per line
in a file.  If both options are used, a unique list of asmbl_ids from both sources will
be fetched.  The user can also specify specific asmbl_id(s) to exclude by a similar 
mechansim (i.e. --exclude_file or --exclude_asmbl_ids flag).  As a short cut, the user
can supply an argument of 'all' to --asmbl_ids or --exclude_asmbl_ids.  This will use all
the asmbl_ids found in the repository.  So if user invoke --asmbl_ids all AND 
--exclude_asmbl_ids all at the same time, the resulting asmbl_ids will be empty as
the two calls completely cancels each other out.

Samples:

1. make a fasta file of assembly sequence A and B called fasta.txt
   make_gene_pep_fasta.pl -b bsml_dir -o fasta.txt -a A,B 

2. make a fasta file of all assemblies except D,E, and F
   make_gene_pep_fasta.pl -b bsml_dir -o fasta.txt -a all -e D,E,F

IMPORTANT:

You must specify at least --asmbl_ids OR --asmbl_file flag.

NOTE:  

Calling the script name with NO flags/options or --help will display the syntax requirement.


=cut


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BSML::BsmlReader;
use BSML::BsmlParserTwig;
use File::Basename;
use File::Path;
use Pod::Usage;

umask(0000);
my %options = ();
my $results = GetOptions (\%options, 'bsml_dir|b=s', 'asmbl_ids|a=s', 'output|o=s', 'asmbl_file|f=s', 
                                     'exclude_file=s', 'exclude_asmbl_ids|e=s', 'help|h', 'man' ) || pod2usage();

###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $asmbl_ids        = $options{'asmbl_ids'};
my $output_file      = $options{'output'};
my $output_dir       = dirname($output_file);
my $BSML_dir         = $options{'bsml_dir'};
$BSML_dir =~ s/\/+$//;         #remove terminating '/'s
my $asmbl_file       = $options{'asmbl_file'};
my $exclude_asmbl_file     = $options{'exclude_file'};   #file containing a list of asmbl_ids to exclude
my $exclude_asmbl_ids = $options{'exclude_asmbl_ids'};   

my $parser = new BSML::BsmlParserTwig;

&cmd_check();
###-------------------------------1------------------------###

my @asm_ids = determine_all_asmbl_id($asmbl_ids, $asmbl_file, $BSML_dir);
my @exclude_asm_ids = determine_all_asmbl_id($exclude_asmbl_ids, $exclude_asmbl_file, $BSML_dir);
my @final_asm_ids = make_final_asmbl_ids(\@asm_ids, \@exclude_asm_ids);

make_assembly_fasta(\@final_asm_ids);


sub make_final_asmbl_ids {
#returns a list of asmbl_ids in @$asm_ids that 
#is not in @$no_asm_ids

    my $asm_ids    = shift;
    my $no_asm_ids = shift;

    my (@final_asmbl_ids, %seen);

    if(@$no_asm_ids == 0) {
	@final_asmbl_ids = @$asm_ids;
    } else {
	@seen{@$no_asm_ids} = ();
	foreach my $item(@$asm_ids) {
	    push(@final_asmbl_ids, $item) unless exists $seen{$item};
	}
    }

    return @final_asmbl_ids;

}

sub make_assembly_fasta {

    my $assembly_ids = shift;

    open(FILE, ">$output_file") || die "Can't open $output_file due to $!";
    foreach my $asmbl_id (@$assembly_ids) {
	my $bsml_file = "$BSML_dir/${asmbl_id}.bsml";
	if (-s $bsml_file) {
	    my $reader = BSML::BsmlReader->new();
	    $parser->parse( \$reader, $bsml_file );
	    my $seq_list = $reader->returnAllSequences();
	    foreach my $seq_obj (@$seq_list) {
		if($seq_obj->returnattr('molecule') eq 'dna') {
		    my $seq_id = $seq_obj->returnattr('id');
		    my $raw_seq = $reader->subSequence($seq_id, -1, -1, '0');
		    my $fastaout = &fasta_out($seq_id, $raw_seq);
		    print FILE $fastaout;
                }
            }
	} else {
	    print STDERR "$bsml_file NOT found!!!! Skipping...\n";
        }
    }
    close FILE;
    unlink $output_file if(-z $output_file);
    chmod 0666, $output_file;
}

sub determine_all_asmbl_id {
#given $asmbl_ids and/or $asmbl_file
#return a list of asmbl_ids
#if both arguments are supplied, will return a list of
#asmbl_ids (duplicates are removed)
    
    my $asmbl_ids  = shift;
    my $asmbl_file = shift;
    my $bsml_dir   = shift;

    my @final_asmbl_ids;
    my (%unique_asmbl_ids, @id_list);

    
    #check asmbl_id passed via command line option
    if(defined($asmbl_ids)) {
	if($asmbl_ids =~ /^all$/) {
	    @id_list = get_all_asmbl_id_via_all($bsml_dir);
	} else{
	    @id_list = split(/,/, $asmbl_ids);
	}
	foreach my $asmbl_id (@id_list) {
	    push(@final_asmbl_ids, $asmbl_id) unless $unique_asmbl_ids{$asmbl_id}++;
	}
    }
       
    #check asmbl_id from a flat file
    if(defined($asmbl_file)) {
	@id_list = read_asmbl_file($asmbl_file);
	foreach my $asmbl_id (@id_list) {
	    push(@final_asmbl_ids, $asmbl_id) unless $unique_asmbl_ids{$asmbl_id}++;
	}
    }

    return @final_asmbl_ids;

}


sub read_asmbl_file {
#parses the flat file containing
#a list of asmbl_ids

    my $file = shift;

    my @asmbl_id_list;

    open (IN, "$file") or die("Unable to read $file due to $!");
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


sub get_all_asmbl_id_via_all {
#shortcut way of getting all asmbl_ids by user
#specifying --asmbl_id all on command line
#The asmbl_ids are obtained by removing the suffix of 
#the bsml file

    my $bsml_dir = shift;

    my @asm_ids;

    opendir(DIR, $bsml_dir) or die "Unable to access $bsml_dir due to $!";
    while( my $file = readdir(DIR)) {
	next if ($file =~ /^\.{1,2}$/);  #skip  "." ,  ".."
	if($file =~ /(.+)\.bsml$/) {
	    push(@asm_ids, $1);
	}
    }

    return  @asm_ids;

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


sub cmd_check {
#quality check

    if( exists($options{'man'})) {
	pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT});
    }   

    if( exists($options{'help'})) {
	pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT});
    }

    if(!$output_file or !$BSML_dir) {
	pod2usage({-exitval => 2,  -message => "$0: All the required options are not specified", -verbose => 1, -output => \*STDERR});
    }
    
    if(!$asmbl_file and ! $asmbl_ids) {
	pod2usage(-exitval =>2, -verbose => 1, -message => "$0: Must specify either --asmbl_ids OR --asmbl_file.", -output => \*STDERR);
    }
    

    #checking BSML repository directory
    if(! -d $BSML_dir) {
	print STDERR "BSML repository directory \"$BSML_dir\" NOT found.  Aborting...\n";
	exit 5;
    }
    #check for user defined output directory
    if(! -d $output_dir) {
	mkpath($output_dir) or die "Unable to create $output_dir.  Aborting...\n";
	#chmod 0777, $output_dir;
    }

}
