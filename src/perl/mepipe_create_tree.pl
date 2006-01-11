#!/usr/local/bin/perl

=head1  NAME

mepipe_create_tree.pl - molecular evolution pipeline script that creates tree files.

=head1 SYNOPSIS

USAGE: mepipe_create_tree.pl
        --input_path=/path/to/directory/containing/input/files/
        --temp_path=/path/to/store/temporary/output
        --output_path=/path/to/store/output
		--paml_bin_dir=/path/to/paml/bin
        --paup_bin_dir=/path/to/paup/bin
		--emboss_bin_dir=/path/to/emboss/bin
        --debug=4

=head1 OPTIONS

B<--input_prefix,-p>
	The file prefix for the gene/cog to be analyzed.
	(NB: This was referred to as 'gene name' in the original script.)

B<--input_path,-i>
    The path to the directory containing required input files.

B<--temp_path,-t>
	The temporary path for execution steps. (eg: $;TEMP_DIR$;)
	
B<--output_path,-o>
	The path to the directory to store output.
	
B<--paup_bin_dir,-b>
	Full path containing PAUP executable.
	
B<--debug,-d>
    Debug level.  Use a large number to turn on verbose debugging.

B<--log,-l>
    Log file

B<--help,-h>
    This help message

=head1   DESCRIPTION

*** TO BE COMPLETED ***

=head1 INPUT

*** TO BE COMPLETED ***

This script takes as input a file with a list of gene names (one per line)
and, for each gene, an amino acid clustal alignment (<geneName.aa>) 
and respective nucleotide seqs (<geneName.seq>) in fasta format,
converts the aa alignment into nucleotide alignment, and creats PAUP 
and PAML input files
Runs PAUP and PAML

NOTE: for each gene, sequence names in the amino acid alignment file
and corresponding nucleotide seq file have to be identical; 
names should contain no spaces.

=head1 OUTPUT

*** TO BE COMPLETED ***

=head1 CONTACT

	Author:
	Joana Silva
	jsilva@tigr.org

	Modifications for running with workflow:
    Brett Whitty
    bwhitty@tigr.org

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use File::Copy;
BEGIN {
    require '/usr/local/devel/ANNOTATION/cas/lib/site_perl/5.8.5/Workflow/Logger.pm';
    import Workflow::Logger;
}

my %options = ();
my $results = GetOptions (\%options,
						  'input_prefix|p=s',
                          'input_path|i=s',
						  'temp_path|t=s',
                          'output_path|o=s',
#						  'paml_bin_dir=s',
						  'paup_bin_dir|b=s',
#						  'emboss_bin_dir=s',
#						  'count_file=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || pod2usage();

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = Workflow::Logger::get_logger();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDERR} );
}

#&check_parameters(\%options);

#my $usage = "usage: $0 <fileList>";
# $0 = program name;  
# <fileList> = file with gene names, one per line; 

#my $infile = $ARGV[0] or die "$usage\n\n";

my $aaAlign;
my $nuclSeqs;
my $nuclAlign;
#my $pathtofiles = "";
my $geneName;

my $outgroup;
my $batchfile;
my @speciesNames;
my $noSpecs;
my $aaAlignLength;
my $PAMLrunmode;

my $err;

#open (INFILE, $infile) or die $!;

#while (<INFILE>) {

#   chomp ( $geneName= $_);

##
## Should be running checks on input/work/output paths that
## A) they exist
foreach my $dir((
#		'paml_bin_dir', 
#		'emboss_bin_dir', 
		'paup_bin_dir',
		'temp_path',
		'input_path',
		'output_path', 
	        )) {
		if (! -d $options{$dir}) {
			die "'$dir' specified as '".$options{$dir}."' does not exist.\n";
		}
		unless ($options{$dir} =~ /\/$/) {
			$options{$dir} .= "/";
		}
}

## B) they do or do not end in / depending on what we want
##

	chdir($options{'temp_path'});
	
	$geneName = $options{'input_prefix'};

#    print "\n$geneName\n";
    
    $aaAlign = "$geneName.aln";
    $nuclSeqs = "$geneName.seq";
    $nuclAlign = "$geneName.nuc.aln";

    $outgroup = "";
    $batchfile = "";
    @speciesNames = ("");
    $noSpecs = 0;

    &readSpeciesNames($options{'temp_path'}."$aaAlign.nex");   #seq names are flanked by single quotes
                                         #list of seq names is writen to array @speciesNames
                                         # checks length of aa alignments and stores in $aaAlignLength
    $noSpecs = scalar @speciesNames;


    if ($noSpecs > 3) {                            #if no. species > 3, PAUP is run to reconstruct tree
		$batchfile = &makePAUPbatch($options{'temp_path'}."$aaAlign");   # makes a batch file to run PAUP
		#chmod(0775, $batchfile);
		chmod(0777, $batchfile);
		#print($options{'paup_bin_dir'}."paup -n -f $batchfile\n");
		$err = system($options{'paup_bin_dir'}."paup -n -f \"$batchfile\"");
		if ($err) {
			die "Paup run failed.\n";
		}
		unlink($batchfile);
		&makeTreeFile($options{'temp_path'}."$geneName", $noSpecs);
		# makes a batch file to run PAUP
        # runs paup resulting in a tree stored in $geneName.aln.PAUPtree
    } else {
		my $treeout = $options{'temp_path'}."$geneName.PAMLtree";
		$" = ",";
		open (TREEOUT, ">$treeout");
		print TREEOUT "(@speciesNames)\n";
		close TREEOUT;
		#print "Processed data and tree files\n";
	}

#print "\nDONE\n";


# end main
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

sub makeTreeFile {
    my $treein = "$_[0].aln.PAUPtree";
    my $treeout = "$_[0].PAMLtree";
    my @line;
    my $seqlen= $_[2];
    my $notaxa= $_[1];
    my @trees;
    my $noTrees=0;
    my @arr;
    my $noElem;
    my $theTree;
    my $openingParenthesis;
    my $closingParenthesis;
    my $parenthesisCounter=0;
    my $lineLen;
			
	## THERE WILL PROBABLY BE SOME PATH PROBLEMS HERE
    open (TREEFILEin, "$treein");
    
    while (<TREEFILEin>) {
		if (/^tree/) {
	    	chomp;
	    	@line= split (' ', $_);
	    	my $noElements = @line;
	    	if ($line[$noElements-1] =~ /^\(/ ){
				$trees[$noTrees++] = "$line[$noElements-1]\n";
	    	}
		}
    }
    close TREEFILEin;
	## THERE WILL PROBABLY BE SOME PATH PROBLEMS HERE
    open (TREEFILEout, ">$treeout");

    print TREEFILEout " $notaxa 1\n";
	#PAUP keeps some names flanked by single quotes in the tree description
    #PAML will have trouble to reconcile these names
    @arr = split (/\'/, $trees[0]); 

	#number of elements in teh string, once split by single quotes
    $noElem=@arr;                   

    for (my $i=0; $i <$noElem; $i++) { 
		#reconstitute string without single quotes
		$theTree .= "$arr[$i]";   
    }
    
# now, to unroot the tree...

    @arr = split (//, $theTree); 

    $lineLen = @arr;
    for (my $i = 0; $i < $lineLen; $i++ ) { 
		if ($arr[$i] eq "(") {
	    	$parenthesisCounter++;
	    
	    	if ($parenthesisCounter == 2) {
				$openingParenthesis = $i;
	    	}
		}
		if ($arr[$i] eq ")") {
	    	$parenthesisCounter--;
	    	if ($parenthesisCounter == 1) { 
				$closingParenthesis = $i;
				last;
	    	}
		}
    }

    $theTree = "";
    for (my $i = 0; $i < $openingParenthesis; $i++ ) { 
		$theTree .= $arr[$i];
    }
    for (my $i = $openingParenthesis + 1; $i < $closingParenthesis; $i++ ) { 
		$theTree .= $arr[$i];
    }
    for (my $i = $closingParenthesis + 1; $i < $lineLen; $i++ ) { 
		$theTree .= $arr[$i];
    }
	#finished unrooting

    print TREEFILEout "$theTree";  
    print TREEFILEout "\n";

    close TREEFILEout;

	## THERE WILL PROBABLY BE SOME PATH PROBLEMS HERE
    chmod(0744, $treeout);
}

sub makePAUPbatch {

    my $batchfile = $options{'temp_path'}."PAUPbatch.tmp";
    my $logfile = "$_[0].PAUPlog";
    my $infile= "$_[0].nex";
    my $treefile = "$_[0].PAUPtree";

	## THERE WILL PROBABLY BE SOME PATH PROBLEMS HERE
	open (BATCHFILE, ">$batchfile") or die $!;

	print BATCHFILE "begin paup;\n";
#	print BATCHFILE "log file=$logfile;\n";
##	Force overwrite of log file if we are re-running
	print BATCHFILE "log file=$logfile replace=yes;\n";
	print BATCHFILE "execute $infile;\n";
	print BATCHFILE "set maxtrees=100 increase=auto criterion=parsimony taxlabels=full;\n";
	print BATCHFILE "bandb;\n";
	print BATCHFILE "savetrees replace=yes format=altnex root=no file=$treefile;\n";
	print BATCHFILE "end;\n";
	print BATCHFILE "quit;\n";
	
	close BATCHFILE;
	return $batchfile;
}


sub readSpeciesNames {   #gets passed a nexus file

	my $infile = shift @_;
	
    my $count=0;
    my @dummyarr;
    my @arr;

	## THERE WILL PROBABLY BE SOME PATH PROBLEMS HERE
	open (NEXUS, "$infile") or die $!;
    while (<NEXUS>) {
		chomp;
		@arr= split(' ', $_);
		if ($arr[0]) {
	    	last if ($arr[0] eq "TAXLABELS");
		}
    }
    while (<NEXUS>) {
		chomp;
		@arr= split(' ', $_);
		last if ($arr[0] eq ";");
		$speciesNames[$count++]= substr($arr[0], 1, (length($arr[0]) - 2) );
    }
    while (<NEXUS>) {
		chomp;
		@arr= split(' ', $_);
		if ($arr[1]) {
		    if ($arr[1] =~ /^NCHAR/) {
				@dummyarr = split(/=/, $arr[1]);
				$aaAlignLength = substr($dummyarr[1], 0, (length($dummyarr[1])-1) );
				last;
	    	}
		}
    }
    close NEXUS;
}
