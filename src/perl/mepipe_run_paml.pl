#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME

mepipe_run_paml.pl - molecular evolution pipeline script that executes PAML

=head1 SYNOPSIS

USAGE: mepipe_run_paml.pl
        --input_path=/path/to/directory/containing/input/files/
        --temp_path=/path/to/store/temporary/output
        --output_path=/path/to/store/output
		--paml_bin_dir=/path/to/paml/bin
		--dn_ds_ratio=0.1
		[--require_more_than_two]
        [--debug=4]

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
	
B<--paml_bin_dir,-b>
	Full path containing PAML executables (specifically codeml).
	
B<--dn_ds_ratio,-w>
	dN/dS ratio estimate (w).

B<--require_more_than_two,-3>
	Will not execute if the number of input sequences is less than three.
	
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
    import Ergatis::Logger;
}


my %options = ();
my $results = GetOptions (\%options,
						  'input_prefix|p=s',
                          'input_path|i=s',
						  'temp_path|t=s',
                          'output_path|o=s',
						  'paml_bin_dir|b=s',
						  'dn_ds_ratio|w=s',
						  'require_more_than_two|3',
#						  'paup_bin_dir=s',
#						  'emboss_bin_dir=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = Ergatis::Logger::get_logger();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDERR} );
}

#&check_parameters(\%options);

my $aaAlign;
my $nuclSeqs;
my $nuclAlign;
my $geneName;

my $outgroup;
my $batchfile;
my @speciesNames;
my $noSpecs;
my $aaAlignLength;
my $PAMLrunmode;

foreach my $dir((
		'paml_bin_dir', 
#		'emboss_bin_dir', 
#		'paup_bin_dir',
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

chdir($options{'temp_path'});
	
$geneName = $options{'input_prefix'};
$aaAlign = "$geneName.aln";
$nuclSeqs = "$geneName.seq";
$nuclAlign = "$geneName.nuc.aln";

$outgroup = "";
$batchfile = "";
@speciesNames = ("");
$noSpecs = 0;

&readSpeciesNames($options{'temp_path'}."$aaAlign.nex");   
#seq names are flanked by single quotes
#list of seq names is writen to array @speciesNames
# checks length of aa alignments and stores in $aaAlignLength
$noSpecs = @speciesNames;

my $wstring = "w".$options{'dn_ds_ratio'};
$wstring =~ s/\./_/g;

if ($noSpecs > 2) { 
	# with 3 species (or more), calculations are tree-based
	# 
	# Trees were either created manually (<= 3 species) 
	# or generated using PAUP -- either way, they have been
	# already created in the temp_dir by 'mepipe_create_tree.pl'
	# 
	$PAMLrunmode = 0; 
#    &makeCtlFile ($options{'temp_path'}."$geneName", 0.1, $PAMLrunmode, "0 1 2 3 7 8", 10); #run PAML with start value for w = 0.1 (w = dN/dS)
    &makeCtlFile ($options{'temp_path'}."$geneName", $options{'dn_ds_ratio'}, $PAMLrunmode, "0 1 2 3 7 8", 10);

#   print "Running first pass of PAML\n";
	## RUNNING PAML/CODEML IN THIS WAY IS PROBABLY GOING TO BE A PROBLEM
	system($options{'paml_bin_dir'}."codeml"); #run paml
	## THERE WILL PROBABLY BE SOME PATH PROBLEMS HERE
	move($options{'temp_path'}."rst", $options{'temp_path'}."$geneName.$wstring.rst");
	## THERE WILL PROBABLY BE SOME PATH PROBLEMS HERE
	move($options{'temp_path'}."lnf", $options{'temp_path'}."$geneName.$wstring.lnf");  
#	print "Finished first pass of PAML\n";
	    
#	&makeCtlFile ($options{'temp_path'}."$geneName", 1.5, $PAMLrunmode, "0 1 2 3 7 8", 10); #run PAML with start value for w = 1.5 (w = dN/dS)

#	print "Running second pass of PAML\n";
	## RUNNING PAML/CODEML IN THIS WAY IS PROBABLY GOING TO BE A PROBLEM
#	system($options{'paml_bin_dir'}."codeml"); #run paml
	## THERE WILL PROBABLY BE SOME PATH PROBLEMS HERE
#	move($options{'temp_path'}."rst", $options{'temp_path'}."$geneName.w01.rst");
	## THERE WILL PROBABLY BE SOME PATH PROBLEMS HERE
#	move($options{'temp_path'}."lnf", $options{'temp_path'}."$geneName.w01.lnf");  
#	print "Finished second pass of PAML\n";

} elsif (! defined($options{'require_more_than_two'})) { 
	$PAMLrunmode = "-2";
	&makeCtlFile ($options{'temp_path'}."$geneName", $options{'dn_ds_ratio'}, $PAMLrunmode, 0, 1); #run PAML with start value for w = 0.1 (w = dN/dS)
	system($options{'paml_bin_dir'}."codeml"); #run paml
	# with 2 species, calculations are pairwise (rumnode = -2)
	## THERE WILL PROBABLY BE SOME PATH PROBLEMS HERE
	move($options{'temp_path'}."rst", $options{'temp_path'}."$geneName.$wstring.rst");
	## THERE WILL PROBABLY BE SOME PATH PROBLEMS HERE
	move($options{'temp_path'}."lnf", $options{'temp_path'}."$geneName.$wstring.lnf");  
}	


## Cleanup steps should be performed as a final step in component

#chdir($options{'temp_path'});
## THERE WILL PROBABLY BE SOME PATH PROBLEMS HERE
#unlink <*.PAUPlog>;
## THERE WILL PROBABLY BE SOME PATH PROBLEMS HERE
#unlink <*.PAUPtree>;
## THERE WILL PROBABLY BE SOME PATH PROBLEMS HERE
#unlink <*.aln.nex>;

### THESE MOVE STEPS SHOULD BE A SEPARATE STEP PROBABLY SO THAT IF THEY
### FAIL WE DON'T HAVE TO RERUN PAML

## MOVE THE .PAMLout FILE TO THE OUTPUT DIRECTORY
move($options{'temp_path'}."$geneName.$wstring.PAMLout", $options{'output_path'}."$geneName.$wstring.PAMLout");

## MOVE THE .rst FILE TO THE OUTPUT DIRECTORY
move($options{'temp_path'}."$geneName.$wstring.rst", $options{'output_path'}."$geneName.$wstring.rst");


# end main
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------


sub makeCtlFile  {   #gets passed geneName, starting w value (dN/dS), runmode, evolution model(s) and # site classes
    my $seqfile = "$_[0].nuc.aln.PAMLseq";
    my $treefile = "$_[0].PAMLtree";
    my $wStart = $_[1];   #starting value for omega
    my $NSsites = $_[3];
    my $ncatG = $_[4]; 

    my $ctlfile = $options{'temp_path'}."codeml.ctl";

	## THERE WILL PROBABLY BE SOME PATH PROBLEMS HERE
    open (CTLout, ">$ctlfile") or die $!;#

    print CTLout "      seqfile = $seqfile\n";
    print CTLout "     treefile = $treefile\n";
	my $wstring = "w".$_[1];
	$wstring =~ s/\./_/g;
    print CTLout "      outfile = $_[0].$wstring.PAMLout\n";
    print CTLout "        noisy = 0 \n"; # 0,1,2,3,9: how much rubbish on the screen
    print CTLout "      verbose = 1 \n"; # 0: concise; 1: detailed, 2: too much
    print CTLout "      runmode = $_[2] \n"; # * 0: user tree;  -2: pairwise\n";

    print CTLout "      seqtype = 1 \n"; # * 1:codons; 2:AAs; 3:codons-->AAs
    print CTLout "    CodonFreq = 2 \n"; # * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
    print CTLout "       aaDist = 0 \n"; # * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
    print CTLout "   aaRatefile = wag.dat\n"; # * only used for aa seqs with model=empirical(_F)
#                                         * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own

    print CTLout "        model = 0\n"; # same dN/dS model for all branches
                                        # models for codons:
                                        #  0:one, 1:b, 2:2 or more dN/dS ratios for branches

    print CTLout "      NSsites = $NSsites\n"; 
#                    0:one w;1:neutral;2:selection; 3:discrete;4:freqs
#                    5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma
#                   10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1; 13:3normal>0

    print CTLout "        icode = 0 \n"; # * 0:universal code; 1:mammalian mt; 2-10:see below
    print CTLout "        Mgene = 0 \n"; # * 0:rates, 1:separate; 

    print CTLout "    fix_kappa = 0 \n"; # * 1: kappa fixed, 0: kappa to be estimated
    print CTLout "        kappa = 2 \n"; # * initial or fixed kappa
    print CTLout "    fix_omega = 0 \n"; #* 1: omega or omega_1 fixed, 0: estimate 
    print CTLout "        omega = $wStart \n"; #* initial or fixed omega, for codons or codon-based AAs
    print CTLout "    fix_alpha = 1 \n"; # * 0: estimate gamma shape parameter; 1: fix it at alpha
    print CTLout "        alpha = 0 \n"; # * initial or fixed alpha, 0:infinity (constant rate)
    print CTLout "       Malpha = 0  \n"; # * different alphas for genes
    print CTLout "        ncatG = 1  \n"; # * no. of categories in dG of NSsites models


    print CTLout "        clock = 0 \n"; # * 0:no clock, 1:global clock; 2:local clock; 3:TipDate
    print CTLout "        getSE = 0  \n"; #* 0: don\'t want them, 1: want S.E.s of estimates
    print CTLout " RateAncestor = 0 \n"; # * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
    print CTLout "    Small_Diff = .5e-6\n";
    print CTLout "*    cleandata = 1 \n"; # * remove sites with ambiguity data (1:yes, 0:no)?
    print CTLout "*        ndata = 1\n";
    print CTLout "*  fix_blength = -1  \n"; #* 0: ignore, -1: random, 1: initial, 2: fixed
    print CTLout "        method = 0  \n"; # * 0: simultaneous; 1: one branch at a time
    close CTLout;
}


sub readSpeciesNames {   #gets passed a nexus file
	my $infile = shift @_;
    my $count=0;
    my @dummyarr;
    my @arr;

	## THERE WILL PROBABLY BE SOME PATH PROBLEMS HERE
	#print $_[0]."\n";
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

