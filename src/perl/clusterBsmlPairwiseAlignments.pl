#!/local/perl/bin/perl
=head1 NAME

clusterBsmlPairwiseAlignments.pl - 

=head1 SYNOPSIS

USAGE:  clusterBsmlPairwiseAlignments.pl  -b bsml_list -m match_list -k linkscore -p percent_identity -u p_value [-l log] [-h help] [-o outdir] [-p percent_identity] [-u p_value]
=head1 OPTIONS

=over 8

=item B<--bsml_list,-b>
    
    List of bsml files containing proteins and assemblies

=item B<--match_list,-m>
    
    List of bsml files containing pairwise matches

=item B<--log,-l>
    
  log file

=item B<--linkscore,-k>

    Link score: <0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9>

=item B<--outdir, -o>
    
    Output directory for the Jaccard cluster files. Defaults to current directory

=item B<--percent_identity, -p>
    
    Percent identity threshold value.  Default is set to 75 

=item B<--p_value, -u>
    
    P_value threshold value.  Default is set to 1e-15

=item B<--help,-h>

    Print this help

=item B<--debug,-d>

    Debug level

=back

=head1 DESCRIPTION

    clusterBsmlPairwiseAlignments.pl - Read BSML pairwise search encoding documents and produce Jaccard cluster file(s)

=cut

#
# chado.analysisfeature fields          ==> BSML attributes
# 
# rawscore double precision NULL        ==> runscore
# normscore double precision NULL       ==> percent_bit_score
# significance double precision NULL    ==> p_value
# pidentity double precision NULL       ==> percent_identity
#
#


use strict;
use Jaccard_coefficient_cluster_resolver;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use BSML::BsmlParserSerialSearch;
use Workflow::Logger;

my %options = ();
my $results = GetOptions (\%options,
			  'bsmlSearchList|m=s',
			  'bsmlModelList|b=s',
			  'linkscore|k=s',
			  'percent_identity|p=s',
			  'p_value|u=s',
			  'outfile|o=s',
			  'log|l=s',
			  'debug|d=s',
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

#MAIN HERE


my $valid_asmbls = &retrieve_asmbls(&get_list_from_file($options{'bsmlModelList'}));

my $pairs = &retrieve_protein_pairs(
				    bsmldoc_list     => &get_list_from_file($options{'bsmlSearchList'}),
				    percent_identity => $options{'percent_identity'},
				    p_value          => $options{'p_value'},
				    valid_asmbls     => $valid_asmbls
				    );


&produce_cluster_output(
			pairs      => $pairs,
			output     => $options{'outfile'},
			linkscore => $options{'linkscore'},
			);




$logger->info("'$0': Finished parsing bsml documents and producing clusters");
$logger->info("Please verify log4perl log file: $options{'log'}");


#---------------------------------------------------------------------------------------------------------------------
#
#                           END OF MAIN  -- SUBROUTINES FOLLOW
#
#---------------------------------------------------------------------------------------------------------------------

sub retrieve_asmbls{
    my ($bsmldoclist) = @_;

    my $valid_asmbls = {};

    foreach my $bsmldoc (@$bsmldoclist){
	$logger->logdie("bsmldoc was not defined") if (!defined($bsmldoc));

	$logger->info("Processing bsml document: $bsmldoc");

	#----------------------------------------------------------
	# New serial parsing manner in which analysis component
	# and seq-pair-alignment components will be serially
	# parsed via callback methods.
	#
	#-----------------------------------------------------------

	print ("Parsing Sequence component for bsml document: $bsmldoc\nAnd extracting assembly-protein pairs\n");
	my $bsml_parser = new BSML::BsmlParserSerialSearch(ReadFeatureTables => 0,
							   SequenceCallBack  => sub {
							       my ($sequence_ref) = @_;

							       my $assembly_id = $sequence_ref->{'BsmlAttr'}->{'ASSEMBLY'} if ((exists $sequence_ref->{'BsmlAttr'}->{'ASSEMBLY'}) and (defined($sequence_ref->{'BsmlAttr'}->{'ASSEMBLY'})));
							       my $protein_id = $sequence_ref->{'attr'}->{'id'} if ((exists $sequence_ref->{'attr'}->{'id'}) and (defined($sequence_ref->{'attr'}->{'id'})));
							       $logger->logdie("assembly_id was not defined") if (!defined($assembly_id));
							       $logger->logdie("protein_id was not defined") if (!defined($protein_id));


									       
							       #
							       # Store only unique protein-assembly pairs
							       #
							       if ((!exists $valid_asmbls->{$protein_id})){
								   $valid_asmbls->{$protein_id} = $assembly_id;
							       }
							       #
							       # Inform the user that something fishy about their BSML search encoding document...
							       #
							       else{
								   $logger->warn("FYI protein:$protein_id was previously encountered with assembly $valid_asmbls->{$protein_id}...") if $logger->is_warn();
							       }
							   }
							   );
	
	$logger->logdie("bsml_parser was not defined") if (!defined($bsml_parser));
	$bsml_parser->parse($bsmldoc);
    }

    return $valid_asmbls;
}


#------------------------------------------------------------
# produce_cluster_output()
#
#------------------------------------------------------------
sub produce_cluster_output {

    $logger->debug("Entered produce_cluster_output") if $logger->is_debug();

    my %param = @_;
    my $paramhash = \%param;


    my $pairs      = $paramhash->{'pairs'}      if ((exists $paramhash->{'pairs'}) and (defined($paramhash->{'pairs'})));
    my $output     = $paramhash->{'output'}     if ((exists $paramhash->{'output'}) and (defined($paramhash->{'output'})));
    my $linkscore = $paramhash->{'linkscore'} if ((exists $paramhash->{'linkscore'}) and (defined($paramhash->{'linkscore'})));


    $logger->logdie("pairs was not defined")      if (!defined($pairs));
    $logger->logdie("output was not defined")     if (!defined($output));
    $logger->logdie("linkscore was not defined") if (!defined($linkscore));


    #
    # BSML document has been parsed and the array of array references "pairs" has been created
    # Now, for each link score, produce clusters and write to file by same name of link score
    #

    $logger->logdie("linkscore was not defined") if (!defined($linkscore));
	
    #
    # 1. Get clustering
    #
    
    my $resolver = new Jaccard_coefficient_cluster_resolver($linkscore);
    $logger->logdie("resolver was not defined") if (!defined($resolver));
    
    my @clusters = $resolver->resolve_clusters(@$pairs);
    
    if (@clusters){
	#
	# 2. Open and write to cluster output file
	#
	open (OUTFILE, ">$output") or $logger->logdie("Could not open $output in write mode");
	
	my $clustercount=0;

	foreach my $cluster (@clusters){
	    $logger->logdie("cluster was not defined") if (!defined($cluster));
	    
	    my $size = scalar(@$cluster);
	    my $clustername = "jaccard_$clustercount";

	    print OUTFILE "COG = $clustername, size $size, connections  = -1, perfect = -1;";
	    
	    foreach my $member (@$cluster){
		print OUTFILE "\t$member\n";
	    }
	    $clustercount++;
	}
	close OUTFILE;
    }
}

#-------------------------------------------------------------------------
# retrieve_protein_pairs()
#
#-------------------------------------------------------------------------
sub retrieve_protein_pairs {

    $logger->debug("Entered retrieve_protein_pairs") if $logger->is_debug();

    my (%param) = @_;
    my $paramref = \%param;
    

    #
    # Receive, parse, extract parameters and then verify whether defined
    #
    my $bsmldoclist      = $paramref->{'bsmldoc_list'}       if ((exists $paramref->{'bsmldoc_list'}) and (defined($paramref->{'bsmldoc_list'})));
    my $percent_identity = $paramref->{'percent_identity'}   if ((exists $paramref->{'percent_identity'}) and (defined($paramref->{'percent_identity'})));
    my $p_value          = $paramref->{'p_value'}            if ((exists $paramref->{'p_value'}) and (defined($paramref->{'p_value'})));
    my $protein2assemblyhash     = $paramref->{'valid_asmbls'}            if ((exists $paramref->{'valid_asmbls'}) and (defined($paramref->{'valid_asmbls'})));

    $logger->logdie("bsmldoclist was not defined")       if (!defined($bsmldoclist));
    $logger->logdie("percent_identity was not defined")  if (!defined($percent_identity));
    $logger->logdie("p_value was not defined")           if (!defined($p_value));
    $logger->logdie("valid_asmbls was not defined")           if (!defined($valid_asmbls));


    my @proteinpairs;

    #
    # for each bsml document we will:
    # 1) verify access permissions
    # 2) validate (if told to)
    # 3) parse the document
    # 4) retrieve all pairs of protein identifiers for which both proteins belong to the same assembly
    #
    foreach my $bsmldoc (@$bsmldoclist){
	$logger->logdie("bsmldoc was not defined") if (!defined($bsmldoc));

	$logger->info("Processing bsml document: $bsmldoc");
    
	print ("Parsing Seq-pair-alignments for bsml document: $bsmldoc\nAnd extract protein-protein pairs\n");
	
	my $bsml_parser = new BSML::BsmlParserSerialSearch(
							AlignmentCallBack  => sub {
							    my ($alignment_ref) = @_;
							    

							    my $compseq = $alignment_ref->{'attr'}->{'compseq'} if ((exists $alignment_ref->{'attr'}->{'compseq'}) and (defined($alignment_ref->{'attr'}->{'compseq'})));
							    my $refseq = $alignment_ref->{'attr'}->{'refseq'} if ((exists $alignment_ref->{'attr'}->{'refseq'}) and (defined($alignment_ref->{'attr'}->{'refseq'})));
							    
							    $logger->logdie("compseq was not defined") if (!defined($compseq));
							    $logger->logdie("refseq was not defined")  if (!defined($refseq));      

							    #
							    # Retrieve the percent_identity and the p_value from the BSML document
							    # 
							    my $pidentity = $alignment_ref->{'BsmlSeqPairRuns'}->[0]->{'BsmlAttr'}->{'percent_identity'} if ((exists $alignment_ref->{'BsmlSeqPairRuns'}->[0]->{'BsmlAttr'}->{'percent_identity'}) and (defined($alignment_ref->{'BsmlSeqPairRuns'}->[0]->{'BsmlAttr'}->{'percent_identity'})));
							    $logger->logdie("pidentity was not defined") if (!defined($pidentity));

							    my $pvalue = $alignment_ref->{'BsmlSeqPairRuns'}->[0]->{'BsmlAttr'}->{'p_value'} if ((exists $alignment_ref->{'BsmlSeqPairRuns'}->[0]->{'BsmlAttr'}->{'p_value'}) and (defined($alignment_ref->{'BsmlSeqPairRuns'}->[0]->{'BsmlAttr'}->{'p_value'})));
							    $logger->logdie("pvalue was not defined") if (!defined($pvalue));

							    
#							    print "compseq:$compseq\t\trefseq:$refseq\n";
							    
							    #
							    # We only keep the protein pairs that meet the following conditions:
							    # 1) the assembly is one in specified assembly list/hash
							    # 2) the percent_identity is above the threshold value
							    # 3) the p_value is above the threshold value
							    #
							    #
							    #
							    if ((exists $protein2assemblyhash->{$compseq}) and (defined($protein2assemblyhash->{$compseq})) and (exists $protein2assemblyhash->{$refseq}) and (defined($protein2assemblyhash->{$refseq}))){
								
								#
								#  pvalue is from blastp bsml file ()
								#  p_value threshold value

								
								if (($pidentity > $percent_identity) and ($pvalue < $p_value)){
								    push (@proteinpairs, [$compseq, $refseq]);# if ($protein2assemblyhash->{$compseq} eq $protein2assemblyhash->{$refseq});
								}
							    }
							    
							}
							);
	
	$logger->logdie("bsml_parser was not defined") if (!defined($bsml_parser));
	$bsml_parser->parse($bsmldoc);


    }

    $logger->debug("Protein pairs to be processed:\n") if $logger->is_debug();

    return \@proteinpairs;

}

sub get_list_from_file{
    my($file) = @_;
    my @lines;
    open( FH, $file ) or $logger->logdie("Could not open $file");
    while( my $line = <FH> ){
	chomp($line);
	push @lines, split(',',$line) if($line =~ /\S+/);
    }
    return \@lines;
}

sub check_parameters{
    my ($options) = @_;
    
    if(0){
	pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});    
    }
}
