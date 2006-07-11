#!/local/perl/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";
=head1 NAME

clusterBsmlPairwiseAlignments.pl - 

=head1 SYNOPSIS

USAGE:  clusterBsmlPairwiseAlignments.pl  -b bsml_list -m match_list -k linkscore -p percent_identity -u p_value [-l log] [-h help] [-o outdir] [-p percent_identity] [-u p_value]
=head1 OPTIONS

=over 8

=item B<--bsml_list,-b>
    
    List of bsml files containing polypeptides and assemblies

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
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Workflow::Logger;
use XML::Parser;
use Jaccard_coefficient_cluster_resolver;
use Data::Dumper;
use DB_File;

my %options = ();
my $results = GetOptions (\%options,
			  'bsmlSearchList|m=s',
			  'asmbl_lookup|a=s',
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


my $valid_asmbls = build_asmbl_lookup($options{'asmbl_lookup'});


my $pairs = &retrieve_polypeptide_pairs(
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
    
    #
    # 2. Open and write to cluster output file
    #
    open (OUTFILE, ">$output") or $logger->logdie("Could not open $output in write mode");
    
    my $clustercount=0;
    
    foreach my $cluster (@clusters){
	$logger->logdie("cluster was not defined") if (!defined($cluster));
	
	my $size = scalar(@$cluster);
	my $clustername = "jaccard_$clustercount";
	
	print OUTFILE "COG = $clustername, size $size, connections  = -1, perfect = -1;\n";
	
	foreach my $member (@$cluster){
	    print OUTFILE "\t$member\n";
	}
	$clustercount++;
    }
    close OUTFILE;
}

sub process_alignment{
    my($polypeptidepairs,$polypeptide2assemblyhash,$compseq,$refseq,$pidentity,$pvalue,$pidentity_cutoff,$pvalue_cutoff) = @_;
    $logger->logdie("compseq was not defined") if (!defined($compseq));
    $logger->logdie("refseq was not defined")  if (!defined($refseq));
    $logger->logdie("pidentity was not defined") if (!defined($pidentity));
    $logger->logdie("pvalue was not defined") if (!defined($pvalue));
    #
    # We only keep the polypeptide pairs that meet the following conditions:
    # 1) the assembly is one in specified assembly list/hash
    # 2) the percent_identity is above the threshold value
    # 3) the p_value is above the threshold value
    #
    #
    #
    if ((exists $polypeptide2assemblyhash->{$compseq}) and (defined($polypeptide2assemblyhash->{$compseq})) and (exists $polypeptide2assemblyhash->{$refseq}) and (defined($polypeptide2assemblyhash->{$refseq}))){
	
	#
	#  pvalue is from blastp bsml file ()
	#  p_value threshold value
	
	
	if (($pidentity > $pidentity_cutoff) and ($pvalue < $pvalue_cutoff)){
	    push (@$polypeptidepairs, [$compseq, $refseq]);
	}
    }
    
}

#-------------------------------------------------------------------------
# retrieve_polypeptide_pairs()
#
#-------------------------------------------------------------------------
sub retrieve_polypeptide_pairs {

    $logger->debug("Entered retrieve_polypeptide_pairs") if $logger->is_debug();

    my (%param) = @_;
    my $paramref = \%param;
    

    #
    # Receive, parse, extract parameters and then verify whether defined
    #
    my $bsmldoclist      = $paramref->{'bsmldoc_list'}       if ((exists $paramref->{'bsmldoc_list'}) and (defined($paramref->{'bsmldoc_list'})));
    my $percent_identity = $paramref->{'percent_identity'}   if ((exists $paramref->{'percent_identity'}) and (defined($paramref->{'percent_identity'})));
    my $p_value          = $paramref->{'p_value'}            if ((exists $paramref->{'p_value'}) and (defined($paramref->{'p_value'})));
    my $polypeptide2assemblyhash     = $paramref->{'valid_asmbls'}            if ((exists $paramref->{'valid_asmbls'}) and (defined($paramref->{'valid_asmbls'})));

    $logger->logdie("bsmldoclist was not defined")       if (!defined($bsmldoclist));
    $logger->logdie("percent_identity was not defined")  if (!defined($percent_identity));
    $logger->logdie("p_value was not defined")           if (!defined($p_value));
    $logger->logdie("valid_asmbls was not defined")           if (!defined($valid_asmbls));


    my @polypeptidepairs;

    #
    # for each bsml document we will:
    # 1) verify access permissions
    # 2) validate (if told to)
    # 3) parse the document
    # 4) retrieve all pairs of polypeptide identifiers for which both polypeptides belong to the same assembly
    #
    my $compseq = undef;
    my $refseq = undef;
    my $pidentity = undef;
    my $pvalue = undef;
    
    my $funcs = {'Seq-pair-alignment'=>
		     sub {
			 my ($expat,$elt,%params) = @_;
			 &process_alignment(\@polypeptidepairs,$polypeptide2assemblyhash,$compseq,$refseq,$pidentity,$pvalue,$percent_identity,$p_value) if(defined $compseq && defined $refseq);
			 $compseq = undef;
			 $refseq = undef;
			 $pidentity = undef;
			 $pvalue = undef;
			 $compseq = $params{'compseq'} if ((exists $params{'compseq'}) and (defined($params{'compseq'})));
			 $refseq = $params{'refseq'} if ((exists $params{'refseq'}) and (defined($params{'refseq'})));
			 
		     },
		 'Attribute'=>
		     sub {
			 my ($expat,$elt,%params) = @_;
			 my $index = scalar(@{$expat->{'Context'}}) - 1;
			 if($expat->{'Context'}->[$index] eq 'Seq-pair-run'){
			     if($params{'name'} eq 'percent_identity'){
				 $pidentity = $params{'content'};
			     }
			     if($params{'name'} eq 'p_value'){
				 $pvalue = $params{'content'};
			     }
			 }
		     }
	     };

    my $x = new XML::Parser(Handlers => 
                            {
                                Start =>
                                    sub {
                                        #$_[1] is the name of the element
                                        if(exists $funcs->{$_[1]}){
                                            $funcs->{$_[1]}(@_);
                                        }
                                    }
                                }
                            );



    foreach my $bsmldoc (@$bsmldoclist){
	$logger->logdie("bsmldoc was not defined") if (!defined($bsmldoc));

	$logger->info("Processing bsml document: $bsmldoc");
	
	$x->parsefile( $bsmldoc );
    }

    &process_alignment(\@polypeptidepairs,$polypeptide2assemblyhash,$compseq,$refseq,$pidentity,$pvalue,$percent_identity,$p_value) if(defined $compseq && defined $refseq);

    $logger->debug("Polypeptide pairs to be processed:\n") if $logger->is_debug();

    return \@polypeptidepairs;

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

sub build_asmbl_lookup{
    my $reader = shift;
    my %lookup;
    $logger->debug("Reading lookup $reader") if($logger->is_debug());
    tie %lookup, 'DB_File', $reader, O_RDWR|O_CREAT, 0660, $DB_BTREE or $logger->logdie("Can't tie $reader $!");
    $logger->debug("Found ".scalar(keys %lookup)." keys") if($logger->is_debug());
    return \%lookup;
}


sub check_parameters{
    my ($options) = @_;
    
    if(0){
	pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});    
    }
}
