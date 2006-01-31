#!/usr/local/bin/perl

=head1  NAME 

berbtab2bsml.pl  - convert info stored in BER btab files into BSML documents 

=head1 SYNOPSIS

USAGE:  berbtab2bsml.pl -b btab_dir -o blastp.bsml -d bsml_dir

=head1 OPTIONS

=over 4

=item *

B<--bsml_dir,-d> [REQUIRED]  Dir containing BSML documents (repository)

=item *

B<--output,-o> [REQUIRED] output BSML file containing bit_score information

=item * 

B<--btab_dir,-b> [REQUIRED] Dir containing btab files

=item *

B<--max_hsp_count,-m> [REQUIRED] Maximum number of HSPs stored per alignment
    
=item *

B<--class,-c> The ref/comp sequence type.  Default is 'assembly'

=item *

B<--help,-h> This help message

=back

=head1   DESCRIPTION

berbtab2bsml.pl is designed to convert information in ber btab files into BSML documents.

Samples:

   btab2bsml.pl -d bsml_dir -b /usr/btab -o blastp.bsml

NOTE:  

Calling the script name with NO flags/options or --help will display the syntax requirement.


=cut


use strict;
use Log::Log4perl qw(get_logger);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use English;
use BSML::BsmlBuilder;
use BSML::BsmlReader;
use BSML::BsmlParserTwig;
use BSML::BsmlRepository;
use File::Basename;
use File::Path;
use Pod::Usage;
use Workflow::Logger;

my %options = ();
my $results = GetOptions (\%options, 
			  'btab_dir|b=s', 
			  'btab_file|f=s',
			  'bsml_dir|d=s', 
			  'output|o=s', 
			  'max_hsp_count|m=s',
			  'pvalue|p=s', 
			  'log|l=s',
			  'debug=s',
			  'class|c=s',
			  'help|h') || pod2usage();

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

if($options{'pvalue'} eq ""){
    $options{'pvalue'} = 10;
}

my $class;
if (!defined($options{'class'})){
    $logger->logdie("class was not defined");
}
else{
    $class = $options{'class'};
}



# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

&check_parameters(\%options);

my $files = &get_btab_files($options{'btab_dir'},$options{'btab_file'});

my $doc = new BSML::BsmlBuilder();

#generate lookups
my ($gene_asmbl_id, $protID_cdsID);
#if($options{'btab_type'} == 1) {
($gene_asmbl_id, $protID_cdsID) = build_id_lookups_ber($options{'bsml_dir'});
$doc->makeCurrentDocument();
parse_ber_btabs($files);
#}
#elsif($options{'btab_type'} == 2) {
#    $gene_asmbl_id = build_id_lookups_blast($options{'bsml_dir'});
#    $doc->makeCurrentDocument();
#    parse_blast_btabs($files);
#} else {
#    print STDERR "Bad btab_type.  Aborting...\n";
#    exit 5;
#}

$doc->write($options{'output'});
exit();   

sub parse_ber_btabs {

    my $btab_files = shift;
    my $num;
    foreach my $file (@$btab_files) {
	$num++; 
	open (BTAB, "$file") or die "Unable to open \"$file\" due to $!";
	$logger->debug("opening $file $num using pvalue cutoff $options{'pvalue'}") if($logger->is_debug());
        my (@btab, $seq, $query_name, $match_name, $query_cds_id, $query_protein_id);
	while(my $line = <BTAB>) {
	    chomp($line);
	    @btab = split("\t", $line);
	    next if($btab[20] > $options{'pvalue'} );
	    next if($btab[3] ne 'praze');
	    next if( !($btab[13] > 0) || !($btab[14] > 0) );
	    next if(!$btab[0] or !$btab[5]);
	    $btab[5] =~ s/\|//g;   #get rid of trailing |
	    $query_protein_id = $btab[0];
	    $btab[0] = $protID_cdsID->{$btab[0]} if $protID_cdsID->{$btab[0]};
	    $query_cds_id = $btab[0];	    
	    $match_name = $btab[5];
	    splice(@btab, 19, 1);
	    
	    for (my $i=0;$i<scalar(@btab);$i++){
		if ($btab[$i] eq 'N/A'){
		    $btab[$i] = undef;
		}
	    }
	    my $align = &createAndAddBtabLine(
					      doc                => $doc,
					      query_name         => $btab[0],
					      date               => $btab[1],
					      query_length       => $btab[2],
					      blast_program      => $btab[3],
					      search_database    => $btab[4],
					      dbmatch_accession  => $btab[5],
					      start_query        => $btab[6],
					      stop_query         => $btab[7],
					      start_hit          => $btab[8],
					      stop_hit           => $btab[9],
					      percent_identity   => $btab[10],
					      percent_similarity => $btab[11],
					      bit_score          => $btab[12],
					      chain_number       => $btab[13],
					      segment_number     => $btab[14],
					      dbmatch_header     => $btab[15],
					      unknown1           => $btab[16],
					      unknown2           => $btab[17],
					      e_value            => $btab[18],
					      p_value            => $btab[19]
					      );

	    $seq = $doc->returnBsmlSequenceByIDR($match_name);
	    my $match_asmbl_id = $gene_asmbl_id->{$match_name};

	}
	close BTAB;
	next if(!defined($query_cds_id));
	$seq = $doc->returnBsmlSequenceByIDR($query_cds_id);
	if($seq) {
	    my $query_asmbl_id = $gene_asmbl_id->{$query_protein_id};
	}
    }
}



    
sub build_id_lookups_ber {

    my $BSML_dir = shift;

    my $bsmlrepo = new BSML::BsmlRepository('BSML_repository'=>$BSML_dir);
 
    my ($files) = $bsmlrepo->list_bsml_files();

    my $parser = new BSML::BsmlParserTwig;

    my $protID_cdsID = {};
    my $gene_asmbl_id ={};
    my $reader;
    foreach my $bsml_doc (@$files) {
	if (-s $bsml_doc) {
	    $logger->debug("parsing $bsml_doc") if($logger->is_debug());
	    $reader = BSML::BsmlReader->new();
	    $parser->parse( \$reader, $bsml_doc );
	    my $hash_ref = $reader->get_all_protein_assemblyId();
	    while(my ($protein_id, $asmbl_id) = each %$hash_ref) {
		$gene_asmbl_id->{$protein_id} = $asmbl_id;
	    }
	    my $rhash = $reader->returnAllIdentifiers();
	    build_protID_cdsID_mapping($rhash, $protID_cdsID); 	
	} else {  
	    $logger->("Empty $bsml_doc...skipping") if($logger->is_debug());
	    }
    }

    return ($gene_asmbl_id, $protID_cdsID);

}


sub build_protID_cdsID_mapping {
#This function builds a mapping between cdsID to proteinID. 
#The returned structure is a hash ref, where key is protID, value is cdsID

    my $rhash = shift;
    my $protID_cdsID=shift;

    foreach my $seqID (keys %$rhash) {
	foreach my $geneID (keys %{ $rhash->{$seqID} }) {
	    foreach my $transcriptID (keys %{ $rhash->{$seqID}->{$geneID} }) {
		my $cdsID = $rhash->{$seqID}->{$geneID}->{$transcriptID}->{'cdsId'};
		my $proteinID = $rhash->{$seqID}->{$geneID}->{$transcriptID}->{'proteinId'};
		$protID_cdsID->{$proteinID} = $cdsID;
	    }
	}
    }

}
    
sub get_btab_files{
    my ($directory,$file) = @_;
    my @files;
    if(-d $directory){
	opendir(DIR, $directory) or die "Unable to access $directory due to $!";
	while( my $filename = readdir(DIR)) {
	    if($filename =~ /(.+)\.btab$/) {
		push (@files, "$directory/$filename");
	    }
	}
    }
    if($file ne ""){
	push @files,$file;
    }
    return \@files;
}


sub check_parameters{
    my ($options) = @_;
    
    if(!$options{'bsml_dir'} or !$options{'output'} or (!$options{'btab_dir'} && !$options{'btab_file'})) {
	pod2usage({-exitval => 2,  -message => "$0: All the required options are not specified", -verbose => 1, -output => \*STDERR});    
    }

    if(!-d $options{'bsml_dir'}) {
	print STDERR "bsml repository directory \"$options{'bsml_dir'}\" cannot be found.  Aborting...\n";
	exit 5;
    }

    if((! -d $options{'btab_dir'}) && (! -e $options{'btab_file'})) {
	print STDERR "btab directory \"$options{'btab_dir'}\" or btab file $options{'btab_file'} cannot be found.  Aborting...\n";
	exit 5;
    } 

    return 1;
}

sub createAndAddBtabLine {

    my %args = @_;
    my $doc = $args{'doc'};

    #determine if the query name and the dbmatch name are a unique pair in the document

    my $alignment_pair_list = BSML::BsmlDoc::BsmlReturnAlignmentLookup( "$args{'query_name'}", "$args{'dbmatch_accession'}" );

    my $alignment_pair = '';
    if( $alignment_pair_list )
    {
	$alignment_pair = $alignment_pair_list->[0];
    }

    if( $alignment_pair  )
	  {
	    #add a new BsmlSeqPairRun to the alignment pair and return
	    my $seq_run = $alignment_pair->returnBsmlSeqPairRunR( $alignment_pair->addBsmlSeqPairRun() );

	    if( $args{'start_query'} > $args{'stop_query'} )
	      {
		$seq_run->setattr( 'refpos', $args{'stop_query'}-1 );
		$seq_run->setattr( 'runlength', $args{'start_query'} - $args{'stop_query'} + 1 );
		$seq_run->setattr( 'refcomplement', 1 );
	      }
	    else
	      {
		$seq_run->setattr( 'refpos', $args{'start_query'}-1 );
		$seq_run->setattr( 'runlength', $args{'stop_query'} - $args{'start_query'} + 1 );
		$seq_run->setattr( 'refcomplement', 0 );
	      }

	    #the database sequence is always 5' to 3'

	    $seq_run->setattr( 'comppos', $args{'start_hit'} -1 )                                   if (defined ($args{'start_hit'}));
	    $seq_run->setattr( 'comprunlength', $args{'stop_hit'} - $args{'start_hit'} + 1 )     if ((defined ($args{'start_hit'})) and (defined ($args{'stop_hit'})));
	    $seq_run->setattr( 'compcomplement', 0 );

	    $seq_run->setattr( 'runscore', $args{'bit_score'} )                                  if (defined ($args{'bit_score'}));
	    $seq_run->setattr( 'runprob', $args{'e_value'} )                                     if (defined ($args{'e_value'}));

	    $seq_run->addBsmlAttr( 'percent_identity', $args{'percent_identity'} )               if (defined ($args{'percent_identity'}));   
	    $seq_run->addBsmlAttr( 'percent_similarity', $args{'percent_similarity'} )           if (defined ($args{'percent_similarity'}));
	    $seq_run->addBsmlAttr( 'chain_number', $args{'chain_number'} )                       if (defined ($args{'chain_number'}));
	    $seq_run->addBsmlAttr( 'segment_number', $args{'segment_number'} )                   if (defined ($args{'segment_number'}));
	    $seq_run->addBsmlAttr( 'p_value', $args{'p_value'} )                                 if (defined ($args{'p_value'}));

	    return $alignment_pair;
	  }

    #no alignment pair matches, add a new alignment pair and sequence run
    #check to see if sequences exist in the BsmlDoc, if not add them with basic attributes
    my $seq;
    
    if( !( $doc->returnBsmlSequenceByIDR( "$args{'query_name'}")) ){
	    $seq = $doc->createAndAddSequence( "$args{'query_name'}", "$args{'query_name'}", $args{'query_length'}, 'aa', $args{'class'} );
    }
    
    if( !( $doc->returnBsmlSequenceByIDR( "$args{'dbmatch_accession'}")) ){
        $seq = $doc->createAndAddSequence( "$args{'dbmatch_accession'}", "$args{'dbmatch_header'}", '', 'aa', $args{'class'} );
    }

    ## see if the dbmatch_header format is recognized.  if so, add some cross-references
    if (defined $args{'dbmatch_header'}) {
        $doc->createAndAddCrossReferencesByParse( sequence => $seq, string => $args{'dbmatch_header'} );
    }

    $alignment_pair = $doc->returnBsmlSeqPairAlignmentR( $doc->addBsmlSeqPairAlignment() );
    

    $alignment_pair->setattr( 'refseq', "$args{'query_name'}" )                                 if (defined ($args{'query_name'}));
    $alignment_pair->setattr( 'compseq', "$args{'dbmatch_accession'}" )                         if (defined ($args{'dbmatch_accession'}));
    
    BSML::BsmlDoc::BsmlSetAlignmentLookup( "$args{'query_name'}", "$args{'dbmatch_accession'}", $alignment_pair );

    $alignment_pair->setattr( 'refxref', ':'.$args{'query_name'})        if (defined ($args{'query_name'}));                     
    $alignment_pair->setattr( 'refstart', 0 );
    $alignment_pair->setattr( 'refend', $args{'query_length'} )      if (defined ($args{'query_length'}));
    $alignment_pair->setattr( 'reflength', $args{'query_length'} )       if (defined ($args{'query_length'}));

    $alignment_pair->setattr( 'method', $args{'blast_program'} )         if (defined ($args{'blast_program'}));

    $alignment_pair->setattr( 'compxref', $args{'search_database'}.':'.$args{'dbmatch_accession'} )  if ((defined ($args{'search_database'})) and (defined ($args{'dbmatch_accession'})));

    my $seq_run = $alignment_pair->returnBsmlSeqPairRunR( $alignment_pair->addBsmlSeqPairRun() );

    if( $args{'start_query'} > $args{'stop_query'} )
      {
	$seq_run->setattr( 'refpos', $args{'stop_query'}-1 );
	$seq_run->setattr( 'runlength', $args{'start_query'} - $args{'stop_query'} + 1 );
	$seq_run->setattr( 'refcomplement', 1 );
      }
    else
      {
	$seq_run->setattr( 'refpos', $args{'start_query'} -1);
	$seq_run->setattr( 'runlength', $args{'stop_query'} - $args{'start_query'} + 1 );
	$seq_run->setattr( 'refcomplement', 0 );
      }

    #the database sequence is always 5' to 3'
    
    $seq_run->setattr( 'comppos', $args{'start_hit'} -1)                                     if (defined  ($args{'start_hit'}));
    $seq_run->setattr( 'comprunlength', ($args{'stop_hit'} - $args{'start_hit'} + 1))      if ((defined ($args{'start_hit'})) and (defined ($args{'stop_hit'})));
    $seq_run->setattr( 'compcomplement', 0 );
    
    $seq_run->setattr( 'runscore', $args{'bit_score'} )                                    if (defined  ($args{'bit_score'}));
    $seq_run->setattr( 'runprob', $args{'e_value'} )                                       if (defined  ($args{'e_value'}));

    $seq_run->addBsmlAttr( 'percent_identity', $args{'percent_identity'} )                 if (defined  ($args{'percent_identity'}));
    $seq_run->addBsmlAttr( 'percent_similarity', $args{'percent_similarity'} )             if (defined  ($args{'percent_similarity'}));
    $seq_run->addBsmlAttr( 'chain_number', $args{'chain_number'} )                         if (defined  ($args{'chain_number'}));
    $seq_run->addBsmlAttr( 'segment_number', $args{'segment_number'} )                     if (defined  ($args{'segment_number'}));
    $seq_run->addBsmlAttr( 'p_value', $args{'p_value'} )                                   if (defined  ($args{'p_value'}));

    return $alignment_pair;


}
