#!/usr/bin/perl

BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

blastbtab2bsml.pl  - convert info stored in btab files into BSML documents

=head1 SYNOPSIS

USAGE:  blastbtab2bsml.pl -b btab_dir -o blastp.bsml

=head1 OPTIONS

=over 4

=item *

B<--output,-o> [REQUIRED] output BSML file containing bit_score information

=item * 

B<--btab_dir,-b> [REQUIRED] Dir containing btab files

=item *

B<--max_hsp_count,-m> Maximum number of HSPs stored per alignment
    
=item *

B<--class,-c> The ref/comp sequence type.  Default is 'assembly'

=item *

B<--gzip,-z> gzip bsml output

=item *

B<--help,-h> This help message

=back

=head1   DESCRIPTION

blastbtab2bsml.pl is designed to convert information in btab files into BSML documents.

NOTE:  

Calling the script name with NO flags/options or --help will display the syntax requirement.

=cut


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use English;
use File::Basename;
use File::Path;
use Pod::Usage;
use Ergatis::Logger;
use BSML::BsmlRepository;
use BSML::BsmlBuilder;
use BSML::BsmlParserTwig;

my %options = ();
my $results = GetOptions (\%options, 
              'analysis_id|a=s',
              'btab_dir|b=s', 
              'btab_file|f=s',
              'query_file_path|q=s',
              'output|o=s', 
              'max_hsp_count|m=s',
              'pvalue|p=s', 
              'log|l=s',
              'debug=s',
              'class|c=s',
              'gzip|z=s',
              'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

if($options{'help'}) {
    pod2usage();
}

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

my $seq_data_import_file = $options{'query_file_path'};

$seq_data_import_file =~ s/\.gz$|\.gzip$//;

my $deflines = get_deflines($options{'query_file_path'});

#generate lookups
$doc->makeCurrentDocument();
parse_blast_btabs($files);

## create sequence stubs for input sequences with null output
foreach my $seq_id(keys(%{$deflines})) {
    if(!($doc->returnBsmlSequenceByIDR($seq_id))){
        my $seq = $doc->createAndAddSequence( 
                                                $seq_id, 
                                                $seq_id, 
                                                '',     # length
                                                '',   # molecule 
                                                '',     # class
                                            );
        
        $seq->addBsmlLink('analysis', '#' . $options{analysis_id}, 'input_of');
        $doc->createAndAddBsmlAttribute( $seq, 'defline', $deflines->{$seq_id} );
        $doc->createAndAddSeqDataImport( $seq, 'fasta', $seq_data_import_file, '', $seq_id);
    }
}

## add the analysis element
my $program = $options{analysis_id};
   $program =~ s/_analysis//;
$doc->createAndAddAnalysis(
                            id => $options{analysis_id},
                            sourcename => $options{'output'},
                            algorithm => $program,
                            program => $program,
                          );

if($options{'gzip'}){
    $doc->write($options{'output'},undef,1);
}
else{
    $doc->write($options{'output'});
}

sub parse_blast_btabs {

    my $btab_files = shift;

    my $num;
    foreach my $file(@$btab_files) {
    $num++;
    open (BTAB, "$file") or $logger->logdie("Unable to open \"$file\" due to $!");
    $logger->debug("opening $file $num using pvalue p_value $options{'pvalue'}, hsp cutoff $options{'max_hsp_count'}") if($logger->is_debug());
    my $hsplookup;
    while(my $line = <BTAB>) {
        chomp($line);
        my @btab = split("\t", $line);
        if(($btab[20] <= $options{'pvalue'}) && ($btab[0] ne "") && ($btab[5] ne "")){
        if(!(exists $hsplookup->{$btab[0]}->{$btab[5]})){
            $hsplookup->{$btab[0]}->{$btab[5]} = [];
        }
        push @{$hsplookup->{$btab[0]}->{$btab[5]}},{'pvalue'=>$btab[20],
                                'line'=>$line};
        }
        else{
        $logger->debug("Skipping $btab[0] $btab[5] pvalue $btab[20] above pvalue cutoff of $options{'pvalue'}") if($logger->is_debug());
        }
    }
    close BTAB;
    
    ## This block prepares an array of HSP data structured so that we can use
    ## existing subs taken from Sybil::Util to calculate coverage/identity info
    my $cov_qual_stats = {};
    my @hsp_ref_array = ();
    my $new_subject = '';
    foreach my $query (keys %$hsplookup) {
        foreach my $subject (keys %{$hsplookup->{$query}}){
            my @hsps = sort {$a->{'pvalue'} <=> $b->{'pvalue'}} @{$hsplookup->{$query}->{$subject}};
            my $maxhsp;
            if($options{'max_hsp_count'} ne ""){
                $maxhsp = ($options{'max_hsp_count'}<scalar(@hsps)) ? $options{'max_hsp_count'} : scalar(@hsps);
            } else {
                $maxhsp = scalar(@hsps);
            }
            my $queryid;
            for(my $i=0;$i<$maxhsp;$i++){
                my $line = $hsps[$i]->{'line'};
                my @btab = split("\t", $line);
                $queryid = $btab[0] if($btab[0] && (!$queryid));

                for (my $i=0;$i<scalar(@btab);$i++){
                    if ($btab[$i] eq 'N/A'){
                        $btab[$i] = undef;
                    }
                }
                
                my $orig_dbmatch_accession = $btab[5];
                $btab[5] =~ s/[^a-z0-9\_\.\-]/_/gi;

                $new_subject = $btab[5];
                
                my $qfmin = $btab[6];
                my $qfmax = $btab[7];
                my $qstrand = 0;
                my $tfmin = $btab[8];
                my $tfmax = $btab[9];
                my $tstrand = 0;
                    
                ## if query positions are on the reverse strand
                if ($btab[6] > $btab[7]) {
                        $qfmin = $btab[7];
                        $qfmax = $btab[6];
                        $qstrand = 1;
                }
                
                ## if target positions are on the reverse strand
                if ($btab[8] > $btab[9]) {
                        $tfmin = $btab[9];
                        $tfmax = $btab[8];
                        $tstrand = 1;
                }
                
                ## transform the start positions to interbase 
                $qfmin = $qfmin - 1;
                $tfmin = $tfmin - 1;    
                    
                my $hsp_ref = { 
                                'query_protein_id'  => $btab[0],
                                'target_protein_id' => $btab[5],
                                'significance'      => $btab[20],
                                'percent_identity'  => $btab[10],
                                'query_seqlen'      => $btab[2],
                                'target_seqlen'     => $btab[18],
                                'query_fmin'        => $qfmin,
                                'query_fmax'        => $qfmax,
                                'query_strand'      => $btab[17],
                                'target_fmin'       => $tfmin,
                                'target_fmax'       => $tfmax,
                                'target_strand'     => $tstrand, ## target strand is not captured in btab 
                              };
                push (@hsp_ref_array, $hsp_ref);
            }
            if (!defined($cov_qual_stats->{$query})) {
               $cov_qual_stats->{$query} = {};
            }
            
            my $coverage_arr_ref = &getAvgBlastPPctCoverage(\@hsp_ref_array);
        
            $cov_qual_stats->{$query}->{$new_subject} = {
                'percent_coverage_refseq'   =>  sprintf("%.1f",$coverage_arr_ref->[0]),
                'percent_coverage_compseq'      =>  sprintf("%.1f",$coverage_arr_ref->[1]),
                                                        };
        }
    }

    foreach my $query (keys %$hsplookup){
        foreach my $subject (keys %{$hsplookup->{$query}}){
        my @hsps = sort {$a->{'pvalue'} <=> $b->{'pvalue'}} @{$hsplookup->{$query}->{$subject}};
        my $maxhsp;
        if($options{'max_hsp_count'} ne ""){
            $maxhsp = ($options{'max_hsp_count'}<scalar(@hsps)) ? $options{'max_hsp_count'} : scalar(@hsps);
        }
        else{
            $maxhsp = scalar(@hsps);
        }
        my $queryid;
        for(my $i=0;$i<$maxhsp;$i++){
            my $line = $hsps[$i]->{'line'};
            my @btab = split("\t", $line);
            $logger->debug("Storing HSP $btab[0] $btab[5] $btab[20]") if($logger->is_debug());
            $queryid = $btab[0] if($btab[0] && (!$queryid));

            for (my $i=0;$i<scalar(@btab);$i++){
            if ($btab[$i] eq 'N/A'){
                $btab[$i] = undef;
            }
            }
            
            ## dbmatch_accession needs to be alphanumeric or _-.
            ##  but the original needs to be passed to createAndAddBtabLine so it can
            ##  be recognized and parsed
            my $orig_dbmatch_accession = $btab[5];
            $btab[5] =~ s/[^a-z0-9\_\.\-]/_/gi;

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
                              raw_score          => $btab[12],
                              bit_score          => $btab[13],
                              dbmatch_header     => $btab[15],
                              blast_frame        => $btab[16],
                              qry_strand         => $btab[17],
                              hit_length         => $btab[18],
                              e_value            => $btab[19],
                              p_value            => $btab[20],
                              percent_coverage_refseq    => $cov_qual_stats->{$btab[0]}->{$btab[5]}->{'percent_coverage_refseq'},
                              percent_coverage_compseq       => $cov_qual_stats->{$btab[0]}->{$btab[5]}->{'percent_coverage_compseq'},
                              class              => $class,
                              orig_dbmatch_accession => $orig_dbmatch_accession
                              );

            my $seq = $doc->returnBsmlSequenceByIDR($btab[5]);



        }
        if($queryid){
            my $seq = $doc->returnBsmlSequenceByIDR($queryid);
        }
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
    
    if(!$options{'output'} or (!$options{'btab_dir'} && !$options{'btab_file'})) {
    pod2usage({-exitval => 2,  -message => "$0: All the required options are not specified", -verbose => 1, -output => \*STDERR});    
    }

    if((! -d $options{'btab_dir'}) && (! -e $options{'btab_file'})) {
    print STDERR "btab directory \"$options{'btab_dir'}\" or btab file $options{'btab_file'} cannot be found.  Aborting...\n";
    exit 5;
    } 

    ## handle some defaults
    if (! $options{analysis_id}) {
        $options{analysis_id} = 'unknown_analysis';
    }

    return 1;
}

sub createAndAddBtabLine {

    my %args = @_;
    my $doc = $args{'doc'};

#    my $orig_dbmatch_accession = $args{dbmatch_accession};
#    $args{dbmatch_accession} =~ s/[^a-z0-9\_]/_/gi;

    #determine if the query name and the dbmatch name are a unique pair in the document
    my $alignment_pair_list = BSML::BsmlDoc::BsmlReturnAlignmentLookup( $args{'query_name'}, $args{'dbmatch_accession'} );

    my $alignment_pair = '';
    if( $alignment_pair_list ) {
        $alignment_pair = $alignment_pair_list->[0];
    }

    if( $alignment_pair ) {
        #add a new BsmlSeqPairRun to the alignment pair and return
        my $seq_run = $alignment_pair->returnBsmlSeqPairRunR( $alignment_pair->addBsmlSeqPairRun() );

        if( $args{'start_query'} > $args{'stop_query'} ) {
            $seq_run->setattr( 'refpos', $args{'stop_query'}-1 );
            $seq_run->setattr( 'runlength', $args{'start_query'} - $args{'stop_query'} + 1 );
            $seq_run->setattr( 'refcomplement', 1 );
        } else {
            $seq_run->setattr( 'refpos', $args{'start_query'}-1 );
            $seq_run->setattr( 'runlength', $args{'stop_query'} - $args{'start_query'} + 1 );
            $seq_run->setattr( 'refcomplement', 0 );
        }

        #the database sequence is always 5' to 3'
        $seq_run->setattr( 'comppos', $args{'start_hit'} -1 )                                if (defined ($args{'start_hit'}));
        $seq_run->setattr( 'comprunlength', $args{'stop_hit'} - $args{'start_hit'} + 1 )     if ((defined ($args{'start_hit'})) and (defined ($args{'stop_hit'})));
        $seq_run->setattr( 'compcomplement', 0 );
        $seq_run->setattr( 'runscore', $args{'bit_score'} )                                  if (defined ($args{'bit_score'}));
        $seq_run->setattr( 'runprob', $args{'e_value'} )                                     if (defined ($args{'e_value'}));
        $seq_run->addBsmlAttr( 'class', 'match_part' );
        $seq_run->addBsmlAttr( 'percent_identity', $args{'percent_identity'} )               if (defined ($args{'percent_identity'}));   
        $seq_run->addBsmlAttr( 'percent_similarity', $args{'percent_similarity'} )           if (defined ($args{'percent_similarity'}));
        $seq_run->addBsmlAttr( 'percent_coverage_refseq', $args{'percent_coverage_refseq'} )                 if (defined ($args{'percent_coverage_refseq'}));   
        $seq_run->addBsmlAttr( 'percent_coverage_compseq', $args{'percent_coverage_compseq'} )                       if (defined ($args{'percent_coverage_compseq'}));   
        $seq_run->addBsmlAttr( 'p_value', $args{'p_value'} )                                 if (defined ($args{'p_value'}));
        $seq_run->addBsmlAttr( 'p_value', $args{'p_value'} )                                 if (defined ($args{'p_value'}));


        return $alignment_pair;
    }

    #no alignment pair matches, add a new alignment pair and sequence run
    #check to see if sequences exist in the BsmlDoc, if not add them with basic attributes
    my $seq;
    
    if( !( $doc->returnBsmlSequenceByIDR( "$args{'query_name'}")) ){
        $seq = $doc->createAndAddSequence( "$args{'query_name'}", "$args{'query_name'}", $args{'query_length'}, '', $args{'class'} );
        $seq->addBsmlLink('analysis', '#' . $options{analysis_id}, 'input_of');
        
        my $defline = $deflines->{$args{'query_name'}};
        $doc->createAndAddBsmlAttribute( $seq, 'defline', $defline );
        
        if ($options{'query_file_path'}) {
            $doc->createAndAddSeqDataImport( $seq, 'fasta', $seq_data_import_file, '', "$args{'query_name'}");
        }
        
    }
    
    if( !( $doc->returnBsmlSequenceByIDR( "$args{'dbmatch_accession'}")) ){
        $seq = $doc->createAndAddSequence( "$args{'dbmatch_accession'}", "$args{'dbmatch_header'}", ($args{'hit_length'} || 0), '', $args{'class'} );
        $doc->createAndAddBsmlAttribute( $seq, 'defline', "$args{orig_dbmatch_accession} $args{dbmatch_header}" );

        ## see if the dbmatch_header format is recognized.  if so, add some cross-references
        if (defined $args{'dbmatch_header'}) {
            $doc->createAndAddCrossReferencesByParse( sequence => $seq, string => $args{orig_dbmatch_accession} );
        }
    }
    
    $alignment_pair = $doc->returnBsmlSeqPairAlignmentR( $doc->addBsmlSeqPairAlignment() );
    
    ## to the alignment pair, add a Link to the analysis
    $alignment_pair->addBsmlLink('analysis', '#' . $options{analysis_id}, 'computed_by');

    $alignment_pair->setattr( 'refseq', "$args{'query_name'}" )                                 if (defined ($args{'query_name'}));
    $alignment_pair->setattr( 'compseq', "$args{'dbmatch_accession'}" )                         if (defined ($args{'dbmatch_accession'}));
    
    BSML::BsmlDoc::BsmlSetAlignmentLookup( "$args{'query_name'}", "$args{'dbmatch_accession'}", $alignment_pair );

    $alignment_pair->setattr( 'refxref', ':'.$args{'query_name'})        if (defined ($args{'query_name'}));                     
    $alignment_pair->setattr( 'refstart', 0 );
    $alignment_pair->setattr( 'refend', $args{'query_length'} )          if (defined ($args{'query_length'}));
    $alignment_pair->setattr( 'reflength', $args{'query_length'} )       if (defined ($args{'query_length'}));
    $alignment_pair->setattr( 'method', $args{'blast_program'} )         if (defined ($args{'blast_program'}));
    $alignment_pair->setattr( 'class', 'match' );
    $alignment_pair->setattr( 'compxref', $args{'search_database'}.':'.$args{'dbmatch_accession'} )  if ((defined ($args{'search_database'})) and (defined ($args{'dbmatch_accession'})));

    my $seq_run = $alignment_pair->returnBsmlSeqPairRunR( $alignment_pair->addBsmlSeqPairRun() );

    if( $args{'start_query'} > $args{'stop_query'} ) {
        $seq_run->setattr( 'refpos', $args{'stop_query'}-1 );
        $seq_run->setattr( 'runlength', $args{'start_query'} - $args{'stop_query'} + 1 );
        $seq_run->setattr( 'refcomplement', 1 );
    } else {
        $seq_run->setattr( 'refpos', $args{'start_query'} -1);
        $seq_run->setattr( 'runlength', $args{'stop_query'} - $args{'start_query'} + 1 );
        $seq_run->setattr( 'refcomplement', 0 );
    }

    #the database sequence is always 5' to 3'
    $seq_run->setattr( 'comppos', $args{'start_hit'} -1)                                   if (defined  ($args{'start_hit'}));
    $seq_run->setattr( 'comprunlength', ($args{'stop_hit'} - $args{'start_hit'} + 1))      if ((defined ($args{'start_hit'})) and (defined ($args{'stop_hit'})));
    $seq_run->setattr( 'compcomplement', 0 );
    $seq_run->setattr( 'runscore', $args{'bit_score'} )                                    if (defined  ($args{'bit_score'}));
    $seq_run->setattr( 'runprob', $args{'e_value'} )                                       if (defined  ($args{'e_value'}));
    $seq_run->addBsmlAttr( 'class', 'match_part' );
    $seq_run->addBsmlAttr( 'percent_identity', $args{'percent_identity'} )                 if (defined  ($args{'percent_identity'}));
    $seq_run->addBsmlAttr( 'percent_similarity', $args{'percent_similarity'} )             if (defined  ($args{'percent_similarity'}));
    $seq_run->addBsmlAttr( 'percent_coverage_refseq', $args{'percent_coverage_refseq'} )                   if (defined ($args{'percent_coverage_refseq'}));   
    $seq_run->addBsmlAttr( 'percent_coverage_compseq', $args{'percent_coverage_compseq'} )                         if (defined ($args{'percent_coverage_compseq'}));   
    $seq_run->addBsmlAttr( 'chain_number', $args{'chain_number'} )                         if (defined  ($args{'chain_number'}));
    $seq_run->addBsmlAttr( 'segment_number', $args{'segment_number'} )                     if (defined  ($args{'segment_number'}));
    $seq_run->addBsmlAttr( 'p_value', $args{'p_value'} )                                   if (defined  ($args{'p_value'}));

    return $alignment_pair;
}


## Returns an array reference where 
##     [0] = query percent coverage, 
## and [1] = target percent coverage
sub getAvgBlastPPctCoverage {
    my($hsps) = @_;
    my $qsum = 0;
    my $tsum = 0;
    my $numHsps = 0;

    # Group by query and target id
    my $hspsByQuery = &groupByMulti($hsps, ['query_protein_id', 'target_protein_id']);

    foreach my $queryId (keys %$hspsByQuery) {
        my $hspsByTarget = $hspsByQuery->{$queryId};

        foreach my $subjId (keys %$hspsByTarget) {
            ++$numHsps;
            my $shsps = $hspsByTarget->{$subjId};
            my $querySeqLen = $shsps->[0]->{'query_seqlen'};
            my $targetSeqLen = $shsps->[0]->{'target_seqlen'};

            my @queryIntervals = map { {'fmin' => $_->{'query_fmin'}, 'fmax' => $_->{'query_fmax'}, 'strand' => $_->{'query_strand'}} } @$shsps;
            my @targetIntervals = map { {'fmin' => $_->{'target_fmin'}, 'fmax' => $_->{'target_fmax'}, 'strand' => $_->{'target_strand'}} } @$shsps;

            my $mergedQueryIntervals = &mergeOverlappingIntervals(\@queryIntervals);
            my $mergedTargetIntervals = &mergeOverlappingIntervals(\@targetIntervals);

            my $queryHitLen = 0;
            my $targetHitLen = 0;

            map { $queryHitLen += ($_->{'fmax'} - $_->{'fmin'}); } @$mergedQueryIntervals;
            map { $targetHitLen += ($_->{'fmax'} - $_->{'fmin'}); } @$mergedTargetIntervals;

            $qsum += $queryHitLen / $querySeqLen;
            $tsum += $targetHitLen / $targetSeqLen;
        }
    }

    if ($numHsps == 0) {
        return undef;
    } else {
        return [($qsum/$numHsps*100.0), ($tsum/$numHsps*100.0)];
    }
    #return ($numHsps > 0) ? ($sum/($numHsps * 2) * 100.0) : undef;
}

# Generalized version of groupBy 
sub groupByMulti {
    my($arrayref, $keyFields) = @_;
    my $nKeys = scalar(@$keyFields);
    my $groups = {};

    foreach my $a (@$arrayref) {
        my @keyValues = map { $a->{$_} } @$keyFields;
        my $hash = $groups;

        for (my $i = 0;$i < $nKeys;++$i) {
            my $kv = $keyValues[$i];

            if ($i < ($nKeys-1)) {
                $hash->{$kv} = {} if (!defined($hash->{$kv}));
                $hash = $hash->{$kv};
            } else {
                $hash->{$kv} = [] if (!defined($hash->{$kv}));
                push(@{$hash->{$kv}}, $a);
            }
        }
    }
    return $groups;
}

# Generate a new set of intervals by merging any that overlap in the original set.
#
sub mergeOverlappingIntervals {
    my($intervals) = @_;

    # result set of intervals
    my $merged = [];

    # sort all intervals by fmin
    my @sorted = sort { $a->{'fmin'} <=> $b->{'fmin'} } @$intervals;
    
    # current interval
    my $current = undef;

    foreach my $i (@sorted) {
        if (!defined($current)) {
            # case 1: no current interval
            $current = $i;
        } else {
            # case 2: compare current interval to interval $i
            if ($i->{'fmin'} > $current->{'fmax'}) {   
                # case 2a: no overlap
                push(@$merged, $current);
                $current = $i;
            } elsif ($i->{'fmax'} > $current->{'fmax'}) {
                # case 2b: overlap, with $i ending to the right of $current
                $current->{'fmax'} = $i->{'fmax'};
            }
        }
    }
    push(@$merged, $current) if (defined($current));

    return $merged;
}


## retrieve deflines from a fasta file
sub get_deflines {
    my ($fasta_file) = @_;

    my $deflines = {};

    my $ifh; 
   
    if (! -e $fasta_file) {
        if (-e $fasta_file.".gz") {
            $fasta_file .= ".gz";
        } elsif (-e $fasta_file.".gzip") {
            $fasta_file .= ".gzip";
        }
    }
    
    if ($fasta_file =~ /\.(gz|gzip)$/) {
        open ($ifh, "<:gzip", $fasta_file)
          || $logger->logdie("can't open input file '$fasta_file': $!");
    } else {
        open ($ifh, $fasta_file)
          || $logger->logdie("Failed opening '$fasta_file' for reading: $!");
      }

    while (<$ifh>) {
        unless (/^>/) {
            next;
        }
        chomp;
        if (/^>((\S+).*)$/) {
            $deflines->{$2} = $1;
        } 
    }
    close $ifh;
    
    if (scalar(keys(%{$deflines})) < 1) {
        $logger->warn("defline lookup failed for '$fasta_file'");
    }

    return $deflines;
}
