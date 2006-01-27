#!/usr/local/packages/perl-5.8.5/bin/perl

eval 'exec /usr/local/packages/perl-5.8.5/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

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
BEGIN {
use Workflow::Logger;
use BSML::BsmlRepository;
use BSML::BsmlBuilder;
use BSML::BsmlParserTwig;
}

my %options = ();
my $results = GetOptions (\%options, 
              'analysis_id|a=s',
              'btab_dir|b=s', 
              'btab_file|f=s',
              'bsml_dir|d=s', ## deprecated.  keeping for backward compat (for now)
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
$doc->makeCurrentDocument();
parse_blast_btabs($files);

## add the analysis element
$doc->createAndAddAnalysis(
                            id => $options{analysis_id},
                            sourcename => $options{'output'},
                          );

$doc->write($options{'output'});

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
        if(($btab[20] < $options{'pvalue'}) && ($btab[0] ne "") && ($btab[5] ne "")){
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
        $seq_run->addBsmlAttr( 'percent_identity', $args{'percent_identity'} )               if (defined ($args{'percent_identity'}));   
        $seq_run->addBsmlAttr( 'percent_similarity', $args{'percent_similarity'} )           if (defined ($args{'percent_similarity'}));
        $seq_run->addBsmlAttr( 'p_value', $args{'p_value'} )                                 if (defined ($args{'p_value'}));

        return $alignment_pair;
    }

    #no alignment pair matches, add a new alignment pair and sequence run
    #check to see if sequences exist in the BsmlDoc, if not add them with basic attributes
    my $seq;
    
    if( !( $doc->returnBsmlSequenceByIDR( "$args{'query_name'}")) ){
        $seq = $doc->createAndAddSequence( "$args{'query_name'}", "$args{'query_name'}", $args{'query_length'}, 'aa', $args{'class'} );
        $seq->addBsmlLink('analysis', '#' . $options{analysis_id}, 'input_of');
    }
    
    if( !( $doc->returnBsmlSequenceByIDR( "$args{'dbmatch_accession'}")) ){
        $seq = $doc->createAndAddSequence( "$args{'dbmatch_accession'}", "$args{'dbmatch_header'}", ($args{'hit_length'} || 0), 'aa', $args{'class'} );
        $seq->addBsmlLink('analysis', '#' . $options{analysis_id}, 'input_of');
    }

    ## see if the dbmatch_header format is recognized.  if so, add some cross-references
    if (defined $args{'dbmatch_header'}) {
        $doc->createAndAddCrossReferencesByParse( sequence => $seq, string => $args{orig_dbmatch_accession} );
    }
    
    $alignment_pair = $doc->returnBsmlSeqPairAlignmentR( $doc->addBsmlSeqPairAlignment() );
    
    ## to the alignment pair, add a Link to the analysis
    $alignment_pair->addBsmlLink('analysis', '#' . $options{analysis_id}, 'computed_by');

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
    $seq_run->addBsmlAttr( 'percent_identity', $args{'percent_identity'} )                 if (defined  ($args{'percent_identity'}));
    $seq_run->addBsmlAttr( 'percent_similarity', $args{'percent_similarity'} )             if (defined  ($args{'percent_similarity'}));
    $seq_run->addBsmlAttr( 'chain_number', $args{'chain_number'} )                         if (defined  ($args{'chain_number'}));
    $seq_run->addBsmlAttr( 'segment_number', $args{'segment_number'} )                     if (defined  ($args{'segment_number'}));
    $seq_run->addBsmlAttr( 'p_value', $args{'p_value'} )                                   if (defined  ($args{'p_value'}));

    return $alignment_pair;
}







