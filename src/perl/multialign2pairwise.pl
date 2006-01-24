#!/usr/local/bin/perl

=head1  NAME 

dummy.pl - do nothing

=head1 SYNOPSIS

USAGE:  dummy.pl --debug debug_level --log log_file

=head1 OPTIONS

=item *

B<--debug,-d> Debug level.  Use a large number to turn on verbose debugging. 

=item *

B<--log,-l> Log file

=item *

B<--help,-h> This help message

=back

=head1   DESCRIPTION

=cut


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

BEGIN {
use Workflow::Logger;
use BSML::BsmlReader;
use BSML::BsmlBuilder;
use BSML::BsmlParserSerialSearch;
}

my %options = ();
my $results = GetOptions (\%options, 
			  'log|l=s',
			  'debug=s',
			  'output|o=s',
			  'directory|d=s',
			  'file|f=s',
			  'filelist|l=s',
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

my @filelist;

if($options{'file'}){
    push @filelist,$options{'file'};
}
if($options{'filelist'}){
    &get_list_from_file(\@filelist,$options{'filelist'});
}
if($options{'directory'}){
    &get_list_from_directory(\@filelist,$options{'directory'});
}


my $clusterlookup = {};
my $count;

my $multiAlnParser = new BSML::BsmlParserSerialSearch( MultipleAlignmentCallBack => \&multipleAlignmentHandler );

my $doc = new BSML::BsmlBuilder();

foreach my $file (@filelist){
    $multiAlnParser->parse( $file );
}

foreach my $alnname (keys %$clusterlookup){
    my @names = keys %{$clusterlookup->{$alnname}};
    for(my $i=0;$i<scalar(@names);$i++){
	for(my $j=0;$j<scalar(@names);$j++){
	    my $refseq = $names[$j];
	    my $compseq = $names[$i];
	    if($refseq ne $compseq){
		my %args = ('refseq' => $refseq,
			    'compseq' => $compseq,
			    'runscore' => 100,
			    'runprob' => 0,
			    'doc' => $doc);
		&generateSeqPairs(%args);
	    }
	}
    }
}


$doc->write($options{'output'});


sub multipleAlignmentHandler
{
    my $aln = shift;
    my $bsml_reader = new BSML::BsmlReader();

    my $maln_ref = $bsml_reader->readMultipleAlignmentTable($aln);
    
    foreach my $alnSum ( @{ $maln_ref->{'AlignmentSummaries'} } )
    {
	$count++;
	foreach my $alnSeq ( @{ $alnSum->{'AlignedSequences'} } )
	{
	    my $name = $alnSeq->{'name'};
	    $name =~ s/:[\d]*//;
	    $clusterlookup->{$count}->{$name} = 1;
	}
    }

    return $clusterlookup;
}

sub generateSeqPairs {

    my %args = @_;
    my $doc = $args{'doc'};

    #determine if the query name and the dbmatch name are a unique pair in the document

    my $alignment_pair_list = BSML::BsmlDoc::BsmlReturnAlignmentLookup( "$args{'query_name'}", "$args{'compseq'}" );

    my $alignment_pair = '';
    if( $alignment_pair_list )
    {
	$alignment_pair = $alignment_pair_list->[0];
    }

    if( $alignment_pair  ){
	    #add a new BsmlSeqPairRun to the alignment pair and return
	    my $seq_run = $alignment_pair->returnBsmlSeqPairRunR( $alignment_pair->addBsmlSeqPairRun() );

	    $seq_run->setattr( 'runscore', $args{'runscore'} );
	    $seq_run->setattr( 'runprob', $args{'runprop'});
	    return $alignment_pair;
	}

    #no alignment pair matches, add a new alignment pair and sequence run

    #check to see if sequences exist in the BsmlDoc, if not add them with basic attributes

    if( !( $doc->returnBsmlSequenceByIDR( "$args{'refseq'}")) ){
	$doc->createAndAddSequence( "$args{'refseq'}", "$args{'refseq'}", '', 'aa' );
    }
    
    if( !( $doc->returnBsmlSequenceByIDR( "$args{'compseq'}")) ){
	$doc->createAndAddSequence( "$args{'compseq'}", "$args{'compseq'}", '', 'aa' );
    }
    
    $alignment_pair = $doc->returnBsmlSeqPairAlignmentR( $doc->addBsmlSeqPairAlignment() );
    

    $alignment_pair->setattr( 'refseq', "$args{'refseq'}" )                                 if (defined ($args{'refseq'}));
    $alignment_pair->setattr( 'compseq', "$args{'compseq'}" )                         if (defined ($args{'compseq'}));

    BSML::BsmlDoc::BsmlSetAlignmentLookup( "$args{'refseq'}", "$args{'compseq'}", $alignment_pair );

    my $seq_run = $alignment_pair->returnBsmlSeqPairRunR( $alignment_pair->addBsmlSeqPairRun() );

    $seq_run->setattr( 'runscore', $args{'runscore'} );
    $seq_run->setattr( 'runprob', $args{'runprob'} );

    return $alignment_pair;

}


sub get_list_from_file{
    my ($filelist,$f) = @_;
    my @elts;
    my @lines;
    my @files = split(',',$f);
    foreach my $file (@files){
	if( $file){
	    open( FH, $file ) or die "Could not open $file";
	    while( my $line = <FH> ){
		chomp($line);
		push @lines,  split(',',$line) if($line =~ /\S+/);
	    }
	    foreach my $line (@lines){
		if($line){
		    my $filename = "$line";
		    push @$filelist,$filename;
		}
	    }
	    close( FH );
	}
    }
}

sub get_list_from_directory{
    my ($filelist, $dir, $glob) = @_;

    my @directories = split(',',$dir);
    foreach my $directory (@directories){
	opendir DIR, "$directory" or $logger->logdie("Can't read directory $directory");
	my @files = grep /\.bsml$/, readdir DIR;
	fisher_yates_shuffle( \@files );    # permutes @array in place
	foreach my $file (@files ){
	    my $filename = "$directory/$file";
	    push @$filelist,$filename;
	}
    }
}



sub check_parameters{
    my ($options) = @_;
    
    if(0){
	pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});    
    }
}


