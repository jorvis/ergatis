#!/usr/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

use strict;
use BSML::BsmlParserSerialSearch;
use Ergatis::Logger;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %options = ();
my $results = GetOptions( \%options, 'cogFile|c=s', 'bsmlModelList|m=s', 'outputDir|o=s', 'outputToken=s', 'maxCogSeqCount|s=s', 'extension|e=s', 'use_feature_ids_in_fasta=i', 'log|l=s', 'debug=s');

my $cogFile = $options{'cogFile'};
my $bsmlModelList = $options{'bsmlModelList'};
my $outDir = $options{'outputDir'};
my $maxCogSeqCount = $options{'maxCogSeqCount'};
my $outputToken = $options{'outputToken'};
$options{'outputToken'} .= "_$$";
if($options{'extension'} eq ''){
    $options{'extension'} = 'fsa';
}
my $use_feature_ids_in_fasta = defined($options{'use_feature_ids_in_fasta'}) && ($options{'use_feature_ids_in_fasta'} > 0);

## logging
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = Ergatis::Logger::get_logger();

# mapping from Seq.id and/or Seq-data-import.identifier to polypeptide sequence
my $Prot = {};
# mapping from polypeptide Feature.id to Seq.id
my $Feat = {};


my %seqParserArgs = ( ReadFeatureTables => 0, SequenceCallBack => \&createPolypeptideLookup );
if ($use_feature_ids_in_fasta) {
    $seqParserArgs{FeatureCallBack} = \&createFeatureLookup;
    $seqParserArgs{ReadFeatureTables} = 1;
}

# set up a serial parser to parse sequence elements, bypassing feature tables for efficiency if possible
my $seqParser = new BSML::BsmlParserSerialSearch(%seqParserArgs);

#Get rid of trailing slashes in directory names

$bsmlModelList =~ s/\/+$//;
$outDir =~ s/\/+$//; 

if( !$bsmlModelList )
{
    $logger->logdie("no Bsml Directory specified");
}
else
{
    if( ! -e $bsmlModelList )
    {
        $logger->logdie("could not open list file: $bsmlModelList");
    }
}

if( !$cogFile )
{
    $logger->logdie("cog file not specified");
}

if( !$outDir )
{
    $logger->logdie("output directory not specified");
}
else
{
    if( ! -d $outDir )
    {
        mkdir( $outDir );
    }
}

open FILE, $bsmlModelList or $logger->logdie("Can't open file $bsmlModelList");
while(my $bsmlFile=<FILE>)
{
    chomp $bsmlFile;
    if(-e $bsmlFile){
        $logger->debug("parsing $bsmlFile");
        $seqParser->parse( $bsmlFile );
    }
    else {
        $logger->logdie("could not open $bsmlFile.");
    }
}

open( INPUTCOGS, "<$cogFile" ) or $logger->logdie("could not open $cogFile.");

my $cog = undef;
my $list = [];

while( my $line = <INPUTCOGS> )
{
    if( $line =~ /^\t([\S]*)/ )
    {
        # A new sequence has been found and added to the current cog.
        $logger->debug("read line $line");
        push( @{$list}, $1 );
    }

    if( $line =~ /^COG\s+=\s+([^,\s]+)/ )
    {
        # A new cog has been encountered, flush the previous to disk if present
        my $newcog = $1;
        $logger->debug("read cog $line $newcog");
        &outputCog($cog, $list) if (defined($cog));
        $cog = $newcog;
        $list = [];
    }
}
outputCog($cog,$list) if (defined($cog));
exit(0);

## subroutines

sub outputCog {
    my($cog, $list) = @_;
    $logger->debug("outputting cog $cog");
    if(scalar(@{$list})>1){
		if(@{$list} <= $maxCogSeqCount){
		    open( OUTFILE, ">$outDir/$cog.$options{'outputToken'}.$options{'extension'}" ) or $logger->logdie("could not open $cog.$options{'extension'}");
		    foreach my $seq ( @{$list} )
		    {
                print OUTFILE ">$seq\n";
                my $residues = undef;
                if ($use_feature_ids_in_fasta) {
                    my $seq_id = $Feat->{$seq};
                    $logger->logdie("no sequence id found for feature with id=$seq") if (!defined($seq_id));
                    $residues = $Prot->{$seq_id};
                } else {
                    $residues = $Prot->{$seq};
                }
                $logger->logdie("no sequence data found for seq=$seq") if (!defined($residues));
                print OUTFILE $residues."\n";
		    }
		    close( OUTFILE );
		}
		else{
		    open( OUTFILE, ">$outDir/$cog.$options{'outputToken'}.$options{'extension'}" ) or $logger->logdie("could not open $cog.$options{'extension'}");
		    foreach my $seq ( @{$list} )
		    {
                print OUTFILE ">$seq\n";
                print OUTFILE "X\n";
		    }
		    close( OUTFILE );
		}
    }
}

sub createPolypeptideLookup
{
    my $seqRef = shift;    

    # We're only interested in polypeptide sequences for all-vs-all and pblast

    if( ($seqRef->returnattr( 'molecule' ) eq 'aa') || ($seqRef->returnattr('class') eq 'polypeptide'))
    {
        my $seq = $seqRef->subSequence(-1,0,0);
        
        my $identifier;
        if((defined $seqRef->{'BsmlSeqDataImport'}) && (!$use_feature_ids_in_fasta)){
            # don't think this is right, but trying to maintain reverse compatibility:
            $identifier = $seqRef->{'BsmlSeqDataImport'}->{'identifier'};
        }
        else{
            $identifier = $seqRef->returnattr( 'id' );
        }

        if( $identifier && $seq )
        {
            $Prot->{$identifier} = $seq;
        }
    }
}

sub createFeatureLookup
{
    my $featRef = shift;
    my $class = $featRef->returnattr('class');
    return unless ($class eq 'polypeptide');
    my $feat_id = $featRef->returnattr('id');
    my $links = $featRef->returnBsmlLinkListR();
    my @seqlinks = grep { $_->{'rel'} eq 'sequence' } @$links;
    my $nsl = scalar(@seqlinks);

    if ($nsl == 1) {
        my $rel = $seqlinks[0]->{'rel'};
        my $seq_id = $seqlinks[0]->{'href'};
        $seq_id =~ s/^\#//;
        if ($rel eq 'sequence') {
            $Feat->{$feat_id} = $seq_id;
        }
    } elsif ($nsl > 1) {
        $logger->logdie("feature $feat_id has $nsl sequence Links");
    }
}

