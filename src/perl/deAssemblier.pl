#! /local/perl/bin/perl
##

use strict;
use warnings;
use Pod::Usage;
use BSML::BsmlParserSerialSearch;
use BSML::BsmlReader;
use BSML::BsmlBuilder;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Data::Dumper;

my %options = ();
my $results = GetOptions( \%options, 'bsmlAssembly|b=s', 'startCoord|s=s', 'windowSize|w=s', 'output_fasta_file|f=s', 'output_bsml_prefix|p=s', 'output_maps_prefix|m=s', 'help|h', 'man' ) || pod2usage();

my $Parser = new BSML::BsmlParserSerialSearch( SequenceCallBack => \&sequenceHandler, GenomeCallBack => \&genomeHandler );
my $Reader = new BSML::BsmlReader;

my $bsmlDoc;
my $startCoord = $options{'startCoord'};
my $endCoord = $options{'windowSize'} + $options{'startCoord'};
my $subAssemblyCount = 0;
my $end = 0;

my $bsmlMap = new BSML::BsmlBuilder;

open( FASTA_FILE, ">>$options{'output_fasta_file'}" ) or die "Could not open fasta file for concatenation.\n";

while( !($end) )
{
    $bsmlDoc = new BSML::BsmlBuilder;
   
    $Parser->parse( $options{'bsmlAssembly'} );

    $bsmlDoc->write( "$options{'output_bsml_prefix'}_$subAssemblyCount.bsml" );

    $startCoord = $startCoord += $options{'windowSize'};
    $endCoord = $startCoord + $options{'windowSize'};
    $subAssemblyCount++;

}

$bsmlMap->makeCurrentDocument();
$bsmlMap->write( "$options{'output_maps_prefix'}.bsml" );

close FASTA_FILE;

sub sequenceHandler
{
    my $seqref = shift;
    
    # determine if the sequence is an assembly sequence or an associated protein sequence
    if( $seqref->returnattr( 'molecule' ) eq 'dna' )
    {
	my $rhash = $Reader->readSequence( $seqref );

	if( $endCoord > $rhash->{'length'} )
	{ 
	    $endCoord = $rhash->{'length'};
	    $end = 1;
	}

	# if the CONTIG has no been added to the chunk map, add it now...

	$bsmlMap->makeCurrentDocument();

	if( !( BSML::BsmlDoc::BsmlReturnDocumentLookup( $rhash->{'id'} ) ) )
	{
	    my $newSeq = $bsmlMap->createAndAddExtendedSequenceN( id => $rhash->{'id'}, 
								  title => $rhash->{'id'}, 
								  molecule => $rhash->{'molecule'},
								  length => $rhash->{'length'} );

	    $bsmlMap->createAndAddSeqDataImportN( seq => $newSeq,
						  format => "BSML",
                                                  sorce =>  $options{'bsmlAssembly'},
                                                  id => "_$rhash->{'id'}" );
	}

        my $newSeq = $bsmlMap->createAndAddExtendedSequenceN( id => "$rhash->{'id'}_$subAssemblyCount", 
							      title => "$rhash->{'id'}_$subAssemblyCount", 
							      molecule => $rhash->{'molecule'},
							      length => $endCoord - $startCoord );

        $bsmlMap->createAndAddNumbering( seq => $newSeq,
                                         seqref => $rhash->{'id'},
                                         refnum => $startCoord,
                                         ascending => '1' );

	$bsmlDoc->makeCurrentDocument();


	# print "DEBUG: $startCoord $endCoord\n";

	my $seqdat = $Reader->subSequence( $seqref, $startCoord, $endCoord, 0 );

	print FASTA_FILE ">$options{'output_bsml_prefix'}_$subAssemblyCount\n";
	print FASTA_FILE "$seqdat\n";

	$newSeq = $bsmlDoc->createAndAddExtendedSequenceN( id => $rhash->{'id'}."_$subAssemblyCount", 
						 title => $rhash->{'id'}."_$subAssemblyCount", 
						 molecule => $rhash->{'molecule'},
						 length => ($endCoord - $startCoord) );
					  
	$newSeq->addBsmlAttr( 'ASSEMBLY', $rhash->{'id'}."_$subAssemblyCount" );

	foreach my $feature_table( @{$Reader->returnAllFeatureTables( $seqref )} )
	{
	    my $newFTable = $bsmlDoc->createAndAddFeatureTableN( seq => $newSeq,
								id => $feature_table->returnattr( 'id' ) );

	    # loop through all the features on the feature table, transforming the coordinates of those in the region of interest
	    # and adding to the bsml document

	    foreach my $feature ( @{$Reader->readFeatures($feature_table)} )
	    {
		if( $feature->{'class'} eq 'GENE' )
		{
		    my $fmin = $feature->{'locations'}->[0]->{'startpos'};
		    my $fmax = $feature->{'locations'}->[0]->{'endpos'};
		    my $complement = $feature->{'locations'}->[0]->{'complement'};

		    if( $fmin > $fmax )
		    { 
			($fmin, $fmax, $complement) = ($fmax, $fmin, 1);
		    }

		    if( ($complement == 0 && ($fmin >= $startCoord && $fmin <= $endCoord )) || ($complement == 1 && ($fmax >= $startCoord && $fmax <= $endCoord ))  )
		    {
			my $mapped_fmin = $fmin - $startCoord;
			my $mapped_fmax = $fmax - $startCoord;

			my $newFeat = $bsmlDoc->createAndAddFeatureWithLocN( FTable => $newFTable,
									     id => $feature->{'id'},
									     title => $feature->{'title'},
									     class => $feature->{'class'},
									     comment => $feature->{'comment'},
									     displayAuto => $feature->{'display-auto'},
									     start => $mapped_fmin,
									     end => $mapped_fmax,
									     complement => $complement );

			foreach my $qualifier ( @{$feature->{'qualifiers'}} )
			{
			    my $qual = $bsmlDoc->createAndAddQualifier( $newFeat, $qualifier->{'key'}, $qualifier->{'value'} );
			}

			foreach my $bsmlattr ( @{$feature->{'bsmlattrs'}} )
			{
			    $newFeat->addBsmlAttr($bsmlattr->{'key'}, $bsmlattr->value );
			}

			$bsmlDoc->makeCurrentDocument();
			BSML::BsmlDoc::BsmlSetDocumentLookup( $feature->{'id'}, $newFeat );
		    }
		}

		if( $feature->{'class'} eq 'TRANSCRIPT' ) 
		{
		    my $fmin = $feature->{'locations'}->[0]->{'startpos'};
		    my $fmax = $feature->{'locations'}->[1]->{'endpos'};
		    my $complement = $feature->{'locations'}->[0]->{'complement'};
		    
		    if( $fmin > $fmax )
		    { 
			($fmin, $fmax, $complement) = ($fmax, $fmin, 1);
		    }

		    if(  ($complement == 0 && ($fmin >= $startCoord && $fmin <= $endCoord )) || ($complement == 1 && ($fmax >= $startCoord && $fmax <= $endCoord )) )
		    { 

		    my $mapped_fmin = $fmin - $startCoord;
		    my $mapped_fmax = $fmax - $startCoord;

		    my $newFeat = $bsmlDoc->createAndAddFeatureN( FTable => $newFTable,
								  id => $feature->{'id'},
								  title => $feature->{'title'},
								  class => $feature->{'class'},
								  comment => $feature->{'comment'},
								  displayAuto => $feature->{'display-auto'});
					
		    $newFeat->addBsmlSiteLoc( $mapped_fmin, $complement, 'START' ); 
		    $newFeat->addBsmlSiteLoc( $mapped_fmax, $complement, 'STOP' ); 

		    foreach my $qualifier ( @{$feature->{'qualifiers'}} )
		    {
			my $qual = $bsmlDoc->createAndAddQualifier( $newFeat, $qualifier->{'key'}, $qualifier->{'value'} );
		    }
		    
		    foreach my $bsmlattr ( @{$feature->{'bsmlattrs'}} )
		    {
			$newFeat->addBsmlAttr($bsmlattr->{'key'}, $bsmlattr->value );
		    }
		    
		    $bsmlDoc->makeCurrentDocument();
		    BSML::BsmlDoc::BsmlSetDocumentLookup( $feature->{'id'}, $newFeat );
		}
		}

		if( $feature->{'class'} eq 'EXON' )
		{
		    my $fmin = $feature->{'locations'}->[0]->{'startpos'};
		    my $fmax = $feature->{'locations'}->[0]->{'endpos'};
		    my $complement = $feature->{'locations'}->[0]->{'complement'};
		    
		    if( $fmin > $fmax )
		    { 
			($fmin, $fmax, $complement) = ($fmax, $fmin, 1);
		    }
		    
		    if(  ($complement == 0 && ($fmin >= $startCoord && $fmin <= $endCoord )) || ($complement == 1 && ($fmax >= $startCoord && $fmax <= $endCoord )) )
		    {
			my $mapped_fmin = $fmin - $startCoord;
			my $mapped_fmax = $fmax - $startCoord;
			
			my $newFeat = $bsmlDoc->createAndAddFeatureWithLocN( FTable => $newFTable,
									     id => $feature->{'id'},
									     title => $feature->{'title'},
									     class => $feature->{'class'},
									     comment => $feature->{'comment'},
									     displayAuto => $feature->{'display-auto'},
									     start => $mapped_fmin,
									     end => $mapped_fmax,
									     complement => $complement );
			
			foreach my $qualifier ( @{$feature->{'qualifiers'}} )
			{
			    my $qual = $bsmlDoc->createAndAddQualifier( $newFeat, $qualifier->{'key'}, $qualifier->{'value'} );
			}
			
			foreach my $bsmlattr ( @{$feature->{'bsmlattrs'}} )
			{
			    $newFeat->addBsmlAttr($bsmlattr->{'key'}, $bsmlattr->value );
			}
			
			$bsmlDoc->makeCurrentDocument();
			BSML::BsmlDoc::BsmlSetDocumentLookup( $feature->{'id'}, $newFeat );
		    }
		}

		if( $feature->{'class'} eq 'CDS' ) 
		{
		    my $fmin = $feature->{'locations'}->[0]->{'startpos'};
		    my $fmax = $feature->{'locations'}->[0]->{'endpos'};
		    my $complement = $feature->{'locations'}->[0]->{'complement'};
		    
		    if( $fmin > $fmax )
		    { 
			($fmin, $fmax, $complement) = ($fmax, $fmin, 1);
		    }
		    
		    if(  ($complement == 0 && ($fmin >= $startCoord && $fmin <= $endCoord )) || ($complement == 1 && ($fmax >= $startCoord && $fmax <= $endCoord )) )
		    { 

		    my $mapped_fmin = $fmin - $startCoord;
		    my $mapped_fmax = $fmax - $startCoord;

		    my $newFeat = $bsmlDoc->createAndAddFeatureN( FTable => $newFTable,
								  id => $feature->{'id'},
								  title => $feature->{'title'},
								  class => $feature->{'class'},
								  comment => $feature->{'comment'},
								  displayAuto => $feature->{'display-auto'});
					
		    $newFeat->addBsmlSiteLoc( $mapped_fmin, $complement, 'START' ); 
		    $newFeat->addBsmlSiteLoc( $mapped_fmax, $complement, 'STOP' ); 

		    foreach my $qualifier ( @{$feature->{'qualifiers'}} )
		    {
			my $qual = $bsmlDoc->createAndAddQualifier( $newFeat, $qualifier->{'key'}, $qualifier->{'value'} );
		    }
		    
		    foreach my $bsmlattr ( @{$feature->{'bsmlattrs'}} )
		    {
			$newFeat->addBsmlAttr($bsmlattr->{'key'}, $bsmlattr->{'value'} );
		    }

		    foreach my $bsmllink ( @{$feature->{'bsmllinks'}} )
		    {
			$newFeat->addBsmlLink($bsmllink->{'rel'}, $bsmllink->{'href'} );

			$bsmllink->{'href'} =~ s/\#//;

			my $protSeq = $bsmlDoc->createAndAddSequence( $bsmllink->{'href'}, '', '', 'aa' );

                        $protSeq->addBsmlAttr( 'ASSEMBLY', "$rhash->{'id'}_$subAssemblyCount" );
			
			$bsmlDoc->makeCurrentDocument();
			BSML::BsmlDoc::BsmlSetDocumentLookup( $bsmllink->{'href'}, $protSeq );
		    }
		    
		    $bsmlDoc->makeCurrentDocument();
		    BSML::BsmlDoc::BsmlSetDocumentLookup( $feature->{'id'}, $newFeat );

		    
		}
		}
	    }
	}

	$bsmlDoc->makeCurrentDocument();

	foreach my $feature_group ( @{$seqref->returnBsmlFeatureGroupListR()} )
	{
	    $rhash = $Reader->readFeatureGroup( $feature_group );

	    if( BSML::BsmlDoc::BsmlReturnDocumentLookup( $rhash->{'group-set'} ) )
	    {
		my $newFGroup = $bsmlDoc->createAndAddFeatureGroup( $newSeq, '', $rhash->{'group-set'} );

		foreach my $member (@{$rhash->{'feature-members'}})
		{
		    $bsmlDoc->createAndAddFeatureGroupMember( $newFGroup, $member->{'featref'}, $member->{'feature-type'}, $member->{'group-type'}, $member->{'cdata'} );
		}
	    }
	}
    }
    else
    {
      if( $seqref->returnattr( 'molecule' ) eq 'aa' )
      {
	  # if the sequence is a protein sequence, determine if it is in the search window, if it is add the sequence to the 
	  # subAssembly document

	  $bsmlDoc->makeCurrentDocument();
	  if( my $protSeq = BSML::BsmlDoc::BsmlReturnDocumentLookup( $seqref->returnattr( 'id' ) ) )
	  {
	      $protSeq->addattr( 'title', $seqref->returnattr( 'title' ) );
	      $protSeq->addattr( 'molecule', $seqref->returnattr( 'molecule' ) );
	      $protSeq->addattr( 'length', $seqref->returnattr( 'length' ));

	      $protSeq->setBsmlSeqData( $seqref->returnSeqData() );				 
	  }
      }
    }

  # if the sequence is an assembly sequence. Replicate the sequence attributes in the subAssembly. 
}

sub genomeHandler
{
  # replicate the genome object in the subAssembly. 
    my $genome = shift;
    my $rhash = $Reader->readGenome( $genome );

    my $organism = $bsmlDoc->createAndAddOrganism( genome => $bsmlDoc->createAndAddGenome(),
						   species => $rhash->{'species'},
						   genus => $rhash->{'genus'} );

    my $strain = $bsmlDoc->createAndAddStrain( organism => $organism,
					       name => $rhash->{'strain'},
					       database => $rhash->{'database'},
					       source_database => $rhash->{'source_database'} );
						   
}

