#!/usr/bin/perl

=head1  NAME 

RNAmmer2bsml.pl - convert RNAmmer output to BSML

=head1 SYNOPSIS

USAGE: RNAmmer2bsml.pl 
        --input=/path/to/rnammerfile 
        --output=/path/to/output.bsml
        --fasta_input=/path/to/fastafile
        --id_repository=/path/to/id_repository
        --project=aa1 
        [ --gzip_output=1       
          --log=/path/to/logfile 
          --debug=3
        ]

=head1 OPTIONS

B<--input,-i> 
    Input file file from a RNAmmer search.

B<--output,-o> 
    Output BSML file (will be created, must not exist)

B<--fasta_input,-a>
    The input file that was used as input for the RNAmmer run

B<--id_repository,-r>
    Path to the project's id repository

B<--gzip_output,-g>
    Optional. A non-zero value will make compressed output.
    
B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--project,-p> 
    Required. Project ID.  Used in creating feature ids. 
    If 'aa1' is passed, repeat feature IDs created are like aa1.repeat_region.4231.1
    If 'parse' is passed, will parse from already present ids.

B<--log,-l> 
    Log file

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used to convert the output from a RNAmmer search into BSML.

=head1 INPUT

RNAmmer can be run using multiple input sequences simultaneously, and this
script supports parsing single or multiple input result sets.  Usual RNAmmer
output looks like:

 seqname    source       feature  start   end     score  +/-  frame  attribute
 ---------  -----------  -------  ------  ------  -----  ---  -----  ---------
 Contig520	RNAmmer-1.2	 rRNA	    273038	273152	97.6	 +	   .	   5s_rRNA	
 Contig520	RNAmmer-1.2	 rRNA	    520669	520783	95.3	 +	   .	   5s_rRNA	
 Contig520	RNAmmer-1.2	 rRNA	     92664	 92778	97.6	 +	   .	   5s_rRNA	
 Contig520	RNAmmer-1.2	 rRNA	    280298	280412	97.6	 +	   . 	   5s_rRNA	
 Contig520	RNAmmer-1.2	 rRNA	    304453	304567	97.6	 +	   .	   5s_rRNA	
 Contig520	RNAmmer-1.2	 rRNA	     33753	 33867	97.6	 +	   .	   5s_rRNA	
 Contig520	RNAmmer-1.2	 rRNA	     87004	 87118	94.6	 +	   .	   5s_rRNA	


You can elimate the headers in the original RNAmmer output file by running 
RNAmmer using the -b option.  If they are present, they should be ignored by 
this script.

You define the input file using the --input option.  This file does not need any
special file extension.

=head1 OUTPUT

After parsing the input file, a file specified by the --output option is created.  This script
will fail if it already exists.  The file is created, and temporary IDs are created for
each result element.  They are only unique to the document, and will need to be replaced
before any database insertion.

#Base positions from the input file are renumbered so that positions start at zero.  The
#current output elements from RNAmmer that are not represented in the BSML file are:

  #  rRNA type
  #  Anti Codon
  #  Cove Score

#These need to be included later.

=head1 CONTACT

    Gaurav Jain
    gjain@som.umaryland.edu

=cut

use warnings;
use strict;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Ergatis::Logger;
use Ergatis::IdGenerator;
use Chado::Gene;
use BSML::GenePredictionBsml;
use BSML::BsmlBuilder;


my $input;      # Input file name
my $output;     # Output file name
my $fasta_input;   # fasta input to RNAmmer
my $idcreator;  # Ergatis::IdGenerator object
my $bsml;       # BSML::BsmlBuilder object
my $data = [];  # parsed RNAmmer data
my $debug = 4;  # debug value.  defaults to 4 (info)
my $parse_project =0;

my %options = ();
my $results = GetOptions (\%options, 
                          'input|i=s',
                          'output|o=s',
                          'fasta_input|a=s',
                          'project|p=s',
                          'log|l=s',
                          'id_repository|r=s',
                          'gzip_output|g=s',
                          'debug=s',
              'help|h') || &_pod;


# Set up the logger
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

## make sure all passed options are peachy
&check_parameters(\%options);

# Use the Gene.pm module if there exists a non-empty input
# otherwise just dump an 'empty' bsml file.
$data = &parse_RNAmmer_input( $input );
$bsml = &generateBsml( $data );
$bsml->writeBsml( $output, '', $options{'gzip_output'});
exit (0);


################ SUBS ######################

sub parse_RNAmmer_input 
{
	# Open the input file and parse the data into an array of gene models.
  my $infile = shift;
  my $genes;
  print("\n\ninput file is : $infile\n");

 
  ## open the input file for parsing
  open (IN,"< $infile") || $logger->logdie("can't open input file for reading");
	
	## Strip out uninteresting lines of data from the file
  #---- will uncomment later -----
  my %rawdata;
  while (<IN>) 
  {
    my @cols = split;
	    unless( $options{'project'} && $options{'project'} ne 'parse' ){
		          my $seq_id = $cols[0];
	            $options{'project'} = $1 if( $seq_id =~ /^([^\.]+)\./ );
	  		}
 
    #check whitespace, no warn
	    next if ( /^\s*$/ );
    
    ## make sure we don't parse the RNAmmer output header lines
	    next if ( /^\#\#gff-version2/i || /^sequence.*bounds.*cove/i || /^name.*end.*score/i || /\-\-.*\-\-/);


	## there should be 9 elements in cols, unless we have an unrecognized format.
	
	#	seqname    source       feature  start   end     score  +/-  frame  attribute
 	#	---------  -----------  -------  ------  ------  -----  ---  -----  ---------
 	#	Contig520	 RNAmmer-1.2	rRNA	   273038	 273152	 97.6	   +	   .	   5s_rRNA	
 	#	Contig520	 RNAmmer-1.2	rRNA	   520669	 520783	 95.3	   +	   .	   5s_rRNA	
 	
        unless (scalar @cols == 9) {
            $logger->error("the following RNAmmer_scan line was not recognized and could not be parsed:\n$_\n") 
                if ($logger->is_error);
            next;
        }
    
        ## add this data row to this sequence
        push( @{$rawdata{shift @cols}}, \@cols );
        
   }
   
   close IN;
   
       ## loop through each of the matches that we found and make a entries in the genes array.
    	 for my $seqid (keys %rawdata) {

        ## loop through each array reference of this key, adding to the data array as necessary.
        foreach my $arr ( @{$rawdata{$seqid}} ) 
        {

						#print "arr[2]= $$arr[2]\n arr[3] = $$arr[3]\n\n";
            ## Determine strandedness
            my $complement = ($$arr[2] > $$arr[3]) ? 1 : 0;

            ## correct for interbase numbering
            if ($complement) {
                $$arr[3]--;     ## End
            } else {
                $$arr[2]--;     ## begin
            }
            
            ## First, create the gene model object
            my $currGene = new Chado::Gene ( $idcreator->next_id( 'type' => 'gene',
                                                           'project' => $options{project} ),
                                      ($complement) ? $$arr[3] : $$arr[2],
                                      ($complement) ? $$arr[2] : $$arr[3],
                                      $complement,
                                      $seqid
                                    );

            ## Next, with the same coords, create the rRNA feature
            my $rrna_feat_id = $currGene->addFeature( $idcreator->next_id ( 'type' => 'rRNA',
                                                         'project' => $options{project} ),
                                   ($complement) ? $$arr[3] : $$arr[2],
                                   ($complement) ? $$arr[2] : $$arr[3],
                                   $complement,
                                   'rRNA'
                                 );
            
            ## add score attribute to feature
            $currGene->addFeatureAttribute(
                                            $rrna_feat_id,
                                            'score',
                                             $$arr[4],
                                          );

            ## add strand attribute to feature
            #$currGene->addFeatureAttribute(
            #                                $rrna_feat_id,
            #                                'strand',
            #                                $$arr[5],
            #                             );
                                          
	    			## add frame attribute to feature
            #$currGene->addFeatureAttribute(
            #                               $rrna_feat_id,
            #                              'frame',
            #                              $$arr[6],
            #                            );                                          
            
            ## add gene_product_name("attribute" in RNAmmer output) attribute to feature
            $currGene->addFeatureAttribute(
                                            $rrna_feat_id,
                                            'gene_product_name',
                                            $$arr[7],
                                          );

	    
            
     				## an exon needs to be added.
                if ($complement) 
                {
                    ## just add the whole thing as an exon
                    &add_exon_and_cds($currGene, $$arr[2], $$arr[3], $complement);
                } 
                else 
                {
                    &add_exon_and_cds($currGene, $$arr[3], $$arr[2], $complement);
                }
          


            # Handle Group now:
            my $count = $currGene->addToGroup( $currGene->getId, { 'all' => 1} );
            &_die("Nothing was added to group") unless ($count);

            push (@{$genes}, $currGene);

        }
    }

    return $genes;
}

sub add_exon_and_cds {
# Add an exon and cds to the gene object.  These features should share 
# start and stop coordinates.
     my ($gm, $start, $stop, $complement) = @_;

     foreach my $type (qw(exon CDS)) {
         
         # Add the exon or CDS to the gene model object
         $gm->addFeature( $idcreator->next_id( 'type' => $type,
                                               'project' => $options{project} ),
                          ($complement) ? $stop : $start,
                          ($complement) ? $start : $stop,
                          $complement,
                          $type
                        );

     }

}

sub generateBsml {
    my $data = shift;
		print "data = $data\n\n";
    #Create the document
    my $gene_pred = new BSML::GenePredictionBsml( 'RNAmmer' );
    my $doc = $gene_pred->{'doc'};

    unless($data && @{$data} > 0) {
        $gene_pred->addSequence( '', $fasta_input );
    }

    foreach my $gene( @{$data} ) {
        $gene_pred->addGene($gene);
        my $addedTo = $gene_pred->addSequence($gene->{'seq'}, $fasta_input);
        $logger->logdie("Could not find identifier ".$gene->{'seq'}." in $fasta_input")
            unless($addedTo || !(-e $fasta_input));
        $logger->logdie("fasta_input file: '$fasta_input' does not exist") unless(-e $fasta_input);
        
    }
    
    return $gene_pred;

}


sub check_parameters {
# Parse options hash, make sure certain ones exist and have meaningful values.
# Also take the time to setup some globally used values.

    my $error = '';

    &_pod if $options{'help'}; 
 
    ## make sure input file was given and exists
    if ($options{'input'}) {

        if (! -e $options{'input'}) {
            $logger->logdie("input file $options{'input'} does not exist")
        }

        $input = $options{'input'};

    } else {
        $error .= "--input_file is a required opttion\n";
    }
    ## make sure output file was given and doesn't exist yet
    if ($options{'output'}) {

        $output = $options{'output'};

    } else {
        $error .= "--output_file is a required option\n";
    }

    ## Make sure we're given the input fasta
    if ($options{'fasta_input'}) {
        $fasta_input = $options{'fasta_input'};
    } else {
        $error .= "--fasta_input is a required option\n";
    } 
    
    ## Now set up the id generator stuff
    if ($options{'id_repository'}) 
    {

        # we're going to generate ids
        $idcreator = new Ergatis::IdGenerator('id_repository' => $options{'id_repository'});

        #---------- Set the pool size ---------------
        #$idcreator->set_pool_size('gene'=>30,'tRNA'=>30,'exon'=>30,'CDS'=>30);

    } 
    else 
    {
        $error .= "--id_repository is a required option\n";
    }

    # The debug option is not required... but let's set it up if it's been given
    if ($options{'debug'}) 
    {
        $debug = $options{'debug'};
    }

    ## set default if project not passed.
    $options{project} ||= 'unknown';

    &_die($error) if $error;

    return 1;

} # END OF check_parameters 

sub _pod 
{
# Used to display the perldoc help.
    pod2usage( {-exitval => 0, -verbose => 2} );
}

sub _die 
{
# Kinda like a hitman, this is called when something needs to DIE.
    my $msg = shift;
    $logger->logdie($msg);
}
