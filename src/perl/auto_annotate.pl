#!/usr/local/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

#auto_annotate --input_file /usr/local/annotation/scratch/kgalens/PROK/output_repository/glimmer3/5559_iter2/glimmer3.bsml.list --hmm_analysis /usr/local/annotation/scratch/kgalens/PROK/output_repository/hmmpfam/5557_protein_ALL_LIB_bin.HMM/ --ber_analysis /usr/local/annotation/scratch/kgalens/PROK/output_repository/ber/5576_AllGroup.niaa/ --output_bsml test.scaff1.output.bsml --hmm_info_db /usr/local/annotation/scratch/kgalens/hmmPandaInfo/hmmInfo.db --panda_header_offsets /usr/local/annotation/scratch/kgalens/hmmPandaInfo/pandaTest --panda_header_file /usr/local/annotation/scratch/kgalens/hmmPandaInfo/pandaTest_header --role_info_db /usr/local/annotation/scratch/kgalens/hmmPandaInfo/grData 

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Workflow::Logger;
use Pod::Usage;
use Data::Dumper;
use File::Find;
use MLDBM 'DB_File';
use XML::Twig;
use BSML::BsmlBuilder;
use Split_DB_File;
require "autoAnnotate.data";
require "/export/prog/autoAnnotate/sharedRoutines.pl";
$|++;


##################### GLOBALS AND CONSTANTS #######################################
my @inputFiles;                             #The list of input files
my $outputBsml;                             #The output bsml file name
my (@hmmBsmlFiles, @berBsmlFiles);          #Analysis file lists
my $gzipOutput = 0;                         #Flag indicating compressed output
my $evidence;                               #Hash ref holding evidence information for models
my $finalAnnotes;                           #The actual annotation decision
my %polypeptideIds;                         #polypeptide ids as key(within one bsml document)
                                            #and gene length as value
my %hmmInfo;                                #Tied to hmm_info_db file (DB_File format)
my $ber_map;                                #Map of CDS (key) to polypeptide and transcript ids.
my %pandaHeaderOffsets;                     #Byte offsets of headers in the headerFile
                                            #Split_DB_File tied hash
my $headerFile;                              
my $per_id_cutoff = 35;                     #Default cutoff for ber hit (percent identity)
my $length_cutoff = 80;                     #Default cutoff for ber hit (length of hit)
my %roleInfo;                               #Holds role guessing data.
my $asmbl;                                  #Contains assembly information
my $outputDir;                              #Directory for output files.
my $geneBoundaries;                         #Hash containing start and end positions keyed by id
my $extension = 300;                        #Extension length for ber (default 300).
use vars qw(%isoType);
my ($acc,$name,$ec_num,$gene_sym,$species,$go_id,$role_id,$exp,$wgp,$cg,$rank) = (0..10);
###################################################################################

###################### OPTION GATHERING AND LOG SETUP ##############################
my %options = ();
my $results = GetOptions (\%options, 
                          'input_list|l=s',
                          'input_file|i=s',
                          'output_bsml|o=s',
                          'output_dir=s',
                          'hmm_analysis|n=s',
                          'ber_analysis|b=s',                          
                          'gzip_output|g=s',
                          'hmm_info_db|m=s',
                          'panda_header_offsets|p=s',
                          'panda_header_file|a=s',
                          'role_info_db|r=s',
                          'percent_id_cutoff=s',
                          'length_cutoff=s',
                          'ber_extension_length|e=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || &_pod;

#Setup the logger
my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

#Check the options
&checkParameters(\%options);
###################################################################################

############################## MAIN ###############################################

#Loop through all of the input files ( bsml with feature: class=polypeptide ).
foreach my $input(@inputFiles) {

    #Make a list of all of the polypeptide id's to be automagically annotated.
    print "Getting polypeptides from $input\n";
    %polypeptideIds = &getPolypeptideIds( $input );

    unless($outputBsml) {
        my $baseName = $1 if($input =~ /^.*\/([^\/]*?)\.bsml/);
        $outputBsml = $outputDir.$baseName.".auto_annotate.bsml";
    }

    unless($outputBsml) {
        $logger->logdie("Output file has not been set");
    }

    my $count = 0;
  
    #Cycle through all these polypeptides
    foreach my $polypeptide ( keys %polypeptideIds ) {       
        next unless($polypeptide eq 'prok.polypeptide.156087.1');
        print "***Annotating protein: $polypeptide***\n";# if($logger->is_debug);

        #Get the best HMM evidence
        $evidence->{$polypeptide}->{'HMM'} = &getHMMEvidence( $polypeptide );

        #Get the best BER evidence unless we've found a great HMM hit.
        #HMM hits are more reliable (if they are equivalog level).
        $evidence->{$polypeptide}->{'BER'} = &getBEREvidence( $polypeptide )
            unless(exists $evidence->{$polypeptide}->{'HMM'}->{'rank'} &&
                   $evidence->{$polypeptide}->{'HMM'}->{'rank'} == 1);
      
        $finalAnnotes->{$polypeptide} = &annotate( $evidence->{$polypeptide} );
        
        $count++;
    }

    print "$count proteins annotated\n";
    &annotation2bsml($outputBsml);
    $finalAnnotes = {};
    $outputBsml = "";
}

########################### SUB-ROUTINES #################################################
sub checkParameters {
    my $opts = shift;
    
    #Make sure an input_file or input_list is supplied
    if( defined($opts->{'input_file'}) ) {
        $logger->logdie("Option input_file ($opts->{'input_file'}) does not exist") unless( -e $opts->{'input_file'} );
        push(@inputFiles, $opts->{'input_file'});
    } elsif( defined($opts->{'input_list'}) ) {
        $logger->logdie("Option input_list ($opts->{'input_list'}) does not exist") unless( -e $opts->{'input_list'} );
        open( FILELIST, "< $opts->{input_list}") or
            $logger->logdie("Unable to open $opts->{input_list} ($!)");
        chomp(@inputFiles = <FILELIST>);
        close(FILELIST);
    } else {
        $logger->logdie("Either option input_list or input_file is required");
    }
    
    #Make sure they specified the output bsml file.
    if($opts->{'output_bsml'}) {
        $outputBsml = $opts->{'output_bsml'};
    } elsif($opts->{'output_dir'}) {
        $outputDir = $opts->{'output_dir'};
        $outputDir .= '/' unless( $outputDir =~ m|/$|);
    } else {
        $logger->logdie("Option output_bsml or output_dir is required");
    }
  

    #Make sure the hmm analysis was passed in
    if( $opts->{'hmm_analysis'} ) {
        @hmmBsmlFiles = &getAnalysisFiles( $opts->{'hmm_analysis'} );
    }

    #Make sure the ber analysis was passed in
    if( $opts->{'ber_analysis'} ) {
        @berBsmlFiles = &getAnalysisFiles( $opts->{'ber_analysis'} );
    }

    #Should the output be compressed?
    if( $opts->{'gzip_output'} ) {
        $gzipOutput = 1;
    }

    #Tie the hash to the HMM info
    if( $opts->{'hmm_info_db'} ) {
        tie(%hmmInfo, 'MLDBM', $opts->{'hmm_info_db'})
            or $logger->logdie("Unable to tie hash to $opts->{'hmm_info_db'}");

    } else {
        $logger->logdie("Option hmm_info_db is required.  See documentation for format".
                        "requirements");
    }

    #Get the panda header db
    if( $opts->{'panda_header_offsets'} ) {
        tie(%pandaHeaderOffsets, 'Split_DB_File', $opts->{'panda_header_offsets'}) 
            or $logger->logdie("Unable to tie hash to $opts->{'panda_header_offsets'}");
    } else {
        $logger->logdie("Option panda_header_offsets is required.  See documentation for format");
    }

    #Get the header text file or set it to the default of whatever the offset db is named plus
    #the string _headers.
    if( $opts->{'panda_header_file'} ) {
        $headerFile = $opts->{'panda_header_file'};
    } else {
        $headerFile = $opts->{'panda_header_offsets'}."_header";
    }

    my $roleFile;
    #Get the role info database
    if($opts->{'role_info'}) {
        $roleFile = $opts->{'role_info'};
    } else {
        $roleFile = "/export/prog/autoAnnotate/grData";
    }
    tie(%roleInfo, 'Split_DB_File', $roleFile) or
        $logger->logdie("Cannot tie hash to $roleFile");

    #Set percent identity cutoff scores for ber matches
    if( $opts->{'percent_id_cutoff'} ) {
        $per_id_cutoff = $opts->{'percent_id_cutoff'};
    }

    #Set length cutoff score for ber matches
    if( $opts->{'length_cutoff'} ) {
        $length_cutoff = $opts->{'length_cutoff'};
    }

    #Set the extension length for the ber run
    if($opts->{'ber_extension_length'}) {
        $extension = $opts->{'ber_extension_length'};
    }
    

}

sub annotate {
    my $evidence = shift;
    my $retval = {};

    print Dumper($evidence) if($logger->is_debug);

    #Check for an HMM hit of level 1
    if(exists $evidence->{'HMM'}->{'rank'} && 
       ($evidence->{'HMM'}->{'rank'} == 1 || $evidence->{'HMM'}->{'rank'} == 2)) {
        print "Taking evidence from equivalog HMM\n" if($logger->is_debug);
        $retval = {
            'annotation_from' => $evidence->{'HMM'}->{'hmm_acc'},
            'com_name'        => $evidence->{'HMM'}->{'hmm_com_name'},
            'gene_sym'        => $evidence->{'HMM'}->{'gene_sym'},
            'ec_num'          => $evidence->{'HMM'}->{'ec_num'}
        };

        $retval->{'role_id'} = $evidence->{'HMM'}->{'role_id'} if(defined $evidence->{'HMM'}->{'role_id'});
        $retval->{'go_id'} = $evidence->{'HMM'}->{'go_term'} if(defined $evidence->{'HMM'}->{'go_term'});

        unless($retval->{'role_id'} && $retval->{'go_id'}) {

            print Dumper($evidence->{'HMM'});
        }


        $logger->warn("No role_id ".$evidence->{'HMM'}->{'hmm_acc'}) unless($retval->{'role_id'});
        $logger->warn("No go_term ".$evidence->{'HMM'}->{'hmm_acc'}) unless($retval->{'go_id'});

        my $autoAnnotate = 0;
        $autoAnnotate = 1 if($evidence->{'HMM'}->{'rank'} == 1);
    } elsif($evidence->{'BER'}->[$rank] && $evidence->{'BER'}->[$rank] != 100) {  #if we have a not crap ber hit
        print "Taking evidence from not crap BER hit\n" if($logger->is_debug);
        $retval->{'annotation_from'} = $evidence->{'BER'}->[$acc];
        $retval->{'com_name'} = $evidence->{'BER'}->[$name];
        $retval->{'gene_sym'} =  $evidence->{'BER'}->[$gene_sym];
        $retval->{'ec_num'} = $evidence->{'BER'}->[$ec_num];
        if ($evidence->{'BER'}->[$role_id]) {
            $retval->{'role_id'} = $evidence->{'BER'}->[$role_id];
            $retval->{'role_from'} = $evidence->{'BER'}->[$acc];
        }
        if (defined $evidence->{'BER'}->[$go_id] && $evidence->{'BER'}->[$go_id]) {
            $retval->{'go_id'} = $evidence->{'BER'}->[$go_id];
            $retval->{'go_from'} = $evidence->{'BER'}->[$acc];
        } elsif (defined $evidence->{'HMM'}->{'go_term'} && $evidence->{'HMM'}->{'goterm'}) {
            $retval->{'go_id'} = $evidence->{'HMM'}->{'go_term'};
            $retval->{'go_from'} = $evidence->{'HMM'}->{'hmm_acc'};
        }

    }
    print "BEFORE IF:".$retval->{'com_name'}."\n";
    if( !$retval->{'com_name'} || $retval->{'com_name'} =~ /hypothetical protein/i) {
        print "The com_name contains hypo...\n";
        if($evidence->{'HMM'}->{'hmm_com_name'}) {
            print "Setting com_name and annotation from from the HMM\n";
            $retval->{'com_name'} = $evidence->{'HMM'}->{'hmm_com_name'};
            $retval->{'annotation_from'} = $evidence->{'HMM'}->{'hmm_acc'};

            if(($evidence->{'HMM'}->{'rank'} == 4 || $evidence->{'HMM'}->{'rank'} == 5) &&
               $retval->{'com_name'} !~ /family/ &&
               $retval->{'com_name'} !~ /protein/ ) {
                
                $retval->{'com_name'} .= " subfamily";
            }
            if($evidence->{'HMM'}->{'rank'} == 3 || $evidence->{'HMM'}->{'rank'} == 5) {
                $retval->{'com_name'} .= ", putative";
            }

            if($evidence->{'HMM'}->{'rank'} == 6 &&
               $retval->{'com_name'} !~ /family/ &&
               $retval->{'com_name'} !~ /protein/) {
                $retval->{'com_name'} .= " superfamily";
            }

            if($evidence->{'HMM'}->{'rank'} == 7 &&
               $retval->{'com_name'} !~ /family/ &&
               $retval->{'com_name'} !~ /protein/) {
                $retval->{'com_name'} .= " family";
            }

            if($evidence->{'HMM'}->{'rank'} == 8) {
                if($retval->{'com_name'} !~ /domain/ &&
                   $retval->{'com_name'} !~ /family/ &&
                   $retval->{'com_name'} !~ /protein/) {
                    $retval->{'com_name'} .= " domain protein";
                } elsif($retval->{'com_name'} !~ /protein/) {
                    $retval->{'com_name'} .= " protein";
                }
            }
        }

        $retval->{'gene_sym'} =  $evidence->{'HMM'}->{'gene_sym'};
        $retval->{'ec_num'} = $evidence->{'HMM'}->{'ec_num'};
        if (defined $evidence->{'HMM'}->{'role_id'}) {
            $retval->{'role_id'} = $evidence->{'HMM'}->{'role_id'};
            $retval->{'role_from'} = $evidence->{'HMM'}->{'acc'};
        }
     
        if ($evidence->{'HMM'}->{'go_id'}) {
            $retval->{'go_id'} = $evidence->{'HMM'}->{'go_term'};
            $retval->{'go_from'} = $evidence->{'HMM'}->{'hmm_acc'};
        }
    }

    unless( defined $retval->{'com_name'} ) {
        $retval->{'com_name'} = 'hypothetical protein';
    }

    if ($retval->{'com_name'} =~ /hypothetical protein|\, putative|\-related|unknown/i) { 
        $retval->{'gene_sym'} = ""; 
    }

    if ($retval->{'com_name'} =~ /^conserved hypothetical protein$/i) { 
        $retval->{'role_id'} = 156;
        $retval->{'go_id'} = "GO:0000004 GO:0005554";
        $retval->{'gene_sym'} = "";
        $retval->{'ec_num'} = "";
    }elsif ($retval->{'com_name'} =~ /^hypothetical protein$/i) {
        $retval->{'gene_sym'} = "";
        $retval->{'role_id'} = "";
        $retval->{'ec_num'} = "";
    } elsif (!$retval->{'role_id'} || $retval->{'role_id'} eq "NULL" || $retval->{'role_id'} == 185) {
        
        print "Trying to find a role id\n" if($logger->is_debug);
        $retval->{'role_id'} = &guessRole($retval->{'com_name'}, $retval->{'gene_sym'}, $retval->{'ec_num'});
        $retval->{'role_from'} = "guess_role";
    }

    print "\n\nI choose you:\n" if($logger->is_debug);
    print Dumper($retval) if($logger->is_debug);
    return $retval;

}

sub guessRole {
    my ($com_name, $gene_sym, $ecs) = @_;
    my @result;

    if($ecs && $ecs ne 'NULL') {
        my @ec_nums = split(/\s+/,$ecs);
        foreach my $ec_num (@ec_nums) {
            print $roleInfo{$ec_num}."\n" if(defined($roleInfo{$ec_num}));
            push(@result, split(/\s+/, $roleInfo{$ec_num})) 
                if(defined($roleInfo{$ec_num}));
        }
    }

    if($gene_sym && $gene_sym ne 'NULL') {
        push(@result, split(/\t/, $roleInfo{$gene_sym})) 
            if(defined($roleInfo{$gene_sym}));
    }

    $com_name =~ s/-related//g;
    $com_name =~ s/,\s*putative//g;
    my $lower_com_name = $com_name;
    $lower_com_name =~ tr/A-Z/a-z/;

    push(@result, split(/\s+/, $roleInfo{$lower_com_name}))
        if(defined($roleInfo{$lower_com_name}));

    my %role;
    foreach my $r ( @result ) {
	my @pairs = split /\s+/, $r;
	foreach my $pair (@pairs) {
        	my ($role_id, $count) = split /\:/, $pair;
        	$role{$role_id}->{'hits'}++;
        	$role{$role_id}->{'count'} += $count;
		if($count eq '4 184') {
			die("$r\n@pairs");
		}
	}
    }

    my $topcount;
    my $most_hits;
    my @role_list;
    foreach my $r_id (sort {$role{$b}->{'hits'} <=> $role{$a}->{'hits'} ||
                               $role{$b}->{'count'} <=> $role{$a}->{'count'}} keys %role) {
        $most_hits = $role{$r_id} unless(defined($most_hits));
        $topcount = $role{$r_id} unless(defined($topcount));

        last if ($role{$r_id}->{'hits'} < $most_hits);

        if ($role{$r_id}->{'count'} > ($topcount * 2/3)) { 
            push(@role_list, $r_id);
        } else { 
            last; 
        }
        $most_hits = $role{$r_id}->{'hits'};
        $topcount = $role{$r_id}->{'count'} if (! $topcount);
    }

    # if we have too many roles, query wasn't specific enough
    if (@role_list > 3) { 
        @role_list = ();
    }
    
    # try finding roles by individual words
    if (@role_list == 0) {
        my %tmpRole;
        my @words = split /\W+/, $com_name;
        my @lower_words = split /\W+/, $lower_com_name;
        my (@condition, $condition);

        for ( my $i=0; $i<@words; $i++ ) {

            # ignore useless words
            next if ($lower_words[$i] =~ /protein|putative|domain|family|hypothetical|conserved|subunit|terminal|dependent|region|product|unknown|function|precursor|predicted|uncharacterized|like|homolog|probable|similar/);

            # ignore short words unless they are potentially gene syms
            next if ($lower_words[$i] eq $words[$i] && length($words[$i]) < 5);
            next if (length($words[$i]) < 3);
            my @result;

            foreach my $key (keys %roleInfo) {
                next unless($key =~ /^lcname(.*)/);
                my $tmpStr = quotemeta($1);
                if($words[$i] =~ /$tmpStr/) {
                    push(@result, split(/\t/, $roleInfo{$key}));
                }
            }

            foreach my $r(@result) {
                my ($role_id, $count) = split /\:/, $r;
                $tmpRole{$role_id}->{'hits'}++;
                $tmpRole{$role_id}->{'count'} += $count;
            }
        }
	

    	# now analyze ROLE hash to find roles that hit the most terms and had the highest counts
        my $topcount;
        my $most_hits;
        foreach my $role_id (sort {$role{$b}->{'hits'} <=> $role{$a}->{'hits'} ||
                                   $role{$b}->{'count'} <=> $role{$a}->{'count'}} keys %role) {
            $most_hits = $role{$role_id} unless(defined($most_hits));
            $topcount = $role{$role_id} unless(defined($topcount));
            last if ($role{$role_id}->{'hits'} < $most_hits);
            if ($role{$role_id}->{'count'} > ($topcount * 2/3)) { 
                push @role_list, $role_id; 
            } else { 
                last; 
            }

            $most_hits = $role{$role_id}->{'hits'};
            $topcount = $role{$role_id}->{'count'} if (! $topcount);
        }

        if (@role_list > 3 or @role_list == 0) {
            if    ($com_name =~ /conserved hypothetical protein/i) { return 156 }
            elsif ($com_name =~ /unknown function/i ||
                   $com_name =~ /family/i) { return 157 }
            return(185);
        } else { 
            return join " ", @role_list;
        }

    } else { 
        return join " ", @role_list; 
    }

}

sub getBEREvidence {
    my $polyId = shift;
    my $berBsmlFile;
    my @berMatches;

    #This will hold a mapping of the bsml attribute 'id' and the real id
    #used by panda.  The bsml attribute 'id' is a mangled version and frankly
    #quite useless (basically, the id without any of the '|', which turns out, 
    #is important)
    my $fakeId2RealId = {};

    #Find the associated file.
    foreach my $bFile ( @berBsmlFiles ) {
        if($bFile =~ /($polyId)/) {
            $berBsmlFile = $bFile;
            last;
        }
    }

    return [] unless($berBsmlFile);

    #Now parse the bsml file.
    my $twig = new XML::Twig( TwigHandlers => 
                              { 
                                  'Sequence'           => 
                                      sub { &mapBsmlId2RealId( @_, $fakeId2RealId ) }, 
                                  'Seq-pair-alignment' =>  sub { &handleBERSpa( @_, \@berMatches )} 
                              }         
                              );
    

    print "\tparsing: $berBsmlFile\n" if($logger->is_debug);
    $twig->parsefile($berBsmlFile);

    return [] if(@berMatches == 0);

    #Now that we have all the berMatches in this array, lets find the best one.
    my $bestMatch = &getBestBERMatch( \@berMatches, $fakeId2RealId );
    
    return $bestMatch;

}

#From getBEREvidence. Twig Handler subroutine.
#spa = Sequence Pair Alignment
#spr = Seq-pair-run
sub handleBERSpa {
    my ($twig, $spa, $berMatches) = @_;
    my $match;

    my $compseq = $spa->att('compseq');
   
    foreach my $spr ( $spa->children('Seq-pair-run') ) {
        my $hitLength = $spr->att('runlength');
        my ($pId, $pVal);
        foreach my $attElem ( $spr->children('Attribute') ) {
            $pId = $attElem->att('content') if($attElem->att('name') eq 'percent_identity');
            $pVal = $attElem->att('content') if($attElem->att('name') eq 'p_value');
        }
        my $id = $spa->att('refseq');

        #Find the length of the extended gene
        my $bounds = $geneBoundaries->{ $ber_map->{$id} };
        my ($start,$end) = ($bounds->{'startpos'} - $extension, 
                            $bounds->{'endpos'} + $extension);
        if($asmbl->{'topology'} && $asmbl->{'topology'} =~ /circular/i) {
            if($end > $asmbl->{'length'}) {
                $logger->logdie("End position is ran off end of molecule");
            }
        } else {
            $start = 0 if($start < 0);
            $end = $asmbl->{'length'} if($end > $asmbl->{'length'});

        }
        
        my $gene_length = $end-$start;

        $match = {
            'hitLength'       => $hitLength,
            'fraction_length' => ($hitLength/$gene_length),
            'bit_score'       => $spr->att('runscore'),
            'percent_id'      => $pId,
            'p_value'         => $pVal,
            'id'              => $id,
            'gene_length'     => $gene_length,
            'compseq'         => $compseq,
            'subject_length'  => $spr->att('runprob')
        };

        #Weed out bad hits.
        if($match->{'percent_id'} > 30 &&  
           $match->{'bit_score'} * $match->{'fraction_length'} > 4) {
            push(@{$berMatches}, $match);
        } elsif( scalar @{$berMatches} == 0 ) {
            $match = 
            { 'compseq'         => $compseq,
              'p_value'         => $pVal,
              'fraction_length' => ($hitLength/$gene_length),
              'percent_id'      => $pId,
              'bit_score'       => $spr->att('runscore')
          };
            push(@{$berMatches}, $match);
        }
        $match = {};
        
    }
    
}

sub getBestBERMatch {
    my ($berMatches, $fakeId2RealId) = @_;

    #Sort the ber matches by p-value
    my @sortedBerMatches = sort { $a->{'p_value'} <=> $b->{'p_value'} } @{$berMatches};
    my $best_pvalue = $sortedBerMatches[0]->{'p_value'};

    #This will store the current best ber match
    my $bestBerMatch = [];
    $bestBerMatch->[$name] = "";
    my $bestBerRank = 100;
    my $bestHitInfo = {};
    
    foreach my $berMatch ( @sortedBerMatches ) {
        
        $best_pvalue = 1e-320 if($best_pvalue == 0);
        $berMatch->{'p_value'} = 1e-320 if($berMatch->{'p_value'} == 0);
 

        next unless ($best_pvalue && $berMatch->{'p_value'} && 
                log(-log($best_pvalue)/log(10))/log(10) - log(-log($berMatch->{'p_value'})/log(10))/log(10) < 0.5);

        #The compseq stored in the ber match have bad accessions (hacked to go inside bsml id
        #attributes (are these strict id requirements helping anyone?)).  So we must pull the real
        #acc and pass it to parseBerHeader so it knows the right accession to access.
        my $realAcc = $fakeId2RealId->{ $berMatch->{'compseq'} };
        
        print "\n$realAcc: Ber match passed pvalue test\n" if($logger->is_debug);

        #If we couldn't parse it die.  I don't know when this would happen, but we don't
        #want it, that's for sure.
        $logger->logdie("Couldn't find the real acc when given id ".$berMatch->{'compseq'})
            unless($realAcc);
        
        my $currentMatch = &getBerHeader( $realAcc );
        next if(@{$currentMatch} < 10);

        $currentMatch->[$name] = &add_annotation_level( $currentMatch->[$name], $berMatch->{'percent_id'},
                                                        $berMatch->{'subject_length'} );

        print "Current Name $currentMatch->[$name]\n" if($logger->is_debug);

        #Get Rank
        my $currentRank = &rankAnnotation($currentMatch);


        #Include ber match information in the ranking.
        if( $berMatch->{'percent_id'} < $per_id_cutoff || 
            $berMatch->{'fraction_length'} * 100 <=$length_cutoff) {
            $currentRank += 6;
        }

        print "Current rank = $currentRank\n" if($logger->is_debug);

        if($currentRank <= $bestBerRank && 
           &swapBestBerHit($bestBerMatch, $bestBerRank, $currentMatch, $currentRank) ) {
            print "Swapping the best hit (Current: $currentMatch->[$name], Best: $bestBerMatch->[$name]\n";
           
            ($bestBerMatch, $bestBerRank, $bestHitInfo) = ($currentMatch, $currentRank, $berMatch); 
        }

        
    }

    #Use strict cutoffs for assigning GO terms.
    unless($bestHitInfo->{'fraction_length'} >= 0.9 &&
           $bestHitInfo->{'percent_id'} >= 40) {  
        print "Getting rid of go_id\n" if($logger->is_debug);
        $bestBerMatch->[$go_id] = "";
    }

    $bestBerMatch->[$rank] = $bestBerRank;

    return $bestBerMatch;
    
    
}

sub getBerHeader {
    my $acc = shift;
    
    my $offset = $pandaHeaderOffsets{$acc};
    return [] unless($offset);
    my ($pos,$len) = split(/\s+/, $offset);
    open(HEAD, "< $headerFile") or $logger->logdie("Unable to open $headerFile ($!)");
    seek(HEAD, $pos, 0);
    my $tabs;
    read(HEAD, $tabs, $len);
    close(HEAD);

    my @cols = split(/\t/,$tabs);

    $logger->logdie("Wrong number of columns") if(@cols < 10);
    
    return \@cols;
}

sub add_annotation_level {
    my($com_name, $per_id, $length) = @_;
    return($com_name) if ($com_name =~ /hypothetical|putative|-related|homolog/i);
##
## want to add the words "putative" or "-related" to those com_names
## that warrant it. 
##
## "high confidence" if the match is greater than $per_id_cutoff spanning
## greater than $length_cutoff
##
## "putative" if the %identity is less than $per_id_cutoff (35) spanning
## greater than $length_cutoff (80).
##
## "-related" if the match is spanning less than $length_cutoff & over
## $per_id_cutoff.
##
## "low confidence" if the match is less than $per_id_cutoff spanning
## under $length_cutoff
##
    if (($per_id < $per_id_cutoff) && ($length >= $length_cutoff)) {
        $com_name .= ", putative";
    } elsif (($per_id >= $per_id_cutoff) && ($length < $length_cutoff)) {
        $com_name = "conserved domain protein";
    } elsif ($per_id < $per_id_cutoff && $length < $length_cutoff) {
        $com_name = "hypothetical protein";
    }
    return($com_name);
}

sub mapBsmlId2RealId {
    my ($twig, $seqElem, $fakeId2RealId) = @_;

    my $fakeId = $seqElem->att('id');
    $logger->logdie("Could not parse id out of sequence element") unless($fakeId);
    
    #Get the attribute element that holds the defline.
    foreach my $attr ( $seqElem->children('Attribute') ) {
        next unless($attr->att('name') eq 'defline');
        
        my $realId = $1 if($attr->att('content') =~ /^(\S+)/);
        $logger->logdie("Could not parse real id out of BSML defline attribute") 
            unless($realId);
        
        $fakeId2RealId->{$fakeId} = $realId;
    }

}

#Searches through a bsml file for feature elements
#with a class of polypeptide and parses out the ids.
#Returns a list of polypeptide ids.
sub getPolypeptideIds {
    my $inputFile = shift;
    my %polypeptides;

    #Now parse the file.
    my $twig = new XML::Twig( TwigHandlers => 
                              { 
                                  'Sequence' => \&getSeqId,
                                  'Feature' => sub { &getPolypeptideHandler(@_, \%polypeptides)},
                                  'Feature-group' => \&featureGroupHandler
                                  }
                              );  
    
    $twig->parsefile($inputFile);

    return %polypeptides;

}

sub getSeqId {
    my ($twig, $seqElem) = @_;
    $asmbl->{'id'} = $seqElem->att('id');
    die("Could not parse id") unless($asmbl->{'id'});

    $asmbl->{'topology'} = $seqElem->att('topology');

    foreach my $sdi ($seqElem->children('Seq-data-import')) {
        $asmbl->{'import'}->{'source'} = $sdi->att('source');
        $asmbl->{'import'}->{'ident'}  = $sdi->att('identifier');
        $asmbl->{'import'}->{'format'} = $sdi->att('id');
    }

    #Find the length of the sequence
    my $seq;
    open(IN, "< $asmbl->{'import'}->{'source'}") or 
        $logger->logdie("Unable to open  $asmbl->{'import'}->{'source'} ($!)");
    while(<IN>) {
        next if(/^>/);
        chomp;
        $seq .= $_;
    }
    close(IN);
    $asmbl->{'length'} = length($seq);
}

sub featureGroupHandler {
    my ($twig, $featGroupElem) = @_;

    my ($poly,$trans,$cds);
    foreach my $fgmElem ( $featGroupElem->children('Feature-group-member') ) {
        $poly  = $fgmElem->att('featref') if($fgmElem->att('feature-type') eq 'polypeptide');
        $trans = $fgmElem->att('featref') if($fgmElem->att('feature-type') eq 'transcript');
        $cds   = $fgmElem->att('featref') if($fgmElem->att('feature-type') eq 'CDS');
    }

    $ber_map->{$cds} = $poly;
    $ber_map->{$poly} = $trans;
}

#A twig handler method for parsing input bsml files from a
#gene prediction output bsml (glimmer3, but should be
#generic).  Parses out polypeptide elements, saves ids
#and gene lengths. 
# NOTE: Parses out start and stop interval loc positions
# of feature class=polypeptide and use it as the gene
# length.
sub getPolypeptideHandler {
    my ($twig, $featElem, $polypeptides) = @_;
    
    return unless( $featElem->att('class') eq 'polypeptide' );

    my $intLoc = $featElem->first_child('Interval-loc');
    $logger->logdie("Interval-loc element for polypeptide $featElem->att('id') does not exist".
                    " or is incomplete") unless($intLoc && $intLoc->att('startpos') &&
                                                $intLoc->att('endpos'));
    my ($id, $start, $end) = ($featElem->att('id'), $intLoc->att('startpos'), $intLoc->att('endpos') );

    $polypeptides->{ $id } =  $end - $start;

    $geneBoundaries->{ $id }->{'startpos'} = $start;
    $geneBoundaries->{ $id }->{'endpos'}   = $end;
    
    
}

#Takes in a list file or a directory and
#returns a list of analysis files.  Will search for
#any bsml files in the directory.
sub getAnalysisFiles {
    my $analysis = shift;
    my @retval;
    if( -d $analysis ) {
        find( sub { &findAnalysisFiles(\@retval) }, $analysis );
    } elsif( -e $analysis ) {
        open(ANFILES, "< $analysis") or
            $logger->logdie("Could not open $analysis ($!)");
        chomp(@retval = <ANFILES>);
        close(ANFILES);
    }
    
    return @retval;
    
}

#The File::Find helper function for sub getAnalysisFiles.
sub findAnalysisFiles {
    my $analysisFiles = shift;
    if($File::Find::name =~ /.*bsml/) {
        push(@$analysisFiles, $File::Find::name);
    }

}

sub getHMMEvidence {
    my $polyId = shift;
    my $hmmBsmlFile;
    my @hmmMatches;

    #Find the associated file.
    foreach my $bFile ( @hmmBsmlFiles ) {
        if($bFile =~ /($polyId)/) {
            $hmmBsmlFile = $bFile;
            last;
        }
    }

    #Now parse the file.
    my $twig = new XML::Twig( TwigHandlers => 
                              { 'Seq-pair-alignment' => sub { &handleHMMSpa(@_, \@hmmMatches)} }         
                              ); 
    print "\tparsing hmmFile: $hmmBsmlFile\n" if($logger->is_debug);;
    $twig->parsefile($hmmBsmlFile);

    return {} if(@hmmMatches == 0);

    #Now that we have all the hmmMatches in this array, lets find the best one.
    my $bestMatch = &getBestHMMMatch( $polyId, @hmmMatches );

    unless($bestMatch) {
        $logger->warn("Couldn't find best match");
    }

    #print "\treturning best hmm match\n" if($logger->is_debug);
    return $bestMatch;
   
}

#From getHMMEvidence.  Twig Handler
#spa = Sequence Pair Alignment
#spr = Seq-pair-run
sub handleHMMSpa {
    my ($twig, $spa, $hmmMatches) = @_;
    my $match;      #Hash ref containing the match

    #Store the accession
    $match->{'hmm_acc'} = $spa->att('compseq');
    
    #Retrieve the total score of the SPA
    my @attElems = $spa->children('Attribute');
    $logger->logdie("spa didn't have attribute elements: $match->{hmm_acc}")
        unless( @attElems > 0 );
    foreach my $attElem ( @attElems ) {
        if( $attElem->att('name') eq 'total_score' ) {
            $match->{'total_score'} =$attElem->att('content');
            last;
        }
    }
    $logger->logdie("Unable to parse total_score information from spa $match->{'hmm_acc'}")
        unless($match->{'total_score'});

    #Get the coordinate information

    #Retrieve the domain scores
    foreach my $spr ( $spa->children("Seq-pair-run") ) {
        my $sprHash = { 
            'domain_score' => $spr->att('runscore'),
            'refpos'       => $spr->att('refpos'),
            'runlength'    => $spr->att('runlength'),
            'comprunlength'=> $spr->att('comprunlength'),
        };

        $logger->logdie("Unable to parse domain score from spr: $match->{'hmm_acc'}")
            unless($sprHash->{'domain_score'});
        push(@{$match->{'domains'}}, $sprHash);                        
    }
    
    #Combines the info stored in the hmm_info_db lookup with the
    #information just parsed from the bsml file.
    my $tmpHash;
    foreach my $hrToCombine ($match, $hmmInfo{$match->{'hmm_acc'}} ) {
        while(my ($k, $v) = each(%$hrToCombine)) {
            $tmpHash->{$k} = $v;
        }
    }

    push(@{$hmmMatches}, $tmpHash);    
    
}

sub getBestHMMMatch {
    my $polyId = shift;
    my @hmmMatches = @_;
    my $winner = {};
    my $beat_trusted = 0;
    
    my $count = 0;

    foreach my $hmmMatch ( @hmmMatches ) {
        foreach my $domain ( @{$hmmMatch->{'domains'}} ) {

            #Make sure that the score hits cutoff:
            #
            # In cases of true equivalogs, the model should only hit the subject protein once.
            # Therefore, domain score should equal total score. There is at least one case
            # of an equivalog model that hits members of another family of proteins two times.
            # The domain scores are lower than trusted, but the total score is higher than 
            # trusted. For this reason, in the case of equivalog hits, domain score should be
            # compared to trusted and not total score. Still don't understand? Talk to Bill Nelson
            # or Dan Haft.
            unless($hmmMatch->{'iso_type'}) {
                print "Doesn't have iso_type\n";
                print Dumper($hmmMatch);
                next;
            }

            unless($hmmMatch->{'trusted_cutoff'}) {
                print "Doesn't have trusted_cutoff\n";
                print Dumper($hmmMatch);
                next;
            }

            unless($domain->{'domain_score'}) {
                print "Doesn't have domain_score\n";
                print Dumper($domain);
                next;
            }

            if($hmmMatch->{'iso_type'} =~ /equivalog/ && $domain->{'domain_score'} >= $hmmMatch->{'trusted_cutoff'}) {
                $winner = $hmmMatch unless($winner);
            } elsif( $hmmMatch->{'total_score'} >= $hmmMatch->{'trusted_cutoff'} ) {
                $winner = $hmmMatch unless($winner);
            } else {
                next;
            }

            print "\tHMM $hmmMatch->{hmm_acc} passed cutoff\n" if($logger->is_debug);

            # test length of HMM hit. If it's longer than 80% of the protein length
            # then we call it a full-length hit, otherwise partial
            # 'full' means full protein hit regardless of hmm hit
            # 'partial' means full hmm hit, partial protein hit
            # 'partial2' means partial hmm hit, partial protien hit
            my ($protLength, $hitType) = ($polypeptideIds{$polyId}/3, "");
            my $protHitLength = (($domain->{'runlength'}/$protLength) * 100);
            my $hmmHitLength  = (($domain->{'comprunlength'}/$hmmMatch->{'hmm_len'}) * 100);
            if ($protHitLength > 70) {
                $hitType = "full";
            } else {
                $hitType = "partial";
            }
            $hitType = "partial2" if ($hitType eq "partial" && $hmmHitLength < 80 );

            #Okay, to determine what the rank is, check autoAnnotate.data.  It holds a fun data structure.
            if(!$winner->{'rank'} || 
               (defined( $isoType{ $hmmMatch->{'iso_type'} }->{ $hitType } ) &&
               $isoType{ $hmmMatch->{'iso_type'} }->{$hitType} < $winner->{'rank'} )) {
                $winner = $hmmMatch;
                $winner->{'rank'} = $isoType{ $hmmMatch->{'iso_type'} }->{$hitType};
                print "\tfound a winning hmm rank ($winner->{'rank'})\n" if($logger->is_debug);
            }
                                   
        }

    }

  
    return $winner;
}

sub annotation2bsml {
    my $output = shift;
    my $doc = new BSML::BsmlBuilder();
    
    #Create the asmbl sequence
    my $seq = $doc->createAndAddSequence($asmbl->{'id'},$asmbl->{'id'},'','dna','assembly');
    my $sdi = $doc->createAndAddSeqDataImport($seq, $asmbl->{'import'}->{'format'}, 
                                              $asmbl->{'import'}->{'source'}, '',
                                              $asmbl->{'import'}->{'identifier'});
    my $link = $doc->createAndAddLink($seq, 'analysis', '#autoAnnotate_analysis', 'input_of');
    my $featTable = $doc->createAndAddFeatureTable($seq);

    foreach my $polyid (keys %{$finalAnnotes}) {
        my $feat = $doc->createAndAddFeature($featTable, $ber_map->{$polyid}, $ber_map->{$polyid}, 'transcript');
        $feat->addBsmlLink('analysis', '#autoAnnotate_analysis', 'computed_by');
        
        foreach my $term (keys %{$finalAnnotes->{$polyid}}) {
            next unless($finalAnnotes->{$polyid}->{$term});

            if($term =~ /(go_id|ec_num|role_id)/) {
                my @vals = split(/\s+/, $finalAnnotes->{$polyid}->{$term});
                foreach my $val(@vals) {
                    $doc->createAndAddBsmlAttribute( $feat, $term, $val );
                }
            } else {
                $doc->createAndAddBsmlAttribute( $feat, $term, $finalAnnotes->{$polyid}->{$term} );
            }
        }
        
    }

    $doc->createAndAddAnalysis( 'id' => 'autoAnnotate_analysis',
                                'sourcename' => $inputFiles[0] );

    print "Writing to $output\n";
    $doc->write($output);

    
}


sub _pod {   
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
