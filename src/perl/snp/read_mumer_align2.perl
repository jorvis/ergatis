#!/usr/bin/perl

# (c) Tim Read July 2000 - all standard disclaimers apply - use at your own risk!
#
#MUMmer can be obtain through license from the TIGR website (www.tigr.org)
#If you have questions about this script email me: tread@tigr.org
#
#usage read_mumer_align.perl <*.align> 
#---> where *.align is output from MUMmer comparison of aligned genomes
#output to STDOUT:
#coord genome1     base genome 1    base genome2      coord genome2
#NOTE: I do not recommend using this script when the MUM parameter in the original
# alignment is set below 10.  This will lead occasionally to double-reporting of
#mutations. (Default parameter is usually 20.)

open (ALIGNFILE, $ARGV[0]);

while (<ALIGNFILE>) {
    @entry = split(/\s+/,$_);

#only reades the first set of alignments (ie not the non-sequential hits)
    if (/^>/ && $linemarker == 0) {$linemarker = 1;}
    elsif (/^>/ && $linemarker == 1) {$linemarker = 0;}

    if ($linemarker == 1) {
	
       #Loads coords at start of alignment
	if ($entry[1] =~ m/\d+/){ 
	    $fasta1coord = $entry[1] - $entry[5];
	    $fasta2coord = $entry[2] - $entry[6];
	   
	    #initializer lets the program know when to start counting the coords
	    $initializer = 0;
	}
	#Loads numbers of errors
	if ($entry[1] =~ /Errors/) {
	    $errors = $entry[3];
	    $TOTerrors = $TOTerrors + $errors;
	}


	#Loads alignment strings
	if ($entry[0] eq "A:") {
	    $Tseq = $entry[1];
	}

	if ($entry[0] eq "B:") {
	    $Sseq = $entry[1];
	}

	#Loads string with position of mutations
	if (/\^/g) {
	    $mutations = $_;


#This is where the action happens ! #####################################
	    while ($mutations =~ /\^/g) {
	  
		$posit = (pos($mutations) - 4);
		if ($initializer == 0) {$start = $posit;} #$start will be zero unless it is the first line of the alignment
		$initializer++;
	        $Tsub = substr($Tseq,$start,($posit - $start));
		$Tcoord = &realsize($Tsub) + $fasta1coord;
		
		$Ssub = substr($Sseq,$start,($posit - $start));
		$Scoord = &realsize($Ssub) + $fasta2coord;
		print $Tcoord."\t".substr($Tseq,$posit,1)."\t".substr($Sseq,$posit,1)."\t".$Scoord."\n";
	    }
###########################################################################
#update genome coords here
	  $fasta1coord = $fasta1coord + &realsize(substr($Tseq,$start,));
	  $fasta2coord = $fasta2coord + &realsize(substr($Sseq,$start,));
	  $start = 0; # because we have already been once through, start is now at the beginning of the alignment
	}#mutaion bracket
	
	

    }# linemarker bracket

}# whole loop bracket

    print "Total number of errors = ".$TOTerrors."\n";


sub realsize {
# finds the length of the substring without spaces introduced by alignmnet program
   local $temp = "";
   $temp = $_[0];
   $temp =~ s/\-//g;
   return length($temp);
}
    
