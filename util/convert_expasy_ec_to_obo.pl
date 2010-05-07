#!/usr/local/bin/perl -w
#-------------------------------------------------------------------------------
#
# script: ecMaker.pl 
# author: sundaram@tigr.org
# date:   2004-09-22
# cvs:    ANNOTATION/chado/util/EnzymeCommission/ecMaker.pl
# 
# purpose: Reads in Enzyme Commission (ec) category file and detailed file
#          and produces ec.obo file.
#
#
# resources:
#
#          Category file: ftp://au.expasy.org/databases/enzyme/release/enzyme.dat
#          Detailed file: ftp://au.expasy.org/databases/enzyme/release/enzclass.dat
#
#          Category file was saved as ANNOTATION/chado/util/EnzymeCommission/enzyme.dat
#          Detailed file was saved as ANNOTATION/chado/util/EnzymeCommission/enclass.dat
#
#          enzyme.dat had to be edited (some upper-level sub-category descriptors were
#          missing) using to the information stored in legacy database: egad..prot_function
#          Edited version: enzyme.edited.dat
#          
#
# invocation:
#
#          ./ecMaker.pl --infile1=enzyme.edited.dat --infile2=enclass.dat --outfile-ec.obo
#
#
#-------------------------------------------------------------------------------
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use FileHandle;
use Log::Log4perl qw(get_logger);


#$/ = "||||";
#my $newline = $/;
#die "newline is '$newline'";


my($infile, $outfile, $help, $man, $log4perl, $verbose, $prefix, $infile2);


&GetOptions("infile|i=s"   => \$infile,
	    "infile2|f=s"  => \$infile2,
	    "outfile|o=s"  => \$outfile,
	    "log4perl|l=s" => \$log4perl,
	    "help"         => \$help,
	    "man"          => \$man,
	    );

pod2usage(1) if $help;
pod2usage({-verbose => 2}) if $man;
pod2usage({-message => "Error:\n     --infile must be specified\n", -exitstatus => 0, -verbose => 0}) if (!$infile);
pod2usage({-message => "Error:\n     --infile2 must be specified\n", -exitstatus => 0, -verbose => 0}) if (!$infile2);
#pod2usage({-message => "Error:\n     --outfile must be specified\n", -exitstatus => 0, -verbose => 0}) if (!$outfile);




#
#  Default log4perl logfile assignment.
#
if (!defined($log4perl)){
    $log4perl = '/tmp/ecMaker.pl.log';
    print "log4perl was not defined, therefore was set to '$log4perl'\n";
}


my $screen_threshold = 'WARN';
if ($verbose){
    $screen_threshold = 'INFO';
}

Log::Log4perl->init(
		    \ qq{
			log4perl.logger                       = INFO, A1, Screen
			log4perl.appender.A1                  = Log::Dispatch::File
			log4perl.appender.A1.filename         = $log4perl
			log4perl.appender.A1.mode             = write
			log4perl.appender.A1.Threshold        = INFO
			log4perl.appender.A1.layout           = Log::Log4perl::Layout::PatternLayout
			log4perl.appender.A1.layout.ConversionPattern = %d %p> %F{1}:%L %M - %m%n 
			log4perl.appender.Screen              = Log::Dispatch::Screen
			log4perl.appender.Screen.layout       = Log::Log4perl::Layout::SimpleLayout
                        #log4perl.appender.Screen.layout.ConversionPattern =%d %p> %F{1}:%L %M - %m%n 
			log4perl.appender.Screen.Threshold    = $screen_threshold
			Log::Log4perl::SimpleLayout
		    }
		    );


my $logger = get_logger("ecMaker.pl");

# --------------------------------------------------------------
# Main program
# --------------------------------------------------------------

$| = 1;


if (!defined($outfile)){
    $outfile = "ec.obo";
    $logger->debug("outfile was not defined, therefore was set to '$outfile'") if $logger->is_debug();
}


if (!defined($prefix)){
    $prefix = 'EC';
    $logger->info("prefix was not defined, therefore was set to '$prefix'");
}



$logger->debug("Opening infile '$infile' and outfile '$outfile'") if $logger->is_debug();

open (INFILE, "<$infile") or $logger->logdie("Could not open infile '$infile'");
open (INFILEE, "<$infile2") or $logger->logdie("Could not open infile '$infile2'");
open (OUTFILE, ">$outfile") or $logger->logdie("Could not open outfile '$outfile'");


my @contents = <INFILE>;
chomp @contents;
close INFILE or $logger->logdie("Could not close infile '$infile'");


my @contents2 = <INFILEE>;
chomp @contents2;
close INFILEE or $logger->logdie("Could not close infile2 '$infile2'");


my $pretree = &pre_build_tree(\@contents);

#print Dumper $pretree;die;

$logger->logdie("pretree was not defined") if (!defined($pretree));


my $forest = &process_detail_file(\@contents2);

#print Dumper $forest;die;

my $merged_trees = &merge_trees($forest, $pretree);


#print Dumper $merged_trees;die;

my $tree = &build_tree($merged_trees);

#print Dumper $tree;die;

&build_obo_file($tree, $outfile);



#--------------------------------------------------------------------------------
# build_obo_file()
#
# This function should write each node
# in the following manner:
#
#
# [Term]
# id: EC:0000001
# name: OXIDOREDUCTASE
# xref_analog: EC:1.-.-.-
#
# [Term]
# id: EC:0000002
# name: ACTING ON THE CH-OH GROUP OF DONOR
# xref_analog: EC:1.1.-.-
# is_a: EC:0000001
#
# [Term]
# id: EC:0000003
# name: WITH NAD(+) OR NADP(+) AS ACCEPTOR
# xref_analog: EC:1.1.1.-
# is_a: EC:0000002
#
# [Term]
# id: EC:0000004
# name: Alcohol dehydrogenase
# def: An alcohol + NAD(+) = an aldehyde or ketone + NADH
# xref_analog: EC:1.1.1.1
# is_a: EC:0000003
#
#-----------------------------------------------------------------------------------
sub build_obo_file {

    my $hash = shift;
    my $outfile = shift;

#    print Dumper $hash; die;


    open (OUTFILE, ">$outfile") or $logger->logdie("Could not open output file '$outfile'");
    


    &write_header();


    my $keycount = 0;
    my $compartment;

    #
    # Lookup hashes
    #
    my $ec2id = {};
    my $id2ec = {};



    $logger->debug("Building the lookup hashes...") if $logger->is_debug();

    foreach my $ec (keys %{$hash}){
	
	$keycount++;

	my $ident = &make_ident($keycount, 7);


	$ec2id->{$ec} = $ident;
	$id2ec->{$ident} = $ec;

    }

    
    #
    # Re-initialize
    #
    $keycount=0;


    $logger->debug("About to start outputting terms to OBO file...") if $logger->is_debug();

    my $unique_names = {};
    


    foreach my $ec (keys %{$hash}){
	
	$keycount++;

	my $ident = &make_ident($keycount, 7);

	$logger->logdie("id was not defined for ec '$ec'") if (!defined($ident));

	my $is_a        = $hash->{$ec}->{'is_a'}; 
	my $name        = $hash->{$ec}->{'full_name'}; 
	my $xref_analog = $ec;
	my $def         = $hash->{$ec}->{'def'};
	my $comment     = $hash->{$ec}->{'comment'};
	my $is_obsolete = $hash->{$ec}->{'is_obsolete'};
	my $synonym     = $hash->{$ec}->{'synonym'};


	$logger->logdie("name was not defined for id '$ident' ec '$ec'") if (!defined($name));

	$name =~ s/\.$//;


	if ((exists $unique_names->{$name}) and (defined($unique_names->{$name}))){

	    #
	    # In order to make this term name truly unique, will suffix with the ec#
	    #

	    if ((defined($is_obsolete)) and ($is_obsolete == 1)){
		$logger->warn("Found obsolete version of a duplicate term name '$name'");
		$unique_names->{$name}++;
	    }
	    else{


		$logger->warn("Had to incorporate ec '$ec' into name '$name' to create truly unique term name'");
		$name .= '||' . $ec;
		$unique_names->{$name} = $name;
	    }
	    
# 	    if ((exists($hash->{$ec}->{'an'})) and (defined($hash->{$ec}->{'an'}))){
		
# 		$name .= '||' . $hash->{$ec}->{'an'};
# 		$name =~ s/\.$//;
		
# 		if ((exists $unique_names->{$name}) and (defined($unique_names->{$name}))){
		    
# 		    $logger->logdie("Extended unique name already exists.");
# 		}
# 		else{
# 		    $logger->warn("Had to incorporate AN line into name '$name' for ec '$ec'");
# 		    $unique_names->{$name} = $name;
# 		}
# 	    }
# 	    else{

# 		if ((defined($is_obsolete)) and ($is_obsolete == 1)){
# 		    $unique_names->{$name}++;
# 		}
# 		else{
# 		    $logger->logdie("The situation is hopeless.  AN line does not exist for ec '$ec' name '$name' is_obsolete '$is_obsolete'.  unique_names:");# . Dumper $unique_names);
# 		}
# 	    }


	}
	else{
	    $unique_names->{$name} = $name;
	}
   


	my $node = "[Term]\n".
	"id: $ident\n".
	"name: $name\n";
	

	if (defined($def)){
	    $def =~ s/\.$//;
	    $node .= "def: \"$def\"\n"; 
	}


	$node .= "comment: \"$comment\"\n" if (defined($comment));
	$node .= "is_obsolete: true\n" if (defined($is_obsolete));
	$node .= "synonym: \"$synonym\" []\n" if (defined($synonym));

	if (defined($xref_analog)){
	    $node .= "xref_analog: " . $prefix . ':' . $xref_analog ."\n";
	}

	if (defined($is_a)){
	    
	    if ((exists ($ec2id->{$is_a})) and (defined($ec2id->{$is_a}))) {

		$node .= "is_a: " . $ec2id->{$is_a} . "\n";

	    }
	    else{
		$logger->logdie("Could not find an EC OBO identifier for the ec# '$is_a'. ec2id:" . Dumper $ec2id);
	    }
	}


#	if ( (defined($def)) and ($def =~ /OBSOLETE/)) {
#	    $node .= "is_obsolete: true\n";
#	}

	print OUTFILE $node . "\n";
	
    }

#    print "ec2id:" . Dumper $ec2id;
#    print "id2ec:" . Dumper $id2ec;




}









#----------------------------------------------------
# process_detail_file()
#
#----------------------------------------------------
sub process_detail_file {
    
    my $contents = shift;
    my $tree = {}; # yet another tree...

    #
    # Want to hear a joke?
    # What did the brother tree say to the sister tree?...
    #


    my $id_encountered = 0;

    my $linectr=0;
    my $recctr=0;


    my $id;
    my $name;
    my $comment;
    my $def;
    my $xref_analog;
    my $an;


    foreach my $line (@{$contents}) {
	
	$linectr++;
	
	if ($line =~ /^\/\//){

	    $logger->debug("record $recctr") if $logger->is_debug();
	    
	    if ($id_encountered == 1){
		

		
		#
		# Un-set the flag as next set of lines belong to the next ID
		#
		$id_encountered = 0;
#	    $tree->{$id}->{'name'} = $name;
#	    $tree->{$id}->{'def'} = $def;
#	    $tree->{$id}->{'xref_analog'} = $xref_analog;
#	    $tree->{$id}->{'comment'} = $comment;
		
		$name = undef;
		$def = undef;
		$xref_analog = undef;
		$comment = undef;
		$an = undef;
	    }
	}
	elsif ($line =~ /^ID/){

#	    $logger->warn("ID $linectr");

	    
	    if ($id_encountered == 0){
		$id_encountered = 1;
		$recctr++;
	    }
	    else{
		$logger->logdie("Corrupted file? Encountered an ID at line $linectr '$line', but was not expecting one...");
	    }
	    
	    if ($line =~ /^ID\s+(\S+)/){
		$id = $1;
		$tree->{$id}->{'full_ec'} = $id;
		

		$logger->info("Found $id");


		my ($ec1, $ec2, $ec3, $ec4) = split(/\./, $id);

		$tree->{$id}->{'e1'} = $ec1;
		$tree->{$id}->{'e2'} = $ec2;
		$tree->{$id}->{'e3'} = $ec3;
		$tree->{$id}->{'e4'} = $ec4;




	    }
	    else{
		$logger->logdie("Could not parse ID at line $linectr '$line'");
	    }
	}
 	elsif ($line =~ /^CC/){
	    

#	    $logger->warn("CC $linectr");


	    if ($id_encountered == 0 ){
		#
		# This is a 
		#
		next;
	    }
	    else{
		#
		# This is a useful comment line
		#
		if (   $line =~ /^CC\s+([\S\s]+)/  ){
		    $comment .= $1;
		    $tree->{$id}->{'comment'} .= $comment;
		}
		else{
		    $logger->logdie("Could not parse CC at line $linectr '$line'");
		}
	    }
	}
 	elsif ($line =~ /^AN/){
	    

#	    $logger->warn("AN $linectr");


	    if ($id_encountered == 0 ){
		#
		# This is a 
		#
		$logger->logdie("Encountered an AN at line $linectr '$line', but did not yet encounter an ID line...");
	    }
	    else{
		#
		# This is a useful comment line
		#
		if (   $line =~ /^AN\s+([\S\s]+)/  ){
		
		    #
		    # Strip all []
		    #
		    my $pre_an = $1;
		    $pre_an =~ s/\[//;
		    $pre_an =~ s/\]//;
		    
		    $an .= $pre_an;
		    $tree->{$id}->{'synonym'} .= $pre_an;
		}
		else{
		    $logger->logdie("Could not parse AN at line $linectr '$line'");
		}
	    }
	}
	elsif ($line =~ /^DE/ ){
	    

#	    $logger->fatal("DE $linectr");


	    if ($id_encountered == 0 ){
		#
		# This is a 
		#
		$logger->logdie("Encountered a DE at line $linectr '$line', but did not yet encounter an ID line...");
	    }
	    else{
		#
		# 
		#
		if ($line =~ /^DE\s+([\S\s]+)/   ){
		    $name = $1;
		    $tree->{$id}->{'partial_name'} .= $name;


		    if ($name =~ /Transferred entry/){

			#
			# Transferred entry
			# Will look like:
			#
			# //
			# ID 1.1.1.249
			# DE Transferred entry: 2.5.1.46.
			# //
			#

			$logger->fatal("Transferred entry discovered at record '$recctr' line '$linectr': $name") if $logger->is_debug();

			$tree->{$id}->{'is_obsolete'} = 1;
			$tree->{$id}->{'def'} .= $name;
		    }
		    elsif ($name =~ /Deleted entry/){


			#
			# Deleted entry
			# Will look like:
			#
			# //
			# ID 1.1.1.249
			# DE Deleted entry.
			# //
			#

			$logger->fatal("Deleted entry discovered at record '$recctr' line '$linectr': $name") if $logger->is_debug();
			
			$tree->{$id}->{'is_obsolete'} = 1;
			$tree->{$id}->{'def'} .= $name;
			
		    }
		}
		else{
		    $logger->logdie("Could not parse DE line $linectr '$line'");
		}
	    }
	}
	elsif($line =~ /^CA/){
	    

#	    $logger->warn("CA $linectr");


	    if ($id_encountered == 0 ){
		#
		#
		#
		$logger->logdie("Encountered a CA line $linectr '$line', but did not yet encounter an ID line...");
	    }
	    else{
		#
		#
		#
		if ($line =~ /^CA\s+([\S\s]+)/){
		    $def .= $1;
		    $tree->{$id}->{'def'} .= $def;
		}
		else{
		    $logger->logdie("Could not parse CA line $linectr '$line'");
		    
		}
	    }
	}
	elsif($line =~ /^PR/){
	    

#	    $logger->warn("PR $linectr");


	    if ($id_encountered == 0 ){
		#
		#
		#
		$logger->logdie("Encountered a PR line $linectr '$line', but did not yet encounter an ID line...");
	    }
	    else{
		#
		#
		#
		if ($line =~ /^PR\s+([\S\s]+)/){
		    $xref_analog .= $1;
		    $tree->{$id}->{'xref_analog'} .= $xref_analog;
		}
		else{
		    $logger->logdie("Could not parse PR line $linectr '$line'");
		    
		}
	    }
	}
	else{
	    if ($line =~ /^(\S+)\s+/){
		$logger->debug("Ignoring this type of line '$1'")# if $logger->is_debug();
	    }
	    else{
		$logger->fatal("Ignoring some kind of line at line $linectr") if $logger->is_debug();
	    }
	}
    }


#    print Dumper $tree;die;

    $logger->info("Found $recctr nodes in the detail file");
#    print "total lines $linectr. Second file:"  . Dumper $tree;die;
    return $tree;
}


#----------------------------------------------------
# build_tree()
#
#----------------------------------------------------
sub build_tree {

    my $pretree =  shift;

    my $tree = {};

    my $missingec=0;
    my $missing = {};

    #
    # Establish relationships among the nodes.
    #

    foreach my $key ( sort keys %{$pretree} ){
	
#	print $key ."\n";
#	next;

	my $parent;


	my $ec2 = $pretree->{$key}->{'e2'};
	if ($ec2 eq '-'){

	    #
	    # This node has no parent.
	    #
#	    $parent = $pretree->{$key}->{'e1'} . '.-.-.-';

	    #
	    # Assign this node to the tree
	    #
	   # $tree->{$key} = $pretree->{$key};

#	    $tree->{$parent}->{'node'} = $pretree->{$parent};
#	    push ( @{$tree->{$parent}->{'children'}}, $pretree->{$key});

	    $pretree->{$key}->{'full_name'} = $pretree->{$key}->{'partial_name'};


	}
	else{
	    my $ec3 = $pretree->{$key}->{'e3'};
	    if ($ec3 eq '-' ) {
		#
		#
		#
		my $g1 = $pretree->{$key}->{'e1'} . '.-.-.-';

		if ((exists ($pretree->{$g1})) and (defined($pretree->{$g1}))){
		    
		    my $fullname = $pretree->{$g1}->{'partial_name'} . '||' . 
		    $pretree->{$key}->{'partial_name'};
		    
#		$tree->{$key}->{'is_a'} = $pretree->{$g1}->{'full_ec'};
		    
		    $pretree->{$key}->{'is_a'} = $pretree->{$g1}->{'full_ec'};
		    $pretree->{$key}->{'full_name'} = $fullname;
		}
		else{
		    $logger->fatal("g1 '$g1' does not exist in pretree!");
		    $missingec++;
		    $missing->{$g1}++;
		    next;
		}
	    }
	    else{
		my $ec4 = $pretree->{$key}->{'e4'};
		if ($ec4 eq '-' ) {
		    #
		    #
		    #

		    my $g2 = $pretree->{$key}->{'e1'} . '.' . $pretree->{$key}->{'e2'} . '.-.-';


		    if (( exists($pretree->{$g2})) and (defined($pretree->{$g2}))){

			my $g1 = $pretree->{$g2}->{'is_a'};
		    
			if (( exists($pretree->{$g1})) and (defined($pretree->{$g1}))){


			    my $fullname = $pretree->{$g1}->{'partial_name'} . '||' .
			    $pretree->{$g2}->{'partial_name'} . '||' .
			    $pretree->{$key}->{'partial_name'};

#		    $tree->{$key}->{'is_a'} = $pretree->{$parent}->{'full_ec'};
			    
			    $pretree->{$key}->{'is_a'} = $pretree->{$g2}->{'full_ec'};
			    
			    $pretree->{$key}->{'full_name'} = $fullname;
			}
			else{
			    $logger->fatal("g1 '$g1' does not exist in pretree!");
			    $missing->{$g1}++;
			    $missingec++;
			    next;
			}
		    }
		    else{
			$logger->fatal("g2 '$g2' does not exist in pretree!");
			$missing->{$g2}++;
			$missingec++;
			next;
		    }
		}
		else{
		    #
		    # This is a leaf node
		    #
		    my $g3 = $pretree->{$key}->{'e1'} . '.' . $pretree->{$key}->{'e2'} . '.'. $pretree->{$key}->{'e3'} . '.-';
	

		    if ( (exists($pretree->{$g3})) and (defined($pretree->{$g3}))){ 

			my $g2 = $pretree->{$g3}->{'is_a'};

			if ( (exists($pretree->{$g2})) and (defined($pretree->{$g2}))){ 
			

			    my $g1 = $pretree->{$g2}->{'is_a'};

			    if ( (exists($pretree->{$g1})) and (defined($pretree->{$g1}))){ 

				my $fullname = $pretree->{$g1}->{'partial_name'} . '||' .  
				$pretree->{$g2}->{'partial_name'} . '||' . 
				$pretree->{$g3}->{'partial_name'} . '||' . 
				$pretree->{$key}->{'partial_name'};
				
#		    $tree->{$key}->{'is_a'} = $pretree->{$parent}->{'full_ec'};
				$pretree->{$key}->{'is_a'} = $pretree->{$g3}->{'full_ec'};
				
				$pretree->{$key}->{'full_name'} = $fullname;
				
#		    die "g1 '$g1' g2 '$g2' g3 '$g3' g4 '$key' fullname '$fullname'";
			    }
			    else{
				$logger->fatal("g1 '$g1' does not exist in pretree!");
				$missingec++;
				$missing->{$g1}++;
				next;
			    }
			}
			else{
			    $logger->fatal("g2 '$g2' does not exist in pretree!");
			    $missingec++;
			    $missing->{$g2}++;
			    next;
			}
		    }
		    else{
			$logger->fatal("g3 '$g3' does not exist in pretree!");
			$missingec++;
			$missing->{$g3}++;
			next;
		    }
		}
	    }
	}
    }


#    print Dumper $pretree;die;


    if ($missingec > 0) {
	$logger->logdie("Upper-level ec subfunctions where missing.  Counted $missingec missing inner-nodes.  Please review source file '$infile' and detailed source file '$infile2'  This script $0 found inconsistencies." . Dumper $missing);
    }




#     #
#     # Build the fully qualified name
#     #
#     foreach my $node (sort keys %{$pretree} ) {
	
# 	if ((exists $pretree->{$node}->{'is_a'} ) and (defined($pretree->{$node}->{'is_a'}))) {

# 	    my $parent = $pretree->{$node}->{'is_a'};
# 	    my $partial_full_name = $pretree->{$node}->{'partial_full_name'};


# 	    if ((exists $pretree->{$parent}->{'is_a'} ) and (defined($pretree->{$parent}->{'is_a'}))) {
		
# 		my $grandparent = $pretree->{$parent}->{'is_a'};
		
		
# 		my $fullname = $pretree->{$grandparent}->{'partial_name'} . '||'.  $partial_full_name;
# 		$pretree->{$node}->{'full_full_name'} = $fullname;
		
# 	    }
# 	    else{
# 		$logger->debug("No grandparent for node '$node'") if $logger->is_debug();
# 	    }
# 	}
#     }

#    print Dumper $pretree;die;
    return $pretree;
}


#----------------------------------------------------
# merge_trees()
#
#----------------------------------------------------
sub merge_trees {
    
    my ($forest, $pretree) = @_;


    my $ctr=0;

    foreach my $id (sort keys %{$forest} ) {

	if ( (exists ($pretree->{$id})) and (defined($pretree->{$id}))) {
	    $logger->logdie("$id should not already exist in the pretree.");
	}
	$pretree->{$id} = $forest->{$id};
	$ctr++;
    }
    
#    print Dumper $pretree;die;


    $logger->info("Merged $ctr nodes from forest to the pretree");
    return $pretree;
}


#----------------------------------------------------
# pre_build_tree()
#
#----------------------------------------------------
sub pre_build_tree {

    #
    # In this section, we are parsing the raw contents of the source file.
    #

    my $contents = shift;
    my $tree = {};
    my $linectr = 0;
    my $ctr=0;


    foreach my $line ( @{$contents} ) {

	$linectr++;


#	die "$line";

	if ($line !~ /^\d/){
	    $logger->debug("Ignoring this line $linectr: '$line'") if $logger->is_debug();
	    next;
	}

	if ($line =~ /^([\d\s]+)\.([\d\s\-]+)\.([\d\s\-]+)\.([\d\s\-]+)\s+([\S\s]+)\./){
	    my ($ec1, $ec2, $ec3, $ec4, $name) = ($1, $2, $3, $4, $5);

	    $ec1 =~ s/\s+//g;
	    $ec2 =~ s/\s+//g;
	    $ec3 =~ s/\s+//g;
	    $ec4 =~ s/\s+//g;

	    $logger->logdie("ec1 '$ec1' line '$line'") if (length($ec1) > 2);
	    $logger->logdie("ec2 '$ec2' line '$line'") if (length($ec2) > 2);
	    $logger->logdie("ec3 '$ec3' line '$line'") if (length($ec3) > 2);
	    $logger->logdie("ec4 '$ec4' line '$line'") if (length($ec4) > 2);


	    $logger->info("ec1 '$ec1' ec2 '$ec2' ec3 '$ec3' ec4 '$ec4' name '$name'");


	    $ctr++;

	    my $qualified_ec = $ec1 . '.' . $ec2 . '.' . $ec3 . '.' . $ec4;
	    
	    $tree->{$qualified_ec}->{'full_ec'} = $qualified_ec;
	    $tree->{$qualified_ec}->{'e1'} = $ec1;
	    $tree->{$qualified_ec}->{'e2'} = $ec2;
	    $tree->{$qualified_ec}->{'e3'} = $ec3;
	    $tree->{$qualified_ec}->{'e4'} = $ec4;
	    $tree->{$qualified_ec}->{'partial_name'} = $name;
#	    $tree->{$qualified_ec}->{'children'} = [];

	}
	else{
	    $logger->logdie("Could not parse line $linectr: '$line'");
	}
    }
#    print Dumper $tree;die;

    $logger->info("Build pretree with $ctr nodes");
    return $tree;
}



#---------------------------------------------------------------------------------
# process_egad_prot_function()
#
#---------------------------------------------------------------------------------
sub process_egad_prot_function {


    my $hash = {};
    my $ctr = {};
    my $name_hash = {};

#
#  Alright, go ahead and iterate over all of the egad.prot_function records.
#
    foreach my $line (@contents){

	my $deleted_entry = 0;
	my $transferred_entry = 0;

	my ($ec, $pfunc, $sub1, $sub2, $sub3, $reaction) = split(/\|\|/, $line);

	$logger->debug("line '$line' produced ec '$ec' pfunc '$pfunc' sub1 '$sub1' sub2 '$sub2' sub3 '$sub3'") if $logger->is_debug();


	if ($ec =~ /^\s+(\S+)\s+/){
	    $ec = $1;
	}
	else{
	    $logger->fatal("Could not parse ec '$ec'.  Line was '$line'. Skipping to next record.");
	    next;
	}


	if ($pfunc =~ /^\s+(\S+)\s*$/){
	    $pfunc = $1;
	    $pfunc =~ s/\s+$//;
	}
	else{
	    $logger->fatal("Could not parse pfunc '$pfunc'.  Line was '$line'. Skipping to next record.");
	    next;
	}


	if ($sub1 =~ /^\s([\S\s]+)\s*$/){
	    $sub1 = $1;
	    $sub1 =~ s/\s+$//;
	}
	else{
	    $logger->fatal("Could not parse sub1 '$sub1'.  Line was '$line'.  Skipping to next record.");
	    next;
	}


	if ($sub2 =~ /^\s([\S\s]+)\s*$/){
	    $sub2 = $1;
	    $sub2 =~ s/\s+$//;
	}
	else{
	    $logger->fatal("Could not parse sub2 '$sub2'.  Line was '$line'.  Skipping to next record.");
	    next;
	}


	if ($sub3 =~ /^\s+([\S\s]+)\s*$/){
	    $sub3 = $1;
	    $sub3 =~ s/\s+$//;
	}
	else{
	    $logger->fatal("Could not parse sub3 '$sub3'.  Line was '$line'.  Skipping to next record.");
	    next;
	}


	if ($reaction =~ /^\s+([\S\s]+)\s*$/){
	    $reaction = $1;
	    $reaction =~ s/\s+$//;
	}
	else{
	    $logger->fatal("Could not parse reaction '$reaction'.  Line was '$line'.  Skipping to next record.");
	    next;
	}




	#
	# Need to split the ec number and store hierarchically
	#
	my ($e1,$e2,$e3,$e4) = split(/\./,$ec);
	$logger->debug("ec '$ec' was parsed: e1 '$e1' e2 '$e2' e3 '$e3' e4 '$e4'") if $logger->is_debug();

	#
	# Here we make sure that the subfunctions 1-3 are defined if the ec sub-category is defined (that is not == '-')
	#
	&check_ec($pfunc,$sub1,$sub2,$sub3,$e1,$e2,$e3,$e4);

	#
	# Sanity checks complete.  Resume our life-saving work here...
	#
	my $ec1 = $e1 . ".-.-.-";

	my $ec2 = $e1 . '.' . $e2 . '.-.-';
	my $ec3 = $e1 . '.' . $e2 . '.' . $e3 . '.-';
	my $ec4 = $e1 . '.' . $e2 . '.' . $e3 . '.' .$e4;
	

	$logger->debug("ec1 '$ec1' ec2 '$ec2' ec3 '$ec3' ec4 '$ec4'") if $logger->is_debug();


	#
	# In some instances the egad.prot_function.sub3 field contains the string 'Deleted entry'.
	# Whenever that occurs we need to insert the appropriate 'is_obsolete: true' into the output OBO file.
	# Also need to prepend OBSOLETE to the 'def:' obo.field.
	#
	if (lc($sub3) eq 'deleted entry'){
	    $deleted_entry = 1;
	    $ctr->{'deleted'}++;

	    #
	    # Preserve the 'Deleted entry' string in the OBO file's defintion field
	    # Note: reaction maps to defintion...
	    #
	    $reaction .= $sub3;
	    $sub3 = undef;
	}


	#
	# Also, in some instances, the egad.prot_function.sub3 field contains the string 'Transferred entry'.
	# Whenever that occurs, we need to insert the appropriate 'is_obsolete: true' into the output OBO file.
	# Also need to prepend OBSOLETE to the 'def:' obo.field.
	# Also need to add to the def field that the new ec# number should be referred to (in GO.obo style)
	#
	my $refer;

	if ( $sub3 =~ /Transferred entry/) {
	    $transferred_entry = 1;
	    $ctr->{'transferred'}++;

	    if ($sub3 =~ /Transferred entry ([\d\.]+)/){
		$refer = $1;
	    }
	    else{
		$logger->warn("Could not parse '$sub3', will attempt to split...");

		my $junk;
		($junk, $refer) = split(/:/, $sub3);
		$refer =~ s/\s+//;
	    }

	    #
	    # Preserve the 'Transferred entry' string in the OBO file's definition field
	    # Remember that reaction maps to definition...
	    $reaction .= $sub3;
	    $sub3 = undef;

	}



	#
	# This is where we build the fully qualified names
	#
	my ($pfunc_name, $sub1_name, $sub2_name, $sub3_name);

	if (defined($pfunc)){
	    $pfunc_name .= $pfunc;
	    $sub1_name  .= $pfunc;
	    $sub2_name  .= $pfunc;
	    $sub3_name  .= $pfunc;
	}
	if ( (defined($sub1)) and ( $sub1 !~ /NULL/) ) {


	    $sub1_name .= '||' . $sub1;
	    $sub2_name .= '||' . $sub1;
	    $sub3_name .= '||' . $sub1;

	}
	if ( (defined($sub2)) and ( $sub2 !~ /NULL/ ) ) {

	    $sub2_name .= '||' . $sub2;
	    $sub3_name .= '||' . $sub2;
	}

	if ( (defined($sub3)) and ($sub3 !~ /NULL/ ) ) {

	    $sub3_name .= '||' . $sub3;
	}

	$name_hash->{$sub3_name}++;




	#
	# Here we build the data structure.
	#

	#
	# Processing function and ec category 1
	#
	if ( (exists $hash->{$ec1}->{'name'}) and (defined($hash->{$ec1}->{'name'}))) {

	    my $old = $hash->{$ec1}->{'name'};
	    if ($old ne $pfunc_name){
		$logger->fatal("Fatal error has occured, though I will permit processing to continue in order to identify all other errors... old '$old' != pfunc_name '$pfunc_name' for ec '$ec' ec1 '$ec1' line '$line'");
		next;
	    }
	}
	else{

	    $hash->{$ec1}->{'name'} = $pfunc_name;
	}

	#
	# Processing sub-function1 and ec category 2
	#
	if ( (exists $hash->{$ec2}->{'name'}) and (defined($hash->{$ec2}->{'name'}))) {

	    my $old = $hash->{$ec2}->{'name'};
	    if ($old ne $sub1_name){
		$logger->fatal("Fatal error has occured, though I will permit processing to continue in order to identify all other errors... old '$old' != sub1_name '$sub1_name' for ec '$ec' ec1 '$ec1' ec2 '$ec2' line '$line'");
		next;
	    }
	}
	else {
#	my $name = $pfunc . "||" . $sub1;
#	my $name = $sub1;
#	$hash->{$ec2}->{'name'} = $name;
	    $hash->{$ec2}->{'name'} = $sub1_name;
#	$hash->{$ec2}->{'is_a'} = $prefix . ':' . $ec1;
	    $hash->{$ec2}->{'is_a'} = $ec1;
	}


	#
	# Processing sub-function2 and ec category 3
	#
	if ( (exists $hash->{$ec3}->{'name'}) and (defined($hash->{$ec3}->{'name'})) ) {

	    my $old = $hash->{$ec3}->{'name'};
	    if ($old ne $sub2_name){
		$logger->fatal("Fatal error has occured, though I will permit processing to continue in order to identify all other errors... old '$old' != sub2_name '$sub2_name' for ec '$ec' ec1 '$ec1' ec2 '$ec2' ec3 '$ec3' line '$line'");
		next;
	    }
	}
	else {
#	my $name = $pfunc . '||' . $sub1 . '||' . $sub2;
#	my $name = $sub2;
#	$hash->{$ec3}->{'name'} = $name;
	    $hash->{$ec3}->{'name'} = $sub2_name;
#	$hash->{$ec3}->{'is_a'} = $prefix . ':' . $ec2;
	    $hash->{$ec3}->{'is_a'} = $ec2;
	}
	


	#
	# Processing sub-function3 and ec category 4
	#
	if ( (exists $hash->{$ec4}->{'name'}) and  (defined($hash->{$ec4}->{'name'})) ) {

	    my $old = $hash->{$ec4}->{'name'};
	    if ($old ne $sub3_name){
		$logger->fatal("Fatal error has occured, though I will permit processing to continue in order to identify all other errors... old '$old' != sub3_name '$sub3_name' for ec '$ec' ec1 '$ec1' ec2 '$ec2' ec3 '$ec3' ec4 '$ec4' line '$line'");
		next;
	    }
	}
	else{
#	my $name = $pfunc . '||' . $sub1 . '||' . $sub2 . '||' . $sub3;
#	my $name = $sub3;
#	$hash->{$ec4}->{'name'} = $name;
	    $hash->{$ec4}->{'name'} = $sub3_name;
#	$hash->{$ec4}->{'is_a'} = $prefix . ':' . $ec3;
	    $hash->{$ec4}->{'is_a'} = $ec3;
	}


	#
	# The egad.prot_function.reaction is stored in the obo.def field.
	#
	if ( (exists $hash->{$ec4}->{'def'}) and (defined($hash->{$ec4}->{'def'})) ) {

	    my $old = $hash->{$ec4}->{'def'};
	    if ($old ne $sub3){
		$logger->fatal("Fatal error has occured, though I will permit processing to continue in order to identify all other errors... old '$old' != reaction '$reaction' for ec '$ec' line '$line'");
		next;
	    }

	}
	else {

	    my $reaction_string;
	    my $obflag=0;

	    if ($reaction !~ /NULL/){
		$logger->warn("Found a NULL reaction for ec '$ec'.  Setting the reaction to undef") ;
		$reaction = undef;
	    }
	    

	    if ($deleted_entry == 1 ){
		#
		# 'Deleted entry' was found in egad.prot_function.sub3 field.
		# 
		$reaction_string = 'OBSOLETE. '. $reaction;
		$obflag=1;
	    }
	    if ($transferred_entry == 1 ){
		#
		# 'Transferred entry' was found in egad.prot_function.sub3 field.
		#
		if ($obflag == 0){
		    $reaction_string = 'OBSOLETE. ' . $reaction;
		}
		$reaction_string .= ' Refer to ' . $refer . ' instead.';
	    }
	    

	    $hash->{$ec4}->{'def'} = $reaction_string;
	}




#    print Dumper $hash;die;

#    print STDERR "ec '$ec' pfunc '$pfunc' sub1 '$sub1' sub2 '$sub2' sub3 '$sub3' reaction '$reaction'\n";
#    die;


#    $hash->{$compartment}->{$mainrole}->{'name'} = $mainrole;
#    push ( @ { $hash->{$compartment}->{$mainrole}->{'array'}}, { subrole => $subrole,
#								 role_id => $role_id
#							     });

#    push ( @ { $hash->{$compartment}->{$mainrole}}, { subrole   => $subrole,
#						      role_id   => $role_id,
#						      roleorder => $roleorder
	#					  });
	
    }



    foreach my $namey (sort keys %{$name_hash} ){

	if ($name_hash->{$namey} > 1 ){
	    $logger->fatal("'$namey' occured '$name_hash->{$namey}'");
	}

    }



    $logger->fatal("name_hash:". Dumper $name_hash);

#die;
    $logger->fatal("Number of transferred entries encountered: '$ctr->{'transferred'}' and the number of deleted entries: '$ctr->{'deleted'}'");
#$logger->logdie("hash:" .  Dumper $hash);

    &write_terms($hash);
}


#--------------------------------------------------------------------------------
# write_terms()
#
# This function should write each node
# in the following manner:
#
#
# [Term]
# id: EC:0000001
# name: OXIDOREDUCTASE
# xref_analog: EC:1.-.-.-
#
# [Term]
# id: EC:0000002
# name: ACTING ON THE CH-OH GROUP OF DONOR
# xref_analog: EC:1.1.-.-
# is_a: EC:0000001
#
# [Term]
# id: EC:0000003
# name: WITH NAD(+) OR NADP(+) AS ACCEPTOR
# xref_analog: EC:1.1.1.-
# is_a: EC:0000002
#
# [Term]
# id: EC:0000004
# name: Alcohol dehydrogenase
# def: An alcohol + NAD(+) = an aldehyde or ketone + NADH
# xref_analog: EC:1.1.1.1
# is_a: EC:0000003
#
#-----------------------------------------------------------------------------------
sub write_terms {

    my $hash = shift;
    my $keycount = 0;
    my $compartment;

    #
    # Lookup hashes
    #
    my $ec2id = {};
    my $id2ec = {};



    $logger->debug("Building the lookup hashes...") if $logger->is_debug();

    foreach my $ec (keys %{$hash}){
	
	$keycount++;

	my $ident = &make_ident($keycount, 7);


	$ec2id->{$ec} = $ident;
	$id2ec->{$ident} = $ec;

    }

    
    #
    # Re-initialize
    #
    $keycount=0;


    $logger->debug("About to start outputting terms to OBO file...") if $logger->is_debug();

    foreach my $ec (keys %{$hash}){
	
	$keycount++;

	my $ident = &make_ident($keycount, 7);



	my $is_a        = $hash->{$ec}->{'is_a'}; 
	my $name        = $hash->{$ec}->{'name'}; 
	my $xref_analog = $ec;
	my $def         = $hash->{$ec}->{'def'};



	my $node = "[Term]\n".
	"id: $ident\n".
	"name: $name\n";
	

	$node .= "def: \"$def\"\n" if (defined($def));

	if (defined($xref_analog)){
	    $node .= "xref_analog: " . $prefix . ':' . $xref_analog ."\n";
	}

	if (defined($is_a)){
	    
	    if ((exists ($ec2id->{$is_a})) and (defined($ec2id->{$is_a}))) {

		$node .= "is_a: " . $ec2id->{$is_a} . "\n";

	    }
	    else{
		$logger->logdie("Could not find an EC OBO identifier for the ec# '$is_a'. ec2id:" . Dumper $ec2id);
	    }
	}


	if ( (defined($def)) and ($def =~ /OBSOLETE/)) {
	    $node .= "is_obsolete: true\n";
	}

	print OUTFILE $node . "\n";
	
    }

#    print "ec2id:" . Dumper $ec2id;
#    print "id2ec:" . Dumper $id2ec;




}


#
# Fancy-schmancy way of assigning properly formatted ontology identifier values.
#
sub make_ident {

    my ($id, $threshold) = @_;


    $threshold++;

    while ( length($id) < $threshold){
	
	my $zero = "0";
	$id = $zero . $id;

#	print "id '$id'\n";
    }

    $logger->debug("id '$id'") if $logger->is_debug();

    my $full_id = $prefix . ':' . $id;
    return $full_id;
}


#
# Simply output some content which is appropriate for the OBO file and will also permit DAGEditor to load.
#
sub write_header {


    $logger->debug("About to write the OBO file header") if $logger->is_debug();


    my $header = "format-version: 1.0\n".
    "date: 24:06:2004 16:32\n".
    "saved-by: sundaram\n".
    "auto-generated-by: DAG-Edit 1.417\n".
    "default-namespace: EC.ontology\n".
    "remark: autogenerated-by\\\:     $0 version 1.1\\nsaved-by\\\:             sundaram\\ndate\\\:                 Thu Feb 26 15\\\:56\\\:55 EST 2004\\nversion\\\:\n";
    

    print OUTFILE $header . "\n\n";



}


#
# In the past, some of the subfunctions were not defined while the ec sub-categories were defined.
# This was incorrect data.  Looks like the egad.prot_function table has since been updated.  However, I do not want to make
# assumptions about the quality of the incoming data...
#
sub check_ec {


    my ($pfunc, $sub1, $sub2, $sub3, $ec1, $ec2, $ec3, $ec4) = @_;


    if ( (defined($ec4)) and (defined($ec3)) and (defined($ec2)) and (defined($ec1)) ) {
	#
	# All four ec categories were defined, therefore we expect all four function and subfunctions to be defined.
	#
	if ( (defined($sub3)) and (defined($sub2)) and (defined($sub1)) and (defined($pfunc)) ){
	    #
	    # Looks good.
	    #
	    $logger->debug("$ec1.$ec2.$ec3.$ec4") if $logger->is_debug();
	}
	else{
	    $logger->fatal("sub-function 3 was not defined for $ec1.$ec2.$ec3.$ec4") if (!defined($sub3));
	    $logger->fatal("sub-function 2 was not defined for $ec1.$ec2.$ec3.$ec4") if (!defined($sub2));
	    $logger->fatal("sub-function 1 was not defined for $ec1.$ec2.$ec3.$ec4") if (!defined($sub1));
	    $logger->fatal("function was not defined for $ec1.$ec2.$ec3.$ec4") if (!defined($pfunc));
	}
    }
    elsif (  (defined($ec3)) and (defined($ec2)) and (defined($ec1)) ) {
	#
	# Only ec categories 1,2,3 were defined, therefore we expect only function, sub1function, sub2function to be defined.
	#
	if ( (defined($sub2)) and (defined($sub1)) and (defined($pfunc)) ){
	    #
	    # Looks good.
	    #
	    $logger->debug("$ec1.$ec2.$ec3") if $logger->is_debug();
	}
	else{
	    $logger->fatal("sub-function 2 was not defined for $ec1.$ec2.$ec3.$ec4") if (!defined($sub2));
	    $logger->fatal("sub-function 1 was not defined for $ec1.$ec2.$ec3.$ec4") if (!defined($sub1));
	    $logger->fatal("function was not defined for $ec1.$ec2.$ec3.$ec4") if (!defined($pfunc));
	}
	if (defined($sub3)){
	    #
	    # We don't expect sub-function3 to be defined.
	    #
	    $logger->warn("Remarkably, sub-function 3 '$sub3' was defined and yet we're processing '$ec1.$ec2.$ec3'");
	}
    }
    elsif (  (defined($ec2)) and (defined($ec1)) ) {
	#
	# Only ec categories 1,2 were defined, therefore we expect only function, sub1function to be defined.
	#
	if ( (defined($sub1)) and (defined($pfunc)) ){
	    #
	    # Looks good.
	    #
	    $logger->debug("$ec1.$ec2") if $logger->is_debug();
	}
	else{
	    $logger->fatal("sub-function 1 was not defined for $ec1.$ec2.$ec3.$ec4") if (!defined($sub1));
	    $logger->fatal("function was not defined for $ec1.$ec2.$ec3.$ec4") if (!defined($pfunc));
	}
	if (defined($sub3)){
	    #
	    # We don't expect sub-function3 to be defined.
	    #
	    $logger->warn("Remarkably, sub-function 3 '$sub3' was defined and yet we're processing '$ec1.$ec2.$ec3'");
	}
	if (defined($sub2)){
	    #
	    # We don't expect sub-function2 to be defined.
	    #
	    $logger->warn("Remarkably, sub-function 2 '$sub2' was defined and yet we're processing '$ec1.$ec2'");
	}
    }
    elsif ( (defined($ec1)) ) {
	#
	# Only ec categories 1 was defined, therefore we expect only function to be defined.
	#
	if ( (defined($pfunc)) ){
	    #
	    # Looks good.
	    #
	    $logger->debug("$ec1") if $logger->is_debug();
	}
	else{
	    $logger->fatal("function was not defined for $ec1.$ec2.$ec3.$ec4") if (!defined($pfunc));
	}
	if (defined($sub3)){
	    #
	    # We don't expect sub-function3 to be defined.
	    #
	    $logger->warn("Remarkably, sub-function 3 '$sub3' was defined and yet we're processing '$ec1.$ec2.$ec3'");
	}
	if (defined($sub2)){
	    #
	    # We don't expect sub-function2 to be defined.
	    #
	    $logger->warn("Remarkably, sub-function 2 '$sub2' was defined and yet we're processing '$ec1.$ec2'");
	}
	if (defined($sub1)){
	    #
	    # We don't expect sub-function1 to be defined.
	    #
	    $logger->warn("Remarkably, sub-function 1 '$sub1' was defined and yet we're processing '$ec1'");
	}
    }
    else{
	$logger->error("Un-expected condition. ec1 '$ec1' ec2 '$ec2' ec3 '$ec3' ec4 '$ec4'");
    }
}

