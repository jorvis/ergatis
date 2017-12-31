#!/usr/local/bin/perl
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Data::Dumper;


## Command-line options
my ($debug, $infile, $infile2, $help, $logfile, $man, $database, $asmbl_id);

my $results = GetOptions (
			  'debug=s'    => \$debug,
			  'help|h'     => \$help,
			  'logfile=s'  => \$logfile,
			  'man|m=s'    => \$man,
			  'infile=s'   => \$infile,
			  'infile2=s'  => \$infile2,
			  'database=s' => \$database,
			  'asmbl_id=s' => \$asmbl_id
			  );

&print_usage() if ($help);
&pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if($man);

my $fatalCtr=0;

if (!$infile){
    print STDERR "infile was not defined\n";
    $fatalCtr++;
}
else {
    if (&checkFile($infile) == 0){
	$fatalCtr++;
    }
}
if ($fatalCtr>0){
    &print_usage();
}

if (!defined($logfile)){
    $logfile = '/tmp/profileBSMLFile2.pl.log';
}

open(LOGFILE, ">$logfile") || die "Could not open log file '$logfile':$!";

my $invocation = "The script invocation was: $0 --infile=$infile --logfile=$logfile";
if ($infile2){
    $invocation .= " --infile2=$infile2";
}
if ($database){
    $invocation .= " --database=$database";
}
if ($asmbl_id){
    $invocation .= " --asmbl_id=$asmbl_id";
}
print LOGFILE $invocation . "\n\n";

my $masterProfileLookup={};
&profileBSMLFile($infile);
&recordFindings($infile);

if (defined($infile2)){
    if (&checkFile($infile2)){
	&profileBSMLFile($infile2);
	&recordFindings($infile2);
	&compareFiles($infile, $infile2);
    }
}

print "$0 program execution has completed\n";
print "Please review log file '$logfile'\n";
exit(0);
  
##--------------------------------------------------------------------------------------------------------------------------------------------
##
##  END OF MAIN -- SUBROUTINES FOLLOW
##
##--------------------------------------------------------------------------------------------------------------------------------------------


##----------------------------------------------------------------------
## addTimeStamp()
##
##----------------------------------------------------------------------
sub addTimeStamp {
    my ($msg) = @_;
    my $date = `date`;
    chomp $date;
    print LOGFILE "$msg $date\n"; 
}

##----------------------------------------------------------------------
## profileBSMLFile()
##
##----------------------------------------------------------------------
sub profileBSMLFile {

    my ($file) = @_;

    &addTimeStamp("Profiling of BSML file '$file' initiated at ");

    open (INFILE, "<$file") || die "Could not open infile '$file' in read mode";

    my $lineCtr=0;

    while (my $line = <INFILE>){

	chomp $line;
	$lineCtr++;
	
	if ($line =~ /^\s+<\/(\S+)/){
	    $masterProfileLookup->{$file}->{'closingTagCtr'}->{$1}++;
	    next;
	}
	if ($line =~/<Sequence /){

	    if ($line =~ /class=\"(\S+)\"/){
		my $class = $1;
		if ($line =~ /id=\"\S+_seq\"/){
		    $masterProfileLookup->{$file}->{'nonLoadedSequenceClassCtr'}->{$class}++;
		    $masterProfileLookup->{$file}->{'nonLoadedSequenceCtr'}++;
		}
		else {
		    $masterProfileLookup->{$file}->{'sequenceClassCtr'}->{$1}++;
		    $masterProfileLookup->{$file}->{'sequenceCtr'}++;
		}
	    }
	    else {
		die "Found sequence with no class at line '$lineCtr' line was '$line'";
	    }
	}
	elsif ($line =~/<Feature /){
	    $masterProfileLookup->{$file}->{'featureCtr'}++;
	    if ($line =~ /class=\"(\S+)\"/){
		$masterProfileLookup->{$file}->{'featureClassCtr'}->{$1}++;
	    }
	    else {
		die "Found feature with no class at line '$lineCtr' line was '$line'";
	    }
	}
	elsif ($line =~ /<Attribute name=\"(.+)\" content/){
	    $masterProfileLookup->{$file}->{'attributeTypeCtr'}->{$1}++;
	    $masterProfileLookup->{$file}->{'attributeCtr'}++;
	}
	elsif ($line =~ /<Cross-reference database=\"(.+)\" identifier=/){
	    $masterProfileLookup->{$file}->{'crossReferenceTypeCtr'}->{$1}++;
	    $masterProfileLookup->{$file}->{'crossReferenceCtr'}++;
	}
	elsif ($line =~ /<Organism genus=\"(.+)\" species=\"(.+)\"/){
	    $masterProfileLookup->{$file}->{'organismTypeCtr'}->{$1,$2}++;
	    $masterProfileLookup->{$file}->{'organismCtr'}++
	}
	elsif ($line =~ /<Feature-group-member /){
	    $masterProfileLookup->{$file}->{'featureGroupMemberCtr'}++;
	    if ($line =~ /feature-type=\"(.+)\" featref/){
		$masterProfileLookup->{$file}->{'featureGroupMemberByClassCtr'}->{$1}++;
	    }
	    if ($line =~ /featref=\"(.+)\"/){
		$masterProfileLookup->{$file}->{'uniqueFeatureGroupMemberCtr'}->{$1}++;
	    }
	}
	elsif ($line =~ /<Feature-group /){
	    $masterProfileLookup->{$file}->{'featureGroupCtr'}++;
	}
	elsif ($line =~ /<Interval-loc /){
	    $masterProfileLookup->{$file}->{'intervalLocCtr'}++;
	}
	elsif ($line =~ /<Attribute-list>/){
	    $masterProfileLookup->{$file}->{'attributeListCtr'}++;
	}
	elsif ($line =~ /<Analysis /){
	    $masterProfileLookup->{$file}->{'analysisCtr'}++;
	}
	elsif ($line =~ /<Seq-pair-alignment /){
	    $masterProfileLookup->{$file}->{'seqPairAlignmentCtr'}++;
	}
	elsif ($line =~ /<Seq-pair-run /){
	    $masterProfileLookup->{$file}->{'seqPairRunCtr'}++;
	}
	elsif ($line =~ /<Link /){
	    $masterProfileLookup->{$file}->{'linkCtr'}++;
	    my $rel;
	    if ($line =~ /rel=\"(.+)\" href/){
		$rel = $1;
		$masterProfileLookup->{$file}->{'linkRelCtr'}->{$rel}++;
	    }
	    if ($line =~ /role=\"(.+)\"/){
		$masterProfileLookup->{$file}->{'linkRoleCtr'}->{$1}++;
	    }
	    if ($rel eq 'analysis'){
		if ($line =~ /href=\"(.+)\" role/){
		    $masterProfileLookup->{$file}->{'linkHrefCtr'}->{$1}++;
		}
	    }
	}
	elsif ($line =~ /<Multiple-alignment-table /){
	    $masterProfileLookup->{$file}->{'multipleAlignmentTableCtr'}++;
	}
	
    }

    $masterProfileLookup->{$file}->{'uniqueFeatureGroupMemberCtr'} = keys %{$masterProfileLookup->{$file}->{'uniqueFeatureGroupMemberCtr'}};

    $masterProfileLookup->{$file}->{'lineCtr'} = $lineCtr;

    &addTimeStamp("Profiling of BSML file '$file' completed at ");
    
}

##----------------------------------------------------------------------
## recordFindings()
##
##----------------------------------------------------------------------
sub recordFindings {

    my ($file) = @_;

    my $profileHeadings = { 'organismCtr' => "<Organism>",
			    'crossReferenceCtr' => "<Cross-reference>",
			    'sequenceCtr' => "<Sequence> (loadable)",
			    'nonLoadedSequenceCtr' => "<Sequence> (non-loadable)",
			    'featureCtr' => "<Feature>",
			    'attributeCtr' => "<Attribute>",
			    'attributeListCtr' => "<Attribute-list>",
			    'featureGroupCtr' => "<Feature-group>",
			    'featureGroupMemberCtr' => "<Feature-group-member>",
			    'uniqueFeatureGroupMemberCtr' => "<Feature-group-member> (unique)",
			    'intervalLocCtr' => "<Interval-loc>",
			    'analysisCtr' => "<Analysis>",
			    'seqPairAlignmentCtr' => "<Seq-pair-alignment>",
			    'seqPairRunCtr' => "<Seq-pair-run>",
			    'linkCtr' => "<Link>",
			    'multipleAlignmentTableCtr' => "<Multiple-alignment-table>"			  
			};

    my $profileLookupHeadings = { 'organismTypeCtr' => "<Organism> elements by genus-species were:",
				  'crossReferenceTypeCtr' => "<Cross-reference> elements by type were:",
				  'sequenceClassCtr' => "loadable <Sequence> elements by type were:",
				  'nonLoadedSequenceClassCtr' => "non-loadable <Sequence> elements by type were:",
				  'featureClassCtr' => "<Feature> elements by type were:",
				  'featureGroupMemberByClassCtr' => "<Feature-group-member> elements by class were:",
				  'attributeTypeCtr' => "<Attribute> elements by type were:",
				  'linkRelCtr' => "<Link> elements by rel were:",
				  'linkHrefCtr' => "<Link> elements by href were:",
				  'linkRoleCtr' => "<Link> elements by role are were:"
			      };
    
    
    print LOGFILE "\n\n$masterProfileLookup->{$file}->{'lineCtr'} lines were processed in BSML file '$file'\n\n";
    print LOGFILE "BSML element counts were:\n";

    foreach my $counterType (keys %{$profileHeadings} ){
	my $count=0;
	if (exists $masterProfileLookup->{$file}->{$counterType}){
	    $count = $masterProfileLookup->{$file}->{$counterType};
	}
	print LOGFILE "$profileHeadings->{$counterType} $count\n";
    }


    foreach my $lookupType (keys %{$profileLookupHeadings} ){
	if (exists $masterProfileLookup->{$file}->{$lookupType}){

	    &report($profileLookupHeadings->{$lookupType}, $masterProfileLookup->{$file}->{$lookupType});
	}
	else {
	    print LOGFILE "One of the files did not contain elements necessary to compute counts for this type '$lookupType'\n";
	}
    }


    print LOGFILE "\n\n";
}


##----------------------------------------------------------------------
## compareFiles()
##
##----------------------------------------------------------------------
sub compareFiles {

    my ($file1, $file2) = @_;

    print LOGFILE "Comparing '$file1' versus '$file2'\n\n";


    my $profileHeadings = { 'organismCtr' => "<Organism>",
			    'crossReferenceCtr' => "<Cross-reference>",
			    'sequenceCtr' => "<Sequence> (loadable)",
			    'nonLoadedSequenceCtr' => "<Sequence> (non-loadable)",
			    'featureCtr' => "<Feature>",
			    'attributeCtr' => "<Attribute>",
			    'attributeListCtr' => "<Attribute-list>",
			    'featureGroupCtr' => "<Feature-group>",
			    'featureGroupMemberCtr' => "<Feature-group-member>",
			    'uniqueFeatureGroupMemberCtr' => "<Feature-group-member> (unique)",
			    'intervalLocCtr' => "<Interval-loc>",
			    'analysisCtr' => "<Analysis>",
			    'seqPairAlignmentCtr' => "<Seq-pair-alignment>",
			    'seqPairRunCtr' => "<Seq-pair-run>",
			    'linkCtr' => "<Link>",
			    'multipleAlignmentTableCtr' => "<Multiple-alignment-table>"			  
			};

    my $profileLookupHeadings = { 'organismTypeCtr' => "<Organism> by genus-species",
				  'crossReferenceTypeCtr' => "<Cross-reference> by database",
				  'sequenceClassCtr' => "<Sequence> (loadable) by class",
				  'nonLoadedSequenceClassCtr' => "<Sequence> (non-loadable) by class",
				  'featureClassCtr' => "<Feature> by class",
				  'featureGroupMemberByClassCtr' => "<Feature-group-member> by class",
				  'attributeTypeCtr' => "<Attribute> by name",
				  'linkRelCtr' => "<Link> by rel",
				  'linkHrefCtr' => "<Link> by href",
				  'linkRoleCtr' => "<Link> by role"
			      };
    
    foreach my $counterType (keys %{$profileHeadings} ){
	my $count1=$masterProfileLookup->{$file1}->{$counterType};
	my $count2=$masterProfileLookup->{$file2}->{$counterType};
	if ($count1 != $count2){
	    print LOGFILE "$profileHeadings->{$counterType}\n".
	    "$file1: $count1\n".
	    "$file2: $count2\n\n";
	}
    }

    foreach my $lookupType (keys %{$profileLookupHeadings} ){

	## Collect all of the different types of keys present in either file
	my $uniqueKeysLookup = {};
 	foreach my $key (keys %{$masterProfileLookup->{$file1}->{$lookupType}} ){
	    $uniqueKeysLookup->{$key}++;
	}
 	foreach my $key (keys %{$masterProfileLookup->{$file2}->{$lookupType}} ){
	    $uniqueKeysLookup->{$key}++;
	}
	
	foreach my $key (keys %{$uniqueKeysLookup}) {
	    my $count1=$masterProfileLookup->{$file1}->{$lookupType}->{$key};
	    my $count2=$masterProfileLookup->{$file2}->{$lookupType}->{$key};
	    if ($count1 != $count2){
		print LOGFILE "$profileLookupHeadings->{$lookupType} '$key'\n".
		"$file1: $count1\n".
		"$file2: $count2\n\n";
	    }
	}
    }

}

##----------------------------------------------------------------------
## checkFile()
##
##----------------------------------------------------------------------
sub checkFile {

    my ($file) = @_;

    my $status=1;

    if (!-e $file){
	print "input BSML file '$file' does not exist\n";
	$status=0;
    }
    if (!-f $file){
	print "input BSML file '$file' is not a file\n";
	$status=0;
    }
    if (!-r $file){
	print "input BSML file '$file' does not have read permissions\n";
	$status=0;
    }
    if (!-s $file){
	print "input BSML file '$file' has no content\n";
	$status=0;
    }
    
    return $status;
}

##----------------------------------------------------------------------
## print_usage()
##
##----------------------------------------------------------------------
sub print_usage {
    print "Sample usage: perl $0 --infile [--infile2] [--logfile]\n";
    exit(1);
}

##----------------------------------------------------------------------
## report()
##
##----------------------------------------------------------------------
sub report {
    my ($msg, $lookup) = @_;
    
    print LOGFILE "\n\n$msg\n";
    foreach my $key ( sort keys %{$lookup}){
	print LOGFILE "$key\t\t$lookup->{$key}\n";
    }
}
