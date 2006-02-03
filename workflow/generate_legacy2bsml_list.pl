#!/usr/local/bin/perl
#--------------------------------------------------------------------------------------------
#
# cvs: ANNOTATION/ergatis/workflow/generate_asmbl_list.pl
# date: 2005-01-25
# contact: sundaram@tigr.org
#
#
# Similar to generate_input_list.pl except instead of iterating on list of BSML files,
# will iterate on list of input asmbl_ids.
#
# Context: Should run in following order:
# 1) generate_asmbl_list.pl
# 2) generate_groups.pl
# 3) generate_subflow.pl
#
#
# $Id$
#
#
#--------------------------------------------------------------------------------------------
=head1  NAME 

generate_asmbl_list.pl - Default output is a workflow iterator that
can be used to iterator over a set of asmbl_ids

=head1 SYNOPSIS

USAGE:  generate_asmbl_list

=head1 OPTIONS

=item *

B<--debug,-d> Debug level.  Use a large number to turn on verbose debugging. 

=item *

B<--log,-l> Log file

=item *

B<--help,-h> This help message

=head1   DESCRIPTION

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

BEGIN {
    require '/usr/local/devel/ANNOTATION/cas/lib/site_perl/5.8.5/Workflow/Logger.pm';
    import Workflow::Logger;
}
use File::Basename;

umask(0000);

my %options = ();


my $results = GetOptions (\%options, 
                          'control_file|f=s', 
			  'output|o=s',
			  'log|l=s',
                          'debug=s', 
                          'help|h' ) || pod2usage();

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Workflow::Logger::get_logger();

my $database_key      = '$;DATABASE$;';
my $asmbl_id_key      = '$;ASMBL_ID$;';
my $sequence_type_key = '$;SEQUENCE_TYPE$;';
my $subflow_key       = '$;SUBFLOW_NAME$;';
my $euk_key           = '$;EUK$;';
my $ntprok_key        = '$;NTPROK$;';
my $include_genefinders_key = '$;INCLUDE_GENEFINDERS$;';
my $exclude_genefinders_key = '$;EXCLUDE_GENEFINDERS$;';


my $valid_organism_types = { 'prok'   => 1,
			     'ntprok' => 1,
			     'euk'    => 1
			 };



if (&verify_control_file($options{'control_file'})){

    if (defined($options{'output'})){

	my $iteratorconf = {
	    $database_key      => [],
	    $asmbl_id_key      => [],
	    $sequence_type_key => [],
	    $subflow_key       => [],
	    $euk_key           => [],
	    $ntprok_key        => [],
	    $include_genefinders_key => [],
	    $exclude_genefinders_key => []
	};



	#
	# Read in the information from asmbl_file OR asmbl_list
	#
	&get_list_from_file($iteratorconf,$options{'control_file'});

	#
	# Output the lists
	#
	&output_lists($iteratorconf, $options{'output'});
    }
    else{
	$logger->logdie("output was not defined");
    }
}


exit;
						     
#---------------------------------------------------------------------------------------------------------
#
#                           END OF MAIN  --  SUBROUTINES FOLLOW
#
#---------------------------------------------------------------------------------------------------------



#-----------------------------------------
# get_list_from_file()
#
#-----------------------------------------
sub get_list_from_file{

    my ($iconf, $f) = @_;

    my $contents = &get_file_contents($f);

    my $orghash = &get_organism_hash($contents);


    foreach my $database (sort keys %{$orghash} ){ 

	my $organism_type = $orghash->{$database}->{'organism_type'};
	my $include_genefinders = $orghash->{$database}->{'include_genefinders'};
	my $exclude_genefinders = $orghash->{$database}->{'exclude_genefinders'};

	foreach my $infohash ( @{$orghash->{$database}->{'infohash'}} ) {

	    my $sequence_type = $infohash->{'sequence_type'};
	    my $asmbl_id = $infohash->{'asmbl_id'};

	    my $subflow_name = $database . "_" . $asmbl_id;

	    &add_entry_to_conf($iconf, $database, $asmbl_id, $subflow_name, $sequence_type, $organism_type, $include_genefinders, $exclude_genefinders);
	}
    }
    
}


#-------------------------------------------
# add_entry_to_conf()
#
#-------------------------------------------
sub add_entry_to_conf{

    my($iteratorconf,$database, $asmbl_id, $subflow_name, $sequence_type, $organism_type, $include_genefinders, $exclude_genefinders) = @_;

    push( @{$iteratorconf->{$database_key}}, $database );
    push( @{$iteratorconf->{$asmbl_id_key}}, $asmbl_id );
    push( @{$iteratorconf->{$subflow_key}}, $subflow_name );
    push( @{$iteratorconf->{$sequence_type_key}}, $sequence_type );
    push( @{$iteratorconf->{$include_genefinders_key}}, $include_genefinders );
    push( @{$iteratorconf->{$exclude_genefinders_key}}, $exclude_genefinders );

    if ($organism_type eq 'ntprok'){
	# legacy2bsml.pl --ntprok=1 --euk=0
	push( @{$iteratorconf->{$ntprok_key}}, 1 );
	push( @{$iteratorconf->{$euk_key}}, 0 );
    }
    elsif ($organism_type eq 'euk'){
	# legacy2bsml.pl --ntprok=0 --euk=1
	push( @{$iteratorconf->{$ntprok_key}}, 0 );
	push( @{$iteratorconf->{$euk_key}}, 1 );
    }
    elsif ( $organism_type eq 'prok'){
	# legacy2bsml.pl --ntprok=0 --euk=0
 	push( @{$iteratorconf->{$ntprok_key}}, 0 );
	push( @{$iteratorconf->{$euk_key}}, 0 );
    }
}

#-------------------------------------------
# verify_control_file()
#
#-------------------------------------------
sub verify_control_file {
    my $file = shift;

    if ((defined($file)) && (-e $file) && (-r $file) && (-s $file)) {
	# control file was defined, exists, has read permissions
	# and has content
	return 1;
    }
    else {
	if (!defined($file)){
	    $logger->logdie("control file '$file' was not defined");
	}
	if (!-e $file){
	    $logger->logdie("control file '$file' does not exist");
	}
	if (!-r $file){
	    $logger->logdie("control file '$file' does not have read permissions");
	}
	if (!-s $file){
	    $logger->logdie("control file '$file' has zero content");
	}
    }
}



#---------------------------------------------
# output_lists()
#
#---------------------------------------------
sub output_lists {

    my ($iteratorconf, $output) = @_;

    open FILE, "+>$output" or $logger->logdie("Can't open output file $output");

    foreach my $key (keys %$iteratorconf){

	print FILE "$key=",join(',',@{$iteratorconf->{$key}}),"\n";

    }

    close FILE;

}



#---------------------------------------------
# get_file_contents()
#
#---------------------------------------------
sub get_file_contents {

    my $f = shift;

    open (CONTROLFILE,  "<$f") or die "Could not open file '$f': $!";
    
    my @lines = <CONTROLFILE>;

    chomp @lines;

    return \@lines;
}

#---------------------------------------------
# get_organism_hash()
#
#---------------------------------------------
sub get_organism_hash {

    my $contents = shift;

    my $hash = {};

    my $database;

    my $unique_asmbl_id_values = {};

    my $linectr=0;

    foreach my $line (@{$contents}){

	$linectr++;

	if ($line =~ /^\s*$/){
	    next; # skip blank lines
	}
	elsif ($line =~ /^\#/){
	    next; # skip comment lines
	}
	elsif ($line =~ /^\-\-/){
	    next; # skip -- lines
	}
	else{

	    if ($line =~ /^database:(\S+)\s+organism_type:(\S+)\s+include_genefinders:(\S*)\s+exclude_genefinders:(\S*)/){
		
		$database      = $1;
		my $organism_type = $2;
		my $include_genefinders = $3;
		my $exclude_genefinders = $4;


#		print "organism_type '$organism_type' include_genefinder '$include_genefinders' exclude_genefinder '$exclude_genefinders'\n";


		if (&verify_organism_type($organism_type, $linectr)){

		    ($include_genefinders, $exclude_genefinders) = &verify_and_set_genefinders($include_genefinders, $exclude_genefinders, $linectr);



		    if (( exists $hash->{$database}) &&
			(defined($hash->{$database}))){
			
			$logger->logdie("This database '$database' was already encountered in the control file!");
		    }
		    else {
		    #
			# Encountered this organism/database for the first time while processing the control file contents
			# therefore go ahead and declare the organism's hash attributes
			#
			$hash->{$database} = { 'asmbl_id_list'       => [],
					       'infohash'            => [],
					       'organism_type'       => $organism_type,
					       'include_genefinders' => $include_genefinders,
					       'exclude_genefinders' => $exclude_genefinders };
			
		    }
		}
	    }
	    elsif ($line =~ /^\s*(\d+)\s*sequence_type:(\S*)/){

		&store_asmbl_id_and_sequence_type($database, $1, $2, $linectr, $unique_asmbl_id_values, $hash);
	    }
	    elsif ($line =~ /^\s*(\d+)\s*/){

		&store_asmbl_id_and_sequence_type($database, $1, undef, $linectr, $unique_asmbl_id_values, $hash);
	    }
	    else {
		$logger->logdie("Could not parse line number '$linectr' - line was '$line'");
	    }
	}
    }

    return $hash;

}



#------------------------------------------------------
# store_asmbl_id_and_sequence_type()
#
#------------------------------------------------------
sub store_asmbl_id_and_sequence_type {
    
    my ($database, $asmbl_id, $sequence_type, $linectr, $unique_asmbl_id_values, $hash) = @_;

    print "database '$database' asmbl_id '$asmbl_id' sequence_type '$sequence_type'\n";


    $sequence_type = &verify_and_set_sequence_type($sequence_type, $linectr);
    
    if (( exists $unique_asmbl_id_values->{$database}->{$asmbl_id}) && 
	(defined($unique_asmbl_id_values->{$database}->{$asmbl_id}))){
	#
	# At the same time, need to ensure that the asmbl_id values are unique
	#
	
	$logger->warn("asmbl_id '$asmbl_id' for organism '$database' was already stored in the lookup");
    }
    else {
	#
	# Update the unique_asmbl_id_values hash to ensure unique assembly identifiers are processed
	#
	$unique_asmbl_id_values->{$database}->{$asmbl_id}++;
	
	#
	# Store the next valid asmbl_id in the organism hash
	#
	my $infohash = { 'asmbl_id' => $asmbl_id,
			 'sequence_type' => $sequence_type };
	
	push( @{$hash->{$database}->{'infohash'}}, $infohash);
    }
}





#-----------------------------------------------------
# verify_organism_type()
#
#-----------------------------------------------------
sub verify_organism_type {

    my ($organism_type, $line) = @_;
    
    if ((exists $valid_organism_types->{$organism_type}) && 
	(defined($valid_organism_types->{$organism_type}))){
	
	return 1;
    }
    else {
	$logger->logdie("Encountered an invalid organism type '$organism_type' at line '$line'");
    }

}

#-----------------------------------------------------
# verify_and_set_sequence_type()
#
#-----------------------------------------------------
sub verify_and_set_sequence_type {

    my ($sequence_type, $linectr) = @_;
    
    
    if ((!defined($sequence_type)) or ($sequence_type =~ /^\s*$/)) {
	$sequence_type = "none";
    }
    
    #
    # Lame. That's all the checking for now.
    #
    return ($sequence_type);
    
    
}

#-----------------------------------------------------
# verify_and_set_genefinders()
#
#-----------------------------------------------------
sub verify_and_set_genefinders {

    my ($include_genefinder, $exclude_genefinder, $linectr) = @_;
    
    if ((!defined($exclude_genefinder)) || 
	($exclude_genefinder =~ /^\s*$/)){

	$exclude_genefinder = 'none';
    } 
    if ((!defined($include_genefinder)) ||
	($include_genefinder =~  /^\s*$/)){

	$include_genefinder = 'all';  
    }

    if ($include_genefinder eq 'none'){
	$exclude_genefinder = 'all';
    }
    elsif ($exclude_genefinder eq 'all'){
	$include_genefinder = 'none';
    }


    return ($include_genefinder, $exclude_genefinder);
    
    
}
