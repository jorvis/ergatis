#!/usr/local/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";
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
    import Ergatis::Logger;
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

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Ergatis::Logger::get_logger();

my $database_key      = '$;SOURCE_DATABASE$;';
my $asmbl_id_key      = '$;ASMBL_ID$;';
my $sequence_type_key = '$;SEQUENCE_TYPE$;';
my $subflow_key       = '$;SUBFLOW_NAME$;';
my $organism_key        = '$;SCHEMA_TYPE$;';
my $include_genefinders_key = '$;INCLUDE_GENEFINDERS$;';
my $exclude_genefinders_key = '$;EXCLUDE_GENEFINDERS$;';
my $alt_database_key = '$;ALT_DATABASE$;';
my $alt_species_key = '$;ALT_SPECIES$;';
my $tu_list_key = '$;TU_LIST_FILE$;';
my $model_list_key = '$;MODEL_LIST_FILE$;';


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
	    $organism_key      => [],
	    $include_genefinders_key => [],
	    $exclude_genefinders_key => [],
	    $alt_database_key => [],
	    $alt_species_key => [],
	    $tu_list_key => [],
	    $model_list_key => []

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

    my $orghash = &get_organism_hash($contents, $f);

    foreach my $database_type (sort keys %{$orghash} ){ 

	my $database      = $orghash->{$database_type}->{'database'};
	my $organism_type = $orghash->{$database_type}->{'organism_type'};
	my $include_genefinders = $orghash->{$database_type}->{'include_genefinders'};
	my $exclude_genefinders = $orghash->{$database_type}->{'exclude_genefinders'};
	my $alt_database = $orghash->{$database_type}->{'alt_database'};
	my $alt_species  = $orghash->{$database_type}->{'alt_species'};

	foreach my $infohash ( @{$orghash->{$database_type}->{'infohash'}} ) {

	    my $sequence_type = $infohash->{'sequence_type'};
	    my $asmbl_id = $infohash->{'asmbl_id'};
	    my $tu_list = $infohash->{'tu_list'};
	    my $model_list = $infohash->{'model_list'};
	    
	    my $subflow_name = $database_type . "_" . $asmbl_id;

	    &add_entry_to_conf($iconf, $database, $asmbl_id, $subflow_name, $sequence_type, $organism_type, $include_genefinders, $exclude_genefinders, $alt_database, $alt_species, $tu_list, $model_list);
	}
    }
    
}


#-------------------------------------------
# add_entry_to_conf()
#
#-------------------------------------------
sub add_entry_to_conf{

    my($iteratorconf,$database, $asmbl_id, $subflow_name, $sequence_type, $organism_type, $include_genefinders, $exclude_genefinders, $alt_database, $alt_species, $tu_list, $model_list) = @_;
    
    push( @{$iteratorconf->{$database_key}}, $database );
    push( @{$iteratorconf->{$asmbl_id_key}}, $asmbl_id );
    push( @{$iteratorconf->{$subflow_key}}, $subflow_name );
    push( @{$iteratorconf->{$sequence_type_key}}, $sequence_type );
    push( @{$iteratorconf->{$include_genefinders_key}}, $include_genefinders );
    push( @{$iteratorconf->{$exclude_genefinders_key}}, $exclude_genefinders );
    push( @{$iteratorconf->{$organism_key}}, $organism_type );
    push( @{$iteratorconf->{$alt_database_key}}, $alt_database );
    push( @{$iteratorconf->{$alt_species_key}}, $alt_species );
    push( @{$iteratorconf->{$tu_list_key}}, $tu_list );
    push( @{$iteratorconf->{$model_list_key}}, $model_list );

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

    my ($contents, $file) = @_;

    my $hash = {};

    my $database_type;

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

	    if ($line =~ /^database:(\S+)\s+organism_type:(\S+)\s+include_genefinders:(\S*)\s+exclude_genefinders:(\S*)\s+alt_database:(\S*)\s+alt_species:(\S*)\s*$/){
		
		my $database      = $1;
		my $organism_type = $2;
		my $include_genefinders = $3;
		my $exclude_genefinders = $4;

		my ($alt_database, $alt_species) = &verify_alt_database_and_species($5, $6);

		$database_type = $database ."_" .$organism_type;

		if (&verify_organism_type($organism_type, $linectr)){
		    
		    ($include_genefinders, $exclude_genefinders) = &verify_and_set_genefinders($include_genefinders, $exclude_genefinders, $linectr);

		    

		    if (( exists $hash->{$database_type}) &&
			(defined($hash->{$database_type}))){
			
			$logger->warn("This database_type '$database_type' was already encountered in the control file!");
		    }
		    else {
			#
			# Encountered this organism/database for the first time while processing the control file contents
			# therefore go ahead and declare the organism's hash attributes
			#
			$hash->{$database_type} = { 'database'            => $database,
						    'asmbl_id_list'       => [],
						    'infohash'            => [],
						    'organism_type'       => $organism_type,
						    'include_genefinders' => $include_genefinders,
						    'exclude_genefinders' => $exclude_genefinders,
						    'alt_database'        => $alt_database,
						    'alt_species'         => $alt_species
						};
			
		    }
		}
	    }
	    elsif ($line =~ /^database:(\S+)\s+organism_type:(\S+)\s+include_genefinders:(\S*)\s+exclude_genefinders:(\S*)/){

		my $database      = $1;
		my $organism_type = $2;
		my $include_genefinders = $3;
		my $exclude_genefinders = $4;

		$database_type = $database ."_" .$organism_type;

		if (&verify_organism_type($organism_type, $linectr)){
		    
		    ($include_genefinders, $exclude_genefinders) = &verify_and_set_genefinders($include_genefinders, $exclude_genefinders, $linectr);

		    

		    if (( exists $hash->{$database_type}) &&
			(defined($hash->{$database_type}))){
			
			$logger->warn("This database_type '$database_type' was already encountered in the control file!");
		    }
		    else {
			#
			# Encountered this organism/database for the first time while processing the control file contents
			# therefore go ahead and declare the organism's hash attributes
			#
			$hash->{$database_type} = { 'database'            => $database,
						    'asmbl_id_list'       => [],
						    'infohash'            => [],
						    'organism_type'       => $organism_type,
						    'include_genefinders' => $include_genefinders,
						    'exclude_genefinders' => $exclude_genefinders,
						    'alt_database'        => 'none',
						    'alt_species'         => 'none'
						};
			
		    }
		}
	    }
	    elsif ($line =~ /^\s*(\d+)\s*/){

		&store_asmbl_id($database_type, $1, $line, $linectr, $unique_asmbl_id_values, $hash, $file);
	    }
	    else {
		$logger->logdie("Could not parse line number '$linectr' - line was '$line'");
	    }
	}
    }

    return $hash;

}



#------------------------------------------------------
# store_asmbl_id()
#
#------------------------------------------------------
sub store_asmbl_id {
    
    my ($database, $asmbl_id, $line, $linectr, $unique_asmbl_id_values, $hash, $file) = @_;

    if (( exists $unique_asmbl_id_values->{$database}->{$asmbl_id}) && 
	(defined($unique_asmbl_id_values->{$database}->{$asmbl_id}))){
	
	$logger->logdie("Already processed information for asmbl_id '$asmbl_id' organism '$database'.  Please review legacy2bsml control file '$file'.");
    }
    else {

	my $sequence_type = 'none';
	
	my $tu_list_file = 'none';

	my $model_list_file = 'none';

	my @attributes = split(/\s+/, $line);
	
	
	foreach my $attribute (@attributes){

	    if ($attribute =~ /:/){

		my ($key, $value) = split(/:/, $attribute);
		
		if ($key eq 'sequence_type'){
		    $sequence_type = &verify_and_set_sequence_type($value, $linectr);
		}
		elsif ($key eq 'tu_list_file'){
		    $tu_list_file = &verify_and_set_tu_list_file($value);
		}
		elsif ($key eq 'model_list_file'){
		    $model_list_file = &verify_and_set_tu_list_file($value);
		}
		else {
		    $logger->warn("Unrecognized attribute");
		}
	    }
	}

	## Store the next valid asmbl_id in the organism hash
	my $infohash = { 'asmbl_id'      => $asmbl_id,
			 'sequence_type' => $sequence_type,
			 'tu_list'       => $tu_list_file,
			 'model_list'    => $model_list_file 
		     };
	
	push( @{$hash->{$database}->{'infohash'}}, $infohash);

	## Update the unique_asmbl_id_values hash to ensure unique assembly identifiers are processed
	$unique_asmbl_id_values->{$database}->{$asmbl_id}++;

    }
}


#-----------------------------------------------------
# verify_and_set_tu_list_file()
#
#-----------------------------------------------------
sub verify_and_set_tu_list_file {

    my ($file) = @_;

    if (-e $file){
	if (-f $file){
	    if (-r $file){
		if (!-z $file){
		    return $file;
		}
		else {
		    $logger->logdie("file '$file' had zero content");
		}
	    }
	    else {
		$logger->logdie("file '$file' does not have read permissions");
	    }
	}
	else {
	    $logger->logdie("file '$file' is not a regular file");
	}
    }
    else {
	$logger->logdie("file '$file' does not exist");
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


sub verify_alt_database_and_species {
    
    my ($alt_database, $alt_species) = @_;
 
    if ((!defined($alt_database)) ||
	($alt_database =~ /^\s*$/)){
	
	$alt_database = "none";
    }

    if ((!defined($alt_species)) ||
	($alt_species =~ /^\s*$/)){
	
	$alt_species = "none";
    }

    return ($alt_database, $alt_species);
}
