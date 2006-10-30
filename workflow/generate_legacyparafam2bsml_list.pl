#!/usr/local/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";
#
#
# $Id$
#
#
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

my $database_key       = '$;SOURCE_DATABASE$;';
my $family_id_key      = '$;FAMILY_ID$;';
my $subflow_key        = '$;SUBFLOW_NAME$;';
my $schema_key         = '$;SCHEMA_TYPE$;';
my $alt_database_key   = '$;ALT_DATABASE$;';
my $alt_species_key    = '$;ALT_SPECIES$;';


my $valid_schema_types = { 'prok'   => 1,
	                   'ntprok' => 1,
			   'euk'    => 1
			 };



if (&verify_control_file($options{'control_file'})){

    if (defined($options{'output'})){

	my $iteratorconf = {
	    $database_key      => [],
	    $family_id_key     => [],
	    $subflow_key       => [],
	    $schema_key        => [],
	    $alt_database_key  => [],
	    $alt_species_key   => []

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

    foreach my $database_type (sort keys %{$orghash} ){ 

	my $database      = $orghash->{$database_type}->{'database'};
	my $schema_type   = $orghash->{$database_type}->{'schema_type'};
	my $alt_database  = $orghash->{$database_type}->{'alt_database'};
	my $alt_species   = $orghash->{$database_type}->{'alt_species'};

	foreach my $infohash ( @{$orghash->{$database_type}->{'infohash'}} ) {

	    my $family_id = $infohash->{'family_id'};

	    my $subflow_name = $database_type . "_" . $family_id;

	    &add_entry_to_conf($iconf, $database, $family_id, $subflow_name, $schema_type, $alt_database, $alt_species);
	}
    }
    
}


#-------------------------------------------
# add_entry_to_conf()
#
#-------------------------------------------
sub add_entry_to_conf{

    my($iteratorconf,$database, $family_id, $subflow_name, $schema_type, $alt_database, $alt_species) = @_;
    
    push( @{$iteratorconf->{$database_key}}, $database );
    push( @{$iteratorconf->{$family_id_key}}, $family_id );
    push( @{$iteratorconf->{$subflow_key}}, $subflow_name );
    push( @{$iteratorconf->{$schema_key}}, $schema_type );
    push( @{$iteratorconf->{$alt_database_key}}, $alt_database );
    push( @{$iteratorconf->{$alt_species_key}}, $alt_species );

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

    my $database_type;

    my $unique_family_id_values = {};

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

	    if ($line =~ /^database:(\S+)\s+schema_type:(\S+)\s+alt_database:(\S*)\s+alt_species:(\S*)\s*$/){
		
		my $database      = $1;
		my $schema_type   = $2;

		my ($alt_database, $alt_species) = &verify_alt_database_and_species($3, $4);

		$database_type = $database ."_" . $schema_type;

		if (&verify_schema_type($schema_type, $linectr)){
		    
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
						    'family_id_list'      => [],
						    'infohash'            => [],
						    'schema_type'         => $schema_type,
						    'alt_database'        => $alt_database,
						    'alt_species'         => $alt_species
						};
			
		    }
		}
	    }
	    elsif ($line =~ /^database:(\S+)\s+schema_type:(\S+)/){

		my $database      = $1;
		my $schema_type   = $2;

		$database_type = $database ."_" .$schema_type;

		if (&verify_schema_type($schema_type, $linectr)){
		    
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
						    'family_id_list'      => [],
						    'infohash'            => [],
						    'schema_type'         => $schema_type,
						    'alt_database'        => 'none',
						    'alt_species'         => 'none'
						};
			
		    }
		}
	    }
	    elsif ($line =~ /^\s*(\d+)\s*/){

		&store_family_id($database_type, $1, $linectr, $unique_family_id_values, $hash);
	    }
	    else {
		$logger->logdie("Could not parse line number '$linectr' - line was '$line'");
	    }
	}
    }

    return $hash;

}



#------------------------------------------------------
# store_family_id()
#
#------------------------------------------------------
sub store_family_id {
    
    my ($database, $family_id, $linectr, $unique_family_id_values, $hash) = @_;


    #
    # Need to ensure that the family identifiers are unique
    #
    if (( exists $unique_family_id_values->{$database}->{$family_id}) && 
	(defined($unique_family_id_values->{$database}->{$family_id}))){
	
	$logger->warn("family_id '$family_id' for database '$database' was already stored in the lookup");
    }
    else {
	#
	# Update the unique_family_id_values hash to ensure that only unique family identifiers are processed
	#
	$unique_family_id_values->{$database}->{$family_id}++;
	
	#
	# Store the next valid family_id in the organism hash
	#
	my $infohash = { 'family_id' => $family_id };
	
	push( @{$hash->{$database}->{'infohash'}}, $infohash);
    }
}





#-----------------------------------------------------
# verify_schema_type()
#
#-----------------------------------------------------
sub verify_schema_type {

    my ($schema_type, $line) = @_;
    
    if ((exists $valid_schema_types->{$schema_type}) && 
	(defined($valid_schema_types->{$schema_type}))){
	
	return 1;
    }
    else {
	$logger->logdie("Encountered an invalid schema type '$schema_type' at line '$line'");
    }

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
