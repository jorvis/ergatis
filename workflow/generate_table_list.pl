#!/usr/local/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";
#--------------------------------------------------------------------------------------------
#
# cvs: ANNOTATION/ergatis/workflow/generate_table_list.pl
# date: 2006-04-16
# contact: sundaram@tigr.org
#
#
# Similar to generate_input_list.pl except instead of iterating on list of BSML files,
# will iterate on list of chado tables.
#
# Context: Should run in following order:
# 1) generate_table_list.pl
# 2) generate_groups.pl
# 3) generate_subflow.pl
#
#
# $Id$
#
#
#--------------------------------------------------------------------------------------------
=head1  NAME 

generate_table_list.pl - Default output is a workflow iterator that
can be used to iterator over a set of chado tables

=head1 SYNOPSIS

USAGE:  generate_table_list

=head1 OPTIONS

=item *

B<--debug,-d> Debug level.  Use a large number to turn on verbose debugging. 

=item *

B<--commit-order> The commit order for chado schema.

=item *

B<--subcomponent> The workflow subcomponent type.

=item *

B<--output-directory> The ouptut directory.

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
                          'commit-order=s', 
			  'index-order=s',
                          'commit-order-file=s', 
			  'index-order-file=s',
			  'output-directory=s',
			  'subcomponent=s',
			  'log|l=s',
                          'debug=s', 
                          'help|h' ) || pod2usage();


#
# Create the Log4perl logger object
#
my $logger = &set_logger(\%options);

#
# Verify the arguments
#
&check_arguments(\%options);

#
#
#
my $outdir = &verify_and_set_outdir(\%options);

#
# Build the iteratorconf
#
my ($iteratorconf) = &build_iteratorconf($options{'subcomponent'}, $options{'commit-order'}, $options{'index-order'}, $options{'commit-order-file'}, $options{'index-order-file'});

#
# Output the lists
#
&output_lists($iteratorconf, $outdir, $options{'subcomponent'});


print "End of program $0\n";
exit(0);



exit;
						     
#---------------------------------------------------------------------------------------------------------
#
#                           END OF MAIN  --  SUBROUTINES FOLLOW
#
#---------------------------------------------------------------------------------------------------------


#---------------------------------------------
# verify_and_set_outdir()
#
#---------------------------------------------
sub verify_and_set_outdir {

    my ($options) = @_;

    my $outdir;

    if (( exists $options{'output-directory'}) &&
	( defined($options{'output-directory'}))) {

	$outdir = $options{'output-directory'};

	if (!-e $outdir){
	    $logger->logdie("output-directory '$outdir' does not exist");
	}
	if (!-w $outdir){
	    $logger->logdie("output-directory '$outdir' does not have write permissions");
	}
    }
    else {
	$logger->warn("output-directory was not specified and therefore was set to current working directory");
	$outdir = ".";
    }

    return $outdir;
}

#---------------------------------------------
# output_lists()
#
#---------------------------------------------
sub output_lists {

    my ($iteratorconf, $outdir, $subcomponent) = @_;

    my $outfile = $outdir . '/' . $subcomponent . '.list';

    open FILE, ">$outfile" or $logger->logdie("Can't open output file '$outfile': $!");

    foreach my $key (keys %$iteratorconf){

	print FILE "$key=",join(',',@{$iteratorconf->{$key}}),"\n";

    }

    close FILE;

}

#-----------------------------------------------------
# set_logger()
#
#-----------------------------------------------------
sub set_logger {

    my ($options) = @_;

    my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
    my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				      'LOG_LEVEL'=>$options{'debug'});

    my $logger = Ergatis::Logger::get_logger();


    return $logger;
}

#-----------------------------------------------------
# check_arguments()
#
#-----------------------------------------------------
sub check_arguments {

    my ($options) = @_;

    if (! exists $options{'subcomponent'}){
	$logger->logdie("subcomponent was not defined");
    }


    if (! exists $options{'commit-order'}){
	$logger->logdie("commit-order was not defined");
    }


    if (! exists $options{'index-order'}){
	$logger->logdie("index-order was not defined");
    }

}    



#-----------------------------------------------------
# build_iteratorconf()
#
#-----------------------------------------------------
sub build_iteratorconf {

    my ($subcomponent, $commit_order, $index_order, $commit_order_file, $index_order_file) = @_;

    my $iteratorconf = {};

    if (( $subcomponent eq 'create_tables') ||
	( $subcomponent eq 'create_views') ||
	( $subcomponent eq 'load_tables') ||
	( $subcomponent eq 'duplicate_detection_elimination') ||
	( $subcomponent eq 'grant_permissions') ||
	( $subcomponent eq 'revoke_permissions') ||
	( $subcomponent eq 'partition_tables') ) {
	
	my $iterator_key = '$;TABLE$;';

	my @commitorder;


	if ((defined($commit_order_file)) &&
	    (-e $commit_order_file) &&
	    (-r $commit_order_file)) {

	    open (INFILE, "<$commit_order_file") or $logger->logdie("Could not open file '$commit_order_file': $!");
	    
	    @commitorder = <INFILE>;
	
	    chomp @commitordert;
	}
	
	@commitorder = split(/,/, $commit_order);
		
	&load_iterator($iterator_key, \@commitorder);
	
	
    }
    elsif (( $subcomponent eq 'drop_tables') ||
	   ( $subcomponent eq 'drop_views') ||
	   ( $subcomponent eq 'dump_tables') ||
	   ( $subcomponent eq 'unpartition_tables') ){

	my $iterator_key = '$;TABLE$;';

	my @commitorder;

	if ((defined($commit_order_file)) &&
	    (-e $commit_order_file) &&
	    (-r $commit_order_file)) {

	    open (INFILE, "<$commit_order_file") or $logger->logdie("Could not open file '$commit_order_file': $!");
	    
	    @commitorder = <INFILE>;
	
	    chomp @commitordert;
	}
	
	@commitorder = split(/,/, $commit_order);
		
	my @reversecommitorder = reverse(@commitorder);

	&load_iterator($iterator_key, \@reversecommitorder);

    }
    elsif ( $subcomponent eq 'create_indices') {
	
	my @commitorder;
	
	my $iterator_key = '$;INDEX;';

	if ((defined($index_order_file)) &&
	    (-e $index_order_file) &&
	    (-r $index_order_file)) {

	    open (INFILE, "<$index_order_file") or $logger->logdie("Could not open file '$index_order_file': $!");
	    
	    @commitorder = <INFILE>;
	
	    chomp @commitordert;
	}
	
	@commitorder = split(/,/, $commit_order);
		
	&load_iterator($iterator_key, \@commitorder);



	
    }
    elsif ( $subcomponent eq 'drop_indices'){
	
	my @commitorder;
	
	my $iterator_key = '$;INDEX;';
	
	if ((defined($index_order_file)) &&
	    (-e $index_order_file) &&
	    (-r $index_order_file)) {

	    open (INFILE, "<$index_order_file") or $logger->logdie("Could not open file '$index_order_file': $!");
	    
	    @commitorder = <INFILE>;
	
	    chomp @commitordert;
	}
	
	@commitorder = split(/,/, $commit_order);

	my @reversecommitorder = reverse(@commitorder);
		
	&load_iterator($iterator_key, \@reversecommitorder);

	
    }
    else {
	$logger->logdie("Unrecognized subcomponent '$subcomponent'");
    }

    return ($iteratorconf);
}



sub load_iterator {
    
    my ($key, $commitorder, $file) = @_;

 
    foreach my $val (@{$commitorder}){
	    
	    next if ($val =~ /^\s*$/);

	    push (@{$iteratorconf->{$iterator_key}}, $val);
	}

}
