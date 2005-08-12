#!/usr/local/bin/perl
#--------------------------------------------------------------------------------------------
#
# cvs: ANNOTATION/cram/workflow/generate_asmbl_list.pl
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

=back

=head1   DESCRIPTION

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

use Workflow::Logger;
use File::Basename;

umask(0000);

my %options = ();

my $results = GetOptions (\%options, 
                          'asmbl_files|f=s', 
                          'asmbl_list|a=s',
			  'output|o=s',
			  'log|l=s',
                          'debug=s', 
                          'help|h' ) || pod2usage();

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Workflow::Logger::get_logger();


my $keyname = '$;ASMBL_ID$;';
my $subflowname = '$;SUBFLOW_NAME$;';

my $iteratorconf = {
    $keyname      => [],
    $subflowname  => [],
};



if ((!defined($options{'asmbl_files'})) and (!defined($options{'asmbl_list'}))){
    $logger->logdie("Either asmbl_files or asmbl_list should be defined");
}

#
# Both asmbl_files and asmbl_list cannot be == none
#
if ( ( lc($options{'asmbl_files'}) eq 'none' )  and ( lc($options{'asmbl_list'}) eq 'none' ) ) {
    $logger->logdie("Either asmbl_files or asmbl_list should be defined");
}

if (!defined($options{'output'})){
    $logger->logdie("output was not defined.");
}



#
# Read in the information from asmbl_file OR asmbl_list
#
if(($options{'asmbl_files'}) && (uc ($options{'asmbl_files'}) ne 'NONE')){
    &get_list_from_file($iteratorconf,$options{'asmbl_files'});
}

#
# editor:  sundaram@tigr.org
# date:    2005-08-04
# bgzcase: 2025
# comment: The script should not process asmbl_list == none
#
if(($options{'asmbl_list'}) && ( uc ($options{'asmbl_list'}) ne 'NONE')){
    &get_list_from_list($iteratorconf,$options{'asmbl_list'});
}


#
# Output the lists
#
if($options{'output'}){
    &output_lists($iteratorconf, $options{'output'});
}
else{
    $logger->logdie("output was not defined");
}

exit;
						     


#---------------------------------------------------------------------------------------------------------
#
#                           END OF MAIN  --  SUBROUTINES FOLLOW
#
#---------------------------------------------------------------------------------------------------------


sub output_lists {

    my ($iteratorconf, $output) = @_;

    open FILE, "+>$output" or $logger->logdie("Can't open output file $output");
    foreach my $key (keys %$iteratorconf){
	print FILE "$key=",join(',',@{$iteratorconf->{$key}}),"\n";
    }
    close FILE;

}



sub get_list_from_file{
    my ($iteratorconf, $f) = @_;

    my @lines;
    my @files = split(',',$f);
    foreach my $file (@files){
	if( $file){
	    open( FH, $file ) or die "Could not open $file";
	    my @contents = <FH>;
	    close( FH );

	    chomp @contents;

	    foreach my $asmbl_id (@contents){
		$asmbl_id =~ s/\s//g;
		push (@lines, $asmbl_id);
	    }

	    fisher_yates_shuffle(\@lines);
	    foreach my $asmbl_id (@lines){
		if($asmbl_id){
		    my $subflow_name = 'asmbl_id_' . $asmbl_id;
		    &add_entry_to_conf($iteratorconf,$asmbl_id, $subflow_name);
		}
	    }

	}
    }
}

sub get_list_from_list{

    my ($iteratorconf, $list) = @_;

    my @array = split(',',$list);

    
    fisher_yates_shuffle(\@array);
    foreach my $asmbl_id (@array){
	if($asmbl_id){
	    my $subflow_name = 'asmbl_id_' . $asmbl_id;
	    &add_entry_to_conf($iteratorconf,$asmbl_id, $subflow_name);
	}
    }
}

sub add_entry_to_conf{
    my($iteratorconf,$filename,$name) = @_;
    push( @{$iteratorconf->{$keyname}}, $filename );
    push( @{$iteratorconf->{'$;SUBFLOW_NAME$;'}}, $name );
}

sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
	my $j = int rand ($i+1);
	next if $i == $j;
	@$array[$i,$j] = @$array[$j,$i];
    }
}
