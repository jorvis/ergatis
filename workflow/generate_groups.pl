#!/usr/local/bin/perl

=head1  NAME 

dummy.pl - do nothing

=head1 SYNOPSIS

USAGE:  dummy.pl --debug debug_level --log log_file

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
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Workflow::Logger;
use IO::File;

my %options = ();
my $results = GetOptions (\%options,
			  'file|f=s',
			  'output_dir|o=s',
			  'prefix|p=s',
			  'log|l=s',
			  'groupsize|g=s',
			  'debug=s',
			  'help|h') || pod2usage();

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Workflow::Logger::get_logger();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDERR} );
}


&check_parameters(\%options);

my $eltshash = &readiteratorconf($options{'file'});

my $groupsconf = {'$;GROUP_FILE$;'=>[],
		  '$;SUBFLOW_NAME$;'=>[]};


foreach my $group (keys %$eltshash){
    my $filename = "$options{'output_dir'}/$options{'prefix'}$group";
    my $fh = IO::File->new("+>$filename") or die "Can't open $filename for writing due to $!\n";
    foreach my $elt (keys %{$eltshash->{$group}}){
	print $fh "$elt=",join(",",@{$eltshash->{$group}->{$elt}}),"\n";
    }
    close $fh;
    push( @{$groupsconf->{'$;GROUP_FILE$;'}}, $filename );
    push( @{$groupsconf->{'$;SUBFLOW_NAME$;'}}, "$options{'prefix'}$group" );
}


my $listfile = "$options{'output_dir'}/$options{'prefix'}.list";
my $lfh = IO::File->new("+>$listfile") or die "Can't open $listfile for writing due to $!\n";

foreach my $key (keys %$groupsconf){
    print $lfh "$key=",join(',',@{$groupsconf->{$key}}),"\n";
}
close $lfh;

sub readiteratorconf{
    my($listfile,$numelts) = @_;

    my $eltshash;

    open FILE, $listfile or $logger->get_logger()->logdie("Can't open file $listfile");

    my $prevnumelts=-1;

    while(my $line=<FILE>){
	chomp $line;
	if($line =~ /=/){
	    my($key,$value)=split(/=/,$line);
	    my @elts = split(/,/,$value);
	    
	    my $numelts;
	    if((scalar(@elts) % $options{'groupsize'})==0){
		$numelts = int(scalar(@elts)/$options{'groupsize'});
	    }
	    else{
		$numelts = int(scalar(@elts)/$options{'groupsize'})+1;
	    }
	    
	    $logger->debug("Size of groups set to $numelts") if($logger->is_debug());
	    if($prevnumelts != -1){
		if($numelts != $prevnumelts){
		    $logger->get_logger->fatal("Number of elements for $key did not match previous count $prevnumelts");
		}
	    }
	    $prevnumelts = $numelts;
	    $logger->debug("Splitting ",scalar(@elts)," elements for key $key") if($logger->is_debug());
	    for(my $i=-1; $i<scalar(@elts); $i+=$numelts){
	    	$logger->debug("Generating group starting at $i") if($logger->is_debug());
		for(my $j=$i+1; $j <= ($i+$numelts); $j++){
		    if($j<scalar(@elts)){
			my $name = ($i+1)."_".($i+$numelts);
			$logger->debug("At elt index $j using group name $name") if($logger->is_debug());
			
			if(!$eltshash->{$name}->{$key}){
			    $eltshash->{$name}->{$key} = [];
			}
			push @{$eltshash->{$name}->{$key}}, $elts[$j];
		    }
		}
	    }
	}
    }
    return $eltshash;
}
	   

sub check_parameters{
    my ($options) = @_;
    
    if(0){
	pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});    
    }
}
