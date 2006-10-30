#!/usr/local/bin/perl

use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

generate_groups.pl

=head1 SYNOPSIS

Reads a .list file to create groups of work

=head1 USAGE EXAMPLE

generate_groups.pl --file=file.list --group_count=150 --output_dir=/some/path --prefix=whatever --debug debug_level --log log_file

=head1 OPTIONS

B<--file,-f> Input file list.  This is usually the output of generate_input_list, which creates files named like "something.list".

B<--group_count,-c> The number of groups to create when splitting the work available.

B<--groupsize,-g> Deprecated (misnomer).  Use --group_count instead.

B<--output_dir,-o> Directory where the output list and groups files should be created.

B<--prefix,-p> Because the root portion of the names of all output files.  See OUTPUT section.

B<--debug,-d> Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> Log file

B<--help,-h> This help message

=head1   DESCRIPTION

This script reads a .list file, usually created by the generate_input_list.pl script,
which contains a list of the files that need to be grouped.  It uses this list, and
the --group_count option to first create another .list file that contains the names
and paths of each of the groups created, and then the group files themselves, named
like $prefixN, where $prefixed is passed with the --prefix option, and N is an 
incrementing group number, depending on how many groups are created. 

=head1   INPUT

The input file, usually a .list created by generate_input_list.pl, must be formatted
like:

    $;BSML_FILE$;=/path/to/file1.bsml,/path/to/file2.bsml,...
    $;SUBFLOW_NAME$;=file1,file2,...

The variable names here, such as BSML_FILE and SUBFLOW_NAME, don't matter.  All 
that is required is that there are variables in this format whose values are 
comma-separated lists.

=head1   OUTPUT

The --prefix option controls a portion of the output file names.  One file will
always be created called:

    $prefix.list

This file is similar to the input list file, but is formatted like:
    
    $;GROUP_NAME$;=1,0,2
    $;SUBFLOW_NAME$;=$prefix1,$prefix2,...
    $;GROUP_FILE$;=/path/to/$prefix1,/path/to/$prefix2,...

The $;GROUP_FILE$; parameter holds a comma-separated list of the paths to each
of the group files created (format explained below). The $;SUBFLOW_NAME$;
parameter holds *names* for each of the group_file elements so that they can
be referred to easily and have nice little names in the workflow viewers.

Each group file created is named like:

    $prefixN

Where 'N' is an incrementing integer for each group file created.  These files 
are formatted like:

    $;BSML_FILE$;=/path/to/file1.bsml,/path/to/file2.bsml,...
    $;SUBFLOW_NAME$;=file1,file2,...

Notice that this looks exactly like the original list file.  This is because it
is meant to be.  The point of this script is to split that original list file
into smaller groups, so the output group files will be exactly input list file
but broken up into sections of work.

=cut


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Ergatis::Logger;
use IO::File;

umask(0000);

my %options = ();
my $results = GetOptions (\%options,
			  'file|f=s',
			  'output_dir|o=s',
			  'prefix|p=s',
			  'log|l=s',
              'group_count|c=i',
			  'groupsize|g=s',
			  'debug=s',
			  'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Ergatis::Logger::get_logger();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDERR} );
}


&check_parameters(\%options);

## map 'groupsize' to 'group_count' if it is still being passed.
##  this can be removed once 'groupsize' is no longer supported
if ($options{groupsize} && (! $options{group_count})) {
    $options{group_count} = $options{groupsize};
} 

my $eltshash = &readiteratorconf($options{'file'});

my $groupsconf = {'$;GROUP_NAME$;' => [],
                  '$;GROUP_FILE$;' => [],
                  '$;SUBFLOW_NAME$;' => []};


foreach my $group (keys %$eltshash){
    #the file generated here is used to produce the subflows within each group
    mkdir("$options{'output_dir'}/$group");
    my $filename = "$options{'output_dir'}/$group/$options{'prefix'}$group";
    my $fh = IO::File->new("+>$filename") or die "Can't open $filename for writing due to $!\n";
    my $size;
    foreach my $elt (keys %{$eltshash->{$group}}){
        $size = scalar(@{$eltshash->{$group}->{$elt}});
        print $fh "$elt=",join(",",@{$eltshash->{$group}->{$elt}}),"\n";
    }
    my $grouplist = "$group,"x$size;
    chop $grouplist;
    print $fh '$;GROUP_NAME$;=' . "$grouplist\n";
    close $fh;
    push( @{$groupsconf->{'$;GROUP_NAME$;'}}, $group );
    push( @{$groupsconf->{'$;GROUP_FILE$;'}}, $filename );
    push( @{$groupsconf->{'$;SUBFLOW_NAME$;'}}, "$options{'prefix'}$group" );
}

#the file $listfile is used to produce the groups.xml file
my $listfile = "$options{'output_dir'}/$options{'prefix'}.list";
my $lfh = IO::File->new("+>$listfile") or die "Can't open $listfile for writing due to $!\n";

foreach my $key (keys %$groupsconf){
    print $lfh "$key=",join(',',@{$groupsconf->{$key}}),"\n";
}
close $lfh;

sub readiteratorconf{
    my($listfile,$numelts) = @_;

    my $eltshash;

    open FILE, $listfile or $logger->logdie("Can't open file $listfile");

    my $prevnumelts=-1;

    while(my $line=<FILE>){
	chomp $line;
	if($line =~ /=/){
	    my($key,$value)=split(/=/,$line);
        my @elts = split(/,/,$value);
	    
	    if(scalar(@elts) == 0){
		$logger->logdie("Invalid value '$value' for key $key in $listfile");
	    }

	    my $numelts;
	    if((scalar(@elts) % $options{'group_count'})==0){
		$numelts = int(scalar(@elts)/$options{'group_count'});
	    }
	    else{
		$numelts = int(scalar(@elts)/$options{'group_count'})+1;
	    }
	    
	    $logger->debug("Size of groups set to $numelts") if($logger->is_debug());
	    if($prevnumelts != -1){
		    if($numelts != $prevnumelts){
		        $logger->fatal("Number of elements for $key did not match previous count $prevnumelts");
		    }
	    }
	    $prevnumelts = $numelts;
	    
        $logger->debug("Splitting ",scalar(@elts)," elements for key $key") if($logger->is_debug());
	    
        ## for each numbered group
        for (my $gnum=0; scalar @elts; $gnum++) {
            ## grab all the members of this group
            for ( 1 .. $numelts ) {
        
                ## initialize this key if needed
                if(!$eltshash->{$gnum}->{$key}) {
                    $eltshash->{$gnum}->{$key} = [];
                }

                ## add this value of the key
                push @{$eltshash->{$gnum}->{$key}}, shift @elts;
                
                ## if groups are uneven, the last group won't be full
                last unless (scalar @elts);
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
