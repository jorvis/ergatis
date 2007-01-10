#!/usr/local/bin/perl -w

=head1  NAME 

create_iterator_groups.pl - generates groups of iterator list files from various input sources.

=head1 SYNOPSIS

USAGE: create_iterator_groups.pl 
        --input_file_list=/path/to/some.list
        --input_file=/path/to/somefile.fsa
        --input_directory=/path/to/some/dir
        --input_directory_extension=fsa
        --group_count=150
        --output_directory=/path/to/some/dir
      [ --group_size_limit=1000
        --log=/path/to/some.log
        --debug=4 ]

=head1 OPTIONS

B<--input_file_list,-l> 
    a plain text file containing the full paths to any number of files, one per line.  
    this can also be a comma-separated list of input file lists.

B<--input_file,-f> 
    the full path to an input file. this can also be a comma-separated list of 
    input files.

B<--input_directory,-d> 
    the full path to an input directory. this can also be a comma-separated list of 
    input directories.

B<--input_directory_extension,-x> 
    to be used in conjuction with the input_directory option, this can be used to
    filter files by extension within any passed input directories.

B<--group_count,-g> 
    all of the input will be separated into, as a target, this many primary groups.
    more groups will be created if the group_size_limit is reached.

B<--group_iter_file,-r> 
    iteration file created for Workflow that summarizes all the groups created.

B<--output_directory,-o> 
    root directory where output groups will be written.

B<--group_size_limit,-z> 
    optional.  all of the input will be separated into, at most, this many primary 
    groups.  (default is 1000)

B<--log> 
    optional.  path to a log file the script should create.  will be overwritten if
    already exists.
    
B<--debug> 
    optional.  the debug level for the logger (an integer)

B<--help,-h> 
    This help message/documentation.

=head1   DESCRIPTION

    This script is used to accept a selection of inputs from either an input list, 
    directory, file, or any combination thereof.  Each of these options can also
    be specified using comma-separated lists.  These inputs will then be distributed
    randomly across a certain number of groups, specified using the --group_count
    option (usually somewhere around 150).
    
    To prevent too many files from being directly in one group, additional groups are 
    automatically created for large input sets.  You can control how many inputs are 
    in a group before another is created using the --group_size_limit option (default 
    1000).  Once the number of inputs in any group reaches this number, another group 
    is created until it is full, and so on.

=head1   OUTPUT

    The location of the output is specified using the --output_directory option.  One
    numbered directory will be created under this for each group needed.
    
    MORE OUTPUT DESCRIPTIONS NEEDED
    
=head1 CONTACT

    Joshua Orvis
    jorvis@tigr.org

=cut

use lib (@INC,$ENV{"PERL_MOD_DIR"});
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev );
use Pod::Usage;
use Ergatis::Logger;

my %options = ();
my $results = GetOptions (\%options, 
              'input_iter_list|i=s',
              'group_count|g=s',
              'output_iter_list|r=s',
              'group_size_limit|z=i',
              'output_directory|o=s',
              'randomize|s=s',
              'log=s',
              'debug=s',
              'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                  'LOG_LEVEL'=>$options{'debug'});
$logger = Ergatis::Logger::get_logger();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

## play nicely
umask(0000);

## make sure all passed options are peachy
&check_parameters(\%options);

## make sure the output directory exists
if (! -d $options{output_directory} ) {
    $logger->debug("attempting to create $options{output_directory}") if $logger->is_debug();
    mkdir($options{output_directory}) || $logger->logdie("failed to create output directory: $!");
}

#skim is list of byte offsets for each line in the file
my $input_skim = [];
#first line of input_iter_list is comma seperated list of key names
my $iterator_names = [];

$logger->debug("gathering input") if $logger->is_debug();

#returns list of byte offsets. don't parse lines to save space
&skim_iterator_list( $options{input_iter_list}, $input_skim, $iterator_names );
my $input_element_count = scalar(@$input_skim);
#shuffle lines
&fisher_yates_shuffle($input_skim) if($options{'randomize'});

## calculate the number of groups we'll need to have.  
my $group_count = 0; 

if($options{group_count} =~ /(\d+)\+/){
    $options{group_count} = $1;
    my $sgestatus = `qstat -q default.q -f | grep lx26`;
    my @lines = split(/\n/,$sgestatus);
    my $totalslots;
    my $totalrunning;
    foreach my $line (@lines){
	  my @x = split(/\s+/,$line);
	  my($used,$tot) = ($x[2] =~ /(\d+)\/(\d+)/);
	  $totalslots+=$tot;
	  $totalrunning+=$used;
      }
    my $free = $totalslots-$totalrunning;
    print STDERR "Free $free\n";
    if($free > $options{group_count}){
	$options{group_count} = $free;
    }
    print STDERR "Using $options{group_count}\n";
}

## is there 1 or 0 elements per group?
if ( $input_element_count <= $options{group_count} ) {
    $group_count = $input_element_count;

## else will the requested group count hold everything, with 1 or more elements per group?
} elsif ( $input_element_count <= ( $options{group_count} * $options{group_size_limit} ) ) {
    $group_count = $options{group_count};
    
} else {    ## there are going to be more groups than the user requested.
    ## find the first group count where the sizes are less than the group_size_limit
    $group_count = $options{group_count} + 1;
    
    while ( 1 ) {
        if ( $options{group_size_limit} >= ( $input_element_count / $group_count ) ) {
            last;
        } else {
            $group_count++;
        }
    }
}

$logger->debug("actual group count was $group_count") if $logger->is_debug;

## make sure we can create the output_iter_list
$logger->debug("attempting to create output_iter_list: $options{output_iter_list}") if $logger->is_debug;

open(my $gif_fh, ">$options{output_iter_list}") || die "can't create output_iter_list: $!";
print $gif_fh '$;GROUP_XML$;',"\t",'$;ITERATOR_LIST$;',"\t",'$;GROUP_NUMBER$;',"\n";
open(my $iter_file, "$options{input_iter_list}") or $logger->logdie("Can't open file $options{input_iter_list}");

##########
## distribute the files across the groups
my @groups;
my $current_group_num = 1;
#assign each line to a group
for my $elem (@$input_skim){
    
    ## add it to the current group and get ready to add to the next group
    push @{ $groups[$current_group_num - 1] }, $elem;
    $current_group_num == $group_count ? $current_group_num = 1 : $current_group_num++;
}
$logger->debug("Created ".scalar(@groups)." group\n") if $logger->is_debug;
##########
## now we can dump the groups to files
my $group_num = 1;
for my $group ( @groups ) {

    ## make sure the group directory exists
    if (! -d "$options{output_directory}/g$group_num" ) {
        $logger->debug("attempting to create $options{output_directory}/g$group_num") if $logger->is_debug();
        mkdir("$options{output_directory}/g$group_num") || $logger->logdie("failed to create group directory: $!");
    }
    
    ## open this group file
    my $group_file = "$options{output_directory}/g$group_num/g$group_num.iter";
    $logger->debug("attempting to create and populate $group_file") if $logger->is_debug();
    open(my $group_fh, ">$group_file") || $logger->logdie( "failed to create group file $group_file\n" );
    
    ## print iterator names
    print $group_fh join("\t",@$iterator_names),"\t",'$;GROUP_NUMBER$;',"\n";
    ## elements of the group
    for my $elt ( @$group ) {
	## seek to byte offset for the line
	seek($iter_file,$elt,0);
	my $line = <$iter_file>;
	chomp $line;
	print $group_fh $line,"\t$group_num\n";
    }
    close $group_fh;
    print $gif_fh "g$group_num/g$group_num\t$group_file\t$group_num\n";
    $group_num++;
}

close $gif_fh;

exit(0);

sub skim_iterator_list{
    my($file,$iterator_elts,$iterator_names) = @_;
    my $linenum = 0;
    open FILE, "$file" or $logger->logdie("Can't open file $file");
    my $currpos=0;
    while(my $line=<FILE>){
	chomp $line;
	if($linenum==0){
	    @$iterator_names = split(/\t/,$line);
	}
	elsif($line =~ /\S/){
	    push @$iterator_elts,$currpos;
	}
	$linenum++;
	$currpos = tell(FILE);
    }
    close FILE;
}

sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    if(@$array){
        for ($i = @$array; --$i; ) {
            my $j = int rand ($i+1);
            next if $i == $j;
            @$array[$i,$j] = @$array[$j,$i];
        }
    }
}

sub check_parameters {
    
    ## handle required arguments
    my @required = qw( output_directory group_count output_iter_list );
    for ( @required ) {
        unless ( $options{$_} ) {
            print STDERR "--$_ is a required option\n\n";
            exit(1);
        }
    }
    
    ## make sure group count is greater than 0
    unless ( $options{group_count} =~ /\d+\+/ || $options{group_count} > 0 ) {
        print STDERR "--group_count must be greater than 0\n\n";
        exit(1);
    }
    $options{group_size_limit} = 1000 if(!$options{group_size_limit});

    ## make sure group count is greater than 0
    unless ( $options{group_size_limit} > 0 ) {
        print STDERR "--group_size_limit must be greater than 0\n\n";
        exit(1);
    }
}

