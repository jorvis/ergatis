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
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through );
use Pod::Usage;
use lib '/usr/local/devel/ANNOTATION/ard/curring/lib/site_perl/5.8.5/';
use Ergatis::Logger;

my %options = ();
my $results = GetOptions (\%options, 
              'input_file_list|l=s',
              'input_file|f=s',
              'input_directory|d=s',
              'input_directory_extension|x=s',
              'group_count|g=i',
              'output_iter_list|r=s',
              'group_size_limit|z=i',
              'output_directory|o=s',
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

## to prevent too much data from being loaded into memory, we store unique
##  directories where with IDs and associate the directory id with the file name
##  this means file names can't be duplicated within an input set (even if
##  they are in separate directories.)
my $id2dir = {};
my $dir2id = {};
my $next_directory_id = 0;

## structure like $$input_elements{'somename.fsa'} = [ directory id, ext ]
my $input_elements = {};

$logger->debug("gathering input") if $logger->is_debug();
my $input_element_count = &gather_input_elements( $input_elements );

## calculate the number of groups we'll need to have.  
my $group_count = 0; 

## is there 1 or 0 elements per group?
if ( $input_element_count <= $options{group_count} ) {
    $group_count = $input_element_count;

## else will the requested group count hold everything, with 1 or more elements per group?
} elsif ( $input_element_count <= ( $options{group_count} * $options{group_size_limit} ) ) {
    $group_count = $options{group_count};
    
} else {
    ## there are going to be more groups than the user requested.
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
$logger->debug("found $input_element_count input elements") if $logger->is_debug;
$logger->debug("there were $next_directory_id distinct directories") if $logger->is_debug;

## make sure we can create the output_iter_list
$logger->debug("attempting to create output_iter_list: $options{output_iter_list}") if $logger->is_debug;
my @group_iter_list;
my @group_output_files;
my @group_numbers;
open(my $gif_fh, ">$options{output_iter_list}") || die "can't create output_iter_list: $!";

##########
## distribute the files across the groups
my @groups;
my $current_group_num = 1;
for my $elem ( keys %$input_elements ) {
    
    ## add it to the current group and get ready to add to the next group
    push @{ $groups[$current_group_num - 1] }, \$elem;
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
    my $group_dir = "g$group_num";
    push @group_iter_list,$group_file;
    push @group_output_files,"$group_dir/g$group_num.xml";
    push @group_numbers,$group_num;
    $logger->debug("attempting to create and populate $group_file") if $logger->is_debug();
    open(my $group_fh, ">$group_file") || $logger->logdie( "failed to create group file $group_file\n" );
    
    ## handle each variable

    ## ITER_FILE_PATH  (the full path to the file)
    print $group_fh '$;ITER_FILE_PATH$;=';
    my $a_pos = 0;
    for my $nameref ( @$group ) {
        print $group_fh $$id2dir{  $$input_elements{ $$nameref }[0]  }, "/$$nameref.$$input_elements{ $$nameref }[1]";
        print $group_fh ',' unless ( ++$a_pos == scalar @$group );
    }

    ## ITER_FILE_NAME (the name of the file including the extension, but without the path)
    print $group_fh "\n", '$;ITER_FILE_NAME$;=';
    $a_pos = 0;
    for my $nameref ( @$group ) {
        print $group_fh "$$nameref.$$input_elements{ $$nameref }[1]";
        print $group_fh ',' unless ( ++$a_pos == scalar @$group );
    }
        
    ## ITER_FILE_BASE (the name of the file without the path or extension)
    print $group_fh "\n", '$;ITER_FILE_BASE$;=';
    $a_pos = 0;
    for my $nameref ( @$group ) {
        print $group_fh "$$nameref";
        print $group_fh ',' unless ( ++$a_pos == scalar @$group );
    }
    
    ## ITER_FILE_EXT (the extension of each file)
    print $group_fh "\n", '$;ITER_FILE_EXT$;=';
    $a_pos = 0;
    for my $nameref ( @$group ) {
        print $group_fh "$$input_elements{ $$nameref }[1]";
        print $group_fh ',' unless ( ++$a_pos == scalar @$group );
    }

    ## ITER_DIR (the path to each file, directory only)
    print $group_fh "\n", '$;ITER_DIR$;=';
    $a_pos = 0;
    for my $nameref ( @$group ) {
        print $group_fh $$id2dir{  $$input_elements{ $$nameref }[0]  };
        print $group_fh ',' unless ( ++$a_pos == scalar @$group );
    }

    ## ITER_DIR (the path to each file, directory only)
    print $group_fh "\n", '$;ITER_XML$;=';
    $a_pos = 0;
    for my $nameref ( @$group ) {
        print $group_fh "$$nameref.xml";
        print $group_fh ',' unless ( ++$a_pos == scalar @$group );
    }

    ## GROUP_NUMBER
    print $group_fh "\n", '$;GROUP_NUMBER$;=';
    $a_pos = 0;
    for my $nameref ( @$group ) {
        print $group_fh "$group_num";
        print $group_fh ',' unless ( ++$a_pos == scalar @$group );
    }

    $group_num++;
    print $group_fh "\n";
}

print $gif_fh '$;ITERATOR_INI$;=';
print $gif_fh join(',',@group_iter_list),"\n";
print $gif_fh '$;GROUP_XML$;=';
print $gif_fh join(',',@group_output_files),"\n";
print $gif_fh '$;GROUP_NUMBER$;=';
print $gif_fh join(',',@group_numbers),"\n";
close $group_fh;
close $gif_fh;

exit(0);

sub add_element {
    my ($elem, $input_elements) = @_;
    
    my $parts = parse_file_parts( \$elem );
    my $directory_id;

    ## we can't have encountered this name already
    if ( exists $$input_elements{ $$parts[2] } ) {
        $logger->logdie("found duplicate basename in input set: $$parts[2]");
    }
    
    ## make sure we have an ID for this directory.  else create one
    if ( exists $$dir2id{ $$parts[0] } ) {
        $directory_id = $$dir2id{ $$parts[0] };
    } else {
        $$dir2id{ $$parts[0] } = $next_directory_id;
        $$id2dir{ $next_directory_id } = $$parts[0];
        $next_directory_id++;
        $directory_id = $$dir2id{ $$parts[0] };
        $logger->debug("keying directory $$parts[0]") if $logger->is_debug();
    }
    
    ## store the file info
    $$input_elements{ $$parts[2] } = [ $directory_id, $$parts[3] ];
}


sub gather_input_elements {
    my $input_elements = shift;
    my $element_count = 0;

    ####
    ##  all three input options can be used in tandem, so we check every one.
    ####

    ## handle any input file lists
    if ( $options{input_file_list} ) {

        for my $input_file_list ( split(',', $options{input_file_list}) ) {

            open(my $ifl_fh, $input_file_list) || $logger->logdie("failed to read $input_file_list");

            while (<$ifl_fh>) {
                next if ( /^\s*$/ );
                chomp;
                &add_element($_, $input_elements);
                $element_count++;
            }
        }
    }

    ## handle any directories
    if ( $options{input_directory} ) {

        for my $input_directory ( split(',', $options{input_directory}) ) {

            opendir(my $idh, $input_directory) || $logger->logdie("failed to read $input_directory");

            while (my $file = readdir $idh) {
                ## only consider files
                next unless ( -f "$input_directory/$file");

                ## if the extension is defined, the file must match it
                if (! $options{input_directory_extension} || $file =~ /\.$options{input_directory_extension}$/ ) {
		    &add_element($input_directory.$file, $input_elements);
                    $element_count++;
                }
            }

        }
    }

    ## handle any individual files
    if ( $options{input_file} ) {
        for my $input_file ( split(',', $options{input_file}) ) {
	    &add_element($input_file, $input_elements);
            $element_count++;
        }
    }

    return $element_count;
}

sub parse_file_parts {
    ## why?  because I couldn't get File::Basename to give me all this
    ##
    ## accepts the path to a file ( such as /path/to/somefile.txt )
    ## returns an array of values in this order:
    ##
    ## iter_dir         ( /path/to )
    ## iter_file_name   ( somefile.txt )
    ## iter_file_base   ( somefile )
    ## iter_file_ext    ( txt )
    ##
    ## if the file doesn't have an extension, the last array element
    ##  is an empty string
    ##
    my $file = shift;
    
    ## match here if the file has an extension
    if ( $$file =~ m|^(.+)/(([^/]+)\.(.+))$| ) {
        return [ $1, $2, $3, $4 ];

    ## match here if the file doesn't have an extension
    } elsif ( $$file =~ m|^(.+)/(([^/]+))$| ) {
        return [ $1, $2, $3, '' ];
        
    ## else die
    } else {
        $logger->logdie("failed to extract file name parts from file: $$file");
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
    unless ( $options{group_count} > 0 ) {
        print STDERR "--group_count must be greater than 0\n\n";
        exit(1);
    }
    
    ## at least one input type is required
    unless ( $options{input_file_list} || 
             $options{input_file} || 
             $options{input_directory} ) {

        print STDERR "no input defined, please read perldoc $0\n\n";
        exit(1);
    }
    
    ## if input_directory passed, directory extension required
#    if ( $options{input_directory} && ! $options{input_directory_extension} ) {
#        print STDERR "--input_directory_extension required when using --input_directory\n\n";
#        exit(1);        
#    }
    
    ## handle defaults
    $options{group_size_limit} = 1000 unless $options{group_size_limit};
    $options{input_directory_extension} = '' unless $options{input_directory_extension};
}






