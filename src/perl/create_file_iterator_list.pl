#!/usr/bin/perl -w

=head1  NAME 

create_file_iterator_list.pl - generates list file from various input sources of filenames.

=head1 SYNOPSIS

USAGE: create_file_iterator_list.pl 
        --input_file_list=/path/to/some.list
        --input_file=/path/to/somefile.fsa
        --input_directory=/path/to/some/dir
        --input_directory_extension=fsa
        [--checksum_filenames=1
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

B<--output_iterator_list,-o>
    output file. Up to 7 lines with comma separated lists of limits

B<--output_iterator_list,-o>
    comma separated list of iterator keys used in the output file.  Each key will be a separate line in the output file.  There are 7 lines in the output file listing FILE, FILE_NAME, FILE_BASE, FILE_EXT, DIRECTORY, XML file

B<--checksum_filenames> 
    use checksums instead of the basename as the iterator name. The checksum will be based on the full path to the file.

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
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Ergatis::Logger;
use Digest::MD5 qw(md5_hex);

my %options = ();
my $results = GetOptions (\%options, 
			  'output_iter_list|o=s',
			  'input_file|f=s',
			  'input_file_list|l=s',
			  'input_directory|d=s',
			  'input_directory_extension=s',
			  'timestamp|t=s',
              'checksum_filenames=s',
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

### to prevent too much data from being loaded into memory, we store unique
##  directories where with IDs and associate the directory id with the file name
##  this means file names can't be duplicated within an input set (even if
##  they are in separate directories.)
my $id2dir = {};
my $dir2id = {};
my $next_directory_id = 0;

# structure like @$input_elements['somename.fsa',directory id, ext ]
# needs to be an array so input order is preserved. 
# previous version used a hash. an array simplies the code a great deal anyway
my $input_elements = [];
#keep track of basename and die if multiple files with same basename
my $elements_check = {};

$logger->debug("gathering input") if $logger->is_debug();
my $input_element_count = &gather_input_elements( $input_elements );

open(my $out_fh, ">$options{output_iter_list}") || die "can't create output_iter_list: $!";

print $out_fh '$;I_FILE_BASE$;',"\t",'$;I_FILE_NAME$;',"\t",'$;I_FILE_PATH$;',"\t",'$;I_FILE_EXT$;',"\t",'$;I_DIR$;';

my $timestamp;
if($options{'timestamp'}){
    print $out_fh "\t",'$;I_TIMESTAMP$;';
    $timestamp = &get_datetime();
}
print $out_fh "\n";

for my $elt (@$input_elements){
    my $path1 = ($elt->[2] eq '') ? 	"$elt->[0]" : "$elt->[0].$elt->[2]";
    my $path2 = ($elt->[2] eq '') ? 	"$id2dir->{$elt->[1]}/$elt->[0]" : "$id2dir->{$elt->[1]}/$elt->[0].$elt->[2]";
    my $field1 = $options{'checksum_filenames'} ? $elt->[3]: $elt->[0];
    print $out_fh 
	"$field1\t",
	"$path1\t",
	"$path2\t",
	"$elt->[2]\t",
	"$id2dir->{$elt->[1]}";
    print $out_fh "\t$timestamp" if($options{'timestamp'});
    print $out_fh "\n";	
}

close $out_fh;

exit;

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
    if ( $$file =~ m|^(.+)/(([^/]+)\.([^\.\d]+))$| ) {
        return [ $1, $2, $3, $4 ];

    ## match here if the file doesn't have an extension
    } elsif ( $$file =~ m|^(.+)/(([^/]+))$| ) {
        return [ $1, $2, $3, '' ];
        
    ## else die
    } else {
        $logger->logdie("failed to extract file name parts from file: $$file");
    }
}

sub add_element {
    my ($elem, $input_elements) = @_;
    
    my $parts = parse_file_parts( \$elem );
    my $checksum = md5_hex $elem;
    my $directory_id;

    ## we can't have encountered this name already
    if($options{'checksum_filenames'}) {
        if(exists $$elements_check{ $checksum }) {
            $logger->logdie("found duplicate file in input set: $elem");
        }
    }
    elsif ( exists $$elements_check{ $$parts[2] } ) {
        $logger->logdie("found duplicate basename in input set: $$parts[2]");
    }
    $elements_check->{$parts->[2]}=1;
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
    push @$input_elements, [$$parts[2],$directory_id, $$parts[3],$checksum ];
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

	    if (($input_file_list =~ /\.gz$/) || ($input_file_list =~ /\.gzip$/)){
		die "Support unavailable for processing gzipped file lists";
	    }

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

sub check_parameters {
    
    ## handle required arguments
    my @required = qw( output_iter_list );
    for ( @required ) {
        unless ( $options{$_} ) {
            print STDERR "--$_ is a required option\n\n";
            exit(1);
        }
    }
    
    ## at least one input type is required
    unless ( $options{input_file_list} || 
             $options{input_file} || 
             $options{input_directory} ) {

        print STDERR "no input defined, please read perldoc $0\n\n";
        exit(1);
    }

    if($options{input_directory} && $options{input_directory} !~ /\/$/){
	$options{input_directory} .= "/";
    }
    
}

#---------------------------------------------
# get_sybase_datetime()
#
#---------------------------------------------
sub get_datetime {


# perl localtime   = Tue Apr  1 18:31:09 2003
# sybase getdate() = Apr  2 2003 10:15AM

    my $datetime = localtime;
    #                  Day of Week                        Month of Year                                       Day of Month  Hour      Mins     Seconds    Year   
    if ($datetime =~ /^(Mon|Tue|Wed|Thu|Fri|Sat|Sun)[\s]+(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\s+([\d]{1,2})\s+([\d]{2}):([\d]{2}):[\d]{2}\s+([\d]{4})$/){
	my $hour = $4;
	my $ampm = "AM";
	if ($4 ge "13"){
	    $hour = $4 - 12;
	    $ampm = "PM";
	}
	$datetime = "$2  $3 $6  $hour:$5$ampm";
    }
    else{
	$logger->logdie("Could not parse datetime");
    }

    return $datetime;
    
}
