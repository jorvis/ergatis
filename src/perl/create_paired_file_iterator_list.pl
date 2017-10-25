#!/usr/bin/env perl -w

=head1  NAME 

create_file_iterator_list.pl - generates list file from various paired input sources of filenames.

=head1 SYNOPSIS

USAGE: create_file_iterator_list.pl 
        --input_pairs_list
        --input_file_list_1
        --input_file_list_2
        --input_file_1
        --input_file_2
        --all_vs_all=[0,1]
        --input_mapping_file=/path/to/file_names.map
        --output_iter_list=/path/to/output
        [--checksum_filenames=1
        --log=/path/to/some.log
        --debug=4 ]

=head1 OPTIONS

B<--input_pairs_list,-p>
    The full path to a list of file pairs. When using this option, will ignore all other 
    input options. Can be a comma separated list of lists. 

B<--input_file_list_1,-l1> 
    a plain text file containing the full paths to any number of files, one per line.  
    this can also be a comma-separated list of input file lists.

B<--input_file_list_2,-l2> 
    a plain text file containing the full paths to any number of files, one per line.  
    this can also be a comma-separated list of input file lists.

B<--input_file_1,-f1> 
    the full path to an input file. this can also be a comma-separated list of 
    input files.

B<--input_file_2,-f2> 
    the full path to an input file. this can also be a comma-separated list of 
    input files.

B<--input_mapping_file,-im>
    If input file lists are used (--input_file_list_1 and --input_file_list_2), and they
    are not in order the order specified in the input_mapping_file will be used. This file
    should be a two column tab file with input 1 in the first column and input 2 in the 
    second. Will look at file names only (no paths). Consequently, this will only work if 
    you have unique file names. If there are different numbers of files in the two input lists
    will fail. Will also fail if some files from the input list cannot be found in the mapping
    file.

B<--all_vs_all,-a>
    Possible values: 0 or non-zero.
    When using many vs many input files, these can be handled two ways. The first way is to take all combination
    of input files. The other way is to assume take them in the order they were provided. For example, the 
    first file of input_list_1 will be paired with the first file of input_list_2. Then the second files
    of the lists would be paired and so on. This program will throw an error if the number of input files for
    each set of inputs does not match and --all_vs_all = 0

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

    This script is used to accept a selection of paired inputs. This is to allow
    on-the-fly creation of n vs m iterations. If both inputs are a single
    file, there will only be one line in the output iterator list. If either (but 
    not both) of the inputs is multiple files, this script will generate n (where n
    is the number of files ) lines in the output iterator list. If both inputs contain
    multiple files, the script will look at the all_vs_all flag to determine how to 
    handle the input. If the value of the option is 1, there will be n * m output lines.
    If the value of all_vs_all = 0, then the script will assume that the inputs are ordered
    and the inputs should contain the same number of files. There will be n = m lines in 
    the output list. 

    Alternatively, if the two sets of inputs used are not in the correct order and the 
    input_mapping_file option is specified, the pairings will be made according to the
    mapping file. The mapping file should contain file names only (no path information) in
    two columns, on representing each input. If this option is used, there must be values 
    specified in the --input_file_[1,2] and/or --input_file_list_[1,2] options and they 
    must contain the same number of files. Also, the --all_vs_all option is ignored in this
    case. If some input files could not be found in the mapping file, the script will die.
    If file names are not unique, the script will also die.

    The --input_pairs_list can also be used. This is expected to be a 2 column list of paired
    input files. If this option is used, all other input options will be ignored. The output list
    will contain the same number of lines as the --input_pairs_list file does.

    
=head1   OUTPUT

    Will print an output tab file (iterator list file) with the following columns:

    $;I_FILE_BASE_1$;	$;I_FILE_NAME_1$;	$;I_FILE_PATH_1$;	$;I_FILE_EXT_1$;	$;I_DIR_1$;  $;I_FILE_BASE_2$;	$;I_FILE_NAME_2$;	$;I_FILE_PATH_2$;	$;I_FILE_EXT_2$;	$;I_DIR_2$;

    The number of lines in the output files will depend on the nature of the input. See DESCRIPTION 
    for more details.
    
=head1 CONTACT

    Kevin Galens
    kevingalens@gmail.com

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Ergatis::Logger;
use Digest::MD5 qw(md5_hex);
use File::Basename;

use Data::Dumper;

##############################
my $logger;
##############################

## play nicely
umask(0000);

my %options = ();
my $results = GetOptions (\%options, 
						  'output_iter_list|o=s',
						  'input_pairs_list|p=s',
						  'input_file_1|f1=s',
						  'input_file_list_1|l1=s',
						  'input_file_2|f2=s',
						  'input_file_list_2|l2=s',
						  'input_mapping_file|im=s',
						  'all_vs_all|a=s',
						  'checksum_filenames|c=s',
						  'log=s',
						  'debug=s',
						  'help|h') || pod2usage();


# Make sure all passed options valid.
&check_parameters(\%options);

open( my $outfh, "> $options{'output_iter_list'}") or 
  die("Can't open $options{'output_iter_list'} for writing: $!");
&create_iter_list( \%options, $outfh );
close( $outfh );

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
	my @file_parts = ('','','','','');
    
    ## match here if the file has an extension
    if ( $file =~ m|^(.+)/(([^/]+)\.([^\.]+))$| ) {
	  @file_parts = ( $3, $2, $file, $4, $1 );
    ## match here if the file doesn't have an extension
    } elsif ( $file =~ m|^(.+)/(([^/]+))$| ) {
	  @file_parts = ( $3, $2, $file, '', $1 );
    } 

	$file_parts[0] = md5_hex $file_parts[2] if( $options{'checksum_filenames'} );

	return @file_parts;
}

sub create_iter_list {
  my ($opts, $outfh) = @_;
  my @pairs;

  # If we have input pairs list, we don't care about anything else
  if( $opts->{'input_pairs_list'} ) {
	open(IN, "< $opts->{'input_pairs_list'}" ) or die("Can't open $opts->{'input_pairs_list'}: $!");
	while( my $line = <IN> ) {
	  chomp( $line );
	  push( @pairs, [ split( /\t/, $line ) ] );
	}
	close(IN);

  } else {

	my @inputs1 = &_get_inputs( $opts->{'input_file_list_1'}, $opts->{'input_file_1'} );
	my @inputs2 = &_get_inputs( $opts->{'input_file_list_2'}, $opts->{'input_file_2'} );

	if( $options{'input_mapping_file'} ) {
	  # Check to make sure we have the same number of inputs
	  die("There were different number of inputs for 1 (".scalar(@inputs1).") and 2 ("
		  .scalar(@inputs2)."). There should be an equal number of files for each input ".
		  "when the --input_mapping_file option is used") unless( @inputs1 == @inputs2 );

	  # Will check for duplicate file names in mapping file.
	  my %mappings = &parse_mapping_file( $options{'input_mapping_file'} );

	  @pairs = &order_inputs_by_mapping( \@inputs1, \@inputs2, \%mappings );


	} else {

	  if ( $opts->{'all_vs_all'} || @inputs1 == 1 || @inputs2 == 1 ) {
		@pairs = &_get_all_combinations( \@inputs1, \@inputs2 );
	  } else {
		die("Different number of inputs from each list") unless( @inputs1 == @inputs2 );
		
		for ( my $i = 0; $i < @inputs1; $i++ ) {
		  push( @pairs, [$inputs1[$i], $inputs2[$i]] );
		}
	  }

	}
	  
  }

  &print_header( $outfh );
  my $seen_pairs = {};
  map { &print_output_pair( $outfh, $_, $seen_pairs ) } @pairs;
}

sub order_inputs_by_mapping {
  my ($inputs1, $inputs2, $mapping) = @_;

  my $list_to_hash = sub {
	my $h = {};
	foreach my $i ( @_ ) {
	  my $b = basename( $i );
	  die("Found duplicate file name in input set ($b)") if( exists( $h->{$b} ) );
	  $h->{$b} = $i;
	}
	return $h;
  };

  my $file_names_1 = $list_to_hash->( @{$inputs1} );
  my $file_names_2 = $list_to_hash->( @{$inputs2} );
  
  print Dumper( $file_names_1 );

  my @pairs;
  while( my ($i1,$i2) = each( %{$mapping} ) ) {
	next unless( exists( $file_names_1->{$i1} ) && exists( $file_names_2->{$i2} ) );
	push(@pairs, [$file_names_1->{$i1},$file_names_2->{$i2}]);
	delete( $file_names_1->{$i1} );
	delete( $file_names_2->{$i2} );
  }
	
  die("There were some input files which were not specified in the mapping file. ".
	  "input 1: [".join(", ", keys %{$file_names_1} )."], input 2: ["
	  .join(", ", keys %{$file_names_2} )."]") 
	if( %{$file_names_1} || %{$file_names_2} );

  return @pairs;
	
}

sub parse_mapping_file {
  my ($mapping_file) = @_;

  ## We'll return a hash with lookups from input 1 files to input 2 files.
  my %retval;

  ## But we also need to make sure that input 2 files are unique.
  my %second_column;

  ## Line number
  my $lno = 0;

  open(IN, "< $mapping_file") or die("Can't open $mapping_file: $!");
  while( my $line = <IN> ) {
	next if( $line =~ /^\#/ || $line =~ /^\s*$/ );
	chomp $line;
	my @c = split(/\s+/, $line);
	die("Mapping file line [$lno] had ".scalar(@c)." columns [@c], expected 2") unless( @c == 2 );
	die("Found duplicate file name in mapping file [$c[0]]") if( exists( $retval{$c[0]} ) );
	die("Found duplicate file name in mapping file [$c[1]]") if( exists( $second_column{$c[1]} ) );

	$retval{$c[0]} = $c[1];
	$second_column{$c[1]} = 1;
	
  }
  close(IN);

  return %retval;
}

sub print_header {
  my ($fh) = @_;
  print $fh join("\t", qw($;I_FILE_BASE_1$; $;I_FILE_NAME_1$; $;I_FILE_PATH_1$; $;I_FILE_EXT_1$; $;I_DIR_1$; $;I_FILE_BASE_2$; $;I_FILE_NAME_2$; $;I_FILE_PATH_2$; $;I_FILE_EXT_2$; $;I_DIR_2$;) )."\n";
}

sub print_output_pair {
  my ($fh, $pair, $seen_pairs) = @_;
  my @file_parts1 = &parse_file_parts( $pair->[0] );
  my @file_parts2 = &parse_file_parts( $pair->[1] );
  my $pair_key = $file_parts1[0].$file_parts2[0];
  die("Duplicate basename combination in input set [$pair_key]") if( exists( $seen_pairs->{$pair_key} ) );
  $seen_pairs->{$pair_key} = 1;
  print $fh join("\t", @file_parts1, @file_parts2 )."\n";
}

sub _get_all_combinations {
  #Umm, yeah. Something like this: http://www.perlmonks.org/?node_id=899125
  #It will only work on strings
  map { [split(/===/)] } glob("{".join(",",@{$_[0]})."}==={".join(",",@{$_[1]})."}");
}

sub _get_inputs {
  my ($list, $file) = @_;
  my @files;
  if( $list ) {
	open(IN, "< $list") or die("Can't open $list: $!");
	while( <IN> ) {
	  chomp;
	  push(@files, $_);
	}
	close(IN);
  } 

  push(@files, $file) if( $file );

  return @files;
}

sub check_parameters {
  
  my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
  $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
								'LOG_LEVEL'=>$options{'debug'});
  $logger = Ergatis::Logger::get_logger();

  # display documentation
  if ( $options{'help'} ) {
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
  }

  ## handle required arguments
  my @required = qw( output_iter_list );
  for ( @required ) {
	unless ( $options{$_} ) {
	  print STDERR "--$_ is a required option\n\n";
	  exit(1);
	}
  }

  $options{'all_vs_all'} = 0 unless( exists( $options{'all_vs_all'} ) );

  unless( $options{'input_pairs_list'} || 
		  ( ($options{'input_file_list_1'} || $options{'input_file_1'} )  && 
			($options{'input_file_list_2'} || $options{'input_file_2'} ) ) ) {
	print STDERR "No input specified. Need something for input 1 and input 2 or input_pairs_list\n";
	exit(1);
  }

}
