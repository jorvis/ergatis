#!/usr/bin/perl

=head1 NAME

create_evidence_file_mapping.pl - Description

=head1 SYNOPSIS

USAGE: create_evidence_file_mapping.pl

    --input_file=/path/to/gene.describing.bsml.list
    --analysis_list=list1,list2,list3...
    --list_of_analysis_lists=/path/to/list.of.lists.file
    --output=/path/to/mapping/file.tab
    --qualified_list_out_dir=/path/to/dir
    --exclude_empty=1
    --empty_file_list=/path/to/list
  [ --log=/path/to/file.log
    --help
  ]

=head1 OPTIONS

B<--input_file,-i>

    REQUIRED. A gene describing BSML file or list file of bsml files (.all.bsml)

B<--analysis_list,-a>

    REQUIRED or list_of_analysis_lists option. List file of output from analysis. Can be 
    comma separated list of list files.

B<--list_of_analysis_lists,-s>

    REQUIRED or analysis_lists option.  List of analysis lists (described above).

B<--output,-o>

    REQUIRED. Output mapping file.

B<--qualified_list_out_dir,-q>

    OPTIONAL. If given a valid directory, will print a list of evidence files found in the
    analysis lists (one file per analysis list) relating to genes described in the input_file.

B<--exclude_empty,-e>

    OPTIONAL.  A non-zero value will exclude the zero size files from qualified lists.

B<--empty_file_list,-m>

    OPTIONAL.  Where analysis files relating to genes described in input_list with a size of
    zero will be printed.

B<--log,-l>

    Logfile.

B<--help,-h>

    Print this message

=head1  DESCRIPTION

    create_evidence_file_mapping.pl will take gene describing bsml files and find evidence
    based on those polypeptide ids. It can print these to a tab delimited file (format described below).
 
=head1  INPUT
    
    The input bsml file must have Feature[@class="polypeptide"] elements with an id.  Currently, 
    the script will use these ids to search file names from the analysis lists to create the mapping.

=head1 OUTPUT

There a couple outputs to the script. Each will be described briefly.

=over 4

=item    * mapping list (--output)

    polypeptide.id     /path/to/analysis_file.1     /path/to/analysis.2
    polypeptide.id ...

=item    * qualified lists (--qualified_list_out_dir)
    
    If, for example, 3 analysis list files were input (say list.a, list.b, list.c) there would
    be three output lists in the qualified_list_out_dir (qualified.list.a, qualified.list.b,
    qualified.list.c)

=item    * list of empty file names (--empty_file_list)

    If the empty_file_list option is specified, this file will create a list of evidence
    files that were size zero.
    
=back

=head1  CONTACT

    Kevin Galens
    kgalens@tigr.org

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use XML::Twig;
use Data::Dumper;
use File::Basename;

my %options;
my $results = GetOptions( \%options,
                          'input_file|i=s',
                          'analysis_list|a=s',
                          'list_of_analysis_lists|s=s',
                          'qualified_list_out_dir|q=s',
                          'output|o=s',
                          'exclude_empty|e=s',
                          'empty_file_list|m=s',
                          'log|l=s',
                          'debug|d=s',
                          'help|h=s' );

################# GLOBALS ##################################
my @input_files;
my @analysis_lists;
my $out_map;
my %polypeptide_ev_map = ();
my $qualified_list_out_dir = undef;
my $exclude_empty = 0;
my $empty_file_list = undef;
my @empty_files;
my $logfh;
############################################################

&check_options( \%options );

my @polypeptide_ids = &get_polypeptides( \@input_files );

%polypeptide_ev_map = &find_evidence_files( \@polypeptide_ids, \@analysis_lists );

&print_mapping_and_list_files( \%polypeptide_ev_map, $out_map );

if( $empty_file_list ) {
    &print_empty_file_list( \@empty_files );
}

print "$0 execution completed\n";
exit(0);

############################################################

## finds evidence files for a given list of polypeptide ids
sub find_evidence_files {
    my ($polypeptide_ids, $analysis_files) = @_;
    my %polypeptide_ev_hash = ();
    
    my $found_ctr=0;
    my $analysis_file_ctr=0;

    foreach my $analysis_file( @{$analysis_files} ) {

	$analysis_file_ctr++;

        my $afh = open_file( $analysis_file );
        chomp( my @list = <$afh> );
        
        foreach my $polypeptide_id( @{$polypeptide_ids} ) {
            
            my @found_files = grep { /$polypeptide_id\./ } @list;

	    if (scalar(@found_files) > 0 ){
		push(@{$polypeptide_ev_hash{$polypeptide_id}->{$analysis_file}}, @found_files);
		$found_ctr += scalar(@found_files) + 1;
	    } 
        }
    }


    if ($analysis_file_ctr == 0){
	die "Did not process any analysis list files!";
    } else {
	print "Processed '$analysis_file_ctr' analysis list files\n";
    }

    if ($found_ctr == 0){
	die "Did not find any matching polypeptide files!";
    } else {
	print "Found '$found_ctr' polypeptide files\n";
    }
    return %polypeptide_ev_hash;
}

## Parse the polypeptide ids from the input bsml files
sub get_polypeptides {
    my $input_files = @_;
    my @polypeptide_ids = ();

    my $twig = new XML::Twig( 'twig_handlers' => {
        'Feature[@class="polypeptide"]' => sub { push(@polypeptide_ids, $_[1]->att('id')) }
        } );

    foreach my $input_file( @input_files ) {
        my $fh = open_file( $input_file );
        $twig->parse( $fh );
    }

    print $logfh "Parsed ".scalar(@polypeptide_ids)." polypeptide ids from BSML\n";

    if( @polypeptide_ids == 0 ){
        print $logfh "ERROR: Could not find any polypeptides in BSML\n";
        die("Could not find any polypeptide ids in the BSML");
    }

    return @polypeptide_ids;

}

## makes the list of empty files
sub print_empty_file_list {
    my ($empty_files) = @_;

    open( my $eh, "> $empty_file_list") or die("Cannot open $empty_file_list for writing ($!)");

    local $" = "\n";
    print $eh "@{$empty_files}\n";
    
    close( $eh );

    print $logfh "Printed empty list file $empty_file_list\n";
    
}

## makes the output mapping file
sub print_mapping_and_list_files {
    my ($polypeptide_ev_map, $out_file) = @_;

    ## open the necesary mapping file and output files
    open( my $out, "> $out_file" ) or die("Unable to open $out_file for writing ($!)");

    my %qualified_lists = ();
    
  POLYPEPTIDE_ID:
    foreach my $polypeptide_id( keys %{$polypeptide_ev_map} ) {
        print $out "$polypeptide_id\t";
        
      ANALYSIS_LIST: 
        foreach my $analysis_list( keys %{$polypeptide_ev_map->{$polypeptide_id}} ) {

          ANALYSIS_FILE:
            foreach my $analysis_file( @{$polypeptide_ev_map->{$polypeptide_id}->{$analysis_list}} ) {
                print $out "$analysis_file\t";

                if( -z $analysis_file ) {
                    push( @empty_files, $analysis_file );
                    
                    if( $exclude_empty ) {
                        next ANALYSIS_FILE;
                    }

                }

                push( @{$qualified_lists{$analysis_list}}, $analysis_file );
                
            }
        }

        print $out "\n";
    }
    
    my %alist_to_qlist_lookup;
    foreach my $analysis_list ( keys %qualified_lists ) {
        
        
        my $base_name = $1 if( $analysis_list =~ m|/([^/]+)$| );
        unless( $base_name ) {
            print $logfh "ERROR: Could not parse base_name from $analysis_list\n";
            die("Could not parse base_name from $analysis_list") unless( $base_name );
        }

	if (!defined($qualified_list_out_dir)){
	    $qualified_list_out_dir = '/tmp';
	}
        my $list_name = "$qualified_list_out_dir/qualified.$base_name";

        $alist_to_qlist_lookup{$analysis_list} = $list_name;
        
    }

    #check list names to make sure they don't exist.  If they do, back them up
    foreach my $list ( values %alist_to_qlist_lookup ) {
        if ( -e $list ) {
            eval {
                print $logfh "Making back up of $list\n";
                system("mv $list $list.bak");
            };
            if( $@ ) {
                print $logfh "Could not back up $list to $list.bak. Please remove list manually\n";
                die("Could not move $list to $list.bak ($@)");
            }
        }
    }

    #Actually print the qualified list files
    foreach my $analysis_list ( keys %qualified_lists ) {
        my $list_name = $alist_to_qlist_lookup{ $analysis_list };

        open( my $fh, ">> $list_name" ) or die("Can't open $list_name for writing ($!)");
        print $logfh "Printing qualified evidence files from $analysis_list to $list_name\n";

        local $" = "\n";
        print $fh "@{$qualified_lists{$analysis_list}}\n";

        close( $fh );

    }


    
}

## checks the options to the program
sub check_options {
    my $opts = shift;

    if( $opts->{'input_file'} ) {
        
        my $in = open_file( $opts->{'input_file'});
        chomp( my $line = <$in> );

        #Determine if it's a bsml file.  If not, see if it's a list
        if( $line =~ /^</ ) {
            push( @input_files, $opts->{'input_file'} );
        } elsif(  &file_exists($line) ) {
            push( @input_files, $line );
            chomp( my @tmp = <$in> );
            push( @input_files, @tmp );
        } else {
            die("Invalid input file");
        }

    } else {
        die("Option input_file is required");
    }

    if (! (( $opts->{'analysis_list'} ) ||( $opts->{'list_of_analysis_lists'} )) ) {
	die "You must specify --analysis_list or --list_of_analysis_lists";
    }

    if( $opts->{'analysis_list'} ) {
        @analysis_lists = split( /,\s*/, $opts->{'analysis_list'} );
    }

    if( $opts->{'list_of_analysis_lists'} ) {
        my $fh = open_file( $opts->{'list_of_analysis_lists'} );
        chomp( my @lists = <$fh> );
        push( @analysis_lists, @lists );
    }
    
    if( @analysis_lists == 0 ) {
        die("No analysis lists were found");
    }

    if( $opts->{'output'} ) {
        $out_map = $opts->{'output'};
    } else {
	$out_map = '/tmp/' . File::Basename::basename($0) . '.mapping-file.txt';
	print STDERR "output was not specified and therefore was set to '$out_map'\n";
    }

    if( $opts->{'exclude_empty'} ) {
        $exclude_empty = 1;
    }

    if( $opts->{'qualified_list_out_dir'} ) {
        $qualified_list_out_dir = $opts->{'qualified_list_out_dir'};
        
        if( !(-d $qualified_list_out_dir) ) {
            die("Value for option qualified_list_out_dir ($qualified_list_out_dir) does".
                " not exist or is not a directory");
        }
    }

    if( $opts->{'log'} ) {
        open( $logfh, "< $opts->{'log'}") or die("Unable to open $opts->{'log'} ($!)");
    } else {
        $logfh = *STDOUT;
    }

    if( $opts->{'exclude_empty'} ) {
        $exclude_empty = 1;
    } else {
        $exclude_empty = 0;
    }

    if( $opts->{'empty_file_list'} ) {
        $empty_file_list = $opts->{'empty_file_list'};

        if( -e $empty_file_list ) {
            print $logfh "file for empty_file_list exists. Overwriting.\n";
        }
    }

    
    
}

## open file wrapper for opening gzip files.
sub open_file {
    my ($file) = @_;
    my $fh;

    if( $file =~ /\.gz$/ ) {
        open( $fh, "<:gzip", $file ) or die("Could not open $file ($!)");
    } elsif( -e $file ) {
        open( $fh, "<$file" ) or die("Could not open $file ($!)");
    } elsif( -e $file.".gz" ) {
        open($fh, "<:gzip", $file ) or die("$file does not exist, but $file.gz does.".
                                           "Failed opening it. ($!)");
    } else {
        die("$file does not exist");
    }

    return $fh;

}

## checks to see if the file or a gzip version of the file exists.
sub file_exists {
    my ($file) = @_;
    my $retval = 0;

    if( -e $file ) {
        $retval = 1;
    } elsif( -e $file.".gz" ) {
        $retval = 1;
    }

    return $retval;
    
}




