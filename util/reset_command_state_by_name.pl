#!/usr/bin/perl

=head1 NAME

reset_command_state_by_name.pl.pl - Modifies workflow XML and sets the status of all commands of a
a given name to 'incomplete'

=head1 SYNOPSIS

USAGE: create_evidence_file_mapping.pl

    --pipeline_xml=/path/to/pipeline.xml
    --name=blast2btab
  [ --status=incomplete
    --log=/path/to/file.log
    --help
  ]

=head1 OPTIONS

B<--pipeline_xml, -p>

    REQUIRED. Any Workflow XML file.  This can be a top-level pipeline.xml, a group file, etc.

B<--name, -n>

    REQUIRED. The name of the command type to reset.

B<--status, -s>

    The status that selected commands should be set to.  Default is 'incomplete'


B<--log,-l>

    Logfile.

B<--help,-h>

    Print this message

=head1  DESCRIPTION

This script parses through any Workflow XML file (and included files) search for command elements
of a given name and sets their status element to a user-defined value (usually incomplete.)  This
can be used to reset just one step of a larger, overall pipeline if it has a unique name.
 
=head1  INPUT

Any pipeline XML that contains elements like this directly or through includes:

    <command>
        <name>hmmpfam</name>
        <id>1440412084</id>
        <startTime>2009-09-24T03:33:31.672-04:00</startTime>
        <endTime>2009-09-24T03:33:38.676-04:00</endTime>
        <state>complete</state>
        <type>RunUnixCommand</type>
        <executable>/usr/local/bin/hmmpfam</executable>
        <retryCount>0</retryCount>
        <retryAttempts>0</retryAttempts>
        <timeOut>0</timeOut>
        <status>
          <retValue>0</retValue>
        </status>
        </param>
        <arg>--acc  /usr/local/projects/db/TIGRFAM_equivalogs/TIGRFAMequivalogs.HMM.LIB.bin /usr/local
        /projects/jorvis/output_repository/split_multifasta/3662_default/i1/g1/RSP_6143.fsa</arg>
    </command>

Of these, only the following are expected:

    <command>
        <name>hmmpfam</name>
        <state>complete</state>
        <type>RunUnixCommand</type>
    </command>

It makes no difference if the input files are compressed via gzip or not, as usual.

=head1 OUTPUT

All files are read, a '.new' version written, then the '.new' version is moved over the original.  The
only changes are the the <state> element values of matching command elements as well as some additional 
whitespace created by the XML::Twig module.

=head1  CONTACT

    Joshua Orvis
    jorvis@users.sf.net

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use XML::Twig;
use File::Copy;

my %options;
my $results = GetOptions( \%options,
                          'pipeline_xml|p=s',
                          'name|n=s',
                          'status|s=s',
                          'log|l=s',
                          'debug|d=s',
                          'help|h=s' );

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

&check_options( \%options );

## recursive call to parse all included files in the pipeline
process_file( $options{pipeline_xml} );


exit(0);

############################################################


sub process_file {
    my $path = shift;
    
    print STDERR "processing $path\n";
    
    my $fh_in = get_conditional_read_fh($path);
    my $fh_out = get_conditional_write_fh($path);

    if ( $fh_in ) {
        my $t = XML::Twig->new(
                           twig_roots => {
                               'command' => sub {
                                   my ($t, $elt) = @_;
                                   
                                   if ( $elt->first_child('name')->text() eq $options{name} ) {
                                        $elt->first_child('state')->set_text( $options{status} );
                                   }
                                   
                                   $elt->print($fh_out);
                               },
                               'commandSet/fileName' => sub {
                                   my ($t, $elt) = @_;
                                   process_file( $elt->text );
                                   
                                   $elt->print($fh_out);
                               },
                           },
                           twig_print_outside_roots => $fh_out,
                           pretty_print => 'indented',
                        );

        $t->parse( $fh_in ); 
        #$t->flush( $fh_out );
        
        ## now move the temp file over the parsed one
        close $fh_in;
        close $fh_out;
        
        move("$path.new", "$path");
        
        
    } else {
        print "WARNING: skipping absent $path\n";
    }
}


sub get_conditional_write_fh {
    my $source_file = shift;
    
    my $is_zipped = 0;
    
    if ( $source_file =~ /\.gz$/ ) {
        $is_zipped = 1;
    }
    
    my $ofh;
    
    if ( $is_zipped ) {
        open($ofh, ">:gzip", "$source_file.new") || die "failed to create output file: $!";
    } else {
        open($ofh, ">$source_file.new") || die "failed to create output file: $!";
    }
    
    return $ofh;
}

sub get_conditional_read_fh {
    my $path = shift;
    my $fh;
    my $found = 0;
    
    if ( -e $path ) {
        $found = 1;
    
    ## auto-handle gzipped files
    } elsif (! -e $path && -e "$path.gz") {
        $path .= '.gz';
        $found = 1;
    } 
    
    if (! $found ) {
        ## we can't find the file.  just return an empty file handle. 
        ##  the process sub is OK with this.
        return $fh;
    }

    if ( $path =~ /\.gz$/ ) {
        open($fh, "<:gzip", $path) || die "can't read file $path: $!";
    } else {
        open($fh, "<$path") || die "can't read file $path: $!";
    }
    
    return $fh;
}

## checks the options to the program
sub check_options {
    my $options = shift;

    ## make sure required arguments were passed
    my @required = qw( pipeline_xml name );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
    
    ##
    ## you can do other things here, such as checking that files exist, etc.
    ##
    
    ## handle some defaults
    $options{status} = 'incomplete'  unless ( defined $options{status});
    
}




