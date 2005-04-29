#!/usr/local/bin/perl

=head1  NAME 

bsml2fasta.pl - convert BSML files to fasta

=head1 SYNOPSIS

USAGE:  bsml2fasta.pl --bsml_input=/path/to/fileORdir | --bsml_list=/path/to/file
        [--class_filter=class --debug=debug_level --log=log_file]

=head1 OPTIONS

B<--bsml_input,-i> 
    Input files or folders.  Can be a comma-separated list of mixed input types.

B<--bsml_list,-s> 
    Text file that is a list of input files and/or folders.

B<--class_filter,-c> 
    Filter output sequences to only include those in the passed class.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--format,-f> 
    Format.  'multi' (default) writes all sequences to a multi-fasta file, and 'single' writes each sequence in a separate file named like $id.fsa

B<--log,-l> 
    Log file

B<--output,-o> 
    Output file (if --format=multi) or directory (if --format=single)

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used to convert BSML to fasta.  The input is meant to be as flexible 
as possible and is described below.  The output can be either a single file with 
multiple entries, or separate files for each entry.

=head1 INPUT

Input can be either a single file, a directory of files, or a list.  A list is
a text file containing a list of either files or directories, one on each line.

Any individual .bsml file can have as many <Sequence> elements as it likes.  They
will all get exported as FASTA unless you use the --class_filter option to include
only certain classes.

The --bsml_input option is used to pass files or directories.  The following
example passes a single bsml file:

    bsml2fasta.pl --bsml_input=/foo/bar/dnaA.bsml

We can pass multiple files in a comma-separated list (use quotes):

    bsml2fasta.pl --bsml_input="/foo/bar/dnaA.bsml, /foo/bar/tonB.bsml"

If /foo/bar had many .bsml files in it and we want to process them all, we could 
just do:

    bsml2fasta.pl --bsml_input=/foo/bar

Also multiple directories are ok:

    bsml2fasta.pl --bsml_input="/foo/bar, /home/you"

You can even mix directories and files:

    bsml2fasta.pl --bsml_input="/foo/bar/dnaA.bsml, /home/you"

Lists are useful if you have a specific set of files or directories to process.
If we have a file named, for example, '/foo/bar/neatstuff.list', which has contents 
like (the .list suffix is completely optional):

    /foo/bar/thingy1.bsml
    /foo/bar/thingy2.bsml
    ...
    /foo/bar/thingy1034.bsml

You can pass this list to the script:

    bsml2fasta.pl --bsml_list=/foo/bar/neatstuff.list

If one of the lines in the list is a path to a directory rather than a file, each
.bsml file in that directory will be processed.  A list can contain paths to both
files and directories such as:

    /foo/bar/thingy1.bsml
    /foo/bar/thingy2.bsml
    ...
    /foo/bar/thingy1034.bsml
    /home/you

Finally, you can be really messy and mix input types and methods to process files,
directories and lists all at the same time:

    bsml2fasta.pl --bsml_input="/foo/bar/dnaA.bsml, /home/you" --bsml_list=/home/you/some.list

Once everything is evaluated down to the file-level, all will be skipped that
don't have a .bsml suffix.

=head1 OUTPUT

Output is specified using the required --output and optional --format options.  By
default the output will be a single file containing multiple FASTA entries.  So:

    bsml2fasta.pl --bsml_input=/foo/bar --output=/home/you/seqs.fsa

This would read all the .bsml files in /foo/bar and write their sequences to the
seqs.fsa file in /home/you in multi-FASTA format.  If you want each sequence to be
output separately, you need to use the --format=single option:

    bsml2fasta.pl --bsml_input=/foo/bar --output=/home/you/data --format=single

This would write each sequence to its own .fsa file into the /home/you/data directory.
Each file will be named using the id attribute of each sequence, like $id.fsa .
--format=multi is the default and does not need to be passed explicitly.  Note that
the only legal characters for the file name are in the set [ a-z A-Z 0-9 - _ . ].  Any
other characters will be replaced with underscores.  These replacements will only occur
in the file name;  the header will still have the original id.  The complete fasta header
for each entry will be composed of the id followed by a single space, and then the
'title' attribute of the sequence, if one exists.  For example:

    >gi28804993 chromosomal DNA replication initiator DnaA [Vibrio parahaemolyticus]

If you are reading in multiple sequences and want to filter which are included in
the output, you can do this based on their class attributes by passing the
--class_filter option:

    bsml2fasta.pl --bsml_input=/foo/bar --output=/home/you/seqs.fsa --class_filter=cds

This will filter the sequences in /foo/bar to only include those that are of the
cds class.  Sequences without class attributes will not be included if --class_filter
is used.  The values of these class attributes can be anything, but should correspond 
to SO terms.

=head1 CONTACT

Joshua Orvis
jorvis@tigr.org

=cut

use strict;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Workflow::Logger;

use BSML::BsmlReader;
use BSML::BsmlParserTwig;

#######
## ubiquitous options parsing and logger creation
my %options = ();
my $results = GetOptions (\%options, 
                            'bsml_input|i=s',
                            'bsml_list|s=s',
                            'bsml_file|b=s',        ## deprecated
                            'bsml_dir|d=s',         ## deprecated
                            'parse_element|p=s',
                            'class_filter|c=s',
                            'type|t=s',             ## deprecated
                            'output_list|u=s',
                            'debug=s',
                            'format|m=s',
                            'log|l=s',
                            'output|o=s',
                            'help|h') || pod2usage();

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = Workflow::Logger::get_logger();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

&check_parameters(\%options);

#######
## gather the list of files we're going to processes.  the user can pass file or
#  directory names using --bsml_input, which can be a single file, a list of 
#  filenames, a directory of files, or a list of directories of files (lists are
#  comma-separated.  it can even be a mixed list of file and dir names.  fancy!
#  a file containing a list of file/dir names should be passed with --bsml_list.  
my @files;
my @list_elements = ();

## did the user pass a list?  if so, get each file/dir names out of it, check it, and add it
if ( $options{bsml_list} ) {
    for my $list ( split(/,/, $options{bsml_list}) ) {
        $list =~ s/\s//g;
        
        if (-f $list && -s $list) {
            open (my $fh, $list) || $logger->logdie("Can't open file $list");
            for ( <$fh> ) {
                chomp;
                if (-e $_ && -s $_) {
                    push @list_elements, $_;
                } else {
                    $logger->warn("Error reading $_ from list $list") if ($logger->is_warn);
                }
            }
        } else {
            $logger->warn("Error reading $list") if ($logger->is_warn);
        }
    }
}

## loop through each input thing passed
for my $thing ( split(/,/, ($options{bsml_input} || '') ),
                split(/,/, ($options{bsml_file}  || '') ),  ## backwards compatibility
                split(/,/, ($options{bsml_dir}   || '') ),  ## backwards compatibility
                @list_elements) {
    $thing =~ s/\s//g;
    next unless ($thing);
    
    ## is this a directory?
    if (-d $thing) {
        opendir(my $dh, $thing) || $logger->warn("Unable to access $thing due to $!");

        for my $file (readdir $dh) {
            &add_file("$thing/$file");
        }

    ## else it's probably a file, make sure it exists and has size and add it, else warn
    } elsif (! &add_file($thing) ) {
        $logger->warn("Error reading $thing") if ($logger->is_warn);
    }
}

## die if no files found
if ( scalar @files == 0 ) {
    $logger->logdie("No files found");
}


## create an output list if passed
my $olist_fh;
if ($options{output_list}) {
    $logger->debug("Creating output list at $options{output_list}") if ($logger->debug);
    open($olist_fh, ">$options{output_list}") || $logger->logdie("couldn't create output list $options{output_list} : $!");
}


## remember which IDs have been found already
my %ids;

#######
## parse out sequences from each file
my $mfh;
my $title;
for my $file ( @files ) {
	my $parser = new BSML::BsmlParserTwig();
	my $reader = new BSML::BsmlReader();
	$logger->debug("Parsing entire file $file") if ($logger->debug);
	$parser->parse( \$reader, $file );

    my $sfh;
    
    if ($options{format} eq "multi") {
        open ($mfh, ">>$options{output}") || $logger->logdie("Unable to write to file $options{output} due to $!");
        $logger->debug("Writing output to multi-fasta file $options{output}") if ($logger->is_debug);
        
        if ($options{output_list}) {
            print $olist_fh "$options{output}\n";
        }
    }
    
    my $seq_list = $reader->returnAllSequences();
    
    ## iterate through BSML::BsmlSequence objects
    my $seq;
    for my $seq_obj (@$seq_list) {
        my $seq_id = $seq_obj->returnattr('id') || $seq_obj->returnattr('identifier');
        $title = $seq_obj->returnattr('title') || '';
        
        ## make sure an identifier was found
        if (! $seq_id) {
            $logger->error("Cowardly refusing to create a fasta entry with no header.  add an id or identifier.") if ($logger->is_error);
            next;
        }
        
        ## if the user passed a class_filter, skip this sequence if it isn't the right class.
        if ($options{parse_element} eq 'sequence') {
            if ($options{class_filter} && $options{class_filter} ne $seq_obj->returnattr('class')) {
                $logger->debug("Skipping $seq_id in $file because it does not match the class_filter passed") if ($logger->is_debug);
                next;
            }
        }
        
        ## check and make sure we have a sequence, else log error and skip it
        ##  (whole sequence is returned with -1 passed as start)
        $seq = $reader->subSequence($seq_obj,-1,0,0);
        $seq =~ s/\s//g;
        
        if (length $seq < 1) {
            $logger->error("sequence $seq_id has no length.  cowardly refusing to export anything.") if ($logger->is_error);
            next;
        }
        
        if ($options{parse_element} eq 'sequence') {
            ## multi or single format?
            write_sequence($seq_id, \$seq);

        } elsif ($options{parse_element} eq 'feature') {
            ## if we are parsing features we need to get all the feature elements within this Sequence.
            my $feat_list = $reader->readFeatures($seq_obj);
            
            ## read each feature.
            foreach my $feat (@$feat_list) {
                my $feat_id = $feat->{id} || die "didn't get an id";
            
                ## are we checking for a specific class?
                if ($options{class_filter} && $options{class_filter} ne $feat->{class}) {
                    $logger->debug("Skipping $feat_id of $seq_id in $file because it does not match the class_filter passed") if ($logger->is_debug);
                    next;
                }

                my $positions = $feat->{locations}->[0];
                
                $logger->debug("attempting to extract $positions->{startpos}, $positions->{endpos}, $positions->{complement}") if ($logger->is_debug);
                my $feat_seq = $reader->extractSequence( \$seq, $positions->{startpos}, $positions->{endpos}, $positions->{complement} );
            
                ## now write the sequence
                write_sequence($feat_id, \$feat_seq);
            }
        }
    }
}

#######
## fin   
exit;

sub add_file {
    ## adds a file to the list of those to process, checking to make sure it
    #  exists and is populated.
    my $file = shift;
    
    ## only do .bsml files
    return 0 unless ( $file =~ /\.bsml$/ );
    
    if (-e $file && -s $file) {
        $logger->debug("Adding file $file for processing") if ($logger->is_debug);
        push @files, $file;
    } else {
        $logger->warn("Error reading file $file") if ($logger->is_warn);
        return 0;
    }
    
    return 1;
}

sub check_parameters {
    my ($options) = @_;
    
    ## they have to pass some form of input
    unless ($options{bsml_input} || $options{bsml_list} ||
            $options{bsml_file}  || $options{bsml_dir}) {
        $logger->logdie("You must specify input with --bsml_input --bsml_list --bsml_file or --bsml_dir");
    }

    ## output is required
    unless ( $options{output} ) {
        $logger->logdie("You must specify an output directory or file with --output");
    }

    ## check the format setting or set a default if it wasn't passed
    if (! $options{format}) {
        $options{format} = 'multi';
    } elsif ( $options{format} ne 'single' && $options{format} ne 'multi' ) {
        $logger->logdie("--format must be either 'single' or 'multi'");
    }
    
    ## if format is single, output must be a directory
    if ($options{format} eq 'single') {
        unless (-d $options{output}) {
            $logger->logdie("if using --format=single then --output must point to a directory");
        }
    ## else if format is multi, output must NOT be a directory
    } elsif ($options{format} eq 'multi') {
        if (-d $options{output}) {
            $logger->logdie("if using --format=multi then --output must NOT point to a directory");
        }
    }

    ## if the user passed a deprecated --type, change it to --class_filter
    if ($options{type}) {
        $options{class_filter} = $options{type};
    }

    ## handle some defaults
    if (defined $options{parse_element}) {
        if ($options{parse_element} ne 'sequence' && $options{parse_element} ne 'feature') {
            $logger->logdie("parse_element must be either 'sequence' or 'feature'");
        }
    } else {
        $options{parse_element} = 'sequence';
    }

    if(0){
        pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});    
    }
}

sub fasta_out {
    ## This subroutine takes a sequence name and its sequence and
    ## outputs a correctly formatted single fasta entry.
    my ($header, $seqref, $fh) = @_;

    print $fh ">$header\n";
    for(my $i=0; $i < length($$seqref); $i+=60){
        print $fh substr($$seqref, $i, 60), "\n";
    }
}

sub write_sequence {
    my ($id, $seqref) = @_;

    ## has this ID already been found?
    if ($ids{$id}) {
        $logger->logdie("duplicate IDs found ($id).  quitting.");
    }
    
    $ids{$id}++;

    if ($options{format} eq 'multi') {
        &fasta_out("$id $title", $seqref, $mfh);

    } elsif ($options{format} eq 'single') {
        ## legal characters only in the output file name:
        my $name = $id;
        $name =~ s/[^a-z0-9\.\-]/_/gi;

        $logger->debug("attempting to create file $options{output}/$name.fsa") if ($logger->is_debug);
        open (my $sfh, ">$options{output}/$name.fsa") || $logger->logdie("Unable to write to file $options{output}/$name due to $!");
        &fasta_out("$id $title", $seqref, $sfh);

        if ($options{output_list}) {
            print $olist_fh "$options{output}/$name.fsa\n";
        }
    }

}

