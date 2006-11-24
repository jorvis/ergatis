#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

fasta2bsml.pl - convert fasta files to BSML

=head1 SYNOPSIS

USAGE:  fasta2bsml.pl 
        --fasta_input=/path/to/fileORdir | --fasta_list=/path/to/file
        --format=multi|single
        --output=/path/to/somefile.fsa   | /path/to/somedir
      [ --output_list=/path/to/somefile.list
        --output_subdir_size=20000
        --output_subdir_prefix='somename'
        --debug=debug_level 
        --log=log_file 
      ]

=head1 OPTIONS

B<--fasta_input,-i> 
    Input files or folders.  Can be a comma-separated list of mixed input types.

B<--fasta_list,-s> 
    Text file that is a list of input files and/or folders.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--format,-f> 
    Format.  'multi' (default) writes all sequences to a multi-entry bsml file, and 'single' writes each sequence in a separate file named like $id.bsml

B<--log,-l> 
    Log file

B<--output,-o> 
    Output file (if --format=multi) or directory (if --format=single)

B<--output_list,-u>
    Optional.  If passed, will create a list file containing the path to each of the
    files created.
    
B<--output_subdir_size,-z>
    Optional.  Number of files to create within each output subdirectory.
    
B<--output_subdir_prefix,-x>
    Optional.  Rather than just plain numberical names (N), output subdirectories will
    be named like "prefixN" if you pass a prefix with this.

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used to convert fasta to BSML.  The input is meant to be as flexible 
as possible and is described below.  The output can be either a single file with 
multiple <Sequence> entries, or separate files for each entry.

=head1 INPUT

Input can be either a single file, a directory of files, or a list.  A list is
a text file containing a list of either files or directories, one on each line.

Any individual fasta file can have a single entry or multiple entries.  This
script will separate them as needed.

The --fasta_input option is used to pass files or directories.  The following
example passes a single fasta file:

    fasta2bsml.pl --fasta_input=/foo/bar/dnaA.fna

We can pass multiple files in a comma-separated list (use quotes):

    fasta2bsml.pl --fasta_input="/foo/bar/dnaA.fna, /foo/bar/tonB.faa"

If /foo/bar had many fasta files in it and we want to process them all, we could 
just do:

    fasta2bsml.pl --fasta_input=/foo/bar

Also multiple directories are ok:

    fasta2bsml.pl --fasta_input="/foo/bar, /home/you"

You can even mix directories and files:

    fasta2bsml.pl --fasta_input="/foo/bar/dnaA.fna, /home/you"

Lists are useful if you have a specific set of files or directories to process.
If we have a file named, for example, '/foo/bar/neatstuff.list', which has contents 
like (the .list suffix is completely optional):

    /foo/bar/thingy1.fna
    /foo/bar/thingy2.faa
    ...
    /foo/bar/thingy1034.fsa

You can pass this list to the script:

    fasta2bsml.pl --fasta_list=/foo/bar/neatstuff.list

If one of the lines in the list is a path to a directory rather than a file, each
fasta file in that directory will be processed.  A list can contain paths to both
files and directories such as:

    /foo/bar/thingy1.bsml
    /foo/bar/thingy2.bsml
    ...
    /foo/bar/thingy1034.bsml
    /home/you

Finally, you can be really messy and mix input types and methods to process files,
directories and lists all at the same time:

    fasta2bsml.pl --fasta_input="/foo/bar/dnaA.bsml, /home/you" --fasta_list=/home/you/some.list

Once everything is evaluated down to the file-level, all will be skipped that
don't have one of the following suffices: .fsa .faa .fna .fasta

=head1 OUTPUT

The fasta format is limited in its data representation in the header.  Because of this,
this script was written very generically to allow any sort of header.  In each <Sequence>
element created, the attributes generated are length, title, and id.  (other attributes
my be handled by default in BsmlBuilder).  The id is assumed to be the first part of
the fasta header up to the initial whitespace.  So with this header:

    >gi46446716  putative chromosomal replication initiator protein, dnaA

The 'id' attribute will be set to 'gi46446716'.  The 'title' attribute is set to
the entire value of the fasta header, including the first word that became the id.

Output is specified using the required --output and optional --format options.  By
default the output will be a single file containing multiple sequences entries.  So:

    fasta2bsml.pl --fasta_input=/foo/bar --output=/home/you/seqs.bsml

This would read all the fasta files in /foo/bar and write their sequences to the
seqs.bsml file in /home/you in multi-entry format.  If you want each sequence to be
output separately, you need to use the --format=single option:

    fasta2bsml.pl --fasta_input=/foo/bar --output=/home/you/data --format=single

This would write each sequence to its own .bsml file into the /home/you/data directory.
Each file will be named using the id attribute of each sequence, like $id.fsa .
--format=multi is the default and does not need to be passed explicitly.  Note that
the only legal characters for the file name are in the set [ a-z A-Z 0-9 - _ . ].  Any
other characters will be replaced with underscores.

You can pass a path to the optional --output_list to create a text file containing the full paths
to each of the BSML files created by this script.

Two other optional arguments, --output_subdir_size and --output_subdir_prefix, can be used
on input sets that are too large to write out to one directory.  This depends on the limitations
of your file system, but you usually don't want 100,000 files written in the same directory.

If you are going to create 95000 sequences, and use the following option:

    --output=/some/path
    --output_subdir_size=30000
    
The following will be created:

    directory              file count
    ---------------------------------
    /some/path/0/          30000
    /some/path/1/          30000
    /some/path/2/          30000
    /some/path/3/           5000

If you choose to create a list file (and you probably want to), it will contain these paths.

You may not want the subdirectories to simply be numbers, as above, so you can use the
--output_subdir_prefix option.  For example:        

    --output=/some/path
    --output_subdir_size=30000
    --output_subdir_prefix=bsml
    
The following will be created:

    directory              file count
    ---------------------------------
    /some/path/bsml0/     30000
    /some/path/bsml1/     30000
    /some/path/bsml2/     30000
    /some/path/bsml3/      5000

=head1 CONTACT

    Joshua Orvis
    jorvis@tigr.org

=cut


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Ergatis::Logger;
use BSML::BsmlBuilder;

my %options = ();
my $results = GetOptions (\%options,
                          'fasta_input|i=s',
                          'fasta_list|s=s',
                          'fasta_dir|d=s',      ## deprecated
                          'fasta_file|f=s',     ## deprecated
                          'output|o=s',
                          'output_list|u=s',
                          'output_subdir_size|z=s',
                          'output_subdir_prefix|x=s',
                          'format|m=s',
			              'log|l=s',
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

#######
## gather the list of files we're going to processes.  the user can pass file or
#  directory names using --fasta_input, which can be a single file, a list of 
#  filenames, a directory of files, or a list of directories of files (lists are
#  comma-separated.  it can even be a mixed list of file and dir names.  fancy!
#  a file containing a list of file/dir names should be passed with --fasta_list.  
my @files;
my @list_elements = ();

## did the user pass a list?  if so, get each file/dir names out of it, check it, and add it
if ( $options{fasta_list} ) {
    for my $list ( split(/,/, $options{fasta_list}) ) {
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
for my $thing ( split(/,/, ($options{fasta_input} || '') ),
                split(/,/, ($options{fasta_file}  || '') ),  ## backwards compatibility
                split(/,/, ($options{fasta_dir}   || '') ),  ## backwards compatibility
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

## variables for handling output sub directory groupings (if needed)
my $sub_dir_num = 0;
my $seq_file_count = 0;


#######
## parse out sequences from each file
my $doc;

## if we're writing out to a multi-entry file, create the doc now:
if ($options{format} eq 'multi') {
    $doc = new BSML::BsmlBuilder;
}

for my $file ( @files ) {
    
    my %seqs = loadMultiSequence( $file );
    
    for my $seqid ( sort {$a<=>$b} keys %seqs ) {
        
        ## capture the first element of the header up to the first whitespace
        my $id;
        if ($seqs{$seqid}{h} =~ /^(\S+)/) {
            $id = $1;
        } else {
            $logger->warn("unrecognized header format: $seqs{$seqid}{h}") if ($logger->is_warn);
        }
        
        ## are we writing each sequence to single files?  If so, start a new doc
        if ($options{format} eq 'single') {
            $doc = new BSML::BsmlBuilder;
        }
         
        my $seq_element = $doc->createAndAddSequence( $id, $seqs{$seqid}{h}, length($seqs{$seqid}{s}) );
        $logger->debug("adding id $id to the bsml doc") if ($logger->is_debug);
        $doc->createAndAddSeqData( $seq_element, $seqs{$seqid}{s} );
        
        ## record the defline
        $doc->createAndAddBsmlAttribute( $seq_element, 'defline', $seqs{$seqid}{h} );
        
        ## write this doc if we're in single mode.  the only thing we can use is the id,
        #  which needs to be made safe first.
        if ($options{format} eq 'single') {
            $id =~ s/[^a-z0-9\.\-]/_/gi;
            
            write_single_doc( $doc, $id );
            
            #$doc->write( "$options{output}/$id.bsml" );
            
            #if ($options{output_list}) {
            #    print $olist_fh "$options{output}/$id.bsml\n";
            #}
        }
    }
}

if ($options{format} eq 'multi') {
    $doc->write( $options{output} );
    
    if ($options{output_list}) {
        print $olist_fh "$options{output}\n";
    }
}

## fin
exit;


sub add_file {
    ## adds a file to the list of those to process, checking to make sure it
    #  exists and is populated.
    my $file = shift;
    
    ## only do .f?a files (.fna .fsa .faa .fasta) unless we're operating on a list
    if (! $options{fasta_list}) {
        return 0 unless ( $file =~ /\.f.a$/ || $file =~ /\.fasta$/);
    }
    
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
    unless ($options{fasta_input} || $options{fasta_list} ||
            $options{fasta_file}  || $options{fasta_dir}) {
        $logger->logdie("You must specify input with --fasta_input --fasta_list --fasta_file or --fasta_dir");
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

    $options{output_subdir_size}   = 0  unless ($options{output_subdir_size});
    $options{output_subdir_prefix} = '' unless ($options{output_subdir_prefix});

    if(0){
        pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});    
    }
}

sub loadMultiSequence {
    #  USAGE:   loadMultiSequence($filepath)
    #  RETURNS: hash
    #
    #  takes a file or path as an argument.  that file should be a multiple-
    #  sequence FASTA file.  It returns a hash with a structure like:
    #      $db{id}{'h'} = header
    #             {'s'} = sequence without whitespace
    #
    #  where id is an incrementing integer that represents that sequence's
    #  order in the file.
    #
    #########################################################################
    my ($file) = @_;
    
    my $seqid = 0;
    my $seq = '';
    my $header;
    my %db;
    
    ## load the sequence file
    open (my $sfh, "<$file") || $logger->logdie("can't open $file because $!");

    for (<$sfh>) {
        ## if we find a header line ...
        if (/^\>(.*)/) {

            $header = $1;

            ## don't do anything if this is the first sequence
            if ($seqid == 0) {
                $seqid++;
                $db{$seqid}{'h'} = $header;
                next;
            } 

            ## remove whitespace
            $seq =~ s/\s//g;
 
            ## record the previous sequence before starting the new one
            $db{$seqid}{'s'} = $seq;

            ## increment the id counter
            $seqid++;

            ## record the new header
            $db{$seqid}{'h'} = $header;

            ## reset the sequence
            $seq = '';

        ## else we've found a sequence line
        } else {
            ## skip it if it is just whitespace
            next if (/^\s*$/);

            ## record this portion of the sequence
            $seq .= $_;
        }
    }
    
    ## don't forget the last sequence
    $seq =~ s/\s//g;
    $db{$seqid}{'s'} = $seq;

    ## close the sequence file
    close $sfh;
    
    return %db;
}

sub write_single_doc {
    my ($doc, $id) = @_;
    my $dirpath = '';
    
    ## the path depends on whether we are using output subdirectories
    if ($options{output_subdir_size}) {
        $dirpath = "$options{'output'}/$options{output_subdir_prefix}$sub_dir_num";
    } else {
        $dirpath = "$options{'output'}";
    }

    ## if the directory doesn't exist, create it.
    mkdir($dirpath) unless (-e $dirpath);
    
    $seq_file_count++;
    $logger->debug("writing file $dirpath/$id.bsml") if ($logger->is_debug);
    $doc->write( "$dirpath/$id.bsml" );

    ## if we are writing multiple subdirectories and have hit our size limit,
    ##  increase the counter to the next one
    if ($options{output_subdir_size} && $options{output_subdir_size} == $seq_file_count) {
        $sub_dir_num++;
        $seq_file_count = 0;
    }

    if ($options{output_list}) {
        print $olist_fh "$dirpath/$id.bsml\n";
    }
}
