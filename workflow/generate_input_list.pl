#!/usr/local/bin/perl

=head1  NAME 

generate_input_list.pl - Default output is a workflow iterator that
can be used to iterator over a set of files

=head1 SYNOPSIS

USAGE:  generate_input_list

=head1 OPTIONS

=item *

B<--debug,-d> Debug level.  Use a large number to turn on verbose debugging. 

=item *

B<--log,-l> Log file

=item *

B<--help,-h> This help message

=head1   DESCRIPTION

Creates the iterator list file to setup up the iterator.  Variables set and
accessible on each iteration are as follows:

    example input: /path/to/somefile.txt
    
    $;ITER_FILE_PATH$; = /path/to/somefile.txt
    $;ITER_DIR$;       = /path/to
    $;ITER_FILE_NAME$; = somefile.txt
    $;ITER_FILE_BASE$; = somefile
    $;ITER_FILE_EXT$;  = txt
    
the following are still available, but have been deprecated and will be
removed at a certain point.

    $;FSA_FILE$;     = /path/to/somefile.txt (variable name changes based on extension)
    $;SUBFLOW_NAME$; = somefile

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Data::Dumper;
BEGIN {
    require '/usr/local/devel/ANNOTATION/cas/lib/site_perl/5.8.5/Workflow/Logger.pm';
    import Workflow::Logger;
}

umask(0000);

my %options = ();

my $results = GetOptions (\%options, 
                          'filelist|l=s', 
                          'file|f=s',
                          'directory|d=s', 
                          'output|o=s',
                          'randomize|r',
                          'extension|x=s',
                          'listfiles',
                          'log|l=s',
                          'debug=s', 
                          'help|h' ) || pod2usage();

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
                  'LOG_LEVEL'=>$options{'debug'});
$logger = Workflow::Logger::get_logger();

my @iteratorelts;

$options{'extension'} = 'bsml' if ($options{'extension'} eq "");

my $keyname = '$;' . uc($options{'extension'}) . '_FILE$;';

my $filehash = {};


my $iteratorconf = {
    '$;ITER_FILE_PATH$;' => [],
    '$;ITER_FILE_BASE$;' => [],
    '$;ITER_FILE_EXT$;'  => [],
    '$;ITER_FILE_NAME$;' => [],
    '$;ITER_DIR$;'       => [],
    
    ## the following two are deprecated and should be removed later:
    $keyname            => [],
    '$;SUBFLOW_NAME$;'  => [],
};

if ($options{file}) {
    if(-e $options{file} && -f $options{file}){
        &add_entry_to_conf($iteratorconf, $options{file});

    } else {
        $logger->logdie("Can't open file $options{file}");
    }
}

if ($options{'filelist'}) {
    &get_list_from_file($iteratorconf,$options{'filelist'});
}

if ($options{'directory'}) {
    &get_list_from_directory($iteratorconf,$options{'directory'});
}

if ($options{'listfiles'}) {
    if($options{'output'}){
        open FILE, "+>$options{'output'}" or $logger->logdie("Can't open output file $options{'output'}");
        print FILE join( "\n", @{$iteratorconf->{$keyname}} ), "\n";
        close FILE;
    } else {
        print STDOUT join( "\n", @{$iteratorconf->{$keyname}} ), "\n";
    }
    
} else {
    if($options{'output'}){
        open FILE, "+>$options{'output'}" or $logger->logdie("Can't open output file $options{'output'}");

        foreach my $key (keys %$iteratorconf){
            print FILE "$key=",join( ',', @{$iteratorconf->{$key}} ), "\n";
        }

        close FILE;
    } else {
        foreach my $key (keys %$iteratorconf){
            print "$key=",join( ',', @{$iteratorconf->{$key}} ), "\n";
        }
    }
}
exit;
                             



sub get_list_from_file{
    my ($iteratorconf, $f) = @_;
    my @elts;



    #
    # Given the following sample invocation of generate_input_list.pl:
    #
    # perl -I ../lib/ generate_input_list.pl --directory='' --file='' --filelist='/usr/local/scratch/annotation/CHADO_TEST2/BSML_repository/legacy2bsml/cea2.bsml.1.list,/usr/local/scratch/annotation/CHADO_TEST2/BSML_repository/legacy2bsml/cea2.bsml.2.list' --output=subflow1.list
    #
    #
    # Sample filelist stored in $f:
    # /usr/local/scratch/annotation/CHADO_TEST2/BSML_repository/legacy2bsml/cea2.bsml.1.list,/usr/local/scratch/annotation/CHADO_TEST2/BSML_repository/legacy2bsml/cea2.bsml.2.list
    #
    my @files = split(',',$f);
    foreach my $file (@files){

        $logger->debug("Processing file '$file'") if $logger->is_debug();

        if ( $file) {
            if (-e $file && -f $file) {
                open( FH, $file ) or $logger->logdie("Could not open $file");


                #
                # editor:     sundaram@tigr.org
                # date:       2005-08-12
                # bgzcase:    2041 
                # URL:        http://serval.tigr.org:8080/bugzilla/show_bug.cgi?id=2041
                # comment:    @lines array should be locally scoped to ensure no files are output more than once!
                #
                my @lines;

                #
                # Contents of sample file /usr/local/scratch/annotation/CHADO_TEST2/BSML_repository/legacy2bsml/cea2.bsml.1.list
                #
                # /usr/local/scratch/annotation/CHADO_TEST2/BSML_repository/legacy2bsml/cea2_2_assembly.bsml
                # /usr/local/scratch/annotation/CHADO_TEST2/BSML_repository/legacy2bsml/cea2_1_assembly.bsml
                # /usr/local/scratch/annotation/CHADO_TEST2/BSML_repository/legacy2bsml/cea2_3_assembly.bsml
                #
                #
                # Contents of sample file /usr/local/scratch/annotation/CHADO_TEST2/BSML_repository/legacy2bsml/cea2.bsml.2.list
                #
                #
                # /usr/local/scratch/annotation/CHADO_TEST2/BSML_repository/legacy2bsml/cea2_4_assembly.bsml
                # /usr/local/scratch/annotation/CHADO_TEST2/BSML_repository/legacy2bsml/cea2_5_assembly.bsml
                # /usr/local/scratch/annotation/CHADO_TEST2/BSML_repository/legacy2bsml/cea2_6_assembly.bsml
                #


                while( my $line = <FH> ){
                    chomp($line);
                    $logger->debug("line '$line'") if $logger->is_debug();
                    push @lines,  split(',', $line) if ($line =~ /\S+/);
                }


                &fisher_yates_shuffle(\@lines);

                foreach my $line (@lines){

                    if($line){
                        my $filename = "$line";
                        $logger->debug("filename '$filename'") if $logger->is_debug();

                        &add_entry_to_conf($iteratorconf, $filename);
                    }
                }
                close( FH );
            } else {
                $logger->logdie("Can't open list file $file");
            }

        }
    }
}

sub get_list_from_directory{
    my ($iteratorconf, $dir, $glob) = @_;

    my @directories = split(',',$dir);
    foreach my $directory (@directories){
        if(-e $directory && -d $directory){
            opendir DIR, "$directory" or $logger->logdie("Can't read directory $directory");
            my @files = grep /\.$options{'extension'}$/, readdir DIR;

            fisher_yates_shuffle( \@files );    # permutes @array in place

            foreach my $file (@files ){
                my $filename = "$directory/$file";
                &add_entry_to_conf($iteratorconf, $filename);
            }
        }
        else{
            $logger->logdie("Can't open directory $directory");
        }
    }
}

sub add_entry_to_conf{
    my ($iteratorconf, $filename) = @_;
    my ($iter_dir, $iter_file_name, $iter_file_base, $iter_file_ext) = &parse_file_parts($filename);
    #
    # editor:     sundaram@tigr.org
    # date:       2005-08-12
    # bgzcase:    2041
    # URL:        http://serval.tigr.org:8080/bugzilla/show_bug.cgi?id=2041
    # comment:    This script should output a unique file list
    #
    if (! exists $filehash->{$filename} ){
    
        push( @{$iteratorconf->{'$;ITER_FILE_PATH$;'}}, $filename );
        push( @{$iteratorconf->{'$;ITER_DIR$;'}},       $iter_dir );
        push( @{$iteratorconf->{'$;ITER_FILE_NAME$;'}}, $iter_file_name );
        push( @{$iteratorconf->{'$;ITER_FILE_BASE$;'}}, $iter_file_base );
        push( @{$iteratorconf->{'$;ITER_FILE_EXT$;'}},  $iter_file_ext );
    
        ## the following two are deprecated and should be removed at some
        ##  later point when all the components are updated.
        push( @{$iteratorconf->{$keyname}}, $filename );
        push( @{$iteratorconf->{'$;SUBFLOW_NAME$;'}}, $iter_file_base );

        $filehash->{$filename}++;
    }
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
    if ( $file =~ m|^(.+)/(([^/]+)\.(.+))$| ) {
        return ( $1, $2, $3, $4 );

    ## match here if the file doesn't have an extension
    } elsif ( $file =~ m|^(.+)/(([^/]+))$| ) {
        return ( $1, $2, $3, '' );
        
    ## else die
    } else {
        die "\tfile $file doesn't match the regex.  can't extract file name parts\n";
    }
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
