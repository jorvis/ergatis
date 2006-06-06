#!/usr/local/bin/perl

use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

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
    
=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Workflow::Logger;

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
    '$;TIMESTAMP$;'      => [],
    
    ## this should be removed later
    $keyname            => [],
    '$;SUBFLOW_NAME$;' => [],
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

    my @files = split(',',$f);
    foreach my $file (@files){

        $logger->debug("Processing file '$file'") if $logger->is_debug();

        if ( $file) {
            if (-e $file && -f $file) {
                open( FH, $file ) or $logger->logdie("Could not open $file");

                my @lines;

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


    my $timestamp = &get_timestamp();

    if (! exists $filehash->{$filename} ){
    
        push( @{$iteratorconf->{'$;ITER_FILE_PATH$;'}}, $filename );
        push( @{$iteratorconf->{'$;ITER_DIR$;'}},       $iter_dir );
        push( @{$iteratorconf->{'$;ITER_FILE_NAME$;'}}, $iter_file_name );
        push( @{$iteratorconf->{'$;ITER_FILE_BASE$;'}}, $iter_file_base );
        push( @{$iteratorconf->{'$;ITER_FILE_EXT$;'}},  $iter_file_ext );
        push( @{$iteratorconf->{'$;TIMESTAMP$;'}},      $timestamp );

        ## this should be removed later
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



#---------------------------------------------
# get_timestamp()
#
#---------------------------------------------
sub get_timestamp {

    # perl localtime   = Tue Apr  1 18:31:09 2003 
    # sybase getdate() = Apr  2 2003 10:15AM
    
    my $timestamp = localtime;
    #                  Day of Week                        Month of Year                                       Day of Month  Hour      Mins     Seconds    Year   
    if ($timestamp =~ /^(Mon|Tue|Wed|Thu|Fri|Sat|Sun)[\s]+(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\s+([\d]{1,2})\s+([\d]{2}):([\d]{2}):[\d]{2}\s+([\d]{4})$/){
	my $hour = $4;
	my $ampm = "AM";
	if ($4 ge "13"){
	    
	    $hour = $4 - 12;
	    $ampm = "PM";
	}
	$timestamp = "$2  $3 $6  $hour:$5$ampm";
    }
    else{
	$logger->logdie("Could not parse timestamp");
    }

    return $timestamp;
    
}#end sub get_timestamp()

