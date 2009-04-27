#!/usr/bin/perl

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

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

use File::Basename;
use Ergatis::Logger;
use XML::Parser;

umask(0000);

my %options = ();

my $results = GetOptions (\%options, 
                          'input_file_list|n=s', 
                          'input_file|f=s',
                          'input_directory|d=s', 
                          'output_dir|o=s',
                          'extension|x=s',
                          'output_extension|u=s',
                          'output_iter_list|l=s',
                          'log|l=s',
                          'debug=s', 
                          'help|h' ) || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                  'LOG_LEVEL'=>$options{'debug'});
$logger = Ergatis::Logger::get_logger();

$options{'extension'} = 'bsml' if($options{'extension'} eq "");
$options{'output_extension'} = 'list' if($options{'output_extension'} eq "");

my $filelist = [];

if ($options{'input_file'}) {
    if (-e $options{'input_file'} && -f $options{'input_file'}) {
        push( @{$filelist}, $options{'input_file'} );
        
    } else {
        $logger->logdie("Can't open file $options{'input_file'}");
    }
}
if ($options{'input_file_list'}) {
    &get_list_from_file($filelist,$options{'input_file_list'});
}
if ($options{'input_directory'}) {
    &get_list_from_directory($filelist,$options{'input_directory'});
}

my $genome;
my $genomelookup = {};

my $funcs = {'Organism'=>
         sub {
                 my ($expat,$elt,%params) = @_;
                 if(exists $params{'genus'} && $params{'species'}){
                 $genome = $params{'genus'}.'_'.$params{'species'};
                 $genome =~ s/\s//g;
                 $genome =~ s/[^\w\.\-\_]/_/g;
             }
         }
     };

foreach my $file (@{$filelist}){
    $logger->info("parsing file $file");
    my $x = new XML::Parser(Handlers => 
                {
                Start =>
                    sub {
                        #$_[1] is the name of the element
                        if(exists $funcs->{$_[1]}){
                            $funcs->{$_[1]}(@_);
                        }
                    }
	    }
			    );
    if (!(-e $file) && -e "$file.gz") {
        $file .= ".gz";
    }
    if(-e $file){
	my $ifh;
	if ($file =~ /\.(gz|gzip)$/) {
	    open ($ifh, "<:gzip", $file) || $logger->logdie("can't read input file $file: $!");
	} else {
	    open ($ifh, "<$file") || $logger->logdie("can't read input file $file: $!");
	}
	$x->parse( $ifh );
	close $ifh;
    }
    else{
	$logger->logdie("Can't read jaccard bsml file $file");
    }		    

    $logger->debug("Found genome $genome in file $file");
    if(! exists $genomelookup->{$genome}){
        $genomelookup->{$genome} = [];
    }
    push @{$genomelookup->{$genome}},$file;
}

open LISTFILE, ">$options{'output_iter_list'}" or $logger->logdie("Can't open $options{'output_iter_list'}");
print LISTFILE '$;I_GENOME$;',"\t",'$;I_GENOME_LIST_FILE$;',"\n";
foreach my $g (keys %$genomelookup){
    $logger->debug("Dumping genome $g");
    my $outfile = $options{'output_dir'}."/$g"."."."$options{'output_extension'}";
    print LISTFILE "$g\t$outfile\n";
    open FILE, "+>$outfile" or $logger->logdie("Can't open output file $outfile");
    print FILE join("\n",@{$genomelookup->{$g}}),"\n";
    close FILE;
}

close LISTFILE;

exit;
    
    
sub get_list_from_file{
    my ($filelist, $f) = @_;
    my @elts;
    my @lines;
    my @files = split(',',$f);
    foreach my $file (@files){
        if( $file){
            if(-e $file && -f $file){
                open( FH, $file ) or $logger->logdie("Could not open $file");
                while( my $line = <FH> ){
                    chomp($line);
                    push @lines,  split(',',$line) if($line =~ /\S+/);
                }

                foreach my $line (@lines){
                    if($line){
                        my $filename = "$line";
                        push( @{$filelist}, $filename );
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
    my ($filelist, $dir, $glob) = @_;

    my @directories = split(',',$dir);
    foreach my $directory (@directories) {
        if (-e $directory && -d $directory) {
            opendir DIR, "$directory" or $logger->logdie("Can't read directory $directory");
            my @files = grep /\.$options{'extension'}$/, readdir DIR;
            foreach my $file (@files ){
                my $filename = "$directory/$file";
                push( @{$filelist}, $filename );
            }
        } else {
            $logger->logdie("Can't open directory $directory");
        }
    }
}


