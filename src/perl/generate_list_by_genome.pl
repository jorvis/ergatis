#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
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

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

use File::Basename;
BEGIN {
use Workflow::Logger;
use BSML::BsmlReader;
use BSML::BsmlParserSerialSearch;
}

umask(0000);

my %options = ();

my $results = GetOptions (\%options, 
                          'filelist|l=s', 
                          'file|f=s',
                          'directory|d=s', 
			  'output_dir|o=s',
			  'extension|x=s',
			  'output_extension|x=s',
			  'log|l=s',
                          'debug=s', 
                          'help|h' ) || pod2usage();

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Workflow::Logger::get_logger();

my @iteratorelts;

$options{'extension'} = 'bsml' if($options{'extension'} eq "");
$options{'output_extension'} = 'list' if($options{'extension'} eq "");

my $keyname = '$;'.uc($options{'extension'}).'_FILE$;';

my $iteratorconf = {
    $keyname            => [],
    '$;SUBFLOW_NAME$;'  => [],
};

if($options{'file'}){
    if(-e $options{'file'} && -f $options{'file'}){
	push( @{$iteratorconf->{$keyname}}, $options{'file'} );
	my $name = &get_name_from_file($options{'file'});
	push( @{$iteratorconf->{'$;SUBFLOW_NAME$;'}}, $name );
    }
    else{
	$logger->logdie("Can't open file $options{'file'}");
    }
}
if($options{'filelist'}){
    &get_list_from_file($iteratorconf,$options{'filelist'});
}
if($options{'directory'}){
    &get_list_from_directory($iteratorconf,$options{'directory'});
}

my $genome;
my $genomelookup = {};
foreach my $file (@{$iteratorconf->{$keyname}}){
    $logger->debug("Parsing file $file");
    my $featParser = new BSML::BsmlParserSerialSearch(GenomeCallBack => \&genomeHandler );
    $featParser->parse($file);
    $logger->debug("Found genome $genome in file $file");
    if(! exists $genomelookup->{$genome}){
	$genomelookup->{$genome} = [];
    }
    push @{$genomelookup->{$genome}},$file;
}

foreach my $g (keys %$genomelookup){
    $logger->debug("Dumping genome $g");
    my $outfile = $options{'output_dir'}."/$g"."."."$options{'output_extension'}";
    open FILE, "+>$outfile"
	or $logger->logdie("Can't open output file $outfile");
    print FILE join("\n",@{$genomelookup->{$g}}),"\n";
    close FILE;
}

exit;
	
    
sub genomeHandler{
    my $bsmlGenome = shift;
    my $reader = new BSML::BsmlReader;
    
    my $rhash = $reader->readGenome( $bsmlGenome );

    if( !defined($rhash->{'strain'}) )
    {
	$rhash->{'strain'} = ' ';
    }

    $genome = $rhash->{'genus'}.'_'.$rhash->{'species'}.'_'.$rhash->{'strain'};
    $genome =~ s/\s//g;
    $genome =~ s/[^\w\.\-\_]/_/g;
}
					     



sub get_list_from_file{
    my ($iteratorconf, $f) = @_;
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
		fisher_yates_shuffle(\@lines);
		foreach my $line (@lines){
		    if($line){
			my $filename = "$line";
			my $name = &get_name_from_file($filename);
			&add_entry_to_conf($iteratorconf,$filename,$name);
		    }
		}
		close( FH );
	    }
	    else{
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
		my $name = &get_name_from_file($filename);
		&add_entry_to_conf($iteratorconf,$filename,$name);
	    }
	}
	else{
	    $logger->logdie("Can't open directory $directory");
	}
    }
}


sub get_name_from_file{
    my($filename) = @_;
    my($name,$path,$suffix) = fileparse($filename,".$options{'extension'}");
    return $name;
}

sub add_entry_to_conf{
    my($iteratorconf,$filename,$name) = @_;
    push( @{$iteratorconf->{$keyname}}, $filename );
    push( @{$iteratorconf->{'$;SUBFLOW_NAME$;'}}, $name );
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
