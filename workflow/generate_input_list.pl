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

=back

=head1   DESCRIPTION

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

use Workflow::Logger;
use File::Basename;

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

$options{'extension'} = 'bsml' if($options{'extension'} eq "");

my $keyname = '$;'.uc($options{'extension'}).'_FILE$;';

my $iteratorconf = {$keyname=>[],
		    '$;SUBFLOW_NAME$;'=>[]};

if($options{'file'}){
    push( @{$iteratorconf->{$keyname}}, $options{'file'} );
    my $name = &get_name_from_file($options{'file'});
    push( @{$iteratorconf->{'$;SUBFLOW_NAME$;'}}, $name );
}
if($options{'filelist'}){
    &get_list_from_file($iteratorconf,$options{'filelist'});
}
if($options{'directory'}){
    &get_list_from_directory($iteratorconf,$options{'directory'});
}

if($options{'listfiles'}){
    if($options{'output'}){
	open FILE, "+>$options{'output'}" or $logger->logdie("Can't open output file $options{'output'}");
	print FILE join("\n",@{$iteratorconf->{$keyname}}),"\n";
	close FILE;
    }
    else{
	print STDOUT join("\n",@{$iteratorconf->{$keyname}}),"\n";
    }
}
else{
    if($options{'output'}){
	open FILE, "+>$options{'output'}" or $logger->logdie("Can't open output file $options{'output'}");
	foreach my $key (keys %$iteratorconf){
	    print FILE "$key=",join(',',@{$iteratorconf->{$key}}),"\n";
	}
	close FILE;
    }
    else{
	foreach my $key (keys %$iteratorconf){
	    print FILE "$key=",join(',',@{$iteratorconf->{$key}}),"\n";
	}
    }
}
exit;
						     



sub get_list_from_file{
    my ($iteratorconf, $f) = @_;
    my @elts;
    my @lines;
    my @files = split(',',$f);
    foreach my $file (@files){
	if( $file){
	    open( FH, $file ) or die "Could not open $file";
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
    }
}

sub get_list_from_directory{
    my ($iteratorconf, $dir, $glob) = @_;

    my @directories = split(',',$dir);
    foreach my $directory (@directories){
	opendir DIR, "$directory" or $logger->logdie("Can't read directory $directory");
	my @files = grep /\.$options{'extension'}$/, readdir DIR;
	fisher_yates_shuffle( \@files );    # permutes @array in place
	foreach my $file (@files ){
	    my $filename = "$directory/$file";
	    my $name = &get_name_from_file($filename);
	    &add_entry_to_conf($iteratorconf,$filename,$name);
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
    for ($i = @$array; --$i; ) {
	my $j = int rand ($i+1);
	next if $i == $j;
	@$array[$i,$j] = @$array[$j,$i];
    }
}
