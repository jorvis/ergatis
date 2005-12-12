#!/usr/local/bin/perl

=head1  NAME

generate_subflow.pl

=head1 SYNOPSIS

Generate a workflow xml instance file for a iterative subflow.

USAGE: generate_subflow.pl -t subflow template xml -i subflow ini file -o subflow instance file -a iterator template xml -b iterator ini -c conf file [--debug level] [--log file]

=head1 OPTIONS

B<--template, -t> Workflow template XML file for the subflow.

B<--inifile, -i> INI file for the subflow.

B<--outfile, -o> Workflow XML instance file to be created.

B<-a> Workflow template XML for the iterator.  This is usually a very
simple template file that specifies whether the iterator is parallel
or serial.  The Papyrus package includes two examples:
default-paralleliterator_template.xml and
default-serialiterator_template.xml.

B<-b> Workflow template INI for the iterator.  

B<--conf, -c> Config file ini INI format for key/value replacement
pairs.  This replacement is applied to the subflow ini file.  See
run_pipeline for more information about the configuration file.

B<--wfname, -w> Workflow name.  This is used in naming files and
directories.

B<--nodistrib, -n> Run all commands locally

B<--log, -l> [OPTIONAL] log file

B<--help, -h> [OPTIONAL]  program help

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
BEGIN {
    require '/usr/local/devel/ANNOTATION/cas/lib/site_perl/5.8.5/Workflow/Builder.pm';
    import Workflow::Builder;
    require '/usr/local/devel/ANNOTATION/cas/lib/site_perl/5.8.5/Workflow/Logger.pm';
    import Workflow::Logger;
    require '/usr/local/devel/ANNOTATION/cas/lib/site_perl/5.8.5/Workflow/IteratorBuilder.pm';
    import Workflow::IteratorBuilder;
}
use File::Basename;
use Config::IniFiles;

umask(0000);

my %options = ();

my $results = GetOptions (\%options, 
                          'template|t=s', 
                          'inifile|i=s',
                          'outputxml|x=s', 
                          'outputdir|o=s' ,
                          'iteratorlist|l=s', 
                          'iteratortemplate=s', 
                          'iteratorini=s', 
			  'conf|c=s',
			  'nodistrib|n=s',
			  'wfname|w=s',
                          'log=s',
                          'debug=s', 
                          'help|h' );

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

my $wfiteratorobj = new Workflow::IteratorBuilder('NAME'=>$options{'wfname'},
						  'OUTPUT_DIR'=>$options{'outputdir'},
						  'NODISTRIB'=>$options{'nodistrib'}
						  );

my $cfg = new Config::IniFiles( -file => $options{'conf'});
$wfiteratorobj->set_config($cfg);
#generate iterative subflow

my $iteratorconf = &builditeratorconf($options{'iteratorlist'});

$wfiteratorobj->set_iterator_template($options{'iteratorini'},$options{'iteratortemplate'});
my $iteratorinstance = $wfiteratorobj->generate_iterator_instance($options{'inifile'}, 
                                                                  $options{'template'}, 
                                                                  $iteratorconf,
								  $options{'outputxml'});
$logger->get_logger()->debug("Created instance file $iteratorinstance") if($logger->get_logger()->is_debug());

exit;

sub builditeratorconf{
    my ($listfile) = @_;
    my $iteratorconf = {};

    open FILE, $listfile or $logger->get_logger()->logdie("Can't open file $listfile");
    while(my $line=<FILE>){
	chomp $line;
	my($key,$value)=split(/=/,$line);
	my @elts = split(/,/,$value);
	if(! (exists $iteratorconf->{$key})){
	    $iteratorconf->{$key} = [];
	}
	foreach my $elt (@elts){
	    push( @{$iteratorconf->{$key}}, $elt );
	}
    }
    return $iteratorconf;
}
