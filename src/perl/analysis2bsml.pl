#!/usr/local/bin/perl

=head1  NAME 

dummy.pl - do nothing

=head1 SYNOPSIS

USAGE:  dummy.pl --debug debug_level --log log_file

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
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Workflow::Logger;
use BSML::BsmlReader;
use BSML::BsmlParserTwig;
use Config::IniFiles;

my %options = ();
my $results = GetOptions (\%options,
			  'bsml_file|b=s',
			  'conf|c=s',
			  'log|l=s',
			  'debug=s',
			  'help|h') || pod2usage();

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Workflow::Logger::get_logger();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDERR} );
}


&check_parameters(\%options);

my $reader = BSML::BsmlReader->new();
my $parser = BSML::BsmlParserTwig->new();

$parser->parse( \$reader, $options{'bsml_file'} );

my $analysis = $reader->returnBsmlAnalysisR( $reader->addBsmlAnalysis() );
my $cfg = new Config::IniFiles( -file => $options{'conf'});
my @sections = $cfg->Sections();
foreach my $section (@sections){
    my @parameters = $cfg->Parameters ($section);
    foreach my $param (@parameters){
	my $value = $cfg->val($section,$param);
	my $delimeter = '$;';
	$param =~ s/\$;//g;
	$param = lc($param);
	if($value ne ""){
	    $analysis->addBsmlAttr( $param, $value );
	}
    }
}

$reader->write( $options{'bsml_file'} );

sub check_parameters{
    my ($options) = @_;
    
    if(! -e $options{'bsml_file'}){
	pod2usage({-exitval => 2,  -message => "Can't read bsml file $options{'bsml_file'}", -verbose => 1, -output => \*STDERR});    
    }
    if(! -e $options{'conf'}){
	pod2usage({-exitval => 2,  -message => "Can't read conf file $options{'conf'}", -verbose => 1, -output => \*STDERR});    
    }
}
