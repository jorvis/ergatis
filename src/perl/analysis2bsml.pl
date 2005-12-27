#!/usr/local/bin/perl

=head1  NAME 

analysis2bsml.pl - add analysis parameters to a BSML result doc

=head1 SYNOPSIS

    USAGE:  analysis2bsml.pl --bsml_file=/path/to/some/output.bsml 
                             --conf=/path/to/some/pipeline.config 
                           [ --componenturl=/path/to/some/pipeline.xml ]
                           [ --pipelineurl=/path/to/some/pipeline0.xml ]
                           [ --debug debug_level  ]
                           [ --log log_file ]
                           [ --analysis_id=some_id ]

=head1 OPTIONS

=item *

B<--debug,-d> Debug level.  Use a large number to turn on verbose debugging. 

=item *

B<--analysis_id,-a> Analysis id (optional).  If given, will be the id attribute for the Analysis
element created or appended to.  If not given, the Analysis element will be created with
no id attribute.

=item *

B<--log,-l> Log file

=item *

B<--help,-h> This help message

=item *

B<--componenturl,-u> Stores component url in <Analysis>

=item *

B<--pipelineurl,-p> Store pipeline.xml url in <Analysis>

=back

=head1   DESCRIPTION

This script is used to record the parameters of an analysis into the BSML file that the
analysis created.  This is done to ensure that each of the run-time settings for the
analysis are recorded for informational purposes or to re-run the analysis.

=cut


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use File::Basename;

BEGIN {
    require '/usr/local/devel/ANNOTATION/cas/lib/site_perl/5.8.5/Workflow/Logger.pm';
    import Workflow::Logger;
}
use XML::Twig;
use Config::IniFiles;

my %options = ();
my $results = GetOptions (\%options,
              'analysis_id|a=s',
			  'bsml_file|b=s',
			  'conf|c=s',
			  'log|l=s',
			  'debug=s',
			  'componenturl|u=s',
			  'pipelineurl|p=s',
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

## build the Twig.  this will toggle the $research_found below if it works
my $research_found = 0;
my $research_text = '';
my $twig = XML::Twig->new(
                            twig_roots => { 'Research' => \&process_research },
                            pretty_print => 'indented',
                         );
$twig->parsefile($options{'bsml_file'});

## if no Research element was found, build one
if (! $research_found) {
    undef $twig;  ## might as well

    my $research = XML::Twig::Elt->new('Research');
    my $analyses = XML::Twig::Elt->new('Analyses');
    my $analysis = XML::Twig::Elt->new('Analysis');
    
    ## add an ID on analysis if it was passed
    if ( defined($options{'analysis_id'}) ) {
        $analysis->set_att('id', $options{'analysis_id'});
    }
    
    ## add the Attributes to the Analysis element
    &add_config_params($analysis);
    
    ## build them together
    $analysis->paste('last_child', $analyses);
    $analyses->paste('last_child', $research);
    
    $research_text = tabbed(\$research->sprint, 1, 2);
}


## open the input file.  if a Research element exists, replace it, else just add it
##  before the closing Bsml tag.  because you can't read and write to the same stream,
##  we'll have to write to a temp file and then move it over the original.
open (my $ifh, "<$options{'bsml_file'}") || $logger->logdie("can't read input BSML file");
open (my $ofh, ">$options{'bsml_file'}.part") || $logger->logdie("can't write output BSML file");

my $replace_mode = 0;
my $research_not_found = 1;

## loop through the file
for (<$ifh>) {
   
    if (m|<Research>|) {
        $replace_mode = 1;
        $research_not_found = 0;
        print $ofh "$research_text\n";
        
    } elsif (m|</Research>|) {
        $replace_mode = 0;
    
    } elsif (m|</Bsml>| && $research_not_found) {
        ## if you hit a </Bsml> and never found <Research>, print the Research
        ##  text before the closing </Bsml>
        print $ofh "$research_text\n";
        print $ofh $_;
        
    } elsif (! $replace_mode) {
        print $ofh $_;
    }
}

## close your files, then move the .part over the original
close $ifh;
close $ofh;
system("mv $options{'bsml_file'}.part $options{'bsml_file'}") && $logger->logdie("failed to move temporary .part file over source");


sub tabbed {
    my ($txtref, $tablevels, $tabspacing) = @_;
    my @txt = split("\n", $$txtref);
    my $result = '';
    
    for (@txt) {
        $result .= (' ' x ($tablevels * $tabspacing)) . "$_\n";
    }
    
    return $result;
}

sub process_research {
    my ($twig, $research) = @_;
    
    $research_found = 1;
    my $analyses_not_found = 0;
    my $analysis_not_found = 1;
    
    ## grab the Analyses element, if there was one
    my $analyses = $research->first_child('Analyses');
    
    ## if there wasn't an Analyses element, create one
    if (! $analyses) {
        $analyses = XML::Twig::Elt->new('Analyses');
        $analyses_not_found = 1;
    }    
    
    ## does Analyses contain any Analysis elements?
    my $analysis;
    if ( $analyses->children_count('Analysis') ) {
        
        ## do any of them match the id we were looking for?
        for ( $analyses->children('Analysis') ) {
            if ( $_->att('id') eq $options{'analysis_id'} ) {
                $analysis = $_;
                $analysis_not_found = 0;
                last;
            }
        }
    }
    
    ## add an Analysis element if we didn't find the one we were looking for
    if ($analysis_not_found) {
        $analysis = XML::Twig::Elt->new('Analysis');
        
        if ( defined($options{'analysis_id'}) ) {
            $analysis->set_att('id', $options{'analysis_id'});
        }
    }
    
    ## add the Attributes to the Analysis element
    &add_config_params($analysis);

    my $sourcename = $analysis->first_child('Attribute',"*[\@name=\"sourcename\"]");
    if($sourcename){
	$sourcename->set_att('content',dirname(dirname($sourcename->att('content'))));
    }
    
    ## add the analysis element to the Analyses element
    $analysis->paste('last_child', $analyses) if ($analysis_not_found);
    
    ## add the analyses element to the Research element
    $analyses->paste('last_child', $research) if ($analyses_not_found);
    
    ## save the text of the Research element
    $research_text = $research->sprint;

}

sub add_config_params {
    my $analysis = shift;
    
    my $cfg = new Config::IniFiles( -file => $options{'conf'});
    my @sections = $cfg->Sections();
    
    for my $section (@sections) {
        my @parameters = $cfg->Parameters($section);
        
        for my $param (@parameters) {
	        my $value = $cfg->val($section,$param);

	        $param =~ s/\$;//g;
	        $param = lc($param);
	        
            if($value ne ""){
                ## create the new Attribute
                my $attribute = XML::Twig::Elt->new('Attribute');
                $attribute->set_att('name', $param);
                $attribute->set_att('content', $value);
                
                ## add it to the Analysis
	            $attribute->paste('last_child', $analysis);
	        }
        }
    }
    
    for my $at ( 'componenturl', 'pipelineurl' ) {
        if ( defined($options{$at}) ){
            my $attribute = XML::Twig::Elt->new('Attribute');
            $attribute->set_att( 'name', $at );
            $attribute->set_att( 'content', $options{$at} );
            $attribute->paste('last_child', $analysis);
        }
    }
}

sub check_parameters{
    my ($options) = @_;
    
    if(! -e $options{'bsml_file'}){
	pod2usage({-exitval => 2,  -message => "Can't read bsml file $options{'bsml_file'}", -verbose => 1, -output => \*STDERR});    
    }
    if(! -e $options{'conf'}){
	pod2usage({-exitval => 2,  -message => "Can't read conf file $options{'conf'}", -verbose => 1, -output => \*STDERR});    
    }
}
