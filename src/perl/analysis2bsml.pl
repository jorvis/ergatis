#!/usr/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

analysis2bsml.pl - add analysis parameters to a BSML result doc

=head1 SYNOPSIS

    USAGE:  analysis2bsml.pl --bsml_file=/path/to/some/output.bsml 
                             --conf=/path/to/some/pipeline.config
                           [ --ergatis_control_file=/path/to/ergatis_install.ini ]
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
element created or appended to.  If the value 'first' is passed, the first Analysis element will
be used.  If not given, the Analysis element will be created with no id attribute.

=item *

B<--log,-l> Log file

=item *

B<--help,-h> This help message

=item *

B<--ergatis_control_file,-e> Control file defining software versions for ergatis install.

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
use Ergatis::Logger;
use XML::Twig;
use Config::IniFiles;

my %options = ();
my $results = GetOptions (\%options,
              'analysis_id|a=s',
			  'bsml_file|b=s',
			  'conf|c=s',
			  'log|l=s',
			  'debug=s',
			  'ergatis_control_file|e=s',
			  'componenturl|u=s',
			  'pipelineurl|p=s',
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

## build the Twig.  this will toggle the $research_found below if it works
my $research_found = 0;
my $research_text = '';
my $twig = XML::Twig->new(
                            twig_roots => { 'Research' => \&process_research },
                            pretty_print => 'indented',
                         );

my $fh;
if (-e $options{'bsml_file'} && $options{'bsml_file'} !~ /\.gz$/){
    print "Opening regular file\n";
    open ($fh, $options{'bsml_file'}) or $logger->logdie("Could not open file '$options{'bsml_file'}: $!");
}
else{
    if ($options{'bsml_file'} !~ /\.gz$/ && -e $options{'bsml_file'}.".gz"){
	$options{'bsml_file'} .= ".gz";
    }
    if($options{'bsml_file'} =~ /\.gz$/ && -e $options{'bsml_file'}){
	print "Opening gzip file\n";
	open ($fh, "<:gzip", $options{'bsml_file'}) || $logger->logdie("Could not open zipped file '$options{'bsml_file'}': $!");
    }
    else{
	$logger->logdie("Could not open file '$options{'bsml_file'}: $!");
    }
}

$twig->parse( $fh );

close $fh;

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
   
	## add software version attributes to the Analysis element
	add_software_versions($analysis);
   	
    ## build them together
    $analysis->paste('last_child', $analyses);
    $analyses->paste('last_child', $research);
    
    $research_text = tabbed(\$research->sprint, 1, 2);
}


## open the input file.  if a Research element exists, replace it, else just add it
##  before the closing Bsml tag.  because you can't read and write to the same stream,
##  we'll have to write to a temp file and then move it over the original.
my $ifh;
my $ofh;
if($options{'bsml_file'} =~ /\.gz$/){
    print "Read/write zip file\n";
    open ($ifh, "<:gzip", $options{'bsml_file'}) || $logger->logdie("can't read input BSML file");
    open ($ofh, ">:gzip", "$options{'bsml_file'}.part") || $logger->logdie("can't write output BSML file");
}
else{
    open ($ifh, "<$options{'bsml_file'}") || $logger->logdie("can't read input BSML file");
    open ($ofh, ">$options{'bsml_file'}.part") || $logger->logdie("can't write output BSML file");
}
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
    }elsif(! defined($options{'analysis_id'}) ) {
	    my $analysis = $analyses->first_child('Analysis');
	    if($analysis){
	        $options{'analysis_id'} = $analysis->att('id');
	    }
    }
    
    ## does Analyses contain any Analysis elements?
    my $analysis;
    if ( $analyses->children_count('Analysis') ) {
        
        ## do any of them match the id we were looking for?
        for ( $analyses->children('Analysis') ) {
            if ( $options{analysis_id} eq 'first' ) {
                $analysis = $_;
                $analysis_not_found = 0;
                last;            
            }
        
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

    ## make sure the sourcename is only the path to the output directory, not the output file itself.
    for my $attribute ( $analysis->children('Attribute') ) {
        if ( $attribute->att('name') eq 'sourcename' ) {
            $attribute->set_att('content',dirname(dirname($attribute->att('content'))));
            last;
        }
    }
    
    ## add the Attributes to the Analysis element
    &add_config_params($analysis);
   
	if (defined($options{'ergatis_control_file'})) {
		## add software version attributes to the Analysis element
		add_software_versions($analysis);
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


sub add_software_versions {
    my $analysis = shift;

    ## this is optional for now.  don't die if the file doesn't exist.  this is only appropriate
    ##  when an organization is doing the same sort of software version tracking as we are.  This
    ##  see bug: http://sourceforge.net/tracker/index.php?func=detail&aid=1831859&group_id=148765&atid=772583
    if (! -e $options{'ergatis_control_file'} ) {
        return;
    }
    
	my %software_version = ();
	my %tag_onto = (
					'prok_prism' 	=> 'prok_prism_version',
					'coati' 		=> 'coati_version',
					'euk_prism' 	=> 'euk_prism_version',
					'chado_prism' 	=> 'chado_prism_version',
					'shared_prism' 	=> 'shared_prism_version',
					'bsml' 			=> 'bsml_version',
					'ergatis' 		=> 'ergatis_version',
					'ontologies' 	=> 'ontologies_version',
					'chado_schema' 	=> 'chado_schema_version',
					'cvdata' 		=> 'cvdata_version',
					'peffect'		=> 'peffect_version',
					'server' 		=> 'database_server_name',
				   );

	open (IN, $options{'ergatis_control_file'}) 
		|| $logger->logdie("Could not open $options{ergatis_control_file} for reading"); 

	while (<IN>) {
		chomp;
		s/\s+//g;
		if ($_ eq '' || /^#/) {next;}
		my ($tag, $version) = split('=');
		if (!defined($tag_onto{$tag})) {
			$logger->logdie("No ontology term defined for software tag '$tag'.");
		}
		if ($version eq '') {
			$logger->logdie("Software version string for tag '$tag' is empty. Check control file for errors.");
		}
		$software_version{$tag_onto{$tag}} = $version;
	}

    foreach my $name(keys(%software_version)) {
		my $content = $software_version{$name};

        ## create the new Attribute
        my $attribute = XML::Twig::Elt->new('Attribute');
        $attribute->set_att('name', $name);
        $attribute->set_att('content', $content);
            
        ## add it to the Analysis
	    $attribute->paste('last_child', $analysis);
	}
}

sub check_parameters{
    my ($options) = @_;
    
    if(! -e $options{'bsml_file'} && ! -e $options{'bsml_file'}.".gz"){
	pod2usage({-exitval => 2,  -message => "Can't read bsml file $options{'bsml_file'} or $options{'bsml_file'}.gz", -verbose => 1, -output => \*STDERR});    
    }
    if(! -e $options{'conf'}){
	pod2usage({-exitval => 2,  -message => "Can't read conf file $options{'conf'}", -verbose => 1, -output => \*STDERR});    
    }
    
    ## this is optional for now.  don't die if the file doesn't exist.  this is only appropriate
    ##  when an organization is doing the same sort of software version tracking as we are.  This
    ##  see bug: http://sourceforge.net/tracker/index.php?func=detail&aid=1831859&group_id=148765&atid=772583
#    if(defined $options{ergatis_control_file} && ! -e $options{'ergatis_control_file'}){
#		pod2usage({-exitval => 2,  -message => "Can't read ergatis install control file $options{'ergatis_control_file'}", -verbose => 1, -output => \*STDERR});    
#	}
}
