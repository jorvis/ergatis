#!/usr/local/bin/perl
#instantiate-pipeline.pl - provide a saved workflow pipeline template directory and a destination project directory
#                          and this will create an instance of the pipeline (but not begin execution)

=head1  NAME

instantiate_pipeline - Creates a runnable Workflow pipeline from a saved template, and executes it.

=head1 SYNOPSIS

USAGE: instantiate_pipeline --template-dir=/directory/containing/template_pipeline_xml --repository-root=/repository/directory/to/run/template [--execute]

=head1 OPTIONS

B<--template-dir,-d>
    Directory containing template 'pipeline.xml' to instantiate.

B<--repository-root,-r>
    Repository root in which to create the pipeline instance.

B<--ergatis-install,-e>
    Path to ergatis install that will be used to run workflow.

B<--execute,-x>
    Execute the pipeline after instantiating (if omitted no pipeline execution will occur).
    
B<--man,-m>
    Display the pod2usage page for this utility.

B<--help,-h>
    Print this help.

=head1   DESCRIPTION

This script is used to do something really magical.

=head1 INPUT

Takes a file of some kind.

=head1 OUTPUT

Outputs another file.

=head1 CONTACT

Aaron Gussman
agussman@tigr.org
Brett Whitty
bwhitty@tigr.org

=cut

BEGIN {
use lib '/home/bwhitty/ergatis_test_suite/lib';
}
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
require Workflow::SavedPipeline;

my %options = ();
my $results = GetOptions (\%options,
              'template-dir|d=s',
              'repository-root|r=s',
			  'ergatis-install|e=s',
	      	  'execute|x',
		  	  'workflow-path|w=s',
              'help|h',
	          'man|m',
                         ) || pod2usage();

# display documentation
if( $options{'man'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 1, -output => \*STDOUT} );
}
if (!defined($options{'template-dir'})) {
    pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
}
if (!defined($options{'repository-root'})) {
    pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
}
if (!defined($options{'ergatis-install'})) {
    pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
}

unless (-d $options{'template-dir'}) {
    print "Template directory '".$options{'repository-root'}."' invalid:\n$!\n";
    pod2usage();
}
unless (-d $options{'repository-root'}) {
    print "Repository root '".$options{'repository-root'}."' invalid:\n$!\n";
    pod2usage();
}

if (!defined($options{'workflow-path'})) {
	$options{'workflow-path'} = '/usr/local/devel/ANNOTATION/workflow';
}

$options{'template-dir'} =~ s/\/$//;
$options{'repository-root'} =~ s/\/$//;
$options{'ergatis-install'} =~ s/\/$//; ##strip off trailing slash if it exists

## free love!
umask(0000);

#make a pipeline from a saved template
my $pipe_test = Workflow::SavedPipeline-> new( 
	template => $options{'template-dir'}."/pipeline.xml"
					     );
$pipe_test->write_pipeline( 
	repository_root => $options{'repository-root'} 
                          );
			  
my $pipeline_id = $pipe_test->pipeline_id();

print "\nPipeline_ID of newly exported WorkFlow: $pipeline_id\n";

#create workflow
chdir($options{'repository-root'}."/Workflow/pipeline/".$pipeline_id);
my $pipeline_ini_path = $options{'repository-root'}."/Workflow/pipeline/".$pipeline_id;
chdir($pipeline_ini_path);
if (defined($options{'execute'})) {
#run pipeline
    $ENV{'WF_ROOT'} = $options{'workflow-path'};
    $ENV{'WF_ROOT_INSTALL'} = $options{'workflow-path'};

my $workflow_root = $options{'workflow-path'};
	
#    doORdie($options{'workflow-path'}.'/CreateWorkflow -t pipeline.xml -c pipeline.xml.ini -i pipeline.xml.instance --autobuild=false');
#    doORdie($options{'workflow-path'}.'/RunWorkflow -i pipeline.xml.instance');
	doORdie($options{'ergatis-install'}."/bin/run_wf --workflow_root $workflow_root --template $pipeline_ini_path/pipeline.xml --ini $pipeline_ini_path/pipeline.xml.ini --instance $pipeline_ini_path/pipeline.xml.instance");
	print ($options{'ergatis-install'}."/bin/run_wf --workflow_root $workflow_root --template $pipeline_ini_path/pipeline.xml --ini $pipeline_ini_path/pipeline.xml.ini --instance $pipeline_ini_path/pipeline.xml.instance\n");
}

exit();

sub doORdie {
    my $debugging = 0;
    my $cmd = shift;
    
    print "$cmd\n" if $debugging;
    
    system($cmd);
    
    if ($? == -1) {
       print STDERR "failed to execute: $!\n";
       exit(1);
    } elsif ($? & 127) {
       printf STDERR "child died with signal %d, %s coredump\n",
           ($? & 127),  ($? & 128) ? 'with' : 'without';
       exit(1);
    }

}
