#!/usr/local/bin/perl
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

=head1  NAME

run_cogs.pl

=head1 SYNOPSIS

Generate cogs search data from a BSML data repository. Each input sequence is searched against all sequences in the repository.

USAGE: run_cogs -d DATABASE_NAME 

=head1 OPTIONS

B<--asmbl_file, -f> [OPTIONAL]   input file containing a list of assembly names, one per line. If no assembly file is specified, run_promer will iterate through all gene model documents in the repository.

B<--database, -d> [OPTIONAL]  database name corresponding to the data repository. If not specified at the command line, this parameter must be included in a custom configuration.

B<--search, -s> [OPTIONAL] search encoding type to be used as input to Cogs. Must be either blastp or ber.

B<--config, -c> [OPTIONAL]  user defined configuration file of the format described using --printconf

B<--printconf, -p> [OPTIONAL]  print the default configuration to STDOUT in the format required for --config 

B<--help, -h> [OPTIONAL]  program help

=cut

my %options = ();
my $results = GetOptions (\%options, 'asmbl_file|f=s', 'database|d=s', 'config|c=s', 'printconf|p', 'search|s=s', 'help|h'  );

# display documentation

if( $options{'help'} )
{
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

# The default COGS configuration

my $WorkflowConfigFileDir = $ENV{'WORKFLOW_DOCS_DIR'};
my $WorkflowSchemaDir = $ENV{'SCHEMA_DOCS_DIR'};
my $WorkflowExecDir = $ENV{'WORKFLOW_WRAPPERS_DIR'};

my $conf = {
    ';DATABASE;' => ' ',
    ';DATABASE_LC;' => ' ',
    ';TMP_DIR;' => '/usr/local/scratch',
    ';REPOSITORY_ROOT;' => '/usr/local/annotation',
    ';INSTALLATION_DIR;' => $WorkflowExecDir,
    ';WORKFLOW_CONFIG_DIR;' => $WorkflowConfigFileDir,
    ';WORKFLOW_SCHEMA_DIR;' => $WorkflowSchemaDir,
    ';MASTER_CONF;' => "$WorkflowConfigFileDir/cogs-master.ini",
    ';MASTER_TEMPLATE;' => "$WorkflowConfigFileDir/cogs-master_template.xml",
    ';BSML_SEARCH;' => 'blastp'
    };

# Override the defaults if the user has specified a config file

if( $options{'config'} )
{
    open( CONFIG, $options{'config'} ) or die "Could not open config file\n";

    while( my $line = <CONFIG> )
    {
	chomp( $line );
	my ($key, $value) = split( '=', $line );
	$conf->{$key} = $value;
    }

    close( CONFIG );
}

if( $options{'database'} )
{
    $conf->{';DATABASE;'} = uc( $options{'database'} );
    $conf->{';DATABASE_LC;'} = lc( $options{'database'} );
}

if( $options{'search'} )
{
    $conf->{';BSML_SEARCH;'} = $options{'search'};
} 

# print out a configuration file to stdout 

if( $options{'printconf'} )
{
    foreach my $key (keys( %{$conf} ))
    {
	print "$key=$conf->{$key}\n";
    } 
    exit(0);
}

# Verify the required configuration options have been specified

if( !( $conf->{';DATABASE;'} && $conf->{';DATABASE_LC;'} ) )
{
    print STDERR "ERROR: User must specify database... \n";
    exit(1);
}

if( !( -e $conf->{';MASTER_CONF;'} ) )
{
    print STDERR "ERROR: Master Workflow Configuration File does not exist ($conf->{'MASTER_CONF'})\n";
    exit(1);
}

if( !( -e $conf->{';MASTER_TEMPLATE;'} ) )
{
    print STDERR "ERROR: Master Workflow Template File does not exist ($conf->{'MASTER_TEMPLATE'})\n";
    exit(1);
}

# Look for the BSML dtd in the specified schema directory

if( !( -e "$conf->{';WORKFLOW_SCHEMA_DIR;'}/bsml3_1.dtd" ))
{
    print STDERR "ERROR: Invalid Workflow Schema Directory ($conf->{';WORKFLOW_SCHEMA_DIR;'})\n";
    exit(1);
}

# Verify that the repository root exists and that Read/Write permissions are enabled

if( !( -w "$conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}" ) || !( -r "$conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}" ) )
{
    print STDERR "ERROR: Repository Root does not exist with read/write permissions ($conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}\n";
    exit(1);
}

# Verify that scratch space exists and that Read/Write permissions are enabled

if( !( -w "$conf->{';TMP_DIR;'}" ) || !( -r "$conf->{';TMP_DIR;'}" ) )
{
    print STDERR "ERROR: Scratch Space does not exist with read/write permissions ($conf->{';TMP_DIR;'})\n";
    exit(1);
}

# Verify that the workflow executables are available

if( !( -x "$conf->{';INSTALLATION_DIR;'}/set_runtime_vars.pl" ) || !( -x "$conf->{';INSTALLATION_DIR;'}/run_wf.sh" ) )
{
    print STDERR "ERROR: Workflow executable directory is not valid ($conf->{';INSTALLATION_DIR;'})\n";
    exit(1);
}

# Verify that the search data type is set to either blastp or ber

if( !( ($conf->{';BSML_SEARCH;'} eq 'blastp') || ($conf->{';BSML_SEARCH;'} eq 'ber')) )
{
    print STDERR "ERROR: Search type must be either blastp or ber ($conf->{';BSML_SEARCH;'})\n";
    exit(1);
}

# Verify that a checksum file exists for the search data

if( !( -r "$conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/BSML_repository/$conf->{';BSML_SEARCH;'}/CHECKSUM" ))
{
    print STDERR "ERROR: Search data not found ($conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/BSML_repository/$conf->{';BSML_SEARCH;'})\n";
    exit(1);
}

# Set the workflow name
$conf->{';WFNAME;'} = 'Cogs';

# Set the workflow id (process id)
$conf->{';WFID;'} = $$;

# Create the workflow archive directory within the BSML repository

if( ! -d "$conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow" )
{
    system "mkdir -m 777 -p $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow";
}

if( ! -d "$conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}" )
{
    system "mkdir -m 777 -p $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}";
}

if( ! -d "$conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}" )
{
    system "mkdir -m 777 -p $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}";
}

# Write the master workflow variable substitution configuration to disk in the repository archive

open( MASTER_SUBS, ">$conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/master.subs" ) or die "Could not open master.subs\n";

foreach my $key (keys( %{$conf} ))
{
    print MASTER_SUBS "$key=$conf->{$key}\n";
}
 
close( MASTER_SUBS );

# Perform variable substitution on the master workflow 

system "$WorkflowExecDir/set_runtime_vars.pl -c $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/master.subs < $conf->{';MASTER_CONF;'} > $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/master.ini";

system "chmod 777 $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/master.ini";
system "chmod 777 $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/master.subs";

system "cp $conf->{';MASTER_TEMPLATE;'} $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/master-template.xml";
system "chmod 777 $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/master-template.xml";

system "$WorkflowExecDir/run_wf.sh -d $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'} -c $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/master.ini -t $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/master-template.xml -i $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/master.xml -l $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/log.txt";

