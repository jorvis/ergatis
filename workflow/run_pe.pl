#!/usr/local/bin/perl
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

=head1  NAME

run_pe.pl

=head1 SYNOPSIS

Generate PEffect search data from a BSML data repository. Each input sequence is searched against all sequences in the repository.

USAGE: run_pe -d DATABASE_NAME 

=head1 OPTIONS

B<--asmbl_file, -f> [OPTIONAL]   input file containing a list of assembly names, one per line. If no assembly file is specified, run_promer will iterate through all gene model documents in the repository.

B<--database, -d> [OPTIONAL]  database name corresponding to the data repository. If not specified at the command line, this parameter must be included in a custom configuration.

B<--search, -s> [OPTIONAL] specifies the input search type. Must be either blastp or ber. 

B<--config, -c> [OPTIONAL]  user defined configuration file of the format described using --printconf

B<--printconf, -p> [OPTIONAL]  print the default configuration to STDOUT in the format required for --config 

B<--help, -h> [OPTIONAL]  program help

=cut

#PEFFECT_WINDOW_SIZE=10
#PEFFECT_GAP_PENALTY=-50
#PEFFECT_ORIENTATION_PENALTY=-100
#PEFFECT_MIN_MATCHES_PER_WINDOW=4

#REGIONS_MIN_CLUSTER=3000
#REGIONS_MIN_GAP=10000
#point to pod (man gene_boundaries_bsml.pl)

my %options = ();
my $results = GetOptions (\%options, 'asmbl_file|f=s', 'database|d=s', 'config|c=s', 'printconf|p', 'search|s=s', 'help|h' );

# display documentation

if( $options{'help'} )
{
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

# The default pe configuration

my $WorkflowConfigFileDir = $ENV{'WORKFLOW_DOCS_DIR'};
my $WorkflowSchemaDir = $ENV{'SCHEMA_DOCS_DIR'};
my $WorkflowExecDir = $ENV{'WORKFLOW_WRAPPERS_DIR'};

my $conf = {
    ';TMP_DIR;' => '/usr/local/scratch',
    ';REPOSITORY_ROOT;' => '/usr/local/annotation',
    ';MASTER_CONF;' => "$WorkflowConfigFileDir/pe-master.ini",
    ';MASTER_TEMPLATE;' => "$WorkflowConfigFileDir/pe-master_template.xml",
    ';SUBFLOW_CONF;' => "$WorkflowConfigFileDir/pe.ini", 
    ';SUBFLOW_TEMPLATE;' => "$WorkflowConfigFileDir/pe_template.xml",
    ';INPUT_SEARCH_TYPE;' => 'ber',
    ';INSTALLATION_DIR;' => $WorkflowExecDir,
    ';WORKFLOW_CONFIG_DIR;' => $WorkflowConfigFileDir,
    ';WORKFLOW_SCHEMA_DIR;' => $WorkflowSchemaDir,
    ';PEFFECT_WINDOW_SIZE;' => 10,
    ';PEFFECT_GAP_PENALTY;' => -50,
    ';PEFFECT_ORIENTATION_PENALTY;' => -100,
    ';PEFFECT_MIN_MATCHES_PER_WINDOW;' => 4,
    ';REGIONS_MIN_CLUSTER;' => 3000,
    ';REGIONS_MIN_GAP;' => 10000
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
    $conf->{';INPUT_SEARCH_TYPE;'} = $options{'search'};
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
    print STDERR "User must specify database... \n";
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

if( !( -e $conf->{';SUBFLOW_CONF;'} ) )
{
    print STDERR "ERROR: Sub-Workflow Configuration File does not exist ($conf->{'SUBFLOW_CONF'})\n";
    exit(1);
}

if( !( -e $conf->{';SUBFLOW_TEMPLATE;'} ) )
{
    print STDERR "ERROR: Sub-Workflow Template File does not exist ($conf->{'SUBFLOW_CONF'})\n";
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

if( !( ($conf->{';INPUT_SEARCH_TYPE;'} eq 'blastp') || ($conf->{';INPUT_SEARCH_TYPE;'} eq 'ber')) )
{
    print STDERR "ERROR: Search type must be either blastp or ber ($conf->{';INPUT_SEARCH_TYPE;'})\n";
    exit(1);
}

# Verify that a checksum file exists for the search data

if( !( -r "$conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/BSML_repository/$conf->{';INPUT_SEARCH_TYPE;'}/CHECKSUM" ))
{
    print STDERR "ERROR: Search data not found ($conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/BSML_repository/$conf->{';INPUT_SEARCH_TYPE;'})\n";
    exit(1);
}

# Set the workflow name
$conf->{';WFNAME;'} = 'PEffect';

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

# Determine the workflow iteration lists and store them in the substitution config file

my @asmbl_list = ();

if( $options{'asmbl_file'} )
{
    open( IDS, $options{'asmbl_file'} ) or die "Could not open $options{'asmbl_file'}";
    while( my $asmbl = <IDS> )
    {
	chomp($asmbl);
	push( @asmbl_list, $asmbl );
    }
    close( IDS );
}
else
{
    my $bsml_path = "$conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/BSML_repository";
    my @asmbl_glob = <$bsml_path/*.bsml>;
    foreach my $asmbl (@asmbl_glob )
    {
	$asmbl =~ s/$bsml_path\///;
	$asmbl =~ s/\.bsml//;
	push( @asmbl_list, $asmbl );
    }
}

foreach my $asmbl (@asmbl_list)
{
    $conf->{';CONFIG_LIST;'} .= "$conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/$asmbl.ini, ";
    $conf->{';INSTANCE_LIST;'} .= "$conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/$asmbl.xml, ";
}

chop $conf->{';CONFIG_LIST;'};
chop $conf->{';INSTANCE_LIST;'};
chop $conf->{';CONFIG_LIST;'};
chop $conf->{';INSTANCE_LIST;'};

# Write the master workflow variable substitution configuration to disk in the repository archive

open( MASTER_SUBS, ">$conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/master.subs" ) or die "Could not open master.subs\n";

foreach my $key (keys( %{$conf} ))
{
    print MASTER_SUBS "$key=$conf->{$key}\n";
}
 
close( MASTER_SUBS );

# Perform variable substitution on the master workflow 

system "$WorkflowExecDir/set_runtime_vars.pl -c $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/master.subs < $conf->{';MASTER_CONF;'} > $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/master.ini";

system "chmod 666 $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/master.ini";

system "$WorkflowExecDir/set_runtime_vars.pl -c $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/master.subs < $conf->{';MASTER_TEMPLATE;'} > $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/master-template.xml";

system "chmod 666 $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/master-template.xml";

system "chmod 666 $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/master.subs";

# iteratively replace the asmbl id in the subflow substitution file and create the workflow configurations 
# for each assembly

foreach my $asmbl (@asmbl_list)
{
    $conf->{';ASMBL;'} = $asmbl;

    open( SUBFLOW_SUBS, ">$conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/subflow.subs" ) or die "Could not open master.subs\n";

    foreach my $key (keys( %{$conf} ))
    {
	print SUBFLOW_SUBS "$key=$conf->{$key}\n";
    }
    
    close( SUBFLOW_SUBS );
    
    system "$WorkflowExecDir/set_runtime_vars.pl -c $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/subflow.subs < $conf->{';SUBFLOW_CONF;'} > $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/$asmbl.ini";

    system "chmod 666 $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/$asmbl.ini";
    system "chmod 666 $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/subflow.subs";
}

system "cp $conf->{';SUBFLOW_TEMPLATE;'} $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/pe_template.xml";
system "chmod 666 $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/pe_template.xml";

system "$WorkflowExecDir/run_wf.sh -d $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'} -c $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/master.ini -t $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/master-template.xml -i $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/master.xml -l $conf->{';REPOSITORY_ROOT;'}/$conf->{';DATABASE;'}/Workflow/$conf->{';WFNAME;'}/$conf->{';WFID;'}/log.txt";

