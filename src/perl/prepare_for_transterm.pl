#!/usr/local/bin/perl


BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

prepare_for_transterm.pl - creates coordinate files for the tranterm program

=head1 SYNOPSIS

USAGE: prepare_for_transterm.pl
    --input_list=/path/to/some/fsa.list
    --output_dir=/directory/for/output/files
    --database=aa1
    --schema=(legacy|chado)
    --pw_file=$EGC_SCRIPTS/password
  [ --log=/path/to/file.log
    --debug=4
    --help ]

=head1 OPTIONS

B<--input_list,-l>
    A list of fasta files to grab model coordinates from

B<--output_dir,-o>
    Directory to place the output .coords files in.
    
B<--database,-d>
    The database in which to look for gene annotation data.
    
B<--schema,-s>
    The schema of above mentioned database (Either legacy or chado)

B<--pw_file,-p>
    File containing user and password for database

B<--log,-l>
    Creates a log file.
    
B<--debug,-d>
    Large number, more verbose.

B<--help,-h>
    Prints out this message.

=head1 Description
    
    This script will take in a list of fasta files, query a database for annotations and use these to create a 
    coordinates file for the program transterm.  

    The input files are assumed to be assembly files following the naming convention
    
       /path/database.assembly.num.fsa
       
       ex.
       /usr/local/annotation/OSA1/asmbls/osa1.assembly.11971.fsa

    For the legacy schema, the models will be only those from the numbered assembly/ies.  For chado, the base
    name is used to query the feature table (i.e. given the above example, osa1.assembly.11971) and all related
    genes are used.

=head1 Input
    
    List of fsa file names.

=head1 Output

    Coordinate file in this format

    gene-name start stop sequenceId

    Where the sequenceId is the fasta header relating to the fasta file (in this case the 'osa1.assembly.11971' part).

=head1 CONTACT

    Kevin Galens (kgalens@tigr.org)

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Workflow::Logger;
use DBI;
use Data::Dumper;

###########CONSTANTS AND GLOBALS###################
use constant SCHEMA => { legacy => 'legacy', chado => 'chado' };
my $sql;
my @inputFiles;
my $output_dir;
my $schema;
my $asmbls = [];
my $coords;
my $dbh;
my $models = [];
my ($database, $user, $pass);

###################################################


my %options = ();
my $results = GetOptions (\%options, 
                          'input_list|i=s',
                          'input_file|f=s',
                          'output_dir|o=s',
                          'schema|s=s',
                          'database|d=s',
                          'pw_file|p=s',
                          'log|l=s',
                          'debug|d=s',
                          'help|h') || &_pod;

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

# display documentation
&_pod if( $options{'help'} );

# make sure all passed options are peachy
&check_parameters(\%options);

#Prepare the database
$dbh = &prepareDB($database, $user, $pass);

#Prepare the sql statement
$sql = &prepareSql($schema);

my $seqIds = &getSeqsFromFileNames(\@inputFiles);

$coords = &getCoords($sql, $schema, $dbh);

$coords = &filterCoords($schema, $seqIds, $coords);

&printCoordFiles( $coords );

#Count genes
my $count = 0;
foreach my $a (keys %{$coords}) {
    $count += scalar (keys %{$coords->{$a}});
}

print "$count records retrieved\n";

$dbh->disconnect();



#######################SUB ROUTINES####################################################

#Connect to the database
sub prepareDB {
    my ($db, $user, $pass) = @_;
    my $retval = DBI->connect("dbi:Sybase:server=SYBTIGR; packetSize=8092",$user, $pass) or
        &_die("Cannot connect to database server");
    $retval->do("use $db");
    print "Just connected to $db\n";
    return $retval;
}

#Prepare the sql statement
sub prepareSql {
    my $schema = shift;
    my $retSql = "";

    if($schema eq 'legacy') {
        
        $retSql = "SELECT a.feat_name, a.end5, a.end3, a.asmbl_id FROM asm_feature a, phys_ev p ".
            "WHERE a.feat_name = p.feat_name  AND a.feat_type = 'model' AND p.ev_type = 'working'";
      
    } elsif($schema eq 'chado') {
        
        $retSql = "SELECT f.name, l.fmin, l.fmax, l.strand, a.name from feature f, featureloc l, ".
            "cvterm c, feature a, cvterm ac WHERE f.type_id=c.cvterm_id AND c.name='gene' ".
            "AND f.is_analysis=0 AND f.is_obsolete=0 AND f.feature_id=l.feature_id AND ".
            "l.srcfeature_id=a.feature_id AND a.type_id=ac.cvterm_id";
        
    } else {
        &_die("Should have gone in legacy or chado\n");
    }

    return $retSql;

}

#From the fasta input list, parse out the basenames from the filenames.
sub getSeqsFromFileNames {
    my $files = shift;
    my @retval;
    
    foreach my $file(@{$files}) {
        chomp $file;

        if($file =~ m|.*/(.*)\.fsa$|) {
            push(@retval, $1);
        } else {
            &_die("Unable to parse the file name $file");
        }
    }
    
    return \@retval;

}

#Query for the coordinates of the models, retain everything.  Takes
# an sql statement, a schema (legacy or chado) and an open database connection handle.
sub getCoords {
    my ($sql, $schema, $dbh) = @_;
    my ($geneName, $start, $stop, $strand, $assembly);
    my $retval;
    my $count = 0;

    my $sth = $dbh->prepare($sql) or &_die("Bad prepare");
    $sth->execute();

    if($schema eq SCHEMA->{'legacy'}) {
        $sth->bind_columns(\$geneName, \$start, \$stop, \$assembly);
    } elsif($schema eq SCHEMA->{'chado'}) {
        $sth->bind_columns(\$geneName, \$start, \$stop, \$strand, \$assembly) or &_die("Bad bind");
    }

    while($sth->fetch()) {
        $count++;
        if($schema eq SCHEMA->{'chado'}) {
            $retval->{$assembly}->{$geneName}->{'strand'} = $strand;
            if($strand == -1) {
                my $tmp = $start;
                $start = $stop;
                $stop = $tmp;
            }
        }
        $retval->{$assembly}->{$geneName}->{'start'} = $start;
        $retval->{$assembly}->{$geneName}->{'stop'} = $stop;
       
    }
    print "Retrieved $count models\n";
    return $retval;
    
    
}

#Filter the coordinates for only the needed sequences.
sub filterCoords {
    my ($schema, $seqIds, $coords) = @_;
    my $retval;

    foreach my $seqId (@{$seqIds}) {
        $seqId = $1 if($seqId =~ /assembly\.(\d+)/ && ($schema eq "legacy"));
        $retval->{$seqId} = $coords->{$seqId};
    }

    return $retval;

}

#Print the coordinates to files.  Each input sequence (fasta) gets its own
#coords file, written to the $output_dir
sub printCoordFiles {
    my $crds = shift;
    
    foreach my $asmbl (keys %{$crds}) {
        my $asmblId = $asmbl;
        $asmblId = "$database.assembly.$asmblId" if($schema eq "legacy");
        print "Making file $asmblId.coords\n";
        open(OUT, "> $output_dir/$asmblId.coords") or 
            &_die("Unable to open $output_dir/$asmblId.coords for writing ($!)");

        foreach my $geneName(keys %{$crds->{$asmbl}}) {
            print OUT $geneName." ".$crds->{$asmbl}->{$geneName}->{'start'}.
                " ".$crds->{$asmbl}->{$geneName}->{'stop'}." ".$asmblId."\n";
        }

        close(OUT);

    }
   
}

#Check the parameters.
sub check_parameters {
    my $options = shift;

    if($options{'input_list'}) {
        &_die("$options{input_list} does not exist") unless(-e $options{'input_list'});
        open(IN, "< $options{'input_list'}") || &_die("Unable to open $options{'input_list'}");
        @inputFiles=<IN>;
        close(IN);
    } elsif($options{'input_file'}) {
        &_die("$options{input_file} does not exist") unless(-e $options{'input_file'});
        push(@inputFiles, $options{'input_file'});
    } else {       
        &_die("Option input_list or input_file is required");
    }
    

    unless(exists($options{'output_dir'})) {
        &_die("option output dir is required");
    } else {
        system("mkdir ".$options{output_dir}) unless(-d $options{'output_dir'});
    }
    $output_dir = $options{'output_dir'};
   
    my @schemas = keys %{&SCHEMA};
    &_die("schema option must be one of the following: @schemas") unless(exists(SCHEMA->{$options->{'schema'}}));
    $schema = $options->{'schema'};

    if($options->{'pw_file'}) {
        &_die("$options->{pw_file} does not exist") unless(-e $options->{'pw_file'});
    } else {
        &_die("option pw_file is required");
    }
    open(UP, "< $options->{pw_file}") or &_die("Cannot open $options->{pw_file} ($!)");
    chomp(($user, $pass) = <UP>);
    close(UP);

    unless($options->{'database'}) {
        &_die("option database is required");
    }
    $database = $options->{'database'};

}

#DIE!!
sub _die {
    my $msg = shift;
    $logger->logdie($msg);
}

#Podding.  Not just for apple anymore.
sub _pod {
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} ); 
}
