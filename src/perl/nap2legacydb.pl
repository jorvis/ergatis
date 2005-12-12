#!/usr/local/bin/perl

=head1  NAME 

nap2legacydb.pl - load nap BSML file(s) into the legacy database schema

=head1 SYNOPSIS

USAGE:  nap2legacydb.pl
            --bsml_file=/path/to/somemapfile.bsml
            --list_file=/path/to/somefile.list
            --database=aa1
          [
            --debug=4
            --log=/path/to/somefile.log
          ]

=head1 OPTIONS

B<--bsml_file,-b> 
    Path to a single BSML nap result file.

B<--list_file,-i> 
    Path to a list of nap BSML output files.

B<--database,-d> 
    Sybase database name for the project in which these nap predictions
    will be loaded.  This is case-sensitive, and most of our project names
    within Sybase are lower-case.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--help,-h> 
    This help message

=head1 DESCRIPTION

This script is used to load gene predictions from the nap program, in BSML format,
into one of the legacy databases.  Such BSML files are most likely generated from
the nap workflow component.

=head1 INPUT

The input may be a single bsml file, defined with --bsml_file, or an input list of
bsml files, defined with --input_list.  To properly extract the assembly IDs from
the file name, the name will need to contain the string ".assembly.N" where N is the 
assembly ID.  

=head1 OUTPUT

This file generates no output on the file system unless a log file is specified.  It
does perform inserts into the project database passed.

=head1 CONTACT

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
BEGIN {
    require '/usr/local/devel/ANNOTATION/cas/lib/site_perl/5.8.5/Workflow/Logger.pm';
    import Workflow::Logger;
}
use XML::Twig;
use DBI;

#######
## ubiquitous options parsing and logger creation
my %options = ();
my $results = GetOptions (\%options, 
                            'bsml_file|b=s',
                            'list_file|i=s',
                            'database|d=s',
                            'debug=s',
                            'log|l=s',
                            'help|h') || pod2usage();

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = Workflow::Logger::get_logger();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

&check_parameters(\%options);

## create the database connection
my $dbh = DBI->connect("dbi:Sybase:server=SYBTIGR; packetSize=8092", 'egc', 'egcpwd', {RaiseError => 1, PrintError => 1});
my $result = $dbh->do("use $options{database}");

## get the query ready
my $qry = "INSERT INTO evidence (feat_name, ev_type, accession, end5, end3, rel_end5, rel_end3, m_lend, m_rend, curated, date, assignby, change_log, save_history, method, per_id, per_sim, score, db, pvalue, domain_score, expect_domain, total_score, expect_whole, chainID) " .
          "VALUES (?, 'nap', ?, ?, ?, -1, -1, ?, ?, 0, getdate(), 'workflow', 0, 0, 'aat_aa', ?, ?, ?, ?, NULL, ?, NULL, ?, NULL, ?)";
my $insert = $dbh->prepare($qry);


## gather the files we are going to process
my @files;
if (defined $options{bsml_file}) {
    $logger->debug("adding $options{bsml_file} to input list") if ($logger->is_debug);
    push @files, $options{bsml_file};
} elsif (defined $options{list_file}) {
    open(my $listfh, "<$options{list_file}") || $logger->logdie("can't read list file $options{list_file} : $!");

    while (<$listfh>) {
        chomp;
        next if (/^\s*$/);
        $logger->debug("adding $_ to input list") if ($logger->is_debug);
        push @files, $_;
    }
}


## process each file
my ($file, $asmbl_id);
for (@files) {
    $file = $_;
    $logger->info("processing $file") if ($logger->is_info);
    
    ## we can parse the assembly ID out of the file name
    if ($file =~ /\.assembly\.(\d+)/) {
        $asmbl_id = $1;
    } else {
        $logger->logdie("couldn't parse asmbl_id from file name!");   
    }
    
    
    ## create the twig
    my $twig = XML::Twig->new(
                                twig_roots  => { 'Seq-pair-alignment' => \&processSeqPairAlignment }
                             );
    $twig->parsefile($file);
    
}

$insert->finish();
$dbh->disconnect();

exit(0);

sub check_parameters {
    my ($options) = @_;

    ## bsml_file or list_file must have been passed
    unless ( $options{bsml_file} || $options{list_file} ) {
        $logger->logdie("bsml_file or list_file must be passed");
    }

    ## check the bsml file, if passed
    if ( $options{bsml_file} && ! -e $options{bsml_file} ) {
        $logger->logdie("bsml_file was passed but does not exist");
    }

    ## check map file
    if ( $options{list_file} && ! -e $options{list_file} ) {
        $logger->logdie("list_file either not passed or does not exist");
    }

    ## make sure a database was passed
    unless ($options{database}) {
        $logger->logdie("database must be passed!");
    }

    if(0){
        pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});    
    }
}

sub get_seqpairrun_attribute {
    my ($spr, $att) = @_;
    my $val = $spr->{att}->{$att};
    if (! defined $val) { $logger->logdie("couldn't find a $att in the seq pair run") }
    return $val;
}


sub qry_db_name {
    my $name = shift;
    my $id;

    my $qry = "SELECT id FROM common..search_dbs WHERE name = '$name'";
    my $dsh = $dbh->prepare($qry);
       $dsh->execute();
    
    while (my $res = $dsh->fetchrow_hashref) {
        $id = $res->{id};
    }
    
    $dsh->finish;
    
    return $id || 0;
}

sub get_db_id {
    my $name = shift;
    my $id;
    
    ## first see if one already exists
    $id = qry_db_name($name);
    
    ## if we didn't get an ID from that, insert this name
    if ($id) {
        $logger->debug("got existing search_dbs id $id for $name") if ($logger->is_debug);
    } else {
        ## we do it this way because the last_insert_id method is unreliable in sybase
        $id = insert_db_name($name);
        $logger->debug("got new search_dbs id $id for $name") if ($logger->is_debug);
    }
    
    ## make sure we got an ID
    if (! $id) { $logger->logdie("failed to get a search_dbs ID for database $name") }
    
    return $id;
}

sub insert_db_name {
    my $name = shift;
    
    my $qry = "INSERT INTO common..search_dbs (name, release_notes, date, assignby, iscurrent, type) " .
              "VALUES ('$name', NULL, getdate(), 'workflow', 1, 'custom')";
    $dbh->do($qry);
    
    return qry_db_name($name);
}

sub processSeqPairAlignment { 
    my ($twig, $feat) = @_;
    
    ## get the total score, which is kept in an Attribute element
    my $total_score;
    for my $attribute ( $feat->children('Attribute') ) {
        if ($attribute->{att}->{'name'} eq 'total_score') {
            $total_score = $attribute->{att}->{'content'};
        }
    }
    
    ## make sure we got the total score
    if (! defined $total_score) {
        $logger->logdie("failed to get total_score from Seq-pair-alignment");
    }
    
    ## get the subject match id
    my $sbj_match_id = $feat->{att}->{compseq} || $logger->logdie("couldn't find a compseq in the seq pair alignment");
    
    ## we need to get the name of the file searched, which is contained within the compxref attribute
    ##  this is usually in the format:
    ##  /usr/local/scratch/annotation/NFA1/search_databases/AFU2_Working_Models_20050414_nt.fasta:afu2_11.m00684_11.t00167
    ##  we need to take off the entry after the ":" and then pull just the file name
    my $db_searched = $feat->{att}->{compxref} || $logger->logdie("couldn't find a compxref in the seq pair alignment");
    if ($db_searched =~ m|.+/(.+):|) {
        $db_searched = $1;
    } else {
        $logger->logdie("compxref ($db_searched) in unrecognized format");
    }    
    
    ## now we need to get a sybase ID for this database.
    my $db_id = get_db_id($db_searched);
    
    my $parent_chain = 0;
    
    for my $seqpairrun ( $feat->children('Seq-pair-run') ) {
        my ($p_sim, $p_ident, $this_chain, $score);
        
        ## get the score
        $score = get_seqpairrun_attribute($seqpairrun, 'runscore');
        
        ## get the positions
        my ($qry_start, $qry_stop, $qry_comp, $qry_run_length, $sbj_start, $sbj_stop, $sbj_comp, $sbj_run_length);
        $qry_start      = get_seqpairrun_attribute($seqpairrun, 'refpos');
        $qry_start++;   ## numbering from one adjustment
        $qry_comp       = get_seqpairrun_attribute($seqpairrun, 'refcomplement');
        $qry_run_length = get_seqpairrun_attribute($seqpairrun, 'runlength');
        $qry_stop       = $qry_start + $qry_run_length;
        
        $sbj_start      = get_seqpairrun_attribute($seqpairrun, 'comppos');
        $sbj_start++;   ## numbering from one adjustment
        $sbj_comp       = get_seqpairrun_attribute($seqpairrun, 'compcomplement');
        $sbj_run_length = get_seqpairrun_attribute($seqpairrun, 'comprunlength');
        $sbj_stop       = $sbj_start + $sbj_run_length;

        ## legacy db stores the coordinates flipped if on complement
        if ($qry_comp == 1) {
            ($qry_start, $qry_stop) = ($qry_stop, $qry_start);
        } elsif ($qry_comp != 0) {
            $logger->logdie("unknown value ($qry_comp) for refcomplement (expected 1 or 0)");
        }
        
        if ($sbj_comp == 1) {
            ($sbj_start, $sbj_stop) = ($sbj_stop, $sbj_start);
        } elsif ($sbj_comp != 0) {
            $logger->logdie("unknown value ($sbj_comp) for compcomplement (expected 1 or 0)");
        }
        
        ## go through the Attribute elements to grab needed values
        for my $attribute ( $seqpairrun->children('Attribute') ) {
            if ($attribute->{att}->{'name'} eq 'chain_number') {
                $this_chain = $attribute->{att}->{'content'};
                
                ## if the parent hadn't been defined yet, this is it
                if (! $parent_chain) {
                    $parent_chain = $this_chain;
                }
                
            } elsif ($attribute->{att}->{'name'} eq 'percent_identity') {
                $p_ident = $attribute->{att}->{'content'};
            
            } elsif ($attribute->{att}->{'name'} eq 'percent_similarity') {
                $p_sim = $attribute->{att}->{'content'};
            }
        }
        
        ## make sure we got everything
        unless ( defined $p_sim && defined $p_ident && defined $this_chain ) {
            $logger->logdie("failed to get one of p_sim, p_ident, or chain from Seq-pair-run");
        }
        
        ## do the insert
        $logger->debug("doing insertion with the following placeholders: (\"$asmbl_id.intergenic\", $sbj_match_id, $qry_start, $qry_stop, $sbj_start, $sbj_stop, $p_ident, $p_sim, $score, $db_id, $score, $total_score, $parent_chain)\;") if ($logger->is_debug);
        $insert->execute("$asmbl_id.intergenic", $sbj_match_id, $qry_start, $qry_stop, $sbj_start, $sbj_stop, $p_ident, $p_sim, $score, $db_id, $score, $total_score, $parent_chain);
    }
}

