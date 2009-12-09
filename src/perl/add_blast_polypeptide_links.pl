#!/usr/bin/perl

=head1 NAME

add_blast_polypeptide_links.pl - Adds a BLAST analysis and associated analysisprops from a BSML file
to a chado database.  Does NOT load match features, only the analysis and attributes.

=head1 SYNOPSIS

USAGE: add_blast_polypeptide_links.pl 
            --bsml_list=/path/to/some_file.list
            --database=entamoeba
            --user=someuser
            --password=somepass
            --analysis_name=ber_analysis
          [ --server=SYBTIGR
            --skip_missing=0
            --program=ber
            --algorithm=ber
            --sourcename=/some/path
            --sequence_class=polypeptide
            --log=/path/to/some.log
          ]

=head1 OPTIONS

B<--bsml_list,-b>
    Input BSML file list containing <Analysis> elements and input_of Sequence features.

B<--database,-d>
    Sybase database name to connect to.

B<--user,-u>
    User account with select, insert, and update privileges on the specified database.

B<--password,-p>
    Password for user account specified.

B<--analysis_name,-a>
    You need to pass the name of the analysis you want searched within the passed BSML.
    This corresponds to the BSML Analysis.id value and is usually something like
    'wu-blastp_analysis'

B<--sequence_class,-e>
    Optional.  Can pass this to limit which BSML sequence element classes are considered.
    Entries in analysisfeature will only be created for these.

B<--sourcename,-o>
    Optional.  For each analysis load the chado analysis.sourcename column needs to be
    loaded with a value.  By default, this will be pulled from the Attribute element
    with a 'name' attribute value of 'sourcename', but you can override this with a
    custom one here.  The sourcename analysisprop will still be loaded.

B<--program,-r>
    Optional.  For each analysis load the chado analysis.program column needs to be
    loaded with a value.  By default, this will be pulled from the Attribute element
    with a 'name' attribute value of 'program', but you can override this with a
    custom one here.  The sourcename analysisprop (if present) will still be loaded.

B<--algorithm,-g>
    Optional.  For each analysis load the chado analysis.algorithm column needs to be
    loaded with a value.  By default, this will be pulled from the Attribute element
    with a 'name' attribute value of 'algorithm', but you can override this with a
    custom one here.  The sourcename analysisprop (if present) will still be loaded.

B<--skip_missing,-k>
    Optional.  By default, the script will fail if any input features are found within 
    the BSML that are not in the database.  Setting this skips these missing features
    and just puts a warning in the logs.  (default = 0)

B<--host,-s>
    Optional.  MySQL server to connect to (default = localhost).

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

We usually use the Ergatis bsml2chado component to load analyses from BSML into a 
chado database instance.  For BLAST results, this involves the creation of feature
rows for qry/sbj/match elements, as well as featureloc coordinate information and
all the rows needed for scoring.  For some tools (like Manatee), most of this isn't
necessary.  They only need the BLAST analysis element loaded and links from each
of the polypeptides as 'input_of' the analysis.  This script does that.

=head1  INPUT

Input can be a list file of BSML that is output from the Ergatis ber component.
Multiple list files can be concatenated together - the script handles simultaneous
input from multiple, separate analyses as long as results aren't mixed within any
single file.

=head1  OUTPUT

This script only produces file output if the --log option is used.  In the database, it 
performs inserts on the analysis, analysisprop and analysisfeature tables.  All row 
identifiers are written to the log.

=head1  CONTACT

    Joshua Orvis
    jorvis@users.sf.net

=cut

use strict;
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use XML::Twig;

my %options = ();
my $results = GetOptions (\%options, 
                          'bsml_list|b=s',
                          'database|d=s',
                          'host|s=s',
                          'user|u=s',
                          'password|p=s',
                          'log|l=s',
                          'program|r=s',
                          'skip_missing|k=i',
                          'analysis_name|a=s',
                          'sourcename|o=s',
                          'algorithm|g=s',
                          'sequence_class|e=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">$options{log}") || die "can't create log file: $!";
}

# Connect to the db
_log("attempting to create database connection");
my $dbh = DBI->connect("DBI:mysql:database=$options{database};host=$options{host};", 
                        $options{user}, $options{password}, 
                        {PrintError=>1, RaiseError=>1});
$dbh->do("use $options{database}");

my %cvterms;

## pull any needed cvterms
$cvterms{input_of} = get_cvtermid('input_of') || die "failed to pull cvterm_id for term 'input_of'";

## structure like:
#   $h{pipelineid.token} = $analysis_id
my %analyses_inserted;

## get the next available IDs
my %next_id;
$next_id{analysis} = $dbh->selectrow_array("SELECT max( analysis_id ) FROM analysis") + 1;
$next_id{analysisprop} = $dbh->selectrow_array("SELECT max( analysisprop_id ) FROM analysisprop") + 1;
$next_id{analysisfeature} = $dbh->selectrow_array("SELECT max( analysisfeature_id ) FROM analysisfeature") + 1;

_log("INFO: first available analysis_id seems to be $next_id{analysis}");
_log("INFO: first available analysisprop_id seems to be $next_id{analysisprop}");

my $qry = qq{
    INSERT INTO analysis ( analysis_id, name, description, program, programversion, algorithm, sourcename, timeexecuted )
         VALUES ( ?, ?, 'none', ?, ?, ?, ?, 'NOW()' )
};
my $analysis_inserter = $dbh->prepare($qry);

$qry = qq{
    INSERT INTO analysisprop ( analysisprop_id, analysis_id, type_id, value )
         VALUES ( ?, ?, ?, ? )
};
my $analysisprop_inserter = $dbh->prepare($qry);

$qry = qq{
    INSERT INTO analysisfeature ( analysisfeature_id, feature_id, analysis_id, type_id )
         VALUES ( ?, ?, ?, ? );
};
my $analysisfeature_inserter = $dbh->prepare($qry);

$qry = qq{
    SELECT analysisfeature_id FROM analysisfeature WHERE analysis_id = ? AND feature_id = ?;
};
my $analysisfeature_checker = $dbh->prepare( $qry );

## read the list file
open (my $ifh, "<$options{bsml_list}") || die "can't read input BSML list: $!";

my $current_analysis_id;

while ( my $line = <$ifh> ) {
    chomp $line;
    next if $line =~ m|^\s*$|;
    
    parse_input_bsml($line);
}


$analysisfeature_inserter->finish();
$analysisprop_inserter->finish();
$analysis_inserter->finish();

$dbh->disconnect;

exit(0);

sub parse_input_bsml {
    my $file = shift;

    
    ## keeps rows to insert - which should happen after each file
    ## structure like: 
    #   $h->{$table} = [ [col1, col2, ...], ... ]
    my $rows = {};
    
    my $twig = XML::Twig->new(
        twig_roots => {
            Analysis => sub { parse_analysis(@_, $rows) },
            'Sequence' => sub { parse_sequence(@_, $rows) },
        },
    );

    _log("INFO: twig parsing file: $file");
    
    ## need to handle whether compressed or not
    my $input_fh;
    if ( ! -e $file && -e "$file.gz" ) {
        $file .= '.gz';
    }
    
    if ( $file =~ /\.gz$/ ) {
        open($input_fh, "<:gzip", $file) || die "can't read compressed file $file: $!";
    } else {
        open($input_fh, "<$file") || die "can't read file $file: $!";
    }
    
    $twig->parse($input_fh);
    
    ## now do the inserts
    for my $row ( @{$$rows{analysis}} ) {
        print "Analysis:\n";
        $" = "|";
        print "\t@$row\n";
        _log("INFO: inserting analysis ID $$row[0]");
        $analysis_inserter->execute( @$row ) or die "failed to insert analysis row: @$row";    
    }
    
    for my $row ( @{$$rows{analysisprop}} ) {
        #print "analysisprop:\n";
        #print "\t@$row\n";
        _log("INFO: inserting cvterm_id $$row[2] analysisprop for analysis $$row[1]");
        $analysisprop_inserter->execute( @$row ) or die "failed to insert analysisprop row: @$row";
    }
    
    for my $row ( @{$$rows{analysisfeature}} ) {
        #print "analysisfeature:\n";
        #print "\t@$row\n";
        _log("INFO: inserting analysisfeature_id $$row[0] for feature $$row[1] as input of analysis_id: $current_analysis_id");
        if( analysisfeature_doesnt_exist( $$row[1], $current_analysis_id ) ) {
            $analysisfeature_inserter->execute( $$row[0], $$row[1], $current_analysis_id, $$row[3] ) or die "failed to insert analysisfeature row: @$row";
        }
    }
    
    ## reset the analysis ID
    $current_analysis_id = undef;
}

sub analysisfeature_doesnt_exist {
    my ($feature_id, $analysis_id) = @_;
    my $retval = 1;
    $analysisfeature_checker->execute( $analysis_id, $feature_id );
    my $results = $analysisfeature_checker->fetchall_arrayref();
    if( @{$results} > 0 ) {
        $retval = 0;
    }
    return $retval;
}

sub parse_sequence {
    my ($t, $elt, $rows) = @_;
    
   
    ## see if this sequence has a Link element with role 'input_of'
    my $sequence_id = undef;
    for my $link ( $elt->children('Link') ) {
        if ( $link->{att}->{role} eq 'input_of' && $link->{att}->{href} eq "\#$options{analysis_name}") {
            $sequence_id = $elt->{att}->{id};
        }
    }
    
    return unless defined $sequence_id;

    ## skip this if it isn't of the right class (if we're even checking)
    if ( $options{sequence_class} && $elt->{att}->{class} ne $options{sequence_class} ) {
        _log( "INFO: skipping sequence $sequence_id because it's not of the right class" );
        return;
    }

    _log("INFO: processing input sequence with id $sequence_id");
    
    my $feature_id = get_featureid($sequence_id);

    if ( $feature_id ) {
        _log("INFO: sequence ID $sequence_id corresponds to feature.id $feature_id");
    
        ## undef here in col 3 is replaced by current analysis id on insert
        push @{$$rows{analysisfeature}}, [ $next_id{analysisfeature}++, $feature_id, 
                                           undef, $cvterms{input_of} ];
    } else {
    
        ## ok, the feature_id wasn't found.  depending on the option passed
        #  we'll either die here or just print a warning
        if ( $options{skip_missing} ) {
            _log("WARN: skipping feature $sequence_id because a feature_id for it could not be found");
            
        } else {
            die("feature_id not found for feature.uniquename $sequence_id");
        }
    }
}

sub parse_analysis {
    my ($t, $elt, $rows) = @_;
    
    my %atts;

    ## we only care about requested analyses
    if ( $elt->{att}->{id} eq $options{analysis_name} ) {

        ## get all the Attributes
        for my $att ( $elt->children('Attribute') ) {
            $atts{$att->{att}->{name}} = $att->{att}->{content};
        }
    } else {
        _log("INFO: skipping Analysis element found with ID: $elt->{att}->{id}");
    }

    ## the following attributes need to have been found
    foreach my $att ( qw( component_name version pipelineid pipeline_token ) ) {
        if (! exists $atts{$att} ) {
            die "Attribute $att not found in $elt->{att}->{id} analysis";
        }
    }

    my $pipeline_id_token = "$atts{pipelineid}.$atts{pipeline_token}";

    ## skip this if we've already found this analysis
    if ( exists $analyses_inserted{$pipeline_id_token} ) {
        $current_analysis_id = $analyses_inserted{$pipeline_id_token};
        return;
    }

    ## sourcename can either be pulled from the BSML or overridden by the user.
    my $sourcename = $options{sourcename} || $atts{sourcename};
    my $program    = $options{program} || $atts{program};
    my $algorithm    = $options{algorithm} || $atts{algorithm};

    push @{$$rows{analysis}}, [ $next_id{analysis}, $atts{component_name}, $program, 
                                $atts{version}, $algorithm, $sourcename ];
    print "Pushed! $atts{component_name}, $program, $atts{version}, $algorithm, $sourcename\n";

    _log("INFO: setting current analysis id to: $next_id{analysis}");
    $current_analysis_id = $next_id{analysis};
    $analyses_inserted{$pipeline_id_token} = $next_id{analysis}++;

    ## now we need to insert the analysisprops
    for my $att ( keys %atts ) {
        my $cvterm_id = get_cvtermid( $att );

        if ( defined $cvterm_id ) {
            push @{$$rows{analysisprop}}, [ $next_id{analysisprop}++, $analyses_inserted{$pipeline_id_token},
                                            $cvterm_id, $atts{$att} ];

        } else {
            _log("WARN: skipping analysisprop '$att' ($atts{$att}) because I couldn't find a cvterm for it");
        }
    }
}

sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
}

sub get_cvtermid {
    my $name = shift;
    my $qry = qq{
        SELECT cvterm_id
        FROM cvterm
        WHERE name = ?
    };
    
    my $cvterm_selector = $dbh->prepare($qry);
       $cvterm_selector->execute($name);
    
    my $cvterm_id = ( $cvterm_selector->fetchrow_array )[0];
    $cvterm_selector->finish();

    _log("INFO: got cvterm_id $cvterm_id for name $name");
    return $cvterm_id;
    
}

sub get_featureid {
    my $uniquename = shift;
    
    my $qry = qq{
        SELECT feature_id
          FROM feature
         WHERE uniquename = ?
    };
    my $feature_id_selector = $dbh->prepare( $qry );
       $feature_id_selector->execute( $uniquename );
    
    my $feature_id = ( $feature_id_selector->fetchrow_array )[0];
    $feature_id_selector->finish;
    
    _log("INFO: got feature_id $feature_id for uniquename $uniquename");
    return $feature_id;
}

sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( bsml_list database user password analysis_name );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
    
    ##
    ## you can do other things here, such as checking that files exist, etc.
    ##
    
    ## handle some defaults
    $options{host} = 'localhost' unless ($options{host});
    $options{skip_missing} = 0 unless ($options{skip_missing});
    $options{sourcename} = undef unless ($options{sourcename});
    $options{program} = undef unless ($options{program});
    $options{sequence_class} = undef unless ($options{sequence_class});
}
