#!/usr/local/bin/perl -w

=head1  NAME 

repeat2legacydb.pl - load repeat data into the legacy database schema.

=head1 SYNOPSIS

USAGE: repeat2legacydb.pl 
        --input_list=/path/to/somefile.bsml.list
        --database=sma1
        --repeat_file=/path/to/some/repeat.db
      [ --log=/path/to/some.log 
        --delete_existing=1
      ]

=head1 OPTIONS

B<--input_list,-i> 
    BSML list file from an repeat-generating workflow component run.  Currently tested
    components are:  repeatmasker, trf.

B<--repeat_file,-r> 
    The full path to the repeat library used.  This will be populated in ORF_attribute.score3  Note:  If loading trf results, enter 'trf' for this parameter.

B<--database,-d> 
    Sybase project database ID.

B<--delete_existing,-x> 
    optional.  will first delete any existing repeats for this repeat type for each
    assembly passed.

B<--log,-d> 
    optional.  will create a log file with summaries of all actions performed.

B<--help,-h> 
    This help message/documentation.

=head1   DESCRIPTION

This script is used to load repeat evidence into the legacy database schema.
It currently supports input from the repeatmasker and trf components as 
specified below.

=head1 INPUT

The file names within the BSML list must be contain the portion 
'assembly.$assemblyid.$program.bsml', like:

    /some/path/sma1.assembly.29106.repeatmasker.bsml
    
So that the assembly number (29106) and program can be extracted from the name.

=head1 OUTPUT

This is a database loading script.  There is no other output unless you use 
the --log option.

=head1 CONTACT

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use XML::Twig;

my %options = ();
my $results = GetOptions (\%options, 
              'input_list|i=s',
              'database|d=s',
              'repeat_file|r=s',
              'delete_existing|x=i',
              'log|l=s',
              'help|h') || pod2usage();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

## make sure all passed options are peachy
&check_parameters(\%options);

## get the current user
my ($user) = getpwuid($<);

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">$options{log}") || die "can't create log file: $!";
}


my %progs_handled = ( trf => 1,
                      repeatmasker => 1
                    );

## connect to the database
my $dbh = DBI->connect("dbi:Sybase:server=SYBTIGR; packetSize=8092", 'egc', 'egcpwd');
$dbh->do("use $options{database}");

## load the file list
my @files;
open(my $listfh, "<$options{input_list}") || die "can't open file list: $!";
    while (<$listfh>) {
        chomp;
        next if (/^\s*$/);
        _log("adding $_ to the process list");
        push @files, $_;
    }
close $listfh;

## prepare a statement for inserting into asm_feature
my $qry = "INSERT INTO asm_feature (feat_type, feat_method, end5, end3, assignby, date, " .
                                   "feat_name, asmbl_id, change_log, save_history, curated) " .
          "VALUES ( 'repeat', ?, ?, ?, '$user', GETDATE(), ?, ?, 1, 0, 0 )";
my $asm_feature_insert = $dbh->prepare($qry);

## prepare a statement for deleting from asm_feature
$qry = "DELETE FROM asm_feature " .
          "WHERE feat_type = 'repeat' " .
          "   AND feat_method = ? " .
          "   AND asmbl_id = ? ";
my $asm_feature_delete = $dbh->prepare($qry);

## prepare a statement for inserting into ORF_attribute
$qry = "INSERT INTO ORF_attribute (feat_name, att_type, curated, method, date, assignby, " .
                           "score, score_desc, score2, score2_desc, score3, score3_desc) " .
       "VALUES ( ?, 'repeat', 0, ?, GETDATE(), '$user', ?, ?, ?, ?, ?, ? )";
my $orf_attribute_insert = $dbh->prepare($qry);

## score3 needs to have name of fasta file searched

## prepare a statement for deleting from ORF_attribute
$qry = "DELETE FROM ORF_attribute " .
          "WHERE feat_name LIKE ? " .
          "  AND att_type = 'repeat' " .
          "  AND method = ?";
my $orf_attribute_delete = $dbh->prepare($qry);

foreach my $file (@files) {
    _log("processing $file");
    
    my ($asmbl_id, $prog_name);
    if ($file =~ /\.assembly\.(\d+)\.(.+)\.bsml/) {
        ($asmbl_id, $prog_name) = ($1, $2);
        _log("extracted asmbl_id $asmbl_id, prog_name $prog_name from $file");
    } else {
        die "unable to parse asmbl_id from $file\n";
    }
    
    ## have we handled this prog?
    if (! $progs_handled{$prog_name}) {
        die "haven't yet handled output from $prog_name.  sorry.\n";
    }
    
    ## are we deleting?
    if ( $options{delete_existing} ) {
        _log("deleting $prog_name entries from asm_feature where asmbl_id = $asmbl_id");
        $asm_feature_delete->execute( $prog_name, $asmbl_id );
        
        _log("deleting $prog_name entries from ORF_attribute where feat_name like $asmbl_id.repeat%");
        $orf_attribute_delete->execute( "$asmbl_id.repeat%", $prog_name );
    }
    
    ## get the latest database ID
    my $nextid = &get_latest_id($asmbl_id, 'repeat') + 1;
    
    my $twig = XML::Twig->new(
                                twig_roots => { 
                                                'Feature' => sub {
                                                                my ($twig, $elt) = @_;
                                                                &process_repeat_feature($twig, $elt, \$asmbl_id, \$prog_name, \$nextid);
                                                             },
                                              },
                             );
    
    my $ifh;
    if (-e "$file.gz") {
        $file .= ".gz";
    } elsif (-e "$file.gzip") {
        $file .= ".gzip";
    }
    
    if ($file =~ /\.(gz|gzip)$/) {
        open ($ifh, "<:gzip", $file);
    } else {
        open ($ifh, "<$file");
    }
    $twig->parse($ifh);
    close $ifh;
}

## clean up
$asm_feature_insert->finish();
$asm_feature_delete->finish();
$orf_attribute_insert->finish();
$orf_attribute_delete->finish();
$dbh->disconnect();

exit;

sub get_latest_id {
    my ($asm_id, $feat_type) = @_;

    my $stmt = $dbh->prepare("SELECT DISTINCT feat_name " .
                             "FROM asm_feature " .
                             "WHERE asmbl_id = ? " .
                             "AND feat_type = ? ");
    $stmt->execute($asm_id, $feat_type);
    my $last_id = "000000";
    while (my $row = $stmt->fetchrow_arrayref) {
        my $id = $& if $$row[0] =~ /\d+$/;
        $last_id = $id if $id > $last_id;
    }

    return $last_id;
}

sub process_repeat_feature {
    my ($twig, $elt, $asmbl_id, $prog_name, $nextid) = @_;
    my ($score, $score_desc, $score2, $score2_desc);
    
    my $feat_id = $elt->{att}->{id} || 'unknown';
    
    ## skip this feature if it isn't a repeat region
    my $class = $elt->{att}->{class};
    if ($class ne 'repeat_region') {
        _log("skipping feature $feat_id because its class is $class");
        return 1;
    }
    
    ## it better have an interval loc
    unless ( $elt->has_child('Interval-loc') ) {
        die "found a repeat_region without an Interval-loc!  panic! ($feat_id)\n";
    }
    
    my $interval_loc = $elt->first_child('Interval-loc');
    my $startpos     = $interval_loc->{att}->{startpos} + 1;
    my $endpos       = $interval_loc->{att}->{endpos};
    
    ## is the repeat on the reverse strand?
    if ( $interval_loc->{att}->{complement} ) {
        ($startpos, $endpos) = ($endpos, $startpos);
    }

    my $feat_name = "$$asmbl_id.repeat" . sprintf("%06d", $$nextid++);
    
    if ($$prog_name eq 'trf') {
        $score       = &get_attribute($elt, 'consensus_text');
        $score_desc  = 'consensus text';   
        $score2      = &get_attribute($elt, 'raw_score');
        $score2_desc = 'score';
    } elsif ($$prog_name eq 'repeatmasker') {
        $score       = &get_attribute($elt, 'matching_repeat');
        $score_desc  = 'class';
        $score2      = &get_attribute($elt, 'sw_score');
        $score2_desc = 'score';
    } else {
        die "unhandled prog name ($$prog_name)\n";
    }
    
    _log("asm_feature_insert->execute($$prog_name, $startpos, $endpos, $feat_name, $$asmbl_id)");
    $asm_feature_insert->execute($$prog_name, $startpos, $endpos, $feat_name, $$asmbl_id);
    
    _log("orf_attribute_insert->execute( $feat_name, $$prog_name, $score, $score_desc, $score2, $score2_desc, $options{repeat_file}, 'repeat database used' )");
    $orf_attribute_insert->execute( $feat_name, $$prog_name, $score, $score_desc, $score2, $score2_desc, $options{repeat_file}, 'repeat database used' );
    
    $twig->purge;
}

sub get_attribute {
    my ($elt, $name) = @_;
    my $att = undef;
    
    for my $attribute ( $elt->children('Attribute') ) {
        if ( $attribute->{att}->{name} eq $name ) {
            $att = $attribute->{att}->{content};
            last;
        }
    }
    
    if (defined $att) {
        return $att;
    } else {
        die "failed to extract $name from attributes\n";
    }
}


sub _log {
    my $msg = shift;
    
    print $logfh "$msg\n" if $logfh;
}

sub check_parameters {
    
    ## database and input_list are required
    unless ( defined $options{database} && $options{input_list} && $options{repeat_file} ) {
        print STDERR "database and input_list options are required\n\n";
        pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
    }
    
    ## make sure input list exists
    unless ( -e $options{input_list} ) {
        print STDERR "\n\ninput_list $options{input_list} does not exist\n\n";
        exit(1);
    }
    
    ## set some defaults
    $options{delete_existing} = 0 unless ( $options{delete_existing} );
}







