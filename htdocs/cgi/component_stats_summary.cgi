#!/usr/bin/perl -w

use strict;
use warnings;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Date::Manip;
use Ergatis::Common;
use Ergatis::ConfigFile;
use Ergatis::Monitor;
use File::stat;
use HTML::Template;
use POSIX;
use XML::Twig;
use File::Basename;
use Statistics::Descriptive;

my $q = new CGI;

#print $q->header( -type => 'text/html' );
print $q->header( -type => 'text/plain' );



## read the ergatis config file
#my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );


my $pipeline_xml = $q->param("pipeline_xml") || die "pass pipeline_xml";
my $pipeline_state = $q->param("pipeline_state") || die "pass pipeline state";

# Print summary info or all the native data rows?
my $is_summary = $q->param("summary");
$is_summary = defined($is_summary) ? $is_summary : 0;

## it may have been compressed
if (! -e $pipeline_xml && -e "$pipeline_xml.gz") {
    $pipeline_xml .= '.gz';
}

my $pipeline_xml_fh;
if ($pipeline_xml =~ /\.gz/) {
    open($pipeline_xml_fh, "<:gzip", "$pipeline_xml") || die "can't read $pipeline_xml: $!";
} else {
    open($pipeline_xml_fh, "<$pipeline_xml") || die "can't read $pipeline_xml: $!";
}

# my $twig = XML::Twig->new( );
my $twig = new XML::Twig( twig_handlers => 
                {
#                 'command/name[string() =~  /create iterato[r].*groupings/ ]'    => \&process_name,     # full path
#                 'command/name'    => \&name_check,     # full path
                }
                       ); 



$twig->parse($pipeline_xml_fh);

my %fstats; # keyed off the file basename
my %gstats; # keyed off the group xml: e.g., g1.xml.gz
# linked in g13.iter.xml.gz: parentFileName=group xml

foreach my $command ( $twig->root->descendants( 'command' ) ) {
  my $name = $command->first_child( 'name' )->text;
  next unless $name =~ /create iterator (.*) groupings/;
  my $iter = $1;

  my $input_iter_list = '';  # e.g., 38_default/i1.list
  my $output_directory = ''; # e.g., 38_default/i1
  my $iterator_xml = '';     # e.g., 38_default/i1/i1.xml.gz

  foreach my $param ( $command->children( 'param' ) ) {
    if ( $param->first_child('key')->text eq '--input_iter_list' ) {
      $input_iter_list = $param->first_child('value')->text;
    }
    elsif ( $param->first_child('key')->text eq '--output_directory' ) {
      $output_directory = $param->first_child('value')->text;
    }
  }

  # get from commandSet sibling
  foreach my $sib ( $command->next_sibling( 'commandSet' ) ) {
    if ( $sib->first_child('name')->text eq "iterator $iter" ) {
      $iterator_xml = $sib->first_child('fileName')->text;
    }
  }

  die "No input_iter_list" if ( $input_iter_list eq '');
  die "No output_directory" if ( $output_directory eq '');
  die "No iterator_xml" if ( $iterator_xml eq '');

  #print "$iterator_xml\n";

  # Parse timing info from iterator_xml: enters qsub/pending, enters prolog, epilog start, epilog end
    if ( -e $iterator_xml ) {

        my $iterator_xml_fh;
        if ($iterator_xml =~ /\.gz/) {
            open($iterator_xml_fh, "<:gzip", "$iterator_xml") || die "can't read $iterator_xml: $!";
        } else {
            open($iterator_xml_fh, "<$iterator_xml") || die "can't read $iterator_xml: $!";
        }

        ## create the twig
        my $ixml = XML::Twig->new( twig_roots => {
                                     'commandSet[@type="remote-serial"]' => \&process_ixml_commandSet, 
                                   }
                                 );
        $ixml->parse($iterator_xml_fh);

    }
    else {
        die "Missing file $iterator_xml";
    }


  # parse input iter files for input file sizes and stats
  open (my $FIN, $input_iter_list) || die "Unable to open input_iter_list ($input_iter_list): $!";
  <$FIN>; # dump header row
  while (my $line = <$FIN>) {
    chomp($line);
    my @row = split( "\t", $line);

    my $base = $row[0];
    my $ifile = $row[2];
    
    $fstats{$base}{iter} = $iter;
    $fstats{$base}{size} = `ls -l $ifile | cut -d ' ' -f 5 `;
    $fstats{$base}{size} =~ s/\s*$//;

# punting on counting non-whitespace characters for now
#    open( my $DIN, $ifile ) || die "Unable to open ifile $ifile: $!";
#      while (my $dline = <$DIN>) {
#      
#      }
#    close($DIN);

    $fstats{$base}{seqs} = `grep '>' $ifile | wc -l`;
    $fstats{$base}{seqs} =~ s/\s*$//;
#    $fstats{$base}{chars} = `grep -v '>' $ifile | wc -m`;

    #print $q->h4( "$base $ifile ".$fstats{$base}{size}." ".$fstats{$base}{seqs} );
  }

  close($FIN);

  # parse timing info from 38_default/i1/ *.iter.xml.gz
  foreach my $filename ( `find $output_directory -name "*.iter.xml.gz"` ) {
    chomp($filename);
    #print $q->h3( $filename);

     if ( -e $filename ) {

        my $filename_fh;
        if ($filename =~ /\.gz/) {
            open($filename_fh, "<:gzip", "$filename") || die "can't read $filename: $!";
        } else {
            open($filename_fh, "<$filename") || die "can't read $filename: $!";
        }

        ## create the twig
        my $twig2 = XML::Twig->new( twig_handlers => {
                                        'commandSet' => \&process_subflowgroup,
                                   }
                                 );
        $twig2->parse($filename_fh);

    }
    else { 
        die "Missing file $filename";
    }
  }

  # Calculate timings
  use Data::Dumper;

  #output timings
  my @times = qw( qsubStart prologStart execStart execEnd epilogStart epilogEnd );
  my $done_header = 0;
  my %negs;
  foreach my $key (sort {$a =~ /_(\d+)$/; my $x=$1; $b =~ /_(\d+)$/; $x <=> $1; } keys %fstats) {
    # Create a hash to hold all the timings info
    my %t = (%{$fstats{$key}}, %{$gstats{ $fstats{$key}{gxml} }});

    my @printme = ($key, @t{ qw( iter group size seqs) });
    my @headings = qw( filebase iter group size seqs);

    # Convert times into Seconds
    foreach ( @times ) {
      $t{$_."S"} = UnixDate($t{$_}, "%s");
      push(@printme, $t{$_}, $t{$_."S"});
      push(@headings, $_, $_."S");
    }
  
    # Calculate deltas, keeping track if any are negative
    for (my $i=1; $i < @times; ++$i) {
       my ($begin, $end) = ($times[$i-1], $times[$i]);

       my $delta_name = "d_"."$begin"."_$end";
 
       $t{$delta_name} = DateCalc($t{$begin}, $t{$end});

       $t{$delta_name."S"} = Delta_Format($t{$delta_name}, 0, "%sh");
 
       push(@printme, $t{$delta_name}, $t{$delta_name."S"});
       push(@headings, $delta_name, $delta_name."S");

       ++$negs{$delta_name} if ( $t{$delta_name."S"} < 0 );
    
    }

    # print row?
    unless ($is_summary) {
      unless ($done_header) {
        print join("\t", @headings)."\n";
        $done_header++;
      }
      print join("\t", @printme)."\n";

    }

    %{$fstats{$key}} = %t;
#    die Dumper(%t);

  }

  if ($is_summary) {

  print "Component: $pipeline_xml\n";
  print "Rows: ".(keys %fstats)."\n";

  print "\nThe following deltas had negatives:\n";
  foreach (keys %negs) {
    print "$_:\t$negs{$_}\n";
  }

  stats_for_column("size");

  print "\n";
  stats_for_column("d_qsubStart_prologStartS", 1);

  print "\n";
  stats_for_column("d_prologStart_execStartS", 1);

  print "\n";
  stats_for_column("d_execStart_execEndS", 1);

  print "\n";
  stats_for_column("d_execStart_execEndS", 1);

  print "\n";
  stats_for_column("d_execEnd_epilogStartS", 1);

  print "\n";
  stats_for_column("d_epilogStart_epilogEndS", 1);
  }

}

sub stats_for_column {
  my ($heading, $is_time) = @_;
 
  # get all values for key2 == $heading, regardless of key1
  my @vals = map { $_->{$heading} } values %fstats;

  my $stat = Statistics::Descriptive::Full->new();

  $stat->add_data( @vals ); 


  my $mean = $stat->mean();
#  my $var  = $stat->variance();
#  my $tm   = $stat->trimmed_mean(.25);
  my $std  = $stat->standard_deviation();
  my $median = $stat->median();
  my $min = $stat->min();
  my $max = $stat->max();
  my $sum = $stat->sum();
  $Statistics::Descriptive::Tolerance = 1e-10;

  if ($is_time) {
  print "\nStats for: $heading\n";
  my @parts = gmtime($mean);
  printf "Mean:\t%12.2f\t%4dd %4dh %4dm %4ds\n", $mean, @parts[7,2,1,0];
  @parts = gmtime($median);
  printf "Median:\t%12.2f\t%4dd %4dh %4dm %4ds\n", $median, @parts[7,2,1,0];
  @parts = gmtime($std);
  printf "StDev:\t%12.2f\t%4dd %4dh %4dm %4ds\n", $std, @parts[7,2,1,0];
  @parts = gmtime($min);
  printf "Min:\t%12.2f\t%4dd %4dh %4dm %4ds\n", $min, @parts[7,2,1,0];
  @parts = gmtime($max);
  printf "Max:\t%12.2f\t%4dd %4dh %4dm %4ds\n", $max, @parts[7,2,1,0];
  @parts = gmtime($sum);
  printf "Sum:\t%12.2f\t%4dd %4dh %4dm %4ds\n", $sum, @parts[7,2,1,0];
  } else {
  print "\nStats for: $heading\n";
  printf "Mean:\t%12.2f\n", $mean;
  printf "Median:\t%12.2f\n", $median;
  printf "StDev:\t%12.2f\n", $std;
  printf "Min:\t%12.2f\n", $min;
  printf "Max:\t%12.2f\n", $max;
  printf "Sum:\t%12.2f\n", $sum;
  }

}


#print $q->b("some text");


sub process_subflowgroup {
    my ($twig, $commandSet) = @_;

    ## within the iN.xml file, each commandSet with a fileName defines a gN.xml file.

    ## do nothing unless this commandSet has a file-based subflow
    return unless $commandSet->has_child('fileName');

    my %sg_props;
    ## get the pointer to the gN file
    $sg_props{file} = $commandSet->first_child('fileName')->text();

    ## pull the group number out of the file name:
#    if ( $sg_props{file} =~ m|(\w+)/g\d+/g(\d+).xml| ) {
#        $sg_props{name} = "$1$2";
#        $sg_props{group_num} = $2;
#    }

    my $bname = basename($sg_props{file}, ('.xml.gz','xml', 'gz', 'xml.gz'));

    defined( $fstats{ $bname } ) || die "No filestats for $bname";

    $fstats{$bname}{group} =  $sg_props{file} =~ m|/i\d+/g(\d+)/| ? $1: '?';
    # get times
#    my $startTime = $commandSet->first_child('startTime')->text();
#    my $endTime = $commandSet->first_child('endTime')->text();

#    my $delta = DateCalc($startTime, $endTime);
#    my $secs = Delta_Format($delta, 0, "%sh");

    # perform date calculations elsewhere, just collect data here
#    $fstats{$bname}{delta} = $delta;
#    $fstats{$bname}{dsecs} = $secs;
    $fstats{$bname}{execStart} = $commandSet->first_child('startTime')->text();
    $fstats{$bname}{execEnd} = $commandSet->first_child('endTime')->text();
#    $fstats{$bname}{startSecs} = UnixDate($startTime, "%s");
#    $fstats{$bname}{endSecs} = UnixDate($endTime, "%s");

    my $gxml = $commandSet->parent->first_child( 'parentFileName' )->text;

    defined( $gstats{ $gxml } ) || die "No group stats for $gxml";

    $fstats{$bname}{gxml} = $gxml;

}

sub process_ixml_commandSet {
  my( $t, $commandSet )= @_;

  my $gxml = $commandSet->first_child('fileName')->text; #group xml

  $gstats{$gxml}{qsubStart} = $commandSet->first_child('startTime')->text;
  $gstats{$gxml}{epilogEnd} = $commandSet->first_child('endTime')->text;

  my $dce = $commandSet->first_child('dceSpec');
  $gstats{$gxml}{prologStart} = $dce->first_child('reqStartTime')->text;

  my $log = $dce->first_child('log')->text;

  # Pull epilog start time:
  # T~~~5087~~~1414~~~Tue Oct 20 22:50:03 UTC 2009~~~job finished on exec.q for www-data~~~ip-10-242-143-1.ec2.internal
  open (my $LOG, $log) || die "Unable to open log $log: $!";
  while (my $line = <$LOG>) {
    next unless $line =~ /^T/;
    my @row = split('~~~', $line);
    $gstats{$gxml}{epilogStart} = $row[3];
    last;
  }
  close($LOG);

#  print $gstats{$gxml}{qsubStart}." ".$gstats{$gxml}{epilogEnd} ." ". $gstats{$gxml}{prologStart}."\n";

}

sub process_name {
  my( $t, $name )= @_;
#  ++$num;
#  print $q->h2( $num ." ".$name->text );
  while ( $name->next_sibling( ) ) {
#	my $key = $param->child( 'key' )->text;
#	print $q->h3( $key );
	print $q->h3( 'test' );
  }
}
