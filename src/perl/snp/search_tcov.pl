#! /usr/local/bin/perl -w

############################################################################
#### Copyright (c) 2002-2003 Jacques Ravel, Pawel Gajer and Timothy Read
###
#############################################################################
### This is UNSUPPORTED code!  I hope you find it useful, but if you can't
### get it to run, you'll have to find someone local to help you.  I'm afraid
### we don't have time to respond to help inquiries.
#############################################################################
###
###                       search_tcov.pl
###
###    This script will display a slice of the coverage information for the 
###    B. anthracis genome sequence released on May 1st, 2003.
###
############################################################################


use strict;
use Getopt::Long;
$| = 1;

### VARIABLE DEFINITION

my $pos;
my $help;
my $newcov;
my @quality;
my $tcov2 = "BA.tcov";  #### IF NECESSARY ADD PATH TO THIS FILE
my $tcov = "index.tcov";#### IF NECESSARY ADD PATH TO THIS FILE
my $add;
my $list;
my @coverage;
my $fh;
my @list;
my $record;
my @location;
my @loca;
my %locat;
my %quali;
my $qualine;

my $refscore;
my $refqual;
my $refcov;
my $range;
my @cov3;


### COMMAND LINE INPUT

Getopt::Long::config("no_ignore_case");

my $result = GetOptions ("h|help"     => \$help,
			 "R=s"        => \$range,
			 "P=s"        => \$pos);

my $HELPTEXT = "
# Copyright (c) 2002-2003 Jacques Ravel, Pawel Gajer and Timothy Read
 
#############################################################################
This is UNSUPPORTED code!  I hope you find it useful, but if you can't
get it to run, you'll have to find someone local to help you.  I'm afraid
we don't have time to respond to help inquiries.
#############################################################################

This script will display a slice of the coverage information for the B. anthracis
genome sequence released on May 1st, 2003.


USAGE -R range -P position 

-R number of bp to display before and after the position

-P position on the genome where the slice will be centered on.

OUTPUT is to STDOUT

See README file on .tcov for OUTPUT explaination.

Example: search_tcov.pl -R 5 -P 5000

  B. anthracis Ames chromosome:
  4995  T       1108    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT       4:34:4:24:45:51:36:51:40:37:28:41:40:37:29:41:36:24:33:4:34:40:41:41:45:33:42:36:16:36:36:34:35
  4996  T       1101    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT       16:39:4:21:45:51:37:51:42:37:29:41:34:35:24:45:35:24:30:4:32:40:41:39:45:33:33:36:17:36:36:34:35
  4997  T       1082    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT       4:34:4:24:51:51:37:51:40:37:32:41:34:35:20:45:32:29:30:4:32:40:41:40:45:34:33:36:14:35:35:34:28
  4998  G       1062    GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG       4:39:4:27:51:51:34:51:37:37:24:41:34:34:25:41:22:30:28:18:29:40:34:39:38:27:33:36:23:36:35:36:24
  4999  T       1072    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT       4:39:21:26:51:51:37:51:44:34:29:30:37:34:18:41:27:26:26:21:27:33:25:39:45:33:33:35:13:36:35:36:35
 *5000  A       1064    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA       4:51:19:28:51:51:37:51:38:34:29:30:37:34:18:41:27:31:24:16:27:30:31:39:45:19:33:36:20:36:34:36:27
  5001  C       1022    CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       4:51:19:25:51:51:37:51:45:34:24:34:37:29:18:41:21:26:22:15:24:25:4:39:45:23:33:36:21:34:34:36:33
  5002  A       1064    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA       4:39:19:32:45:40:37:45:45:34:24:41:33:32:32:37:21:31:28:20:32:26:27:41:45:21:40:36:17:34:34:36:36
  5003  T       908     TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT        4:32:20:21:37:39:36:45:37:29:24:31:30:24:32:30:4:21:18:15:20:25:4:39:45:31:37:36:35:35:36:36
  5004  C       929     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        4:30:28:25:37:42:37:45:37:34:32:30:30:18:32:30:18:22:18:4:28:26:4:39:38:31:37:36:35:35:36:31
  5005  G       870     GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG        4:30:20:4:37:40:36:45:40:34:32:30:30:21:28:30:4:17:4:4:22:34:4:51:38:22:37:34:34:34:34:36
  COV: 33 CBQUAL: 1064 QUAL: 32.24

NECESSARY FILES:

BA.tcov  and index.tcov
Both files can be dwnloaded from at ftp://tigr.org/pub/data/b_anthracis_ames/
Make sure these two files are in the same directory as this script or indicate the path in the 
variable definition section above.
The offset (index.tcov file )is used in this script to retrieve information from the .tcov file.

HELP: -h (This text)


";


if ($result == 0) {
    die ("Command line parsing failed\n")
    }


if (! defined $pos) {
    print STDERR "You must specify a position (-P), use option -h for help\n";
    exit(0);
}

if (! defined $range) {
    print STDERR "You must specify a range (-R), use option -h for help\n";
    exit(0);
}


### OPEN  THE INDEX FILE (index.tcov) Store it into an array

open(QUAL, $tcov);
	
@coverage = <QUAL>;

close (QUAL);

### OPEN THE UNGAPPED.tcov file (just it into filehandle)

open($fh, $tcov2);


######################################################################
if ($pos < 1013547) {
    $pos = $pos + 576235;
}elsif ( $pos > 1013546 && $pos < 1288255) {
    $pos = $pos + 576234;
}elsif ( $pos > 1288254 && $pos < 3390649) {
    $pos = $pos + 576235;
}elsif ( $pos > 3390648 && $pos < 4651059) {
    $pos = $pos + 576236;
}elsif ( $pos > 4651058) {
    $pos = $pos - 4651058;
}
#######################################################################
     
#### CREATE AN ARRAY WITH ONLY $range  LINE BEFORE AND AFTER

my $i;

push(@cov3, "B. anthracis Ames chromosome:\n");

for ($i = $pos-($range+1) ; $i < $pos+$range ; ++$i) {

    ### grab the offset from index.tcov
    my $posoff = $coverage[$i];
	
    ### position the cursor to the position offset
    seek($fh, $posoff, 0);
    
    ### Grab the line with subroutine get_line
    
    $record = get_line($fh);
    
    ### format for output (put a * in front of right position)
    
    my @temp4 = split(" ", $record);
    
    my $pos_t = pos_trans2($temp4[0]);
    
    $record = "$pos_t\t$temp4[1]\t$temp4[2]\t$temp4[3]\t$temp4[4]\n";
    
    if ($i == $pos-1) {
	
	$refcov = length($temp4[3]);
	$refscore = $temp4[2]/$refcov;
	$refscore = int($refscore*100)/100;
	$refqual = $temp4[2];
	$record = "*" . $record;
	
    }else{
	
	$record = " " . $record;
    }
    
    push (@cov3, $record);
    
}

push (@cov3, "COV: $refcov CBQUAL: $refqual QUAL: $refscore\n");

### print to STDOUT

print "@cov3";

close ($fh);

exit;


##############################################################################
###               SUBROUTINES
##############################################################################


sub get_line {


    my ($fh) = @_;

    my $offset;
    my $record = '';
    my ($save_input_separator) = $/;

    $/ = "\n";

    $record = <$fh>;
    $/ = $save_input_separator;

    return $record;

}

sub get_next_record {

    my ($fh) = @_;

    my $offset;
    my $record = '';
    my ($save_input_separator) = $/;

    $/ = "\/\/\n";

    $record = <$fh>;
    $/ = $save_input_separator;

    return $record;

}


sub pos_trans2 {

 my $pos = shift;

    my $pos2;

    if ($pos > 576235 && $pos < 1589782) {
	$pos2 = $pos - 576235;
    }elsif ( $pos > 1589781 && $pos < 1864489) {
	$pos2 = $pos - 576234;
    }elsif ( $pos > 1864488 && $pos < 3966884) {
	$pos2 = $pos - 576235;
    }elsif ( $pos > 3966883 && $pos < 5227295) {
	$pos2 = $pos - 576236;
    }elsif ( $pos > 0 && $pos < 576236 ) {
	$pos2 = $pos + 4651058;
    }

    return $pos2;

}



