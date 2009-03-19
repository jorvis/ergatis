package File::OpenFile;

use Carp;
require Exporter;

@ISA = qw(Exporter);
@EXPORT_OK = qw (open_file);

sub open_file {
    my $file = shift;
    my $direction = shift;
    my $fh;

    if( $direction eq 'out' ) {
        if( $file =~ /\.gz$/ ) {
            open( $fh, ">:gzip", $file ) or confess("can't open $file ($!)");
            print "using gzip\n";
        } else {
            open( $fh, "> $file" ) or confess("Can't open $file ($!)");
        }
    } elsif( $direction eq 'concat' ) {
        if( $file =~ /\.gz$/ ) {
            open( $fh, ">>:gzip", $file ) or confess("can't open $file ($!)");
            print "using gzip\n";
        } else {
            open( $fh, ">> $file" ) or confess("Can't open $file ($!)");
        }
    } elsif( $direction eq 'in' ) {

        if( -e $file ) {
            
            if( $file =~ /\.gz$/ ) {
                open( $fh, "<:gzip", $file ) or confess("can't open $file ($!)");
            } else {
                open( $fh, "< $file") or confess("can't open $file ($!)");
            } 
        } elsif( -e $file.".gz" ) {
            my $tmp = $file.".gz";
            open( $fh, "<:gzip", $tmp ) or confess("Can't open $tmp ($!)");
        } else {
            confess("Could not find $file or a gz version");
        }

    } else {
        confess("Please specifiy a direction.  'in', 'out', or 'concat'");
    }

    return $fh;
}

1;
