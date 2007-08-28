package File::OpenFile;

require Exporter;

@ISA = qw(Exporter);
@EXPORT_OK = qw (open_file);

sub open_file {
    my $file = shift;
    my $direction = shift;
    my $fh;

    if( $direction eq 'out' ) {
        if( $file =~ /\.gz$/ ) {
            open( $fh, ">:gzip", $file ) or die("can't open $file ($!)");
        } else {
            open( $fh, "> $file" ) or die("Can't open $file ($!)");
        }
    } elsif( $direction eq 'in' ) {

        if( -e $file ) {
            
            if( $file =~ /\.gz$/ ) {
                open( $fh, "<:gzip", $file ) or die("can't open $file ($!)");
            } else {
                open( $fh, "< $file") or die("can't open $file ($!)");
            } 
        } elsif( -e $file.".gz" ) {
            my $tmp = $file.".gz";
            open( $fh, "<:gzip", $tmp ) or die("Can't open $tmp ($!)");
        } else {
            die("Could not find $file or a gz version");
        }

    } else {
        die("Please specifiy a direction.  'in' or 'out'");
    }

    return $fh;
}

1;
