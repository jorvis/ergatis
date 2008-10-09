#!/usr/bin/perl -w

use strict;
use Date::Manip;
use XML::Twig;

my %cmds;

my $source_file = shift || die "pass an XML file to parse\n";

process_file( $source_file );

print "command summary:\n";

my $total_time = 0;
my $total_count = 0;

foreach my $cmd ( keys %cmds ) {
    $total_time += $cmds{$cmd}{time};
    print "\t$cmd - count: $cmds{$cmd}{count} - time: $cmds{$cmd}{time}s\n";
}

print "\ntotal CPU time: ${total_time}s\n";
print "total command count: $total_count\n\n";

exit();

sub process_file {
    my $path = shift;
    
    print "processing $path\n";
    
    my $fh = get_conditional_read_fh($path);

    if ( $fh ) {
        my $t = XML::Twig->new(
                           twig_roots => {
                               'command' => sub {
                                   my ($t, $elt) = @_;
                                   $cmds{ $elt->first_child('name')->text() }{time} += time_info( $elt );
                                   $cmds{ $elt->first_child('name')->text() }{count}++;
                               },
                               'commandSet/fileName' => sub {
                                   my ($t, $elt) = @_;
                                   process_file( $elt->text );
                               },
                           },);

        $t->parse( $fh ); 
    } else {
        print "WARNING: skipping absent $path\n";
    }
}

sub get_conditional_read_fh {
    my $path = shift;
    my $fh;
    my $found = 0;
    
    if ( -e $path ) {
        $found = 1;
    
    ## auto-handle gzipped files
    } elsif (! -e $path && -e "$path.gz") {
        $path .= '.gz';
        $found = 1;
    } 
    
    if (! $found ) {
        ## we can't find the file.  just return an empty file handle. 
        ##  the process sub is OK with this.
        return $fh;
    }

    if ( $path =~ /\.gz$/ ) {
        open($fh, "<:gzip", $path) || die "can't read file $path: $!";
    } else {
        open($fh, "<$path") || die "can't read file $path: $!";
    }
    
    return $fh;
}


sub time_info {
    my $command = shift;
    
    ## make sure we can at least get start time
    if (! $command->first_child('startTime') ) {
        return 0;
    }
    
    my $state = $command->first_child('state')->text;
    my $start_time_obj = ParseDate( $command->first_child('startTime')->text );
    my $start_time     = UnixDate( $start_time_obj, "%c" );
    
    my ($end_time_obj, $end_time);
    ## end time may not exist (if running, for example)
    if ( $command->first_child('endTime') ) {
        $end_time_obj   = ParseDate( $command->first_child('endTime')->text );
        $end_time       = UnixDate( $end_time_obj, "%c" );
    }

    ## doing it here manually because strftime was behaving badly (or I coded it badly)
    if ($start_time_obj) {
        my $diffstring;
        
        if ($end_time_obj) {
            $diffstring = DateCalc($start_time_obj, $end_time_obj, undef, 1);
        } else {
            $diffstring = DateCalc($start_time_obj, "now", undef, 1);
        }
        
        ## take out any non \d: characters
        $diffstring =~ s/[^0-9\:]//g;

        my @parts = split(/:/, $diffstring);
        
        ## years + months + weeks + days
        my $days = ($parts[0] * 365) + ($parts[1] * 30) + ($parts[2] * 7) + ($parts[3]);
        
        # ( days + hours + min + sec )
        return ($days * 86400) + ($parts[4] * 3600) + ($parts[5] * 60) + $parts[6];
    } else {
        return 0;
    }
}
