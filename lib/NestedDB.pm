package NestedDB;

use strict;
use warnings;
use Carp;
use Data::Dumper;
use MLDBM "DB_File";
use File::Copy;

{
sub whowasi { (caller(1))[3]  }
sub usage { croak("usage: @{[&whowasi]}(dbtype, filename)") }
sub _die { croak("@{[&whowasi]}: $_[0]") }
sub _warn { carp("@{[&whowasi]}: $_[0]") }
my $DEBUG = 0;
sub debug { $DEBUG = @_ ? shift : 1 }

## User provies path to a file
## Will use 2 files:
##     /path/to/file
##     /path/to/file.tci 
sub TIEHASH {
    _warn &whowasi if $DEBUG;
    my ($class, @args) = @_;
    my $f = $args[0] || usage;
    
    croak( "$f is not a valid ".__PACKAGE__." [no Tchar index file (expected $f.tci)]")
        if( -e $f && !-e "$f.tci");
    
    my %tci;
    tie(%tci, "MLDBM", "$f.tci") or _die("Can't tie $f.tci to hash");

    my $data;
    open($data, "+>> $f") || croak("Could not open file for writing: $f [$!]");

    my $self = {
        'tci' => \%tci,
        'data' => $data,
        'data_file' => $f,
        };
    
    #&debug;

    return bless( $self, $class );
}

sub DESTROY {
    _warn &whowasi if $DEBUG;
    my ($self) = @_;
    my $fh = $self->{'data'};
    close($fh);
    untie( %{$self->{'tci'}} );
}

sub FETCH {
    _warn &whowasi if $DEBUG;

    my ($self, $key) = @_;
    unless( exists( $self->{'tci'}->{$key} ) ) {
        _warn( "$key does not exist" );
        return undef;
    }

    my $loc = $self->{'tci'}->{$key};
    my $d;
    seek( $self->{'data'}, $loc->[0], 0 ) or _die("Problem seeking in data file");
    read( $self->{'data'}, $d, $loc->[1] );
    
    my $retval;
    eval( "\$retval = $d" );
    _die("Could not eval $d") unless( defined( $retval ) );
    return $retval;
}

sub STORE {
    my ($self, $key, $value) = @_;

    local $Data::Dumper::Indent = 0;
    local $Data::Dumper::Purity = 1;
    local $Data::Dumper::Terse = 1;

    my $s = Dumper( $value );

    _warn(&whowasi.": $key, $s") if $DEBUG;
    my $data = $self->{'data'};
    my $d = $s;
    my $p = tell($self->{'data'});
    my $l = length($d);
    
    #Seek to end of file and print
    seek( $data, 0, 2 );
    print $data $d;

    $self->{'tci'}->{$key} = [$p,$l];
 
}

sub DELETE {
    _warn &whowasi if $DEBUG;

    my ($self, $key) = @_;
    my $retval;
    $retval = $self->{'tci'}->{$key} if( exists( $self->{'tci'}->{$key} ) );
    delete( $self->{'tci'}->{$key} );
    return $retval;
}

sub CLEAR {
    _warn &whowasi if $DEBUG;

    my ($self) = @_;
    my $fh = $self->{'data'};
    close($fh);

    unlink( $self->{'data_file'} ) or _die("Can't unlink file $self->{'data_file'}");
    open($fh, "+>> $self->{'data_file'}") or _die("Can't open file $self->{'data_file'}: $!");
    $self->{'data'} = $fh;
    %{$self->{'tci'}} = ();
    1;
}

sub EXISTS {
    _warn &whowasi if $DEBUG;
    my ($self, $key) = @_;
    return exists $self->{'tci'}->{$key};
}

sub FIRSTKEY {
    _warn &whowasi if $DEBUG;
    my ($self) = @_;

    ## resets the internal iterator with no overhead
    keys %{$self->{'tci'}};
    my ($k,$v) = each %{$self->{'tci'}};
    $k;
}

sub NEXTKEY {
    _warn &whowasi if $DEBUG;
    my ($self) = @_;
    my ($k,$v) = each %{$self->{'tci'}};
    return () if(!defined( $k ) );
    $k;
}

sub OPTIMIZE {
   _warn &whowasi if $DEBUG;
   my ($self) = @_;
   
   #make a tmp file
   my $tmp_db = "/tmp/tchar.$$.db";
   my %new_data;
   tie(%new_data, "TcharFile", $tmp_db);
   
   foreach my $k ( keys %{$self->{'tci'}} ) {
       $new_data{$k} = $self->FETCH($k);
   }
   untie(%new_data);

   copy( $tmp_db, $self->{'data_file'} );
   copy( "$tmp_db.tci", "$self->{'data_file'}.tci" );

   unlink( $tmp_db );
   unlink( $tmp_db.".tci" );
}

}
1==1;
