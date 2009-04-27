#!/usr/bin/perl

=head1 NAME

go_mapping_file_to_mldbm.pl - Will include information from a go_mapping file (pfam2go or tigrfams2go)
    to an mldbm file.

=head1 SYNOPSIS

 USAGE: go_mapping_file_to_mldbm.pl
       --mldbm_file=/path/to/db
       --go_mapping_file=/path/to/pfam2go
    [  --handle_different_values=(die|concat|replace|skip) 
    ]

=head1 OPTIONS

B<--mldbm_file,-m>
    Path to an mldbm file (created by hmmlib_to_mldbm). 

B<--go_mapping_file,-g>
    File which maps hmm (Pfam or TIGRfam) to GO terms. 
    http://www.geneontology.org/GO.indices.shtml

    Can be a comma separated list of mapping files.

B<--handle_different_values,-a>
    If there are already GO terms associated with a certain HMM and 
    the terms in this lookup file differ from those already in the mldbm, how should these be handled?

    die = die if differences are found
    concat = add new go terms to the list
    replace = replace all go terms with the new ones
    skip = skip the new terms and leave the old.

    Default = die

=head1  DESCRIPTION

 This file is used to add GO information to the HMM lookup file used in applying
 annotation based on HMMs.  The go mapping file can be retrieved from the URL above.

 Will map the go terms to the hmms in the lookup association file.  If the hmm is not in
 the mldbm file it will be skipped (printing a warning). 
 
=head1  INPUT
    Describe the input

=head1 OUTPUT
    Describe the output

=head1  CONTACT

    Kevin Galens
    kgalens@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use MLDBM 'DB_File';
use Data::Dumper;
use Fcntl qw(O_RDWR);

my ($DIE, $CONCAT, $REPLACE, $SKIP) = (0..3);
my $diff_lookup = { "die" => $DIE,
                    'concat' => $CONCAT,
                    'replace' => $REPLACE,
                    'skip' => $SKIP };

my %options;
my $results = GetOptions (\%options,
                          'mldbm_file|m=s',
                          'go_mapping_file|g=s',
                          'handle_different_values|a=s',
                          'help|h',
                          );

my %hmm_lookup;
my $hpv;
my @mapping_files;

&check_options(\%options);

# we need to make a lookup for hmm accessions with
# .## format.  Example: PF04889.4 -> PF04889
my %additional_lookup = &create_additional_lookup( %hmm_lookup );

foreach my $map_file ( @mapping_files ) {
    #parse the mapping file
    my %to_go = &parse_mapping_file( $map_file );

    foreach my $hmm ( keys %to_go ) {
        my $acc = $hmm;
        if( !exists( $hmm_lookup{$acc} ) && exists( $additional_lookup{$acc} ) ) {
            $acc = $additional_lookup{$acc};
        }

        print "$acc\t$hmm\n";
        next unless( exists( $hmm_lookup{$acc} ) );

        if( @{$hmm_lookup{$acc}->{'go'}} == 0 ) {
            my $tmp = $hmm_lookup{$acc};
            $tmp->{'go'} = $to_go{$hmm};
            $hmm_lookup{$acc} = $tmp;
        } elsif( my $concatted = &is_different( $hmm_lookup{$acc}->{'go'}, $to_go{$hmm} ) ) {
            if( $hpv eq $DIE ) {
                die("Found different go terms for hmm ($acc) in $map_file and look up");
            } elsif( $hpv eq $SKIP ) {
                next;
            } elsif( $hpv eq $CONCAT ) {
                my $tmp = $hmm_lookup{$acc};
                $tmp->{'go'} = $concatted;
                $hmm_lookup{$acc} = $tmp;
            } elsif( $hpv eq $REPLACE ) {
                my $tmp = $hmm_lookup{$acc};
                $tmp->{'go'} = $to_go{$hmm};
                $hmm_lookup{$acc} = $tmp;
            } 
        }
    }
}

sub is_different {
    my ($arr1, $arr2) = @_;
    my $retval;

    #make them into hashes
    my (%h1,%h2);
    map { $h1{$_} = 1 } @{$arr1};
    map { $h2{$_} = 1 } @{$arr2};

    #combine them
    my %h3;
    map { $h3{$_} = 1 } (keys %h1, keys %h2);

    #check the sizes
    my $count = scalar(keys %h3);
    unless( $count == scalar(@{$arr1}) && $count == scalar(@{$arr2})) {
        map { push(@{$retval}, $_) } (keys %h3);
    } else {
        $retval = 0;
    }
    return $retval;
}

sub parse_mapping_file {
    my ($mapping_file) = @_;
    my %retval;
    open(IN, "< $mapping_file" ) or die("Can't open $mapping_file: $!");
    while(<IN>) {
        next if( /^\!/ );
        my @c = split(/\s+/);
        my (undef,$acc) = split(/\:/,$c[0]);
        my $go_term = pop @c;
        push( @{$retval{$acc}}, $go_term );
    }
    close(IN);

    return %retval;
}

sub create_additional_lookup {
    my %hmm_lookup = @_;
    my %retval;
    map { if( /^(\S+)\.\d+$/ ) { $retval{$1} = $_; } } keys %hmm_lookup;
    return %retval;
}

sub check_options {
   my $opts = shift;

   if( $opts->{'help'} ) {
       &_pod;
   }

   foreach my $req ( qw(mldbm_file go_mapping_file) ) {
       die("Option --$req is required") unless( $opts->{$req} );
   }

   tie(%hmm_lookup, "MLDBM", $opts->{'mldbm_file'}, O_RDWR ) or die("Could not tie $opts->{'mldbm_file'}");
   my $tmp = $opts->{'go_mapping_file'};
   @mapping_files = split(/,\s*/, $tmp );
   
   if( $opts->{'handle_different_values'} ) {
       if( exists( $diff_lookup->{$opts->{'handle_different_values'}} ) ) {
           $hpv = $diff_lookup->{$opts->{'handle_different_values'}};
       }
   }
   
}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => *STDERR} );
}
