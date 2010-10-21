#!/usr/local/bin/perl

=head1 NAME

bsml2interevidence_fasta.pl - Searches intergenic and evidence-less regions
    
=head1 SYNOPSIS

 USAGE: interevidence_search.pl
       --input_bsml=/path/to/input.bsml,input.bsml2
       --input_list=/path/to/input.list,input.list2
       --evidence_list=/path/to/some/blast_results.list,/path/to/other.bsml.list
       --output=/path/to/dir
       --length_cutoff=100
     [ --log=/path/to/file.log
       --debug=4
       --help
     ]

=head1 OPTIONS

B<--input_bsml,-b>
    -b or -i are required. Path to input bsml file (gene describing). Can be a comma
    separated list of files.

B<--input_list,-i>
    -b or -i are required.  Path to input list of gene describing bsml files. Can be
    a comma separated list of files.

B<--evidence_list,-e>
    Required. Path to evidence list (or lists) of alignment data (i.e. HMM bsml
    or BER bsml)

B<--output,-o>
    Required. Path to output directory

B<--length_cutoff,-c>
    Optional [Default = 100].  Interevidence regions smaller than reported cutoff will not be 
    written to fasta files.

B<--log,-l>
    Logfile.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

 DESCRIPTION
 
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
use XML::Twig;
use File::OpenFile qw(open_file);
use Data::Dumper;

$|++;

my @input_files;
my @evidence_files;
my %sequences;
my $output_dir;
my $length_cutoff = 100;

my %options;
my $results = GetOptions (\%options,
                          'input_bsml|b=s',
                          'input_list|i=s',
                          'evidence_list|e=s',
                          'output|o=s',
                          'length_cutoff|c=s',
                          'log|l=s',
                          'help|h',
                          );

&check_options(\%options);

#create lookup of sequences
my %feature_lookup = &parse_input_files( @input_files );

#parse the evidence files
&parse_evidence_files( \@evidence_files, \%feature_lookup );

#remove all the records in the feature_lookup{'polypeptides'} hashref
#that have a value of zero for 'has_evidence'
map { delete( $feature_lookup{'polypeptides'}->{$_} ) } 
grep( $feature_lookup{'polypeptides'}->{$_}->{'has_evidence'} == 0, 
      keys %{$feature_lookup{'polypeptides'}} );

&print_intergenic_fasta( \%feature_lookup, \%sequences, $output_dir );


sub print_intergenic_fasta {
    my ($f_lookup, $seqs, $odir) = @_;
    my $polys = $f_lookup->{'polypeptides'};

    foreach my $mol ( keys %{$f_lookup->{'molecules'}} ) {

        my $ofh = open_file( "$odir/$mol.interevidence.fsa", "out" );

        if( @{$f_lookup->{'molecules'}->{$mol}} > 0 ) {
            print $ofh ">${mol}\n";
            print $ofh $1."\n" while( $seqs->{$mol}->{'residues'} =~ /(\w{1,60})/g );
        } else {

            my @sorted_ids = sort { 
                if( !exists( $polys->{$a} ) ) {
                    die("Could not find poly $a in hash\n");
                }
                $polys->{$a}->{'startpos'} <=> $polys->{$b}->{'startpos'};
            } grep( exists( $polys->{$_} ), @{$f_lookup->{'molecules'}->{$mol}} );

            #using interbase numbering (from bsml)
            my $left = 0;
            my $right;
            foreach my $pid ( @sorted_ids ) {
                my $p = $polys->{$pid};
                $right = $p->{'startpos'};
                
                if( ($right - $left) >= $length_cutoff ) {
                    my $subseq = substr( $seqs->{$mol}->{'residues'},
                                         $left, $right - $left );
                    print $ofh ">${mol}_$left-$right\n";
                    print $ofh $1."\n" while( $subseq =~ /(\w{1,60})/g );
                }
                $left = $p->{'endpos'};
            }

            if( ($right - $left) >= $length_cutoff ) {
                my $subseq = substr( $seqs->{$mol}->{'residues'},
                                     $left, $right - $left );
                print $ofh ">${mol}_$left-".length($seqs->{$mol}->{'residues'})."\n";
                print $ofh $1."\n" while( $subseq =~ /(\w{1,60})/g );
            }
        }
        close($ofh);
        print "$odir/$mol.interevidence.fsa\n";
    }        
    
}

sub parse_evidence_files {
    my ($ev_files, $f_lookup) = @_;
    
    foreach my $e ( @{$ev_files} ) {
        my $in = open_file( $e, "in" );
        chomp( my @spas = grep( $_ =~ /\<Seq-pair-alignment/, <$in> ) );
        close($in);
        map {
            my $cid = $1 if( $_ =~ /refseq=\"(\S+)\"/ );
            next unless( exists( $f_lookup->{'CDS'}->{$cid} ) );
            my $pid = $f_lookup->{'CDS'}->{$cid};
            $f_lookup->{'polypeptides'}->{$pid}->{'has_evidence'} = 1;
        } @spas;
    }
}
sub parse_input_files {
    my @input_files = @_;
    my %retval;

    my $twig = new XML::Twig( 'twig_roots' => {
        'Sequence' => sub {
            my ($t, $seq) = @_;
            my $sid = $seq->att('id');
            map {
                my $pid = $_->att('id');
                my $il = $_->first_child('Interval-loc');
                my ($left, $right) = ( $il->att('startpos'), $il->att('endpos') );
                $retval{'polypeptides'}->{$pid} = {
                    'startpos'   => $left,
                    'endpos'     => $right,
                    'has_evidence' => 0
                    };
                push( @{$retval{'molecules'}->{$sid}}, $pid );
            } $seq->find_nodes('//Feature[@class="polypeptide"]');

            #parse sequence info
            my $sdi = $seq->first_child("Seq-data-import");
            $sequences{$sid}->{'residues'} = &parse_sequence( $sdi->att('identifier'), $sdi->att('source') );
        },
        'Feature-group' => sub {
            my ($t, $fg) = @_;
            my $pid = ($fg->find_nodes('Feature-group-member[@feature-type="polypeptide"]'))[0]->att('featref');
            my $cid = ($fg->find_nodes('Feature-group-member[@feature-type="CDS"]'))[0]->att('featref');
            $retval{'polypeptides'}->{$pid}->{'CDS'} = $cid;
            $retval{'CDS'}->{$cid} = $pid;
        }
    } );

    map { $twig->parsefile( $_ ) } @input_files;
                              
    return %retval;
}

sub parse_sequence {
    my ($id, $file) = @_;
    my $flag = 0;
    my $seq;
    my $in = open_file( $file, "in");
    for( <$in> ) {
        chomp;
        if( /^>(\S+)/ ) {
            last if( $flag );
            $flag = 1;
        } elsif( $flag ) {
            $seq .= $_;
        }
    }
    return $seq;
}

sub check_options {

   my $opts = shift;

   if( $opts->{'help'} ) {
       &_pod;
   }
   
   my @reqs = qw( evidence_list output );
   foreach my $req ( @reqs ) {
       unless( $opts->{$req} ) {
           die("Option $req is required");
       }
   }

   $output_dir = $opts->{'output'};

   my @evidence_lists = split(/[,\s]+/, $opts->{'evidence_list'});
   map {
       my $in = open_file( $_, 'in' );
       chomp( my @tmp = <$in> );
       close($in);
       push(@evidence_files, @tmp);
   } @evidence_lists;

   if( $opts->{'input_list'} ) {
       map {
           my $in = open_file( $_, 'in' );
           chomp( my @tmp = <$in> );
           close($in);
           push(@input_files, @tmp );
       } split(/[,\s]+/, $opts->{'input_list'} );

   } elsif( $opts->{'input_bsml'} ) {
       map { push(@input_files, $_) } split(/[,\s]+/, $opts->{'input_bsml'} );
   } else {
       die("Either input_list or input_bsml option is required");
   }

}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => *STDERR} );
}
