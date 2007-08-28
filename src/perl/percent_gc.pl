#!/usr/local/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

percent_gc.pl - calculates the percent GC content of nucleotide sequences

=head1 SYNOPSIS

USAGE: 
    outer_membrane_domain.pl 
        --input=/path/to/input_file
        --output_bsml=/path/to/output.bsml
        --sourcename=/path/to
      [ --class=assembly 
        --log
        --debug
        --help ]
      

=head1 OPTIONS

B<--input,-i>
    REQUIRED. Can be a single or multi fasta file, or bsml file.

B<--output_bsml,-o>
    REQUIRED. The output bsml file.

B<--class,-c>
    OPTIONAL.  The class of sequence to calculated gc percentage for.
    Used with bsml only.

B<--sourcename,-s>
    REQUIRED. The value used in the analysis section of the bsml

B<--log,-l>
    Logfile.

B<--debug,-d>
    Higher number = more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

    Calculates the GC percentage of nucleotides sequences and prints the results to
    a bsml file.

=head1  INPUT

    Single or multi fasta file.

=head1 OUTPUT

    Bsml.

=head1  CONTACT

    Kevin Galens
    kgalens@jcvi.org

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use BSML::BsmlBuilder;
use Ergatis::Logger;
use XML::Twig;

####### GLOBALS AND CONSTANTS ###########
my $input;                              #Holds the input file
my $input_is_fsa = 0;                   #Flag.  True if input is fasta
my $output;                             #Output file
my $class = ['all'];                    #The class of element to calculate percent_gc
my $gc_content;                         #Holds percent gc for sequences
my $debug;                              #The debug variable
my $sourcename;                         #For the analysis section of the bsml
#########################################

my %options = ();
my $results = GetOptions (\%options, 
                          'input|i=s',
                          'class|c=s',
                          'output_bsml|o=s',
                          'sourcename|s=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || &_pod;

#Setup the logger
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

# Check the options.
&check_parameters(\%options);

if( $input_is_fsa ) {
    $gc_content = &get_gc_from_fsa_input( $input );
} else {
    $gc_content = &get_gc_from_bsml_input( $input, $class );
}

&generate_bsml( $gc_content, $output );
print "Fin!\n" if( $logger->is_debug );


######################## SUB ROUTINES #######################################
sub get_gc_from_fsa_input {
    my $input = shift;
    my $retval;
    
    my $fh = &open_file( $input, 'in' );
    my $seq;

    while( my $line = <$fh> ) {
        
        if( $line =~ /^>((\S+)\s+.*)/ ) {
            if( $seq ) {
                my $seqData = $retval->{'Sequence'}->{$seq}->{'seq-data'};
                $retval->{'Sequence'}->{$seq}->{'total_length'} = 
                    length($seqData );
                $seqData =~ s/[^GCgc]//g;
                $retval->{'Sequence'}->{$seq}->{'percent_gc'} = 
                    (int((length( $seqData )/$retval->{'Sequence'}->{$seq}->{'total_length'})*10000))/100;
            }
            $seq = $2;
            $retval->{'Sequence'}->{$seq}->{'defline'} = $1;
        } else {
            chomp $line;
            $retval->{'Sequence'}->{$seq}->{'seq-data'} .= $line;
        }

    }
    close( $fh );

    my $seqData = $retval->{'Sequence'}->{$seq}->{'seq-data'};
    $retval->{'Sequence'}->{$seq}->{'total_length'} = length( $seqData );
    $seqData =~ s/[^GCgc]//g;
    $retval->{'Sequence'}->{$seq}->{'percent_gc'} = 
        (int((length( $seqData )/$retval->{'Sequence'}->{$seq}->{'total_length'})*10000))/100;
    $retval->{'Sequence'}->{$seq}->{'source'} = $input;
    $retval->{'Sequence'}->{$seq}->{'identifier'} = $seq;
    
    return $retval;
}

sub get_gc_from_bsml_input {
    my ($input, $class) = @_;

    my $retval = {};

    my $twig = new XML::Twig( 'twig_handlers' => {
        'Sequence' => sub { 
            my $seq_id = $_[1]->att('id');
            my $molecule = $_[1]->att('molecule');
            my $seq_class = $_[1]->att('class');
            $logger->warn( "Skipping $seq_id because it is an amino acid sequence")
                if( $molecule && $molecule eq 'aa' );
            $logger->warn( "Skipping $seq_id because it does not have a class" )
                unless( $seq_class );

            return unless( map(/($seq_class|all)/, @{$class}) );

            my $seq_data_import = $_[1]->first_child('Seq-data-import');
            my $ident = $seq_data_import->att('identifier');
            my $source = $seq_data_import->att('source');

            my $bh = &open_file( $source, "in" );

            my $seq = "";
            my $flag = 0;
            while( my $line = <$bh> ) {
                chomp $line;

                if( $line =~ /^>/ ) {
                    $flag = 0;
                    $flag = 1 if( $line =~ /^>$ident/ );
                } elsif( $flag ) {
                    $seq .= $line;
                }

            }
            close( $bh );

            $retval->{'Sequence'}->{$seq_id}->{'total_length'} = length( $seq );
            $seq =~ s/[^GCgc]//g;
            $retval->{'Sequence'}->{$seq_id}->{'percent_gc'} = 
                (int((length( $seq )/$retval->{'Sequence'}->{$seq_id}->{'total_length'})*10000))/100;
            $retval->{'Sequence'}->{$seq_id}->{'source'} = $source;
            $retval->{'Sequence'}->{$seq_id}->{'identifier'} = $ident;
            my ($defline) = $_[1]->find_nodes( '//Attribute[@name="defline"]' );
            $retval->{'Sequence'}->{$seq_id}->{'defline'} = $defline->att('content')
                if( $defline );
            $retval->{'Sequence'}->{$seq_id}->{'class'} = $seq_class;
            
        }

    });

    my $fh = &open_file( $input, 'in' );
    $twig->parse( $fh );
    close( $fh );

    return $retval;
    
    
}

sub generate_bsml {
    my ( $data, $outfile ) = @_;

    my $bsml = new BSML::BsmlBuilder;

    foreach my $seq_id ( keys %{$data->{'Sequence'}} ) {
        
        #Add the sequence.
        my $bsml_seq = $bsml->createAndAddSequence( "$seq_id", "$seq_id", 
                                                    $data->{'Sequence'}->{$seq_id}->{'total_length'},
                                                    'na', $data->{'Sequence'}->{$seq_id}->{'class'});
        
        #Add the Seq-data-import
        my $sdi = $bsml->createAndAddSeqDataImport( $bsml_seq, 'fasta',
                                                 $data->{'Sequence'}->{$seq_id}->{'source'}, '',
                                                 $data->{'Sequence'}->{$seq_id}->{'identifier'} );

        #Add the percent gc
        $bsml->createAndAddBsmlAttribute( $bsml_seq, 'percent_gc', 
                                          $data->{'Sequence'}->{$seq_id}->{'percent_gc'} );

        #Add the defline attribute
        $bsml->createAndAddBsmlAttribute( $bsml_seq, 'defline',
                                          $data->{'Sequence'}->{$seq_id}->{'defline'} );

        #Add the link
        $bsml->createAndAddLink( $bsml_seq, 'analysis', '#percent_gc_analysis', 
                                     'input_of' );
        

    }

    $bsml->createAndAddAnalysis( ( 
        'id'             => 'percent_gc_analysis',
        'programversion' => 'current',
        'sourcename'     => $sourcename,
        'algorithm'      => 'percent_gc' ) );

    $bsml->write( $outfile );
    print "Wrote to $outfile\n" if( $logger->is_debug );
    
}

sub is_fsa {
    my $file = shift;
    my $retval;

    my $fh = &open_file( $file, 'in' );
    my ($first_line) = <$fh>;
    close( $fh );
    
    if( $first_line =~ /^\s*</ ) {
        $retval = 0;
    } elsif( $first_line =~ /^>/ ) {
        $retval = 1;
    } else {
        $logger->logdie("Cannot determine the format of file $file");
    }
    
    return $retval;
}

sub check_parameters {
    my $opts = shift;

    &_pod if($opts->{'help'});

    if($opts->{'input'}) {
        $logger->logdie("Option input_file ($opts->{'input'}) does not exist") 
            unless( -e $opts->{'input'} );
        $input = $opts->{'input'};
        $input_is_fsa = &is_fsa( $input );
    } else {
        $logger->logdie("Option --input is required");
    }

    print "Using $input as input\n" if( $logger->is_debug );

    unless($opts->{'output_bsml'}) {
        $logger->logdie("Option --output_bsml is required");
    } else {
        $output = $opts->{'output_bsml'};
    }

    print "Using $output as output\n" if( $logger->is_debug );

    unless( $input_is_fsa ) {
        
        if( $opts->{'class'} ) {
            my @tmp = split( /\,/, $opts->{'class'} );
            $class = \@tmp;
        }        

        print "\tparsing @{$class} classes\n" if( $logger->is_debug );
        $"=" ";
    } else {
        print "Ignoring the class option because using fsa input\n" 
            if($logger->is_debug);
    }

    unless( $opts->{'sourcename'} ) {
        $logger->logdie("Option --sourcename is required");
    } else {
        $sourcename = $opts->{'sourcename'};
    }

    print "Sourcename for this analysis is $sourcename\n" if( $logger->is_debug );
    
    if($opts->{'debug'}) {
        $debug = $opts->{'debug'};
    }
    
}

sub open_file {
    my ($file, $direction) = @_;
    
    my $arrow = "<";
    $arrow = ">" if( $direction eq 'out' );

    if( $file =~ /\.gz/ ) {
        $arrow .= ":gzip $file";
    } else {
        $arrow .= " $file";
    }

    my $fh;
    open( $fh, $arrow ) or $logger->logdie("Could not open $file ($!)");

    return $fh;
}

sub _pod {   
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
