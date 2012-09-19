#!/usr/bin/env perl

=head1 NAME

mugsy_concat_fasta.pl - Description

=head1 SYNOPSIS

 USAGE: mugsy_concat_fasta.pl
       --input_fasta_list=/path/to/fasta.list
       --input_bsml_list=/path/to/bsml.list
       --output_directory=/path/to/output_dir
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1  DESCRIPTION

    Script will create a individual fasta file for all sequences of the same organism. This is if the genus, species and strain
    match exactly. The organism is parsed from the BSML input and the ids should be reference in the BSML files.
 
=head1  CONTACT

    Kevin Galens
    kgalens@gmail.com

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use XML::Parser;
use Pod::Usage;
use File::OpenFile qw(open_file);

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
####################################################

my %options;
my $results = GetOptions (\%options,
			  "input_fasta_list|f=s",
			  "input_bsml_list|b=s",
			  "output_directory|o=s",
			  "log|l=s",
			  "debug|d=s",
			  "help|h"
                          );

&check_options(\%options);

my $id_lookup = &parse_bsml($options{'input_bsml_list'});
&concat_fsa( $options{'input_fasta_list'}, $id_lookup, $options{'output_directory'} );

sub concat_fsa {
    my ($fsa_list, $id_lookup, $outdir) = @_;
    open( IN, "< $fsa_list") or die("Can't open $fsa_list: $!");
    chomp( my @fsa_files = <IN> );
    close(IN);
    
    my %outfhs;
    
    foreach my $fsa_file ( @fsa_files ) {
	open(FSA, "< $fsa_file") or die("Could not open $fsa_file: $!");

	my ($seq_id, $seq);
	while( my $line = <FSA> ) {
	    next if( $line =~ /^\s+$/ );
	    chomp $line;
	    if( $line =~ /^>(\S+)/ ) {
		my $new_seq_id = $1;
		print "Starting sequence $new_seq_id\n";
		if( defined( $seq ) && $seq ne "" ) {
		    print "Writing sequence $seq_id\n";
		    &write_seq( $seq_id, $seq, \%outfhs, $id_lookup, $outdir );
		    undef $seq;
		}
		$seq_id = $new_seq_id;
	    } else {
		$seq .= $line;
	    }
	}

	&write_seq( $seq_id, $seq, \%outfhs, $id_lookup, $outdir ) if( defined( $seq ) && $seq ne "" );

	close(FSA);
    }

    # Close all the output file handles;
    close for( values %outfhs );
    
}

sub write_seq {
    my ($id, $seq, $outfhs, $id_lookup, $outdir) = @_;

    die("Could not find organism for sequence $id. This sequence was not found ".
	"in the input bsml") unless( exists( $id_lookup->{$id} ) );
    my $org = $id_lookup->{$id};
    my $filename = &org_to_filename( $org );
    my $fh;
    unless( exists( $outfhs->{$filename} ) ) {
	my $outfile = $outdir."/$filename.fsa";
	open($fh, "> $outfile") or die("Can't open $outfile for writing: $!");
	$outfhs->{$filename} = $fh;
    } else {
	$fh = $outfhs->{$filename};
    }
    print $fh ">$id\n";
    print $fh $1."\n" while( $seq =~ /(\w{1,60})/g );
    return 1;
    
}

# Remove all spaces and replace them with underscores.
# Organism name should not have a period. This will create
# issues when trying to parse the Mugsy MAF file because
# it will use the organism name (parsed from the fasta file
# name) concatenated with the molecule name by a period. For 
# example: PY_76.PY_76.assembly.1. Later we parse out everything
# to the left of the first period as the organism name.
sub org_to_filename {
    my ($org) = @_;
    $org =~ s/[\s\/]+/_/g;
    $org =~ s/\.//g;
    $org;
}

# Will break if there are nested Sequence elements. 
# Shouldn't happen though?
sub parse_bsml {
    my ($list) = @_;
    
    my ($genus, $species, $strain);
    my @sequence_ids;
    my $strain_flag = 0; # Should let us know when we are in a strain element. 
    my $start_subs = {
	'Organism' => sub {
	    my ($expat,$elt,%params) = @_;
	    undef $strain;
	    $genus = $params{'genus'} or die("There was no genus attribute for organism element");
	    $species = $params{'species'} or die("There was no species attribute for organism element");
	},
	'Strain' => sub { $strain_flag = 1 },
	'Attribute' => sub { 
	    my ($expat,$elt,%params) = @_;
	    return unless( $strain_flag && $params{'name'} eq 'name' );
	    $strain = $params{'content'} or die("Could not parse strain from attribute");
	},
	'Sequence' => sub {
	    my ($expat,$elt,%params) = @_;
	    return unless( $params{'class'} eq 'assembly' );
	    push(@sequence_ids, $params{'id'});
	}
    };
    my $end_subs = {
	'Strain' => sub { $strain_flag = 0 }
    };

    my $parser = new XML::Parser('Handlers' => {
	'Start' => sub { $start_subs->{$_[1]}(@_) if( exists($start_subs->{$_[1]}) ) },
	'End' => sub { $end_subs->{$_[1]}(@_) if( exists($end_subs->{$_[1]}) ) }
    });

    open(BSML, "< $list") or die("Could not open $list: $!");
    chomp( my @files = <BSML> );
    close(BSML);

    my $retval = {};
    foreach my $file ( @files ) {
	my $fh = open_file($file, "in");
	$parser->parse( $fh );
	close($fh);
	my $org = "$genus $species";
	$org .= " $strain" if( $strain );
	foreach my $seq_id ( @sequence_ids ) {
	    die("Already found sequence with id of $seq_id ($file)")
		if( exists( $retval->{$seq_id} ) );
	    $retval->{$seq_id} = $org;
	}
	@sequence_ids = ();
	$org = "";
    }

    return $retval;
}

sub check_options {

   my $opts = shift;

   if( $opts->{'help'} ) {
       &_pod;
   }

   if( $opts->{'log'} ) {
       open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
   }

   $debug = $opts->{'debug'} if( $opts->{'debug'} );

   foreach my $req ( qw(input_fasta_list output_directory input_bsml_list) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }
}

sub _log {
   my ($level, $msg) = @_;
   if( $level <= $debug ) {
      print STDOUT "$msg\n";
   }
   print $logfh "$msg\n" if( defined( $logfh ) );
   exit(1) if( $level == $ERROR );
}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
