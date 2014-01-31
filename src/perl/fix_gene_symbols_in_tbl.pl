#!/usr/bin/env perl

=head1 NAME

fix_gene_symbols_in_tbl.pl - Removes bad genes and duplicate genes from a TBL file
and makes product name substitutions also.

=head1 SYNOPSIS

 USAGE: fix_gene_symbols_in_tbl.pl
       --input_file=/path/to/file.tbl
       --output_file=/path/to/file_new.tbl
       --rules_file=/path/to/rules.conf
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>
    TBL file

B<--output_file -o>
    Location to store new TBL file
    
B<--rules_file, -r>
    Path to a file listing additional bad gene symbols.  They can be comma- or tab-delimited
        
B<--log,-l>
    Logfile.

B<--debug,-D>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

    Uses an existing TBL file and writes a new TBL file from it, removing duplicate gene 
    symbols and bad gene symbols, and subtituting certain product 
    keywords encountered with more relevant values.
 
=head1  INPUT
    A TBL file

=head1 OUTPUT
    A new TBL file

=head1  CONTACT

    Shaun Adkins
    sadkins@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Data::Dumper;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;

my $rules_file;
my @bad_syms;
####################################################

my %options;
my $results = GetOptions (\%options,
                         "input_file|i=s",
                         "output_file|o=s",
                         "rules_file|r=s",
                         "log|l=s",
                         "debug|d=s",
                         "help|h"
                          );

&check_options(\%options);
&parse_rules_list($rules_file);

open(OUT, "> $options{'output_file'}") or 
    &_log($ERROR, "Could not open $options{'output_file'} for writing: $!");
open( IN, "< $options{'input_file'}") or 
    &_log($ERROR, "Could not open $options{'input_file'}: $!");
chomp(my @lines = <IN>);
close(IN);

my $flag = 0;
my $sym;
my %all_syms;	# Keeps track of all gene symbols encountered
my %copy_syms;	# Keeps track of gene syms encountered more than once
foreach my $line( @lines ) {
    if( $line =~ /\d+\s+\d+\s+gene/ ) {
        if( $sym ) {
            if( exists( $all_syms{$sym} ) ) {
                $all_syms{$sym}++;
                if( exists( $copy_syms{$sym} ) ) {
                    $copy_syms{$sym}++;
                } else {
                    $copy_syms{$sym} = 2;
                }
            } else {
                $all_syms{$sym} = 1;
            }
        }
        undef $sym;
    } elsif( $line =~ /^\s+gene\s+(\w+)/ ) {
        $sym = $1;
    }
};

my %is_bad_sym;

undef %is_bad_sym;
for (@bad_syms) { $is_bad_sym{$_} = 1; }	#Making keys accessible for quick searching

my $bad =0;

foreach my $line (@lines) {
    if( $line =~ /^\s+gene\s+(\S+)/ ) {
        next if( exists( $copy_syms{$1} ) );
        my $gene = $1;
        if($gene =~ /^[A-Z]/){	#Prok genes start with lower-case letters
	    #print $line, "\t$gene\n";
	    $bad++;
	    next;
	}
    }
    
    if($line =~ m/conserved hypothetical protein|conserved domain protein|conserved protein/ ) {
    	$line =~ s/conserved hypothetical protein|conserved domain protein|conserved protein/hypothetical protein/gi;
    }
    if($line =~ m/tRNA-Pseudo/) {
	print "tRNA-Pseudo issue\n";
    }
    print OUT "$line\n";
}
print "Removed ".scalar(keys %copy_syms)." repeated gene symbols\n";
print "Removed $bad bad gene symbols\n";
print "$options{'output_file'}\n";
close(OUT);

sub parse_rules_list {
    my $file = shift;
    
    open RULE, $file or die "Cannot open file $file for reading: $!\n";
     while (<RULE>) {
         my $line = $_;
         chomp $line;
	 push @bad_syms, split (/[,\t]/, $line);
     }
     
    close RULE;
}

sub check_options {

   my $opts = shift;

   if( $opts->{'help'} ) {
       &_pod;
   }

   if( $opts->{'log'} ) {
       open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
   }

   foreach my $req ( qw(input_file output_file) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }
   
   $rules_file = $opts->{'rules_file'} if defined $opts->{'rules_file'};
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
