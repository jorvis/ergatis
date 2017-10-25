#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use IO::File;

use Blast::BlastHitDataType;

my %p_values	=	();
my %map_hash = ();
my $praze_btab	=	*STDIN;
my $out		=	*STDOUT;

my $IDENTITY_CUTOFF = 0.5;
my $SIMILARITY_CUTOFF = 0.5;
my $MATCH_CUTOFF = 0.1;	# Match % isn't stored in btab... only in BER raw output

&parse_options;
&add_pvalues;

sub print_usage {
	die << "END";
usage: add_pvalue_to_praze_output.pl -blast_btab|b <blast_btab_data>
	-praze_btab|p <praze_brab_data> [-output|o <output>]
	[-help|h]
END
}

sub parse_options {
	my %opts = ();
	GetOptions(\%opts, "blast_btab|b=s", "praze_btab|p=s", "mapping_file|m=s",
		   "output|o=s", "help|h");
	print_usage if $opts{help};
	
	if ($opts{mapping_file}) {
		print "Mapping polypeptide features to corresponding CDS features...\n";
		map_polypeptide_to_CDS($opts{mapping_file});
	}
	initialize_pvalues($opts{blast_btab}, \%map_hash) if ($opts{blast_btab});
    
    ## if the input file doesn't exist, touch it.  This is a dirty hack used
    ##  only because praze doesn't create an empty output file if there are no
    ##  matches, which it should.
    if ( ! -e $opts{praze_btab} ) {
        open(my $infh, ">$opts{praze_btab}") || die "can't create empty btab output file from praze: $!";
        close $infh;
    }
    
    ## sigh.  also make sure the regular (non-btab) file is there
    if ( $opts{praze_btab} =~ /(.+)\.btab$/ ) {
        my $non_btab = $1;    
    
        if ( ! -e $non_btab ) {
            open(my $infh, ">$non_btab") || die "can't create empty non-btab output file from praze: $!";
            close $infh;
        }

    }    
    
	$praze_btab = new IO::File($opts{praze_btab}) or
		die "Error reading praze btab data $opts{praze_btab}: $!"
		if $opts{praze_btab};
	$out = new IO::File($opts{output}, "w") or
		die "Error writing to $opts{output}: $!"
		if $opts{output};
}

sub map_polypeptide_to_CDS {
	my $map_list = shift;
	open MAP, $map_list or die "Cannot open mapping list for reading: $!\n";
	while (<MAP>) {
		chomp;
		my $line = $_;
		my @l = split(/\t/, $line);
		$map_hash{$l[0]} = $l[1];
	}
	close MAP
}

sub initialize_pvalues {
	my $fname = shift;
	my $map_h = shift;
	my $name;
	my $fh = new IO::File($fname) or
		die "Error reading BLAST btab data $fname: $!";
	while (my $line = <$fh>) {
		chomp $line;
		my $hit = new Blast::BlastHitDataType($line);
		my $query_name = $hit->GetQueryName;	#get polypeptide ID
        my $subject_name = $hit->GetSubjectName;	#get UNIREF ID
        $subject_name =~ s/^lcl\|//;
        print "Subject Name = $subject_name, Polypeptide Name = $query_name\n";

        # if we used a mapping list, assign key as CDS, otherwise assign as polypeptide.
        if (scalar keys %{$map_h} > 0) {
        	$name = $map_h->{$query_name};
		} else {
			$name = $query_name;
		}
		$p_values{$name}{$subject_name} = $hit->GetPValue;
	}
}

sub add_pvalues {
	#Add dummy p-value to praze results that did not have BLAST hits to ensure it never ranks among the best results
	#my $dummy = "1.00e+100";	

	while (my $line = <$praze_btab>) {
		chomp $line;
		my @tokens = split /\t/, $line;
		my $query_name = $tokens[0];
		my $subject_name = $tokens[5];
		my $p_ident = $tokens[10];
		my $p_sim = $tokens[11];
        $subject_name =~ s/^lcl\|//;
        my $p_value;
        if (exists $p_values{$query_name}{$subject_name}){
			#next if $p_ident < $IDENTITY_CUTOFF || $p_sim < $SIMILARITY_CUTOFF;
			$p_value = $p_values{$query_name}{$subject_name};
			print $out "$line\t$p_value\n";
		} #else {
			#$p_value = $dummy;
			#print STDERR "No p-value found for $subject_name on query $query_name.  Inserting dummy p-value\n";
		#}
		#print $out "$line\t$p_value\n";
	}
}
