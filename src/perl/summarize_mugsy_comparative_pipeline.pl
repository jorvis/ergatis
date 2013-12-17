#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use IntervalTree;
use Pod::Usage;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
my %i_trees;
####################################################

my %options;
my $results = GetOptions (\%options,
			  "mugsy_map_file|m=s",
			  "mugsycog_raw_file|c=s",
			  "mugsy_snps_file|s=s",
			  "output_dir|o=s",
			  "log|l=s",
			  "debug|d=s",
			  "help|h"
                          );

&check_options(\%options);

my ($genomes, $molecules) = &parse_map_file( $options{'mugsy_map_file'} );
my $clusters = &parse_cogs( $options{'mugsycog_raw_file'} );
my $snps = &parse_snps( $options{'mugsy_snps_file'} );

foreach my $genome ( keys %{$genomes} ) {
    print "Printing summary file for $genome\n";
    unless( exists( $options{'gene_summary'} ) && $options{'gene_summary'} == 0 ) {
	my $outfile = $options{'output_dir'}."/$genome.gene_summary.txt";
	my $partial_outfile = $options{'output_dir'}."/$genome.partial_gene_summary.txt";
	open( my $ofh, "> $outfile") or die("Could not open $outfile for writing: $!");
	open( my $pofh, "> $partial_outfile") or die("Could not open $partial_outfile for writing: $!");
	&print_gene_summary_info( $genomes->{$genome}, $clusters, $snps, $molecules, $genome, $ofh, $pofh);
	print "\t$outfile\n$partial_outfile\n";
	close($ofh);
	close($pofh);
    }
    #unless( exists( $options{'gene_summary'} ) && $options{'gene_summary'} == 0 ) {
#	my $outfile = $options{'output_dir'}."/$genome.genome_summary.txt";
#	open( my $ofh, "> $outfile") or die("Could not open $outfile for writing: $!");
#	&print_genome_summary_info( $genomes->{$genome}, $clusters, $snps, $molecules, $genome, $ofh );
#	print "\t$outfile\n";
#	close($ofh);
#    }
}

sub print_genome_summary_info {
    my ($gene_lookup, $clusters, $snps, $molecules, $org, $ofh) = @_;

    ## Count some snps
    print $ofh "$org\n";
    my %snp_counts;
    foreach my $molecule( keys %{$snps->{$org}} ) {
    foreach my $gene( keys %{$snps->{$org}->{$molecule}} ) {
    my $c = scalar(keys %{$snps->{$org}->{$molecule}->{$gene}});
    $snp_counts{"by_molecule"}->{$molecule} += $c;
	    $snp_counts{"by_gene"}->{$gene} = $c;
	    $snp_counts{"total"} += $c;
	    my $uniques = &count_unique( $snps->{$org}->{$molecule}->{$gene}, $org );
	    $snp_counts{'unique'} += $uniques;
}
}
    print $ofh "Total SNPs: $snp_counts{'total'} [unique: $snp_counts{'unique'}]\n";
    foreach my $mol ( keys %{$snp_counts{'by_molecule'}} ) {
    print $ofh "\t$mol: $snp_counts{'by_molecule'}->{$mol}\n";
    }		
    print $ofh "\n";

    ## Gene cluster stuff
    my @ordered_genomes = &get_ordered_genomes( $molecules );
    my %cluster_info = ("shared" => [], "core" => [], "unique" => []);
    my $total_genes = 0;
    foreach my $mol ( keys %{$gene_lookup} ) {
    my $genes = $gene_lookup->{$mol};
    $total_genes += scalar(keys %{$genes});
	foreach my $gene ( sort { $genes->{$a}->[0] <=> $genes->{$b}->[0] } keys %{$genes} ) {
    if( exists( $clusters->{$gene} ) ) {
    my $orgs = {};
    foreach my $mol ( keys %{$clusters->{$gene}} ) {
    my $o = $molecules->{$mol};
    $orgs->{$o} = 1;
}

		my $count_orgs_in_cluster = scalar(keys %{$orgs});
		if( $count_orgs_in_cluster == @ordered_genomes ) {
    push(@{$cluster_info{'core'}}, $gene);
		}
		    
} else {
    push(@{$cluster_info{'unique'}}, $gene);
	    }
}
}

    print $ofh "Total Genes: $total_genes\n";

    print $ofh "Core Genes: ".scalar(@{$cluster_info{'core'}})."\n";
    foreach my $gene ( @{$cluster_info{'core'}} ) {
    print $ofh "\t$gene\n";
    }

    print $ofh "Unique Genes: ".scalar(@{$cluster_info{'unique'}})."\n";
    foreach my $gene ( @{$cluster_info{'unique'}} ) {
    print $ofh "\t$gene\n";
    }
}

sub count_unique {
    my ($snp, $org) = @_;
    my $c = 0;
    foreach my $pos ( keys %{$snp} ) {
	my $refbase = $snp->{$pos}->{$org};
	my $same = 0;
	next if( $refbase eq '-' || $refbase eq 'N' );
	foreach my $qorg ( keys %{$snp->{$pos}} ) {
	    next if( $qorg eq $org );
	    if( $snp->{$pos}->{$qorg} eq $refbase ) {
		$same++;
		last;
	    }
	}
	$c++ unless( $same );
    }
    return $c;
}

## Or parsnips maybe?
sub parse_snps {
    my ($snps_file, $molecules) = @_;
    open(IN, "< $snps_file") or die("Couldn't open $snps_file: $!");

    chomp( my $header = <IN> );
    my @h = split(/\t/, $header );
    my %orgi;
    for( my $i = 0; $i < @h; $i++ ) {
	my ($org, $col) = split(/:/, $h[$i]);
	$orgi{$org}->{$col} = $i;
    }

    my $retval = {};
    while( my $line = <IN> ) {
	next if( $line =~ /^\#/ || $line =~ /^\s*$/ );
	chomp($line);

	my $tmp  = {};
	my $snps = [];

	my @c = split(/\t/, $line);
	foreach my $org ( keys %orgi ) {
	    my $mol   = $c[$orgi{$org}->{'mol'}];
	    my $pos   = $c[$orgi{$org}->{'pos'}];
	    my $base  = $c[$orgi{$org}->{'base'}];

	    my $report = (!defined($base) || $base eq '') ? "NA" : $base;
	    $tmp->{$org} = [$report,$pos];

	    if( $mol ) {
		my $itree = $i_trees{$org}->{$mol};
		unless( $itree ) {
			print STDERR "Couldn't find tree for $org, $mol\n";
			next;
		}
		
		# do we overlap with any genes?
		my @overlaps = $itree->searchInterval( $pos, $pos );
		
		if( @overlaps ) {
		    foreach my $ol ( @overlaps ) {
			push(@{$snps}, [$org, $mol, $ol->[2], $pos, $base]);
		    }
		} else {
		    #print "No overlaps!\n";
		    # this is where we would handle intergenic snps, but we don't right now.
		}
	    }

	}
	
	foreach my $snp ( @{$snps} ) {
	    $retval->{$snp->[0]}->{$snp->[1]}->{$snp->[2]}->{$snp->[3]} = $tmp;
	}
	
    }
    close(IN);

    return $retval;
}

sub parse_map_file {
    my ($map_file) = @_;
    my $genomes = {};
    my $molecules = {};

    # 0 - Gene ID
    # 1 - Molecule name
    # 2 - Left
    # 3 - Right
    # 4 - Strand (-1, 1)
    # 5 - BSML polypeptide ID
    # 6 - BSML gene ID
    # 7 - Genome Name
    # 8 - Gene name
    my $process = sub {
	my ($gid, $mol, $le, $ri, $str, $poly, $bgid, $genome, $name) = @{$_[0]};

	# Mugsy changes the organism identifier. So store the genome name as 
	# the same as mugsy. Totally black-boxing how the identifiers are modified
	# so this is by no means exhaustive.
#        print STDERR "$genome is a genome\n";
	#$genome =~ s/[\+\-\:\s]//g;
     #   print STDERR "$genome is a genome with molecule $mol\n";
	$molecules->{$mol} = $genome;
	$genomes->{$genome}->{$mol}->{$gid} = [$le,$ri,$str,$name];
	$i_trees{$genome}->{$mol} = new IntervalTree() unless( exists($i_trees{$genome}->{$mol}) );
	$i_trees{$genome}->{$mol}->addInterval( $gid, $le, $ri );
    };    
    &parse_tab_file( $map_file, $process );

    # build all the interval trees
    for my $genome ( keys %i_trees ) {
	map { $_->buildTree() } values %{$i_trees{$genome}};
    }

    return ($genomes, $molecules);
}

sub parse_tab_file {
    my ($file, $run) = @_;
    open(IN, "< $file") or die("Unable to open $file: $!");
    while( my $line = <IN> ) {
	chomp($line);
	my @cols = split(/\t/, $line);
	$run->(\@cols);
    }
    close(IN);
}

sub parse_cogs {
    my ($cog_raw) = @_;
    my $clusters = {};

    open(IN, "< $cog_raw") or die("Unable to open $cog_raw: $!");

    my $store_tmp = sub {
	my $tmp = shift;
	map { 
	    for my $gene ( @{$tmp->{$_}} ) {
		$clusters->{$gene} = $tmp;
	    }
	} keys %{$tmp};
    };

    my $tmp = {};
    while( my $line = <IN> ) {
	next if( $line =~ /^\#/ || $line =~ /^\s*$/ );
	last if( $line =~ /^Class legend/ );
	chomp( $line );

	if( $line =~ /^>(\S+)/ ) {
	    if( keys %{$tmp} ) {
		$store_tmp->($tmp);
	    }
	    $tmp = {};
	} else {
	    my @c = split(/\t/, $line);
	    my @genes = split(/,/, $c[0]);
	    $tmp->{$c[2]} = \@genes;
	}
	
    }
    $store_tmp->($tmp);
    close(IN);

    return $clusters;
}

sub print_gene_summary_info {
    my ($gene_lookup, $clusters, $snps, $molecules, $org, $ofh, $pofh) = @_;
    my $partial = 0;
    my @ordered_genomes = &get_ordered_genomes( $molecules );

    # print the header
    my $preheader = "Gene Info\t\t\t\t\t\tGenes In Cluster\t".("\t"x (scalar(@ordered_genomes)-1));
    $preheader .= "SNPs";
    my @header = qw(molecule gene start stop strand name);
    push(@header, @ordered_genomes);
    my @other_genomes = grep { !($_ eq $org) } @ordered_genomes;
    push(@header, "Current Genome");
    push(@header, @other_genomes);
    print $ofh $preheader."\n";
    print $ofh join("\t", @header )."\n";
    print $pofh $preheader."\n";
    print $pofh join("\t", @header )."\n";

    foreach my $molecule ( sort keys %{$gene_lookup} ) {

	my $genes = $gene_lookup->{$molecule};
	foreach my $gene ( sort { $genes->{$a}->[0] <=> $genes->{$b}->[0] } keys %{$genes} ) {
	    my @row = ($molecule, $gene, @{$genes->{$gene}});
	    $partial = 0;
	    if( !exists( $clusters->{$gene} ) ) {
		# if it's not in a cluster, push onto the row the correct number
		# of dashes (total genomes minus one)
		push(@row, split(//, "-"x (scalar(@ordered_genomes)) ) );
	    } else {
		my %tmp;

		## Add all the genes to tmp that are in the cluster and
		## key them by organism and molecule.
		foreach my $qmol ( keys %{$clusters->{$gene}} ) {
		    my $qgenome = $molecules->{$qmol};
		    
		    # we should filter out the key gene here.
		    my @other_genes = grep { !($_ eq $gene) } @{$clusters->{$gene}->{$qmol}};
		    if(@{$clusters->{$gene}->{$qmol}} > 1) {
			$partial = 1;
		    }

		    $tmp{$qgenome}->{$qmol} = \@other_genes unless( @other_genes == 0 );
		}
		## For each of the genomes, add genes in cluster
		foreach my $g ( @ordered_genomes ) {
		    if( exists( $tmp{$g} ) ) {
			my @qrow;
			foreach my $qmol ( keys %{$tmp{$g}} ) {
			    push(@qrow, $qmol.":".join("|",@{$tmp{$g}->{$qmol}}) );
			}
			push(@row, join(",",@qrow));
		    } else {
			push(@row, "-");
		    }
		}

		## Are there any SNPs in the gene?
		if( exists( $snps->{$org}->{$molecule}->{$gene} ) ) {
		    my $gene_snps = $snps->{$org}->{$molecule}->{$gene};
		    my %curr_bases;
		    map { $curr_bases{$_} = [] } @ordered_genomes;
		    foreach my $pos ( sort keys %{$gene_snps} ) {
			my $this_snp = $gene_snps->{$pos};
			my $ref_base = $this_snp->{$org}->[0];
			push(@{$curr_bases{$org}}, $pos);

			foreach my $og ( @other_genomes ) {
			    my $qstring = "NA";
			    my $qbase = $this_snp->{$og}->[0];
			    if( $qbase ne "NA" ) {
				my $qpos  = $this_snp->{$og}->[1];
				$qstring = "$qpos:$qbase>$ref_base";
			    }
			    push(@{$curr_bases{$og}}, $qstring);
			}

		    }
		    map { push(@row, join(",", @{$curr_bases{$_}}) ) } ($org, @other_genomes);
		}
	    }
# Separating printing of partial genes in a new file
	    if($partial == 0) {
	    	print $ofh join("\t", @row)."\n";
	    } else {
	    	print $pofh join("\t", @row)."\n"; 
	    }
	}
    }

}

sub get_ordered_genomes {
    my ($molecules) = @_;
    my %seen;
    sort grep { ! $seen{$_}++ } values %{$molecules};
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
    
    foreach my $req ( qw(mugsy_map_file mugsycog_raw_file mugsy_snps_file output_dir) ) {
	&_log($ERROR, "Option $req is required") unless( $opts->{$req} );
    }
}

sub _log {
    my ($level, $msg) = @_;
    if( $level <= $debug ) {
	print STDOUT "$msg\n";
	print $logfh "$msg\n" if( defined( $logfh ) );
    }
    exit(1) if( $level == $ERROR );
}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}


=head1 NAME

create_genome_comparative_summary.pl - Will create a summary file 

=head1 SYNOPSIS

 USAGE: create_genome_comparative_summary.pl
       --input_file=/path/to/some/input.file
       --output=/path/to/transterm.file
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>

B<--output_file,-o>

B<--log,-l>
    Logfile.

B<--debug,-d>
    1,2 or 3. Higher values more verbose.

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
    kgalens@gmail.com

=cut
