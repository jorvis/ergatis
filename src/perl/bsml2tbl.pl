#!/usr/local/bin/perl

use strict;
use warnings;

use XML::Twig;
use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;
use IO::File;

my $in                  = "/dev/stdin";
my $out                 = *STDOUT;
my $output_mrna_feats   = 0;
my $extract_all_ec      = 0;
my $percent_n_cutoff    = 10;
my $asmbl_seq = '';
my $asmbl_length;
## genes with equal or more than this % of Ns will be skiped (warning printed)


parse_options();
convert();

sub print_usage
{
    my $progname = basename($0);
    die << "END";
usage: $progname [--input|-i <input_bsml>]
        [--output|-o <output_tbl>] [--mrna|-m] [--ec_all|-e] [-h]

        -m: export mRNA features [ default - 0 (false) ]
        -e: extract all ec numbers (including incomplete e.g. #.#.#.-)
            [ default - false ]
END
}

sub parse_options
{
    my %opts = ();
    GetOptions(\%opts, "input|i=s", "output|o=s", "mrna|m=i", "percent_n_cutoff|n=i",
            "ec_all|e=i", "help|h");
    print_usage() if $opts{help};
    $in = $opts{input} if $opts{input};
    $out = new IO::File($opts{output}, "w") or
        die "Error writing tbl to $opts{output}: $!"
        if $opts{output};
    $output_mrna_feats = 0 unless $opts{mrna};
    $extract_all_ec = 1 if $opts{ec_all};
    if ( defined $opts{percent_n_cutoff} ) {
        $percent_n_cutoff = $opts{percent_n_cutoff};
    }
}

sub convert
{
    my $twig = new XML::Twig();
    my %feats = ();
    $twig->setTwigRoots({'Feature' => sub { process_feature(\%feats, @_); },
                         'Sequence' => \&process_seq });
    $twig->parsefile($in);
    $twig = new XML::Twig();
    $twig->setTwigRoots({'Feature-group' =>
            sub { process_feature_group(\%feats, @_); } });
    $twig->parsefile($in);
    process_misc_feats(\%feats);
}

sub process_seq
{
    my ($twig, $elt) = @_;
    return if $elt->att('class') ne "assembly";
    print $out ">Features ", $elt->att('id'), "\n";
    
    ## get the sequence length
    if ( $elt->has_child('Seq-data-import') ) {
        open(my $ifh, "<" . $elt->first_child('Seq-data-import')->att('source')) || die "can't open seq-data-import: $!";
        while (<$ifh>) {
            if (! /^\>/ ) {
                $asmbl_seq .= $_;
            }
        }
        $asmbl_seq =~ s/\s+//g;
        $asmbl_length = length $asmbl_seq;
    } else {
        die "error: sequence element has no seq-data-import\n";
    }
}

sub process_feature
{
    my ($feats, $twig, $elt) = @_;
    my $id = $elt->att('id');
    $feats->{$id} = $elt;
    $twig->purge();
}

sub process_feature_group
{
    my ($feats, $twig, $elt) = @_;
    my @exons = ();
    my $cds = undef;
    my $class = undef;
    my $gene = undef;
    my $polypeptide = undef;
    my $transcript = undef;
    my @repeats = ();
    foreach my $feat_member ($elt->children('Feature-group-member')) {
        my $feat_type = $feat_member->att('feature-type');
        my $featref = $feat_member->att('featref');
        my $feat = $feats->{$featref};
        if ($feat_type eq "exon") {
            push @exons, $feat;
        }
        elsif ($feat_type eq "transcript") {
            $class = "transcript";
            $transcript = $feat;
        }
        elsif ($feat_type eq "tRNA") {
            $class = "tRNA";
            $transcript = $feat;
        }
        elsif ($feat_type eq "ncRNA") {
            $class = "ncRNA";
        }
        elsif ($feat_type eq "gene") {
            $gene = $feat;
        }
        elsif ($feat_type eq "CDS") {
            $cds = $feat;
        }
        elsif ($feat_type eq "polypeptide") {
            $polypeptide = $feat;
        }
        elsif ($feat_type eq "repeat_region") {
            push @repeats, $feat;
        }
    }

    process_gene($gene, $transcript) if $gene;
    print_feat(\@repeats, "repeat_region") if scalar(@repeats);

    if ($class eq "transcript") {
        my $product = $transcript->first_child('Attribute[@name="gene_product_name"]');

	my $productContent;
	if (defined($product)){
	    $productContent = $product->att('content');
	}

	my $protein_id;

	if (defined($polypeptide)){

	    $protein_id = $polypeptide->att('id');

	    if (!defined($product)){
		## The product was not associated with the transcript.
		## Check the polypeptide.
		my $product = $polypeptide->first_child('Attribute[@name="gene_product_name"]');
		if (defined($product)){
		    $productContent = $product->att('content');
		}
	    }
	}

	if (!defined($productContent)){
	    die "productContent was not defined";
	}

        my %attrs = (   transcript_id  => $transcript->att('id'),
                        protein_id => $protein_id,
                        product => $productContent
                    );

        my @ec_numbers = ();
        foreach my $att_list( $transcript->children('Attribute-list') ) {
            my $ec_number = $att_list->first_child('Attribute[@name="EC"]');
            push @ec_numbers, $ec_number if $ec_number;
        }
        foreach my $ec_number (@ec_numbers) {
            my $ec = $ec_number->att('content');
            if ($extract_all_ec || $ec !~ /-/) {
                push @{$attrs{EC_number}}, $ec_number->att('content');
            }
        }
        print_feat(\@exons, "mRNA", \%attrs) if $output_mrna_feats;
        process_cds(\@exons, $cds, \%attrs);
    }
    elsif ($class eq "tRNA") {
        my %attrs = ();
        
        if ( $transcript->has_child('Attribute[@name="gene_product_name"]') ) {
            my $product = $transcript->first_child('Attribute[@name="gene_product_name"]');
            %attrs = ( product => $product->att('content') );
        }
        
        print_feat(\@exons, "tRNA", \%attrs);
    }
    elsif ($class eq "ncRNA") {
        print_feat(\@exons, "ncRNA");
    }
}

sub process_gene
{
    my ($gene, $transcript) = @_;
    my $locus = $gene->first_child('Cross-reference[@database="TIGR_moore"]');
    my $locus_tag = $gene->first_child('Attribute[@name="locus_tag"]');
    my $cross_reference = $gene->first_child('Cross-reference[@identifier-type="locus"]');
    my %attrs = ();
    $attrs{'locus_tag'} = $locus_tag->att('content') if $locus_tag;
    $attrs{'locus_tag'} = $locus->att('identifier') if $locus;
    $attrs{'locus_tag'} = $cross_reference->att('identifier') if $cross_reference;
    if ($transcript) {
        my $symbol =
            $transcript->first_child('Attribute[@name="gene_symbol"]');
        $attrs{'gene'} = $symbol->att('content') if $symbol;
    }
    print_feat([$gene], "gene", \%attrs);
}

sub process_cds
{
    my ($exons, $cds, $attrs) = @_;
    my @cdss;

    if (!defined($cds)){
	die "cds was not defined";
    }

    my $cds_loc = $cds->first_child('Interval-loc');
    my $cds_from = $cds_loc->att('startpos');
    my $cds_to = $cds_loc->att('endpos');
    foreach my $exon (@{$exons}) {
        my $exon_loc = $exon->first_child('Interval-loc');
        my $exon_from = $exon_loc->att('startpos');
        my $exon_to = $exon_loc->att('endpos');
        next if ($exon_to < $cds_from || $exon_from > $cds_to);
        my $elt = $exon->copy();
        if ($exon_from < $cds_from) {
            $elt->first_child('Interval-loc')->set_att('startpos', $cds_from);
        }
        if ($exon_to > $cds_to) {
            $elt->first_child('Interval-loc')->set_att('endpos', $cds_to);
        }
        push @cdss, $elt;
    }
    my $cog = $cds->first_child('Attribute[@name="top_cog_hit"]');
    if ($cog) {
        $attrs->{note} = $cog->att('content');
    }
    print_feat(\@cdss, "CDS", $attrs);
}

sub process_misc_feats
{
    my ($feats) = @_;
    while (my ($key, $val) = each %{$feats}) {
        my $class = $val->att('class');
        if ($class eq "repeat_region") {
            print_feat([$val], "repeat_region");
        }
    }
}

sub print_feat
{
    my ($feats, $class, $attrs) = @_;
    for (my $i = 0; $i < scalar(@{$feats}); ++$i) {
        my $feat = $feats->[$i];
        my $loc = $feat->first_child('Interval-loc');
        my $from = $loc->att('startpos') + 1;
        my $to = $loc->att('endpos');
        my $codon_start = 0;
        if ($loc->att('complement') == 1) {
            ($from, $to) = ($to, $from);
        }
        my $tmp_from = $from;
        my $tmp_to = $to;
        
        ## check to see if we're off the end of the sequence here
        #   I'm addressing the other cases as I see them.
        if ( $from > $asmbl_length && $from > $to ) {
            $tmp_from = "<$asmbl_length";
            $attrs->{codon_start} = 4 - ( $from - $asmbl_length );
        } elsif( $to > $asmbl_length && $to > $from ) {
            $tmp_to = ">$asmbl_length";
        }

        ## Deal with negative coordinates
        if ($from <= 0 && $from < $to) {
            if ($class eq "CDS") {
                my $frame = abs($from - 1) % 3;
                if ($frame == 1) {
                    $attrs->{codon_start} = 3;
                }
                elsif ($frame == 2) {
                    $attrs->{codon_start} = 2;
                }
            }
            $tmp_from = "<1";
        }

        if($to <= 0) {
            $tmp_to = ">1";
        }
            
       
        $from = $tmp_from;
        $to = $tmp_to;
        
        if ( $class eq 'gene' || $class eq 'CDS' || $class eq 'mRNA' ) {
            ## bail if this sequence has too many Ns
            my $feat_seq = substr($asmbl_seq, $loc->att('startpos'), $loc->att('endpos') - $loc->att('startpos') );
            my $n_count = ( $feat_seq =~ tr/Nn// );
            my $n_perc = ($n_count / length($feat_seq)) * 100;

            if ( $n_perc >= $percent_n_cutoff ) {
                print STDERR "warning: skipping $class because it has too many Ns (", sprintf("%.1f",$n_perc), "\%)\n";
                return;
            }
        }
        
        print $out "$from\t$to";
        print $out "\t$class" if $i == 0;
        print $out "\n";
    }

    while (my ($key, $val) = each %{$attrs}) {
        if (ref($val) eq "ARRAY") {
            foreach my $v (@{$val}) {
                print_attribute($out, $key, $v);
            }
        }
        else {
            print_attribute($out, $key, $val);
        }
    }
}

sub print_attribute
{
    my ($out, $key, $val) = @_;
    print $out "\t\t\t$key\t$val\n";
}
